!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 593 $)
!
!  MODULE: m_constraints, m_velocity_verlet, m_routines, m_mtrandom,
!          m_constraints_data,
!  SUBROUINE: constrained_dynamics_initialize, constrained_dynamics_dump,
!             constrained_dynamics_finalize
!  FUNCTION: f_getBoolValue
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
!!!!BRANCH_P ORG_Parallel
#ifndef DISABLE_CONSTRAINTS
! ==============================================================================
module  m_variables_for_atoms
  use m_Const_Parameters, only: DP,PAI, INITIAL
  use m_Control_Parameters, only : icond
  use m_Ionic_System, only : natm,altv,cps,cpd_l,imdtyp,amion,ityp,lattice_system_from_m_CS_SG
  implicit none

  integer :: Nfree
  integer :: Nreservoir

  real(DP)                      :: Cell(3,3), Cell_inv(3,3), Reciprocal(3,3)
  real(DP), allocatable         :: x_nopbc(:), y_nopbc(:), z_nopbc(:)

  real(DP)                      :: a,b,c,alpha,beta,gamma
  logical,  allocatable         :: mobile(:),dummy_atom(:)

  contains

  subroutine init()
    integer :: i
    mobile = .true.
    Nfree = 0
    do i=1,natm
      if (imdtyp(i)==0) then
        mobile(i)=.false.
      else
        Nfree = Nfree+3
      endif
    enddo
    if(Nfree>0) Nfree = Nfree-3 ! remove overall translation
    dummy_atom  = .false.
    Cell = altv
    call m_vfa_init_nopbc()
  end subroutine init

  subroutine m_vfa_init_nopbc()
    integer :: i
    do i=1,natm
      x_nopbc(i)=cps(i,1)
      y_nopbc(i)=cps(i,2)
      z_nopbc(i)=cps(i,3)
    enddo
  end subroutine m_vfa_init_nopbc

  subroutine m_vfa_alloc_atoms()
    Nfree = 3*(natm-1)
    Nreservoir = Nfree+1
    allocate(mobile(natm))
    allocate(dummy_atom(natm))
    allocate( x_nopbc(natm),y_nopbc(natm),z_nopbc(natm) )
    call init()
  end subroutine m_vfa_alloc_atoms

  subroutine m_vfa_dealloc()
    deallocate(mobile)
    deallocate(dummy_atom)
    deallocate(x_nopbc,y_nopbc,z_nopbc)
  end subroutine m_vfa_dealloc

  real(DP) function get_total_mass()
    integer :: i
    real(DP) :: tm
    tm=0.d0
    do i=1,natm
      tm = tm + amion(ityp(i))
    enddo
    get_total_mass = tm
  end function get_total_mass

end module m_variables_for_atoms

module m_variables_for_dynamics
  use m_Const_Parameters, only : CONST_kB,VERLET,T_CONTROL,QUENCHED_MD,DP,DAMPED_MD,VELOCITY_SCALING
  use m_Control_Parameters, only : printable,imdalg,forccr
  use m_Ionic_System, only : tkb,qmass,natm,t_ctrl_method
  use m_Files, only : nfout
  implicit none

  integer :: curr_md_step = 0
  real(DP)              :: max_force  = 1.0d-3
  real(DP)              :: Temperature, Temperature_atom
  real(DP)              :: mass_thermo

  integer, private :: allocatestatus
  integer :: mdstep_for_this_run
  real(DP) :: forcmx

  real(DP), allocatable, dimension(:,:) :: v_old

  real(DP) :: H_NVT
  integer :: n_md_step

  contains

  logical function is_md()
    is_md=.false.
    is_md=imdalg==VERLET .or. imdalg==T_CONTROL .or. imdalg==VELOCITY_SCALING
  end function is_md

  subroutine m_vfd_init()
    implicit none
    integer :: i
    max_force = forccr
    curr_md_step=0
    mdstep_for_this_run=0
    Temperature_atom = 0.d0
    if(allocated(tkb))   Temperature_atom = tkb(1)/CONST_kB
    if(allocated(qmass)) mass_thermo = qmass(1)
    if(t_ctrl_method==VELOCITY_SCALING) imdalg=VELOCITY_SCALING
  end subroutine m_vfd_init

  real(DP) function get_mass_thermo_default()
    get_mass_thermo_default = 10.d0*dble(natm)/8.d0 !8-ko de 10; toriaezu
  end function get_mass_thermo_default

end module m_variables_for_dynamics

module m_routines

use m_Const_Parameters,   only : DP,PAI,INITIAL,CONTINUATION
use m_Control_Parameters, only : icond, printable
use m_Files, only : nfout
use m_Parallelization, only : mype,sw_wdir,workdir

implicit none

contains

!character(len=256) function to_string(i,keta)
!  integer, intent(in) :: i
!  integer, intent(in) :: keta
!  character :: string*256
!  character :: ch*256
!  integer :: iketa
!  integer :: itmp
!  iketa = keta
!  itmp=int(log10(real(mype_conf)))+1
!  if(itmp>keta)iketa=itmp
!  write(ch,*) iketa
!  write(string,'(i'//trim(adjustl(ch))//'.'//trim(adjustl(ch))//')') i
!  to_string = trim(adjustl(string))
!end function to_string

subroutine cross_product(a,b,c)
  implicit none
  real(DP), dimension(3), intent(in) :: a,b
  real(DP), dimension(3), intent(out)::  c
  c=0.d0
  c(1) = a(2)*b(3)-a(3)*b(2)
  c(2) = a(3)*b(1)-a(1)*b(3)
  c(3) = a(1)*b(2)-a(2)*b(1)
end subroutine cross_product

real(DP) function get_norm(vec)
  implicit none
  real(DP), dimension(3),intent(in) :: vec
  get_norm = dsqrt(vec(1)**2+vec(2)**2+vec(3)**2)
end function get_norm

real(DP) function get_distance_nopbc(j,i,exclude)
  integer, intent(in) :: j,i
  integer,dimension(3),intent(in),optional :: exclude
  integer :: k
  real(DP), dimension(3) :: tmprealx
  real(DP) :: tmpdistance
  integer,dimension(3) :: excl
  excl = 0
  if(present(exclude)) excl = exclude
  tmpdistance = 0.d0
  tmprealx    = 0.d0
  call get_vector_nopbc(j,i,tmprealx)
  do k=1,3
    if(excl(k)==1) cycle
    tmpdistance=tmpdistance+tmprealx(k)**2
  enddo
  get_distance_nopbc = dsqrt(tmpdistance)
end function get_distance_nopbc

subroutine get_vector_nopbc(j,i,ret)
  use m_variables_for_atoms, only : x_nopbc,y_nopbc,z_nopbc
  integer, intent(in) :: j,i
  real(DP), dimension(3),intent(out) :: ret
  integer :: k,l

  ret(1) = x_nopbc(j)-x_nopbc(i)
  ret(2) = y_nopbc(j)-y_nopbc(i)
  ret(3) = z_nopbc(j)-z_nopbc(i)
end subroutine get_vector_nopbc

integer function get_unused_unitnumber()
  integer :: i
  logical :: op
  do i=500,9999
    inquire(unit=i,opened=op)
    if( .not.op ) then
      get_unused_unitnumber = i
      return
    endif
  enddo
  call phase_error_with_msg(nfout,'  insufficient file handle',__LINE__,__FILE__)
end function get_unused_unitnumber

subroutine close_all_opened_files()
  integer :: i
  logical :: op
  if(mype==0)then
    do i=500,9999
      inquire(unit=i,opened=op)
      if ( op ) then
         call flush(i)
         close(i)
      endif
    enddo
  endif
end subroutine close_all_opened_files

subroutine smart_output_open(unit_id,fname,force_replace,force_append)
  integer, intent(in) :: unit_id
  character(len=*), intent(in) :: fname
  logical, optional, intent(in) :: force_replace
  logical, optional, intent(in) :: force_append
  logical :: exi
  logical :: repl
  logical :: appe
  logical :: op
  integer :: num
  if (mype/=0) return
  repl=.false.
  appe=.false.
  if(present(force_replace))repl=force_replace
  if(present(force_append)) appe=force_append

  inquire(file=fname,exist=exi,opened=op,number=num)
  if (op) close(unit_id)

  if(appe) then
    if(sw_wdir==1)then
      open(unit_id,&
    & file=trim(workdir)//fname,status='old',position='append')
    else
      open(unit_id,&
    & file=fname,status='old',position='append')
    endif
    if(printable) write(nfout,'(a,a32,a)') ' opened output file: ',fname,  ', status="old",     position="append".'
  endif
  if(repl) then
    if(sw_wdir==1)then
      open(unit_id,&
    & file=trim(workdir)//fname,status='replace',position='rewind')
    else
      open(unit_id,&
    & file=fname,status='replace',position='rewind')
    endif
    if(printable) write(nfout,'(a,a32,a)') ' opened output file: ',fname,  ', status="replace", position="rewind".'
  endif
  if(appe .or. repl) return

  if(icond==INITIAL) then
    if(sw_wdir==1)then
      open(unit_id,&
    & file=trim(workdir)//fname,status='replace',position='rewind')
    else
      open(unit_id,&
    & file=fname,status='replace',position='rewind')
    endif
    if(printable) write(nfout,'(a,a32,a)') ' opened output file: ',fname,  ', status="replace", position="rewind".'
  else if(icond==CONTINUATION)then
    if(exi)then
      if(sw_wdir==1)then
        open(unit_id,&
    &   file=trim(workdir)//fname,status='old',position='append')
      else
        open(unit_id,&
    &   file=fname,status='old',position='append')
      endif
      if(printable) write(nfout,'(a,a32,a)') ' opened output file: ',fname,', status="old",     position="append".'
    else
      if(sw_wdir==1)then
        open(unit_id,&
    &   file=trim(workdir)//fname,status='replace',position='rewind')
      else
        open(unit_id,&
    &   file=fname,status='replace',position='rewind')
      endif
      if(printable) write(nfout,'(a,a32,a)') ' opened output file: ',fname,', status="replace", position="rewind".'
    endif
  endif
end subroutine smart_output_open

subroutine smart_output_open0(unit_id,fname)
  integer, intent(in) :: unit_id
  character(len=*), intent(in) :: fname
  logical :: exi
  logical :: op
  integer :: num
  inquire(file=fname,exist=exi,opened=op,number=num)
  if (op) close(num)
  if(exi)then
    open(unit_id,&
  &   file=fname,status='old',position='append')
    if(printable) write(nfout,'(a,a32,a)') ' opened output file: ',fname,', status="old",     position="append".'
  else
    open(unit_id,&
  &   file=fname,status='replace',position='rewind')
    if(printable) write(nfout,'(a,a32,a)') ' opened output file: ',fname,', status="replace", position="rewind".'
  endif
end subroutine smart_output_open0

end module m_routines

!
! random number generator using MT algorithm by MATSUMOTO Makoto.
!
module m_mtrandom
  implicit none
  integer, parameter :: defaultsd = 4357
  integer, parameter :: N = 624, MT_N1 = N + 1
  integer, save, dimension(0:N-1) :: mt
  integer, save                   :: mti = MT_N1
  integer, private :: sd = defaultsd

contains
  subroutine init_seed()
    call sgrnd(sd)
  end subroutine init_seed

  subroutine set_seed(ssd)
    integer, intent(in) :: ssd
    sd = ssd
  end subroutine set_seed

  integer function get_seed()
    get_seed = sd
  end function get_seed

  subroutine sgrnd(seed)
    implicit none
    integer, intent(in) :: seed

    mt(0) = iand(seed,-1)
    do mti=1,N-1
       mt(mti) = iand(69069 * mt(mti-1),-1)
    enddo
    return
  end subroutine sgrnd

  real(8) function grnd()
    implicit none
    integer, parameter :: M = 397, MATA  = -1727483681
    integer, parameter :: LMASK =  2147483647
    integer, parameter :: UMASK = -LMASK - 1
    integer, parameter :: TMASKB= -1658038656, TMASKC= -272236544
    integer :: y, kk
    integer, dimension(0:1) :: mag01
    data mag01/0, MATA/
    save mag01

    if(mti.ge.N) then
       if(mti.eq.N+1) then
          call sgrnd( sd )
       endif
       do kk=0,N-M-1
          y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
          mt(kk)=ieor(ieor(mt(kk+M),ishft(y,-1)),mag01(iand(y,1)))
       enddo
       do kk=N-M,N-2
          y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
          mt(kk)=ieor(ieor(mt(kk+(M-N)),ishft(y,-1)),mag01(iand(y,1)))
       enddo
       y=ior(iand(mt(N-1),UMASK),iand(mt(0),LMASK))
       mt(N-1)=ieor(ieor(mt(M-1),ishft(y,-1)),mag01(iand(y,1)))
       mti = 0
    endif

    y=mt(mti)
    mti = mti + 1
    y=ieor(y,tshftu(y))
    y=ieor(y,iand(tshfts(y),TMASKB))
    y=ieor(y,iand(tshftt(y),TMASKC))
    y=ieor(y,tshftl(y))

    if(y .lt. 0) then
       grnd=(dble(y)+2.0d0**32)/(2.0d0**32-1.0d0)
    else
       grnd=dble(y)/(2.0d0**32-1.0d0)
    endif
    return
  end function grnd

  integer function tshftu(i)
    integer, intent(in) :: i
    tshftu = ishft(i,-11)
    return
  end function tshftu

  integer function tshfts(i)
    integer, intent(in) :: i
    tshfts = ishft(i,7)
    return
  end function tshfts

  integer function tshftt(i)
    integer, intent(in) :: i
    tshftt = ishft(i,15)
    return
  end function tshftt

  integer function tshftl(i)
    integer, intent(in) :: i
    tshftl = ishft(i,-18)
    return
  end function tshftl
end module m_mtrandom

!! TODO
!! implement pos-vec constraint
!! implement type 2 coordination number constraint
module m_constraints_data
use m_Const_Parameters, only : DP, PAI, CARTS, PUCV, FMAXVALLEN, LOWER, FMAXUNITLEN
use m_Files, only : nfout, F_CNTN, F_CNTN_bak, m_Files_reopen_nfcntn
use m_Control_Parameters, only : printable,ipri,iprimd
use m_variables_for_atoms, only : x_nopbc,y_nopbc,z_nopbc,Cell
use m_Ionic_System, only : cps,natm,amion,ityp,pos,cps_in,pos_in
!use m_routines, only : get_distance_nopbc, get_vector_nopbc, get_norm, cross_product
use m_routines
implicit none

integer, parameter :: NO_UNIT = -1
integer, parameter :: LENGTH  =  0
integer, parameter :: ANGLE   =  1

! specify constrainable coordinate type
integer, parameter :: DISTANCE_FROM_POS  = -1
integer, parameter :: BOND_LENGTH        =  0
integer, parameter :: BOND_ANGLE         =  1
integer, parameter :: DIHEDRAL_ANGLE     =  2
integer, parameter :: PLANE              =  3
integer, parameter :: CENTER_OF_MASS     =  4
integer, parameter :: BOND_LENGTH_DIFF   =  5
integer, parameter :: COORDNATION_NUMBER =  6
#ifdef USER_DEFINED_CONSTRAINT
integer, parameter :: USER_DEFINED       =  7
#endif
integer, parameter :: BOND_ANGLE_DIFF    =  8
integer, parameter :: BOND_LENGTH_SUM    =  9
integer, parameter :: DISTANCE_FROM_COM  = 10
integer, parameter :: DISTANCE_FROM_REF  = 11

integer, parameter :: REAC_COORD_VIA_FILE  = 0
integer, parameter :: REAC_COORD_VIA_INPUT = 1

real(DP), parameter :: very_small = 1.d-14

!tags for parsing the input file
character(len("constrainable")), parameter :: tag_constrainable = "constrainable"
character(len("reaction_coordinate")), parameter :: tag_reaction_coordinate = "reaction_coordinate"
character(len("sw_reaction_coordinate")),parameter :: tag_sw_reaction_coordinate = "sw_reaction_coordinate"
character(len("init_value")), parameter :: tag_init_value = "init_value"
character(len("final_value")), parameter :: tag_final_value = "final_value"
character(len("increment")), parameter :: tag_increment = "increment"
character(len("mobile")), parameter :: tag_mobile = "mobile"
character(len("monitor")), parameter :: tag_monitor = "monitor"
character(len("atom")), parameter :: tag_atom ="atom"
character(len("atom_group")), parameter :: tag_atom_group ="atom_group"
character(len("id")), parameter :: tag_id ="id"
character(len("no")), parameter :: tag_no ="no"
character(len("type")), parameter :: tag_type ="type"

character(len("normalization_factor")), parameter :: tag_normalization_factor='normalization_factor'

character(len("pos_vec")), parameter :: tag_pos_vec = "pos_vec"
character(len("distance_from_pos")), parameter :: tag_distance_from_pos = "distance_from_pos"
character(len("distance_from_com")), parameter :: tag_distance_from_com = "distance_from_com"
character(len("distance_from_ref")), parameter :: tag_distance_from_ref = "distance_from_ref"
character(len("bond_length")), parameter :: tag_bond_length = "bond_length"
character(len("bond_angle")), parameter :: tag_bond_angle = "bond_angle"
character(len("dihedral_angle")), parameter :: tag_dihedral_angle = "dihedral_angle"
character(len("plane")), parameter :: tag_plane = "plane"
character(len("center_of_mass")), parameter :: tag_center_of_mass = "center_of_mass"
character(len("bond_length_diff")), parameter :: tag_bond_length_diff = "bond_length_diff"
character(len("bond_length_sum")), parameter :: tag_bond_length_sum = "bond_length_sum"
character(len("coordination_number")), parameter :: tag_coordination_number = "coordination_number"
#ifdef USER_DEFINED_CONSTRAINT
character(len("user_defined")), parameter :: tag_user_defined = "user_defined"
#endif
character(len("bond_angle_diff")), parameter :: tag_bond_angle_diff = "bond_angle_diff"

character(len("reac_coord_generation")), parameter :: tag_reac_coord_generation="reac_coord_generation"
character(len("via_file")), parameter :: tag_via_file="via_file"
character(len("via_input")), parameter :: tag_via_input="via_input"

character(len("reference_structure")), parameter :: tag_reference_structure="reference_structure"
character(len("rx")), parameter :: tag_rx="rx"
character(len("ry")), parameter :: tag_ry="ry"
character(len("rz")), parameter :: tag_rz="rz"
character(len("coordinate_system")), parameter :: tag_coordinate_system="coordinate_system"
character(len("cartesian")), parameter :: tag_cartesian="cartesian"
character(len("xyz")), parameter :: tag_xyz="xyz"
character(len("atoms")), parameter :: tag_atoms="atoms"

!type for 'constrainable coordinates'
type constrainable_coords_t
  character(30)  :: nam
  character(256) :: descri
  character(10)  :: unit_name
  integer :: typ
  integer, pointer :: associated_atoms(:)
  integer :: n_associated_atoms
  integer :: nfree_decrement

  integer :: ndim
  real(DP), pointer :: value(:)
  real(DP) :: normalization_factor
  real(DP), pointer :: sigma(:)
  real(DP), pointer :: dsigma(:,:,:)
  real(DP), pointer :: dsigma_old(:,:,:)
  logical :: mobile
  logical :: monitor

  real(DP) :: det_metric

  real(DP), pointer, dimension(:) :: lambda

  real(DP) :: rmix_coords
  real(DP) :: rmix_velocity

  logical :: is_reaction_coordinate
  integer :: n_reaction_coords
  real(DP), pointer, dimension(:,:) :: reaction_coords
  logical,  pointer, dimension(:) :: finished

  integer :: reac_id
  integer :: reac_coord_count
  real(DP), pointer, dimension(:) :: origvalue
  integer :: reac_direction

  logical :: pbc

  real(DP), pointer :: refx(:)
  real(DP), pointer :: refy(:)
  real(DP), pointer :: refz(:)

  integer           :: ngroup
  integer, pointer  :: atom_groups(:,:)
  integer, pointer  :: natm_per_group(:)
end type constrainable_coords_t

integer :: nvariable
integer, allocatable, dimension(:) :: nvariable_map
logical, allocatable, dimension(:) :: reac_set_finished
integer :: reac_coord_generation=REAC_COORD_VIA_INPUT

contains
  subroutine alloc_constrainable_coord(constrainable_coord, allocate_map)
    type(constrainable_coords_t) :: constrainable_coord
    logical, intent(in), optional :: allocate_map
    integer :: nassoc,ndi
    integer :: i,j,icount
    logical :: amap
    amap = .true.
    if(present(allocate_map)) amap = allocate_map
    nassoc = constrainable_coord%n_associated_atoms
    ndi = constrainable_coord%ndim
    allocate(constrainable_coord%value(ndi));constrainable_coord%value=0.d0
    allocate(constrainable_coord%sigma(ndi));constrainable_coord%sigma=0.d0
    allocate(constrainable_coord%lambda(ndi));constrainable_coord%lambda=0.d0
    if(amap) then
      allocate(constrainable_coord%associated_atoms(nassoc));constrainable_coord%associated_atoms=0
    endif
    allocate(constrainable_coord%dsigma(ndi,nassoc,3));constrainable_coord%dsigma=0.d0
    allocate(constrainable_coord%dsigma_old(ndi,nassoc,3));constrainable_coord%dsigma_old=0.d0
    if(constrainable_coord%mobile) constrainable_coord%nfree_decrement=0
    if(constrainable_coord%ngroup>0) then
      icount=1
      do i=1,constrainable_coord%ngroup
        do j=1,constrainable_coord%natm_per_group(i)
          constrainable_coord%associated_atoms(icount) = constrainable_coord%atom_groups(j,i)
          icount = icount+1
        enddo
      enddo
    endif
  end subroutine alloc_constrainable_coord

  subroutine dealloc_constrainable_coord(constrainable_coord)
    type(constrainable_coords_t) :: constrainable_coord
    deallocate(constrainable_coord%value)
    deallocate(constrainable_coord%sigma)
    deallocate(constrainable_coord%lambda)
    deallocate(constrainable_coord%dsigma)
    deallocate(constrainable_coord%dsigma_old)
    deallocate(constrainable_coord%associated_atoms)
    deallocate(constrainable_coord%origvalue)
    if(constrainable_coord%ngroup>0) then
      deallocate(constrainable_coord%natm_per_group)
      deallocate(constrainable_coord%atom_groups)
    endif
    !if(constrainable_coord%is_reaction_coordinate) then
    !  deallocate(constrainable_coord%reaction_coords)
    !  deallocate(constrainable_coord%finished)
    !endif
  end subroutine dealloc_constrainable_coord

  subroutine parse_atom_groups(constrainable_coord)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    character(256) countstr,tmpstr
    integer :: i,j,icount,iret,iiret,nmax,ii, ntot
    integer :: f_selectNextTableLine, f_selectFirstTableLine, f_selectBlock, f_getIntValue, f_selectParentBlock

    icount = 1

    do while(.true.)
      write(countstr,*) icount
      tmpstr = trim(adjustl(tag_atom//trim(adjustl(countstr))))
      iiret = f_selectBlock(tmpstr)
      if(iiret/=0)then
        tmpstr = trim(adjustl(tag_atom_group//trim(adjustl(countstr))))
        iiret = f_selectBlock(tmpstr)
      endif
      if(iiret/=0)then
        tmpstr = trim(adjustl(tag_atom//trim(adjustl(countstr))))
        ii = f_getIntValue(tmpstr,iret)
      endif
      if(iiret/=0 .and. ii/=0) then
        constrainable_coord%ngroup = icount-1
        exit
      endif
      icount = icount+1
      if(iiret==0) iiret = f_selectParentBlock()
    enddo
    if(constrainable_coord%ngroup==0) return
    allocate(constrainable_coord%natm_per_group(constrainable_coord%ngroup))
    do i=1,constrainable_coord%ngroup
      write(countstr,*) i
      tmpstr = trim(adjustl(tag_atom//trim(adjustl(countstr))))
      iiret = f_selectBlock(tmpstr)
      if(iiret/=0) then
        tmpstr = trim(adjustl(tag_atom_group//trim(adjustl(countstr))))
        iiret = f_selectBlock(tmpstr)
      endif
      icount = 1
      if(iiret==0) then
        do while(.true.)
          if (icount == 1) then
             if(f_selectFirstTableLine() /= 0) then
                exit
             end if
          else
             if(f_selectNextTableLine() /= 0) then
                constrainable_coord%natm_per_group(i) = icount-1
                exit
             end if
          end if
          icount = icount+1
        enddo
        iiret = f_selectParentBlock()
      else
        tmpstr = trim(adjustl(tag_atom//trim(adjustl(countstr))))
        if(f_getIntValue(tmpstr,iret)==0) constrainable_coord%natm_per_group(i) = 1
      endif
    enddo
    nmax = maxval(constrainable_coord%natm_per_group)
    allocate(constrainable_coord%atom_groups(nmax,constrainable_coord%ngroup))
    do i=1,constrainable_coord%ngroup
      write(countstr,*) i
      tmpstr = trim(adjustl(tag_atom//trim(adjustl(countstr))))
      iiret = f_selectBlock(tmpstr)
      if(iiret/=0) then
        tmpstr = trim(adjustl(tag_atom_group//trim(adjustl(countstr))))
        iiret = f_selectBlock(tmpstr)
      endif
      icount = 1
      if(iiret==0) then
        do while(.true.)
          if (icount == 1) then
             if(f_selectFirstTableLine() /= 0) then
                exit
             end if
          else
             if(f_selectNextTableLine() /= 0) then
                exit
             end if
          end if
          if(f_getIntValue(tag_id,iret)==0) then
            constrainable_coord%atom_groups(icount,i) = iret
          endif
          if(f_getIntValue(tag_no,iret)==0) then
            constrainable_coord%atom_groups(icount,i) = iret
          endif
          icount = icount+1
        enddo
        iiret = f_selectParentBlock()
      else
        tmpstr = trim(adjustl(tag_atom//trim(adjustl(countstr))))
        if(f_getIntValue(tmpstr,iret)==0) constrainable_coord%atom_groups(1,i) = iret
      endif
    enddo

    ntot = 0
    do i=1,constrainable_coord%ngroup
      ntot = ntot+constrainable_coord%natm_per_group(i)
    enddo
    constrainable_coord%n_associated_atoms = ntot

    if(printable) then
      write(nfout,'(a,i8)') '!** number of groups defined ',constrainable_coord%ngroup
      do i=1,constrainable_coord%ngroup
        write(nfout,'(a,i8,a,i8)')  '!** number of atoms defined in group ',i,':',constrainable_coord%natm_per_group(i)
        write(nfout,'(a,i8)')       '!**           atoms defined in group ',i
        write(nfout,'(8i8)')        (constrainable_coord%atom_groups(j,i),j=1,constrainable_coord%natm_per_group(i))
      enddo
    endif
  end subroutine parse_atom_groups

  subroutine parse_associated_atoms(constrainable_coord)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    integer :: f_getIntValue
    character(256) tmpstr
    integer :: iret,iiret,i
    call parse_atom_groups(constrainable_coord)
!    if (constrainable_coord%n_associated_atoms==0)then
!     i=1
!     do while(.true.)
!       write(tmpstr,*) i
!       tmpstr = trim(adjustl(tag_atom//trim(adjustl(tmpstr))))
!       iiret = f_getIntValue(tmpstr,iret)
!       if(iiret/=0)then
!         constrainable_coord%n_associated_atoms = i-1
!         exit
!       endif
!       i=i+1
!     enddo
!     if (constrainable_coord%n_associated_atoms==0) then
!       write(0,*) 'you must specify at least one atom in order to define a '//trim(constrainable_coord%nam)//' constraint'
!       call phase_error_with_msg(nfout, &
!       'you must specify at least one atom in order to define a '//trim(constrainable_coord%nam)//' constraint' &
!       ,__LINE__,__FILE__)
!     endif
!     if(constrainable_coord%ndim==0) constrainable_coord%ndim = constrainable_coord%n_associated_atoms
!     call alloc_constrainable_coord(constrainable_coord)
!   endif
!   do i=1,constrainable_coord%n_associated_atoms
!     write(tmpstr,*) i
!     tmpstr = trim(adjustl(tag_atom//trim(adjustl(tmpstr))))
!     iiret = f_getIntValue(tmpstr,iret)
!     if(iiret/=0) then
!       if(printable)write(nfout,'(a)') &
!   &   'you must specify '//trim(tmpstr)//' in order to define a '//trim(constrainable_coord%nam)//' constraint'
!       call phase_error_with_msg(nfout,&
!       'you must specify '//trim(tmpstr)//' in order to define a '//trim(constrainable_coord%nam)//' constraint' &
!        ,__LINE__,__FILE__)
!     endif
!     !if(iret<=0 .or. iret>natm)then
!     if(abs(iret)>natm .or. iret==0) then
!       if(printable) write(nfout,'(a,i5)') 'invalid specification for '&
!      &  //trim(tag_atom)//trim(tmpstr)//' : ' ,iret
!       call phase_error_with_msg(nfout,'invalid specification for '//trim(tag_atom)//trim(tmpstr)&
!       ,__LINE__,__FILE__)
!     endif
!     constrainable_coord%associated_atoms(i) = iret
!   enddo
  end subroutine parse_associated_atoms

  subroutine parse_mobile_and_monitor(constrainable_coord)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    integer :: f_getIntValue
    integer :: iret
    if( f_getIntValue(tag_mobile,iret)==0 ) then
      if(iret==1) then
        constrainable_coord%mobile = .true.
      else
        constrainable_coord%mobile = .false.
      endif
    endif

    if( f_getIntValue(tag_monitor,iret)==0 ) then
      if(iret==1) then
        constrainable_coord%monitor = .true.
      else
        constrainable_coord%monitor = .false.
      endif
    endif
  end subroutine parse_mobile_and_monitor

  subroutine prep_reac_coords_1D(constrainable_coord,uni,min_val)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    integer, intent(in) :: uni
    real(DP), intent(in), optional :: min_val
    character(len=256) :: tag
    integer :: f_getIntValue, f_getRealValue,f_selectBlock,f_selectParentBlock
    integer :: i
    integer :: iret
    real(DP) :: dret
    real(DP) :: factor
    real(DP) :: fval,incre,val
    logical  :: has_minval
    real(DP) :: mval
    has_minval = .false.
    if(present(min_val)) then
      mval = min_val
      has_minval = .true.
    endif
    tag=''
    if( f_selectBlock(tag_reaction_coordinate)==0 ) then
      factor = 1
      if(uni.eq.LENGTH .or. uni.eq.DISTANCE_FROM_POS)then
        factor = 1.0d0
        tag = 'bohr'
      else if(uni.eq.ANGLE)then
        tag = 'degree'
        factor = PAI/180.0d0
      endif

      if(f_getIntValue(tag_sw_reaction_coordinate,iret)==0)then
        if(iret==0)then
          constrainable_coord%is_reaction_coordinate = .false.
        else
          constrainable_coord%is_reaction_coordinate = .true.
        endif
      endif
      if(f_getRealValue(tag_init_value,dret,trim(tag))==0)then
        constrainable_coord%value(1) = dret * factor
      endif
      if(constrainable_coord%is_reaction_coordinate)then
        if(f_getRealValue(tag_final_value,dret,trim(tag))==0)then
          fval = dret * factor
        endif
        incre = 0.1d0
        if(f_getRealValue(tag_increment,dret,trim(tag))==0)then
          incre = dabs(dret) * factor * sign(1.0d0,fval-constrainable_coord%value(1))
        endif
        constrainable_coord%n_reaction_coords = abs(nint((fval-constrainable_coord%value(1))/dabs(incre)))
        allocate(constrainable_coord%reaction_coords(constrainable_coord%n_reaction_coords &
      &        , constrainable_coord%ndim))
        allocate(constrainable_coord%finished(constrainable_coord%n_reaction_coords+1))
        constrainable_coord%finished=.false.
        do i=1,constrainable_coord%n_reaction_coords
          val = constrainable_coord%value(1) + incre*i
          if(has_minval) then
            if (val<min_val) then
              constrainable_coord%n_reaction_coords = constrainable_coord%n_reaction_coords-1
              cycle
            endif
          endif
         constrainable_coord%reaction_coords(i,1) = constrainable_coord%value(1) &
      & + incre * i
        enddo
      endif
      iret = f_selectParentBlock()
    endif
  end subroutine prep_reac_coords_1D

  subroutine update_det_metric(constrainable_coord)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    integer :: i,j
    constrainable_coord%det_metric=0.0d0
    do i=1,constrainable_coord%n_associated_atoms
      constrainable_coord%det_metric = constrainable_coord%det_metric + &
     & dot_product(constrainable_coord%dsigma(1,i,:),constrainable_coord%dsigma(1,i,:)) &
     & /(amion(ityp(constrainable_coord%associated_atoms(i))))
    enddo
  end subroutine update_det_metric

  subroutine init_constrainable_coord(constrainable_coord)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    constrainable_coord%mobile                 = .false.
    constrainable_coord%monitor                = .false.
    constrainable_coord%nam                    = ''
    constrainable_coord%descri                 = ''
    constrainable_coord%unit_name              = ''
    constrainable_coord%ndim                   = 0
    constrainable_coord%n_associated_atoms     = 0
    constrainable_coord%nfree_decrement        = 0
    constrainable_coord%det_metric             = 0.d0
    constrainable_coord%rmix_coords            = 1.d0
    constrainable_coord%rmix_velocity          = 1.d0
    constrainable_coord%n_reaction_coords      = 0
    constrainable_coord%is_reaction_coordinate = .false.
    constrainable_coord%reac_id              = -1
    constrainable_coord%reac_coord_count     = 1
    constrainable_coord%reac_direction       = 1
    constrainable_coord%normalization_factor = 1.d0
    constrainable_coord%ngroup               = 0
  end subroutine init_constrainable_coord

  logical function is_vector_non_zero(vect,n)
    integer, intent(in) :: n
    real(DP), intent(in), dimension(n) :: vect
    real(DP) :: eps
    integer :: i
    is_vector_non_zero=.false.
    eps=1.d-12
    do i=1,n
      if(dabs(vect(i)).gt.eps) then
         is_vector_non_zero = .true.
      endif
    enddo
  end function is_vector_non_zero

  subroutine get_curr_com(constrainable_coord,com,sum_mass,ignore)
    type(constrainable_coords_t), intent(in) :: constrainable_coord
    real(DP), dimension(3),intent(out) :: com
    real(DP), intent(out) :: sum_mass
    integer,intent(in) :: ignore
    real(DP), dimension(3) :: res
    integer :: i,i1
    integer :: ig=-1
    ig = ignore
    res=0.d0
    sum_mass = 0.d0
    do i=1,constrainable_coord%n_associated_atoms
      i1 = constrainable_coord%associated_atoms(i)
      if(i1==ig) cycle
      res(1) = res(1) + amion(ityp(i1))*x_nopbc(i1)
      res(2) = res(2) + amion(ityp(i1))*y_nopbc(i1)
      res(3) = res(3) + amion(ityp(i1))*z_nopbc(i1)
      sum_mass = sum_mass + amion(ityp(i1))
    enddo
    res(:) = res(:)/sum_mass
    com = res
  end subroutine get_curr_com

  subroutine get_curr_com_group(ig,constrainable_coord,com,sum_mass)
    integer, intent(in) :: ig
    type(constrainable_coords_t), intent(in) :: constrainable_coord
    real(DP), dimension(3),intent(out) :: com
    real(DP), intent(out) :: sum_mass
    real(DP), dimension(3) :: res
    integer :: i,i1
    res=0.d0
    sum_mass = 0.d0
    do i=1,constrainable_coord%natm_per_group(ig)
      i1 = constrainable_coord%atom_groups(i,ig)
      res(1) = res(1) + amion(ityp(i1))*x_nopbc(i1)
      res(2) = res(2) + amion(ityp(i1))*y_nopbc(i1)
      res(3) = res(3) + amion(ityp(i1))*z_nopbc(i1)
      sum_mass = sum_mass + amion(ityp(i1))
    enddo
    res(:) = res(:)/sum_mass
    com = res
  end subroutine get_curr_com_group

  real(DP) function get_distance_between_com(i,j,constrainable_coord)
    integer, intent(in) :: i,j
    type(constrainable_coords_t), intent(in) :: constrainable_coord
    real(DP), dimension(3) :: comi, comj, comij
    real(DP) :: sum_mass
    call get_curr_com_group(i,constrainable_coord,comi,sum_mass)
    call get_curr_com_group(j,constrainable_coord,comj,sum_mass)
    comij = comi-comj
    get_distance_between_com = sqrt(dot_product(comij,comij))
  end function get_distance_between_com

  integer function get_atom_id_constrainable_coord(i,constrainable_coord)
    integer, intent(in) :: i
    type(constrainable_coords_t), intent(in) :: constrainable_coord
    integer ::j
    get_atom_id_constrainable_coord = 0
    do j=1,constrainable_coord%n_associated_atoms
      if(constrainable_coord%associated_atoms(j)==i) then
        get_atom_id_constrainable_coord = j
        return
      endif
    enddo
    return
  end function get_atom_id_constrainable_coord

end module m_constraints_data

module distance_from_com_constraint
use m_constraints_data

implicit none

real(kind=DP), private :: sum_mass
integer, private :: target_atom

contains
  subroutine update_sigma_dfc(constrainable_coord)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    constrainable_coord%sigma(1) = &
    & get_curr_dfc(constrainable_coord)-constrainable_coord%value(1)
  end subroutine update_sigma_dfc

  subroutine update_dsigma_dfc(constrainable_coord)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    real(kind=DP) :: dfc,idfc
    real(kind=DP),dimension(3) :: com,diffcom
    integer :: itarg,i
    constrainable_coord%dsigma = 0.d0
    dfc = get_curr_dfc(constrainable_coord)
    if(dfc<1e-15) then
      return
    endif
    idfc = 1.0d0/dfc
    itarg = constrainable_coord%associated_atoms(target_atom)
    call get_curr_com(constrainable_coord,com,sum_mass,itarg)
    diffcom(1) = x_nopbc(itarg)-com(1)
    diffcom(2) = y_nopbc(itarg)-com(2)
    diffcom(3) = z_nopbc(itarg)-com(3)
    do i=1,constrainable_coord%n_associated_atoms
      if(i==target_atom)then
        constrainable_coord%dsigma(1,i,:) = idfc*diffcom(:)
      else
        constrainable_coord%dsigma(1,i,:) = idfc*diffcom(:) &
     &  * amion(ityp(constrainable_coord%associated_atoms(i)))/sum_mass
      endif
    enddo
  end subroutine update_dsigma_dfc

  subroutine monitor_dfc(constrainable_coord)
    type(constrainable_coords_t), intent(in) :: constrainable_coord
    real(kind=DP) :: dfc
    dfc = get_curr_dfc(constrainable_coord)
    write(nfout,'(a,i8,a,f10.5)') 'distance between atom ',constrainable_coord%associated_atoms(target_atom), &
    & ' and the center of mass for the specified atoms : ',dfc
  end subroutine monitor_dfc

  subroutine read_dfc_constraint(constrainable_coord)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    integer :: iret
    integer :: f_selectBlock, f_getIntValue, f_selectParentBlock
    real(kind=DP), dimension(3) :: com
    constrainable_coord%nam                = 'distance_from_com'
    constrainable_coord%descri             = &
  & 'constrain the distance from a target atom and the center of mass for a set of atoms'
    constrainable_coord%unit_name          = 'bohr'
    constrainable_coord%typ                = DISTANCE_FROM_COM
    constrainable_coord%nfree_decrement    = 4
    constrainable_coord%ndim               = 1
    target_atom = 1
    call parse_associated_atoms(constrainable_coord)
    call alloc_constrainable_coord(constrainable_coord)
    call parse_mobile_and_monitor(constrainable_coord)
    if(f_selectBlock(tag_distance_from_com)==0)then
      if(f_getIntValue('target_atom',iret)==0) target_atom = iret
      iret = f_selectParentBlock()
    endif
    constrainable_coord%value(1) = get_curr_dfc(constrainable_coord)
    call prep_reac_coords_1D(constrainable_coord,LENGTH,very_small)
  end subroutine read_dfc_constraint

  real(kind=DP) function get_curr_dfc(constrainable_coord)
    type(constrainable_coords_t), intent(in) :: constrainable_coord
    real(kind=DP),dimension(3) :: com
    integer :: itarg
    itarg = constrainable_coord%associated_atoms(target_atom)
    call get_curr_com(constrainable_coord,com,sum_mass,itarg)
    get_curr_dfc = &
    & dsqrt((x_nopbc(itarg)-com(1))**2+ &
    &       (y_nopbc(itarg)-com(2))**2+ &
    &       (z_nopbc(itarg)-com(3))**2)
  end function get_curr_dfc

end module distance_from_com_constraint

module distance_from_pos_constraint
use m_constraints_data
implicit none
real(DP), private, dimension(3) :: from_pos
contains
  subroutine update_sigma_dfp(constrainable_coord)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    constrainable_coord%sigma(1) = &
   & get_curr_dfp(constrainable_coord)-constrainable_coord%value(1)
  end subroutine update_sigma_dfp

  real(DP) function get_curr_dfp(constrainable_coord)
    type(constrainable_coords_t), intent(in) :: constrainable_coord
    integer :: i1
    i1 = constrainable_coord%associated_atoms(1)
    get_curr_dfp = dsqrt((x_nopbc(i1)-from_pos(1))**2  &
  &                            +       (y_nopbc(i1)-from_pos(2))**2  &
  &                            +       (z_nopbc(i1)-from_pos(3))**2)
  end function get_curr_dfp

  subroutine update_dsigma_dfp(constrainable_coord)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    integer :: i1
    real(DP) :: tmpx,tmpy,tmpz,tmpdist
    real(DP),dimension(3) :: tmpvec
    tmpdist = get_curr_dfp(constrainable_coord)
    tmpvec = 0.d0
    i1 = constrainable_coord%associated_atoms(1)
    tmpvec(1) = x_nopbc(i1)-from_pos(1)
    tmpvec(2) = y_nopbc(i1)-from_pos(2)
    tmpvec(3) = z_nopbc(i1)-from_pos(3)
    tmpx = tmpvec(1)/tmpdist
    tmpy = tmpvec(2)/tmpdist
    tmpz = tmpvec(3)/tmpdist
    constrainable_coord%dsigma(1,1,1) = +tmpx
    constrainable_coord%dsigma(1,1,2) = +tmpy
    constrainable_coord%dsigma(1,1,3) = +tmpz
  end subroutine update_dsigma_dfp

  subroutine monitor_dfp(constrainable_coord)
    type(constrainable_coords_t), intent(in) :: constrainable_coord
    integer :: i1
    real(DP) blen
    i1=constrainable_coord%associated_atoms(1)
    blen = get_curr_dfp(constrainable_coord)
    if(printable)write(nfout,'(a,3f10.5,a,i8,a,f10.5,a)') &
  & '    distance from position (',from_pos(1),from_pos(2),from_pos(3),') to atom ',i1,':',blen,' bohr'
  end subroutine monitor_dfp

  subroutine read_dfp_constraint(constrainable_coord)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    integer :: iret,iiret,i
    integer :: f_getIntValue,f_selectBlock,f_selectParentBlock,f_getRealValue
    real(DP) :: dret
    real(DP) :: dr,fval
    integer :: number_of_reac_coords
    character(256) tmpstr

    constrainable_coord%nam                = 'distance_from_pos'
    constrainable_coord%descri             = 'constrain the distance from some point in space to the specified atom'
    constrainable_coord%unit_name          = 'bohr'
    constrainable_coord%typ                = DISTANCE_FROM_POS
    constrainable_coord%nfree_decrement    = 1
    constrainable_coord%ndim               = 1
    constrainable_coord%n_associated_atoms = 1

    call parse_associated_atoms(constrainable_coord)
    call alloc_constrainable_coord(constrainable_coord)
    call parse_mobile_and_monitor(constrainable_coord)
    from_pos=0
    if(f_selectBlock(tag_distance_from_pos)==0)then
      if (f_getRealValue('posx',dret,'bohr')==0) from_pos(1) = dret
      if (f_getRealValue('posy',dret,'bohr')==0) from_pos(2) = dret
      if (f_getRealValue('posz',dret,'bohr')==0) from_pos(3) = dret
      iret = f_selectParentBlock()
    else
      write(0,'(a)') 'the coordinate of some point in space is needed in order to use this constraint.'
      call phase_error_with_msg(nfout,'the coordinate of some point in space is needed in order to use this constraint.'&
                               ,__LINE__,__FILE__)
    endif
    constrainable_coord%value(1) = get_curr_dfp(constrainable_coord)

    call prep_reac_coords_1D(constrainable_coord,DISTANCE_FROM_POS,very_small)

  end subroutine read_dfp_constraint
end module distance_from_pos_constraint

module bond_length_constraint
use m_constraints_data
implicit none
integer,private,dimension(3) :: exclude
contains
  subroutine update_sigma_bl(constrainable_coord)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    constrainable_coord%sigma(1) = &
    get_curr_bl(constrainable_coord)-constrainable_coord%value(1)
  end subroutine update_sigma_bl

  real(DP) function get_curr_bl(constrainable_coord)
    type(constrainable_coords_t), intent(in) :: constrainable_coord
    integer :: i1,i2
    get_curr_bl = get_distance_between_com(1,2,constrainable_coord)
  end function get_curr_bl

  subroutine update_dsigma_bl(constrainable_coord)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    integer :: i,j
    integer :: k,l
    real(DP) :: sumi,sumj,isumi,isumj
    real(DP), dimension(3) :: comi, comj, comij
    real(DP) :: sigma,sigmai
    i = 1
    j = 2
    sigma = get_distance_between_com(i,j,constrainable_coord)
    if(sigma<very_small) then
        if(printable) write(nfout,'(a,2i8,a,f20.15)') 'distance between atom groups ',i,j,' is too short',sigma
        call phase_error_with_msg(nfout,&
        'distance between atom groups is too short '&
        ,__LINE__,__FILE__)
    endif
    sigmai = 1.d0/sigma
    call get_curr_com_group(i,constrainable_coord,comi,sumi)
    call get_curr_com_group(j,constrainable_coord,comj,sumj)
    comij = comi-comj
    isumi = 1.d0/sumi
    isumj = 1.d0/sumj
    constrainable_coord%dsigma = 0.d0
    do k=1,constrainable_coord%natm_per_group(i)
      l = get_atom_id_constrainable_coord( &
          constrainable_coord%atom_groups(k,i),constrainable_coord)
      if(l==0) then
        if(printable) write(nfout,'(a,i8)') 'invalid atom index in atom group ',i
        call phase_error_with_msg(nfout,&
        'invalid atom index in atom group'&
        ,__LINE__,__FILE__)
      endif
      constrainable_coord%dsigma(1,l,:) = &
      sigmai*comij(:)*amion(ityp(constrainable_coord%atom_groups(k,i)))*isumi
    enddo
    do k=1,constrainable_coord%natm_per_group(j)
      l = get_atom_id_constrainable_coord( &
          constrainable_coord%atom_groups(k,j),constrainable_coord)
      if(l==0) then
        if(printable) write(nfout,'(a,i8)') 'invalid atom index in atom group ',j
        call phase_error_with_msg(nfout,&
        'invalid atom index in atom group'&
        ,__LINE__,__FILE__)
      endif
      constrainable_coord%dsigma(1,l,:) = constrainable_coord%dsigma(1,l,:) &
     -sigmai*comij(:)*amion(ityp(constrainable_coord%atom_groups(k,j)))*isumj
    enddo
  end subroutine update_dsigma_bl

  subroutine monitor_bl(constrainable_coord)
    type(constrainable_coords_t), intent(in) :: constrainable_coord
    integer :: i1,i2
    real(DP) :: sigma
    sigma = get_distance_between_com(1,2,constrainable_coord)
    if(constrainable_coord%natm_per_group(1) == 1 .and. &
       constrainable_coord%natm_per_group(2) == 1) then
       i1=constrainable_coord%associated_atoms(1)
       i2=constrainable_coord%associated_atoms(2)
       if(printable)write(nfout,'(a,i5,a,i5,a,f10.5,a)') &
    & '    bond-length between atoms ',i1,' and ',i2,' : ',sigma,' bohr'
    else
      if(printable) write(nfout,'(a,f10.5,a)') &
      'distance between atom group 1 and 2 : ',sigma,' bohr'
    endif
  end subroutine monitor_bl

  subroutine read_bl_constraint(constrainable_coord)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    integer :: iret,iiret,i
    integer :: f_getIntValue,f_selectBlock,f_selectParentBlock,f_getRealValue
    real(DP) :: dret
    real(DP) :: dr,fval
    integer :: number_of_reac_coords,ntot
    character(10) strx,stry,strz
    character(256) tmpstr

    constrainable_coord%nam                = 'bond_length'
    constrainable_coord%descri             = 'constrain the length between two atoms/group of atoms'
    constrainable_coord%unit_name          = 'bohr'
    constrainable_coord%typ                = BOND_LENGTH
    constrainable_coord%nfree_decrement    = 1
    constrainable_coord%ndim               = 1

    call parse_atom_groups(constrainable_coord)
    if(constrainable_coord%ngroup<2) then
      call phase_error_with_msg(nfout, &
      'you must specify at least two atoms in order to define a '//trim(constrainable_coord%nam)//' constraint' &
      ,__LINE__,__FILE__)
    endif
    call alloc_constrainable_coord(constrainable_coord)
    call parse_mobile_and_monitor(constrainable_coord)
    exclude = 0
    if(f_selectBlock(tag_bond_length)==0) then
      if(f_getIntValue("exclude_x",iret)==0) exclude(1) = iret
      if(f_getIntValue("exclude_y",iret)==0) exclude(2) = iret
      if(f_getIntValue("exclude_z",iret)==0) exclude(3) = iret
      if(exclude(1)==1.or.exclude(2)==1.or.exclude(3)==1)then
        if(printable)then
          strx=''
          stry=''
          strz=''
          if(exclude(1)==1) strx = 'x'
          if(exclude(2)==1) stry = 'y'
          if(exclude(3)==1) strz = 'z'
          write(nfout,'(a)') 'components excluded for bond-length calculation : '&
          & //trim(strx)//' '//trim(stry)//' '//trim(strz)
        endif
      endif
      iret = f_selectParentBlock()
    endif
    constrainable_coord%value(1) = get_distance_between_com(1,2,constrainable_coord)

    call prep_reac_coords_1D(constrainable_coord,LENGTH,very_small)

  end subroutine read_bl_constraint

end module bond_length_constraint

module bond_angle_constraint
use m_constraints_data
implicit none
contains

  subroutine read_bang_constraint(constrainable_coord)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    integer :: i,ntot
    constrainable_coord%nam                = 'bond_angle'
    constrainable_coord%descri             = 'constrain the bond-angle among three atoms/group of atoms'
    constrainable_coord%unit_name          = 'radian'
    constrainable_coord%typ                = BOND_ANGLE
    constrainable_coord%nfree_decrement    = 1
    constrainable_coord%ndim               = 1
    call parse_atom_groups(constrainable_coord)
    if(constrainable_coord%ngroup<3) then
      call phase_error_with_msg(nfout, &
      'you must specify at least three atoms in order to define a '//trim(constrainable_coord%nam)//' constraint' &
      ,__LINE__,__FILE__)
    endif
    call alloc_constrainable_coord(constrainable_coord)
    call parse_mobile_and_monitor(constrainable_coord)
    if(constrainable_coord%mobile) constrainable_coord%nfree_decrement=0
    constrainable_coord%value(1) = get_curr_bang(constrainable_coord)
    call prep_reac_coords_1D(constrainable_coord,ANGLE)
  end subroutine read_bang_constraint

  subroutine monitor_bang(constrainable_coord)
    type(constrainable_coords_t), intent(in) :: constrainable_coord
    integer :: i1, i2, i3
    real(DP) bang
    bang = get_curr_bang(constrainable_coord)
    if(constrainable_coord%natm_per_group(1)==1 .and. &
       constrainable_coord%natm_per_group(1)==1 .and. &
       constrainable_coord%natm_per_group(1)==1) then
       i1 = constrainable_coord%associated_atoms(1)
       i2 = constrainable_coord%associated_atoms(2)
       i3 = constrainable_coord%associated_atoms(3)
       if(printable)write(nfout,'(a,i5,a,i5,a,i5,a,f10.5,a)') &
    & '    bond-angle among atoms ',i1,', ',i2,' and ',i3,' : ',180.0d0*bang/PAI,' degrees'
    else
      if(printable)write(nfout,'(a,f10.5,a)') &
    & '    bond-angle among atom group 1 2 3 : ',180.0d0*bang/PAI,' degrees'
    endif
  end subroutine monitor_bang

  real(DP) function get_curr_bang(constrainable_coord)
    type(constrainable_coords_t), intent(in) :: constrainable_coord
    integer :: i1,i2,i3
    real(DP) :: rji,rki,cosine,sum_mass
    real(DP), dimension(3) :: vecji,vecki,com1,com2,com3
    real(DP) :: eps=1.d-12
    rji = get_distance_between_com(1,2,constrainable_coord)
    rki = get_distance_between_com(3,2,constrainable_coord)
    call get_curr_com_group(1,constrainable_coord,com1,sum_mass)
    call get_curr_com_group(2,constrainable_coord,com2,sum_mass)
    call get_curr_com_group(3,constrainable_coord,com3,sum_mass)
    vecji = com1-com2
    vecki = com3-com2
    cosine = dot_product(vecji,vecki)/(rji*rki)
    if(cosine.gt.1)  cosine =  1.d0
    if(cosine.lt.-1) cosine = -1.d0
    get_curr_bang = dacos(cosine)
  end function get_curr_bang

  subroutine update_sigma_bang(constrainable_coord)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    real(DP) :: currbang
    currbang = get_curr_bang(constrainable_coord)
    constrainable_coord%sigma(1) = currbang-constrainable_coord%value(1)
  end subroutine update_sigma_bang

  subroutine update_dsigma_bang(constrainable_coord)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    real(DP) :: cosine,sine,rji,rki,theta,isum1,isum2,isum3,sm
    real(DP),dimension(3) :: vecji,vecki,com1,com2,com3
    real(DP), dimension(3) :: der1,der3
    integer :: i,j,k
    constrainable_coord%dsigma = 0.0d0

    rji = get_distance_between_com(1,2,constrainable_coord)
    rki = get_distance_between_com(3,2,constrainable_coord)
    call get_curr_com_group(1,constrainable_coord,com1,sm)
    isum1 = 1.d0/sm
    call get_curr_com_group(2,constrainable_coord,com2,sm)
    isum2 = 1.d0/sm
    call get_curr_com_group(3,constrainable_coord,com3,sm)
    isum3 = 1.d0/sm
    vecji = com1-com2
    vecki = com3-com2

    theta = get_curr_bang(constrainable_coord)
    cosine = dcos(theta)
    sine   = dsin(theta)
    der1(:) = ((vecji(:)/rji)*cosine - (vecki(:)/rki))/(rji*sine)
    der3(:) = ((vecki(:)/rki)*cosine - (vecji(:)/rji))/(rki*sine)

    do j=1,constrainable_coord%natm_per_group(1)
      k = get_atom_id_constrainable_coord( &
          constrainable_coord%atom_groups(j,1),constrainable_coord)
      if(k==0) then
        if(printable) write(nfout,'(a,i8)') 'invalid atom index in atom group ',i
        call phase_error_with_msg(nfout,&
        'invalid atom index in atom group'&
        ,__LINE__,__FILE__)
      endif
      constrainable_coord%dsigma(1,k,:) = der1(:)*isum1*amion(ityp(constrainable_coord%atom_groups(j,1)))
    enddo
    do j=1,constrainable_coord%natm_per_group(3)
      k = get_atom_id_constrainable_coord( &
          constrainable_coord%atom_groups(j,3),constrainable_coord)
      if(k==0) then
        if(printable) write(nfout,'(a,i8)') 'invalid atom index in atom group ',i
        call phase_error_with_msg(nfout,&
        'invalid atom index in atom group'&
        ,__LINE__,__FILE__)
      endif
      constrainable_coord%dsigma(1,k,:) = constrainable_coord%dsigma(1,k,:) &
                                        + der3(:)*isum3*amion(ityp(constrainable_coord%atom_groups(j,3)))
    enddo
    do j=1,constrainable_coord%natm_per_group(2)
      k = get_atom_id_constrainable_coord( &
          constrainable_coord%atom_groups(j,2),constrainable_coord)
      if(k==0) then
        if(printable) write(nfout,'(a,i8)') 'invalid atom index in atom group ',i
        call phase_error_with_msg(nfout,&
        'invalid atom index in atom group'&
        ,__LINE__,__FILE__)
      endif
      constrainable_coord%dsigma(1,k,:) = constrainable_coord%dsigma(1,k,:) &
                                        - (der1(:)+der3(:))*isum2*amion(ityp(constrainable_coord%atom_groups(j,2)))
    enddo
  end subroutine update_dsigma_bang
end module bond_angle_constraint

module dihedral_angle_constraint
use m_constraints_data
implicit none
character(len("abs")), private, parameter :: tag_abs = "abs"
logical, private :: babs
contains
  subroutine read_dang_constraint(constrainable_coord)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    integer :: f_selectBlock, f_getBoolValue, f_selectParentBlock, iret
    logical :: bret
    integer :: i,ntot

    babs                                   = .false.
    constrainable_coord%nam                = 'dihedral_angle'
    constrainable_coord%descri             = 'constrain the dihedral angle among four atoms/group of atomss'
    constrainable_coord%unit_name          = 'radian'
    constrainable_coord%typ                = DIHEDRAL_ANGLE
    constrainable_coord%nfree_decrement    = 1
    constrainable_coord%ndim               = 1
    call parse_atom_groups(constrainable_coord)
    if(constrainable_coord%ngroup<4) then
      call phase_error_with_msg(nfout, &
      'you must specify at least four atoms in order to define a '//trim(constrainable_coord%nam)//' constraint' &
      ,__LINE__,__FILE__)
    endif
    call alloc_constrainable_coord(constrainable_coord)
    call parse_mobile_and_monitor(constrainable_coord)
    if(f_selectBlock('dihedral')==0)then
      if(f_getBoolValue(tag_abs,bret)==0) babs=bret
      if(babs)then
        constrainable_coord%descri         = &
        & 'constrain the absolute value of the dihedral angle among four atoms'
      endif
      iret = f_selectParentBlock()
    endif

    constrainable_coord%value(1) = get_curr_dang(constrainable_coord)
    call prep_reac_coords_1D(constrainable_coord,ANGLE)
  end subroutine read_dang_constraint

  real(DP) function get_curr_dang(constrainable_coord,sig)
    type(constrainable_coords_t), intent(in) :: constrainable_coord
    real(DP), intent(out), optional :: sig
    integer :: ii,ij,ik,il
    real(DP), dimension(3) :: mvec,nvec,mnvec
    real(DP), dimension(3) :: com1,com2,com3,com4
    real(DP), dimension(3) :: rij,rkj,rkl
    real(DP) :: m,n
    real(DP) :: cosine
    real(DP) :: eps
    real(DP) :: hugo
    real(DP) :: sum_mass
    eps=1.d-12

    call get_curr_com_group(1,constrainable_coord,com1,sum_mass)
    call get_curr_com_group(2,constrainable_coord,com2,sum_mass)
    call get_curr_com_group(3,constrainable_coord,com3,sum_mass)
    call get_curr_com_group(4,constrainable_coord,com4,sum_mass)
    rij = com1-com2
    rkj = com3-com2
    rkl = com3-com4
    call cross_product(rij,rkj,mvec)
    call cross_product(rkj,rkl,nvec)
    m = get_norm(mvec)
    n = get_norm(nvec)
    if (m.lt.eps.or.n.lt.eps) then
      get_curr_dang = 0.d0
      return
    endif
    cosine = dot_product(mvec,nvec)/(m*n)
    if(cosine.ge.1)  cosine =  1.d0
    if(cosine.le.-1) cosine = -1.d0
    call cross_product(mvec,nvec,mnvec)
    hugo = sign(1.0d0,dot_product(rkj,mnvec))
    if(.not.babs) then
      get_curr_dang = dabs(dacos(cosine))*hugo
    else
      get_curr_dang = dabs(dacos(cosine))
    endif
    !get_curr_dang = dacos(cosine)
    if(printable.and.ipri>=2) write(nfout,'(9f15.10)') mvec(1),mvec(2),mvec(3),nvec(1),nvec(2),nvec(3),m,n,cosine
    if(present(sig))sig=hugo
  end function get_curr_dang

  subroutine monitor_dang(constrainable_coord)
    type(constrainable_coords_t), intent(in) :: constrainable_coord
    integer :: i1,i2,i3,i4
    real(DP) :: dihed,hugo
    dihed = get_curr_dang(constrainable_coord,hugo)
    if(constrainable_coord%natm_per_group(1)==1 .and. &
       constrainable_coord%natm_per_group(2)==1 .and. &
       constrainable_coord%natm_per_group(3)==1 .and. &
       constrainable_coord%natm_per_group(4)==1) then
       i1 = constrainable_coord%associated_atoms(1)
       i2 = constrainable_coord%associated_atoms(2)
       i3 = constrainable_coord%associated_atoms(3)
       i4 = constrainable_coord%associated_atoms(4)
       if(printable)write(nfout,'(a,i5,a,i5,a,i5,a,i5,a,f10.5,a)') &
   &   '    dihedral-angle among atoms ',i1,', ',i2,', ',i3,' and ',i4,' : ',180.d0*dihed/PAI,' degrees'
    else
      if(printable)write(nfout,'(a,f10.5,a)') &
     & '    dihedral-angle among atom groups 1 2 3 4 ',180.d0*dihed/PAI,' degrees'
    endif
    if(printable.and.ipri>=2) write(nfout,'(a,f10.1)') '      its sign:',hugo
  end subroutine monitor_dang

  subroutine update_sigma_dang(constrainable_coord)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    real(DP) :: currdihed
    currdihed = get_curr_dang(constrainable_coord)
    constrainable_coord%sigma(1) = currdihed-constrainable_coord%value(1)
  end subroutine update_sigma_dang

  subroutine update_dsigma_dang(constrainable_coord)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    integer :: ii,ij,ik,il
    real(DP) :: dang
    real(DP), dimension(3) :: mvec,nvec,mnvec
    real(DP), dimension(3) :: rij,rkj,rkl
    real(DP) :: m,n,rijl,rkjl,rkll
    real(DP),dimension(3) :: di,dl
    real(DP) :: cosine
    real(DP),dimension(3) :: tmpvec,tmpvec1,tmpvec2
    real(DP) :: eps
    real(DP) :: hugo
    real(DP) :: sum_mass,imass1,imass2,imass3,imass4
    real(DP), dimension(3) :: com1,com2,com3,com4
    integer  :: i,j,iatm
    eps=1.d-14
    !dang = get_curr_dang(constrainable_coord)
    constrainable_coord%dsigma=0.0d0

    call get_curr_com_group(1,constrainable_coord,com1,sum_mass)
    imass1 = 1.d0/sum_mass
    call get_curr_com_group(2,constrainable_coord,com2,sum_mass)
    imass2 = 1.d0/sum_mass
    call get_curr_com_group(3,constrainable_coord,com3,sum_mass)
    imass3 = 1.d0/sum_mass
    call get_curr_com_group(4,constrainable_coord,com4,sum_mass)
    imass4 = 1.d0/sum_mass
    rij = com1-com2
    rkj = com3-com2
    rkl = com3-com4

    call cross_product(rij,rkj,mvec)
    call cross_product(rkj,rkl,nvec)
    m = get_norm(mvec)
    n = get_norm(nvec)
    if (m.lt.eps.or.n.lt.eps) then
      return
    endif
    call cross_product(mvec,nvec,mnvec)
    hugo = 1.d0
    if(babs) hugo = sign(1.0d0,dot_product(rkj,mnvec))
    rijl = get_norm(rij)
    rkjl = get_norm(rkj)
    rkll = get_norm(rkl)
    di(:) = +hugo*rkjl*mvec(:)/(m*m)
    dl(:) = -hugo*rkjl*nvec(:)/(n*n)

    do i=1,constrainable_coord%natm_per_group(1)
      iatm = constrainable_coord%atom_groups(i,1)
      j = get_atom_id_constrainable_coord(iatm,constrainable_coord)
      constrainable_coord%dsigma(1,j,:) =  di(:)*amion(ityp(iatm))*imass1
    enddo

    do i=1,constrainable_coord%natm_per_group(2)
      iatm = constrainable_coord%atom_groups(i,2)
      j = get_atom_id_constrainable_coord(iatm,constrainable_coord)
      constrainable_coord%dsigma(1,j,:) = constrainable_coord%dsigma(1,j,:) + &
     & (-di(:) &
     & + di(:)*dot_product(rij,rkj)/(rkjl*rkjl)  &
     & - dl(:)*dot_product(rkl,rkj)/(rkjl*rkjl)) * amion(ityp(iatm))*imass2
    enddo

    do i=1,constrainable_coord%natm_per_group(3)
      iatm = constrainable_coord%atom_groups(i,3)
      j = get_atom_id_constrainable_coord(iatm,constrainable_coord)
      constrainable_coord%dsigma(1,j,:) = constrainable_coord%dsigma(1,j,:) + &
     & (-dl(:) &
     & - di(:)*dot_product(rij,rkj)/(rkjl*rkjl)  &
     & + dl(:)*dot_product(rkl,rkj)/(rkjl*rkjl)) * amion(ityp(iatm))*imass3
    enddo

    do i=1,constrainable_coord%natm_per_group(4)
      iatm = constrainable_coord%atom_groups(i,4)
      j = get_atom_id_constrainable_coord(iatm,constrainable_coord)
      constrainable_coord%dsigma(1,j,:) =  dl(:)*amion(ityp(iatm))*imass4
    enddo

  end subroutine update_dsigma_dang

end module dihedral_angle_constraint

module plane_constraint
use m_constraints_data
implicit none
real(DP), dimension(3), private :: norm
real(DP), dimension(3), private :: origin
integer :: iret
character(len("normx")), private, parameter :: tag_normx = "normx"
character(len("normy")), private, parameter :: tag_normy = "normy"
character(len("normz")), private, parameter :: tag_normz = "normz"
contains
  subroutine read_plane_constraint(constrainable_coord)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    integer :: iret,iiret,i,j,k
    real(DP) :: dret, dop
    integer :: f_getIntValue,f_selectBlock,f_selectParentBlock,f_getRealValue
    constrainable_coord%nam                = 'plane'
    constrainable_coord%descri             = 'constrain the specified atoms within the specified plane.'
    constrainable_coord%unit_name          = ''
    constrainable_coord%typ                = PLANE
    constrainable_coord%ndim               = 1
    call parse_associated_atoms(constrainable_coord)
    constrainable_coord%nfree_decrement    = constrainable_coord%n_associated_atoms
    call alloc_constrainable_coord(constrainable_coord)
    call parse_mobile_and_monitor(constrainable_coord)

    norm=0.d0
    if(f_getRealValue(tag_normx,dret,'')==0) norm(1) = dret
    if(f_getRealValue(tag_normy,dret,'')==0) norm(2) = dret
    if(f_getRealValue(tag_normz,dret,'')==0) norm(3) = dret
    if(f_selectBlock(tag_plane)==0)then
      if(f_getRealValue(tag_normx,dret,'')==0) norm(1) = dret
      if(f_getRealValue(tag_normy,dret,'')==0) norm(2) = dret
      if(f_getRealValue(tag_normz,dret,'')==0) norm(3) = dret
      iret = f_selectParentBlock()
    endif

    if(.not.is_vector_non_zero(norm,3)) then
      write(0,*) 'norm is undefined. you must define a norm in order to constrain an atom within a plane.'
      call phase_error_with_msg(nfout,'norm is undefined. you must define a norm in order to constrain an atom within a plane.' &
      ,__LINE__,__FILE__)
    endif

    !call generate_origin(constrainable_coord)
    call generate_origin2(constrainable_coord)

    dop = dot_product(norm,origin)
    do i=1,constrainable_coord%n_associated_atoms
      constrainable_coord%value(i) = dop
      if(printable)write(nfout,'(a)') '  norm dot (currxyz-origin)'
      if(printable)write(nfout,'(e15.5)') get_xyz_dot_norm(constrainable_coord,i)-constrainable_coord%value(i)
    enddo
    call prep_reac_coord_plane(constrainable_coord)
  end subroutine read_plane_constraint

  subroutine prep_reac_coord_plane(constrainable_coord)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    integer :: f_getIntValue, f_getRealValue,f_selectBlock,f_selectParentBlock
    integer :: iret,i,j
    real(DP) :: dret,fval,ival,incre
    real(DP),dimension(3) :: tmporig
    if(f_selectBlock(tag_reaction_coordinate)==0)then
      if(f_getIntValue(tag_sw_reaction_coordinate,iret)==0)then
        if(iret==0)then
          constrainable_coord%is_reaction_coordinate = .false.
        else
          constrainable_coord%is_reaction_coordinate = .true.
        endif
      endif
      if(constrainable_coord%is_reaction_coordinate)then
        ival=0
        if(f_getRealValue(tag_init_value,dret,'bohr')==0)then
          ival = dret
        endif
        fval=1
        if(f_getRealValue(tag_final_value,dret,'bohr')==0)then
          fval = dret
        endif
        incre = 0.1d0
        if(f_getRealValue(tag_increment,dret,'bohr')==0)then
          incre = dret
        endif
        constrainable_coord%n_reaction_coords = abs(nint((fval-ival)/incre))
        allocate(constrainable_coord%reaction_coords(constrainable_coord%n_reaction_coords &
      &        , constrainable_coord%ndim))
        allocate(constrainable_coord%finished(constrainable_coord%n_reaction_coords+1))
        constrainable_coord%finished=.false.
        tmporig(:) = ival*norm(:)+origin(:)
        constrainable_coord%value(:) = dot_product(norm,tmporig)
        do i=1,constrainable_coord%n_reaction_coords
          tmporig(:) = origin(:) + norm(:)*(ival+i*incre)
          constrainable_coord%reaction_coords(i,:) = dot_product(norm,tmporig)
        enddo
      endif
      iret = f_selectParentBlock()
    endif
  end subroutine prep_reac_coord_plane

  subroutine generate_origin2(constrainable_coord)
    type(constrainable_coords_t), intent(in) :: constrainable_coord
    integer :: i
    real(DP) :: eps = 1.e-10
    real(DP) :: delt = 0.1d0
    real(DP) :: tmpx,tmpy,tmpz
    real(DP) :: normnorm
    tmpx = cps(constrainable_coord%associated_atoms(1),1)-delt
    tmpy = cps(constrainable_coord%associated_atoms(1),2)-delt
    tmpz = cps(constrainable_coord%associated_atoms(1),3)
    normnorm=0.d0
    do i=1,3
      normnorm=normnorm+norm(i)**2
    enddo
    norm(:) = norm(:)/dsqrt(normnorm)
    if(dabs(norm(3)).gt.eps)then
       origin(1) = tmpx
       origin(2) = tmpy
       origin(3) = tmpz + (delt*norm(1)+delt*norm(2))/norm(3)
    else if (dabs(norm(2)).gt.eps)then
       origin(1) = tmpx
       origin(2) = cps(constrainable_coord%associated_atoms(1),2)+delt*norm(1)/norm(2)
       origin(3) = 0.d0
    else
       origin(1) = cps(constrainable_coord%associated_atoms(1),1)
       origin(2) = tmpy
       origin(3) = 0.d0
    endif
  end subroutine generate_origin2

  subroutine generate_origin(constrainable_coord)
    type(constrainable_coords_t), intent(in) :: constrainable_coord
    integer :: i,j,k
    real(DP) :: eps,theta,phi
    real(DP) :: tmpx,tmpy,tmpz
    real(DP) :: normnorm
    real(DP), dimension(3,3) :: rotmat,axis
    real(DP),dimension(3) :: tmpvec
    eps=1.d-12
    tmpx = cps(constrainable_coord%associated_atoms(1),1)
    tmpy = cps(constrainable_coord%associated_atoms(1),2)
    tmpz = cps(constrainable_coord%associated_atoms(1),3)
    normnorm=0.d0
    do i=1,3
      normnorm=normnorm+norm(i)**2
    enddo
    norm(:) = norm(:)/dsqrt(normnorm)

    axis = 0.0d0
    do i=1,3
      axis(i,i) = 1.0d0
    enddo

    if(dabs(norm(3)-1).gt.eps)then
      ! we first rotate the axis about the z-axis, then rotate the axis about the y'-axis.
      ! the resultant z''-axis has the same direction as the norm, and accordingly,
      ! the perpendicular direction to the norm is either the x''-axis or the y''-axis.
      theta = dacos(norm(3))
      tmpvec=0.0d0;tmpvec(1)=norm(1);tmpvec(2)=norm(2)
      normnorm=0.d0
      do i=1,3
        normnorm=normnorm+tmpvec(i)**2
      enddo
      tmpvec(:) = tmpvec(:)/dsqrt(normnorm)
      if(tmpvec(1)-1.gt.eps)  tmpvec(1)=tmpvec(1)-eps
      if(tmpvec(1)+1.lt.-eps) tmpvec(1)=tmpvec(1)+eps
      phi = dacos(tmpvec(1))

      rotmat = 0.0d0
      rotmat(1,1) =  dcos(phi)
      rotmat(2,2) =  dcos(phi)
      rotmat(2,1) = -dsin(phi)
      rotmat(1,2) = +dsin(phi)
      rotmat(3,3) = 1.0d0

      do k=1,3
        tmpvec = 0.0d0
        do i=1,3
          do j=1,3
            tmpvec(i) = tmpvec(i)+rotmat(i,j)*axis(j,k)
          enddo
        enddo
        do i=1,3
          axis(i,k) = tmpvec(i)
        enddo
      enddo

      rotmat = 0.0d0
      rotmat(1,1) =  dcos(theta)
      rotmat(2,2) =  1.0d0
      rotmat(3,1) =  dsin(theta)
      rotmat(1,3) = -dsin(theta)
      rotmat(3,3) =  dcos(theta)

      do k=1,3
        tmpvec = 0.0d0
        do i=1,3
          do j=1,3
            tmpvec(i) = tmpvec(i)+rotmat(i,j)*axis(j,k)
          enddo
        enddo
        do i=1,3
          axis(i,k) = tmpvec(i)
        enddo
      enddo
    endif

    origin(1) = axis(1,1)+tmpx
    origin(2) = axis(1,2)+tmpy
    origin(3) = axis(1,3)+tmpz

    if(printable)write(nfout,'(a)') 'plane statistics'
    if(printable)write(nfout,'(a)') '  norm'
    if(printable)write(nfout,'(3f15.10)') norm(1),norm(2),norm(3)
    if(printable)write(nfout,'(a)') '  rotated axis'
    do i=1,3
      if(printable)write(nfout,'(3f15.10)') axis(i,1),axis(i,2),axis(i,3)
    enddo
    if(printable)write(nfout,'(a)') '  origin'
    if(printable)write(nfout,'(3f15.10)') origin(1),origin(2),origin(3)
    if(printable)write(nfout,'(a)') '  currxyz'
    if(printable)write(nfout,'(3f15.10)') tmpx,tmpy,tmpz

  end subroutine generate_origin

  real(DP) function get_rvec_dot_norm(constrainable_coord,iatomind)
    type(constrainable_coords_t), intent(in) :: constrainable_coord
    integer, intent(in) :: iatomind
    integer :: i1
    real(DP), dimension(3) :: xyzdiff
    real(DP) :: dop
    i1=constrainable_coord%associated_atoms(iatomind)
    xyzdiff(1) = x_nopbc(i1)-origin(1)
    xyzdiff(2) = y_nopbc(i1)-origin(2)
    xyzdiff(3) = z_nopbc(i1)-origin(3)
    dop = dot_product(norm,xyzdiff)
    get_rvec_dot_norm = dop
  end function get_rvec_dot_norm

  real(DP) function get_xyz_dot_norm(constrainable_coord,iatomind)
    type(constrainable_coords_t), intent(in) :: constrainable_coord
    integer, intent(in) :: iatomind
    real(DP), dimension(3) :: xyz
    integer :: i1
    i1 = constrainable_coord%associated_atoms(iatomind)
    xyz(1) = x_nopbc(i1)
    xyz(2) = y_nopbc(i1)
    xyz(3) = z_nopbc(i1)
    get_xyz_dot_norm = dot_product(norm,xyz)
  end function get_xyz_dot_norm

  subroutine update_sigma_plane(constrainable_coord,iidim)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    integer, intent(in) :: iidim
    integer  :: i1
    real(DP), dimension(3) :: xyz
    real(DP) :: dop
    i1 = constrainable_coord%associated_atoms(iidim)
    xyz(1) = x_nopbc(i1)
    xyz(2) = y_nopbc(i1)
    xyz(3) = z_nopbc(i1)
    dop = dot_product(norm,xyz)
    constrainable_coord%sigma(iidim) = dop-constrainable_coord%value(iidim)
  end subroutine update_sigma_plane

  subroutine update_dsigma_plane(constrainable_coord,iidim)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    integer, intent(in) :: iidim
    integer :: i
    constrainable_coord%dsigma(iidim,iidim,:) = norm(:)
  end subroutine update_dsigma_plane

  subroutine monitor_plane(constrainable_coord)
    type(constrainable_coords_t), intent(in) :: constrainable_coord
    integer :: i
    do i=1,constrainable_coord%n_associated_atoms
      if(printable)write(nfout,'(a,i7,a,f10.5)')                  &
     & '    the dot product of the norm and atom ', constrainable_coord%associated_atoms(i),':',&
     & constrainable_coord%value(i)
    enddo
  end subroutine monitor_plane

end module plane_constraint

module center_of_mass_constraint
use m_constraints_data

implicit none
real(DP), private :: sum_mass
real(DP), private, dimension(3) :: direction
character(len("directionx")), private, parameter :: tag_directionx = "directionx"
character(len("directiony")), private, parameter :: tag_directiony = "directiony"
character(len("directionz")), private, parameter :: tag_directionz = "directionz"

contains
  subroutine read_com_constraint(constrainable_coord)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    integer :: i,j,i1
    real(DP),dimension(3) :: com
    constrainable_coord%nam                = 'center_of_mass'
    constrainable_coord%descri             = 'constrain the center of mass of the specified atoms'
    constrainable_coord%unit_name          = ''
    constrainable_coord%typ                = CENTER_OF_MASS
    constrainable_coord%ndim               = 3
    constrainable_coord%nfree_decrement    = 3
    call parse_associated_atoms(constrainable_coord)
    call alloc_constrainable_coord(constrainable_coord)
    call parse_mobile_and_monitor(constrainable_coord)
    call get_curr_com(constrainable_coord,com,sum_mass,-1)
    constrainable_coord%value(:) = com(:)

    call prep_reac_coord_com(constrainable_coord)
  end subroutine read_com_constraint

  subroutine prep_reac_coord_com(constrainable_coord)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    integer :: f_getIntValue, f_getRealValue,f_selectBlock,f_selectParentBlock
    integer :: iret,i,j
    real(DP) :: dret,fval,ival,incre
    real(DP),dimension(3) :: tmporig
    direction=0;direction(1)=1
    if(f_selectBlock(tag_center_of_mass)==0)then
      if(f_getRealValue(tag_directionx,dret,'')==0) direction(1) = dret
      if(f_getRealValue(tag_directiony,dret,'')==0) direction(2) = dret
      if(f_getRealValue(tag_directionz,dret,'')==0) direction(3) = dret
      iret = f_selectParentBlock()
    endif
    direction = direction/get_norm(direction)
    if(f_selectBlock(tag_reaction_coordinate)==0)then
      if(f_getIntValue(tag_sw_reaction_coordinate,iret)==0)then
        if(iret==0)then
          constrainable_coord%is_reaction_coordinate = .false.
        else
          constrainable_coord%is_reaction_coordinate = .true.
        endif
      endif
      if(constrainable_coord%is_reaction_coordinate)then
        if(.not.is_vector_non_zero(direction,3))then
          write(0,*) 'direction is undefined. you must define a direction in order to use the center &
        & of mass as a reaction coordinate.'
          call phase_error_with_msg(nfout,&
          'direction is undefined. you must define a direction in order to use the center of mass as a reaction coordinate.' &
          ,__LINE__,__FILE__)
        endif
        if(printable)write(nfout,'(a,3f15.10)') 'the direction to which the center of mass will change:',&
       & direction(1),direction(2),direction(3)
        ival=0
        if(f_getRealValue(tag_init_value,dret,'bohr')==0)then
          ival = dret
        endif
        fval=1
        if(f_getRealValue(tag_final_value,dret,'bohr')==0)then
          fval = dret
        endif
        incre = 0.1d0
        if(f_getRealValue(tag_increment,dret,'bohr')==0)then
          incre = dret
        endif
        constrainable_coord%n_reaction_coords = abs(nint((fval-ival)/incre))
        allocate(constrainable_coord%reaction_coords(constrainable_coord%n_reaction_coords &
      &        , constrainable_coord%ndim))
        allocate(constrainable_coord%finished(constrainable_coord%n_reaction_coords+1))
        constrainable_coord%finished=.false.
        constrainable_coord%value(:) = ival*direction(:)+constrainable_coord%value(:)
        do i=1,constrainable_coord%n_reaction_coords
          constrainable_coord%reaction_coords(i,:) = constrainable_coord%value(:) + direction(:)*(ival+i*incre)
        enddo
      endif
      iret = f_selectParentBlock()
    endif
  end subroutine prep_reac_coord_com

!  function get_curr_com(constrainable_coord)
!    type(constrainable_coords_t), intent(in) :: constrainable_coord
!    real(DP), dimension(3) :: get_curr_com
!    real(DP), dimension(3) :: res
!    integer :: i,i1
!    res=0.d0
!    sum_mass = 0.d0
!    do i=1,constrainable_coord%n_associated_atoms
!      i1 = constrainable_coord%associated_atoms(i)
!      res(1) = res(1) + amion(ityp(i1))*x_nopbc(i1)
!      res(2) = res(2) + amion(ityp(i1))*y_nopbc(i1)
!      res(3) = res(3) + amion(ityp(i1))*z_nopbc(i1)
!      sum_mass = sum_mass + amion(ityp(i1))
!    enddo
!    res(:) = res(:)/sum_mass
!    get_curr_com = res
!  end function get_curr_com

  subroutine update_sigma_com(constrainable_coord,iidim)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    integer, intent(in) :: iidim
    real(DP),dimension(3) :: curr_com
    call get_curr_com(constrainable_coord,curr_com,sum_mass,-1)
    constrainable_coord%sigma(iidim) = curr_com(iidim)-constrainable_coord%value(iidim)
  end subroutine update_sigma_com

  subroutine update_dsigma_com(constrainable_coord,iidim)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    integer, intent(in) :: iidim
    integer :: i
    character(256) :: str
    do i=1,constrainable_coord%n_associated_atoms
      constrainable_coord%dsigma(iidim,i,iidim) &
    & = amion(ityp(constrainable_coord%associated_atoms(i)))/sum_mass
    enddo
  end subroutine update_dsigma_com

  subroutine monitor_com(constrainable_coord)
    type(constrainable_coords_t), intent(in) :: constrainable_coord
    real(DP),dimension(3) :: curr_com
    character(256) :: str
    integer :: i,iiatom
    write(str,*) constrainable_coord%n_associated_atoms
    call get_curr_com(constrainable_coord,curr_com,sum_mass,-1)
    if(printable)write(nfout,'(a,'//trim(adjustl(str))//'i5,a,3f10.5)') &
    & '    the center of mass for atoms ',                              &
    & (constrainable_coord%associated_atoms(iiatom),                    &
    & iiatom=1,constrainable_coord%n_associated_atoms),                 &
    & ' : ',curr_com(1),curr_com(2),curr_com(3)
  end subroutine monitor_com

end module center_of_mass_constraint

module coord_num_constraint
use m_constraints_data
implicit none

real(DP), private :: kappa
real(DP), private :: rcut
integer :: iiatom
character(len("kappa")), private, parameter     :: tag_kappa = "kappa"
character(len("kappa_inv")), private, parameter :: tag_kappa_inv = "kappa_inv"
character(len("rcut")),  private, parameter :: tag_rcut  = "rcut"
contains
  subroutine read_coord_num_constraint(constrainable_coord)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    integer :: f_getRealValue, f_selectBlock, f_selectParentBlock, f_getIntValue
    integer :: iret
    integer :: i,iat
    real(DP) :: dret
    real(DP) :: fval,dr
    integer :: number_of_reac_coords
    kappa = 5.0d0
    rcut  = 2.0d0
    if(f_selectBlock(tag_coordination_number)==0)then
      if(f_getRealValue(tag_kappa_inv,dret,'bohr')==0)then
        kappa = 1.0d0/(dret)
      endif
      if(f_getRealValue(tag_kappa,dret,'')==0)then
        kappa = dret
      endif
      if(f_getRealValue(tag_rcut,dret,'bohr')==0)then
        rcut = dret
      endif
      iret = f_selectParentBlock()
    endif

    constrainable_coord%nam                = 'coordination number'
    constrainable_coord%descri             = 'constrain the coordination number defined by &
   & the Fermi-distribution function'
    constrainable_coord%unit_name          = ''
    constrainable_coord%typ                = COORDNATION_NUMBER
    constrainable_coord%nfree_decrement    = 1
    constrainable_coord%ndim               = 1
    iiatom = 0
    if(f_getIntValue("atom1",iret)==0)then
      iiatom = iret
      constrainable_coord%n_associated_atoms = natm
    else
      write(0,*) 'you must specify atom1 in order to define a coordination-number constraint.'
      call phase_error_with_msg(nfout,'you must specify atom1 in order to define a coordination-number constraint.' &
                               ,__LINE__,__FILE__)
    endif

    call alloc_constrainable_coord(constrainable_coord)
    do i=1,natm
      constrainable_coord%associated_atoms(i) = i
    enddo
    call parse_mobile_and_monitor(constrainable_coord)

    constrainable_coord%value(1) = get_curr_coord_num(constrainable_coord)

    call prep_reac_coords_1D(constrainable_coord,NO_UNIT,very_small)

  end subroutine read_coord_num_constraint

  real(DP) function get_curr_coord_num(constrainable_coord)
    type(constrainable_coords_t), intent(in) :: constrainable_coord
    integer :: i,iat
    real(DP) :: S,dr
    real(DP) :: tmpval
    S = 0.d0
    do i=1,constrainable_coord%n_associated_atoms
      iat = constrainable_coord%associated_atoms(i)
      if(iat.eq.iiatom)cycle
      tmpval = 1.d0/(dexp(kappa*(get_distance_nopbc(iat,iiatom)-rcut))+1)
      S = S+tmpval
    enddo
    get_curr_coord_num = S
  end function get_curr_coord_num

  subroutine update_sigma_coord_num(constrainable_coord)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    integer :: i,iat
    constrainable_coord%sigma(1) = get_curr_coord_num(constrainable_coord) &
  & - constrainable_coord%value(1)
  end subroutine update_sigma_coord_num

  subroutine update_dsigma_coord_num(constrainable_coord)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    integer :: i,iat
    real(DP), dimension(3) :: vec
    real(DP) :: rlen
    constrainable_coord%dsigma = 0.0d0
    do i=1,constrainable_coord%n_associated_atoms
      iat = constrainable_coord%associated_atoms(i)
      if(iat.eq.iiatom)cycle
      call get_vector_nopbc(iat,iiatom,vec)
      rlen = get_distance_nopbc(iat,iiatom)
      constrainable_coord%dsigma(1,i,:) = -0.5d0*kappa/(dcosh(kappa*(rlen-rcut))+1) *vec(:)/rlen
      constrainable_coord%dsigma(1,iiatom,:) = constrainable_coord%dsigma(1,iiatom,:) &
     & + 0.5d0*kappa/(dcosh(kappa*(rlen-rcut))+1) *vec(:)/rlen
    enddo
  end subroutine update_dsigma_coord_num

  subroutine monitor_coord_num(constrainable_coord)
    type(constrainable_coords_t), intent(in) :: constrainable_coord
    real(DP) :: coordi
    coordi = get_curr_coord_num(constrainable_coord)
    if(printable)write(nfout,'(a,i5,a,f10.5)') '    coordination number of atom ',iiatom,' : ',coordi
  end subroutine monitor_coord_num

end module coord_num_constraint

module bond_length_diff_constraint
use m_constraints_data
implicit none
contains
  subroutine read_bl_diff_constraint(constrainable_coord)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    constrainable_coord%nam                = 'bond_length_diff'
    constrainable_coord%descri             = 'constrain the difference of the bond-length of either three or four atoms'
    constrainable_coord%unit_name          = 'bohr'
    constrainable_coord%typ                = BOND_LENGTH_DIFF
    constrainable_coord%nfree_decrement    = 1
    constrainable_coord%ndim               = 1

    call parse_associated_atoms(constrainable_coord)
    if(constrainable_coord%n_associated_atoms.lt.3)then
      write(0,'(a)') 'you must define at least 3 atoms in order to define a bond-length-diff constraint.'
      call phase_error_with_msg(nfout,'you must define at least 3 atoms in order to define a bond-length-diff constraint.'&
                               ,__LINE__,__FILE__)
    endif
    call alloc_constrainable_coord(constrainable_coord)
    call parse_mobile_and_monitor(constrainable_coord)

    constrainable_coord%value(1) = get_curr_bl_diff(constrainable_coord)
    call prep_reac_coords_1D(constrainable_coord,LENGTH)
  end subroutine read_bl_diff_constraint

  function get_indices(constrainable_coord)
    type(constrainable_coords_t), intent(in) :: constrainable_coord
    real(DP), dimension(4) :: get_indices
    get_indices(1)=constrainable_coord%associated_atoms(1)
    get_indices(2)=constrainable_coord%associated_atoms(2)
    if(constrainable_coord%n_associated_atoms.eq.3)then
      get_indices(3) = get_indices(2)
      get_indices(4) = constrainable_coord%associated_atoms(3)
    else
      get_indices(3) = constrainable_coord%associated_atoms(3)
      get_indices(4) = constrainable_coord%associated_atoms(4)
    endif
  end function get_indices

  real(DP) function get_curr_bl_diff(constrainable_coord)
    type(constrainable_coords_t), intent(in) :: constrainable_coord
    integer, dimension(4) :: ii
    real(DP) :: r1,r2
    ii = get_indices(constrainable_coord)
    r1 = get_distance_nopbc(ii(1),ii(2))
    r2 = get_distance_nopbc(ii(3),ii(4))
    get_curr_bl_diff = r1-r2
  end function get_curr_bl_diff

  subroutine update_sigma_bl_diff(constrainable_coord)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    real(DP) :: currvalue
    currvalue = get_curr_bl_diff(constrainable_coord)
    constrainable_coord%sigma(1) = currvalue - constrainable_coord%value(1)
  end subroutine update_sigma_bl_diff

  subroutine update_dsigma_bl_diff(constrainable_coord)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    real(DP), dimension(3) :: r12,r34
    real(DP) :: r12l,r34l
    integer, dimension(4) :: indices
    indices  = get_indices(constrainable_coord)
    call get_vector_nopbc(indices(1),indices(2),r12)
    call get_vector_nopbc(indices(3),indices(4),r34)
    r12l = get_distance_nopbc(indices(1),indices(2))
    r34l = get_distance_nopbc(indices(3),indices(4))
    r12(:) = r12(:)/r12l
    r34(:) = r34(:)/r34l
    constrainable_coord%dsigma(1,1,:) =  r12(:)
    constrainable_coord%dsigma(1,2,:) = -r12(:)
    if(constrainable_coord%n_associated_atoms.eq.3)then
      constrainable_coord%dsigma(1,2,:) = constrainable_coord%dsigma(1,2,:)-r34(:)
      constrainable_coord%dsigma(1,3,:) = r34(:)
    else
      constrainable_coord%dsigma(1,3,:) = -r34(:)
      constrainable_coord%dsigma(1,4,:) = +r34(:)
    endif
  end subroutine update_dsigma_bl_diff

  subroutine monitor_bl_diff(constrainable_coord)
    type(constrainable_coords_t), intent(in) :: constrainable_coord
    integer, dimension(4) :: ind
    real(DP) :: value
    ind = get_indices(constrainable_coord)
    value = get_curr_bl_diff(constrainable_coord)
    if(printable)write(nfout,'(a,i5,i5,a,i5,i5,a,f10.5,a)') &
   & '    the difference of the bond-length between atom ',ind(1),ind(2), &
   & ' and between atom ',ind(3),ind(4),' : ',value,' bohr'
  end subroutine monitor_bl_diff

end module bond_length_diff_constraint

module bond_length_sum_constraint
use m_constraints_data
implicit none
contains
  subroutine read_bl_sum_constraint(constrainable_coord)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    constrainable_coord%nam                = 'bond_length_sum'
    constrainable_coord%descri             = 'constrain the sum of the bond-length of either three or four atoms'
    constrainable_coord%unit_name          = 'bohr'
    constrainable_coord%typ                = BOND_LENGTH_SUM
    constrainable_coord%nfree_decrement    = 1
    constrainable_coord%ndim               = 1

    call parse_associated_atoms(constrainable_coord)
    if(constrainable_coord%n_associated_atoms.lt.3)then
      write(0,'(a)') 'you must define at least 3 atoms in order to define a bond-length-sum constraint.'
      call phase_error_with_msg(nfout,'you must define at least 3 atoms in order to define a bond-length-sum constraint.'&
                               ,__LINE__,__FILE__)
    endif
    call alloc_constrainable_coord(constrainable_coord)
    call parse_mobile_and_monitor(constrainable_coord)

    constrainable_coord%value(1) = get_curr_bl_sum(constrainable_coord)
    call prep_reac_coords_1D(constrainable_coord,LENGTH,very_small)
  end subroutine read_bl_sum_constraint

  function get_indices(constrainable_coord)
    type(constrainable_coords_t), intent(in) :: constrainable_coord
    real(DP), dimension(4) :: get_indices
    get_indices(1)=constrainable_coord%associated_atoms(1)
    get_indices(2)=constrainable_coord%associated_atoms(2)
    if(constrainable_coord%n_associated_atoms.eq.3)then
      get_indices(3) = get_indices(2)
      get_indices(4) = constrainable_coord%associated_atoms(3)
    else
      get_indices(3) = constrainable_coord%associated_atoms(3)
      get_indices(4) = constrainable_coord%associated_atoms(4)
    endif
  end function get_indices

  real(DP) function get_curr_bl_sum(constrainable_coord)
    type(constrainable_coords_t), intent(in) :: constrainable_coord
    integer, dimension(4) :: ii
    real(DP) :: r1,r2
    ii = get_indices(constrainable_coord)
    r1 = get_distance_nopbc(ii(1),ii(2))
    r2 = get_distance_nopbc(ii(3),ii(4))
    get_curr_bl_sum = r1+r2
  end function get_curr_bl_sum

  subroutine update_sigma_bl_sum(constrainable_coord)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    real(DP) :: currvalue
    currvalue = get_curr_bl_sum(constrainable_coord)
    constrainable_coord%sigma(1) = currvalue - constrainable_coord%value(1)
  end subroutine update_sigma_bl_sum

  subroutine update_dsigma_bl_sum(constrainable_coord)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    real(DP), dimension(3) :: r12,r34
    real(DP) :: r12l,r34l
    integer, dimension(4) :: indices
    indices  = get_indices(constrainable_coord)
    call get_vector_nopbc(indices(1),indices(2),r12)
    call get_vector_nopbc(indices(3),indices(4),r34)
    r12l = get_distance_nopbc(indices(1),indices(2))
    r34l = get_distance_nopbc(indices(3),indices(4))
    r12(:) = r12(:)/r12l
    r34(:) = r34(:)/r34l
    constrainable_coord%dsigma(1,1,:) =  r12(:)
    constrainable_coord%dsigma(1,2,:) = -r12(:)
    if(constrainable_coord%n_associated_atoms.eq.3)then
      constrainable_coord%dsigma(1,2,:) = constrainable_coord%dsigma(1,2,:)+r34(:)
      constrainable_coord%dsigma(1,3,:) = -r34(:)
    else
      constrainable_coord%dsigma(1,3,:) = +r34(:)
      constrainable_coord%dsigma(1,4,:) = -r34(:)
    endif
  end subroutine update_dsigma_bl_sum

  subroutine monitor_bl_sum(constrainable_coord)
    type(constrainable_coords_t), intent(in) :: constrainable_coord
    integer, dimension(4) :: ind
    real(DP) :: value
    ind = get_indices(constrainable_coord)
    value = get_curr_bl_sum(constrainable_coord)
    if(printable)write(nfout,'(a,i5,i5,a,i5,i5,a,f10.5,a)') &
   & '    the sum of the bond-length between atom ',ind(1),ind(2), &
   & ' and between atom ',ind(3),ind(4),' : ',value,' bohr'
  end subroutine monitor_bl_sum

end module bond_length_sum_constraint

#ifdef USER_DEFINED_CONSTRAINT
module user_defined_constraint
  use m_constraints_data
  implicit none

  contains

  subroutine read_udef_constraint(constrainable_coord)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    integer :: udefc_get_nfree_decrement, udefc_get_ndim, udefc_get_n_associated
    character(30) :: udefc_get_name
    character(10) :: udefc_get_unit_name
    character(256) :: udefc_get_description
    real(DP), allocatable, dimension(:) :: valtmp
    integer, allocatable, dimension(:) :: atomtmp
    integer :: i,iunit
    integer :: f_getIntValue, f_getRealValue, f_getStringValue, f_selectBlock, f_selectTop, &
   & f_selectParentBlock, f_selectNextTableLine, f_selectFirstTableLine
    external f_getIntValue, f_getRealValue, f_getStringValue, f_selectBlock, f_selectTop, &
   & f_selectParentBlock, f_selectNextTableLine, f_selectFirstTableLine

    constrainable_coord%nam                = udefc_get_name()
    constrainable_coord%descri             = udefc_get_description()
    constrainable_coord%unit_name          = udefc_get_unit_name()
    constrainable_coord%typ                = USER_DEFINED
    constrainable_coord%nfree_decrement    = udefc_get_nfree_decrement()
    constrainable_coord%ndim               = udefc_get_ndim()
    constrainable_coord%n_associated_atoms = udefc_get_n_associated()

    call parse_associated_atoms(constrainable_coord)
    call alloc_constrainable_coord(constrainable_coord)
    call parse_mobile_and_monitor(constrainable_coord)

    allocate(valtmp(constrainable_coord%ndim));valtmp=0.d0
    allocate(atomtmp(constrainable_coord%n_associated_atoms));atomtmp=constrainable_coord%associated_atoms

    call udefc_get_curr_value( x_nopbc,y_nopbc,z_nopbc,natm, &
    & atomtmp, constrainable_coord%n_associated_atoms, valtmp, constrainable_coord%ndim )
    do i=1,constrainable_coord%ndim
      constrainable_coord%value(i) = valtmp(i)
    enddo
    deallocate(valtmp)
    deallocate(atomtmp)

    if(constrainable_coord%ndim==1)then
      iunit = NO_UNIT
      if(constrainable_coord%unit_name=='bohr') iunit=LENGTH
      if(constrainable_coord%unit_name=='radian') iunit=ANGLE
      call prep_reac_coords_1D(constrainable_coord,iunit)
    endif

    if(.not.associated(constrainable_coord%reaction_coords)) then
      allocate(constrainable_coord%reaction_coords(1,1)) !dummy allocation
      allocate(constrainable_coord%finished(1)) !dummy allocation
    endif

    call udefc_parse_input(f_getIntValue, f_getRealValue, f_getStringValue,  &
   & f_selectBlock, f_selectParentBlock, f_selectTop, f_selectNextTableLine, &
   !!& f_selectFirstTableLine, constrainable_coord%reaction_coords)
   & f_selectFirstTableLine)

  end subroutine read_udef_constraint

  subroutine update_sigma_udef_constraint(constrainable_coord, iidim)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    integer, intent(in) :: iidim
    integer :: i
    real(DP), allocatable, dimension(:) :: valtmp
    integer, allocatable, dimension(:) :: attmp
    allocate(valtmp(constrainable_coord%ndim));valtmp=0.d0
    allocate(attmp(constrainable_coord%n_associated_atoms))
    attmp=constrainable_coord%associated_atoms
    call udefc_get_curr_value( x_nopbc,y_nopbc,z_nopbc,natm,attmp,&
   & constrainable_coord%n_associated_atoms,valtmp,constrainable_coord%ndim)

    do i=1,constrainable_coord%ndim
      constrainable_coord%sigma(i) = &
    & valtmp(i)-constrainable_coord%value(i)
    enddo
    deallocate(valtmp)
    deallocate(attmp)
  end subroutine update_sigma_udef_constraint

  subroutine update_dsigma_udef_constraint(constrainable_coord, iidim)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    integer, intent(in) :: iidim
    integer :: i,j
    real(DP), allocatable, dimension(:,:,:) :: dsigtmp
    integer, allocatable, dimension(:) :: attmp
    allocate(dsigtmp(constrainable_coord%ndim,constrainable_coord%n_associated_atoms,3));dsigtmp=0.d0
    allocate(attmp(constrainable_coord%n_associated_atoms));attmp=constrainable_coord%associated_atoms
    constrainable_coord%dsigma=0.d0
    call udefc_get_dsigma(                           &
   &  x_nopbc,y_nopbc,z_nopbc,natm,                  &
   &  attmp, constrainable_coord%n_associated_atoms, &
   &  constrainable_coord%ndim,dsigtmp)
    do i=1,constrainable_coord%n_associated_atoms
      do j=1,constrainable_coord%ndim
        constrainable_coord%dsigma(j,i,:) = dsigtmp(j,i,:)
      enddo
    enddo
    deallocate(dsigtmp)
    deallocate(attmp)
  end subroutine update_dsigma_udef_constraint

  subroutine monitor_udef_constraint(constrainable_coord)
    type(constrainable_coords_t), intent(in) :: constrainable_coord
    integer, allocatable, dimension(:) :: atoms
    allocate(atoms(constrainable_coord%n_associated_atoms))
    atoms = constrainable_coord%associated_atoms
    call udefc_monitor(nfout,x_nopbc,y_nopbc,z_nopbc,natm,atoms,constrainable_coord%n_associated_atoms,&
   & constrainable_coord%ndim,printable)
    deallocate(atoms)
  end subroutine monitor_udef_constraint

end module user_defined_constraint
#endif

module bond_angle_diff_constraint
use m_constraints_data
implicit none
contains

  subroutine read_bang_diff_constraint(constrainable_coord)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    constrainable_coord%nam                = 'bond_angle_diff'
    constrainable_coord%descri             = 'constrain the difference between two bond-angles'
    constrainable_coord%unit_name          = 'radian'
    constrainable_coord%typ                = BOND_ANGLE_DIFF
    constrainable_coord%nfree_decrement    = 1
    constrainable_coord%ndim               = 1
    constrainable_coord%n_associated_atoms = 6
    call parse_associated_atoms(constrainable_coord)
    call alloc_constrainable_coord(constrainable_coord)
    call parse_mobile_and_monitor(constrainable_coord)
    if(constrainable_coord%mobile) constrainable_coord%nfree_decrement=0
    constrainable_coord%value(1) = get_curr_bang_diff(constrainable_coord)
    call prep_reac_coords_1D(constrainable_coord,ANGLE)
  end subroutine read_bang_diff_constraint

  subroutine monitor_bang_diff(constrainable_coord)
    type(constrainable_coords_t), intent(in) :: constrainable_coord
    integer :: i1, i2, i3, i4, i5, i6
    real(DP) bang
    bang = get_curr_bang_diff(constrainable_coord)
    i1 = constrainable_coord%associated_atoms(1)
    i2 = constrainable_coord%associated_atoms(2)
    i3 = constrainable_coord%associated_atoms(3)
    i4 = constrainable_coord%associated_atoms(4)
    i5 = constrainable_coord%associated_atoms(5)
    i6 = constrainable_coord%associated_atoms(6)
    if(printable)write(nfout,'(a,i5,a,i5,a,i5,a,i5,a,i5,a,i5,a,f10.5,a)') &
  & '    the difference of the bond-angle for atoms ',i1,', ',i2,',',i3,  &
  & ' and atoms ',i4,',',i5,',',i6,&
  & ' : ',180.0d0*bang/PAI,' degrees'
  end subroutine monitor_bang_diff

  real(DP) function get_curr_bang_diff(constrainable_coord)
    type(constrainable_coords_t), intent(in) :: constrainable_coord
    integer :: i1,i2,i3
    real(DP) :: ang1, ang2
    ang1 = get_bang(1,2,3,constrainable_coord)
    ang2 = get_bang(4,5,6,constrainable_coord)
    get_curr_bang_diff = ang1-ang2
  end function get_curr_bang_diff

  subroutine update_sigma_bang_diff(constrainable_coord)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    real(DP) :: currbangdiff
    currbangdiff = get_curr_bang_diff(constrainable_coord)
    constrainable_coord%sigma(1) = currbangdiff-constrainable_coord%value(1)
  end subroutine update_sigma_bang_diff

  subroutine update_dsigma_bang_diff(constrainable_coord)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    constrainable_coord%dsigma = 0.d0
    call dsigma_bang_diff(1,2,3, 1.0d0,constrainable_coord)
    call dsigma_bang_diff(4,5,6,-1.0d0,constrainable_coord)
  end subroutine update_dsigma_bang_diff

  real(DP) function get_bang(at1,at2,at3,constrainable_coord)
    integer, intent(in) :: at1, at2, at3
    type(constrainable_coords_t), intent(in) :: constrainable_coord
    integer :: i1,i2,i3
    real(DP) :: rji,rki,cosine
    real(DP), dimension(3) :: vecji,vecki
    real(DP) :: eps=1.d-12
    i1 = constrainable_coord%associated_atoms(at1)
    i2 = constrainable_coord%associated_atoms(at2)
    i3 = constrainable_coord%associated_atoms(at3)
    rji = get_distance_nopbc(i1,i2)
    rki = get_distance_nopbc(i3,i2)
    call get_vector_nopbc(i1,i2,vecji)
    call get_vector_nopbc(i3,i2,vecki)
    cosine = dot_product(vecji,vecki)/(rji*rki)
    if(cosine.gt.1)  cosine =  1.d0
    if(cosine.lt.-1) cosine = -1.d0
    get_bang = dacos(cosine)
  end function get_bang

  subroutine dsigma_bang_diff(at1,at2,at3,fac,constrainable_coord)
    integer, intent(in) :: at1, at2, at3
    real(DP), intent(in) :: fac
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    real(DP) :: cosine,sine,rji,rki,theta
    real(DP),dimension(3) :: vecji,vecki
    real(DP), dimension(3) :: der1,der3
    integer :: i1,i2,i3
    i1 = constrainable_coord%associated_atoms(at1)
    i2 = constrainable_coord%associated_atoms(at2)
    i3 = constrainable_coord%associated_atoms(at3)
    rji = get_distance_nopbc(i1,i2)
    rki = get_distance_nopbc(i3,i2)
    call get_vector_nopbc(i1,i2,vecji)
    call get_vector_nopbc(i3,i2,vecki)
    theta = get_bang(at1,at2,at3,constrainable_coord)
    cosine = dcos(theta)
    sine   = dsin(theta)
    der1(:) = ((vecji(:)/rji)*cosine - (vecki(:)/rki))/(rji*sine)
    der3(:) = ((vecki(:)/rki)*cosine - (vecji(:)/rji))/(rki*sine)
    constrainable_coord%dsigma(1,at1,:) = constrainable_coord%dsigma(1,at1,:) &
   & + der1(:)*fac
    constrainable_coord%dsigma(1,at3,:) = constrainable_coord%dsigma(1,at3,:) &
   & + der3(:)*fac
    constrainable_coord%dsigma(1,at2,:) = constrainable_coord%dsigma(1,at2,:) &
   & -(der1(:)+der3(:))*fac
  end subroutine dsigma_bang_diff

end module bond_angle_diff_constraint

module distance_from_ref_constraint
use m_constraints_data
use m_Ionic_System, only : altv
implicit none
integer,private,dimension(3) :: exclude
contains
  subroutine update_sigma_dref(constrainable_coord)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    real(DP) :: currdistance
    constrainable_coord%sigma(1) = &
    get_curr_distance_from_ref(constrainable_coord)-constrainable_coord%value(1)
  end subroutine update_sigma_dref

  subroutine update_dsigma_dref(constrainable_coord)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    integer :: i,it
    real(DP) :: sigma,isigma
    sigma = get_curr_distance_from_ref(constrainable_coord)
    if(sigma<very_small) then
      call phase_error_with_msg(nfout, &
      'current structure is too close to the reference structure',  &
      __LINE__, __FILE__)
    endif
    isigma = 1.d0/sigma
    do i=1,constrainable_coord%n_associated_atoms
      it = constrainable_coord%associated_atoms(i)
      constrainable_coord%dsigma(1,i,1) = isigma*(cps(it,1)-constrainable_coord%refx(i))
      constrainable_coord%dsigma(1,i,2) = isigma*(cps(it,2)-constrainable_coord%refy(i))
      constrainable_coord%dsigma(1,i,3) = isigma*(cps(it,3)-constrainable_coord%refz(i))
    enddo
  end subroutine update_dsigma_dref

  subroutine monitor_dref(constrainable_coord)
    type(constrainable_coords_t), intent(in) :: constrainable_coord
    integer :: i1,i2
    real(DP) dref
    dref = get_curr_distance_from_ref(constrainable_coord)
    write(nfout,'(a,f10.5)') 'distance from the reference structure ', dref
  end subroutine monitor_dref

  real(DP) function get_curr_distance_from_ref(constrainable_coord)
    type(constrainable_coords_t), intent(in) :: constrainable_coord
    integer :: i,it
    real(DP) :: ret
    ret = 0.d0
    do i=1,constrainable_coord%n_associated_atoms
      it = constrainable_coord%associated_atoms(i)
      ret = ret+(constrainable_coord%refx(i)-cps(it,1))**2 &
               +     (constrainable_coord%refy(i)-cps(it,2))**2 &
               +     (constrainable_coord%refz(i)-cps(it,3))**2
    enddo
    get_curr_distance_from_ref = sqrt(ret)
  end function get_curr_distance_from_ref

  subroutine read_dref_constraint(constrainable_coord)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    integer :: iret,iiret,i
    integer :: f_getIntValue,f_selectBlock,f_selectParentBlock,f_getRealValue
    real(DP) :: dret
    real(DP) :: dr,fval
    integer :: number_of_reac_coords
    character(10) strx,stry,strz

    constrainable_coord%nam                = 'distance_from_ref'
    constrainable_coord%descri             = 'constrain the distance from the specified reference structure'
    constrainable_coord%unit_name          = 'bohr'
    constrainable_coord%typ                = DISTANCE_FROM_REF
    constrainable_coord%nfree_decrement    = 1
    constrainable_coord%ndim               = 1

    call rd_reference_structure(constrainable_coord)
    call alloc_constrainable_coord(constrainable_coord,allocate_map = .false.)
    call parse_mobile_and_monitor(constrainable_coord)
    constrainable_coord%value(1) = get_curr_distance_from_ref(constrainable_coord)
    call prep_reac_coords_1D(constrainable_coord,LENGTH,very_small)

  end subroutine read_dref_constraint

  subroutine rd_reference_structure(constrainable_coord)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    integer :: f_getIntValue,f_selectBlock,f_selectParentBlock,f_getRealValue, f_getStringValue
    integer :: f_selectFirstTableLine, f_selectNextTableLine
    integer :: iret
    character(len=FMAXVALLEN) :: rstr
    logical :: tf
    integer :: coordsys,icount, iat, nat
    character(len=FMAXUNITLEN) :: unit_f
    real(DP), allocatable, dimension(:,:) :: tmpcoord,tmpcoord2
    real(DP) :: rx,ry,rz
    coordsys = PUCV
    if(f_selectBlock(tag_reference_structure)==0) then
      if(f_getStringValue(tag_coordinate_system,rstr,LOWER)==0) then
        call strncmp0(trim(rstr),tag_cartesian,tf)
        if(.not.tf) call strncmp0(trim(rstr),tag_XYZ,tf)
        if(tf) coordsys = CARTS
      endif
      if(f_selectBlock(tag_atoms)==0) then
        unit_f = ''
        if(coordsys == CARTS) unit_f = 'bohr'
        icount = 1
        do while(.true.)
           if (icount == 1) then
              if(f_selectFirstTableLine() /= 0) then
                 exit
              end if
           else
              if(f_selectNextTableLine() /= 0) then
                 exit
              end if
           end if
           icount = icount+1
        enddo
        constrainable_coord%n_associated_atoms = icount-1
        allocate(constrainable_coord%refx(constrainable_coord%n_associated_atoms))
        allocate(constrainable_coord%refy(constrainable_coord%n_associated_atoms))
        allocate(constrainable_coord%refz(constrainable_coord%n_associated_atoms))
        allocate(tmpcoord(constrainable_coord%n_associated_atoms,3))
        allocate(constrainable_coord%associated_atoms(constrainable_coord%n_associated_atoms))
        icount = 1
        do while(.true.)
           if (icount == 1) then
              if(f_selectFirstTableLine() /= 0) then
                 exit
              end if
           else
              if(f_selectNextTableLine() /= 0) then
                 exit
              end if
           end if
           if (icount>natm) then
             call phase_error_with_msg(nfout, &
             'the number of atoms for the reference and base structures must be equal',  &
             __LINE__, __FILE__)
           endif
           iat = icount
           if(f_getIntValue(tag_id,iret)==0) iat = iret
           if(f_getIntValue(tag_no,iret)==0) iat = iret
           constrainable_coord%associated_atoms(icount) = iat
           if(f_getRealValue(tag_rx, rx, unit_f)==0) tmpcoord(icount,1) = rx
           if(f_getRealValue(tag_ry, ry, unit_f)==0) tmpcoord(icount,2) = ry
           if(f_getRealValue(tag_rz, rz, unit_f)==0) tmpcoord(icount,3) = rz
           icount = icount+1
        enddo
        iret = f_selectParentBlock()
!        if(icount<=natm) then
!          call phase_error_with_msg(nfout, &
!         'the number of atoms for the reference and base structures must be equal',  &
!          __LINE__, __FILE__)
!        endif
      endif
      iret = f_selectParentBlock()
    endif
    if(coordsys==PUCV) then
      nat = constrainable_coord%n_associated_atoms
      allocate(tmpcoord2(nat,3))
      call change_of_coordinate_system(altv,tmpcoord,nat,nat,tmpcoord2)!-(b_I.S.) pos -> cps
      constrainable_coord%refx(:) = tmpcoord2(:,1)
      constrainable_coord%refy(:) = tmpcoord2(:,2)
      constrainable_coord%refz(:) = tmpcoord2(:,3)
      deallocate(tmpcoord2)
    else
      constrainable_coord%refx(:) = tmpcoord(:,1)
      constrainable_coord%refy(:) = tmpcoord(:,2)
      constrainable_coord%refz(:) = tmpcoord(:,3)
    endif
    if(printable)then
      write(nfout,'(a)') '!** reference atomic coordinates'
      do iat=1,constrainable_coord%n_associated_atoms
        write(nfout,'(i8,3f20.10)') constrainable_coord%associated_atoms(iat), &
                                    constrainable_coord%refx(iat), &
                                    constrainable_coord%refy(iat), &
                                    constrainable_coord%refz(iat)
      enddo
    endif
    deallocate(tmpcoord)
  end subroutine rd_reference_structure

end module distance_from_ref_constraint

module m_constraints

use m_Const_Parameters,   only : LOWER, UPPER, NOCONV, DP, WITHOUTTAG, T_CONTROL, VERLET, VELOCITY_SCALING, INITIAL, Hartree &
                               , CONST_EV, CONST_NA, CONST_CALth
use m_Control_Parameters, only : printable, ipri, af, m_CtrlP_wd_isolver, icond, icond_org
use m_Parallelization, only : mype,npes,MPI_CommGroup,nrank_conf,mype_conf,conf_para
use m_Ionic_System, only : amion,ityp,natm,natm2,iwei,dtio,m_IS_wd_speciesname_etc,imdalg,napt,m_IS_wd_pos_and_v,cpd_l &
&                        , lattice_system_from_m_CS_SG, nopr_supercell, napt_supercell, iop_supercell
use m_Crystal_Structure, only : nopr, op, rltv, b2pmat, tau
use m_Files, only : F_ENF, F_DYNM, F_ZAJ, F_CNTN, F_CNTN_BIN, F_CHGT &
  &, F_ZAJ_in,F_CNTN_in,F_CNTN_BIN_in,F_CHGT_in,nfcntn, m_Files_reopen_nfcntn
use m_variables_for_atoms, only : Nfree,mobile,Nreservoir
use m_variables_for_dynamics, only : curr_md_step,mdstep_for_this_run,forcmx,is_md
use m_Total_Energy, only : etotal
use m_Force, only : m_Force_wd_force_and_cps
use m_IterationNumbers, only : m_Iter_wd_iteration_numbers, iteration_ionic, iteration
use m_Total_Energy, only : m_TE_wd_total_energy

use m_mtrandom

use distance_from_pos_constraint
use distance_from_com_constraint
use bond_length_constraint
use bond_angle_constraint
use dihedral_angle_constraint
use plane_constraint
use center_of_mass_constraint
use coord_num_constraint
use bond_length_diff_constraint
use bond_length_sum_constraint
#ifdef USER_DEFINED_CONSTRAINT
use user_defined_constraint
#endif
use bond_angle_diff_constraint
use distance_from_ref_constraint
use mpi

implicit none

!include 'mpif.h'

character(len('constrainable_coords')), parameter, private :: tag_constrainable_coords='constrainable_coords'
type(constrainable_coords_t), private, allocatable, dimension(:) &
 & :: constrainable_coords
integer, private :: n_constrainable_coords

integer, private :: n_cnstr_src=0
integer, private :: n_cnstr_coords_inp

real(DP), private :: eps_rattle = 1.d-10
integer,  private :: max_iteration_rattle = 10000

integer, private :: n_total_reac_coords = 1
integer, private :: reac_coord_id = 0
integer, private :: ifile_handle
integer :: ifile_handle_structure
integer :: ifile_handle_energy
integer, private :: ifile_handle_structure_reac
integer, private :: ifile_handle_energy_reac
integer, private :: ifile_handle_reac
logical, private :: convfile_opened=.false.

logical, allocatable, dimension(:) :: is_constrained_atom

real(DP), allocatable, dimension(:) :: force_of_constraintx
real(DP), allocatable, dimension(:) :: force_of_constrainty
real(DP), allocatable, dimension(:) :: force_of_constraintz

real(DP), parameter :: Hartree2eV  = Hartree
!real(DP), parameter :: autime2fs   = AU_TIME * 1.0d15  ! 2.418884327d-2
real(DP), parameter :: eV2kJ_mol   = CONST_EV*CONST_NA  ! 2625.5d0/Hartree2eV  ~96.4853
real(DP), parameter :: eV2kcal_mol = eV2kJ_mol / CONST_CALth  ! 627.509391d0/Hartree2eV  ~23.0605

logical :: reset_md_step

integer, private :: istart_reac_coord=1

character(len=260) :: F_CNTN_org,F_CNTN_BIN_org,F_ZAJ_org,F_CHGT_org
character(len=260) :: F_REAC='./reac_coords.data'

logical :: contfile_read=.false.
logical :: skip_coords=.false.

contains

  integer function m_cnstr_get_n_reac_coords()
    m_cnstr_get_n_reac_coords=n_total_reac_coords
  end function m_cnstr_get_n_reac_coords

  integer function m_cnstr_get_reac_coords_s()
    m_cnstr_get_reac_coords_s=istart_reac_coord
  end function m_cnstr_get_reac_coords_s

  subroutine m_cnstr_set_ncnstr_src(nn)
    integer, intent(in) :: nn
    n_cnstr_src = nn
  end subroutine m_cnstr_set_ncnstr_src

  logical function m_cnstr_reac_coords_variable()
    if(n_total_reac_coords.gt.1)then
      m_cnstr_reac_coords_variable = .true.
    else
      m_cnstr_reac_coords_variable = .false.
    endif
  end function m_cnstr_reac_coords_variable

  character(len=10) function get_id_char(imax,i)
    integer, intent(in) :: imax,i
    integer :: iketa=2
    integer :: itmp
    character(len=256) :: ch
    character(len=256) :: string
    itmp=int(log10(real(imax)))+1
    if(itmp>2)iketa=itmp
    write(ch,*) iketa
    write(get_id_char,'(i'//trim(adjustl(ch))//'.'//trim(adjustl(ch))//')') i
  end function get_id_char

  character(len=10) function m_cnstr_get_reac_coord_char()
    m_cnstr_get_reac_coord_char = get_id_char(n_total_reac_coords,reac_coord_id)
  end function m_cnstr_get_reac_coord_char

  integer function m_cnstr_get_id()
    m_cnstr_get_id = reac_coord_id
  end function m_cnstr_get_id

!TODO file no atukaikata nado wo hannyousei no aru keisikini
  ! will return .true. if the corresponding reaction coordinate should be 'skipped'.
  logical function m_cnstr_set_reac_coord(ireac)
    integer, intent(in) :: ireac
    character(len=10) :: cid,cid0
    logical :: ex,logi
    integer :: i

    reac_coord_id=ireac

    if(printable.and.m_cnstr_reac_coords_variable()) &
    &    write(nfout,'(a,i2.2)') "entering reaction coordinate no. ",reac_coord_id

    m_cnstr_set_reac_coord=.false.
    if (conf_para) then !preparation for the configuration-parallelization scheme.
      cid = m_cnstr_get_reac_coord_char()
      F_ZAJ=trim(F_ZAJ_org)//'_reac'//trim(adjustl(cid))
      F_CNTN=trim(F_CNTN_org)//'_reac'//trim(adjustl(cid))
      F_CNTN_BIN=trim(F_CNTN_BIN_org)//'_reac'//trim(adjustl(cid))
      F_CHGT=trim(F_CHGT_org)//'_reac'//trim(adjustl(cid))
      F_ZAJ_in=trim(F_ZAJ_org)//'_reac'//trim(adjustl(cid))
      F_CNTN_in=trim(F_CNTN_org)//'_reac'//trim(adjustl(cid))
      F_CNTN_BIN_in=trim(F_CNTN_BIN_org)//'_reac'//trim(adjustl(cid))
      F_CHGT_in=trim(F_CHGT_org)//'_reac'//trim(adjustl(cid))
      if (icond==CONTINUATION) then
        inquire(file=trim(F_CNTN),exist=ex)
        if (ex) then !! the continuation file for the current reaction coordinate
          if(printable) &
        & write(nfout,'(a)') 'found continuation file : '//trim(F_CNTN)//' for reaction coordinate '//&
        & trim(adjustl(cid))//'; loading data.'
         call prepare_continuation_of_reac(nfcntn,.true.)
         if(m_cnstr_finished(reac_coord_id))then
           m_cnstr_set_reac_coord=.true.
           contfile_read=.false.
         else
           call Ewald_and_Structure_Factor()
           call scf_rd_wf_and_chg(logi)
           contfile_read=.true.
         endif
        else !! look for the 'closest' continuation file
          if(.not.contfile_read)then !! the operation will be done only for the first reac-coord of the run.
            do i=ireac-1,1,-1
              cid0 = get_id_char(n_total_reac_coords,i)
              inquire(file=trim(F_CNTN_org)//'_reac'//trim(adjustl(cid0)),exist=ex)
              if(ex)then
                if(printable) &
                &   write(nfout,'(a)') &
                &  'found (the closest) continuation file : '//trim(F_CNTN_org)//'_reac'//trim(adjustl(cid0)) &
                &  //' for reaction coordinate '//trim(adjustl(cid))//'; loading data.'
                F_CNTN=trim(F_CNTN_org)//'_reac'//trim(adjustl(cid0))
                F_ZAJ=trim(F_ZAJ_org)//'_reac'//trim(adjustl(cid0))
                F_CHGT=trim(F_CHGT_org)//'_reac'//trim(adjustl(cid0))
                call prepare_continuation_of_reac(nfcntn,.false.)
                call Ewald_and_Structure_Factor()
                call scf_rd_wf_and_chg(logi)
                F_CNTN=trim(F_CNTN_org)//'_reac'//trim(adjustl(cid))
                F_ZAJ=trim(F_ZAJ_org)//'_reac'//trim(adjustl(cid))
                F_CHGT=trim(F_CHGT_org)//'_reac'//trim(adjustl(cid))
                iteration = 0
                reset_md_step=.true.
                contfile_read=.true.
                skip_coords=.false.
                exit
              endif
            enddo
          endif
        endif
      endif
    else
!      if (icond==CONTINUATION) then
!        inquire(file=trim(F_CNTN),exist=ex)
!        if (ex) then !! the continuation file for the current reaction coordinate
!           call prepare_continuation_of_reac(nfcntn,.true.)
!           call Ewald_and_Structure_Factor()
!           call scf_rd_wf_and_chg(logi)
!        endif
!      endif
    endif

    if(m_cnstr_set_reac_coord)return

    if(.not.reset_md_step)then
      reset_md_step = .true.
    else
      curr_md_step = 0
      iteration_ionic=curr_md_step+1
    endif
    call m_cnstr_update_reac_coords(reac_coord_id)
  end function m_cnstr_set_reac_coord

  logical function m_cnstr_finished(ireac)
    integer, intent(in) :: ireac
    integer :: i
    m_cnstr_finished=.false.
    if(nvariable.eq.1.and.reac_coord_generation==REAC_COORD_VIA_INPUT)then
      do i=1,n_constrainable_coords
        if(constrainable_coords(i)%is_reaction_coordinate)then
          m_cnstr_finished = constrainable_coords(i)%finished(ireac)
          if(m_cnstr_finished.and.printable) write(nfout,'(a,i5,a)') '  reaction coordinate no. ' &
      & ,ireac,' has already been taken into account.'
          return
        endif
      enddo
    else
      m_cnstr_finished = reac_set_finished(ireac)
    endif
  end function m_cnstr_finished

  subroutine m_cnstr_update_reac_coords(ireac)
    use m_routines, only : get_unused_unitnumber,smart_output_open
    integer, intent(in) :: ireac
    integer :: i,j,iidim,ir,icount
    integer :: isum
    integer :: tmpid
    character(256) :: str,suffix
    character(256) :: fname
    real(DP), allocatable, dimension(:) :: tmpval
    integer :: tmpdim
    if(mype==0)then
      suffix=''
      if(m_cnstr_reac_coords_variable()) suffix = '.reac'//trim(adjustl(m_cnstr_get_reac_coord_char()))
      if(m_cnstr_reac_coords_variable().and.(imdalg==T_CONTROL.or.imdalg==VELOCITY_SCALING))then
        ifile_handle = get_unused_unitnumber()
        fname = './nfbluemoon.data'//trim(suffix)
        call smart_output_open(ifile_handle,trim(fname))
        if(icond==INITIAL) write(ifile_handle,'(a)') '#nsteps value1,value2, ... lambda1, lambda2, ...'
      endif
      ifile_handle_structure = get_unused_unitnumber()
      call smart_output_open(ifile_handle_structure,trim(F_DYNM)//trim(suffix))
      call m_IS_wd_speciesname_etc(ifile_handle_structure)

      ifile_handle_energy = get_unused_unitnumber()
      call smart_output_open(ifile_handle_energy,trim(F_ENF)//trim(suffix))
      if(.not.is_md())then
         write(ifile_handle_energy,'(a)') ' iter_ion, iter_total, etotal, forcmx, torqmx, transmx'
      else if (imdalg==VERLET.or.imdalg==T_CONTROL.or.imdalg==VELOCITY_SCALING) then
        write(ifile_handle_energy,'(a)') 'iter_ion, iter, ekina, temperature, etotal, econst'
      endif
    endif
!    if (ireac.eq.1) return
    if( .not. m_cnstr_reac_coords_variable()) return
    if(nvariable.eq.1.and.reac_coord_generation==REAC_COORD_VIA_INPUT)then
      if(ireac .gt. 1) then
      do i=1,n_constrainable_coords
        if(constrainable_coords(i)%is_reaction_coordinate)then
          constrainable_coords(i)%value(:) = constrainable_coords(i)%reaction_coords(ireac-1,:)
          write(str,*) constrainable_coords(i)%ndim
          if(printable)write(nfout,'(a,'//trim(str)//'f20.10)')    &
         & ' new value for the reaction coordinate : '&
         & ,(constrainable_coords(i)%value(iidim),iidim=1,constrainable_coords(i)%ndim)
        endif
      enddo
      endif
    else
      tmpdim=0
      do i=1,nvariable
        tmpdim = tmpdim + constrainable_coords(nvariable_map(i))%ndim
      enddo
      allocate(tmpval(tmpdim))
      write(str,*) tmpdim
      ifile_handle_reac=get_unused_unitnumber()
      open(unit=ifile_handle_reac,file=trim(adjustl(F_REAC)),position='rewind')
      do i=1,ireac-1
        read(unit=ifile_handle_reac,fmt=*,end=10,err=10)
      enddo
      read(unit=ifile_handle_reac,fmt=*,err=10) ir,(tmpval(iidim),iidim=1,tmpdim)
      if(printable) write(nfout,'(a,i8,'//trim(adjustl(str))//'f20.10)') &
     & 'new value for the reaction coordinate:',i,(tmpval(iidim),iidim=1,tmpdim)
      icount=0
      do i=1,nvariable
        tmpid=nvariable_map(i)
        do j=1,constrainable_coords(tmpid)%ndim
          icount = icount+1
          constrainable_coords(tmpid)%value(j) = tmpval(icount)
        enddo
      enddo
      deallocate(tmpval)
      close(ifile_handle_reac)
    endif
    return

10  continue
    if(printable) write(nfout,'(a)') 'encountered error while reading '//F_REAC
    call phase_error_with_msg(nfout,'encountered error while reading '//F_REAC,__LINE__,__FILE__)
  end subroutine m_cnstr_update_reac_coords

  subroutine m_cnstr_pp_reac_coords(stat)
    integer, intent(in) :: stat
    character(len=256) :: str
    integer :: i,j,ii
    real(DP), allocatable, dimension(:) :: tmpval
    integer :: tmpdim,tmpid,icount

    if(mype==0)then
      close(ifile_handle)
      close(ifile_handle_structure)
      close(ifile_handle_energy)
    endif
    if(m_cnstr_reac_coords_variable() .and. .not.is_md().and.mype==0.and.(stat==0.or.stat==1))then
      write(ifile_handle_structure_reac,'(" cps and forc at (iter_ion, iter_total = " &
           & ,i5,i8," )")') curr_md_step, iteration
      call m_Force_wd_force_and_cps(ifile_handle_structure_reac,WITHOUTTAG,cps,natm)
      if(nvariable.le.1)then
        do i=1,n_constrainable_coords
          if(constrainable_coords(i)%is_reaction_coordinate)then
            write(str,*) constrainable_coords(i)%ndim
            write(ifile_handle_energy_reac,'('//trim(str)//'f20.10,i6,4e20.10)')   &
         &  (constrainable_coords(i)%value(ii),ii=1,constrainable_coords(i)%ndim), &
         &   curr_md_step, etotal,etotal*Hartree2eV,etotal*Hartree2eV*eV2kcal_mol,etotal*Hartree2eV*eV2kJ_mol
          endif
        enddo
      else
        tmpdim=0
        icount=0
        do i=1,nvariable
          tmpdim = tmpdim + constrainable_coords(nvariable_map(i))%ndim
        enddo
        allocate(tmpval(tmpdim));tmpval=0.0d0
        do i=1,nvariable
          tmpid=nvariable_map(i)
          do j=1,constrainable_coords(tmpid)%ndim
            icount = icount+1
            tmpval(icount) = constrainable_coords(tmpid)%value(j)
          enddo
        enddo
        write(str,*) tmpdim
        write(ifile_handle_energy_reac,'('//trim(str)//'f20.10,i6,e20.10)') (tmpval(ii),ii=1,tmpdim),curr_md_step,etotal
        deallocate(tmpval)
      endif
    endif

    if(m_cnstr_reac_coords_variable())then
      if(nvariable.eq.1.and.reac_coord_generation==REAC_COORD_VIA_INPUT)then
        do i=1,n_constrainable_coords
          if(constrainable_coords(i)%is_reaction_coordinate)then
            if(stat==0.or.stat==1)then
              constrainable_coords(i)%finished(reac_coord_id)=.true.
            else
              constrainable_coords(i)%finished(reac_coord_id)=.false.
            endif
          endif
        enddo
      else
        if(stat==0.or.stat==1)then
          reac_set_finished(m_cnstr_get_id()) = .true.
        else
          reac_set_finished(m_cnstr_get_id()) = .false.
        endif
      endif
    endif

!    if(mype==0)then
!      call m_Files_reopen_nfcntn()
!      call m_TE_wd_total_energy(nfcntn)
!      call m_CtrlP_wd_isolver(nfcntn)
!      call m_Iter_wd_iteration_numbers(nfcntn)
!      call m_IS_wd_pos_and_v(nfcntn)
!      call constrained_dynamics_dump()
!      close(nfcntn)
!    endif

    if(stat==0.or.stat==1) iteration = 0

  end subroutine m_cnstr_pp_reac_coords

  subroutine m_cnstr_initialize()
    integer :: i,icount,j,k,iidim
    integer :: total_reac_coords
    integer :: ierr
    logical :: force_append
    logical :: exi
    integer :: ireac
    character(len=10) :: suf
    integer :: nreac_max
    character(len=256) :: str
    integer :: tmpid
    real(DP), allocatable, dimension(:) :: tmpval
    integer :: tmpdim,icond_tmp
    icond_tmp = icond
    icond = icond_org

    n_total_reac_coords = 1
    nvariable = 0
    nreac_max = 0
    do i=1,n_constrainable_coords
      if(constrainable_coords(i)%is_reaction_coordinate)then
        n_total_reac_coords = n_total_reac_coords * (constrainable_coords(i)%n_reaction_coords+1)
        nvariable=nvariable+1
        if(constrainable_coords(i)%n_reaction_coords>nreac_max) nreac_max=constrainable_coords(i)%n_reaction_coords
      endif
      allocate(constrainable_coords(i)%origvalue(constrainable_coords(i)%ndim))
      constrainable_coords(i)%origvalue=constrainable_coords(i)%value
    enddo

    if(reac_coord_generation==REAC_COORD_VIA_FILE)then
      inquire(file=trim(adjustl(F_REAC)),exist=exi)
      if(.not.exi)then
        if(printable) write(nfout,'(a)') 'the '//trim(adjustl(F_REAC)) &
      & //' file must exist in order to use reac_coord_generation=via_file'
        call phase_error_with_msg(nfout,'the '//trim(adjustl(F_REAC)) &
        & //' file must exist in order to use reac_coord_generation=via_file' &
        ,__LINE__,__FILE__)
      endif
      ifile_handle_reac=get_unused_unitnumber()
      n_total_reac_coords=0
      open(ifile_handle_reac,file=trim(adjustl(F_REAC)),position='rewind')
      do
        read(ifile_handle_reac,fmt=*,end=10) str
        if(trim(adjustl(str)).eq.'') exit
        n_total_reac_coords=n_total_reac_coords+1
      enddo
10    continue
      if(printable)then
        if(nvariable.eq.1) write(nfout,'(a,i8,a)')       &
       & 'number of reaction coordinate defined in the ' &
       & //trim(adjustl(F_REAC))//' file : ',n_total_reac_coords
        if(nvariable.gt.1) write(nfout,'(a,i8,a)')           &
       & 'number of reaction coordinate set defined in the ' &
       & //trim(adjustl(F_REAC))//' file : ',n_total_reac_coords
      endif
      close(ifile_handle_reac)
    endif

    if(n_total_reac_coords.gt.1)then
      allocate(reac_set_finished(n_total_reac_coords));reac_set_finished=.false.
      allocate(nvariable_map(nvariable))
      icount=1
      do i=1,n_constrainable_coords
        if(constrainable_coords(i)%is_reaction_coordinate)then
          nvariable_map(icount)=i
          constrainable_coords(i)%reac_id=icount
          icount = icount+1
        endif
      enddo
    endif

    if(n_total_reac_coords>1.and..not.is_md().and.mype==0)then
      force_append = convfile_opened
      if(.not.convfile_opened)then
        convfile_opened=.true.
      endif
      suf=''
      if(conf_para.and.mype_conf.gt.0)then
        suf = '_'//trim(adjustl(get_id_char(nrank_conf,mype_conf)))
      endif
      ifile_handle_structure_reac = get_unused_unitnumber()
      call smart_output_open(ifile_handle_structure_reac,trim(F_DYNM)//'.converged'//trim(adjustl(suf)))
      ifile_handle_energy_reac = get_unused_unitnumber()
      call smart_output_open(ifile_handle_energy_reac,trim(F_ENF)//'.converged'//trim(adjustl(suf)))

      call m_IS_wd_speciesname_etc(ifile_handle_structure_reac)
      if(nvariable.le.1)then
        write(ifile_handle_energy_reac,'(a)') &
       & 'reaction coordinate, iter_ion, etotal (hartree), etotal (eV), etotal (kcal/mol), etotal (kJ/mol)'
      else
        write(ifile_handle_energy_reac,'(a)') &
       & 'reaction coordinate 1, reaction_coordinate 2, ... etotal1 (hartree), etotal2 (hartree), ...'
      endif
    endif

    if(nvariable>1.and.icond==INITIAL.and.reac_coord_generation==REAC_COORD_VIA_INPUT)then
      tmpdim=0
      do i=1,nvariable
        tmpdim = tmpdim + constrainable_coords(nvariable_map(i))%ndim
      enddo
      allocate(tmpval(tmpdim))
      write(str,*) tmpdim

      if(mype==0.and.mype_conf==0) then
        ifile_handle_reac=get_unused_unitnumber()
        call smart_output_open(ifile_handle_reac,trim(F_REAC))
      endif

      do i=1,nvariable
        constrainable_coords(nvariable_map(i))%reac_coord_count=1
      enddo

      do i=1,m_cnstr_get_n_reac_coords()
        if(i.eq.1)then
          icount=0
          do j=1,nvariable
            tmpid=nvariable_map(j)
            do k=1,constrainable_coords(tmpid)%ndim
              icount = icount+1
              tmpval(icount)= constrainable_coords(tmpid)%value(k)
            enddo
          enddo
          if(mype==0.and.mype_conf==0)then
            write(ifile_handle_reac,'(i8,'//trim(str)//'f20.10)') i,(tmpval(j),j=1,tmpdim)
          endif
          cycle
        endif
        call incre(constrainable_coords(1))
        icount=0
        tmpval=0.d0
        do j=1,nvariable
          tmpid=nvariable_map(j)
          if(constrainable_coords(tmpid)%reac_coord_count.le.1)then
            do k=1,constrainable_coords(tmpid)%ndim
              icount = icount+1
              tmpval(icount) = constrainable_coords(tmpid)%origvalue(k)
            enddo
          else
            do k=1,constrainable_coords(tmpid)%ndim
              icount = icount+1
              tmpval(icount) = &
            & constrainable_coords(tmpid)%reaction_coords(constrainable_coords(tmpid)%reac_coord_count-1,k)
            enddo
          endif
        enddo
        if(mype==0.and.mype_conf==0)then
          write(ifile_handle_reac,'(i8,'//trim(str)//'f20.10)') i,(tmpval(j),j=1,tmpdim)
        endif
      enddo
      deallocate(tmpval)
    endif
    close(ifile_handle_reac)

    icond = icond_tmp
    reset_md_step = .true.
    if(icond==INITIAL) call m_cnstr_print_cnstr(nfout)
    istart_reac_coord=1
  end subroutine m_cnstr_initialize

  recursive subroutine incre(constrainable_coord)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    integer :: i,j
    character(len=256) :: ch
    if(.not.constrainable_coord%is_reaction_coordinate)return
    if(constrainable_coord%reac_direction==1)then
      if(constrainable_coord%reac_coord_count.eq.constrainable_coord%n_reaction_coords+1)then
        constrainable_coord%reac_coord_count=constrainable_coord%n_reaction_coords+2
        constrainable_coord%reac_direction=constrainable_coord%reac_direction * (-1)
        if(constrainable_coord%reac_id+1.le.nvariable)then
          call incre(constrainable_coords(nvariable_map(constrainable_coord%reac_id+1)))
        endif
      endif
    else
      if(constrainable_coord%reac_coord_count.eq.1)then
        constrainable_coord%reac_coord_count=0
        constrainable_coord%reac_direction=constrainable_coord%reac_direction * (-1)
        if(constrainable_coord%reac_id+1.le.nvariable)then
          call incre(constrainable_coords(nvariable_map(constrainable_coord%reac_id+1)))
        endif
      endif
    endif
    constrainable_coord%reac_coord_count=constrainable_coord%reac_coord_count+constrainable_coord%reac_direction
  end subroutine incre

  subroutine m_cnstr_set_constraint(i,constrainable_coord)
    integer, intent(in) :: i
    type(constrainable_coords_t), intent(in) :: constrainable_coord
    constrainable_coords(n_cnstr_coords_inp+i) = constrainable_coord
  end subroutine m_cnstr_set_constraint

  subroutine m_cnstr_set_Nfree(nf)
    integer, intent(in) :: nf
    Nfree = nf
  end subroutine m_cnstr_set_Nfree

  subroutine m_cnstr_init_cnstr_src()
    integer :: i,j
    do i=n_cnstr_coords_inp+1,n_constrainable_coords
      Nfree = Nfree - constrainable_coords(i)%nfree_decrement
      if(.not.constrainable_coords(i)%mobile) then
        do j=1,constrainable_coords(i)%n_associated_atoms
          is_constrained_atom(constrainable_coords(i)%associated_atoms(j)) = .true.
        enddo
      endif
    enddo
  end subroutine m_cnstr_init_cnstr_src

  subroutine m_cnstr_init_and_read_cnstr(i,constrainable_coord,tag)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    integer, intent(in) :: i
    character(len=*), intent(in) :: tag
    call init_constrainable_coord(constrainable_coord)
    call read_constrainable_coord(i,constrainable_coord,tag)
  end subroutine m_cnstr_init_and_read_cnstr

  ! parse the input
  ! it is assumed that the 'structure' block is being selected, and that the
  ! atomic coordinates have already been read ...
  subroutine m_cnstr_parse_input()
    integer :: f_selectBlock,f_selectParentBlock, &
   &           f_getIntValue,f_getRealValue,f_getStringValue
    integer :: i,iret,j
    character(len=256) :: str,strstr
    character(len=256) :: typ

    ! first we must check the number of constrainable coordinates...
    i=0
    do while(.true.)
      i = i+1
      write(str,*) i
      strstr = trim(adjustl(str))
      if(f_selectBlock(tag_constrainable//strstr)==0)then
        iret = f_selectParentBlock()
      else
        exit
      endif
    enddo

    n_cnstr_coords_inp = i-1
    n_constrainable_coords = n_cnstr_coords_inp+n_cnstr_src
    if (n_constrainable_coords==0) return

    if(f_getStringValue(tag_reac_coord_generation,str,LOWER)==0)then
      if(trim(adjustl(str)).eq.tag_via_file)  reac_coord_generation=REAC_COORD_VIA_FILE
      if(trim(adjustl(str)).eq.tag_via_input) reac_coord_generation=REAC_COORD_VIA_INPUT
    endif

    allocate(constrainable_coords(n_constrainable_coords))
    allocate(force_of_constraintx(natm));force_of_constraintx=0.d0
    allocate(force_of_constrainty(natm));force_of_constrainty=0.d0
    allocate(force_of_constraintz(natm));force_of_constraintz=0.d0
    allocate(is_constrained_atom(natm));is_constrained_atom=.false.

    do i=1,n_cnstr_coords_inp
      call init_constrainable_coord(constrainable_coords(i))
      call read_constrainable_coord(i,constrainable_coords(i),tag_constrainable)
      !!call read_constrainable_coord(i)
    enddo

    do i=1,n_cnstr_coords_inp
      Nfree = Nfree - constrainable_coords(i)%nfree_decrement
      if(.not.constrainable_coords(i)%mobile) then
        do j=1,constrainable_coords(i)%n_associated_atoms
          is_constrained_atom(constrainable_coords(i)%associated_atoms(j)) = .true.
        enddo
      endif
    enddo

    Nreservoir = Nfree+1
!    call m_cnstr_print_cnstr(nfout)
!    contains

!    subroutine read_constrainable_coord(i)
!      integer, intent(in) :: i
!      integer :: f_selectBlock,f_getStringValue,f_selectParentBlock
!      character(len=256)    :: string
!      character(len=256) :: str,strstr
!      character(len=256) :: typ
!      real(DP) :: dret
!      integer :: iret
!      write(str,*) i
!      strstr = trim(adjustl(str))
!      iret = f_selectBlock(tag_constrainable//strstr) ! we already know that this will succeed.
!
!      if(f_getStringValue(tag_type,string,LOWER)==0)then
!        if(adjustl(trim(string)).eq.tag_bond_length) then
!          call read_bl_constraint(constrainable_coords(i))
!        else if(adjustl(trim(string)).eq.tag_distance_from_pos)then
!          call read_dfp_constraint(constrainable_coords(i))
!        else if(adjustl(trim(string)).eq.tag_bond_angle) then
!          call read_bang_constraint(constrainable_coords(i))
!        else if(adjustl(trim(string)).eq.tag_dihedral_angle) then
!          call read_dang_constraint(constrainable_coords(i))
!        else if(adjustl(trim(string)).eq.tag_plane) then
!          call read_plane_constraint(constrainable_coords(i))
!        else if(adjustl(trim(string)).eq.tag_center_of_mass) then
!          call read_com_constraint(constrainable_coords(i))
!        else if(adjustl(trim(string)).eq.tag_coordination_number) then
!          call read_coord_num_constraint(constrainable_coords(i))
!        else if(adjustl(trim(string)).eq.tag_bond_length_diff) then
!          call read_bl_diff_constraint(constrainable_coords(i))
!#ifdef USER_DEFINED_CONSTRAINT
!        else if(adjustl(trim(string)).eq.tag_user_defined) then
!          call read_udef_constraint(constrainable_coords(i))
!#endif
!        else if(adjustl(trim(string)).eq.tag_bond_angle_diff) then
!          call read_bang_diff_constraint(constrainable_coords(i))
!        else
!          write(0,*) "unknown constraint type : "//trim(string)
!          stop
!        endif
!      else
!        stop " 'type' is a necessary field for the specification of constrainable coordinates."
!      endif
!      iret = f_selectParentBlock()
!
!    end subroutine read_constrainable_coord

  end subroutine m_cnstr_parse_input

  subroutine read_constrainable_coord(i,constrainable_coord,tag)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    integer, intent(in) :: i
    character(len=*), intent(in) :: tag
    integer :: f_selectBlock,f_getStringValue,f_selectParentBlock, f_getRealValue
    character(len=256)    :: string
    character(len=256) :: str,strstr
    character(len=256) :: typ
    real(DP) :: dret
    integer :: iret
    write(str,*) i
    strstr = trim(adjustl(str))
    if( f_selectBlock(tag//strstr)==0 ) then! we already know that this will succeed.

    if(f_getStringValue(tag_type,string,LOWER)==0)then
      if(adjustl(trim(string)).eq.tag_bond_length) then
        call read_bl_constraint(constrainable_coord)
      else if(adjustl(trim(string)).eq.tag_bond_angle)then
        call read_bang_constraint(constrainable_coord)
      else if(adjustl(trim(string)).eq.tag_dihedral_angle)then
        call read_dang_constraint(constrainable_coord)
      else if(adjustl(trim(string)).eq.tag_plane) then
        call read_plane_constraint(constrainable_coord)
      else if(adjustl(trim(string)).eq.tag_center_of_mass) then
        call read_com_constraint(constrainable_coord)
      else if(adjustl(trim(string)).eq.tag_coordination_number) then
        call read_coord_num_constraint(constrainable_coord)
      else if(adjustl(trim(string)).eq.tag_bond_length_diff) then
        call read_bl_diff_constraint(constrainable_coord)
      else if(adjustl(trim(string)).eq.tag_bond_length_sum) then
        call read_bl_sum_constraint(constrainable_coord)
#ifdef USER_DEFINED_CONSTRAINT
      else if(adjustl(trim(string)).eq.tag_user_defined) then
        call read_udef_constraint(constrainable_coord)
#endif
      else if(adjustl(trim(string)).eq.tag_bond_angle_diff) then
        call read_bang_diff_constraint(constrainable_coords(i))
      else if(adjustl(trim(string)).eq.tag_distance_from_pos) then
        call read_dfp_constraint(constrainable_coords(i))
      else if(adjustl(trim(string)).eq.tag_distance_from_com) then
        call read_dfc_constraint(constrainable_coords(i))
      else if(adjustl(trim(string)).eq.tag_distance_from_ref) then
        call read_dref_constraint(constrainable_coords(i))
      else
        write(0,*) "unknown 'type' : "//trim(string)
        call phase_error_with_msg(nfout,"unknown 'type' : "//trim(string),__LINE__,__FILE__)
      endif
    else
      call phase_error_with_msg(nfout," 'type' is a necessary field for the specification of constrainable coordinates."&
                               ,__LINE__,__FILE__)
    endif
    call read_rattle_parameters(constrainable_coords(i))
    iret = f_selectParentBlock()

    endif

  end subroutine read_constrainable_coord

  subroutine read_rattle_parameters(constrainable_coord)
    type(constrainable_coords_t),intent(out) :: constrainable_coord
    integer :: f_selectBlock,f_selectParentBlock,f_getRealValue,f_getIntValue
    real(kind=DP) :: dret
    integer :: iret
    if(f_selectBlock("rattle")==0)then
      if(f_getRealValue("tol",dret,'')==0)eps_rattle = dret
      if(f_getRealValue("factor",dret,'')==0) then
        constrainable_coord%rmix_coords = dret
        constrainable_coord%rmix_velocity = dret
      endif
      if(f_getIntValue("max_iteration",iret)==0) max_iteration_rattle = iret
      if(printable .and. ipri>=1) then
        write(nfout,'(a)') ''
        write(nfout,'(a)') ' --- parameters for the RATTLE solver ---'
        write(nfout,'(a,f20.15)') ' tolerance     : ',eps_rattle
        write(nfout,'(a,f20.5)')  ' factor        : ',constrainable_coord%rmix_coords
        write(nfout,'(a,i20)')    ' max_iteration : ',max_iteration_rattle
      endif
      iret = f_selectParentBlock()
    endif
  end subroutine read_rattle_parameters

  ! returns true if any constraints are specified in the input file.
  logical function m_cnstr_constraints_exist()
    integer :: i
    m_cnstr_constraints_exist = allocated(constrainable_coords)
    if(m_cnstr_constraints_exist)then
      do i=1,n_constrainable_coords
        if(.not.constrainable_coords(i)%mobile)then
          return
        endif
      enddo
    endif
    m_cnstr_constraints_exist = .false.
  end function m_cnstr_constraints_exist

  ! print status of the defined constrainable coordinates to the specified file handler
  subroutine m_cnstr_print_cnstr(fil)
    integer, intent(in) :: fil
    real(DP), allocatable, dimension(:) :: tmpval
    integer :: i,j,ii,tmpdim,ir
    character(len=256) str
    if (.not.printable) return
    write(nfout,'(a)') ''
    if (.not.m_cnstr_constraints_exist()) then
      write(fil,'(a)') '  constrainable coordinates do not exist'
      return
    endif
    write(fil,'(a)')    '--- status for the constrainable coordinates ---'
    write(fil,'(a,i5)') '  number of freedom : ',Nfree
    do i=1,n_constrainable_coords
      write(fil,'(a,i5)') '  constrainable coordinate no. ',i
      write(fil,'(a,a)')  '    name        : ', trim(constrainable_coords(i)%nam)
      write(fil,'(a,a)')  '    description : ', trim(constrainable_coords(i)%descri)
      write(fil,'(a,1f20.10,a)') '    value  : ',constrainable_coords(i)%value(1),&
     & ' '//trim(constrainable_coords(i)%unit_name)
      write(fil,'(a,l5)') '    mobile  : ',constrainable_coords(i)%mobile
      write(fil,'(a,l5)') '    monitor : ',constrainable_coords(i)%monitor
      write(str,*) constrainable_coords(i)%n_associated_atoms
      write(fil,'(a,'//trim(adjustl(str))//'i5)')  '    associated atoms : ',&
     &  (constrainable_coords(i)%associated_atoms(j), &
     &   j=1, constrainable_coords(i)%n_associated_atoms)
      write(fil,'(a,l5)')  '    reaction coordinate : ',  &
     &  constrainable_coords(i)%is_reaction_coordinate
      if (constrainable_coords(i)%is_reaction_coordinate.and.nvariable.eq.1) then
        write(str,*) constrainable_coords(i)%ndim
        write(fil,'(a)')  '    the following values will be taken into account by the current rank.'
        do j=1,constrainable_coords(i)%n_reaction_coords+1
          if(mod(j-1,nrank_conf)/=mype_conf) cycle
          if(j==1)then
            write(fil,'('//trim(str)//'f20.10,a,l5)') (constrainable_coords(i)%value(ii),&
         &  ii=1,constrainable_coords(i)%ndim),&
         &  ' '//trim(constrainable_coords(i)%unit_name) &
         & , constrainable_coords(i)%finished(1)
          else
            write(fil,'('//trim(str)//'f20.10,a,l5)') (constrainable_coords(i)%reaction_coords(j-1,ii),&
         &   ii=1,constrainable_coords(i)%ndim),&
         &   ' '//trim(constrainable_coords(i)%unit_name) &
         &   , constrainable_coords(i)%finished(j)
          endif
        enddo
      endif
    enddo

    if (n_total_reac_coords.gt.1) then
      write(nfout,'(a,i5)') '  total number of reaction coordinates : ',n_total_reac_coords
      if(reac_coord_generation==REAC_COORD_VIA_FILE)then
        write(fil,'(a)')    '  reaction coordinate generation scheme : '//tag_via_file
      else
        write(fil,'(a)')    '  reaction coordinate generation scheme : '//tag_via_input
      endif
    endif

    if(nvariable.gt.1.or.reac_coord_generation==REAC_COORD_VIA_FILE)then
      write(fil,'(a,i8)') '  number of variable reaction coordinates : ',nvariable
      write(fil,'(a)')    '  reaction coordinates taken into account by the current rank:'
      tmpdim=0
      do i=1,nvariable
        tmpdim = tmpdim + constrainable_coords(nvariable_map(i))%ndim
      enddo
      allocate(tmpval(tmpdim))
      write(str,*) tmpdim
      ifile_handle_reac=get_unused_unitnumber()
      open(unit=ifile_handle_reac,file=trim(adjustl(F_REAC)),position='rewind')
      do i=1,n_total_reac_coords
        read(unit=ifile_handle_reac,fmt=*,end=10,err=10) ir,(tmpval(ii),ii=1,tmpdim)
        if(mod(i-1,nrank_conf)/=mype_conf) cycle
        write(fil,'(a,i8,'//trim(adjustl(str))//'f20.10,l5)') &
        & '    ',ir,(tmpval(ii),ii=1,tmpdim),reac_set_finished(i)
      enddo
10    continue
      deallocate(tmpval)
      close(ifile_handle_reac)
    endif
    write(fil,'(a)')
  end subroutine m_cnstr_print_cnstr

  subroutine update_sigma(constrainable_coord,iidim)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    integer, intent(in) :: iidim
    if(constrainable_coord%typ==BOND_LENGTH)then
      call update_sigma_bl(constrainable_coord)
    else if (constrainable_coord%typ==DISTANCE_FROM_POS)then
      call update_sigma_dfp(constrainable_coord)
    else if (constrainable_coord%typ==DISTANCE_FROM_COM)then
      call update_sigma_dfc(constrainable_coord)
    else if (constrainable_coord%typ==BOND_ANGLE)then
      call update_sigma_bang(constrainable_coord)
    else if (constrainable_coord%typ==DIHEDRAL_ANGLE)then
      call update_sigma_dang(constrainable_coord)
    else if (constrainable_coord%typ==PLANE)then
      call update_sigma_plane(constrainable_coord,iidim)
    else if (constrainable_coord%typ==BOND_LENGTH_DIFF)then
      call update_sigma_bl_diff(constrainable_coord)
    else if (constrainable_coord%typ==BOND_LENGTH_SUM)then
      call update_sigma_bl_sum(constrainable_coord)
    else if (constrainable_coord%typ==CENTER_OF_MASS)then
      call update_sigma_com(constrainable_coord,iidim)
    else if (constrainable_coord%typ==COORDNATION_NUMBER)then
      call update_sigma_coord_num(constrainable_coord)
#ifdef USER_DEFINED_CONSTRAINT
    else if (constrainable_coord%typ==USER_DEFINED)then
      call update_sigma_udef_constraint(constrainable_coord,iidim)
#endif
    else if (constrainable_coord%typ==BOND_ANGLE_DIFF)then
      call update_sigma_bang_diff(constrainable_coord)
    else if (constrainable_coord%typ==DISTANCE_FROM_REF)then
      call update_sigma_dref(constrainable_coord)
    else
      if(printable)write(nfout,'(a)') 'unimplemented constraint type : '//constrainable_coord%nam
    endif
  end subroutine update_sigma

  subroutine update_dsigma(constrainable_coord,iidim)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    integer, intent(in) :: iidim
    integer :: i,j
    real(DP), allocatable, dimension(:,:) :: forc_l_tmp
    real(DP), pointer, dimension(:,:) :: forcw
    integer, pointer, dimension(:) :: ipfrc
    if(constrainable_coord%typ==BOND_LENGTH)then
      call update_dsigma_bl(constrainable_coord)
    else if (constrainable_coord%typ==DISTANCE_FROM_POS)then
      call update_dsigma_dfp(constrainable_coord)
    else if (constrainable_coord%typ==DISTANCE_FROM_COM)then
      call update_dsigma_dfc(constrainable_coord)
    else if (constrainable_coord%typ==BOND_ANGLE)then
      call update_dsigma_bang(constrainable_coord)
    else if (constrainable_coord%typ==DIHEDRAL_ANGLE)then
      call update_dsigma_dang(constrainable_coord)
    else if (constrainable_coord%typ==PLANE)then
      call update_dsigma_plane(constrainable_coord,iidim)
    else if (constrainable_coord%typ==BOND_LENGTH_DIFF)then
      call update_dsigma_bl_diff(constrainable_coord)
    else if (constrainable_coord%typ==BOND_LENGTH_SUM)then
      call update_dsigma_bl_sum(constrainable_coord)
    else if (constrainable_coord%typ==CENTER_OF_MASS)then
      call update_dsigma_com(constrainable_coord,iidim)
    else if (constrainable_coord%typ==COORDNATION_NUMBER)then
      call update_dsigma_coord_num(constrainable_coord)
#ifdef USER_DEFINED_CONSTRAINT
    else if (constrainable_coord%typ==USER_DEFINED)then
      call update_dsigma_udef_constraint(constrainable_coord,iidim)
#endif
    else if (constrainable_coord%typ==BOND_ANGLE_DIFF)then
      call update_dsigma_bang_diff(constrainable_coord)
    else if (constrainable_coord%typ==DISTANCE_FROM_REF)then
      call update_dsigma_dref(constrainable_coord)
    else
      if(printable)write(nfout,'(a)') 'unimplemented constraint type : '//constrainable_coord%nam
      return
    endif

    if(nopr>=2)then
      do i=1,constrainable_coord%ndim
        allocate(forc_l_tmp(natm,3));forc_l_tmp=0.d0
        allocate(forcw(natm,3));allocate(ipfrc(natm2))
        do j=1,constrainable_coord%n_associated_atoms
          forc_l_tmp(constrainable_coord%associated_atoms(j),:) = constrainable_coord%dsigma(i,j,:)
        enddo
        if(nopr_supercell <= nopr) then
           call fd_symmetrize(natm2,natm,natm,napt,nopr+af,nopr,op,iwei &
                &, forc_l_tmp, forcw, ipfrc) ! -(bottom_Subroutines_para)
        else
           Call fd_supercell_symmetrize(natm2,natm,natm,napt_supercell,nopr_supercell,op,nopr+af,iwei &
                & , forc_l_tmp, forcw, ipfrc, iop_supercell)
        end if
        do j=1,constrainable_coord%n_associated_atoms
          constrainable_coord%dsigma(i,j,:) = forc_l_tmp(constrainable_coord%associated_atoms(j),:)
          constrainable_coord%dsigma_old(i,j,:) = forc_l_tmp(constrainable_coord%associated_atoms(j),:)
        enddo
        deallocate(forc_l_tmp);deallocate(forcw);deallocate(ipfrc)
      enddo
    endif

  end subroutine update_dsigma

  subroutine monitor_constrainable_coord(constrainable_coord)
    type(constrainable_coords_t), intent(in) :: constrainable_coord
    if(constrainable_coord%typ==BOND_LENGTH)then
      call monitor_bl(constrainable_coord)
    else if (constrainable_coord%typ==DISTANCE_FROM_POS)then
      call monitor_dfp(constrainable_coord)
    else if (constrainable_coord%typ==DISTANCE_FROM_COM)then
      call monitor_dfc(constrainable_coord)
    else if (constrainable_coord%typ==BOND_ANGLE)then
      call monitor_bang(constrainable_coord)
    else if (constrainable_coord%typ==DIHEDRAL_ANGLE)then
      call monitor_dang(constrainable_coord)
    else if (constrainable_coord%typ==BOND_LENGTH_DIFF)then
      call monitor_bl_diff(constrainable_coord)
    else if (constrainable_coord%typ==BOND_LENGTH_SUM)then
      call monitor_bl_sum(constrainable_coord)
    else if (constrainable_coord%typ==PLANE)then
      call monitor_plane(constrainable_coord)
    else if (constrainable_coord%typ==CENTER_OF_MASS)then
      call monitor_com(constrainable_coord)
    else if (constrainable_coord%typ==COORDNATION_NUMBER)then
      call monitor_coord_num(constrainable_coord)
#ifdef USER_DEFINED_CONSTRAINT
    else if (constrainable_coord%typ==USER_DEFINED)then
      call monitor_udef_constraint(constrainable_coord)
#endif
    else if (constrainable_coord%typ==BOND_ANGLE_DIFF)then
      call monitor_bang_diff(constrainable_coord)
    else if (constrainable_coord%typ==DISTANCE_FROM_REF)then
      call monitor_dref(constrainable_coord)
    else
      if(printable)write(nfout,'(a)') 'monitor unimplemented for constraint type : '//trim(constrainable_coord%nam)
    endif
  end subroutine monitor_constrainable_coord

  subroutine m_cnstr_update_dsigma_old()
    integer :: i
    do i=1,n_constrainable_coords
      constrainable_coords(i)%dsigma_old(:,:,:) = 0.0d0
      constrainable_coords(i)%dsigma_old(:,:,:) = constrainable_coords(i)%dsigma(:,:,:)
    enddo
  end subroutine m_cnstr_update_dsigma_old

  subroutine print_lambda()
    integer :: i,iidim
    character(256) :: str
    do i=1,n_constrainable_coords
      if(constrainable_coords(i)%mobile)cycle
      write(str,*) constrainable_coords(i)%ndim
      if(printable) &
     & write(nfout,'(a,i5,a,'//trim(adjustl(str))//'f20.10)')    &
     & '    lambda for constraint no. ',i,':',                   &
     & (constrainable_coords(i)%lambda(iidim),iidim=1,constrainable_coords(i)%ndim)

      if(constrainable_coords(i)%is_reaction_coordinate.and.mype==0.and. &
       &(imdalg==T_CONTROL.or.imdalg==VELOCITY_SCALING))then
        write(str,*) 2*constrainable_coords(i)%ndim+1
        write(ifile_handle,'(i8,'//trim(str)//'f20.10)') curr_md_step, &
     & (constrainable_coords(i)%value(iidim), iidim=1,constrainable_coords(i)%ndim), &
     & (constrainable_coords(i)%lambda(iidim),iidim=1,constrainable_coords(i)%ndim), &
     &  constrainable_coords(i)%det_metric
      endif

    enddo
  end subroutine print_lambda

  subroutine m_cnstr_print_status()
    integer :: i
    call print_lambda()
    do i=1,n_constrainable_coords
      if(constrainable_coords(i)%monitor.and.printable)then
        call monitor_constrainable_coord(constrainable_coords(i))
      endif
    enddo
  end subroutine m_cnstr_print_status

  subroutine m_cnstr_get_cnstr_coords(xx)
    real(DP), dimension(:,:), intent(inout) :: xx
    integer  :: i,j,k,l,i1,iidim
    real(DP) :: dt2,dt2inv
    real(DP) :: factor
    real(DP) :: denom
    real(DP) :: currlambda
    real(DP) :: large_force=100.d0
    logical  :: all_constraints_converged
    dt2 = 0.5d0*dtio*dtio
    do i=1,n_constrainable_coords
      constrainable_coords(i)%lambda = 0.d0
    enddo

    if(printable.and.ipri>=2)write(nfout,'(a)') 'fixing the coordinates according to the RATTLE algorithm'
    dt2inv = 1.d0/dtio/dtio
    do i=1,max_iteration_rattle
      cloop:do j=1,n_constrainable_coords
        if(constrainable_coords(j)%mobile) cycle
        do iidim=1,constrainable_coords(j)%ndim
          call update_sigma(constrainable_coords(j),iidim)
          call update_dsigma(constrainable_coords(j),iidim)
          denom = 0.d0
          do k=1,constrainable_coords(j)%n_associated_atoms
            i1 = constrainable_coords(j)%associated_atoms(k)
            if(i1.eq.0)cycle
            factor = 0.5d0/(amion(ityp(i1)))
            denom  = denom + factor * dot_product(         &
          & constrainable_coords(j)%dsigma_old(iidim,k,:), &
          & constrainable_coords(j)%dsigma(iidim,k,:))
          enddo

          if(dabs(denom).lt.very_small)then
            write(0,*) 'at get_constrained_coords : denom is too small'
            write(0,*) 'iteration no.',i
            call phase_error_with_msg(nfout,'at get_constrained_coords : denom is too small',__LINE__,__FILE__)
!             exit cloop
          endif
          currlambda = constrainable_coords(j)%rmix_coords * dt2inv * constrainable_coords(j)%sigma(iidim)/denom
          constrainable_coords(j)%lambda(iidim) = constrainable_coords(j)%lambda(iidim)+currlambda
          do k=1,constrainable_coords(j)%n_associated_atoms
            i1 = constrainable_coords(j)%associated_atoms(k)
            if (.not.mobile(i1)) cycle
            factor = dt2/(amion(ityp(i1)))
            xx(i1,1) = xx(i1,1) - factor * currlambda * constrainable_coords(j)%dsigma_old(iidim,k,1)
            xx(i1,2) = xx(i1,2) - factor * currlambda * constrainable_coords(j)%dsigma_old(iidim,k,2)
            xx(i1,3) = xx(i1,3) - factor * currlambda * constrainable_coords(j)%dsigma_old(iidim,k,3)
            x_nopbc(i1) = x_nopbc(i1) - factor * currlambda * constrainable_coords(j)%dsigma_old(iidim,k,1)
            y_nopbc(i1) = y_nopbc(i1) - factor * currlambda * constrainable_coords(j)%dsigma_old(iidim,k,2)
            z_nopbc(i1) = z_nopbc(i1) - factor * currlambda * constrainable_coords(j)%dsigma_old(iidim,k,3)
          enddo
        enddo
      enddo cloop

      all_constraints_converged = .true.
      do j=1,n_constrainable_coords
        if(constrainable_coords(j)%mobile) cycle
        do iidim=1,constrainable_coords(j)%ndim
          !!$if(printable) write(0,*) i,constrainable_coords(j)%sigma(iidim)
          if(dabs(constrainable_coords(j)%sigma(iidim)).gt.eps_rattle) then
            all_constraints_converged = .false.
            exit
          endif
        enddo
      enddo
      if (all_constraints_converged) exit

      if(i.eq.max_iteration_rattle)then
        write(0,*) 'RATTLE iteration unconverged (coordinates)'
        do j=1,n_constrainable_coords
          write(0,*) 'lambda1 : ',constrainable_coords(j)%lambda(1)
        enddo
        call phase_error_with_msg(nfout,'RATTLE iteration unconverged (coordinates)',__LINE__,__FILE__)
      endif

    enddo
!    if(printable.and.ipri>=2)write(nfout,'(a,i5)') 'done RATTLE (coords). iterations needed for convergence: ',i
    if(printable.and.iprimd>=2)write(nfout,'(a,i5)') 'done RATTLE (coords). iterations needed for convergence: ',i
    do i=1,n_constrainable_coords
      call update_det_metric(constrainable_coords(i))
    enddo
  end subroutine m_cnstr_get_cnstr_coords

  subroutine m_cnstr_udate_for_cnstr()
    integer :: i,j,i1,iidim
    force_of_constraintx=0.d0;force_of_constrainty=0.d0;force_of_constraintz=0.0d0
    ! store the 'force of constraint'
    do i=1,n_constrainable_coords
      if(constrainable_coords(i)%mobile)cycle
      do iidim=1,constrainable_coords(i)%ndim
        do j=1,constrainable_coords(i)%n_associated_atoms
          i1 = constrainable_coords(i)%associated_atoms(j)
          if (.not.mobile(i1)) cycle
          force_of_constraintx(i1) = force_of_constraintx(i1) - &
        & constrainable_coords(i)%lambda(iidim) * constrainable_coords(i)%dsigma(iidim,j,1)
          force_of_constrainty(i1) = force_of_constrainty(i1) - &
        & constrainable_coords(i)%lambda(iidim) * constrainable_coords(i)%dsigma(iidim,j,2)
          force_of_constraintz(i1) = force_of_constraintz(i1) - &
        & constrainable_coords(i)%lambda(iidim) * constrainable_coords(i)%dsigma(iidim,j,3)
        enddo
      enddo
    enddo
  end subroutine m_cnstr_udate_for_cnstr

  subroutine m_cnstr_get_cnstr_vel(veloc)
    real(DP), dimension(:,:), intent(inout) :: veloc
    integer :: i,j,k,l,i1,iidim
    real(DP) :: denom,numera,data1,denom1,numera1,eta1
    real(DP), dimension(3) :: tmpv
    real(DP) :: tmpeta,tmpsum,massinv
    real(DP) :: largest_value
    logical :: all_constraints_converged
    if(printable.and.ipri>=2)write(nfout,'(a)') 'fixing the velocities according to the RATTLE algorithm'
    do i=1,max_iteration_rattle
      cloop:do j=1,n_constrainable_coords
        if(constrainable_coords(j)%mobile) cycle
        do iidim=1,constrainable_coords(j)%ndim
          !calculate eta
          numera = 0.d0
          denom  = 0.d0
          do k=1,constrainable_coords(j)%n_associated_atoms
            i1 = constrainable_coords(j)%associated_atoms(k)
            numera = numera + veloc(i1,1)*constrainable_coords(j)%dsigma(iidim,k,1) &
          &                 + veloc(i1,2)*constrainable_coords(j)%dsigma(iidim,k,2) &
          &                 + veloc(i1,3)*constrainable_coords(j)%dsigma(iidim,k,3)
            denom  = denom  + dot_product(             &
          & constrainable_coords(j)%dsigma(iidim,k,:), &
          & constrainable_coords(j)%dsigma(iidim,k,:))/amion(ityp(i1))
          enddo
          if(dabs(denom).lt.very_small)then
            write(0,*) 'at get_constrained_velocities : denom is too small'
            call phase_error_with_msg(nfout,'at get_constrained_velocities : denom is too small',__LINE__,__FILE__)
!             exit cloop
          endif
          tmpeta = constrainable_coords(j)%rmix_velocity*numera/denom
          do k=1,constrainable_coords(j)%n_associated_atoms
            i1 = constrainable_coords(j)%associated_atoms(k)
            if(.not.mobile(i1))cycle
            massinv = 1.0d0/amion(ityp(i1))
            veloc(i1,1) = veloc(i1,1) - massinv*tmpeta*constrainable_coords(j)%dsigma(iidim,k,1)
            veloc(i1,2) = veloc(i1,2) - massinv*tmpeta*constrainable_coords(j)%dsigma(iidim,k,2)
            veloc(i1,3) = veloc(i1,3) - massinv*tmpeta*constrainable_coords(j)%dsigma(iidim,k,3)
          enddo
        enddo
      enddo cloop
      all_constraints_converged = .true.
      largest_value=0.d0
      do j=1,n_constrainable_coords
        if (constrainable_coords(j)%mobile) cycle
        do iidim=1,constrainable_coords(j)%ndim
          tmpsum = 0.d0
          do k=1,constrainable_coords(j)%n_associated_atoms
            i1 = constrainable_coords(j)%associated_atoms(k)
            tmpsum = tmpsum + (veloc(i1,1)*constrainable_coords(j)%dsigma(iidim,k,1) &
          &                 +  veloc(i1,2)*constrainable_coords(j)%dsigma(iidim,k,2) &
          &                 +  veloc(i1,3)*constrainable_coords(j)%dsigma(iidim,k,3))
          enddo
          if(dabs(tmpsum).gt.largest_value)largest_value=tmpsum
          if(dabs(tmpsum).gt.eps_rattle)then
            all_constraints_converged = .false.
            exit
          endif
        enddo
      enddo
      if(i.eq.1)then
        data1=largest_value
        denom1=denom
        numera1=numera
        eta1=tmpeta
      endif
      if (all_constraints_converged) exit
      if (i.eq.max_iteration_rattle) then
        write(0,*) 'RATTLE iteration unconverged (velocities)'
        call phase_error_with_msg(nfout,'RATTLE iteration unconverged (velocities)',__LINE__,__FILE__)
      endif
    enddo
!    if(printable.and.ipri>=2) write(nfout,'(a,i5)') 'done RATTLE (velocities). iteration needed for convergence: ',i
    if(printable.and.iprimd>=2) write(nfout,'(a,i5)') 'done RATTLE (velocities). iteration needed for convergence: ',i
  end subroutine m_cnstr_get_cnstr_vel

  ! deallocate arrrays, and do any finalizations deemed necessary
  subroutine m_cnstr_finalize()
    integer :: i
    if(m_cnstr_constraints_exist()) then
      if(m_cnstr_reac_coords_variable().and..not.is_md().and.mype==0)then
        close(ifile_handle_structure_reac)
        close(ifile_handle_energy_reac)
      endif
      deallocate(force_of_constraintx)
      deallocate(force_of_constrainty)
      deallocate(force_of_constraintz)
      do i=1,n_constrainable_coords
        call dealloc_constrainable_coord(constrainable_coords(i))
      enddo
      deallocate(constrainable_coords)
      deallocate(is_constrained_atom)
      if(allocated(nvariable_map)) deallocate(nvariable_map)
      if(allocated(reac_set_finished)) deallocate(reac_set_finished)
    endif
  end subroutine m_cnstr_finalize

  subroutine m_cnstr_first_step()
    integer :: i,iidim
    do i=1,n_constrainable_coords
      do iidim=1,constrainable_coords(i)%ndim
        call update_dsigma(constrainable_coords(i),iidim)
      enddo
    enddo
    call m_cnstr_update_dsigma_old()
  end subroutine m_cnstr_first_step

  subroutine m_cnstr_wd_constraints(nfcntn)
    integer, intent(in) :: nfcntn
    integer :: i0,i,j,ii
    if(.not.m_cnstr_constraints_exist()) return
    if(mype/=0)return
    write(nfcntn,'(a)') tag_constrainable_coords
    write(nfcntn,'(a)') '(n_constrainable_coords)'
    write(nfcntn,'(i5)') n_constrainable_coords
    if(m_cnstr_reac_coords_variable())then
      write(nfcntn,'(a)') '(reac_coord_id)'
      write(nfcntn,'(i5)') reac_coord_id
    endif
    write(nfcntn,'(a)') '(force of constraint x,y,z)'
    do i=1,natm
      write(nfcntn,'(3e20.10)') force_of_constraintx(i),force_of_constrainty(i),force_of_constraintz(i)
    enddo
    do i0=1,n_constrainable_coords
      write(nfcntn,'(a,i5,a)') '(constrainable coordinate no.',i0,')'
      call m_cnstr_dump(constrainable_coords(i0),nfcntn)
    enddo

    if(nvariable.gt.1.or.reac_coord_generation==REAC_COORD_VIA_FILE)then
      write(nfcntn,'(a)') '(reaction coordinate set taken into account)'
      do i=1,n_total_reac_coords
        write(nfcntn,'(i8,l5)') i,reac_set_finished(i)
      enddo
    endif
  end subroutine m_cnstr_wd_constraints

  subroutine m_cnstr_dump(constrainable_coord,nfcntn)
    type(constrainable_coords_t), intent(in) :: constrainable_coord
    integer, intent(in)                      :: nfcntn
    integer :: ii,i,j
    character(len=256) :: str
    write(nfcntn,'(a)') '(name)'
    write(nfcntn,'(a)') trim(constrainable_coord%nam)
    write(nfcntn,'(a)') '(description)'
    write(nfcntn,'(a)') trim(constrainable_coord%descri)
    write(nfcntn,'(a)') '(unit)'
    write(nfcntn,'(a)') trim(constrainable_coord%unit_name)
    write(nfcntn,'(a)') '(type)'
    write(nfcntn,'(i2)') constrainable_coord%typ
    write(nfcntn,'(a)') '(number of associated atoms)'
    write(nfcntn,'(i5)') constrainable_coord%n_associated_atoms
    write(nfcntn,'(a)') '(associated atoms)'
    write(str,*) constrainable_coord%n_associated_atoms
    write(nfcntn,'('//trim(adjustl(str))//'i5)') (constrainable_coord%associated_atoms(ii), &
   & ii=1,constrainable_coord%n_associated_atoms)
    write(nfcntn,'(a)') '(dimension)'
    write(nfcntn,'(i5)') constrainable_coord%ndim
    write(nfcntn,'(a)') '(value)'
    write(str,*) constrainable_coord%ndim
    write(nfcntn,'('//trim(adjustl(str))//'f20.10)') (constrainable_coord%value(ii), &
   & ii=1,constrainable_coord%ndim)
    write(nfcntn,'(a)') '(sigma)'
    write(nfcntn,'('//trim(adjustl(str))//'f20.10)') (constrainable_coord%sigma(ii), &
   & ii=1,constrainable_coord%ndim)
    write(nfcntn,'(a)') '(dsigma)'
    do i=1,constrainable_coord%ndim
      write(nfcntn,'(a,i5,a)') '  (for dimension ',i,' )'
      do j=1,constrainable_coord%n_associated_atoms
        write(nfcntn,'(a,i5,a)') '    (for atom ',j,' )'
        write(nfcntn,'(3e20.10)') (constrainable_coord%dsigma(i,j,ii),ii=1,3)
      enddo
    enddo
    write(nfcntn,'(a)') '(dsigma_old)'
    do i=1,constrainable_coord%ndim
      write(nfcntn,'(a,i5,a)') '  (for dimension ',i,' )'
      do j=1,constrainable_coord%n_associated_atoms
        write(nfcntn,'(a,i5,a)') '    (for atom ',j,' )'
        write(nfcntn,'(3e20.10)') (constrainable_coord%dsigma_old(i,j,ii),ii=1,3)
      enddo
    enddo
    write(nfcntn,'(a)') '(mobile?)'
    write(nfcntn,'(l5)') constrainable_coord%mobile
    write(nfcntn,'(a)') '(monitor?)'
    write(nfcntn,'(l5)') constrainable_coord%monitor
    write(nfcntn,'(a)') '(determinant of the metric tensor)'
    write(nfcntn,'(1e20.10)') constrainable_coord%det_metric
    write(nfcntn,'(a)') '(lambda)'
    write(str,*) constrainable_coord%ndim
    write(nfcntn,'('//trim(adjustl(str))//'e20.10)') (constrainable_coord%lambda(ii), &
   & ii=1,constrainable_coord%ndim)
    write(nfcntn,'(a)') '(reaction coordinate?)'
    write(nfcntn,'(l5)') constrainable_coord%is_reaction_coordinate
    if(constrainable_coord%is_reaction_coordinate)then
      write(nfcntn,'(a)') '(number of reaction coordinates)'
      write(nfcntn,'(i5)') constrainable_coord%n_reaction_coords
      write(nfcntn,'(a)') '(reaction coordinates taken into account)'
      write(str,*) constrainable_coord%ndim
      do i=1,constrainable_coord%n_reaction_coords
        write(nfcntn,'('//trim(adjustl(str))//'e20.10)') &
       & (constrainable_coord%reaction_coords(i,ii),ii=1,constrainable_coord%ndim)
      enddo
      if(conf_para.and.constrainable_coord%is_reaction_coordinate.and. &
      & nvariable.eq.1.and.reac_coord_generation==REAC_COORD_VIA_INPUT) then
        write(nfcntn,'(a)') '(status for each reaction coordinate)'
        do i=1,m_cnstr_get_n_reac_coords()
          write(nfcntn,'(l5)') constrainable_coord%finished(i)
        enddo
      endif
    endif
  end subroutine m_cnstr_dump

  subroutine m_cnstr_rd_constraints(nfcntn)
    integer, intent(in)    :: nfcntn
    integer, parameter     :: len_str = 132
    character(len=len_str) :: str
    integer :: i0,i,j,ierr,itmp,ii,ir
    logical :: tag_is_found, EOF_reach
    if(mype==0)then
      call rewind_to_tag0(nfcntn,len(tag_constrainable_coords),tag_constrainable_coords &
     & , EOF_reach, tag_is_found, str, len_str)
    endif

    call mpi_bcast(tag_is_found,1,mpi_logical,0,MPI_CommGroup,ierr)

    if(.not.tag_is_found) then
      if(printable)write(nfout,'(a)') 'tag: '//tag_constrainable_coords//' not found in the nfcntn file.'
      return
    else
      if(printable) write(nfout,'(a)') 'tag: '//tag_constrainable_coords//&
     & ' found; loading constraint(s) from the nfcntn file.'
    endif

    if(mype==0)then
      read(nfcntn,*)
      read(nfcntn,*)
      if(m_cnstr_reac_coords_variable())then
        read(nfcntn,*)
        read(nfcntn,*) itmp
        write(nfout,'(a,i5)') 'last run terminated at reaction coordinate no. ',itmp
        istart_reac_coord = itmp
        if(.not.conf_para) reac_coord_id = itmp-1
      else
        if(.not.conf_para) reac_coord_id = 0
      endif

      read(nfcntn,*)
      do i=1,natm
        read(nfcntn,*) force_of_constraintx(i),force_of_constrainty(i),force_of_constraintz(i)
      enddo
      write(nfout,'(a)') 'read force of constraint'
      do i=1,natm
        write(nfout,'(3f20.10)') force_of_constraintx(i),force_of_constrainty(i),force_of_constraintz(i)
      enddo

    endif

    if(npes>1)then
      if(.not.conf_para) call mpi_bcast(reac_coord_id,1,mpi_integer,0,MPI_CommGroup,ierr)
      call mpi_bcast(istart_reac_coord,1,mpi_integer,0,MPI_CommGroup,ierr)
      call mpi_bcast(force_of_constraintx,natm,mpi_double_precision,0,MPI_CommGroup,ierr)
      call mpi_bcast(force_of_constrainty,natm,mpi_double_precision,0,MPI_CommGroup,ierr)
      call mpi_bcast(force_of_constraintz,natm,mpi_double_precision,0,MPI_CommGroup,ierr)
    endif

    do i0=1,n_constrainable_coords
      if(printable) write(nfout,'(a,i5)') 'loading constraint no. ',i0
      call m_cnstr_load(constrainable_coords(i0),nfcntn)
    enddo

    if((nvariable.gt.1 .or. reac_coord_generation==REAC_COORD_VIA_FILE))then
      if(mype==0)then
        read(nfcntn,*)
        do i=1,n_total_reac_coords
          read(nfcntn,*) ir,reac_set_finished(i)
        enddo
      endif
      if(npes>1)then
        call mpi_bcast(reac_set_finished,n_total_reac_coords,mpi_logical,0,MPI_CommGroup,ierr)
      endif
      if(printable)then
        do i=1,n_total_reac_coords
          write(nfout,'(a,i8,a,l5)') ' reaction coordinate set no. ',i,' has already been taken into account:',&
         & reac_set_finished(i)
        enddo
      endif
    endif

    if(printable) write(nfout,'(a)') 'successfully read constraints from file.'

!!$    call m_cnstr_print_cnstr(nfout)
  end subroutine m_cnstr_rd_constraints

  subroutine m_cnstr_load(constrainable_coord,nfcntn)
    type(constrainable_coords_t), intent(inout) :: constrainable_coord
    integer, intent(in)                      :: nfcntn
    real(DP), allocatable, dimension(:) :: tmpval
    real(DP), allocatable, dimension(:,:) :: tmpval2
    integer :: ntmp,ii,i,j,ierr
    integer, parameter     :: len_str = 132
    character(len=len_str) :: str
    if(mype==0)then
      read(nfcntn,*)
      read(nfcntn,*)
      read(nfcntn,*)
      read(nfcntn,*)
      read(nfcntn,*)
      read(nfcntn,*)
      read(nfcntn,*)
      read(nfcntn,*)
      read(nfcntn,*)
      read(nfcntn,*)
      read(nfcntn,*)
      read(nfcntn,*)
      read(nfcntn,*)
      read(nfcntn,*)
      read(nfcntn,*)
      read(nfcntn,*)
      allocate(tmpval(constrainable_coord%ndim))
      read(nfcntn,*) (tmpval(ii),ii=1,constrainable_coord%ndim)
      write(str,*) constrainable_coord%ndim
      write(nfout,'(a,'//trim(adjustl(str))//'f20.10)') 'read value: ', &
    & (tmpval(ii),ii=1,constrainable_coord%ndim)
      deallocate(tmpval)
      read(nfcntn,*)
      read(nfcntn,*) (constrainable_coord%sigma(ii),ii=1,constrainable_coord%ndim)
      write(nfout,'(a,'//trim(adjustl(str))//'f20.10)') 'read sigma: ', &
    &   (constrainable_coord%sigma(ii),ii=1,constrainable_coord%ndim)

      read(nfcntn,*)
      do i=1,constrainable_coord%ndim
        write(nfout,'(a,i5)') 'read dsigma for dimension ',i
        read(nfcntn,*)
        do j=1,constrainable_coord%n_associated_atoms
          read(nfcntn,*)
          read(nfcntn,*) (constrainable_coord%dsigma(i,j,ii),ii=1,3)
          write(nfout,'(3e20.10)') (constrainable_coord%dsigma(i,j,ii),ii=1,3)
        enddo
      enddo

      read(nfcntn,*)
      do i=1,constrainable_coord%ndim
        read(nfcntn,*)
        write(nfout,'(a,i5)') 'read dsigma_old for dimension ',i
        do j=1,constrainable_coord%n_associated_atoms
          read(nfcntn,*)
          read(nfcntn,*) (constrainable_coord%dsigma_old(i,j,ii),ii=1,3)
          write(nfout,'(3e20.10)') (constrainable_coord%dsigma_old(i,j,ii),ii=1,3)
        enddo
      enddo

      read(nfcntn,*)
      read(nfcntn,*)
      read(nfcntn,*)
      read(nfcntn,*)
      read(nfcntn,*)
      read(nfcntn,*) constrainable_coord%det_metric
      write(nfout,'(a)') 'read the determinant of the metric'
      write(nfout,'(e20.10)') constrainable_coord%det_metric

      read(nfcntn,*)
      read(nfcntn,*) (constrainable_coord%lambda(ii),ii=1,constrainable_coord%ndim)
      write(nfout,'(a,'//trim(adjustl(str))//'e20.10)') 'read lambda: ', &
    &   (constrainable_coord%lambda(ii),ii=1,constrainable_coord%ndim)

      read(nfcntn,*)
      read(nfcntn,*)
      if(constrainable_coord%is_reaction_coordinate)then
        read(nfcntn,*)
        read(nfcntn,*) ntmp
        if(ntmp/=constrainable_coord%n_reaction_coords)then
          constrainable_coord%n_reaction_coords = ntmp
          if(associated(constrainable_coord%reaction_coords))then
            deallocate(constrainable_coord%reaction_coords)
          endif
          allocate(constrainable_coord%reaction_coords(   &
    &              constrainable_coord%n_reaction_coords, &
    &              constrainable_coord%ndim))
        endif
        read(nfcntn,*)
        allocate(tmpval2(constrainable_coord%n_reaction_coords,constrainable_coord%ndim));tmpval2=0.d0
        do i=1,constrainable_coord%n_reaction_coords
          read(nfcntn,*) (tmpval2(i,ii),ii=1,constrainable_coord%ndim)
          write(nfout,'(a,'//trim(adjustl(str))//'e20.10)') 'read reaction coordinate: ', &
    &     (tmpval2(i,ii),ii=1,constrainable_coord%ndim)
        enddo
        if(conf_para.and.constrainable_coord%is_reaction_coordinate.and. &
    &      nvariable.eq.1.and.reac_coord_generation==REAC_COORD_VIA_INPUT)then
          read(nfcntn,fmt=*,err=10)
          do i=1,m_cnstr_get_n_reac_coords()
            read(nfcntn,fmt=*,err=10) constrainable_coord%finished(i)
            write(nfout,'(a,l5)') 'status of each reaction coordinate: ', &
    &       constrainable_coord%finished(i)
          enddo
        endif
      endif

    endif

    ntmp = constrainable_coord%ndim*constrainable_coord%n_associated_atoms*3

    if(npes>1)then
      call mpi_bcast(constrainable_coord%sigma,constrainable_coord%ndim,&
    & mpi_double_precision,0,MPI_CommGroup,ierr)
      call mpi_bcast(constrainable_coord%dsigma,ntmp,mpi_double_precision,0,MPI_CommGroup,ierr)
      call mpi_bcast(constrainable_coord%dsigma_old,ntmp,mpi_double_precision,0,MPI_CommGroup,ierr)
      call mpi_bcast(constrainable_coord%det_metric,1,mpi_double_precision,0,MPI_CommGroup,ierr)
      call mpi_bcast(&
    &   constrainable_coord%lambda,constrainable_coord%ndim,mpi_double_precision,0,MPI_CommGroup,ierr)
      if(conf_para.and.constrainable_coord%is_reaction_coordinate.and. &
    &    nvariable.eq.1.and.reac_coord_generation==REAC_COORD_VIA_INPUT)then
        call mpi_bcast(&
    &   constrainable_coord%finished,m_cnstr_get_n_reac_coords(),mpi_logical,0,MPI_CommGroup,ierr)
      endif
    endif

    return
10  continue
    if(printable) write(nfout,'(a)') 'encountered error while reading '//trim(adjustl(F_CNTN))
    call phase_error_with_msg(nfout,'encountered error while reading '//trim(adjustl(F_CNTN)),__LINE__,__FILE__)
  end subroutine m_cnstr_load

end module m_constraints

module m_meta_dynamics

  use m_variables_for_dynamics
  use m_constraints
  use m_Force, only : forc_l, m_Force_wd_force_cps_cpd
  use mpi

  implicit none

  ! state constants
  integer, private :: EVERY_STEP=0
  integer, private :: PER_NGAP=1

  integer, private :: BIAS_AND_FICTITIOUS=0
  integer, private :: BIAS_ONLY=1
  integer, private :: BIAS_GENERATION=2

  integer, private :: FROM_HISTOGRAM=0
  integer, private :: FROM_CONSTRUCTION=1
  integer, private :: NO_BIAS_OUTPUT=2

  ! type for a 'collective variable'
  type collective_variable_t
    integer  :: id
    real(DP) :: mass
    real(DP) :: k
    real(DP),pointer,dimension(:) :: s
    real(DP) :: delta_s_buf
    real(DP) :: delta_s
    real(DP) :: delta_s_at_e
    real(DP),pointer,dimension(:) :: vs
    real(DP) :: vs_old
    real(DP),pointer,dimension(:) :: fs
    real(DP) :: fsold
    real(DP) :: fsspring,fsbias
    real(DP) :: smin,smax
    real(DP) :: ds
    integer  :: ns
    integer  :: icount
    logical  :: control_vs
    real(DP),pointer,dimension(:) :: s_thermo
    real(DP) :: v_thermo, v_thermo_old
    real(DP) :: F_thermo, F_thermo_old
    real(DP) :: m_thermo
    real(DP) :: target_KE
    real(DP),pointer,dimension(:) :: KE
    type(constrainable_coords_t),pointer,dimension(:) :: collective_variable
    logical :: take_abs_on_output
  end type collective_variable_t

  ! input identifiers
  character(len('collective_variable')), parameter, private :: tag_collective_variable='collective_variable'
  character(len('meta_dynamics')), parameter, private :: tag_meta_dynamics='meta_dynamics'
  character(len('meta_dynamics_type')), parameter, private :: tag_meta_dynamics_type='meta_dynamics_type'
  character(len('bias_buildup_scheme')), parameter, private :: tag_bias_buildup_scheme='bias_buildup_scheme'
  character(len('every_step')), parameter, private :: tag_every_step='every_step'
  character(len('update_frequency')), parameter, private :: tag_update_frequency='update_frequency'
  character(len('per_ngap')), parameter, private :: tag_per_ngap='per_ngap'
  character(len('bias_output_scheme')), parameter, private :: tag_bias_output_scheme='bias_output_scheme'
  character(len('from_histogram')), parameter, private :: tag_from_histogram='from_histogram'
  character(len('from_construction')), parameter, private :: tag_from_construction='from_construction'
  character(len('no_bias_output')), parameter, private :: tag_no_bias_output='no_bias_output'
  character(len('max_bias_update')), parameter, private :: tag_max_bias='max_bias_update'
  character(len('bias_only')), parameter, private :: tag_bias_only='bias_only'
  character(len('bias_generation')), parameter, private :: tag_bias_generation='bias_generation'
  character(len('bias_and_fictitious')), parameter, private :: tag_bias_and_fictitious='bias_and_fictitious'

  character(len('bias_potential')), parameter, private :: tag_bias_potential='bias_potential'
  character(len('bias_id')), parameter, private :: tag_bias_id='bias_id'
  character(len('bias_output_frequency')), parameter, private :: tag_bias_output_frequency='bias_output_frequency'
  character(len('output_frequency')), parameter, private :: tag_output_frequency='output_frequency'
  character(len('mass')), parameter, private :: tag_mass='mass'
  character(len('k')), parameter, private :: tag_k='k'
  character(len('height')), parameter, private :: tag_height='height'
  character(len('add_large_bias_at_edge')), parameter, private :: tag_add_large_bias_at_edge='add_large_bias_at_edge'
  character(len('height_at_edge')), parameter, private :: tag_height_at_edge='height_at_edge'
  character(len('delta_s')), parameter, private :: tag_delta_s='delta_s'
  character(len('delta_s_at_edge')), parameter, private :: tag_delta_s_at_edge='delta_s_at_edge'
  character(len('ds')), parameter, private :: tag_ds='ds'

  character(len('smin')), parameter, private :: tag_smin   ='smin'
  character(len('smax')), parameter, private :: tag_smax   ='smax'
  character(len('snorm')), parameter, private :: tag_snorm ='snorm'
  character(len('take_abs_cvar')), parameter, private :: tag_take_abs_cvar='take_abs_cvar'

  character(len('control_velocity')), parameter, private :: tag_control_velocity='control_velocity'
  character(len('mass_thermo')), parameter, private :: tag_mass_thermo='mass_thermo'
  character(len('target_KE')), parameter, private :: tag_target_KE='target_KE'

  character(len('extensive_output')), parameter, private :: tag_extensive_output='extensive_output'
  character(len('output_per_rank')), parameter, private :: tag_output_per_rank='output_per_rank'

  character(len('continuation_strategy')), parameter, private :: tag_continuation_strategy='continuation_strategy'
  character(len('randomize_velocity')), parameter, private :: tag_randomize_velocity='randomize_velocity'
  character(len('scale_velocity')), parameter, private :: tag_scale_velocity='scale_velocity'
  character(len('velocity_scaling_factor')), parameter, private :: tag_velocity_scaling_factor='velocity_scaling_factor'
  character(len('configuration_from_input')), parameter, private :: tag_configuration_from_input='configuration_from_input'

  character(len('output_cvar_every_step')), parameter, private :: tag_output_cvar_every_step= &
  & 'output_cvar_every_step'

  ! for meta-dynamics control
  logical, private :: mtd_enabled=.false.
  integer, private :: meta_dynamics_type
  integer, private :: max_bias
  integer, private :: update_frequency=20
  integer, private :: curr_bias=0
  integer, private :: curr_step=0

  integer,  private :: bias_buildup_scheme
  integer,  private :: bias_output_scheme
  real(DP), private :: height
  real(DP), private :: height_at_edge
  integer,  private :: n_collective_variables
  integer,  private :: maxns
  real(DP), allocatable, dimension(:), private :: Vs
  real(DP), allocatable, dimension(:,:), private :: dVsds
  logical, private :: add_large_bias_at_e = .false.
  type(collective_variable_t), allocatable, dimension(:) :: collective_variables

  ! files
  character(len('bias_potential.data')), parameter, private :: bias_potential_data = 'bias_potential.data'
  character(len('curr_bias_potential.data')), parameter, private :: curr_bias_potential_data = 'curr_bias_potential.data'
  character(len('fparticle_coordinates.data')), parameter, private :: &
 & fparticle_coordinates_data = 'fparticle_coordinates.data'
  character(len('collective_variables.data')), parameter, private :: &
 & collective_variables_data = 'collective_variables.data'
  character(len('bias_potential_parameters.data')), parameter, private :: &
 & bias_potential_parameters_data = 'bias_potential_parameters.data'

  integer, private :: nfbias
  integer, private :: nfcurrbias
  integer, allocatable, dimension(:), private :: nffparticle
  integer, allocatable, dimension(:), private :: nfcvar
  integer, private :: nffparticle_es
  integer, private :: nfcvar_es
  integer, private :: nfbpotparam
  integer, allocatable, dimension(:), private :: nfcoord
  integer, allocatable, dimension(:), private :: nfenergy

  ! the 'true' total energy
  real(DP), private :: true_total_energy

  logical, private :: bias_generation_only=.false.
  integer, private :: bias_output_freq=10

  logical, private :: extensive_output=.false.

  logical, private :: output_per_rank=.false.

  !!real(DP), private, allocatable, dimension(:,:) :: x_gathered,y_gathered,z_gathered
  real(DP), private, allocatable, dimension(:,:,:) :: cps_gathered,cpd_gathered
  real(DP), private, allocatable, dimension(:,:) :: x_gathered,y_gathered,z_gathered
  real(DP), private, allocatable, dimension(:) :: temperature_gathered,energy_gathered,true_energy_gathered

  logical, private :: randomize_velocity = .true.
  logical, private :: scale_velocity = .false.
  real(DP), private :: velocity_scaling_factor = 1.d0
  logical, private :: configuration_from_input = .false.

  logical, private :: skip_bias_update=.false.

  integer, private :: md_step_for_curr_bias=0

  real(DP), private :: etotal_mtd=0.d0

  logical, private :: output_cvar_every_step=.false.

  real(DP), private, allocatable, dimension(:,:) :: cvarbuf
  real(DP), private, allocatable, dimension(:)   :: heightbuf
  real(DP), private, allocatable, dimension(:,:) :: sbuf

  integer, private :: nrank_prevrun = -1

  interface m_mtd_get_fictitious_KE
    module procedure get_fictitious_KE0
    module procedure get_fictitious_KE1
  end interface m_mtd_get_fictitious_KE

  interface m_mtd_get_bias_id_character
    module procedure get_id_char0
    module procedure get_id_char1
    module procedure get_id_char2
  end interface m_mtd_get_bias_id_character

  contains

  subroutine m_mtd_set_bias_gen_only(bgen_only)
    logical, intent(in) :: bgen_only
    bias_generation_only = bgen_only
  end subroutine m_mtd_set_bias_gen_only

  subroutine m_mtd_open_mtd_files()
    logical :: ex
    integer :: nf
    character(len=256) :: tmpstr
    !if(mype/=0) return
    allocate(nffparticle(nrank_conf_()+1))
    allocate(nfcvar(nrank_conf_()+1))
    allocate(nfcoord(nrank_conf_()+1))
    allocate(nfenergy(nrank_conf_()+1))
    if(.not.bias_generation_only) then
      call open_mtd_files()
    else
      nfcvar(1) = get_unused_unitnumber()
      nffparticle(1) = get_unused_unitnumber()
      if(meta_dynamics_type==BIAS_AND_FICTITIOUS)then
        tmpstr=fparticle_coordinates_data
        nf=nffparticle(1)
      else
        tmpstr=collective_variables_data
        nf=nfcvar(1)
      endif
      inquire(file=trim(adjustl(tmpstr)),exist=ex)
      if(.not.ex)then
        if(printable)write(nfout,'(a)') 'the '//trim(adjustl(tmpstr))//' file does not exist.'
        call phase_error_with_msg(nfout,'the '//trim(adjustl(tmpstr))//' file does not exist.',__LINE__,__FILE__)
      endif
      open(nf,file=trim(adjustl(tmpstr)),status='old',position='rewind')
      inquire(file=trim(adjustl(bias_potential_parameters_data)),exist=ex)
      if(.not.ex)then
        if(printable)write(nfout,'(a)') 'the '//trim(adjustl(bias_potential_parameters_data))//' file does not exist.'
        call phase_error_with_msg(nfout, 'the '//trim(adjustl(bias_potential_parameters_data))//' file does not exist.' &
                                 ,__LINE__,__FILE__)
      endif
      nfbpotparam = get_unused_unitnumber()
      open(nfbpotparam,file=trim(adjustl(bias_potential_parameters_data)),status='old',position='rewind')
    endif
  end subroutine m_mtd_open_mtd_files

  subroutine open_mtd_files()
    integer :: i,nf
    character(len=256) :: str, suffix

    if(mype==0)then

      nf=1
      if(output_per_rank) nf=nrank_conf_()

    ! open files
      if(mype==0.and.bias_output_scheme/=NO_BIAS_OUTPUT)then
        nfcurrbias = get_unused_unitnumber()
        call smart_output_open(nfcurrbias,trim(curr_bias_potential_data))
      endif

      if(meta_dynamics_type==BIAS_AND_FICTITIOUS)then
        if(output_cvar_every_step)then
          nffparticle_es = get_unused_unitnumber()
          call smart_output_open(nffparticle_es,trim(fparticle_coordinates_data)//'_every_step'//trim(adjustl(str)))
        endif
        if(mype_conf==0)then
          nffparticle(1) = get_unused_unitnumber()
          call smart_output_open(nffparticle(1),trim(fparticle_coordinates_data))
          if(output_per_rank.and.nrank_conf_()>1)then
            do i=2,nf+1
              nffparticle(i) = get_unused_unitnumber()
              str=''
              write(str,*) i-2
              call smart_output_open(nffparticle(i),trim(fparticle_coordinates_data)//trim(adjustl(str)))
            enddo
          endif
        endif
      endif

      if(mype_conf==0)then
        nfcvar(1) = get_unused_unitnumber()
        call smart_output_open(nfcvar(1),trim(collective_variables_data))
        if(output_per_rank.and.nrank_conf_()>1)then
          do i=2,nf+1
            nfcvar(i) = get_unused_unitnumber()
            str=''
            write(str,*) i-2
            call smart_output_open(nfcvar(i),trim(collective_variables_data)//trim(adjustl(str)))
          enddo
        endif

        nfbpotparam = get_unused_unitnumber()
        call smart_output_open(nfbpotparam,trim(bias_potential_parameters_data))
      endif

      str=''
      if(mype_conf>0)then
        write(str,*) mype_conf
      endif
      if(output_cvar_every_step)then
        nfcvar_es = get_unused_unitnumber()
        call smart_output_open(nfcvar_es,trim(collective_variables_data)//'_every_step'//trim(adjustl(str)))
      endif
      ifile_handle_structure = get_unused_unitnumber()
      call smart_output_open(ifile_handle_structure,trim(F_DYNM)//trim(adjustl(str)))
      call m_IS_wd_speciesname_etc(ifile_handle_structure)

      ifile_handle_energy = get_unused_unitnumber()
      call smart_output_open(ifile_handle_energy,trim(F_ENF)//trim(adjustl(str)))
      write(ifile_handle_energy,'(a)') 'iter_ion, iter, temperature, etotal, emtd'

      if(mype_conf==0)then
        nfcoord(1)=get_unused_unitnumber()
        call smart_output_open(nfcoord(1),trim(F_DYNM)//'_at_bias')
        call m_IS_wd_speciesname_etc(nfcoord(1))
        call flush(nfcoord(1))
        if(output_per_rank.and.nrank_conf_()>1)then
          do i=2,nf+1
            nfcoord(i)=get_unused_unitnumber()
            str=''
            write(str,*) i-2
            call smart_output_open(nfcoord(i),trim(F_DYNM)//'_at_bias'//trim(adjustl(str)))
            call m_IS_wd_speciesname_etc(nfcoord(i))
            call flush(nfcoord(i))
          enddo
        endif

        nfenergy(1)=get_unused_unitnumber()
        call smart_output_open(nfenergy(1),trim(F_ENF)//'_at_bias')
        write(nfenergy(1),'(a)') 'iter_bias, iter, temperature, etotal, etotal_mtd'
        call flush(nfenergy(1))
        if(output_per_rank.and.nrank_conf_()>1)then
          do i=2,nf+1
            str=''
            nfenergy(i)=get_unused_unitnumber()
            write(str,*) i-2
            call smart_output_open(nfenergy(i),trim(F_ENF)//'_at_bias'//trim(adjustl(str)))
            write(nfenergy(i),'(a)') 'iter_bias, iter, temperature, etotal, etotal_mtd'
            call flush(nfenergy(i))
          enddo
        endif
      endif
    endif

  end subroutine open_mtd_files

  integer function m_mtd_get_curr_bias()
    m_mtd_get_curr_bias=curr_bias
  end function m_mtd_get_curr_bias

  subroutine init_fparticle_coords()
    integer :: i
    do i=1,n_collective_variables
      collective_variables(i)%s = 0.d0
      collective_variables(i)%vs = 0.d0
      collective_variables(i)%fs = 0.d0
      collective_variables(i)%fsold = 0.d0
      call update_sigma(collective_variables(i)%collective_variable(mype_conf+1),1)
    enddo
    do i=1,n_collective_variables
      collective_variables(i)%s(mype_conf+1) = collective_variables(i)%collective_variable(mype_conf+1)%sigma(1)
    enddo
  end subroutine init_fparticle_coords

  subroutine m_mtd_parse_input()
    integer :: f_selectBlock, f_selectParentBlock, f_getIntValue
    integer :: iret
    mtd_enabled=.true.
    call m_mtd_init()
    if (f_selectBlock(tag_meta_dynamics)==0) then
      call read_collective_variables()
      call read_mtd_cntrl()
      call read_bias_pot()
      call read_cont_strat()
      iret=f_selectParentBlock()
    else
      if(printable)write(nfout,'(a)') 'tag: '//tag_meta_dynamics//' was not found.'
      if(.not.bias_generation_only)then
      if(printable)write(nfout,'(a)') 'you must define the '//tag_meta_dynamics// &
    & ' block in order to perform a meta-dynamics simulation.'
      else
      if(printable)write(nfout,'(a)') 'you must define the '//tag_meta_dynamics// &
    & ' block in order generate the bias potential'
      endif
      call phase_error_with_msg(nfout,'tag: '//tag_meta_dynamics//' was not found.',__LINE__,__FILE__)
    endif

    call m_mtd_open_mtd_files()

    call alloc_bias_potential()
    call alloc_auxil()

    if(meta_dynamics_type==BIAS_AND_FICTITIOUS) call init_fparticle_coords()

    contains

    subroutine read_cont_strat()
      integer :: f_getRealValue,f_selectBlock,f_selectParentBlock,f_getIntValue
      integer :: iret
      real(DP) :: dret
      character(len=256)::str
      if(f_selectBlock(tag_continuation_strategy)==0)then
        if(f_getIntValue(tag_randomize_velocity,iret)==0)then
          if(iret==0)then
            randomize_velocity=.false.
          else
            randomize_velocity=.true.
          endif
        endif
        if(f_getIntValue(tag_scale_velocity,iret)==0)then
          if(iret==0)then
            scale_velocity=.false.
          else
            scale_velocity=.true.
          endif
        endif
        if(scale_velocity)then
          if(f_getRealValue(tag_velocity_scaling_factor,dret,'')==0) velocity_scaling_factor=dret
        endif
        write(str,*) mype_conf
        if(f_getIntValue(tag_configuration_from_input//trim(adjustl(str)),iret)==0)then
          if(iret==0)then
            configuration_from_input=.false.
          else
            configuration_from_input=.true.
          endif
        endif
        iret=f_selectParentBlock()
      endif
    end subroutine read_cont_strat

    subroutine read_collective_variables()
      integer :: f_getRealValue,f_selectBlock,f_selectParentBlock
      integer :: iret
      integer :: i,j
      type(collective_variable_t) :: tmpcol0,tmpcol1
      character(len=256) :: str
      real(DP) :: dret
      real(DP) :: minv,maxv,delta,maxn

!!      if(.not.bias_generation_only)then
!!        call preprocess_coords()
!!        call calc_atom_configuration()
!!      endif

      ! first we check the number of collective variables defined ...
      i=0
      do while(.true.)
        i = i+1
        write(str,*) i
        if(f_selectBlock(tag_collective_variable//trim(adjustl(str)))==0)then
          iret = f_selectParentBlock()
        else
          exit
        endif
      enddo

      n_collective_variables=i-1
      if(n_collective_variables==0)then
        if(printable)write(nfout,'(a)') 'collective variable undefined.'
        if(printable)write(nfout,'(a)') &
       & 'you must define atleast one collective variable in order to perform a meta-dynamics dimulation.'
        call phase_error_with_msg(nfout,'collective variable undefined.',__LINE__,__FILE__)
      endif

      allocate(collective_variables(n_collective_variables))
      do i=1,n_collective_variables
        allocate(collective_variables(i)%s(nrank_conf_()));collective_variables(i)%s=0.d0
        allocate(collective_variables(i)%vs(nrank_conf_()));collective_variables(i)%vs=0.d0
        allocate(collective_variables(i)%fs(nrank_conf_()));collective_variables(i)%fs=0.d0
        allocate(collective_variables(i)%s_thermo(nrank_conf_()));collective_variables(i)%s_thermo=0.d0
        allocate(collective_variables(i)%collective_variable(nrank_conf_()))
        allocate(collective_variables(i)%KE(nrank_conf_()));collective_variables(i)%KE=0.d0
      enddo

      do i=1,n_collective_variables
        if(.not.bias_generation_only)then
          do j=1,nrank_conf_()
            call m_cnstr_init_and_read_cnstr(&
         &  i,collective_variables(i)%collective_variable(j),tag_collective_variable)
            collective_variables(i)%collective_variable(j)%value=0.d0
            collective_variables(i)%collective_variable(j)%pbc=.true.
            call update_sigma(collective_variables(i)%collective_variable(j),1)
          enddo
        endif
        call set_cval_defaults(collective_variables(i))
      enddo

      if(f_selectBlock(tag_collective_variable)==0)then
        do i=1,n_collective_variables
          call read_collective_variables_sub(i)
        enddo
        iret=f_selectParentBlock()
      endif

      do i=1,n_collective_variables
        write(str,*) i
        if(f_selectBlock(tag_collective_variable//trim(adjustl(str)))==0)then
          call read_collective_variables_sub(i)
          iret = f_selectParentBlock()
        endif
      enddo

      maxn=0
      do i=1,n_collective_variables
        minv  = collective_variables(i)%smin
        maxv  = collective_variables(i)%smax
        if(minv.ge.maxv)then
          if(printable)write(nfout,'(a,i5)') 'detected smin>=smax for collective variable no. ',i
          call phase_error_with_msg(nfout,'detected smin>=smax in the collective variable',__LINE__,__FILE__)
        endif
        delta = collective_variables(i)%ds
        collective_variables(i)%ns = floor((maxv-minv)/delta)+1
        if(collective_variables(i)%ns>maxn) then
          maxn=collective_variables(i)%ns
        endif
      enddo
      maxns = maxn

      ! simple sort according to the number of mesh for each col. var.
      do i=1,n_collective_variables-1
        do j=i+1,n_collective_variables
          if( &
          & collective_variables(i)%ns > &
          & collective_variables(j)%ns) then
            tmpcol0 = collective_variables(i)
            tmpcol1 = collective_variables(j)
            collective_variables(j) = tmpcol0
            collective_variables(i) = tmpcol1
          endif
        enddo
      enddo
      do i=1,n_collective_variables
        collective_variables(i)%id=i
      enddo

      !if(.not.bias_generation_only)then
        !do i=1,nrank_conf_()
        !  call m_constraints_dealloc_constrainable_coord(tmpcol0%collective_variable(i))
        !  call m_constraints_dealloc_constrainable_coord(tmpcol1%collective_variable(i))
        !enddo
      !endif
    end subroutine read_collective_variables

    subroutine set_cval_defaults(collective_variable)
      type(collective_variable_t), intent(inout) :: collective_variable
      collective_variable%mass         = 1.d0
      collective_variable%k            = 0.1d0
      collective_variable%delta_s      = 0.1d0
      collective_variable%smin         = 0.d0
      collective_variable%smax         = 2.d0
      collective_variable%ds           = 0.1d0
      collective_variable%fs           = 0.d0
      collective_variable%fsold        = 0.d0
      collective_variable%fsspring     = 0.d0
      collective_variable%fsbias       = 0.d0
      if(.not.bias_generation_only) &
     & collective_variable%s(mype_conf+1) = collective_variable%collective_variable(mype_conf+1)%sigma(1)
      collective_variable%control_vs   = .false.
      collective_variable%icount       = 0
      collective_variable%s_thermo     = 0.d0
      collective_variable%v_thermo     = 0.d0
      collective_variable%v_thermo_old = 0.d0
      collective_variable%F_thermo     = 0.d0
      collective_variable%F_thermo_old = 0.d0
      collective_variable%m_thermo     = 0.1d0
      collective_variable%target_KE    = 0.2d0
      collective_variable%vs           = 0.d0
      collective_variable%take_abs_on_output = .false.
      if(.not.bias_generation_only) collective_variable%vs_old = collective_variable%vs(mype_conf+1)
    end subroutine set_cval_defaults

    subroutine read_collective_variables_sub(i)
      integer, intent(in) :: i
      integer :: iret
      real(DP) :: dret, factor
      character(len=10) :: tag,uname
      integer :: f_getRealValue,f_getIntValue,f_getBoolValue
      logical :: bret

      tag='';factor=1.d0
      uname = collective_variables(i)%collective_variable(1)%unit_name
      if(uname.eq.'bohr') then
        tag='bohr'
      else if (uname.eq.'radian')then
        tag='degree'
        factor=PAI/180.0d0
      endif
      if(f_getRealValue(tag_mass,dret,'au_mass')==0) collective_variables(i)%mass=dret
      if(f_getRealValue(tag_k,dret,'')==0) collective_variables(i)%k=dret
      if(f_getRealValue(tag_delta_s,dret,tag)==0) then
        collective_variables(i)%delta_s      = dret * factor
        collective_variables(i)%delta_s_buf  = dret * factor
        collective_variables(i)%delta_s_at_e = dret * factor
      endif
      if(f_getRealValue(tag_delta_s_at_edge,dret,tag)==0) collective_variables(i)%delta_s_at_e = dret * factor
      if(f_getRealValue(tag_smin,dret,tag)==0) collective_variables(i)%smin                    = dret * factor
      if(f_getRealValue(tag_smax,dret,tag)==0) collective_variables(i)%smax                    = dret * factor
      if(f_getRealValue(tag_ds,dret,tag)==0) collective_variables(i)%ds                        = dret * factor
      if(f_getRealValue(tag_normalization_factor,dret,'')==0)then
        collective_variables(i)%collective_variable%normalization_factor = dret
      endif
      if(f_getIntValue(tag_control_velocity,iret)==0) then
        if(iret==0)then
          collective_variables(i)%control_vs=.false.
        else
          collective_variables(i)%control_vs=.true.
        endif
      endif
      if(f_getRealValue(tag_mass_thermo,dret,'')==0) then
          collective_variables(i)%m_thermo=dret
      endif
      if(f_getRealValue(tag_target_KE,dret,'')==0) then
          collective_variables(i)%target_KE=dret
      endif
      collective_variables(i)%vs(mype_conf+1) = &
     & dsqrt(collective_variables(i)%target_KE*2.0d0/collective_variables(i)%m_thermo)
      collective_variables(i)%vs_old = collective_variables(i)%vs(mype_conf+1)
      if(f_getBoolValue(tag_take_abs_cvar,bret)==0)then
        collective_variables(i)%take_abs_on_output = bret
      endif
    end subroutine read_collective_variables_sub

    subroutine read_mtd_cntrl()
      integer :: f_getStringValue, f_getIntValue, f_getBoolValue
      integer :: iret
      character(len=256) :: string
      logical :: bret
      if (f_getStringValue(tag_meta_dynamics_type,string,LOWER)==0) then
        if(string==tag_bias_only) meta_dynamics_type=BIAS_ONLY
        if(string==tag_bias_and_fictitious) meta_dynamics_type=BIAS_AND_FICTITIOUS
        if(string==tag_bias_generation) then
          meta_dynamics_type=BIAS_GENERATION
          call m_mtd_set_bias_gen_only(bgen_only = .true.)
        endif
      endif
      if (f_getIntValue(tag_max_bias,iret)==0) max_bias=iret
      if (f_getIntValue(tag_extensive_output,iret)==0)then
        if(iret==1) extensive_output=.true.
      endif
      if (f_getIntValue(tag_output_per_rank,iret)==0)then
        if(iret==1) output_per_rank=.true.
      endif
      if (f_getBoolValue(tag_output_cvar_every_step,bret)==0)then
        output_cvar_every_step=bret
      endif
    end subroutine read_mtd_cntrl

    subroutine read_bias_pot()
      integer :: f_getRealValue, f_selectParentBlock, f_getIntValue, f_selectBlock, f_getStringValue
      integer :: iret
      real(DP) :: dret
      character(len=256) :: string
      if (f_selectBlock(tag_bias_potential)==0) then
        if(f_getIntValue(tag_update_frequency,iret)==0)then
          update_frequency = iret
        endif
        if (f_getStringValue(tag_bias_buildup_scheme,string,LOWER)==0) then
          if(string==tag_per_ngap) bias_buildup_scheme=PER_NGAP
        endif
        if (f_getStringValue(tag_bias_output_scheme,string,LOWER)==0) then
          if(string==tag_from_histogram) bias_output_scheme=FROM_HISTOGRAM
          if(string==tag_no_bias_output) bias_output_scheme=NO_BIAS_OUTPUT
        endif
        if(f_getRealValue(tag_height,dret,'hartree')==0) height=dret
        height_at_edge = height * 100.d0
        if(f_getRealValue(tag_height_at_edge,dret,'hartree')==0) height_at_edge = dret
        if(f_getIntValue(tag_add_large_bias_at_edge,iret)==0)then
          if(iret==0)add_large_bias_at_e=.false.
          if(iret==1)add_large_bias_at_e=.true.
        endif
        if(bias_generation_only) then
          if(f_getIntValue(tag_output_frequency,iret)==0) bias_output_freq=iret
          if(f_getIntValue(tag_bias_output_frequency,iret)==0) bias_output_freq=iret
        endif
        iret=f_selectParentBlock()
      endif
    end subroutine read_bias_pot

    subroutine alloc_bias_potential()
      integer :: i,ntot
      ntot = 1
      do i=1,n_collective_variables
        ntot = ntot * collective_variables(i)%ns
      enddo
      allocate(Vs(ntot))
      Vs=0.d0
      if(bias_buildup_scheme==PER_NGAP)then
        allocate(dVsds(ntot,n_collective_variables))
        dVsds=0.d0
      endif
    end subroutine alloc_bias_potential

    subroutine alloc_auxil()
      !!allocate(x_gathered(natm,nrank_conf_()),y_gathered(natm,nrank_conf_()),z_gathered(natm,nrank_conf_()))
      !!x_gathered=0.d0;y_gathered=0.d0;z_gathered=0.d0
      allocate(cps_gathered(natm,3,nrank_conf_()));cps_gathered=0.d0
      allocate(cpd_gathered(natm,3,nrank_conf_()));cpd_gathered=0.d0
      allocate(x_gathered(natm,nrank_conf_()));x_gathered=0.d0
      allocate(y_gathered(natm,nrank_conf_()));y_gathered=0.d0
      allocate(z_gathered(natm,nrank_conf_()));z_gathered=0.d0
      allocate(temperature_gathered(nrank_conf_()),energy_gathered(nrank_conf_()),true_energy_gathered(nrank_conf_()))
      temperature_gathered=0.d0;energy_gathered=0.d0;true_energy_gathered=0.d0
    end subroutine alloc_auxil

  end subroutine m_mtd_parse_input

  subroutine m_mtd_do_mtd()
    integer :: i,ierr

    etotal_mtd = etotal

    if(meta_dynamics_type==BIAS_AND_FICTITIOUS) call do_fictitious_coords()

    do i=1,n_collective_variables
      call update_sigma(collective_variables(i)%collective_variable(mype_conf+1),1)
      call update_dsigma(collective_variables(i)%collective_variable(mype_conf+1),1)
    enddo

    if(curr_bias.gt.1)then
      if(printable) write(nfout,'(a)') 'building the bias potential for this step...'
      if (bias_buildup_scheme==EVERY_STEP) then
        call build_bias_pot_per_step()
      else if (bias_buildup_scheme==PER_NGAP) then
        call intpl_bias_pot()
      endif
      if(printable) write(nfout,'(a)') '...done.'
    endif


    if(meta_dynamics_type==BIAS_AND_FICTITIOUS) then
      call do_fictitious_energy_and_force()
      call do_fictitious_velocity()
    endif
  end subroutine m_mtd_do_mtd

  subroutine m_mtd_io_per_step()
    character(len=256) :: str
    integer :: j

    if(mype==0)then
       write(ifile_handle_energy,'(i6,i10,5f20.10)') &
     & curr_md_step,iteration,temperature, etotal,H_NVT,m_mtd_get_etotal_mtd(),m_mtd_get_econst_mtd()

      if(output_cvar_every_step)then
        write(str,*) n_collective_variables
        write(nfcvar_es,'(2i8,'//trim(adjustl(str))//'f20.10)') curr_bias+mype_conf, curr_md_step, &
        & (collective_variables(j)%collective_variable(mype_conf+1)%sigma(1),j=1,n_collective_variables)
        call flush(nfcvar_es)

        if(meta_dynamics_type==BIAS_AND_FICTITIOUS)then
          if(extensive_output)then
            write(str,*) n_collective_variables*4+1
            write(nffparticle_es,'(2i8,'//trim(adjustl(str))//'f20.10)') &
          & curr_bias+mype_conf, curr_md_step,                                &
          & (collective_variables(j)%s(mype_conf+1),j=1,n_collective_variables),        &
          & (collective_variables(j)%vs(mype_conf+1),j=1,n_collective_variables),       &
          & (collective_variables(j)%fs(mype_conf+1),j=1,n_collective_variables),       &
          & (collective_variables(j)%s_thermo(mype_conf+1),j=1,n_collective_variables), &
          & (collective_variables(j)%KE(mype_conf+1),j=1,n_collective_variables)
          else
            write(nffparticle_es,'(2i8,'//trim(adjustl(str))//'f20.10)') &
          & curr_bias+mype_conf, curr_md_step, (collective_variables(j)%s(mype_conf+1),j=1,n_collective_variables)
          endif
          call flush(nffparticle_es)
        endif
      endif
    endif
  end subroutine m_mtd_io_per_step

  subroutine do_fictitious_coords()
    integer :: i
    real(DP) :: factor,dt2,massinv
    dt2 = 0.5d0
    do i=1,n_collective_variables
      factor = dt2
      massinv = 1.0d0/(collective_variables(i)%mass)
      if(collective_variables(i)%control_vs)then
        collective_variables(i)%s(mype_conf+1) = collective_variables(i)%s(mype_conf+1) &
   &   +collective_variables(i)%vs(mype_conf+1)+factor*(massinv*collective_variables(i)%fs(mype_conf+1) &
   &   -collective_variables(i)%v_thermo*collective_variables(i)%vs(mype_conf+1))

        collective_variables(i)%s_thermo(mype_conf+1)=collective_variables(i)%s_thermo(mype_conf+1) &
   &   +collective_variables(i)%v_thermo &
   &   +factor*collective_variables(i)%F_thermo/collective_variables(i)%m_thermo
      else
        collective_variables(i)%s(mype_conf+1) = collective_variables(i)%s(mype_conf+1) &
   &   +collective_variables(i)%vs(mype_conf+1)+factor*massinv*collective_variables(i)%fs(mype_conf+1)
      endif
      collective_variables(i)%fsold = collective_variables(i)%fs(mype_conf+1)
    enddo
  end subroutine do_fictitious_coords

  subroutine do_fictitious_energy_and_force()
    integer :: i,j,k
    real(DP) :: vspring
    real(DP),allocatable,dimension(:,:) :: f
    allocate(f(natm,3));f=0.d0
    vspring=0.d0

    do i=1,n_collective_variables
      vspring = vspring &
    &         + 0.5d0*collective_variables(i)%k*( &
    &           collective_variables(i)%collective_variable(mype_conf+1)%sigma(1) &
    &         - collective_variables(i)%s(mype_conf+1))**2
    enddo
    etotal_mtd = etotal_mtd+vspring

    do i=1,n_collective_variables
      do j=1,collective_variables(i)%collective_variable(mype_conf+1)%n_associated_atoms
        k=collective_variables(i)%collective_variable(mype_conf+1)%associated_atoms(j)
        f(k,:) = f(k,:) &
    &          + collective_variables(i)%k*( &
    &            collective_variables(i)%collective_variable(mype_conf+1)%sigma(1)-collective_variables(i)%s(mype_conf+1) &
    &          )*collective_variables(i)%collective_variable(mype_conf+1)%dsigma(1,j,:)
      enddo
    enddo
    forc_l(:,:) = forc_l(:,:) - f(:,:)
    deallocate(f)

    do i=1,n_collective_variables
      collective_variables(i)%fsspring = &
   & -collective_variables(i)%k*( &
   &  collective_variables(i)%s(mype_conf+1)-collective_variables(i)%collective_variable(mype_conf+1)%sigma(1))
      collective_variables(i)%fs(mype_conf+1) = collective_variables(i)%fsspring+collective_variables(i)%fsbias
    enddo

  end subroutine do_fictitious_energy_and_force

  subroutine do_fictitious_velocity()
    integer :: i,j
    real(DP) :: factor,massinv
    real(DP) :: v_thermo_tmp
    real(DP) :: eps = 1.d-12
    do i=1,n_collective_variables
      if(collective_variables(i)%control_vs)then
        do j=1,1000
          factor = 0.5d0
          massinv = 1.0d0/(collective_variables(i)%mass)
          collective_variables(i)%v_thermo = collective_variables(i)%v_thermo_old + &
   &      (collective_variables(i)%F_thermo_old+collective_variables(i)%F_thermo) &
   &      * 0.5d0/collective_variables(i)%m_thermo

          collective_variables(i)%vs(mype_conf+1) = collective_variables(i)%vs_old &
   &     + factor * (massinv*(collective_variables(i)%fs(mype_conf+1)+collective_variables(i)%fsold) &
         - collective_variables(i)%vs(mype_conf+1)*collective_variables(i)%v_thermo &
         - collective_variables(i)%vs_old*collective_variables(i)%v_thermo_old)

          collective_variables(i)%F_thermo =  2.0d0*( &
   &      m_mtd_get_fictitious_KE(collective_variables(i))-collective_variables(i)%target_KE)

          v_thermo_tmp = collective_variables(i)%v_thermo_old+ &
   &     (collective_variables(i)%F_thermo+collective_variables(i)%F_thermo_old)* &
   &     0.5d0/collective_variables(i)%m_thermo

         if(dabs(v_thermo_tmp-collective_variables(i)%v_thermo).lt.eps) exit
        enddo
        collective_variables(i)%vs_old       = collective_variables(i)%vs(mype_conf+1)
        collective_variables(i)%v_thermo_old = collective_variables(i)%v_thermo
        collective_variables(i)%F_thermo_old = collective_variables(i)%F_thermo
      else
        factor = 0.5d0
        massinv = 1.0d0/(collective_variables(i)%mass)
        collective_variables(i)%vs(mype_conf+1) = collective_variables(i)%vs(mype_conf+1) &
   &    + factor * massinv*(collective_variables(i)%fs(mype_conf+1)+collective_variables(i)%fsold)
      endif
    enddo
  end subroutine do_fictitious_velocity

  real(DP) function get_fictitious_KE0()
    integer :: i
    get_fictitious_KE0=0.d0
    do i=1,n_collective_variables
      get_fictitious_KE0 = get_fictitious_KE0 &
   & +0.5d0*collective_variables(i)%mass*collective_variables(i)%vs(mype_conf+1)**2
    enddo
  end function get_fictitious_KE0

  real(DP) function get_fictitious_KE1(collective_variable)
    type(collective_variable_t), intent(in) :: collective_variable
    get_fictitious_KE1 = 0.5d0 * collective_variable%mass * collective_variable%vs(mype_conf+1)**2
  end function get_fictitious_KE1

  subroutine large_bias_at_edge(cvar,v,f)
    real(DP), intent(in), dimension(:) :: cvar
    real(DP), intent(out) :: v
    real(DP), intent(out), dimension(:) :: f
    integer :: i
    real(DP) :: tmpds,tmpsmin,tmpsmax
    v=0;f=0
    do i=1,n_collective_variables
      tmpds = collective_variables(i)%delta_s_at_e
      tmpsmin = collective_variables(i)%smin
      tmpsmax = tmpsmin+collective_variables(i)%ds*(collective_variables(i)%ns-1)
      v = v + dexp(-0.5d0*((cvar(i)-tmpsmin)/tmpds)**2) + dexp(-0.5d0*((cvar(i)-tmpsmax)/tmpds)**2)
      f(i) = -(-1.0d0/(tmpds**2))*(cvar(i)-tmpsmin)+(-1.0d0/(tmpds**2))*(cvar(i)-tmpsmax)
    enddo
    f = f*v
  end subroutine large_bias_at_edge

  subroutine intpl_bias_pot()
    integer :: ipos,ipostmp
    integer :: i,j,k
    real(DP) :: tmpsmin, tmpds, nstmp, ntmp
    real(DP), allocatable, dimension(:) :: tmpcvar
    allocate(tmpcvar(n_collective_variables))
    if(meta_dynamics_type==BIAS_ONLY)then
      ipos=1
      do i=1,n_collective_variables
        tmpcvar(i) = collective_variables(i)%collective_variable(mype_conf+1)%sigma(1)
      enddo
      do i=1,n_collective_variables
        tmpsmin = collective_variables(i)%smin
        tmpds   = collective_variables(i)%ds
        ipostmp= floor((tmpcvar(i)-tmpsmin)/tmpds)+1
        if(ipostmp.le.0.or.ipostmp.gt.collective_variables(i)%ns)then
          if(printable)write(nfout,'(a,i5)')     'WARN: index out of range for collective variable ',i
          if(printable)write(nfout,'(a,f15.10)') 'WARN: value: ',tmpcvar(i)
          if(printable)write(nfout,'(a,i5)')     'WARN: the value is deemed to be at the edge.'
          if(ipostmp.le.0)then
            ipostmp=1
          else
            ipostmp=collective_variables(i)%ns
          endif
        endif
        nstmp = 1
        do j=1,i-1
          nstmp = nstmp * collective_variables(j)%ns
        enddo
        ipos = ipos + ipostmp*nstmp
      enddo
      etotal_mtd = etotal_mtd+Vs(ipos)
      do i=1,n_collective_variables
        do j=1,collective_variables(i)%collective_variable(mype_conf+1)%n_associated_atoms
          k = collective_variables(i)%collective_variable(mype_conf+1)%associated_atoms(j)
          forc_l(k,:) = forc_l(k,:) - dVsds(ipos,i) * &
         & collective_variables(i)%collective_variable(mype_conf+1)%dsigma(1,j,:)
        enddo
      enddo
    endif
    deallocate(tmpcvar)
  end subroutine intpl_bias_pot

  subroutine m_mtd_reset_bias_id()
    curr_bias = 0
  end subroutine m_mtd_reset_bias_id

  logical function m_mtd_increment_step()
    integer :: ntmp
    m_mtd_increment_step=.true.
    md_step_for_curr_bias = md_step_for_curr_bias+1
    if(md_step_for_curr_bias.gt.update_frequency) then
      m_mtd_increment_step=.false.
      md_step_for_curr_bias=0
    else
      ntmp = curr_bias+mype_conf
      if(printable) write(nfout,'(a,i8,a,i8)') 'in bias potential no. ', ntmp,', md step no. ',md_step_for_curr_bias
    endif
  end function m_mtd_increment_step

  logical function m_mtd_update_bias()
    integer :: i,ierr
    real(DP) :: icount0
    if(skip_bias_update)then
      skip_bias_update=.false.
      m_mtd_update_bias=.true.
      return
    endif
    if(curr_bias.le.1)then
      curr_bias = curr_bias+1
    else
      if(nrank_prevrun>0)then
        curr_bias = curr_bias+nrank_prevrun
        nrank_prevrun = -1
      else
        curr_bias = curr_bias+nrank_conf_()
      endif
    endif
    if(printable)then
      if(printable)write(nfout,'(a,i8)') 'bias potential update no. ',curr_bias
      if(printable)write(nfout,*)
    endif
    if(curr_bias.gt.max_bias.and.max_bias.gt.0) then
      m_mtd_update_bias=.false.
      curr_bias=0
      return
    endif

    if (curr_bias.gt.1) then
      if(printable) write(nfout,'(a)') 'outputting meta-dynamics related data at bias update ...'
      call do_gather()
      call wd_bpot_parameters()
      call wd_collective_variable(nfcvar,nffparticle,'')
      call wd_coords_and_energy()
      call cnstrct_and_wd_bias_potential()
      if(nrank_conf_()>1) call mpi_barrier(mpi_comm_world,ierr)
      call rd_and_send_colvars()
      if(printable) write(nfout,'(a)') '...done'
      if(printable) write(nfout,'(a)')
    endif

    m_mtd_update_bias = .true.
    md_step_for_curr_bias = 0
    !!$curr_md_step = 0
    !!$iteration_ionic = 1
    !!$call set_step_id(0)
  end function m_mtd_update_bias

  subroutine do_gather()
    integer :: i,j,ierr
    real(DP) :: icount
    real(DP) :: tmpke
    real(DP), allocatable, dimension(:) :: buftmp
    real(DP), allocatable, dimension(:) :: buf1d
    real(DP), allocatable, dimension(:,:):: bufx,bufy,bufz

    if(meta_dynamics_type==BIAS_AND_FICTITIOUS)then
      tmpke = m_mtd_get_fictitious_KE()
      do i=1,n_collective_variables
        collective_variables(i)%KE(mype_conf+1) = tmpke
      enddo
    endif

    if(nrank_conf_()>1)then
      allocate(buf1d(nrank_conf_()));allocate(buftmp(nrank_conf_()))
      do i=1,n_collective_variables
        do j=0,nrank_conf_()-1
          if(j==mype_conf)cycle
          collective_variables(i)%s(j+1)=0.d0
          if(extensive_output)then
            collective_variables(i)%vs(j+1)=0.d0
            collective_variables(i)%fs(j+1)=0.d0
            collective_variables(i)%s_thermo(j+1)=0.d0
            collective_variables(i)%KE(j+1)=0.d0
          endif
        enddo
        call mpi_allreduce(collective_variables(i)%s,buf1d,nrank_conf_(),mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
        collective_variables(i)%s=buf1d/dble(npes)
        buftmp=0.d0
        buftmp(mype_conf+1) = collective_variables(i)%collective_variable(mype_conf+1)%sigma(1)
        call mpi_allreduce(buftmp,buf1d,nrank_conf_(),mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
        do j=1,nrank_conf_()
          collective_variables(i)%collective_variable(j)%sigma(1)=buf1d(j)/dble(npes)
        enddo
        if(extensive_output)then
          buf1d=0.d0
          call mpi_allreduce(collective_variables(i)%vs,buf1d,nrank_conf_(),mpi_double_precision, &
          & mpi_sum,mpi_comm_world,ierr)
          collective_variables(i)%vs=buf1d/dble(npes)

          buf1d=0.d0
          call mpi_allreduce(collective_variables(i)%fs,buf1d,nrank_conf_(),mpi_double_precision, &
          & mpi_sum,mpi_comm_world,ierr)
          collective_variables(i)%fs=buf1d/dble(npes)

          buf1d=0.d0
          call mpi_allreduce(collective_variables(i)%s_thermo,buf1d,nrank_conf_(),mpi_double_precision, &
          & mpi_sum,mpi_comm_world,ierr)
          collective_variables(i)%s_thermo=buf1d/dble(npes)

          buf1d=0.d0
          call mpi_allreduce(collective_variables(i)%KE,buf1d,nrank_conf_(),mpi_double_precision, &
          & mpi_sum,mpi_comm_world,ierr)
          collective_variables(i)%KE=buf1d/dble(npes)
        endif
      enddo

      allocate(bufx(natm,nrank_conf_()));bufx=0.d0
      allocate(bufy(natm,nrank_conf_()));bufy=0.d0
      allocate(bufz(natm,nrank_conf_()));bufz=0.d0
      x_gathered=0.d0;y_gathered=0.d0;z_gathered=0.d0
      x_gathered(:,mype_conf+1) = cps(:,1)
      y_gathered(:,mype_conf+1) = cps(:,2)
      z_gathered(:,mype_conf+1) = cps(:,3)
      call mpi_allreduce(x_gathered,bufx,natm*nrank_conf_(),mpi_double_precision, mpi_sum,mpi_comm_world,ierr)
      call mpi_allreduce(y_gathered,bufy,natm*nrank_conf_(),mpi_double_precision, mpi_sum,mpi_comm_world,ierr)
      call mpi_allreduce(z_gathered,bufz,natm*nrank_conf_(),mpi_double_precision, mpi_sum,mpi_comm_world,ierr)

      cpd_gathered = 0.d0
      cpd_gathered(:,:,mype_conf+1) = cpd_l(:,:)
      call mpi_allreduce(mpi_in_place,cpd_gathered,natm*nrank_conf_()*3,mpi_double_precision, mpi_sum,mpi_comm_world,ierr)
      cpd_gathered = cpd_gathered/dble(npes)

      buf1d=0.d0;buftmp=0.d0;buftmp(mype_conf+1)=temperature
      call mpi_allreduce(buftmp,buf1d,nrank_conf_(),mpi_double_precision, mpi_sum, mpi_comm_world,ierr)
      temperature_gathered=buf1d/dble(npes)

      buf1d=0.d0;buftmp=0.d0;buftmp(mype_conf+1)=etotal_mtd
      call mpi_allreduce(buftmp,buf1d,nrank_conf_(),mpi_double_precision, mpi_sum, mpi_comm_world,ierr)
      energy_gathered=buf1d/dble(npes)

      buf1d=0.d0;buftmp=0.d0;buftmp(mype_conf+1)=etotal
      call mpi_allreduce(buftmp,buf1d,nrank_conf_(),mpi_double_precision, mpi_sum, mpi_comm_world,ierr)
      true_energy_gathered=buf1d/dble(npes)
      do i=1,natm
        cps_gathered(i,1,:) = bufx(i,:)/dble(npes)
        cps_gathered(i,2,:) = bufy(i,:)/dble(npes)
        cps_gathered(i,3,:) = bufz(i,:)/dble(npes)
      enddo
      deallocate(buf1d);deallocate(buftmp)
      deallocate(bufx,bufy,bufz)
    else
      cps_gathered(:,:,1)=cps(:,:)
      cpd_gathered(:,:,1)=cpd_l(:,:)
      temperature_gathered(1)=temperature
      energy_gathered(1)=etotal_mtd
      true_energy_gathered(1) = etotal
    endif

  end subroutine do_gather

  subroutine build_bias_pot_per_step()
    integer :: i,j,k,l,ierr
    integer :: itmp
    real(DP) :: v,vtmp
    real(DP) :: vspring
    real(DP), allocatable, dimension(:,:) :: f
    real(DP), allocatable, dimension(:) :: fftmp,ff
    real(DP), allocatable, dimension(:) :: tmpcvar,tmpcvar0
    real(DP), allocatable, dimension(:,:) :: tmptmpcvar
    real(DP), allocatable, dimension(:)   :: tmptmpheight,tmpsbuf
    real(DP) :: icount0
    integer :: nf
    integer :: ntmp
    real(DP) :: tmpheight
    allocate(tmpcvar(n_collective_variables));tmpcvar=0.d0
    allocate(tmpcvar0(n_collective_variables));tmpcvar0=0.d0
    allocate(fftmp(n_collective_variables));fftmp=0.d0
    allocate(ff(n_collective_variables));ff=0.d0
    v=0.d0
    ff=0.d0
    do i=1,n_collective_variables
      if(meta_dynamics_type==BIAS_ONLY)then
        tmpcvar0(i) = collective_variables(i)%collective_variable(mype_conf+1)%sigma(1)
      else
        tmpcvar0(i) = collective_variables(i)%s(mype_conf+1)
      endif
    enddo

    !nf=get_unused_unitnumber()
    !if(meta_dynamics_type==BIAS_ONLY)then
      !nf=nfcvar(1)
    !  open(unit=nf,status='old',position='rewind',file=collective_variables_data)
    !else
    !  open(unit=nf,status='old',position='rewind',file=fparticle_coordinates_data)
    !  !nf=nffparticle(1)
    !endif

    !rewind(unit=nf)
    !rewind(unit=nfbpotparam)

    ntmp=curr_bias-1
    if(ntmp.eq.1)ntmp=1+nrank_conf_()-1
    do i=1,ntmp
      tmpcvar(:) = cvarbuf(i,:)
      tmpheight = heightbuf(i)
      do j=1,n_collective_variables
        collective_variables(j)%delta_s_buf=sbuf(i,j)
      enddo
      call multiply_exponent(vtmp,fftmp,tmpcvar0,tmpcvar)
      v = v+vtmp*tmpheight
      ff = ff+fftmp*tmpheight
    enddo

!    do i=2,ntmp
!      read(unit=nf,fmt=*,err=100,end=100) itmp,itmp,(tmpcvar(j),j=1,n_collective_variables)
!      read(unit=nfbpotparam,fmt=*,err=100,end=100) itmp,tmpheight, &
!     & (collective_variables(j)%delta_s_buf,j=1,n_collective_variables)
!      call multiply_exponent(vtmp,fftmp,tmpcvar0,tmpcvar)
!      v = v+vtmp*tmpheight
!      ff = ff+fftmp*tmpheight
!    enddo
!
!    close(nf)

    if(add_large_bias_at_e)then
      call large_bias_at_edge(tmpcvar0,vtmp,fftmp)
      v  = v+vtmp*height_at_edge
      ff = ff+fftmp*height_at_edge
    endif

    if(meta_dynamics_type==BIAS_ONLY)then
      allocate(f(natm,3));f=0.d0
      do i=1,n_collective_variables
        do j=1,collective_variables(i)%collective_variable(mype_conf+1)%n_associated_atoms
          k=collective_variables(i)%collective_variable(mype_conf+1)%associated_atoms(j)
          f(k,:) = f(k,:) + ff(i)*collective_variables(i)%collective_variable(mype_conf+1)%dsigma(1,j,:)
        enddo
      enddo
      forc_l(:,:) = forc_l(:,:) - f(:,:)
      deallocate(f)
    else
      do i=1,n_collective_variables
        collective_variables(i)%fsbias = -ff(i)
      enddo
    endif

    etotal_mtd = etotal_mtd+v

    deallocate(tmpcvar)
    deallocate(tmpcvar0)
    deallocate(fftmp)
    deallocate(ff)

!    return
!
!100 continue
!    if(printable) write(nfout,'(a)') 'failed to read collective variables from file'
!    write(0,'(a,2i8)') 'mype, mype_conf : ',mype,mype_conf
!    stop
  end subroutine build_bias_pot_per_step

  subroutine m_mtd_wd_bias_potential()
    integer :: j,itmp,nf
    character(len=256) :: fname
    real(DP), allocatable, dimension(:) :: var
    real(DP) :: tmpheight
    integer :: nfnf
    integer :: nb,iketa
    allocate(var(n_collective_variables))
    nf=nfcvar(1)
    if(meta_dynamics_type==BIAS_AND_FICTITIOUS) nf=nffparticle(1)
    nfnf=get_unused_unitnumber()
    nb=1
    rewind(unit=nf)
    do
      read(unit=nf,fmt=*,err=1,end=1) itmp,itmp,(var(j),j=1,n_collective_variables)
      nb=nb+1
    enddo
1   continue
    if(printable) write(nfout,'(a,i10)') 'total number of bias potential update : ',nb
    iketa=int(log10(real(nb)))+1
    rewind(unit=nf);rewind(unit=nfbpotparam)
    do
      read(unit=nf,fmt=*,err=2,end=2) itmp,itmp,(var(j),j=1,n_collective_variables)
      do j=1,n_collective_variables
        if(collective_variables(j)%take_abs_on_output) var(j) = abs(var(j))
      enddo
      read(unit=nfbpotparam,fmt=*,end=2) itmp,tmpheight,(collective_variables(j)%delta_s_buf,j=1,n_collective_variables)
      if(max_bias>0.and.itmp.gt.max_bias) exit
      if(printable) write(nfout,'(a,i8)') ' doing bias update no. ',itmp
      call cnstrct_bias_potential(var,tmpheight)
      if(mod(itmp,bias_output_freq)==0)then
        fname = trim(adjustl(bias_potential_data))//trim(adjustl(m_mtd_get_bias_id_character(itmp,iketa)))
        open(nfnf,file=trim(adjustl(fname)),position='rewind')
        call wd_bias_potential(nfnf)
        if(printable) write(nfout,'(a)')  '   written down bias potential to : '//trim(adjustl(fname))
        close(nfnf)
      endif
    enddo
2   continue
    fname = trim(adjustl(curr_bias_potential_data))
    open(nfnf,file=trim(adjustl(fname)),position='rewind')
    call wd_bias_potential(nfnf)
    if(printable) write(nfout,'(a)')  '   written down the most recent bias potential to : '//trim(adjustl(fname))
    close(nfnf)
    deallocate(var)
  end subroutine m_mtd_wd_bias_potential

  subroutine cnstrct_bias_potential(vars,height0)
    real(DP), dimension(:), intent(in) :: vars
    real(DP), intent(in) :: height0
    real(DP) :: vtmp,tmpds
    real(DP) :: ex
    integer  :: i,j,k,l
    integer  :: ntot,ipos,ipostmp
    integer  :: nstmp
    real(DP), allocatable, dimension(:) :: tmpcvar0
    real(DP), allocatable, dimension(:) :: ftmp
    real(DP) :: tmpsmin
    integer :: nf
    if(bias_output_scheme==FROM_HISTOGRAM.and.bias_buildup_scheme/=PER_NGAP)then
      vtmp=0.d0
      ipos=1
      do i=1,n_collective_variables
        tmpsmin = collective_variables(i)%smin
        tmpds   = collective_variables(i)%ds
        nstmp=1
        do j=1,i-1
          nstmp = nstmp * collective_variables(i)%ns
        enddo
        ipostmp = floor((vars(i)-tmpsmin)/tmpds)+1
        if(ipostmp.le.0.or.ipostmp.gt.collective_variables(i)%ns)then
          if(printable) write(nfout,'(a,i5)') 'WARN: index out of range for collective variable ',i
          if(printable) write(nfout,'(a,i5)') 'WARN: the value is deemed to be at the edge.'
          if(ipostmp.le.0)ipostmp=1
          if(ipostmp.gt.collective_variables(i)%ns)ipostmp=collective_variables(i)%ns
        endif
        vtmp = vtmp+height0
        ipos = ipos+nstmp*ipostmp
      enddo
      Vs(ipos) = Vs(ipos)+vtmp
    else
      allocate(ftmp(n_collective_variables));ftmp=0.d0
      allocate(tmpcvar0(n_collective_variables));tmpcvar0=0.d0
      ntot = 1
      do i=1,n_collective_variables
        ntot = ntot * collective_variables(i)%ns
      enddo
      do i=1,ntot
        ipos=1
        do j=1,n_collective_variables
          tmpsmin = collective_variables(j)%smin
          tmpds   = collective_variables(j)%ds
          tmpcvar0(j) = tmpsmin+(collective_variables(j)%icount)*tmpds
          ipostmp=1
          do k=1,j-1
            ipostmp = ipostmp*collective_variables(k)%ns
          enddo
          ipos = ipos+collective_variables(j)%icount*ipostmp
        enddo
        call multiply_exponent(vtmp,ftmp,tmpcvar0,vars)
        Vs(ipos) = Vs(ipos)+vtmp*height0
        if (bias_buildup_scheme==PER_NGAP) then
          dVsds(ipos,:) = dVsds(ipos,:) + ftmp(:)*height0
        endif
        if(add_large_bias_at_e.and.curr_bias.eq.2)then
          call large_bias_at_edge(tmpcvar0,vtmp,ftmp)
          Vs(ipos)  = Vs(ipos)+vtmp*height_at_edge
          if(bias_buildup_scheme==PER_NGAP)then
            dVsds(ipos,:) = dVsds(ipos,:)+ftmp(:)*height_at_edge
          endif
        endif
        nf=-1
        call increment_cvar(collective_variables(1),nf)
      enddo
      deallocate(ftmp)
      deallocate(tmpcvar0)
    endif

    call reset_icount()

  end subroutine cnstrct_bias_potential

  subroutine rd_bias_potential()
    integer :: ns,ntot
    integer :: i,j,k
    real(DP) :: tmpsmin,tmpds
    real(DP), allocatable, dimension(:) :: tmpcvar
    character(len=256) :: str
    integer :: nf, ierr
    logical :: ex

    if(printable) write(nfout,'(a)') '  reading bias potential from file'
    if(mype==0)then
      inquire(file=trim(adjustl(curr_bias_potential_data)),exist=ex)
      if(.not.ex)then
        if(printable) write(nfout,'(a)') 'file: '//curr_bias_potential_data//' does not exist.'
        return
      endif
      ntot = 1
      do i=1,n_collective_variables
        ntot = ntot * collective_variables(i)%ns
      enddo

      rewind(nfcurrbias)
      nf=-1
      write(str,*) n_collective_variables+1
      allocate(tmpcvar(n_collective_variables))
      do i=1,ntot
        do j=1,n_collective_variables
          tmpsmin = collective_variables(j)%smin
          tmpds = collective_variables(j)%ds
          tmpcvar(j) = tmpsmin+(collective_variables(j)%icount)*tmpds
        enddo
        read(nfcurrbias,*,err=10,end=10) (tmpcvar(k),k=1,n_collective_variables),Vs(i)
        call increment_cvar(collective_variables(1),nf)
      enddo
      deallocate(tmpcvar)
10    continue
      call reset_icount()
    endif

    if(npes>1)then
      call mpi_bcast(ntot,1,mpi_integer,0,MPI_CommGroup,ierr)
      call mpi_bcast(Vs,ntot,mpi_double_precision,0,MPI_CommGroup,ierr)
    endif
    if(printable) write(nfout,'(a)') '  ...done'
  end subroutine rd_bias_potential

  subroutine wd_bias_potential(nf)
    integer, intent(in) :: nf
    integer :: ntot
    integer :: i,j,k
    real(DP) :: tmpsmin,tmpds
    real(DP), allocatable, dimension(:) :: tmpcvar
    character(len=256) :: str
    if(.not.bias_generation_only.and.printable) &
   & write(nfout,'(a)') '  outputting the bias potential to : '//trim(curr_bias_potential_data)
    write(str,*) n_collective_variables+1
    allocate(tmpcvar(n_collective_variables))
    ntot = 1
    do i=1,n_collective_variables
      ntot = ntot * collective_variables(i)%ns
    enddo
    do i=1,ntot
      do j=1,n_collective_variables
        tmpsmin = collective_variables(j)%smin
        tmpds = collective_variables(j)%ds
        tmpcvar(j) = tmpsmin+(collective_variables(j)%icount)*tmpds
      enddo
      write(nf,'('//trim(adjustl(str))//'f20.10)') (tmpcvar(k),k=1,n_collective_variables),Vs(i)
      call increment_cvar(collective_variables(1),nf)
    enddo
    deallocate(tmpcvar)
    call reset_icount()
  end subroutine wd_bias_potential

  subroutine cnstrct_and_wd_bias_potential()
    integer  :: ipe
    integer  :: i
    real(DP), allocatable, dimension(:) :: tmpcvar

    if(bias_output_scheme==NO_BIAS_OUTPUT.and.bias_buildup_scheme/=PER_NGAP) return

    allocate(tmpcvar(n_collective_variables));tmpcvar=0.d0
    do ipe=1,nrank_conf_()
      do i=1,n_collective_variables
        if(meta_dynamics_type==BIAS_ONLY)then
          tmpcvar(i) = collective_variables(i)%collective_variable(ipe)%sigma(1)
        else
          tmpcvar(i) = collective_variables(i)%s(ipe)
        endif
      enddo
      call cnstrct_bias_potential(tmpcvar,height)
    enddo

    deallocate(tmpcvar)

    if(mype==0.and.mype_conf==0) then
      nfcurrbias = get_unused_unitnumber()
      call smart_output_open(nfcurrbias,trim(curr_bias_potential_data),.true.)
      !!write(nfbias,'(a)') '# bias potential for bias update no. '//trim(adjustl(m_mtd_get_bias_id_character()))
      call wd_bias_potential(nfcurrbias)
      call flush(nfcurrbias)
      close(nfcurrbias)
    endif
  end subroutine cnstrct_and_wd_bias_potential

  subroutine multiply_exponent(v,f,cvar0,cvar1)
    real(DP), intent(out) :: v
    real(DP), dimension(:), intent(out) :: f
    real(DP), dimension(:), intent(in) :: cvar0,cvar1
    integer :: i,j
    real(DP) :: tmpds,ex
    v=1
    f=0
    do i=1,n_collective_variables
      tmpds = collective_variables(i)%delta_s_buf
      ex = dexp(-0.5d0*((cvar0(i)-cvar1(i))/tmpds)**2)
      v = v*ex
      f(i) = (-1.0d0/(tmpds**2))*(cvar0(i)-cvar1(i))
    enddo
    f=f*v
  end subroutine multiply_exponent

  recursive subroutine increment_cvar(collective_variable,nf)
    type(collective_variable_t), intent(inout) :: collective_variable
    integer, intent(in) :: nf
    integer :: i
    if (collective_variable%icount.eq.collective_variable%ns-1) then
      collective_variable%icount=0
      if (collective_variable%id+1.le.n_collective_variables) then
        call increment_cvar(collective_variables(collective_variable%id+1),nf)
        if(nf>0) write(nf,*) ''
      endif
    else
      collective_variable%icount = collective_variable%icount+1
    endif
  end subroutine

  subroutine reset_icount()
    integer :: i
    do i=1,n_collective_variables
      collective_variables(i)%icount=0
    enddo
  end subroutine reset_icount

  subroutine wd_bpot_parameters()
    character(len=256) :: str
    integer :: i,j,ierr
    if(.not.(mype==0.and.mype_conf==0))return
    if(printable) write(nfout,'(a)') '  outputting bias potential parameters to : '//bias_potential_data
    write(str,*) n_collective_variables+1
    do i=1,nrank_conf_()
      write(nfbpotparam,'(i8,'//trim(adjustl(str))//'f20.10)') &
     & curr_bias+i-1, &
     & height, (collective_variables(j)%delta_s_buf,j=1,n_collective_variables)
    enddo
    call flush(nfbpotparam)
  end subroutine wd_bpot_parameters

  subroutine wd_collective_variable(nfcolvar,nfparticle,tagtag)
    integer, dimension(:), intent(in) :: nfcolvar,nfparticle
    character(len=*), intent(in) :: tagtag
    integer :: i,j,ierr,nf
    character(len=256)::str
    real(DP) :: tmpke
    if(.not.(mype==0.and.mype_conf==0))return

    nf=nrank_conf_()
    !!if(output_per_rank) nf=nrank_conf_()

    write(str,*) n_collective_variables

    if(printable) write(nfout,'(a)') '  outputting collective variables to : '//collective_variables_data//trim(tagtag)
    do i=1,nf
      write(nfcolvar(1),'(2i8,'//trim(adjustl(str))//'f20.10)') curr_bias+i-1, curr_md_step, &
    & (collective_variables(j)%collective_variable(i)%sigma(1),j=1,n_collective_variables)
      if(output_per_rank.and.nrank_conf_()>1)then
        write(nfcolvar(i+1),'(2i8,'//trim(adjustl(str))//'f20.10)') curr_bias+i-1, curr_md_step, &
      & (collective_variables(j)%collective_variable(i)%sigma(1),j=1,n_collective_variables)
      endif
    enddo
    call flush(nfcolvar(1))
    if(output_per_rank.and.nrank_conf_()>1)then
      do i=2,nrank_conf_()+1
        call flush(nfcolvar(i))
      enddo
    endif

    if(meta_dynamics_type==BIAS_AND_FICTITIOUS)then
      if(printable) write(nfout,'(a)') '  outputting the coordinates of the fictious particles to : '// &
      & fparticle_coordinates_data//trim(tagtag)
      if(extensive_output) write(str,*) n_collective_variables*4+1
      do i=1,nf
        if(extensive_output)then
          write(nfparticle(1),'(2i8,'//trim(adjustl(str))//'f20.10)') &
      &   curr_bias+i-1, curr_md_step,                                &
      &   (collective_variables(j)%s(i),j=1,n_collective_variables),        &
      &   (collective_variables(j)%vs(i),j=1,n_collective_variables),       &
      &   (collective_variables(j)%fs(i),j=1,n_collective_variables),       &
      &   (collective_variables(j)%s_thermo(i),j=1,n_collective_variables), &
      &   (collective_variables(j)%KE(i),j=1,n_collective_variables)
          if(output_per_rank.and.nrank_conf_()>1)then
          write(nfparticle(i+1),'(2i8,'//trim(adjustl(str))//'f20.10)') &
      &   curr_bias+i-1, curr_md_step,                                &
      &   (collective_variables(j)%s(i),j=1,n_collective_variables),        &
      &   (collective_variables(j)%vs(i),j=1,n_collective_variables),       &
      &   (collective_variables(j)%fs(i),j=1,n_collective_variables),       &
      &   (collective_variables(j)%s_thermo(i),j=1,n_collective_variables), &
      &   (collective_variables(j)%KE(i),j=1,n_collective_variables)
          endif
        else
          write(nfparticle(1),'(2i8,'//trim(adjustl(str))//'f20.10)') &
        & curr_bias+i-1, curr_md_step, (collective_variables(j)%s(i),j=1,n_collective_variables)
          if(output_per_rank.and.nrank_conf_()>1)then
          write(nfparticle(i+1),'(2i8,'//trim(adjustl(str))//'f20.10)') &
        & curr_bias+i-1, curr_md_step, (collective_variables(j)%s(i),j=1,n_collective_variables)
          endif
        endif
      enddo
      call flush(nfparticle(1))
      if(output_per_rank.and.nrank_conf_()>1)then
        do i=2,nrank_conf_()+1
          call flush(nfparticle(i))
        enddo
      endif
    endif
  end subroutine wd_collective_variable

  subroutine wd_coords_and_energy()
    integer :: i,nn
    if(.not.(mype==0.and.mype_conf==0))return
    if(printable) write(nfout,'(a)') '  outputting coordinates to : '//trim(F_DYNM)//'_at_bias'
    if(printable) write(nfout,'(a)') '  outputting energy to : '//trim(F_ENF)//'_at_bias'
    nn=nrank_conf_()
    do i=1,nn
      write(nfcoord(1),'(" cps and forc at (iter_bias, iter_total = ",i5,i8," )")') curr_bias+i-1, iteration
      call m_Force_wd_force_cps_cpd(nfcoord(1),WITHOUTTAG,cps_gathered(:,:,i),cpd_gathered(:,:,i),natm)
      write(nfenergy(1),'(i6,i10,4f20.10)') curr_bias+i-1,iteration,temperature_gathered(i), &
      & true_energy_gathered(i), energy_gathered(i)
      if(output_per_rank.and.nrank_conf_()>1)then
        write(nfcoord(i+1),'(" cps and forc at (iter_bias, iter_total = ",i5,i8," )")') curr_bias+i-1, iteration
        call m_Force_wd_force_cps_cpd(nfcoord(i+1),WITHOUTTAG,cps_gathered(:,:,i),cpd_gathered(:,:,i),natm)
        write(nfenergy(i+1),'(i6,i10,4f20.10)') curr_bias+i-1,iteration,temperature_gathered(i), &
        & true_energy_gathered(i), energy_gathered(i)
      endif
    enddo
  end subroutine wd_coords_and_energy

  subroutine m_mtd_finalize()
    integer :: i,j
    if(.not.bias_generation_only)then
      do i=1,n_collective_variables
        deallocate(collective_variables(i)%s)
        deallocate(collective_variables(i)%vs)
        deallocate(collective_variables(i)%fs)
        deallocate(collective_variables(i)%s_thermo)
        deallocate(collective_variables(i)%KE)
        do j=1,nrank_conf_()
          call dealloc_constrainable_coord(collective_variables(i)%collective_variable(j))
        enddo
        deallocate(collective_variables(i)%collective_variable)
      enddo
      deallocate(collective_variables)
    endif
    deallocate(Vs)
    if(bias_buildup_scheme==PER_NGAP)then
      deallocate(dVsds)
    endif
    deallocate(x_gathered,y_gathered,z_gathered)
    deallocate(cps_gathered)
    deallocate(cpd_gathered)
    deallocate(temperature_gathered,energy_gathered,true_energy_gathered)
    deallocate(nffparticle)
    deallocate(nfcvar)
    deallocate(nfcoord)
    if(allocated(cvarbuf))   deallocate(cvarbuf)
    if(allocated(heightbuf)) deallocate(heightbuf)
    if(allocated(sbuf))      deallocate(sbuf)
  end subroutine m_mtd_finalize

  subroutine m_mtd_wd_mtd(nfcntn)
    integer, intent(in) :: nfcntn
    integer :: i
    write(nfcntn,'(a)') 'meta_dynamics_variables'
    write(nfcntn,'(a)') '(curr_bias)'
    write(nfcntn,'(i8)') curr_bias
    !write(nfcntn,'(a)') '(curr_step)'
    !write(nfcntn,'(i8)') get_curr_step()
    write(nfcntn,'(a)') '(md_step_for_curr_bias)'
    write(nfcntn,'(i8)') md_step_for_curr_bias
    write(nfcntn,'(a)') '(n_collective_variables)'
    write(nfcntn,'(i5)') n_collective_variables
    write(nfcntn,'(a)') '(coordinates, velocity and force of the fictitious particles)'
    do i=1,n_collective_variables
      write(nfcntn,'(3f20.10)') collective_variables(i)%s(mype_conf+1),collective_variables(i)%vs(mype_conf+1), &
     & collective_variables(i)%fs(mype_conf+1)
    enddo
    write(nfcntn,'(a)') '(coodinates, velocity and force of the thermostat associated with the fictitious particles)'
    do i=1,n_collective_variables
      write(nfcntn,'(3f20.10)') &
    & collective_variables(i)%s_thermo(mype_conf+1), collective_variables(i)%v_thermo, collective_variables(i)%F_thermo
    enddo
    write(nfcntn,'(a)') '(nrank_conf)'
    write(nfcntn,'(i5)') nrank_conf_()
    do i=1,n_collective_variables
      write(nfcntn,'(a,i5,a)') '(collective variable no. ',i,')'
      call m_cnstr_dump(collective_variables(i)%collective_variable(mype_conf+1),nfcntn)
    enddo
  end subroutine m_mtd_wd_mtd

  character(len=256) function get_char_from(ii)
    integer, intent(in) :: ii
    character :: string*256
    character :: ch*256
    integer :: iketa=2
    integer :: itmp
    itmp=int(log10(real(ii)))+1
    if(itmp>4)iketa=itmp
    write(ch,*) iketa
    write(string,'(i'//trim(adjustl(ch))//'.'//trim(adjustl(ch))//')') ii
    get_char_from = trim(adjustl(string))
  end function get_char_from

  subroutine m_mtd_rd_mtd(nfcntn,contfile_matches_rank)
    integer, intent(in) :: nfcntn
    logical, intent(in) :: contfile_matches_rank
    logical :: tag_is_found, EOF_reach
    integer, parameter     :: len_str = 132
    character(len=len_str) :: str
    integer :: ntmp,i, ierr

    if(mype==0) &
    call rewind_to_tag0(nfcntn,len('meta_dynamics_variables'),'meta_dynamics_variables' &
   & , EOF_reach, tag_is_found, str, len_str)
    call mpi_bcast(tag_is_found,1,mpi_logical,0,MPI_CommGroup,ierr)
    if(.not.tag_is_found) then
      if(printable) write(nfout,'(a)') 'tag: meta_dynamics_variables not found in the nfcntn file.'
      return
    else
      if(printable) write(nfout,&
     & '(a)') 'tag: meta_dynamics_variables found; loading meta-dynamics variables from the nfcntn file.'
    endif
    if(mype==0)then
      read(nfcntn,*)
      read(nfcntn,*) curr_bias
      if(printable) write(nfout,'(a,i8)') 'last run terminated at bias-update no.',curr_bias
!!$    if(get_curr_step()==0) curr_bias=curr_bias+1
      read(nfcntn,*)
      read(nfcntn,*) md_step_for_curr_bias
      md_step_for_curr_bias=md_step_for_curr_bias-1
      read(nfcntn,*)
      read(nfcntn,*) ntmp
      if(printable) write(nfout,'(a,i5)') 'number of collective variables read from the nfcntn file:',ntmp
      if(ntmp.ne.n_collective_variables) then
        if(printable) write(nfout,'(a,2i5)') 'ERROR : the number of collective coordinates is inconsistent ',&
        & ntmp,n_collective_variables
        call phase_error_with_msg(nfout,'the number of collective coordinates is inconsistent',__LINE__,__FILE__)
      endif
      read(nfcntn,*)
      do i=1,ntmp
        read(nfcntn,*) collective_variables(i)%s(mype_conf+1),collective_variables(i)%vs(mype_conf+1), &
      &  collective_variables(i)%fs(mype_conf+1)
        if(printable) &
      & write(nfout,'(i5,3f20.10)') i,collective_variables(i)%s(mype_conf+1),collective_variables(i)%vs(mype_conf+1), &
      & collective_variables(i)%fs(mype_conf+1)
      enddo
      read(nfcntn,*)
      do i=1,ntmp
        read(nfcntn,*) collective_variables(i)%s_thermo(mype_conf+1),collective_variables(i)%v_thermo, &
       & collective_variables(i)%F_thermo
          if(printable) write(nfout,'(i5,3f20.10)') &
       & i,collective_variables(i)%s_thermo(mype_conf+1),collective_variables(i)%v_thermo,collective_variables(i)%F_thermo
      enddo
      read(nfcntn,*)
      read(nfcntn,*) nrank_prevrun
    endif

    if(npes>1)then
      call mpi_bcast(curr_bias,1,mpi_integer,0,MPI_CommGroup,ierr)
      call mpi_bcast(md_step_for_curr_bias,1,mpi_integer,0,MPI_CommGroup,ierr)
      do i=1,n_collective_variables
        call mpi_bcast(collective_variables(i)%s(mype_conf+1),1,mpi_double_precision,0,MPI_CommGroup,ierr)
        call mpi_bcast(collective_variables(i)%vs(mype_conf+1),1,mpi_double_precision,0,MPI_CommGroup,ierr)
        call mpi_bcast(collective_variables(i)%fs(mype_conf+1),1,mpi_double_precision,0,MPI_CommGroup,ierr)
        call mpi_bcast(collective_variables(i)%v_thermo,1,mpi_double_precision,0,MPI_CommGroup,ierr)
        call mpi_bcast(collective_variables(i)%F_thermo,1,mpi_double_precision,0,MPI_CommGroup,ierr)
      enddo
      call mpi_bcast(nrank_prevrun,1,mpi_integer,0,MPI_CommGroup, ierr)
    endif

    do i=1,n_collective_variables
      call m_cnstr_load(collective_variables(i)%collective_variable(mype_conf+1),nfcntn)
    enddo

    if(bias_output_scheme/=NO_BIAS_OUTPUT) call rd_bias_potential()

    if(curr_bias.ge.2) call rd_and_send_colvars()

    if (.not.contfile_matches_rank) call adjust_contfile_data()

    skip_bias_update=.true.

  end subroutine m_mtd_rd_mtd

  subroutine adjust_contfile_data()
    real(DP) :: tmpke,tmpke0
    integer :: i,j
    real(DP) :: vstmp
    if(randomize_velocity)then
      tmpke=0.d0
!!$      call init_seed()
      call set_seed(mype_conf)
      if(printable) write(nfout,'(a)') ' randomizing the velocities read from F_CNTN file...'
      if(printable)then
        write(nfout,'(a)') ' --- original velocities ---'
        do i=1,natm
        write(nfout,'(i8,3f20.10)') i,v_old(i,1),v_old(i,2),v_old(i,3)
        enddo
      endif
      do i=1,natm
          if(.not.mobile(i))cycle
          tmpke=tmpke+v_old(i,1)**2+v_old(i,2)**2+v_old(i,3)**2
      enddo
      tmpke=dsqrt(tmpke)
      do i=1,3
        do j=1,natm
          if(.not.mobile(j))cycle
          v_old(j,i)=sqrt(-2.0_DP*log(grnd()))*cos(2.0_DP*PAI*grnd())
        enddo
      enddo
      tmpke0=0
      do i=1,natm
        if(.not.mobile(i))cycle
        tmpke0=tmpke0+v_old(i,1)**2+v_old(i,2)**2+v_old(i,3)**2
      enddo
      if(tmpke0.gt.1e-10) v_old(:,:) = (tmpke/tmpke0) * v_old(:,:)
      if(meta_dynamics_type==BIAS_AND_FICTITIOUS)then
        do i=1,n_collective_variables
          vstmp=sqrt(-2.0_DP*log(grnd()))*cos(2.0_DP*PAI*grnd())
          collective_variables(i)%vs(mype_conf+1) = collective_variables(i)%vs(mype_conf+1)*vstmp
          collective_variables(i)%vs_old    = collective_variables(i)%vs(mype_conf+1)
        enddo
      endif
      if(printable)then
        write(nfout,'(a)') ' --- randomized velocities ---'
        do i=1,natm
        write(nfout,'(i8,3f20.10)') i,v_old(i,1),v_old(i,2),v_old(i,3)
        enddo
      endif
    endif
    if(scale_velocity)then
      if(printable) write(nfout,'(a)') ' scaling the velocities read from F_CNTN file.'
      v_old(:,:) = velocity_scaling_factor * v_old(:,:)
      if(meta_dynamics_type==BIAS_AND_FICTITIOUS)then
        collective_variables(i)%vs(mype_conf+1) = collective_variables(i)%vs(mype_conf+1)*velocity_scaling_factor
        collective_variables(i)%vs_old    = collective_variables(i)%vs(mype_conf+1)
      endif
      if(printable)then
        write(nfout,'(a)') ' --- scaled velocities ---'
        do i=1,natm
        write(nfout,'(i8,3f20.10)') i,v_old(i,1),v_old(i,2),v_old(i,3)
        enddo
      endif
    endif

    if(configuration_from_input)then
      cps=cps_in;pos=pos_in
      write(nfout,'(a)') ' resetted the coordinates read from F_CNTN file to those defined in the input.'
    endif
  end subroutine adjust_contfile_data

  logical function m_mtd_enabled()
    m_mtd_enabled=mtd_enabled
  end function m_mtd_enabled

  logical function m_mtd_gen_bias_mode()
    m_mtd_gen_bias_mode = meta_dynamics_type == BIAS_GENERATION
  end function m_mtd_gen_bias_mode

  subroutine m_mtd_init()
    curr_bias           =  0
    max_bias            = -1
    meta_dynamics_type  =  BIAS_ONLY
    bias_buildup_scheme =  EVERY_STEP
    bias_output_scheme  =  FROM_CONSTRUCTION
    height               =  0.1d0
  end subroutine m_mtd_init

  subroutine m_mtd_print_status()
    character(len=256) :: tmpstr,uname
    integer :: i
    if(.not.bias_generation_only)then

    if(printable) then

    write(nfout,'(a)') ''
    write(nfout,'(a)')          ' -- meta-dynamics related parameters and status --'
    tmpstr = tag_bias_and_fictitious
    if(meta_dynamics_type==BIAS_ONLY) tmpstr=tag_bias_only
    if(meta_dynamics_type==BIAS_GENERATION) tmpstr=tag_bias_generation
    write(nfout,'(a)')          '    dynamics type '//trim(adjustl(tmpstr))
    if(max_bias>1) then
    write(nfout,'(a,i5)')       '    max. number of bias-potential update  ',max_bias
    else
    write(nfout,'(a,i5,a)')     '    max. number of bias-potential update  ',max_bias, &
   & ' (remark : negative value means that the program will not terminate by this criteria)'
    endif
    tmpstr = tag_every_step
    if(bias_buildup_scheme==PER_NGAP) tmpstr=tag_per_ngap
    write(nfout,'(a)')          '    bias potential'
    write(nfout,'(a)')          '      potential build-up scheme '//trim(adjustl(tmpstr))
    tmpstr = tag_from_construction
    if(bias_output_scheme==FROM_HISTOGRAM) tmpstr=tag_from_histogram
    if(bias_output_scheme==NO_BIAS_OUTPUT) tmpstr=tag_no_bias_output
    write(nfout,'(a)')          '      potential output scheme   '//trim(adjustl(tmpstr))
    write(nfout,'(a,f15.10,a)') '      height                      ',height,' [hartree]'
    write(nfout,'(a,l)')        '      add large bias at edge    ',add_large_bias_at_e
    if(add_large_bias_at_e)then
    write(nfout,'(a,f15.10,a)') '      height of the large bias  ',height_at_edge,' [hartree]'
    endif
    write(nfout,'(a,i6)')       '      bias potential update frequency ',update_frequency
    write(nfout,'(a,i5)')       '    number of collective variables',n_collective_variables
    write(nfout,'(a)')          '    fictitious particle and associated collective variables'
    do i=1,n_collective_variables
    uname = ' ['//trim(collective_variables(i)%collective_variable(1)%unit_name)//']'
    write(nfout,'(a,i5)')       '       no. ',i
    write(nfout,'(a)')          '         type                      '// &
   & trim(adjustl(collective_variables(i)%collective_variable(mype_conf+1)%nam))
    write(nfout,'(a)')          '         description            '//&
  &  trim(adjustl(collective_variables(i)%collective_variable(mype_conf+1)%descri))
    write(nfout,'(a,f15.10,a)') '         mass                   ',collective_variables(i)%mass,' [atomic mass]'
    write(nfout,'(a,f15.10)')   '         k                      ',collective_variables(i)%k
    write(nfout,'(a,f15.10,a)') '         delta_s                ',collective_variables(i)%delta_s,trim(uname)
    if(add_large_bias_at_e)then
    write(nfout,'(a,f15.10,a)') '         delta_s_at_e           ',collective_variables(i)%delta_s_at_e,trim(uname)
    endif
    write(nfout,'(a,f15.10,a)') '         smin                   ',collective_variables(i)%smin,trim(uname)
    write(nfout,'(a,f15.10,a)') '         smax                   ',collective_variables(i)%smax,trim(uname)
    write(nfout,'(a,f15.10,a)') '         ds                     ',collective_variables(i)%ds,trim(uname)
    write(nfout,'(a,i8)')       '         ns                     ',collective_variables(i)%ns
    write(nfout,'(a,f15.10,a)') '         sigma                  ', &
   & collective_variables(i)%collective_variable(mype_conf+1)%sigma(1),trim(uname)
    write(nfout,'(a,f15.10)')   '         normalization factor   ', &
   & collective_variables(i)%collective_variable(mype_conf+1)%normalization_factor
    if(meta_dynamics_type==BIAS_AND_FICTITIOUS)then
    write(nfout,'(a,l)')        '         velocity control       ',collective_variables(i)%control_vs
    if(collective_variables(i)%control_vs)then
    write(nfout,'(a,f15.10)')   '           target KE            ',collective_variables(i)%target_KE
    write(nfout,'(a,f15.10)')   '           mass of thermostat   ',collective_variables(i)%m_thermo
    endif
    endif
    enddo
    write(nfout,'(a)')          '    output'
    write(nfout,'(a,l)')        '      output per rank ', output_per_rank
    write(nfout,'(a,l)')        '      extensive output ', extensive_output
    write(nfout,'(a,l)')        '      output collective variables at every MD step ', output_cvar_every_step
    write(nfout,'(a)')          '    continuation strategy (meaningful when nrank_conf>number of continuation files)'
    write(nfout,'(a,l)')        '      randomize velocity ',randomize_velocity
    write(nfout,'(a,l)')        '      scale velocity     ',scale_velocity
    if(scale_velocity)then
    write(nfout,'(a,f10.5)')    '        velocity scalingfactor     ',velocity_scaling_factor
    endif
    write(nfout,'(a,l)')        '      use configuration defined in the input file ',configuration_from_input
    write(nfout,'(a)') ''

    endif

    else

    if(printable)then
    write(nfout,'(a)')          ' -- parameters --'
    tmpstr = tag_bias_and_fictitious
    if(meta_dynamics_type==BIAS_ONLY) tmpstr=tag_bias_only
    if(meta_dynamics_type==BIAS_GENERATION) tmpstr=tag_bias_generation
    write(nfout,'(a)')          '    dynamics type '//trim(adjustl(tmpstr))
    write(nfout,'(a)')          '    bias potential'
    write(nfout,'(a)')          '      potential output scheme    '//trim(adjustl(tmpstr))
    write(nfout,'(a,f15.10,a)') '      height                      ',height,' [hartree]'
    write(nfout,'(a,l)')        '      add large bias at edge     ',add_large_bias_at_e
    if(add_large_bias_at_e)then
    write(nfout,'(a,f15.10)')   '      height of the large bias    ',height_at_edge
    endif
    write(nfout,'(a,i8)')       '      bias output frequency      ',bias_output_freq
    write(nfout,'(a,i5)')       '    number of collective variables',n_collective_variables
    do i=1,n_collective_variables
    write(nfout,'(a,i5)')       '       no. ',i
    write(nfout,'(a)')          '         type                      '// &
   & trim(adjustl(collective_variables(i)%collective_variable(mype_conf+1)%nam))
    write(nfout,'(a)')          '         description            '//&
  &  trim(adjustl(collective_variables(i)%collective_variable(mype_conf+1)%descri))
    write(nfout,'(a,f15.10)')   '         delta_s                ',collective_variables(i)%delta_s
    if(add_large_bias_at_e)then
    write(nfout,'(a,f15.10)')   '         delta_s_at_e           ',collective_variables(i)%delta_s_at_e
    endif
    write(nfout,'(a,f15.10)')   '         smin                   ',collective_variables(i)%smin
    write(nfout,'(a,f15.10)')   '         smax                   ',collective_variables(i)%smax
    write(nfout,'(a,f15.10)')   '         ds                     ',collective_variables(i)%ds
    write(nfout,'(a,i8)')       '         ns                     ',collective_variables(i)%ns
    enddo

    endif

    endif
  end subroutine m_mtd_print_status

  function m_mtd_get_etotal_mtd() result(res)
    real(DP) :: res
    res = etotal_mtd
  end function m_mtd_get_etotal_mtd

  function m_mtd_get_econst_mtd() result(res)
    real(DP) :: res
    res = etotal_mtd + H_NVT - etotal
    if(meta_dynamics_type==BIAS_AND_FICTITIOUS) res = res + m_mtd_get_fictitious_KE()
  end function m_mtd_get_econst_mtd

  integer function m_mtd_get_n_col_val()
    m_mtd_get_n_col_val = n_collective_variables
  end function m_mtd_get_n_col_val

  character(len=256) function get_id_char0()
    get_id_char0=get_id_char1(curr_bias)
  end function get_id_char0

  character(len=256) function get_id_char1(i)
    integer, intent(in) :: i
!!$    character :: string*256
!!$    character :: ch*256
!!$    integer :: iketa=4
!!$    integer :: itmp
    get_id_char1=get_id_char2(i,4)
!!$    itmp=int(log10(real(i)))+1
!!$    if(itmp>4)iketa=itmp
!!$    write(ch,*) iketa
!!$    write(string,'(i'//trim(adjustl(ch))//'.'//trim(adjustl(ch))//')') i
!!$    get_id_char1 = trim(adjustl(string))
  end function get_id_char1

  character(len=256) function get_id_char2(i,iketa)
    integer, intent(in) :: i
    integer, intent(in) :: iketa
    integer :: ik
    character :: string*256
    character :: ch*256
    integer :: itmp
    itmp=int(log10(real(i)))+1
    ik=iketa
    if(itmp>iketa)ik=itmp
    write(ch,*) ik
    write(string,'(i'//trim(adjustl(ch))//'.'//trim(adjustl(ch))//')') i
    get_id_char2 = trim(adjustl(string))
  end function get_id_char2

  integer function nrank_conf_()
    if(nrank_conf<1)then
      nrank_conf_ = 1
      return
    endif
    nrank_conf_ = nrank_conf
  end function nrank_conf_

  subroutine rd_and_send_colvars()
    integer :: ntmp,nf,itmp,i,j,ierr,ntmpbuf
    real(DP), allocatable, dimension(:)   :: work1
    real(DP), allocatable, dimension(:,:) :: work2
    character(len=256) :: str, suffix

    ntmp=0
    if(mype==0.and.mype_conf==0)then
      if(meta_dynamics_type==BIAS_ONLY)then
        nf=nfcvar(1)
      else
        nf=nffparticle(1)
      endif
      rewind(unit=nf)
      do while(.true.)
        read(unit=nf,fmt=*,err=100,end=80) itmp,itmp
        ntmp = ntmp+1
      enddo
80    continue
    endif

    ntmpbuf=0
    call mpi_allreduce(ntmp,ntmpbuf,1,mpi_integer,mpi_sum,mpi_comm_world,ierr)
    ntmp=ntmpbuf

    if(printable) write(nfout,'(a,i8)') 'total number of bias potential update at this point : ',ntmp

    if(allocated(cvarbuf))   deallocate(cvarbuf)
    if(allocated(heightbuf)) deallocate(heightbuf)
    if(allocated(sbuf))      deallocate(sbuf)
    allocate(cvarbuf(ntmp,n_collective_variables));cvarbuf=0.d0
    allocate(heightbuf(ntmp));heightbuf=0.d0
    allocate(sbuf(ntmp,n_collective_variables));sbuf=0.d0
    allocate(work1(ntmp));work1=0.d0
    allocate(work2(ntmp,n_collective_variables));work2=0.d0

    if(mype==0.and.mype_conf==0)then
      if(printable) write(nfout,'(a)') '  reading collective variables from file...'
      if(meta_dynamics_type==BIAS_ONLY)then
        nf=nfcvar(1)
      else
        nf=nffparticle(1)
      endif

      rewind(unit=nf)
      rewind(unit=nfbpotparam)

      do i=1,ntmp
        read(unit=nf,fmt=*,err=100,end=100) itmp,itmp,(cvarbuf(i,j),j=1,n_collective_variables)
        read(unit=nfbpotparam,fmt=*,err=100,end=100) itmp,heightbuf(i),(sbuf(i,j),j=1,n_collective_variables)
      enddo

      call flush(nf)
      call flush(nfbpotparam)
!      close(nf)
!      close(nfbpotparam)
!      str=''
!      if(mype_conf>0)then
!        write(str,*) mype_conf
!      endif
!      if(meta_dynamics_type==BIAS_ONLY)then
!        open(unit=nf,file=trim(collective_variables_data)//trim(adjustl(str)),status='old',position='append')
!      else
!        open(unit=nf,file=trim(fparticle_coordinates_data),status='old',position='append')
!      endif
!      open(unit=nfbpotparam,file=trim(bias_potential_parameters_data),status='old',position='append')
    endif

    if(printable) write(nfout,'(a)') '  sending data to all processes...'

    call mpi_allreduce(heightbuf,work1,ntmp,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
    heightbuf = work1
    call mpi_allreduce(cvarbuf,work2,ntmp*n_collective_variables,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
    cvarbuf = work2
    call mpi_allreduce(sbuf,work2,ntmp*n_collective_variables,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
    sbuf = work2

    if(printable) write(nfout,'(a)') '  ...done'
    deallocate(work1);deallocate(work2)
    return

100 continue
    if(printable) write(nfout,'(a)') 'failed to read collective variables from file'
    call phase_error_with_msg(nfout,'failed to read collective variables from file',__LINE__,__FILE__)
  end subroutine rd_and_send_colvars

end module m_meta_dynamics

module m_velocity_verlet
  use m_Control_Parameters, only : ipri,dtio, icond, iprimd
  use m_Const_Parameters, only : LOWER, UPPER, NOCONV, CONST_kB, WITHOUTTAG, ON, PAI, VERLET, &
 &    T_CONTROL, QUENCHED_MD, DAMPED_MD, INITIAL, VELOCITY_SCALING, AU_TIME
  use m_Ionic_System, only : cpd_l, cps, amion, ityp, m_IS_cps_to_pos, m_IS_cp_cps2cpo,natm2,natm, &
 &    m_IS_pack_all_ions_in_uc,ntyp, imdalg
  use m_Force, only : forc_l, m_Force_wd_force_and_cps, m_Force_wd_force_cps_cpd
  use m_Total_Energy, only : etotal
  use m_variables_for_atoms
  use m_variables_for_dynamics
  use m_constraints, only : ifile_handle_structure,ifile_handle_energy,m_cnstr_get_cnstr_coords,   &
 &    m_cnstr_get_cnstr_vel, m_cnstr_first_step,     &
 &    m_cnstr_constraints_exist, m_cnstr_get_id, m_cnstr_udate_for_cnstr,         &
 &    m_cnstr_update_dsigma_old, m_cnstr_reac_coords_variable,force_of_constraintx, &
 &    force_of_constrainty, force_of_constraintz, reset_md_step, skip_coords
  use m_mtrandom
  use m_Parallelization, only : npes,mype,MPI_CommGroup,mype_conf,conf_para
  use m_IterationNumbers, only : iteration,iteration_ionic
  use m_meta_dynamics, only : m_mtd_enabled, m_mtd_do_mtd, m_mtd_get_curr_bias, m_mtd_increment_step, &
 &    m_mtd_get_etotal_mtd, m_mtd_io_per_step
  use mpi

  implicit none

!  include 'mpif.h'

  real(DP), private, allocatable, dimension(:) :: foldx
  real(DP), private, allocatable, dimension(:) :: foldy
  real(DP), private, allocatable, dimension(:) :: foldz

  real(DP), private :: coord_thermo
  real(DP), private :: coord_thermo_dot
  real(DP), private :: coord_thermo_dot_old
  real(DP), private :: F_thermo
  real(DP), private :: F_thermo_old
  real(DP), private :: H00
  real(DP), private :: KE
  real(DP), private :: constant_of_motion

  real(DP), private :: eps_reservor = 1.d-10
  integer,  private :: max_iteration_reservor = 1000

  real(DP), private, parameter :: fs2au=1.0d-15 / AU_TIME  ! 100.d0/2.418884327
  real(DP), private, parameter :: tau_nitrogen  = 14.1d0*fs2au
  real(DP), private, parameter :: mass_nitrogen = 25532.4958d0
  real(DP), allocatable, private, dimension(:) :: damping_factor
  real(DP), allocatable, private, dimension(:,:) :: f_tilda
  real(DP), allocatable, private, dimension(:) :: f_sd_1st_step
  real(DP), allocatable, private, dimension(:,:) :: cpd_l_sd
  real(DP), private :: n_div_mode
  real(DP), private :: dtio_max

  integer, private :: n_resample
  integer, private :: ncount

  logical,private :: automatic_dt
  logical,private :: resample_damping_parameters
  logical,private :: upper_limit_for_dt

  contains

  real(DP) function get_ke(include_molecule)
    implicit none
    logical, intent(in) :: include_molecule
    logical :: mol
    integer :: i
    real(DP) :: KE
    mol=include_molecule
    KE=0.0d0
    do i=1,natm
      if(.not.mobile(i))cycle
      KE = KE+  amion(ityp(i)) &
     & * (cpd_l(i,1)*cpd_l(i,1)+cpd_l(i,2)*cpd_l(i,2)+cpd_l(i,3)*cpd_l(i,3))
    enddo
    get_ke = KE*0.5d0
  end function get_ke

  real(DP) function get_heat_reservoir_force()
    implicit none
    integer :: i
    real(DP) :: nkbt,KE
    nkbt = Temperature_atom*0.5d0*(Nfree+1)*CONST_kB
    KE = get_ke(.false.)
    get_heat_reservoir_force = 2.d0 * (KE-nkbt)
  end function get_heat_reservoir_force

  real(DP) function get_curr_temperature(KE)
    implicit none
    real(DP), intent(in) :: KE
    get_curr_temperature = KE/(0.5*real(Nfree)*CONST_kB)
  end function get_curr_temperature

  real(DP) function get_constant_of_motion()
    real(DP) :: H
    H = get_ke(.true.)+etotal
    if(imdalg==T_CONTROL)then
      H = H + &
   &     0.5d0 * mass_thermo * coord_thermo_dot * coord_thermo_dot &
   &   + coord_thermo * (Nfree+1)*Temperature_atom*CONST_kB
    endif
    get_constant_of_motion = H
  end function get_constant_of_motion

  subroutine m_vv_parse_input()
    integer :: f_selectBlock,f_selectParentBlock,f_selectTop,f_getIntValue,f_getBoolValue,f_getRealValue
    integer :: iret,ival,i
    logical :: bret
    real(DP) :: minmass,dret
    logical :: matom_exists
    iret = f_selectTop();iret = f_selectBlock('structure_evolution')
    if(f_getIntValue('n_md_step',iret)==0)then
      n_md_step = iret
    endif

    automatic_dt=.false.
    resample_damping_parameters = .false.
    upper_limit_for_dt = .true.
    n_div_mode = 7.0d0
    n_resample = 30
    dtio_max = -1
    if(f_selectBlock('damp')==0)then
      if(f_getBoolValue('resample_damping_parameters',bret)==0)then
        resample_damping_parameters=bret
      endif

      if(resample_damping_parameters) automatic_dt=.true.
      if(f_getBoolValue('automatic_dt',bret)==0)then
        automatic_dt = bret
      endif

      if(f_getRealValue('div',dret,'')==0)then
        n_div_mode = dret
      endif

      if(f_getRealValue('dt_max',dret,'au_time')==0)then
        dtio_max = dret
      endif

      if(f_getIntValue('resample_period',iret)==0)then
        n_resample = iret
      endif

      if(f_getBoolValue('upper_limit_for_dt',bret)==0)then
        upper_limit_for_dt = bret
      endif
      iret = f_selectParentBlock()
    endif
    if(printable) write(nfout,'(a,i8)') 'number of MD steps : ',n_md_step
    iret = f_selectTop()

    if(imdalg==DAMPED_MD)then
      allocate(damping_factor(natm))
      allocate(f_tilda(natm,3));f_tilda=0.d0
      allocate(f_sd_1st_step(natm))
      allocate(cpd_l_sd(natm,3));cpd_l_sd=0.d0

      matom_exists=.false.
      do i=1,natm
        if(imdtyp(i)/=0) then
          matom_exists=.true.
          exit
        endif
      enddo
      if(.not.matom_exists) then
        if(printable) write(nfout,'(a)') '** mobile atoms do not exist'
        return
      endif

      minmass = 1.d+10
      do i=1,natm
        if(amion(ityp(i)).lt.minmass.and.imdtyp(i)/=0) minmass=amion(ityp(i))
      enddo
      if(printable) write(nfout,'(a,f15.5)') 'mass of lightest element:',minmass
      if(automatic_dt) then
        dtio =  2.0d0*(tau_nitrogen/n_div_mode)*dsqrt(minmass/mass_nitrogen)
        if(printable) write(nfout,'(a,f10.5)') &
        & 'automatically determined dt from the lightest element: ',dtio
      endif

      if(dtio_max<0) dtio_max = dtio
      damping_factor = 2*PAI/dtio/n_div_mode
      if(printable) write(nfout,'(a)')       '--- damped md stats ---'
      if(printable) write(nfout,'(a,f10.5)') '  damping factor estimated from dt:',damping_factor(1)
      if(printable) write(nfout,'(a,l5)')     '  resample damping parameters: ',resample_damping_parameters
      if(printable .and. resample_damping_parameters) &
                 &  write(nfout,'(a,i5)')    '  resampling period: ',n_resample
      if(printable .and. resample_damping_parameters) &
                 &  write(nfout,'(a,l5)')     '  impose upper limit for dt:', upper_limit_for_dt
      if(printable.and.upper_limit_for_dt.and.resample_damping_parameters)  &
                 &  write(nfout,'(a,f10.5)') '  max. value of dt : ',dtio_max
    endif
  end subroutine m_vv_parse_input

  subroutine m_vv_io_per_step()
    if (mype.ne.0) return
    write(ifile_handle_structure,'(" cps and forc at (iter_ion, iter_total = " &
         & ,i5,i8," )")') curr_md_step, iteration
    if(is_md()) then
      call m_Force_wd_force_cps_cpd(ifile_handle_structure,WITHOUTTAG,cps,cpd_l,natm)
    else
      call m_Force_wd_force_and_cps(ifile_handle_structure,WITHOUTTAG,cps,natm)
    endif
    if(.not.is_md())then
        write(ifile_handle_energy,'(i6,i10,2f20.10)') curr_md_step,iteration,etotal,forcmx
    else if(is_md().and. .not.m_mtd_enabled())then
      write(ifile_handle_energy,'(i6,i10,4f20.10)') curr_md_step,iteration,KE,temperature, &
     & etotal,constant_of_motion
    else if(is_md().and. m_mtd_enabled())then
      call m_mtd_io_per_step()
!      write(ifile_handle_energy,'(i6,i10,5f20.10)') curr_md_step,iteration,KE,temperature, &
!     & etotal,m_mtd_get_etotal_mtd(),constant_of_motion
    endif
  end subroutine m_vv_io_per_step

  integer function m_vv_get_curr_md_step()
    m_vv_get_curr_md_step = curr_md_step
  end function m_vv_get_curr_md_step

  logical function m_vv_increment_md_step()
    logical :: ret
    curr_md_step = curr_md_step+1
    ret = .true.
    if(curr_md_step>n_md_step.and.n_md_step>0) ret = .false.
    if(.not.ret)then
      if(printable) write(nfout,'(a)') 'reached max. MD steps'
    endif
    if(m_mtd_enabled().and.ret)then
      ret = m_mtd_increment_step()
      if(.not.ret) curr_md_step=curr_md_step-1
    endif
    mdstep_for_this_run = mdstep_for_this_run+1
    m_vv_increment_md_step=ret
  end function m_vv_increment_md_step

  subroutine m_vv_init()
    integer :: i
    real(DP) :: minmass
    allocate(foldx(natm));foldx=0.d0
    allocate(foldy(natm));foldy=0.d0
    allocate(foldz(natm));foldz=0.d0
    allocate(v_old(natm,3))
    v_old = 0.d0

    coord_thermo_dot = 0.d0;coord_thermo_dot_old=0.0d0
    coord_thermo = 0.d0
    H00 = 0.0d0
    F_thermo=0.0d0;F_thermo_old=0.0d0
    curr_md_step = 0
    mdstep_for_this_run = 0
    if(.not. m_mtd_enabled())then
      n_md_step = 1000
    else
      n_md_step = -1
    endif

    if(is_md()) call init_velocities()
    skip_coords = .false.
  end subroutine m_vv_init

  logical function is_resampling()
    is_resampling = (mod(mdstep_for_this_run,n_resample)==1 .or. mod(mdstep_for_this_run,n_resample)==2) &
   & .and. resample_damping_parameters .and. imdalg==DAMPED_MD
  end function is_resampling

  subroutine resamp_damping_factor()
    integer  :: i
    real(DP) :: max_damp
    max_damp = -1000000
    if(.not.resample_damping_parameters)return
    if(mod(mdstep_for_this_run,n_resample)==1)then
      if(printable) write(nfout,'(a)') '  resampling damping parameters, step 1.'
      damping_factor = 1.0d0/dtio
      cpd_l_sd = 0.d0
    else if (mod(mdstep_for_this_run,n_resample)==2)then
      if(printable) write(nfout,'(a)') '  resampling damping parameters, step 2.'
      do i=1,natm
        f_sd_1st_step(i) = dsqrt(dot_product(forc_l(i,:),forc_l(i,:)))
      enddo
    else if(mod(mdstep_for_this_run,n_resample)==3)then
      if(printable) write(nfout,'(a)') '  resampling damping parameters, step 3.'
      do i=1,natm
        if(f_sd_1st_step(i).lt.1.d-12)then
          if(printable) write(nfout,'(a)') 'invalid value for the force obtained from the steepest-descent step;'
          if(printable) write(nfout,'(a)') 'damping parameters will not be changed.'
          damping_factor = 2*PAI/dtio/n_div_mode
          call print_damping_parameters()
          return
        endif
        damping_factor(i) = dsqrt(dabs(dlog((f_sd_1st_step(i)/dsqrt(dot_product(forc_l(i,:),forc_l(i,:)))))))/dtio
        if(damping_factor(i).gt.max_damp) max_damp = damping_factor(i)
      enddo
      dtio = 2.0d0*PAI/max_damp/n_div_mode
      call print_damping_parameters()

      if(dtio.gt.dtio_max.and.upper_limit_for_dt)then
        dtio=dtio_max
        damping_factor = 2*PAI/dtio/n_div_mode
        if(printable) write(nfout,'(a)') 'dt>dt_max; setting dt=dt_max'
        call print_damping_parameters()
      endif
    endif

    contains

    subroutine print_damping_parameters()
      integer :: i
      do i=1,natm
        if(printable) write(nfout,'(a,i6,a,f20.10)') &
        & '  new damping factor for atom ',i,' : ',damping_factor(i)
      enddo
      if(printable) write(nfout,'(a,f10.5)') '  new dt : ',dtio
    end subroutine print_damping_parameters

  end subroutine resamp_damping_factor

  subroutine init_velocities()
    real(DP),dimension(3) :: sumv
    real(DP) :: tmp
    integer :: i,j,ia
    real(DP) :: KE
    real(kind=DP) :: sumt,mcom
    real(kind=DP), dimension(3) :: pcom
    !cpd_l=0.d0

    v_old = cpd_l

    F_thermo = get_heat_reservoir_force()
    F_thermo_old = F_thermo
    KE = get_ke(.true.)
    temperature = get_curr_temperature(KE)
    if(printable) write(nfout,'(a,f10.5)') 'initial temperature : ',temperature
  end subroutine init_velocities

  subroutine correct_momentum()
    integer :: i
    integer :: nn
    real(DP),dimension(3) :: sumv
    nn=0
    do i=1,natm
      if(.not.mobile(i))cycle
      sumv(:) = sumv(:) + cpd_l(i,:)*amion(ityp(i))
      nn=nn+1
    enddo
    do i=1,natm
      if(.not.mobile(i))cycle
      cpd_l(i,:) = cpd_l(i,:)-sumv(:)/real(nn)
      cpd_l(i,:) = cpd_l(i,:)/amion(ityp(i))
    enddo
  end subroutine correct_momentum

  subroutine m_vv_do_dynamics(stat)
    integer, intent(out) :: stat
    real(DP) :: H,factor
    logical :: break,skip_ewald
    integer :: i

!!$    if(curr_md_step==1 .and. m_cnstr_get_id()==1 .and. m_cnstr_constraints_exist())then
    if(curr_md_step==1 .and. m_cnstr_constraints_exist())then
      call m_cnstr_first_step()
    endif

    if(imdalg==DAMPED_MD) call resamp_damping_factor()

    skip_ewald = .false.
    if(.not.skip_coords) then
      if(printable) write(nfout,'(a)') 'updating the coordinates.'
      call do_coords()
    else
      if(printable) write(nfout,'(a)') 'skipped update of coordinates.'
      skip_coords = .false.
      skip_ewald = .true.
    endif
    if(printable) write(nfout,'(a)') 'calculating energy and the force acting on each atom.'
    call scf_do_scf_and_force(break,skip_ewald)

    !! meta-dynamics
    if(m_mtd_enabled()) call m_mtd_do_mtd()

    if(imdalg==DAMPED_MD)then
      do i=1,natm
        factor = 0.5d0*dtio/(amion(ityp(i)))
        f_tilda(i,1) = (1/(1+damping_factor(i)))*(forc_l(i,1)-damping_factor(i)*(cpd_l(i,1)+foldx(i)*factor)/factor)
        f_tilda(i,2) = (1/(1+damping_factor(i)))*(forc_l(i,2)-damping_factor(i)*(cpd_l(i,2)+foldy(i)*factor)/factor)
        f_tilda(i,3) = (1/(1+damping_factor(i)))*(forc_l(i,3)-damping_factor(i)*(cpd_l(i,3)+foldz(i)*factor)/factor)
      enddo
      foldx(:) = f_tilda(:,1)
      foldy(:) = f_tilda(:,2)
      foldz(:) = f_tilda(:,3)
    endif

    stat = 0

    if(break)then
      stat = 2
      return
    endif

    if(imdalg==QUENCHED_MD.or.imdalg==DAMPED_MD)then
      call quench_velocities()
      if(check_convergence()) then
        stat = 1
        return
      endif
    endif

    if(printable) write(nfout,'(a)') 'updating the velocities.'
    if(curr_md_step/=1.or.(m_mtd_enabled().and.m_mtd_get_curr_bias().gt.1)) call do_velocities()


    KE = get_ke(.true.)
    H  = get_constant_of_motion()
    temperature = get_curr_temperature(KE)
    if(curr_md_step==1) then
      H00=H
    endif
    constant_of_motion = H-H00
    H_NVT = H
  end subroutine m_vv_do_dynamics

  subroutine m_vv_finalize()
    deallocate(foldx)
    deallocate(foldy)
    deallocate(foldz)
    deallocate(v_old)
    if(imdalg==DAMPED_MD)then
      deallocate(f_tilda)
      deallocate(damping_factor)
      deallocate(cpd_l_sd)
    endif
  end subroutine m_vv_finalize

  subroutine quench_velocities()
    integer :: i,j
    real(DP) :: forc
    real(DP),dimension(3) :: tmpforc,tmpveloc
    if(imdalg==DAMPED_MD.and.is_resampling())return
    do i=1,natm
      tmpforc(:)  = forc_l(i,:)
      tmpveloc(:) = cpd_l(i,:)
      if (m_cnstr_constraints_exist()) then
        tmpforc(1)=tmpforc(1)+force_of_constraintx(i)
        tmpforc(2)=tmpforc(2)+force_of_constrainty(i)
        tmpforc(3)=tmpforc(3)+force_of_constraintz(i)
      endif
      forc=0.d0
      do j=1,3
        forc = tmpforc(j)**2
      enddo
      if(forc.eq.0)then
        cycle
      endif
      if(dot_product(tmpveloc,tmpforc).lt.0)then
        cpd_l(i,:)=0.d0
        if(printable.and.ipri>=1)write(nfout,'(a,i8)') '    quenched velocity of atom ',i
      endif
    enddo
  end subroutine quench_velocities

  logical function check_convergence()
    integer :: i,j
    real(DP) :: max_constraint_force
    real(DP),dimension(3) :: tmpforc,tmpconstraintforc
    real(DP) :: tmpval
    forcmx = 0.d0
    check_convergence = .false.
    do i=1,natm
      if(.not.mobile(i))cycle
      tmpforc(1)=forc_l(i,1);tmpforc(2)=forc_l(i,2);tmpforc(3)=forc_l(i,3)
      if (m_cnstr_constraints_exist()) then
        tmpforc(1)=tmpforc(1)+force_of_constraintx(i)
        tmpforc(2)=tmpforc(2)+force_of_constrainty(i)
        tmpforc(3)=tmpforc(3)+force_of_constraintz(i)
      endif
      tmpval = 0.d0
      do j=1,3
        tmpval = tmpval+tmpforc(j)**2
      enddo
      if(tmpval.eq.0)cycle
      tmpval = dsqrt(tmpval)
      if(forcmx.lt.tmpval) forcmx=tmpval
    enddo

    if(printable) then
      if(m_cnstr_constraints_exist())then
        write(nfout,'(a,f20.10)') &
    &   ' max. force (including the force of constraint) [hartree/bohr]:',forcmx
      else
        write(nfout,'(a,f20.10)') &
    &   ' max. force [hartree/bohr]:',forcmx
      endif
    endif

    if ( forcmx<max_force ) then
      if(printable)then
        write(nfout,'(a)') ' relaxation converged.'
      endif
      check_convergence = .true.
    endif

  end function check_convergence

  subroutine do_coords()
    integer  :: i,iidim,j
    real(DP) :: dt2
    real(DP) :: factor
    real(DP), dimension(3) :: tmpv
    real(DP),allocatable,dimension(:,:) :: tmpcps
    integer,allocatable,dimension(:) :: tmpityp
    dt2 = 0.5d0*dtio*dtio

    if(curr_md_step/=1.or.(m_mtd_enabled().and.m_mtd_get_curr_bias().gt.1))then
      do i=1,natm
        if ( .not.mobile(i) ) cycle
        factor = dt2/(amion(ityp(i)))
        tmpv(:) = cpd_l(i,:)
        if(is_resampling()) tmpv = cpd_l_sd(i,:)
        cps(i,1)   = cps(i,1)   + dtio * tmpv(1) + forc_l(i,1)*factor
        cps(i,2)   = cps(i,2)   + dtio * tmpv(2) + forc_l(i,2)*factor
        cps(i,3)   = cps(i,3)   + dtio * tmpv(3) + forc_l(i,3)*factor
        x_nopbc(i) = x_nopbc(i) + dtio * tmpv(1) + forc_l(i,1)*factor
        y_nopbc(i) = y_nopbc(i) + dtio * tmpv(2) + forc_l(i,2)*factor
        z_nopbc(i) = z_nopbc(i) + dtio * tmpv(3) + forc_l(i,3)*factor
      enddo
      foldx(:) = forc_l(:,1)
      foldy(:) = forc_l(:,2)
      foldz(:) = forc_l(:,3)

    ! NoseNoseNose
      if(imdalg==T_CONTROL)then
        do i=1,natm
          if( .not. mobile(i) ) cycle
          cps(i,:)   = cps(i,:)   - dt2*cpd_l(i,:)*coord_thermo_dot
          x_nopbc(i) = x_nopbc(i) - dt2*cpd_l(i,1)*coord_thermo_dot
          y_nopbc(i) = y_nopbc(i) - dt2*cpd_l(i,2)*coord_thermo_dot
          z_nopbc(i) = z_nopbc(i) - dt2*cpd_l(i,3)*coord_thermo_dot
        enddo
        coord_thermo = coord_thermo + dtio * coord_thermo_dot + &
     &                 F_thermo*(dt2/mass_thermo)
      endif
    endif

    if(printable .and. iprimd>=2)then
       write(nfout,'(a)') 'coordinates before enforcing constraints'
       do i=1,natm
          write(nfout,'(3f25.15)') cps(i,1),cps(i,2),cps(i,3)
       enddo
    endif
    ! NoseNoseNose
    if(m_cnstr_constraints_exist()) then
      ! RATTLE RATTLE RATTLE
      call m_cnstr_get_cnstr_coords(cps)
      call m_cnstr_udate_for_cnstr()
      call m_cnstr_update_dsigma_old()
      ! RATTLE RATTLE RATTLE
    endif
    if(printable .and. iprimd>=2)then
       write(nfout,'(a)') 'coordinates after update'
       do i=1,natm
          write(nfout,'(3f25.15)') cps(i,1),cps(i,2),cps(i,3)
       enddo
    endif
    call m_IS_cps_to_pos()
  end subroutine do_coords

  subroutine do_velocities()
    integer :: i,j,i1
    real(DP) :: factor
    real(DP) :: massinv
    real(DP), dimension(3) :: tmpv
    real(DP) :: tmpsum
    real(DP) :: coord_thermo_dot_tmp

    if (imdalg==VERLET.or.imdalg==QUENCHED_MD.or.imdalg==VELOCITY_SCALING)then
      do i=1,natm
        if ( .not.mobile(i) ) cycle
        factor = dtio*0.5d0/(amion(ityp(i)))
        cpd_l(i,1) = cpd_l(i,1) + factor*(forc_l(i,1)+foldx(i))
        cpd_l(i,2) = cpd_l(i,2) + factor*(forc_l(i,2)+foldy(i))
        cpd_l(i,3) = cpd_l(i,3) + factor*(forc_l(i,3)+foldz(i))
      enddo
    endif
    if(imdalg==DAMPED_MD)then
      do i=1,natm
        if ( .not.mobile(i) ) cycle
        factor = dtio*0.5d0/(amion(ityp(i)))
        if(.not.is_resampling())then
          cpd_l(i,1) = cpd_l(i,1) + factor*(f_tilda(i,1)+foldx(i))
          cpd_l(i,2) = cpd_l(i,2) + factor*(f_tilda(i,2)+foldy(i))
          cpd_l(i,3) = cpd_l(i,3) + factor*(f_tilda(i,3)+foldz(i))
        else
          cpd_l_sd(i,1) = cpd_l_sd(i,1) + factor*(f_tilda(i,1)+foldx(i))
          cpd_l_sd(i,2) = cpd_l_sd(i,2) + factor*(f_tilda(i,2)+foldy(i))
          cpd_l_sd(i,3) = cpd_l_sd(i,3) + factor*(f_tilda(i,3)+foldz(i))
        endif
      enddo
    endif

    if(m_cnstr_constraints_exist())then
      !RATTLE RATTLE RATTLE
      if(imdalg==VERLET.or.imdalg==QUENCHED_MD.or.imdalg==DAMPED_MD.or.imdalg==VELOCITY_SCALING) &
      &  call m_cnstr_get_cnstr_vel(cpd_l)
      if(imdalg==T_CONTROL) call m_cnstr_get_cnstr_vel(v_old)
      !RATTLE RATTLE RATTLE
    endif

    if(imdalg==VELOCITY_SCALING)then
      call scale_velocity()
    endif

    if(imdalg==T_CONTROL)then
      !NoseNoseNose
      do i=1,max_iteration_reservor
        coord_thermo_dot = coord_thermo_dot_old + &
      &   dtio*(F_thermo_old+F_thermo)*0.5d0/mass_thermo
        do j=1,natm
          if ( .not.mobile(j) ) cycle
          factor = dtio*0.5d0
          massinv = 1.d0/amion(ityp(j))
          cpd_l(j,1) = v_old(j,1) + factor * (massinv*forc_l(j,1)-cpd_l(j,1)*coord_thermo_dot &
        &     + massinv*foldx(j)-v_old(j,1)*coord_thermo_dot_old)
          cpd_l(j,2) = v_old(j,2) + factor * (massinv*forc_l(j,2)-cpd_l(j,2)*coord_thermo_dot &
        &     + massinv*foldy(j)-v_old(j,2)*coord_thermo_dot_old)
          cpd_l(j,3) = v_old(j,3) + factor * (massinv*forc_l(j,3)-cpd_l(j,3)*coord_thermo_dot &
        &     + massinv*foldz(j)-v_old(j,3)*coord_thermo_dot_old)
        enddo

        !RATTLE RATTLE RATTLE
        if(m_cnstr_constraints_exist()) then
          call m_cnstr_get_cnstr_vel(cpd_l)
        endif
        !RATTLE RATTLE RATTLE

        F_thermo = get_heat_reservoir_force()
        coord_thermo_dot_tmp = coord_thermo_dot_old + &
     &   dtio*(F_thermo+F_thermo_old)*0.5d0/mass_thermo
        if (dabs(coord_thermo_dot_tmp-coord_thermo_dot).lt.eps_reservor) exit

        if(i.eq.max_iteration_reservor)then
          write(0,*) 'reservoir iteration unconverged'
          write(0,*) dabs(coord_thermo_dot_tmp-coord_thermo_dot)
          call phase_error_with_msg(nfout,'reservoir iteration unconverged',__LINE__,__FILE__)
        endif

      enddo
      if(printable.and.iprimd>=2)write(nfout,'(a,i5)') 'done reservoir iteration. iteration needed for convergence: ',i

      v_old = cpd_l
      coord_thermo_dot_old = coord_thermo_dot
      F_thermo_old = F_thermo
      !NoseNoseNose
    endif

  end subroutine do_velocities

  subroutine scale_velocity()
    real(DP) :: tempera,kina
    tempera=Temperature_atom*0.5d0*(Nfree)*CONST_kB
    kina=get_ke(.false.)
    if(kina.gt.1.d-12) cpd_l(:,:) = cpd_l(:,:) * dsqrt(tempera/kina)
  end subroutine scale_velocity

  subroutine m_vv_wd_vv_variables(nfcntn)
    integer, intent(in) :: nfcntn
    integer :: i
    if(mype/=0)return
    write(nfcntn,'(a)') 'velocity_verlet_variables'
    write(nfcntn,'(a)') '(curr_md_step)'
    write(nfcntn,'(i8)') curr_md_step

    !write(nfcntn,'(a)') '(f)'
    !do i=1,natm
    !  write(nfcntn,'(3e20.10)') forc_l(i,1),forc_l(i,2),forc_l(i,3)
    !enddo
    write(nfcntn,'(a)') '(fold)'
    do i=1,natm
      write(nfcntn,'(3e20.10)') foldx(i),foldy(i),foldz(i)
    enddo
    write(nfcntn,'(a)') '(coord_thermo)'
    write(nfcntn,'(e20.10)') coord_thermo
    write(nfcntn,'(a)') '(coord_thermo_dot)'
    write(nfcntn,'(e20.10)') coord_thermo_dot
    write(nfcntn,'(a)') '(coord_thermo_dot_old)'
    write(nfcntn,'(e20.10)') coord_thermo_dot_old
    write(nfcntn,'(a)') '(F_thermo)'
    write(nfcntn,'(e20.10)') F_thermo
    write(nfcntn,'(a)') '(F_thermo_old)'
    write(nfcntn,'(e20.10)') F_thermo_old
    write(nfcntn,'(a)') '(H0)'
    write(nfcntn,'(e20.10)') H00
    !write(nfcntn,'(a)') '(v)'
    !do i=1,natm
    !  write(nfcntn,'(3e20.10)') cpd_l(i,1:3)
    !enddo
    write(nfcntn,'(a)') '(v_old)'
    do i=1,natm
      write(nfcntn,'(3e20.10)') v_old(i,1),v_old(i,2),v_old(i,3)
    enddo
    write(nfcntn,'(a)') '(pos_nopbc)'
    do i=1,natm
      write(nfcntn,'(3e20.10)') x_nopbc(i),y_nopbc(i),z_nopbc(i)
    enddo

  end subroutine m_vv_wd_vv_variables

  subroutine m_vv_rd_vv_variables(nfcntn,skip)
    integer, intent(in) :: nfcntn
    logical, intent(in) :: skip
    logical :: tag_is_found, EOF_reach
    integer, parameter     :: len_str = 132
    character(len=len_str) :: str
    integer :: i,ierr

    skip_coords = skip ! this should be done regardless of the status of the continuation file
    curr_md_step = iteration_ionic-1

    if(mype==0)then
      call rewind_to_tag0(nfcntn,len('velocity_verlet_variables'),'velocity_verlet_variables' &
     & , EOF_reach, tag_is_found, str, len_str)
    endif

    call mpi_bcast(tag_is_found,1,mpi_logical,0,MPI_CommGroup,ierr)

    if(.not.tag_is_found) then
      if(printable) write(nfout,'(a)') 'tag: velocity_verlet_variables not found in the nfcntn file.'
      return
    else
      if(printable) write(nfout,'(a)') 'tag: velocity_verlet_variables found; &
    &  loading velocity-verlet variables from the nfcntn file.'
    endif

    if(mype==0)then
      read(nfcntn,*)
      read(nfcntn,*) curr_md_step
      write(nfout,'(a,i8)') 'last run terminated at md step no. ',curr_md_step
      curr_md_step = curr_md_step-1

      !write(nfout,'(a)') 'force'
      !read(nfcntn,*)
      !do i=1,natm
      !  read(nfcntn,*) forc_l(i,1),forc_l(i,2),forc_l(i,3)
      !  write(nfout,'(3e20.10)') forc_l(i,1),forc_l(i,2),forc_l(i,3)
      !enddo

      write(nfout,'(a)') 'force (old)'
      read(nfcntn,*)
      do i=1,natm
        read(nfcntn,*) foldx(i),foldy(i),foldz(i)
        write(nfout,'(3e20.10)') foldx(i),foldy(i),foldz(i)
      enddo

      read(nfcntn,*)
      read(nfcntn,*) coord_thermo
      write(nfout,'(a)') 'coordinate of the thermostat'
      write(nfout,'(e20.10)') coord_thermo

      read(nfcntn,*)
      read(nfcntn,*) coord_thermo_dot
      write(nfout,'(a)') 'the time-derivative for the coordinate of the thermostat'
      write(nfout,'(e20.10)') coord_thermo_dot

      read(nfcntn,*)
      read(nfcntn,'(e20.10)') coord_thermo_dot_old
      write(nfout,'(a)') 'the time-derivative for the coordinate of the thermostat (old)'
      write(nfout,'(e20.10)') coord_thermo_dot_old

      write(nfout,'(a)') 'force acting on the thermostat'
      read(nfcntn,*)
      read(nfcntn,*) F_thermo
      write(nfout,'(e20.10)') F_thermo

      write(nfout,'(a)') 'force acting on the thermostat (old)'
      read(nfcntn,*)
      read(nfcntn,*) F_thermo_old
      write(nfout,'(e20.10)') F_thermo_old

      write(nfout,'(a)') 'H0'
      read(nfcntn,*)
      read(nfcntn,*) H00
      write(nfout,'(e20.10)') H00

      !write(nfout,'(a)') 'velocity'
      !read(nfcntn,*)
      !do i=1,natm
      !  read(nfcntn,*) cpd_l(i,1:3)
      !  write(nfout,'(3e20.10)') cpd_l(i,1:3)
      !enddo

      write(nfout,'(a)') 'velocity (old)'
      read(nfcntn,*)
      do i=1,natm
        read(nfcntn,*) v_old(i,1),v_old(i,2),v_old(i,3)
        write(nfout,'(3e20.10)') v_old(i,1),v_old(i,2),v_old(i,3)
      enddo

      write(nfout,'(a)') 'coordinates (pbc not taken into account)'
      read(nfcntn,*)
      do i=1,natm
        read(nfcntn,*) x_nopbc(i),y_nopbc(i),z_nopbc(i)
        write(nfout,'(3e20.10)') x_nopbc(i),y_nopbc(i),z_nopbc(i)
      enddo

    endif

    if(npes>1) then
      call mpi_bcast(curr_md_step,1,mpi_integer,0,MPI_CommGroup,ierr)
      !call mpi_bcast(forc_l,natm*3,mpi_double_precision,0,MPI_CommGroup,ierr)
      call mpi_bcast(foldx,natm,mpi_double_precision,0,MPI_CommGroup,ierr)
      call mpi_bcast(foldy,natm,mpi_double_precision,0,MPI_CommGroup,ierr)
      call mpi_bcast(foldz,natm,mpi_double_precision,0,MPI_CommGroup,ierr)
      call mpi_bcast(coord_thermo,1,mpi_double_precision,0,MPI_CommGroup,ierr)
      call mpi_bcast(coord_thermo_dot,1,mpi_double_precision,0,MPI_CommGroup,ierr)
      call mpi_bcast(coord_thermo_dot_old,1,mpi_double_precision,0,MPI_CommGroup,ierr)
      call mpi_bcast(F_thermo,1,mpi_double_precision,0,MPI_CommGroup,ierr)
      call mpi_bcast(F_thermo_old,1,mpi_double_precision,0,MPI_CommGroup,ierr)
      call mpi_bcast(H00,1,mpi_double_precision,0,MPI_CommGroup,ierr)
      !call mpi_bcast(cpd_l,natm*3,mpi_double_precision,0,MPI_CommGroup,ierr)
      call mpi_bcast(v_old,natm*3,mpi_double_precision,0,MPI_CommGroup,ierr)
      call mpi_bcast(x_nopbc,natm,mpi_double_precision,0,MPI_CommGroup,ierr)
      call mpi_bcast(y_nopbc,natm,mpi_double_precision,0,MPI_CommGroup,ierr)
      call mpi_bcast(z_nopbc,natm,mpi_double_precision,0,MPI_CommGroup,ierr)
    endif

    reset_md_step = .false.
    iteration_ionic = curr_md_step+1

  end subroutine m_vv_rd_vv_variables

  subroutine m_vv_set_skip_coords(skip)
    logical, intent(in) :: skip
     skip_coords = skip
  end subroutine m_vv_set_skip_coords

end module m_velocity_verlet

integer function f_getBoolValue(tag,ret)
  implicit none
  character(*),intent(in) :: tag
  logical,intent(out)     :: ret
  integer :: integ,iret,f_getIntValue
  iret = f_getIntValue(tag,integ)
  if ( integ==0 ) then
    ret = .false.
  else
    ret = .true.
  endif
  f_getBoolValue = iret
end function f_getBoolValue
#endif
! === "phase" repository merge: To make 3D_Parallel by tkato 2011/11/10 ========
!!!!BRANCH_P_END ORG_Parallel
! ==============================================================================
