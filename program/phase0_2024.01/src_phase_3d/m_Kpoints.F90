!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 629 $)
!
!  MODULE: m_Kpoint
!
!  AUTHORS: T. Yamasaki, T. Yamamoto,  August/20/2003
!
!  Contact address :  Phase System Consortium
!                     E-mail: phase_system@nims.go.jp URL https://azuma.nims.go.jp
!
!
!
!=======================================================================
!  patch 10.1 by K. Tagami @adv    2011/06/18
!
!  patch 10.1 :  k_sample_mesh is changed from private to public
!=======================================================================
!
!     The original version of this set of the computer programs "PHASE"
!  was developOPed by the members of the Theory Group of Joint Research
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
module m_Kpoints
!     (m_Kp)
!  $Id: m_Kpoints.F90 629 2020-07-07 02:21:33Z ktagami $
!
  use m_Crystal_Structure,   only : il,imag,inv,ngen,igen,jgen, a,b,c,ca,cb,cc &
       &                          , altv, rltv, nbztyp , nbztyp_spg, n1_sc, n2_sc, n3_sc
  use m_Files,               only : nfout,nfkpoint, nfmatbp, m_Files_open_kpoint_files
  use m_Timing,              only : tstatc0_begin, tstatc0_end
  use m_Control_Parameters,  only : nspin, ipri, tag_accuracy, icond &
       &                          , ipriinputfile, ipri_kp, sw_berry_phase, ekmode, printable, kimg &
       &                          , fixed_charge_k_parallel, sw_hybrid_functional &
       &                          , m_CtrlP_way_of_smearing, m_CtrlP_set_way_ksample,way_of_smearing, way_ksample &
       &                          , sw_modified_kpoint_increment
  use m_Const_Parameters,    only : CARTS, BUCS, CRDTYP, NONAME &
       &                          , MESH, MONKHORST_PACK &
       &                          , SKPS_DIRECT_IN, GAMMA, FILE, DP, PAI2 &
       &                          , SURFACE, WHOLE_BZ, SIMPLE_CUBIC &
       &                          , BCC, FCC, DIAMOND, HEXAGONAL, ORTHORHOMBIC &
       &                          , RUTILE, C2v_SURFACE, Csy_SURFACE &
       &                          , GENERAL, GENERAL_LARGER &
       &                          , HEX1fold, HEX2fold, HEX3fold &
       &                          , LOWER, NODATA, FMAXVALLEN &
       &                          , TETRAHEDRON, ON, OFF, GRID, DELTA &
       &                          , GAMMA_base_symmetrization, PARA, ONE_BY_ONE
  use m_Parallelization,     only : MPI_CommGroup,npes,ierr,nrank_k,mype, myrank_k, ista_kv3_ek, nis_kv3_ek

! ================================= added by. Tagami =========== 11.0
  use m_Control_Parameters,    only : ndim_spinor, noncol, SpinOrbit_mode
  use m_Const_Parameters,      only : Neglected
! ============================================================= 11.0

! ================================= added by. Tagami =========== 12.0A
  use m_Crystal_Structure,   only : gen_name_in_carts, use_altv_rltv
! ============================================================== 12.0A

! === KT_add ==== 2014/08/14+
  use m_CS_Magnetic,  only : magmom_dir_inversion_opr_flag, sw_neglect_magmom
! =============== 2014/08/14+

! === KT_add ==== 2014/09/30, 2018/02/26
  use m_Const_Parameters,  only : DELTA07, USE_OPR_IN_STAR_OF_K, &
       &                          chg_symm_level1, chg_symm_level2
  use m_Crystal_Structure, only : nopr, op
! =============== 2014/09/30, 2018/02/26
! === KT_add ==== 13.2S
  use m_Crystal_Structure, only : sw_use_magnetic_symmetry
! =============== 13.2S

  use m_Control_Parameters,  only : sw_excitation, sw_band_unfolding, charge_symm_mode, &
       &                            sw_write_bxsf_file, use_metagga, vtau_exists
  use m_Crystal_Structure,  only : rltv_refcell, altv_refcell
  use m_IterationNumbers,   only : nkgroup
  use mpi

  implicit none
!  include 'mpif.h'

  integer, target                                      :: kv3=0 ! #sampling k-points
  real(kind=DP), allocatable, target, dimension(:,:,:) :: vkxyz ! d(kv3,3,CRDTYP), sampling k-points
  real(kind=DP), allocatable, target, dimension(:)     :: qwgt  ! d(kv3), weight factors
  integer, allocatable, dimension(:)           :: k_symmetry ! d(kv3), symmetries of k-points, Gamma or not
  integer, target                              :: kv3_ek ! k-points with charge-fixed
  integer                                      :: ek_group = 1 ! #k-point groups with charge-fixed
  real(kind=DP), allocatable, target, dimension(:,:,:) :: vkxyz_ek ! k-points
  real(kind=DP), allocatable, target, dimension(:)     :: qwgt_ek ! weight factors

! ========================= modified by K. Tagami ==================== 10.1
!  integer,          private,  dimension(3,2)   :: k_sample_mesh = 1
  integer,                     dimension(3,2)   :: k_sample_mesh = 1
! =================================================================== 10.1

  integer                                      :: kv3_previous = 0
  integer, private                             :: ipri_kp_count = 0, ipri_kp_t = 0
  integer                                      :: nkek  ! #kpoints for ek curve
  integer, private                             :: base_reduction_for_GAMMA = ON ! default=ON if possible
  integer, private                             :: base_symmetrization_for_GAMMA = ON ! default=ON if possible
  integer, parameter, private :: len_str = 132
  character(len=len_str),private :: str
  real(kind=DP),private,pointer,dimension(:,:) :: work

  integer                                      :: np0,np1,np2,nx1,ny1,nz1,nd
  integer,allocatable, dimension(:)            :: ip20,ip10,ip02,ip12
  integer,allocatable, dimension(:)            :: ip01,ip21,iu21,iv21,iwt,ip2cub
  integer,allocatable, dimension(:)            :: nxyz_tetra
  real(kind=DP), dimension(3,3)                :: trmat

  ! --- Ksampling ---
  character(len("ksampling")),private,parameter :: tag_ksampling = "ksampling"
  character(len("method")),private,parameter ::    tag_method    = "method"
  character(len("mesh")),private,parameter ::      tag_mesh      = "mesh"
  character(len("nx")),private,parameter ::        tag_nx        = "nx"
  character(len("ny")),private,parameter ::        tag_ny        = "ny"
  character(len("nz")),private,parameter ::        tag_nz        = "nz"
  character(len("monk")),private,parameter ::      tag_monkhorst_pack = "monk"
  character(len("mp_index")),private,parameter ::  tag_mp_index  = "mp_index"
  character(len("n1")),private,parameter ::        tag_n1        = "n1"
  character(len("n2")),private,parameter ::        tag_n2        = "n2"
  character(len("n3")),private,parameter ::        tag_n3        = "n3"
  character(len("kshift")),private,parameter ::    tag_kshift    = "kshift"
  character(len("k1")),private,parameter ::        tag_k1        = "k1"
  character(len("k2")),private,parameter ::        tag_k2        = "k2"
  character(len("k3")),private,parameter ::        tag_k3        = "k3"
  character(len("file")),private,parameter ::      tag_file      = "file"
  character(len("gamma")),private,parameter ::     tag_gamma     = "gamma"
  character(len("directin")),private,parameter ::  tag_directin  = "directin"
  character(len("num_kpoints")),private,parameter ::tag_num_kpoints = "num_kpoints"
  character(len("sum_weight")),private,parameter:: tag_sum_weight = "sum_weight"
  character(len("kpoints")),private,parameter ::   tag_kpoints   = "kpoints"
  character(len("kx")),private,parameter ::        tag_kx        = "kx"
  character(len("ky")),private,parameter ::        tag_ky        = "ky"
  character(len("kz")),private,parameter ::        tag_kz        = "kz"
  character(len("denom")),private,parameter ::     tag_denom     = "denom"
  character(len("weight")),private,parameter ::    tag_weight    = "weight"
  character(len("base_reduction_for_GAMMA")),private,parameter :: &
       &                       tag_base_reduction_for_GAMMA = "base_reduction_for_GAMMA"
  character(len("base_symmetrization_for_GAMMA")),private,parameter :: &
       &                       tag_base_symmetrization_for_G = "base_symmetrization_for_GAMMA"
  character(len("kv3")),private,parameter ::       tag_kv3       = "kv3"
  character(len("use_trs")),private,parameter ::  tag_use_trs  = "use_trs"
  character(len("density")),private,parameter :: tag_density = "density"
  ! --- temporary ---
  real(kind=DP), allocatable, dimension(:) :: kx_t,ky_t,kz_t
  real(kind=DP), allocatable, dimension(:) :: weight_t
  integer ::                                  kv3_t = -1

  real(kind=DP) :: k_density = 4.d0

  ! --- Monknorst-Pack method---
  !!integer, parameter, private :: ng=45,ngrid=(2*ng+1)**3,nshell=10000
  integer, private :: ng,ngrid,nshell
  integer, parameter, private :: nface_max = 125
  type shell
    integer :: ndim
    real(kind=DP) :: length
    real(kind=DP), pointer :: rindex(:,:) !(3,ndim)
  end type
  integer :: mp_index(3)
  real(kind=DP) :: kshift(3)

  integer :: itrs = ON

! ==================================== added by K. Tagami ==================== 12.0A
  character(len("gen_tetramesh_mode")), private,parameter :: &
       &                 tag_gen_tetramesh_mode    = "gen_tetramesh_mode"
  integer :: gen_tetramesh_mode = 0
!
!!  character(len("use_altv_rltv")), private,parameter :: &
 !!      &                 tag_use_altv_rltv    = "use_altv_rltv"
!!!  integer :: sw_use_altv_rltv = ON
!!!  logical :: use_altv_rltv = .true.
! ============================================================================ 12.0A

! ==================================== added by K. Tagami ==================== 11.0P
  character(len("use_op_before_sym_reduction")),private,parameter ::  &
       &         tag_use_op_before_sym_reduction = "use_op_before_sym_reduction"
  logical :: use_op_before_sym_reduction = .false.
  integer :: sw_use_op_before_sym_reduction = OFF
! ============================================================================ 11.0P

! === KT_add ==== 2014/09/30
  integer :: kv3_fbz
  integer, allocatable :: to_ibz_from_fbz_for_kpoint(:)
  real(kind=DP), allocatable :: vkxyz_fbz(:,:,:) ! d(kv3_fbz,3,CRDTYP)
!
  integer :: max_num_star_of_k
  integer, allocatable :: iopr_k_fbz_to_ibz(:)
  integer, allocatable :: trev_k_fbz_to_ibz(:)
  integer, allocatable :: num_star_of_k(:)
  integer, allocatable :: star_of_k(:,:)
! =============== 2014/09/30

! ==== ASMS ==== 2018/02/19
  integer :: nopr_from_fbz_to_ibz
  integer, allocatable :: flg_opr_from_fbz_to_ibz(:)
! ==== ASMS ==== 2018/02/19

  character(len("fix_ibz_on_fbz_mesh")),private,parameter :: &
       &  tag_fix_ibz_on_fbz_mesh  = "fix_ibz_on_fbz_mesh"
  integer :: sw_fix_ibz_on_fbz_mesh = OFF

! === unfolding ===
  real(kind=DP), allocatable :: vkxyz_refcell(:,:,:) ! d(kv3,3,CRDTYP)
! === unfolding ===

  character(len("sw_force_kpt_inside_bz")),private,parameter :: &
       &  tag_sw_force_kpt_inside_bz = "sw_force_kpt_inside_bz"
  integer :: sw_force_kpt_inside_bz = OFF
  integer, allocatable :: GvecTrans_kpt(:,:)
  integer, allocatable :: GvecTrans_kpt_ek(:,:)

! === Fast Monk-horst Pack
  character(len("segmentation")),private,parameter :: &
       &  tag_segmentation = "segmentation"
  character(len("sw_kpt_segmentation")),private,parameter :: &
       &  tag_sw_kpt_segmentation = "sw_kpt_segmentation"
  character(len("kpt_save_memory_mode")),private,parameter :: &
       &  tag_kpt_save_memory_mode = "kpt_save_memory_mode"
  integer :: sw_kpt_segmentation = OFF
  integer :: ndiv_segment(3)
  integer :: kpt_save_memory_mode = 1

contains

  subroutine alloc_kxyzweight(na)
    integer, intent(in) :: na
    allocate(kx_t(na)); allocate(ky_t(na));allocate(kz_t(na));allocate(weight_t(na))
  end subroutine alloc_kxyzweight

  subroutine dealloc_kxyzweight()
    if(allocated(kx_t)) deallocate(kx_t)
    if(allocated(ky_t)) deallocate(ky_t)
    if(allocated(kz_t)) deallocate(kz_t)
    if(allocated(weight_t)) deallocate(weight_t)
  end subroutine dealloc_kxyzweight

  subroutine m_Kp_sample_mesh(kmesh)
    integer, intent(out), dimension(3) :: kmesh
    kmesh(1:3) = k_sample_mesh(1:3,1)
  end subroutine m_Kp_sample_mesh

  subroutine m_Kp_wd_kv3(nfcntn)
    integer, intent(in) :: nfcntn
    if(mype==0) then
       write(nfcntn,'(a11)') tag_kv3
       write(nfcntn,'(i10)') kv3
    end if
  end subroutine m_Kp_wd_kv3

  subroutine m_Kp_rd_kv3(nfcntn,nfout)
    integer, intent(in) :: nfcntn,nfout
    logical :: EOF_reach, tag_is_found
    if(mype==0) then
       call rewind_to_tag0(nfcntn,len(tag_kv3),tag_kv3 &
            & , EOF_reach, tag_is_found, str,len_str)
       if(.not.tag_is_found) then
          write(nfout,'(" tag_kv3 is not found")')
       else
          read(nfcntn,*) kv3_previous
       end if
    end if
    if(npes > 1) call mpi_bcast(kv3_previous,1,mpi_integer,0,MPI_CommGroup,ierr)
    if(printable) write(nfout,'(i8, " : kv3_previous_job")') kv3_previous
  end subroutine m_Kp_rd_kv3

  subroutine m_Kp_rd_n(nfout)
    integer, intent(in) :: nfout
    character(len=FMAXVALLEN) :: rstr
    integer :: iret, f_selectBlock, f_getStringValue, f_getIntValue &
            & ,f_getRealValue, f_selectFirstTableLine
    integer :: f_selectParentBlock, f_selectTop
    integer :: sum_weight
    integer :: i,j
    real(kind=DP) :: dret
    logical :: prealloc = .false.
    real(kind=DP), dimension(3) :: rltv_len
    logical :: mesh_block_exists

    if(ipriinputfile >= 2) write(nfout,'(" !**  << m_Kp_rd_n >>")')
    ! --- accuracy ---
    iret = f_selectTop()
    if(icond > 1) then
       base_reduction_for_GAMMA = 0
       base_symmetrization_for_GAMMA = 0
    endif
    if(way_of_smearing == TETRAHEDRON) then
       call m_CtrlP_set_way_ksample(MESH)
    endif
    if( f_selectBlock( tag_accuracy) == 0) then
       if(ipriinputfile >= 2) write(nfout,'(" !** -- tag_accuracy --")')
       if( f_selectBlock(tag_ksampling) == 0) then
         if(f_getRealValue(tag_density,dret,'bohr') == 0) then
            k_density = dret
            write(nfout, '(a,f10.5)') ' !** density of the sampling kpoints',k_density
         endif
         iret = f_selectParentBlock()
       endif
       if(k_density>0)then
          do i=1,3
             rltv_len(i) = 0.d0
             do j=1,3
                rltv_len(i) = rltv_len(i)+rltv(i,j)*rltv(i,j)
             enddo
             rltv_len(i) = dsqrt(rltv_len(i))
             k_sample_mesh(i,1) = int(rltv_len(i)*k_density)
             if(k_sample_mesh(i,1)<1) k_sample_mesh(i,1) = 1
             k_sample_mesh(i,2) = k_sample_mesh(i,1)
          enddo
          mesh_block_exists = .false.
          if( f_selectBlock( tag_ksampling) == 0) then
             if(f_selectBlock(tag_mesh) ==0)then
                iret = f_selectParentBlock()
                mesh_block_exists = .true.
             endif
             iret = f_selectParentBlock()
          endif
          if(.not.mesh_block_exists)then
            if(ipriinputfile>=1) write(nfout,'(a,f10.5,3i8)') &
           & ' !** kmesh density and sampling kmesh resolved from it : ' &
           & , k_density, k_sample_mesh(1:3,1)
          endif
          if(way_ksample == MONKHORST_PACK)then
            do i=1,3
               mp_index(i) = k_sample_mesh(i,1)
            enddo
            kshift(1:3) = 0.5d0   ! default value for cubic system
            if(il == 0) then   ! default value for hexagonal system
               kshift(1:2) = 0.0d0;     kshift(3) = 0.5d0
            end if
          endif
       else
          if(mype==0) write(nfout,'(a,f10.5)') ' !** WARNING : invalid density specification ',k_density
       endif
       if( f_selectBlock( tag_ksampling) == 0) then

          if(icond > 1) base_reduction_for_GAMMA = 0
          if( f_getIntValue(tag_base_reduction_for_GAMMA, iret)==0) then
             base_reduction_for_GAMMA = iret
             if(ipri_kp >= 1) then
                write(nfout,'(" !** base_reduction_for_GAMMA = ",i3 &
                     & ," in tag_ksampling <<m_Kp_rd_n>>")') base_reduction_for_GAMMA
             end if
             if(base_reduction_for_GAMMA < 0 .or. base_reduction_for_GAMMA >1) then
                base_reduction_for_GAMMA = 1
                if(ipri_kp >= 1) write(nfout,'(" !** base_reduction_for_GAMMA = ",i3 &
                     & ," (reset) in tag_ksampling <<m_Kp_rd_n>>")') base_reduction_for_GAMMA
             end if
          else
             if(ipri_kp >= 1) then
                write(nfout,'(" !** base_reduction_for_GAMMA = ",i3 &
                     & ," (=default value) <<m_Kp_rd_n>>")') base_reduction_for_GAMMA
             end if
          end if

          if(icond > 1) base_symmetrization_for_GAMMA = 0
          if( f_getIntValue(tag_base_symmetrization_for_G, iret)==0) then
             base_symmetrization_for_GAMMA = iret
             if(ipri_kp >= 1) then
                write(nfout,'(" !** base_symmetrization_for_GAMMA = ",i3 &
                     & ," in tag_ksampling <<m_Kp_rd_n>>")') base_symmetrization_for_GAMMA
             end if
             if(base_symmetrization_for_GAMMA < 0 .or. base_symmetrization_for_GAMMA >1) then
                base_symmetrization_for_GAMMA = 1
                write(nfout,'(" !** base_symmetrization_for_GAMMA = ",i3 &
                     & ," (reset) in tag_ksampling <<m_Kp_rd_n>>")') base_symmetrization_for_GAMMA
             end if
          else
             if(ipri_kp >= 1) then
                write(nfout,'(" !** base_symmetrization_for_GAMMA = ",i3 &
                     & ," (=default value) <<m_Kp_rd_n>>")') base_symmetrization_for_GAMMA
             end if
          end if

          if ( sw_excitation == ON .or. ( use_metagga .and. vtau_exists ) ) then
             base_reduction_for_GAMMA = off
             base_symmetrization_for_GAMMA = off
             write(nfout,*) "!** base_reduction_for_GAMMA is forced to OFF"
             write(nfout,*) "!** base_symmetrization_for_GAMMA is forced to OFF"
          endif

          if( f_getStringValue( tag_method, rstr, LOWER) == 0) then
             call set_ksamplingmethod(rstr)
             if(way_ksample == MONKHORST_PACK .and. way_of_smearing == TETRAHEDRON) then
               !if(printable) write(nfout,'(a)') ' !** tetrahedral smearing cannot be used in &
               !                       & conjunction with monk sampling. specify mesh instead.'
               call phase_error_with_msg(nfout,' !** tetrahedral smearing cannot be used in &
                                      & conjunction with monk sampling. specify mesh instead.'&
                                      ,__LINE__,__FILE__)
             endif
          else
             if(ipri_kp >= 1) then
                write(nfout,'(" !** a method of sampling is not defined in the accuracy block <<m_Kp_rd_n>>")')
             end if
          end if
          if( way_ksample == MESH ) then
             if( f_selectBlock( tag_mesh) == 0) then
                if( f_getIntValue( tag_nx, iret) == 0) k_sample_mesh(1,1) = iret
                if( f_getIntValue( tag_ny, iret) == 0) k_sample_mesh(2,1) = iret
                if( f_getIntValue( tag_nz, iret) == 0) k_sample_mesh(3,1) = iret
                k_sample_mesh(:,2) = k_sample_mesh(:,1)
                iret = f_selectParentBlock()
             else
!!                k_sample_mesh(1:3,1) = 4   ! default value
             end if

! ================================ added by K. Tagami ======================== 12.0A
!             if( f_getIntValue( tag_use_altv_rltv, iret) == 0 ) then
!                sw_use_altv_rltv = iret
!             endif
!             if ( sw_use_altv_rltv == ON ) then
!                use_altv_rltv = .true.
!             else
!                use_altv_rltv = .false.
!             endif             ! ---- judgement is moved to m_Crystal_Structre.F90

             if( f_getIntValue( tag_gen_tetramesh_mode, iret) == 0 ) then
                if ( iret == 0 .or. iret == 1  ) then
                   gen_tetramesh_mode = iret
                endif
             endif

! === KT 2018/03/07
             if ( use_altv_rltv ) gen_tetramesh_mode = 1
! === KT 2018/03/07

             write(nfout,*) '!************************************ '
             write(nfout,*) '!** gen_tetramesh_mode is set to ', gen_tetramesh_mode
!!             write(nfout,*) '!** use_altv_rltv is set to ', use_altv_rltv
             write(nfout,*) '!************************************ '

!
             if ( noncol ) then
                write(*,*) '***** Caution *****'
                write(*,*) 'In noncollinear systems "method = mesh" is experimental '
             endif
!
             if( f_getIntValue( tag_use_trs, iret) == 0) then
                itrs = iret
             else
                itrs = OFF       ! for consistency with previous programs
             endif
! ----------
             if ( noncol ) then
                if ( sw_neglect_magmom == OFF ) itrs = OFF
!                itrs = OFF
             else
                if ( imag /= PARA ) itrs = OFF
             endif
! ----------
! ================================================================== 12.0A
             if ( sw_write_bxsf_file == ON ) then
                if ( nopr > 1 ) call phase_error_with_msg(nfout," kpoint = mesh is not supported when nopr > 1",__LINE__,__FILE__)
             endif

           else if( way_ksample == MONKHORST_PACK ) then
              if( f_selectBlock( tag_mp_index) == 0) then
                 if( f_getIntValue( tag_n1, iret) == 0) mp_index(1) = iret
                 if( f_getIntValue( tag_n2, iret) == 0) mp_index(2) = iret
                 if( f_getIntValue( tag_n3, iret) == 0) mp_index(3) = iret
                 iret = f_selectParentBlock()
              else if( f_selectBlock( tag_mesh) == 0) then
                 if( f_getIntValue( tag_nx, iret) == 0) mp_index(1) = iret
                 if( f_getIntValue( tag_ny, iret) == 0) mp_index(2) = iret
                 if( f_getIntValue( tag_nz, iret) == 0) mp_index(3) = iret
                 iret = f_selectParentBlock()
              else
                 !!$ mp_index(1:3) = 4   ! default value
                 continue
              end if

              kshift(1:3) = 0.5d0   ! default value for cubic system
              if(il == 0) then   ! default value for hexagonal system
                 kshift(1:2) = 0.0d0;     kshift(3) = 0.5d0
              end if
              if( f_selectBlock( tag_kshift) == 0) then
                 if( f_getRealValue( tag_k1, dret, '' ) == 0) kshift(1) = dret
                 if( f_getRealValue( tag_k2, dret, '' ) == 0) kshift(2) = dret
                 if( f_getRealValue( tag_k3, dret, '' ) == 0) kshift(3) = dret
                 iret = f_selectParentBlock()
              end if

              if ( sw_excitation == ON .or. sw_write_bxsf_file == ON ) then
                 kshift = 0.0d0
                 write(nfout,*) "!** kshift is forced to 0.0"
              endif

              if( f_getIntValue( tag_use_trs, iret) == 0) then
                 itrs = iret
              else
                 if ( noncol ) itrs = OFF
              endif
! ----------
              if ( noncol ) then
                 if ( sw_neglect_magmom == OFF ) itrs = OFF
!                 itrs = OFF
              else
                 if ( imag /= PARA ) itrs = OFF
              endif
! ----------
! ==== ASMS ===
              if( f_getIntValue( tag_fix_ibz_on_fbz_mesh, iret) == 0) then
                 sw_fix_ibz_on_fbz_mesh = iret
                 write(nfout,*) "!** sw_fix_ibz_on_fbz_mesh is set to ", iret
              endif
! ==== ASMS ===

              if( f_selectBlock( tag_segmentation) == 0) then
                 if( f_getIntValue( tag_sw_kpt_segmentation, iret) == 0 ) then
                    sw_kpt_segmentation = iret
                    write(nfout,*) "!** sw_kpt_segmentation is set to ", iret
                 endif
                 if( f_getIntValue( tag_kpt_save_memory_mode, iret) == 0 ) then
                    kpt_save_memory_mode = iret
                    write(nfout,*) "!** kpt_save_memory_mode set to ", iret
                 endif
                 ndiv_segment(1:3) = mp_index(1:3)
                 if( f_getRealValue( tag_nx, dret, '' ) == 0) ndiv_segment(1) = dret
                 if( f_getRealValue( tag_ny, dret, '' ) == 0) ndiv_segment(2) = dret
                 if( f_getRealValue( tag_nz, dret, '' ) == 0) ndiv_segment(3) = dret
                 iret = f_selectParentBlock()
              end if

              if( f_getIntValue( tag_use_op_before_sym_reduction, iret) == 0 ) then
                 sw_use_op_before_sym_reduction = iret
              endif
              if ( sw_use_op_before_sym_reduction == ON ) then
                 use_op_before_sym_reduction = .true.
              else
                 use_op_before_sym_reduction = .false.
              endif
              if(ipri_kp >= 1) then
                 write(nfout,'("<< Monkhorst-Pack scheme")')
                 write(nfout,'(" MP index:",3(1x,i3))') mp_index(1:3)
                 write(nfout,'(" kp shift:",3(1x,f10.5))') kshift(1:3)
                 write(nfout,'(" il = ",i3)') il
                 write(nfout,'(" use_trs = ",i3)') itrs

                 if ( noncol ) then
                    write(nfout,*) "use_op_before_sym_reduction = ", &
                         &            use_op_before_sym_reduction
                 endif

                 write(nfout,'("   Monkhorst-Pack scheme >>")')
              end if

              do i=1,3
                if(mp_index(i) > 100) then
                   if(printable) then
                      write(nfout,*)  &
                           & 'Monkhorst-Pack: found a bad parameter'
                      write(nfout,*)  &
                           & '#### Use Monknorst-Pack indeces less than 101. ###'
                   end if
                   call phase_error_with_msg(nfout,'Monkhorst-Pack: found a bad parameter',__LINE__,__FILE__)
                end if
              end do

           else if( way_ksample == FILE ) then
              if( f_getIntValue( tag_sw_force_kpt_inside_bz, iret) == 0) then
                 sw_force_kpt_inside_bz = iret
                 write(nfout,*) "!** sw_force_kpt_inside_bz is set to ", iret
              endif

           else if( way_ksample == SKPS_DIRECT_IN ) then
             sum_weight = 0
             if( f_getIntValue( tag_sum_weight,  iret) == 0) sum_weight = iret
             if( f_getIntValue( tag_num_kpoints, iret) == 0) then
                kv3_t = iret
                if(kv3_t <= 0) call phase_error_with_msg(nfout,' kv3 is not positive value << m_Kp_rd_n >>',__LINE__,__FILE__)
             else
                prealloc = .true.
                call set_kxyz(prealloc,iret)
                kv3_t = iret
                if(ipriinputfile >= 3) write(nfout,'(" !* kv3_t = ",i6," <<m_Kp_rd_n>>")') kv3_t
             end if
             prealloc = .false.
             call alloc_kxyzweight(kv3_t)
             call set_kxyz(prealloc,iret)
          end if
          iret = f_selectParentBlock()
       end if
       iret = f_selectParentBlock()
    end if
    if(way_ksample == MESH) then
!!$       if(ipriinputfile >= 1) write(nfout,'(" !** k-point sampling method = ",i6)') way_ksample
       if(ipriinputfile >= 1) write(nfout,'(" !** k_sample_mesh = ",3i6)') k_sample_mesh(1:3,1)
    end if
  contains
    subroutine set_kxyz(prealloc,iret)
      logical, intent(in) ::  prealloc
      integer, intent(out)::  iret
      integer :: i, f_readKPoints, f_selectFirstTableLine, f_selectNextTableLine
      real(kind=DP),dimension(3) :: kvec
      integer ::              weight, sum
      if( f_selectBlock(tag_kpoints) == 0) then
         sum = 0
         i = 1
         do while(.true.)
            if( i == 1 ) then
               if( f_selectFirstTableLine() /= 0 ) then
                  exit
               end if
            else
               if( f_selectNextTableLine() /= 0 ) then
                  exit
               end if
            end if
            iret = f_readKPoints( tag_kx,tag_ky,tag_kz,tag_denom,tag_weight &
                 & ,kvec, weight )
!!$            print '(f8.4, f8.4, f8.4, I3)', kvec(1),  kvec(2), kvec(3), weight
            if(.not.prealloc) then
               if(i > kv3_t) exit
               if(weight < 1.d-13) weight = 1
               sum = sum + weight
               kx_t(i) = kvec(1); ky_t(i) = kvec(2); kz_t(i) = kvec(3)
               weight_t(i) = weight
            end if
            i = i+1
         end do
         if(ipriinputfile >= 3) write(nfout,'(" !*  weight-sum = ", i6)') sum
         if(.not.prealloc .and. sum_weight /= 0) then
            if(sum_weight /= sum) then
               if(printable) write(nfout,'(" !* Given sum_weight is not equal to summed weight")')
            end if
         end if
         if(.not.prealloc) weight_t = weight_t/dble(sum)
         iret = f_selectParentBlock()
      else
         call phase_error_with_msg(nfout,' ! No kpoints is given in the inputfile <<m_Kp_rd_n>>',__LINE__,__FILE__)
      end if
      if(prealloc) iret = i-1
      if(.not.prealloc) then
         if(ipriinputfile >= 2) then
            do i = 1, kv3_t
               write(nfout,'(" !* i = ",i6," kvxyz,weight = ",4f8.4)') i,kx_t(i),ky_t(i),kz_t(i),weight_t(i)
            end do
         end if
      end if
    end subroutine set_kxyz

    subroutine set_ksamplingmethod(rstr)
      character(len=FMAXVALLEN),intent(in) :: rstr
      logical :: tf
      integer :: ksamp
      ksamp = MONKHORST_PACK
      if(way_of_smearing == TETRAHEDRON) then
         call m_CtrlP_set_way_ksample(MESH)
      endif
      call strncmp0(tag_monkhorst_pack, trim(rstr), tf)
      if(tf) then
         ksamp = MONKHORST_PACK
         goto 1001
      end if
      call strncmp0(tag_mesh, trim(rstr), tf)
      if(tf) then
         ksamp = MESH
         goto 1001
      end if
      call strncmp0(tag_file, trim(rstr), tf)
      if(tf) then
         ksamp = FILE
         goto 1001
      end if
      call strncmp0(tag_gamma,trim(rstr),tf)
      if(tf) then
         ksamp = GAMMA
         goto 1001
      end if
      call strncmp0(tag_directin, trim(rstr),tf)
      if(tf) then
         ksamp = SKPS_DIRECT_IN
         goto 1001
      end if
      call phase_error_with_msg(nfout,' ! tag for ksampling is invalid <<m_Kp_rd_n.set_ksamplingmethod>>',__LINE__,__FILE__)
1001  continue
      if(ipri_kp>= 2) write(nfout,'(" !** ksampling method = ",a14)') trim(rstr)
      call m_CtrlP_set_way_ksample(ksamp)
    end subroutine set_ksamplingmethod
  end subroutine m_Kp_rd_n

  subroutine m_Kp_alloc_kpoints()
    if(.not.allocated(vkxyz)) then
       allocate(vkxyz(kv3,3,CRDTYP)); vkxyz = 0.d0
       if(ipri_kp >= 2) then
          write(nfout,'(" !kp vkxyz is allocated now <<m_Kp_alloc_kpoints>>")')
          write(nfout,'(" !kp kv3 = ", i8," CRDTYP = ",i3)') kv3,CRDTYP
       end if
    else
       if(ipri_kp >= 2) write(nfout,'(" !kp vkxyz is already allocated <<m_Kp_alloc_kpoints>>")')
    end if
    if(.not.allocated(qwgt)) then
       allocate(qwgt(kv3))
       if(ipri_kp >= 2) write(nfout,'(" !kp qwgt is allocated now <<m_Kp_alloc_kpoints>>")')
       qwgt = 1.d0
    else
       if(ipri_kp >= 2) write(nfout,'(" !kp qwgt is already allocated <<m_Kp_alloc_kpoints>>")')
    end if
    if(.not.allocated(k_symmetry)) then
       if(ekmode==OFF .or. kv3_ek==0) then
           allocate(k_symmetry(kv3)); k_symmetry = 0
       else
           allocate(k_symmetry(kv3_ek)); k_symmetry = 0
       endif
       if(ipri_kp >= 2) then
          write(nfout,'(" !kp k_symmetry is allocated now <<m_Kp_alloc_symmetry>>")')
       end if
    else
       if(ipri_kp >= 2) write(nfout,'(" !kp k_symmetry is already allocated <<m_Kp_alloc_kpoints>>")')
    end if

! === KT_add === 2014/09/30
!    write(*,*) "kv3_fbz = ", kv3_fbz, kv3
    if ( .not. allocated(vkxyz_fbz) ) allocate(vkxyz_fbz(kv3_fbz,3,CRDTYP))
    if ( .not. allocated(to_ibz_from_fbz_for_kpoint) ) then
       allocate( to_ibz_from_fbz_for_kpoint(kv3_fbz) )
    endif
! ============== 2014/09/30
    if ( sw_band_unfolding == ON ) then
       if ( .not. allocated(vkxyz_refcell) ) allocate(vkxyz_refcell(kv3,3,CRDTYP) )
    endif
    if ( sw_force_kpt_inside_bz == ON ) then
       if ( .not. allocated(GvecTrans_kpt) ) allocate(GvecTrans_kpt(kv3,3) )
    endif

  end subroutine m_Kp_alloc_kpoints

  subroutine m_Kp_alloc_kpoints_ek()
    kv3_ek = kv3
    if(.not.allocated(vkxyz_ek)) allocate(vkxyz_ek(kv3_ek,3,CRDTYP))
    if(.not.allocated(qwgt_ek))  allocate(qwgt_ek(kv3_ek))
    if ( sw_force_kpt_inside_bz == ON ) then
       if ( .not. allocated(GvecTrans_kpt_ek) ) then
          allocate(GvecTrans_kpt_ek(kv3_ek,3) )
       endif
    endif
    if(allocated(k_symmetry)) then
      deallocate(k_symmetry)
    endif
    allocate(k_symmetry(kv3_ek));k_symmetry = 0
  end subroutine m_Kp_alloc_kpoints_ek

  subroutine m_Kp_cp_vkxyz_to_vkxyz_ek()
    vkxyz_ek = vkxyz
    qwgt_ek  = qwgt
    if ( sw_force_kpt_inside_bz == ON ) then
       GvecTrans_kpt_ek = GvecTrans_kpt
    endif
  end subroutine m_Kp_cp_vkxyz_to_vkxyz_ek

  subroutine m_Kp_cp_vkxyz_ek_to_vkxyz(nk)
    integer, intent(in) :: nk
    integer :: i, kvt, nks, is, kv3t, ikt
    if(ipri_kp >= 1) write(nfout,'(" !kp nk, kv3, kv3_ek = ",3i6)') nk, kv3, kv3_ek
    if(nk+kv3-1 > kv3_ek) then
       if(fixed_charge_k_parallel /= ONE_BY_ONE) &
            & call phase_error_with_msg(nfout,' nk is illegal (m_Kp_cp_vkxyz_ek_to_vkxyz)',__LINE__,__FILE__)
    end if

    if(sw_modified_kpoint_increment == ON)then
       do i=0,nrank_k-1
          do is=1,nspin
             ikt = nspin*(nis_kv3_ek(i)-1)+(nkgroup-1)*nspin+is
             if(ikt>kv3_ek) ikt = kv3_ek
             vkxyz(i*nspin+is,1:3,1:CRDTYP) = vkxyz_ek(ikt,1:3,1:CRDTYP)
             qwgt(i*nspin+is) = qwgt_ek(ikt)
          enddo
       enddo
    else
       if(fixed_charge_k_parallel == ONE_BY_ONE .and. nk+kv3-1 > kv3_ek) then
          kvt = kv3_ek - nk + 1
          vkxyz(1:kvt,1:3,1:CRDTYP) = vkxyz_ek(nk:kv3_ek,1:3,1:CRDTYP)
          qwgt(1:kvt) = qwgt_ek(nk:kv3_ek)
          
          if ( sw_force_kpt_inside_bz == ON ) then
             GVecTrans_kpt(1:kvt,1:3) = GvecTrans_kpt_ek(nk:kv3_ek,1:3)
          endif
          do i = kvt+1, kv3
             vkxyz(i,1:3,1:CRDTYP) = vkxyz_ek(kv3_ek,1:3,1:CRDTYP)
             qwgt(i) = qwgt_ek(kv3_ek)
             if ( sw_force_kpt_inside_bz == ON ) then
                GVecTrans_kpt(i,1:3) = GvecTrans_kpt_ek(kv3_ek,1:3)
             endif
          end do
       else
          vkxyz(1:kv3,1:3,1:CRDTYP) = vkxyz_ek(nk:nk+kv3-1,1:3,1:CRDTYP)
          qwgt(1:kv3) = qwgt_ek(nk:nk+kv3-1)
          if ( sw_force_kpt_inside_bz == ON ) then
             GVecTrans_kpt(i,1:3) = GvecTrans_kpt_ek(kv3_ek,1:3)
          endif
       endif
    end if

    if(kimg == 2) then
       if(base_reduction_for_GAMMA == ON) then
          do i = 1, kv3
             if(vkxyz(i,1,1)**2 + vkxyz(i,2,1)**2 + vkxyz(i,3,1)**2 < DELTA) then
                k_symmetry(i) = GAMMA
             else
                k_symmetry(i) = 0
             end if
          end do
       else if(base_reduction_for_GAMMA == OFF .and. base_symmetrization_for_GAMMA == ON) then
          do i = 1, kv3
             if(vkxyz(i,1,1)**2 + vkxyz(i,2,1)**2 + vkxyz(i,3,1)**2 < DELTA) then
                k_symmetry(i) = GAMMA_base_symmetrization
             else
                k_symmetry(i) = 0
             end if
          end do
       end if
    else
       do i = 1, kv3
          k_symmetry(i) = 0
       end do
    end if

! ================================== added by K. Tagami =============== 11.0
    if ( noncol ) then
       write(nfout,*) '!----- '
       write(nfout,*) '! k_symmetry is forced to 0 in the non-collinear system'
       write(nfout,*) '!----- '
       k_symmetry(:) = 0
    endif
! ===================================================================== 11.0

    if(ipri_kp >= 1) then
       do i = 1, kv3
          if(k_symmetry(i) == GAMMA) then
             write(nfout,'(" !kp ",i6, 6f8.4, " GAMMA")') &
                  & i, vkxyz(i,1:3,CARTS), vkxyz(i,1:3,BUCS)
          else
             write(nfout,'(" !kp ",i6, 6f8.4)') &
                  & i, vkxyz(i,1:3,CARTS), vkxyz(i,1:3,BUCS)
          end if
       end do
    end if
  end subroutine m_Kp_cp_vkxyz_ek_to_vkxyz

  subroutine m_Kp_set_ek_group()
! ============================== modified by K. Tagami =========== 11.0
!    if(nspin == 1) then
!       ek_group = ceiling(dble(kv3_ek)/nrank_k)
!    else
!       ek_group = ceiling(dble(kv3_ek/nspin)/nrank_k)*nspin
!    end if
!    kv3 = nrank_k*nspin
!
    if ( noncol ) then
      ek_group = ceiling(dble(kv3_ek/ndim_spinor)/nrank_k)*ndim_spinor
      kv3 = nrank_k *ndim_spinor
    else
      if(nspin == 1) then
         ek_group = ceiling(dble(kv3_ek)/nrank_k)
      else
         ek_group = ceiling(dble(kv3_ek/nspin)/nrank_k)*nspin
      end if
      kv3 = nrank_k *nspin
    endif
! ================================================================ 11.0

    if(ipri_kp >=1 ) then
       write(nfout,'(" !kp ek_group = ",i8)') ek_group
       write(nfout,'(" !kp kv3_ek   = ",i8)') kv3_ek
       write(nfout,'(" !kp kv3      = ",i8)') kv3
    end if
  end subroutine m_Kp_set_ek_group

  subroutine m_Kp_realloc_kpoints2()
    if(allocated(vkxyz)) deallocate(vkxyz)
    if(allocated(qwgt))  deallocate(qwgt)
    if(allocated(k_symmetry)) deallocate(k_symmetry)
    call m_Kp_alloc_kpoints()
  end subroutine m_Kp_realloc_kpoints2

  subroutine m_Kp_realloc_kpoints()
! ================================ modified by K. Tagami ========== 11.0
!    kv3 = nspin
!
    if ( noncol ) then
       kv3 = ndim_spinor
    else
       kv3 = nspin
    endif
! ================================================================= 11.0

    if(allocated(vkxyz)) deallocate(vkxyz)
    if(allocated(qwgt))  deallocate(qwgt)
    if(allocated(k_symmetry)) deallocate(k_symmetry)

! === KT_add === 2014/09/30
    if ( allocated(vkxyz_fbz) ) deallocate(vkxyz_fbz)
    if ( allocated(to_ibz_from_fbz_for_kpoint) ) deallocate(to_ibz_from_fbz_for_kpoint)
! ============== 2014/09/30
    if ( sw_force_kpt_inside_bz == ON ) then
       if ( allocated(GvecTrans_kpt) ) deallocate(GvecTrans_kpt)
    endif
    if ( sw_band_unfolding == ON ) then
       if ( allocated(vkxyz_refcell) ) deallocate(vkxyz_refcell)
    endif

    call m_Kp_alloc_kpoints()

  end subroutine m_Kp_realloc_kpoints

  subroutine m_Kp_rd_way_ksample_etc(nfinp,nlines)
    integer, intent(in) :: nfinp,nlines
    integer, parameter  :: NWK = 6
    integer             :: natm0, ntyp0
    allocate(work(NWK,1))

    if(npes >= 2) call mpi_barrier(MPI_CommGroup,ierr)
    rewind nfinp
    if(ipri_kp >= 2) write(nfout,'(" <<< m_Kp_rd_way_ksample_etc >>>")')
!       ==========================
    call read_natm_ntyp          ! -> natm0, ntyp0
    if(ipri_kp >= 2) then
       write(nfout,'(" !kp nlines+14+natm0+ntyp0 = ",i8)') nlines+14+natm0+ntyp0
       write(nfout,'(" !kp ntyp0, natm0 = ",2i8)') ntyp0, natm0
    end if
    call skip_lines(nfinp,nlines+14+natm0+ntyp0)
    call read_way_ksample_etc    ! -> way_ksample, k_sample_mesh
!       ==========================
    deallocate(work)
  contains
    subroutine read_natm_ntyp
      integer nsize
      read(nfinp,'(a132)') str
      call chnnm(str,len_str,4,work,nsize)
      ntyp0 = nint(work(3,1))
      natm0 = nint(work(4,1))
    end subroutine read_natm_ntyp

    subroutine read_way_ksample_etc
      integer :: i, ikpnt, sw_k_coord_system
      integer :: ksamp
      read(nfinp,'(a132)') str
      if(ipri_kp >= 2) write(nfout,'(a132)') str
      call ch4ksample(str,len_str,ksamp)
      if(ksamp == NONAME .or. ksamp == MESH) then
         if(ksamp == NONAME) then
            read(str,*) (k_sample_mesh(i,1), i=1,3)
            ksamp = MESH
         else
            read(nfinp,*) (k_sample_mesh(i,1),i=1,3)
         endif
         if(ipri_kp>=1) write(nfout,680) (k_sample_mesh(i,1),i=1,3)
         read(nfinp,*) (k_sample_mesh(i,2),i=1,3)
         if(ipri_kp>=1) write(nfout,690) (k_sample_mesh(i,2),i=1,3)
         if(k_sample_mesh(1,1) == 0  &
              & .or.(     k_sample_mesh(1,2) == 0 &
              &     .and. k_sample_mesh(2,2) == 0 &
              &     .and. k_sample_mesh(3,2) == 0)) &
              &  ksamp = GAMMA
      else if(ksamp.eq.SKPS_DIRECT_IN) then
         if(ipri_kp>=1) write(nfout,'(" ! way_ksample = SKPS_DIRECT_IN")')
         read(nfinp,'(a132)') str
         if(ipri_kp>=2) write(nfout,'(" ! str : ",a50)') str(1:50)
         call read_nkpnt(str,len_str,ikpnt) !-(b_Kpoints)
         if(ipri_kp>=2) write(nfout,'(" ! ikpnt(#skipped lines) = ",i8)') ikpnt
         read(nfinp,'(a132)') str
         call read_coordsystem(str,len_str,sw_k_coord_system) !-(b_Kpoints)
         if(sw_k_coord_system == NODATA) ikpnt = ikpnt - 1
         call skip_lines(nfinp,ikpnt)      !-(bottom_Subroutines)
!         call read_konly
!     >        (nfinp,nfout,knv3,altv,rltv,nspin
!     <        ,kv3,vkxyz,qwgt,wka)
      else if(ksamp == FILE) then
         if(ipri_kp>=1) write(nfout,'(" ! way_ksample = FILE")')
         call skip_lines(nfinp,2)
      else if(ksamp == MONKHORST_PACK) then
         if(ipri_kp>=1) write(nfout,'(" ! way_ksample = MONKHORST_PACK")')
         call read_mpindex(nfinp) !-(b_SpeicalKpints)
      endif

      call m_CtrlP_set_way_ksample(ksamp)
680   format(' ',3i6,'   : nkx, nky, nkz ')
690   format(' ',3i6,'   : nkx2,nky2,nkz2')
    end subroutine read_way_ksample_etc
  end subroutine m_Kp_rd_way_ksample_etc

  subroutine m_Kp_gnrt_or_rd_k_points(nfinp,preallocation)
    ! -> kv3(when preallocation == .true.), kv3,vkxyz,qwgt(when preallocation == .flase.)
    external kpmsf0, kpmwbz0, sccm0, bccm0, fccm0, hexm0
    integer, intent(in)      :: nfinp
    logical, intent(in)      :: preallocation
    integer                  :: knv3_t, id_sname = -1, i, kv3_adj, j
    integer :: knv3_t_fbz
    call tstatc0_begin('m_Kp_gnrt_or_rd_k_points ', id_sname)

    if(ipri_kp >= 2) then
       write(nfout,'(" !kpgnrt way_ksample = ",i6)') way_ksample
       write(nfout,'(" !kpgnrt nbztyp      = ",i6)') nbztyp
       write(nfout,'(" !kpgnrt nbztyp_spg  = ",i6)') nbztyp_spg
    end if

    if(preallocation) then
       allocate(vkxyz(1,1,1))
       allocate(qwgt(1))
       allocate(vkxyz_fbz(1,1,1))
       allocate(to_ibz_from_fbz_for_kpoint(1))
    end if

    knv3_t = kv3;  knv3_t_fbz = kv3_fbz

    if(sw_berry_phase == 1) then
       continue

    else if(way_ksample == GAMMA) then
       call gen_gamma_p0(preallocation,knv3_t,kv3,vkxyz,qwgt)  !-(b_Kpoints)

       if ( preallocation ) then
          kv3_fbz = kv3
       else
          kv3_fbz = kv3_fbz /nspin
          vkxyz_fbz = vkxyz
          to_ibz_from_fbz_for_kpoint = 1
       endif

    else if(way_ksample == FILE) then
       call m_Files_open_kpoint_files(way_ksample,nbztyp_spg)
       if ( sw_band_unfolding == OFF ) then
          call readk0(preallocation,nfout,ipri_kp,knv3_t,rltv,nfkpoint,nfmatbp&
               &,kv3,vkxyz,qwgt) !-(b_Kpoints)
       else
          call readk0(preallocation,nfout,ipri_kp,knv3_t,rltv_refcell,nfkpoint,nfmatbp&
               &,kv3,vkxyz_refcell,qwgt) !-(b_Kpoints)
          call readk0_for_band_unfolding( preallocation, nfout, ipri_kp, knv3_t, &
               &                          rltv_refcell, altv,&
               &                          nfkpoint, nfmatbp, kv3, vkxyz, qwgt )
       endif

       if(ipriinputfile >= 2 .and. .not.preallocation) then
          do i=1, kv3
             write(nfout,'(" !** id, vkxyz, qwgt = ",i6,4f8.4," << m_Kp_gnrt_or_rd_k_points>>")') &
                  & i,vkxyz(i,1,CARTS),vkxyz(i,2,CARTS),vkxyz(i,3,CARTS),qwgt(i)
          end do
       end if
    else if(way_ksample == SKPS_DIRECT_IN) then
       if(ipri_kp>=2) write(nfout,'(" !kp kv3_t = ", i6)') kv3_t
       if(kv3_t > 0) then  ! directin from inputfile in new format style
          kv3 = kv3_t
          if(preallocation) then
             if(ipriinputfile >= 2) then
                do i=1, kv3_t
                   write(nfout,'(" !** id,kx_t,ky_t,kz_t,weight_t = " &
                        &    ,i6,4f8.4," <<m_Kp_gnrt_or_rd_k_points>>")') &
                        & i,kx_t(i),ky_t(i),kz_t(i),weight_t(i)
                end do
             end if
          else
             kv3_adj = ubound(vkxyz,1)-lbound(vkxyz,1) + 1
             call set_kfromtemporaryarray(nfout,ipri,nfmatbp,kv3_adj,kv3_t,kx_t,ky_t,kz_t,weight_t &
                  &   ,rltv,vkxyz,qwgt)
!!$             call dealloc_kxyzweight()
             if(ipriinputfile >= 2) then
                do i=1, kv3
                   write(nfout,'(" !** id, vkxyz, qwgt = ",i6,4f8.4," << m_Kp_gnrt_or_rd_k_points>>")') &
                        & i,vkxyz(i,1,CARTS),vkxyz(i,2,CARTS),vkxyz(i,3,CARTS),qwgt(i)
                end do
             end if
          end if
       else
          call read_konly(preallocation,nfinp,nfout,knv3_t &
               & ,altv,rltv,kv3,vkxyz,qwgt,str,len_str)
       end if

    else if(way_ksample == MONKHORST_PACK) then

       call gen_special_kpoints( preallocation, nfout, ipri, knv3_t, &
            &                    nfmatbp, kv3, vkxyz, qwgt, imag, itrs, &
            &                    kv3_fbz, vkxyz_fbz, to_ibz_from_fbz_for_kpoint )

    else if( nbztyp == GENERAL .or. nbztyp == GENERAL_LARGER) then
       ipri_kp_t = ipri_kp*(1-ipri_kp_count)
       if(ipri_kp_count == 0) ipri_kp_count = 1


! ======================================= modified by K. Tagami =========== 12.0A
!       call gnrt_k0_n( il, ngen, inv, igen, jgen, a, b, c, ca, cb, cc, &
!            &          nbztyp_spg, &
!            &          k_sample_mesh(1,1),k_sample_mesh(2,1),k_sample_mesh(3,1), &
!            &          preallocation, nfout, ipri, knv3_t, rltv, &
!            &          kv3, vkxyz, qwgt, ipri_kp_t, trmat )
!                                                               !-(b_Kpoints)
!
       call gnrt_k0_n( il, ngen, inv, igen, jgen, a, b, c, ca, cb, cc, &
            &          nbztyp_spg, &
            &          k_sample_mesh(1,1),k_sample_mesh(2,1),k_sample_mesh(3,1), &
            &          preallocation, nfout, ipri, knv3_t, rltv, &
            &          kv3, vkxyz, qwgt, ipri_kp_t, trmat, &
            &          gen_tetramesh_mode, use_altv_rltv, altv, itrs, &
            &          gen_name_in_carts, &
            &          knv3_t_fbz, kv3_fbz, vkxyz_fbz, to_ibz_from_fbz_for_kpoint )
                                                         !-(b_Kpoints)
       if(ipri_kp_t >= 1) then
          write(nfout,'(" << trmat >>")')
          write(nfout,'(3f8.4)') ((trmat(i,j),j=1,3),i=1,3)
       end if
!        === modified by T. Y. 2015/01/16 ===
!!$       call gen_kpoint_in_full_BZ( preallocation, kv3_fbz, vkxyz_fbz )
       call gen_kpoint_in_full_BZ( preallocation )
!        ====================================
! ========================================================================= 12.0A

    else if(way_ksample == MESH) then

       if(ipri_kp>=1) call wd_k_sample_mesh          !-(contained here) k_sample_mesh->
       BZtype: select case(nbztyp)
       case (SURFACE, C2v_SURFACE, Csy_SURFACE)
          call call_k_generator(kpmsf0)
       case (WHOLE_BZ, ORTHORHOMBIC, RUTILE)
          call call_k_generator(kpmwbz0)
       case (SIMPLE_CUBIC)
          call call_k_generator(sccm0)
       case (BCC)
          call call_k_generator(bccm0)
       case (FCC, DIAMOND)
          call call_k_generator(fccm0)
       case (HEXAGONAL, 9)
          call call_k_generator(hexm0)
!!$          case (HEX1fold,HEX2fold,HEX3fold)
!!$          call kccrdf0(paramset,kv3,nfout,nbztyp,nfkpoint,nfmatbp &
!!$               &,vkxyz,qwgt,kv3)
       end select BZtype

    endif

! ==================================== modified by K. Tagami ======= 11.0
!!!    kv3 = kv3*nspin
!
    if ( noncol ) then             ! non-collinear case
       kv3 = kv3 *ndim_spinor;   kv3_fbz = kv3_fbz *ndim_spinor
    else
       kv3 = kv3 *nspin;         kv3_fbz = kv3_fbz *nspin
    endif
! =================================================================== 11.0

    if(.not.preallocation) then
       if(ipri_kp >= 2 ) write(nfout,*) ' << generate_or_read_k_sampling_points >>'

       if(way_ksample == MESH .and.  m_CtrlP_way_of_smearing() == TETRAHEDRON &
            & .and. nbztyp == WHOLE_BZ )   call change_kpoints_order

! ======================= modified by K. Tagami ================ 11.0
!!!       if(nspin == 2)  call double_k_points
!
       if ( noncol ) then
	  call gen_k_points_noncol
       else
          if ( nspin == 2 )  call double_k_points
       endif
! ========================================================= 11.0

       call k_points_in_BUCS_from_CARTS

       if ( way_ksample == FILE ) then
          if ( sw_force_kpt_inside_bz == ON ) then
             call m_Kp_force_kpoint_into_BZ
          endif
       endif

       call set_k_symmetry()
       if(ipri_kp>=1) call wd_kpoints()
! ==== KT_add === 2014/09/30
       if ( way_ksample == MONKHORST_PACK .or. way_ksample == MESH ) then
          if(ipri_kp>=2) call wd_kpoints_fbz()
       endif
! =============== 2014/09/30
! ==== ASMS ==== 2018/02/26
       if ( charge_symm_mode >= chg_symm_level1 ) then
          call m_Kp_chk_opr_from_fbz_to_ibz
       endif
       if ( charge_symm_mode >= chg_symm_level2 ) then
          call m_Kp_set_star_of_k
       endif
! ==== ASMS ==== 2018/02/26

       if ( way_ksample == FILE .and. sw_band_unfolding == ON ) then
          call wd_kpoints_refcell()
       endif

    else
       if(ipri_kp>=1) write(nfout,'(" !kp kv3 = ",i8, " nspin = ", i5)') kv3,nspin

! ======================= added by K. Tagami ================== 11.0
       if(ipri_kp>=1) write(nfout,*) 'ndim_spinor = ', ndim_spinor
! ============================================================= 11.0
    end if

    if(preallocation) then
       deallocate(vkxyz)
       deallocate(qwgt)
       deallocate(vkxyz_fbz)
       deallocate(to_ibz_from_fbz_for_kpoint)
       if ( sw_force_kpt_inside_bz == ON ) then
          if ( allocated(GvecTrans_kpt) ) deallocate(GvecTrans_kpt)
       endif
    end if

    call tstatc0_end(id_sname)
  contains
    subroutine set_k_symmetry()
      integer :: ik
      if(kimg == 2) then
         if(base_reduction_for_GAMMA == ON) then
            do ik = 1, kv3
               if(vkxyz(ik,1,CARTS)**2 + vkxyz(ik,2,CARTS)**2 &
                    & + vkxyz(ik,3,CARTS)**2 < DELTA) then
                  k_symmetry(ik) = GAMMA
               else
                  k_symmetry(ik) = 0
               end if
            end do
         else if(base_reduction_for_GAMMA == OFF .and. base_symmetrization_for_GAMMA == ON) then
            do ik = 1, kv3
               if(vkxyz(ik,1,CARTS)**2 + vkxyz(ik,2,CARTS)**2 &
                    & + vkxyz(ik,3,CARTS)**2 < DELTA) then
                  k_symmetry(ik) = GAMMA_base_symmetrization
               else
                  k_symmetry(ik) = 0
               end if
            end do
         else
            do ik = 1, kv3
               k_symmetry(ik) = 0
            end do
         end if
      else
         do ik = 1, kv3
            k_symmetry(ik) = 0
         end do
      end if

! ================================== added by K. Tagami =============== 11.0
      if ( noncol ) then
         write(nfout,*) '!----- '
         write(nfout,*) '! k_symmetry is forced to 0 in the non-collinear system'
         write(nfout,*) '!----- '
         k_symmetry(:) = 0
      endif
! ===================================================================== 11.0
    end subroutine set_k_symmetry

    subroutine wd_kpoints
      integer ik,i
      write(nfout,*)
      write(nfout,*) ' === k-points generated or read ==='
      write(nfout,'(" ik",8x,"CARTS",22x,"BUCS",19x,"QITG")')
      do ik = 1, kv3
         if(k_symmetry(ik) == GAMMA) then
            write(nfout,9002) ik,(vkxyz(ik,i,CARTS),i=1,3) &
                 &        ,           (vkxyz(ik,i,BUCS) ,i=1,3),qwgt(ik)
         else
            write(nfout,9001) ik,(vkxyz(ik,i,CARTS),i=1,3) &
                 &        ,           (vkxyz(ik,i,BUCS) ,i=1,3),qwgt(ik)
         end if
9001     format(i6,3f8.4,3x,3f8.4,3x,f8.4)
9002     format(i6,3f8.4,3x,3f8.4,3x,f8.4," GAMMA")
      enddo
    end subroutine wd_kpoints

! === KT_add === 2014/09/30
    subroutine wd_kpoints_fbz
      integer ik,i
      write(nfout,*)
      write(nfout,*) ' === k-points in full BZ ==='
      write(nfout,'(" ik",8x,"CARTS",22x,"BUCS",19x,"in IBZ")')
      do ik = 1, kv3_fbz
         write(nfout,9001) ik, (vkxyz_fbz(ik,i,CARTS),i=1,3), &
              &                (vkxyz_fbz(ik,i,BUCS) ,i=1,3), &
              &                to_ibz_from_fbz_for_kpoint(ik)
9001     format(i6,3f8.4,3x,3f8.4,3x,i6)
      enddo
    end subroutine wd_kpoints_fbz

    subroutine wd_kpoints_refcell
      integer ik,i
      write(nfout,*)
      write(nfout,*) ' === k-points in ref cell ==='
      write(nfout,'(" ik",8x,"CARTS",22x,"BUCS",19x,"in refcell")')
      do ik = 1, kv3
         write(nfout,9001) ik, (vkxyz_refcell(ik,i,CARTS),i=1,3), &
              &                (vkxyz_refcell(ik,i,BUCS) ,i=1,3)
9001     format(i6,3f8.4,3x,3f8.4)
      enddo
    end subroutine wd_kpoints_refcell

!!$    subroutine wd_kpoints(kv3,vkxyz)
!!$      integer, intent(in) :: kv3
!!$      real(kind=DP),intent(inout),dimension(kv3,3,CRDTYP) :: vkxyz
!!$
!!$      integer :: ik,i
!!$
!!$      write(nfout,*)
!!$      write(nfout,*) ' === k-points generated ==='
!!$      write(nfout,'(" ik",8x,"CARTS",22x,"BUCS",19x,"QITG")')
!!$      do ik = 1, kv3
!!$         write(nfout,9001) ik,(vkxyz(ik,i,CARTS),i=1,3) &
!!$              &     ,         (vkxyz(ik,i,BUCS) ,i=1,3),qwgt(ik)
!!$9001     format(i3,3f8.4,3x,3f8.4,3x,f8.4)
!!$      enddo
!!$    end subroutine wd_kpoints

    subroutine call_k_generator(external_k_gen)
      external   external_k_gen
      call external_k_gen(preallocation,nfout,altv,rltv,k_sample_mesh &
           &,knv3_t,kv3,vkxyz,qwgt)
    end subroutine call_k_generator

    subroutine double_k_points
      integer :: ik, is, k1
      integer, allocatable :: ilist(:)
      real(kind=DP), allocatable :: vkxyz_fbz_org(:,:,:)

      do ik = kv3/nspin, 1, -1
         vkxyz(ik*2,  1:3,CARTS) = vkxyz(ik,1:3,CARTS)
         vkxyz(ik*2-1,1:3,CARTS) = vkxyz(ik,1:3,CARTS)
         qwgt(ik*2)            = qwgt(ik)*0.5d0
         qwgt(ik*2-1)          = qwgt(ik)*0.5d0
      enddo
      if ( sw_band_unfolding == ON ) then
         do ik = kv3/nspin, 1, -1
            vkxyz_refcell(ik*2,  1:3,CARTS) = vkxyz_refcell(ik,1:3,CARTS)
            vkxyz_refcell(ik*2-1,1:3,CARTS) = vkxyz_refcell(ik,1:3,CARTS)
         enddo
      endif

      allocate( vkxyz_fbz_org(kv3_fbz,3,2 ) );
      vkxyz_fbz_org(1:kv3_fbz/nspin,:,:) = vkxyz_fbz(1:kv3_fbz/nspin,:,:)

      allocate( ilist(kv3_fbz) )
      ilist(1:kv3_fbz/nspin) = to_ibz_from_fbz_for_kpoint(1:kv3_fbz/nspin)

      do ik=1, kv3_fbz/nspin
         Do is=1, nspin
            k1 = ( ik -1 )*nspin + is
            if ( is == 1 ) then
               vkxyz_fbz( k1,:,CARTS) = vkxyz_fbz_org(ik,:,CARTS)
               to_ibz_from_fbz_for_kpoint(k1) = 2*ilist(ik) -1
            else
               vkxyz_fbz( k1,:,CARTS) = vkxyz_fbz_org(ik,:,CARTS)
               to_ibz_from_fbz_for_kpoint(k1) = 2*ilist(ik)
            endif
         End do
      End do

      deallocate( vkxyz_fbz_org ); deallocate( ilist )

    end subroutine double_k_points

! =========================== added by K. Tagami =========== 11.0
    subroutine gen_k_points_noncol
      integer :: ik, is, k1
      integer, allocatable :: ilist(:)
      real(kind=DP), allocatable :: vkxyz_org(:,:,:), qwgt_org(:)
      real(kind=DP), allocatable :: vkxyz_fbz_org(:,:,:)
!
!      allocate( vkxyz_org(kv3/ndim_spinor,3,2 ) );
      allocate( vkxyz_org(kv3,3,2 ) );
      vkxyz_org(1:kv3/ndim_spinor,:,:) = vkxyz(1:kv3/ndim_spinor,:,:)

!      allocate( qwgt_org(kv3/ndim_spinor) );
      allocate( qwgt_org(kv3) );
      qwgt_org(1:kv3/ndim_spinor) =  qwgt(1:kv3/ndim_spinor)
!
      do ik=1, kv3/ndim_spinor
         Do is=1, ndim_spinor
            k1 = ( ik -1 )*ndim_spinor + is

	    if ( is == 1 ) then
              vkxyz( k1,:,CARTS) = vkxyz_org(ik,:,CARTS)
              qwgt( k1 ) = qwgt_org(ik)
            else
              vkxyz( k1,:,CARTS) = vkxyz_org(ik,:,CARTS)
              qwgt( k1 ) = 0.0d0
           endif
         End do
      enddo
!
      if ( sw_band_unfolding == ON ) then
         do ik = kv3/nspin, 1, -1
            vkxyz_refcell(ik*2,  1:3,CARTS) = vkxyz_refcell(ik,1:3,CARTS)
            vkxyz_refcell(ik*2-1,1:3,CARTS) = vkxyz_refcell(ik,1:3,CARTS)
         enddo
      endif

      deallocate( vkxyz_org, qwgt_org )

!      allocate( vkxyz_fbz_org(kv3_fbz/ndim_spinor,3,2 ) );
      allocate( vkxyz_fbz_org(kv3_fbz,3,2 ) );
      vkxyz_fbz_org(1:kv3_fbz/ndim_spinor,:,:) = vkxyz_fbz(1:kv3_fbz/ndim_spinor,:,:)

!      allocate( ilist(kv3_fbz/ndim_spinor) )
      allocate( ilist(kv3_fbz) )
      ilist(1:kv3_fbz/ndim_spinor) = to_ibz_from_fbz_for_kpoint(1:kv3_fbz/ndim_spinor)

      do ik=1, kv3_fbz/ndim_spinor
         Do is=1, ndim_spinor
            k1 = ( ik -1 )*ndim_spinor + is
	    if ( is == 1 ) then
               vkxyz_fbz( k1,:,CARTS) = vkxyz_fbz_org(ik,:,CARTS)
               to_ibz_from_fbz_for_kpoint(k1) = 2*ilist(ik) -1
            else
               vkxyz_fbz( k1,:,CARTS) = vkxyz_fbz_org(ik,:,CARTS)
               to_ibz_from_fbz_for_kpoint(k1) = 2*ilist(ik) -1
            endif
         End do
      End do
      deallocate( vkxyz_fbz_org ); deallocate( ilist )

    end subroutine gen_k_points_noncol
! ========================================================== 11.0

    subroutine k_points_in_BUCS_from_CARTS
      integer       :: i, ik
      real(kind=DP) :: fact
      do i = 1, 3
         do ik = 1, kv3
            vkxyz(ik,i,BUCS) = (altv(1,i)*vkxyz(ik,1,CARTS)&
                 &           +  altv(2,i)*vkxyz(ik,2,CARTS)&
                 &           +  altv(3,i)*vkxyz(ik,3,CARTS))/PAI2
         enddo
      enddo

      ! -- >      check of BUCS system
      do i = 1, 3
         do ik = 1, kv3
            fact =   rltv(i,1)*vkxyz(ik,1,BUCS)&
                 & + rltv(i,2)*vkxyz(ik,2,BUCS)&
                 & + rltv(i,3)*vkxyz(ik,3,BUCS)
            if(dabs(fact - vkxyz(ik,i,CARTS)) > 1.d-5) then
               if(ipri_kp>=1) then
                  write(nfout,*) ' !D (i = ',i,',ik = ',ik,')'
                  write(nfout,*) vkxyz(ik,i,BUCS),fact, vkxyz(ik,i,CARTS)
               end if
               call phase_error_with_msg(nfout,' Coordinate transformation is invalid',__LINE__,__FILE__)
            endif
         enddo
      enddo

! === KT_add === 2014/09/30
      do i = 1, 3
         do ik = 1, kv3_fbz
            vkxyz_fbz(ik,i,BUCS) = ( altv(1,i) *vkxyz_fbz(ik,1,CARTS) &
                 &                    +  altv(2,i) *vkxyz_fbz(ik,2,CARTS) &
                 &                    +  altv(3,i) *vkxyz_fbz(ik,3,CARTS) )/PAI2
         enddo
      enddo
! ============= 2014/09/30

      if ( sw_band_unfolding == ON ) then
         do i = 1, 3
            do ik = 1, kv3
               vkxyz_refcell(ik,i,BUCS) &
                    &        = ( altv_refcell(1,i) *vkxyz_refcell(ik,1,CARTS) &
                    &          + altv_refcell(2,i) *vkxyz_refcell(ik,2,CARTS) &
                    &          + altv_refcell(3,i) *vkxyz_refcell(ik,3,CARTS) )/PAI2
            enddo
         enddo
      endif
    end subroutine k_points_in_BUCS_from_CARTS

    subroutine wd_k_sample_mesh
      integer :: i
      write(nfout,680) (k_sample_mesh(i,1),i=1,3)
      write(nfout,690) (k_sample_mesh(i,2),i=1,3)
680   format(' ',3i6,'   : nkx, nky, nkz ')
690   format(' ',3i6,'   : nkx2, nky2, nkz2 ')
    end subroutine wd_k_sample_mesh
  end subroutine m_Kp_gnrt_or_rd_k_points



!!$  subroutine m_Kp_rd_kpoints_table
!!$!    for tetrahedron method
!!$!     read k-table from the file 'f.kp0 '
!!$!           Tsuyoshi Miyazaki  '94.8.10
!!$#ifndef NO_TETRAHEDRON
!!$    integer :: nx1,ny1,nz1
!!$    integer,       pointer, dimension(:,:) :: ip2cub_wk
!!$    real(kind=DP), pointer, dimension(:,:) :: pa0,pb0
!!$
!!$    write(nfout,*) '-- read k-table from <f.kp0> --'
!!$    call m_Files_open_nfkindex
!!$    read(nfkindex,*) np0
!!$    read(nfkindex,*) np0,np1,np2
!!$    write(nfout,'(" np0,np1,np2 = ",3i6)') np0,np1,np2
!!$    if(np2 /= kv3/nspin) then
!!$       write(nfout,*) ' np0,np1,np2 ',np0,np1,np2
!!$       write(nfout,*) ' kv3/nspin ',kv3/nspin
!!$       write(nfout,*) ' np2  /=  (kv3/nspin)'
!!$       stop
!!$    endif
!!$    allocate(ip20(np0)) ; ip20 = 0
!!$    allocate(ip10(np0)) ; ip10 = 0
!!$    allocate(ip01(np1)) ; ip01 = 0
!!$    allocate(ip21(np1)) ; ip21 = 0
!!$    allocate(iu21(np1)) ; iu21 = 0
!!$    allocate(iv21(np1)) ; iv21 = 0
!!$    allocate(ip02(np2)) ; ip02 = 0
!!$    allocate(ip12(np2)) ; ip12 = 0
!!$    allocate(iwt(np2)) ; iwt = 0
!!$    allocate(ip2cub(np1)) ; ip2cub = 0
!!$    allocate(nxyz_tetra(3)) ; nxyz_tetra = 0
!!$    allocate(pa0(3,np0)) ; pa0 = 0.d0
!!$    allocate(pb0(3,np0)) ; pb0 = 0.d0
!!$    rewind nfkindex
!!$    write(nfout,'(" nfkindex = ",i9)') nfkindex
!!$    call nskpr0(nfkindex,nxyz_tetra(1),nxyz_tetra(2),nxyz_tetra(3) &
!!$           &    ,nx1,ny1,nz1,np0,np1,np2 &
!!$           &    ,pa0,pb0 &
!!$           &    ,ip10,ip20,ip01,ip21,ip02,ip12,iu21,iv21)
!!$    deallocate(pa0)
!!$    deallocate(pb0)
!!$    allocate(ip2cub_wk(9,nxyz_tetra(1)*nxyz_tetra(2)*nxyz_tetra(3)))
!!$    ip2cub_wk = 0.d0
!!$    call wtetra &
!!$      &  (nxyz_tetra(1),nxyz_tetra(2),nxyz_tetra(3),np0,np2,ip20 &
!!$      &  ,iwt,ip2cub &
!!$      &  ,ip2cub_wk)
!!$    deallocate(ip2cub_wk)
!!$#endif
!!$  end subroutine m_Kp_rd_kpoints_table

!! added by mizouchi@adv !!
!!!!!!!! modified by mizouchi@adv 2003.2.25 !!!!!
  subroutine m_Kp_cr_kpoints_table
!    for tetrahedron method
#ifndef NO_TETRAHEDRON
!!$    integer :: nx1,ny1,nz1,nd
    integer :: nx,ny,nz
    integer :: nxx,nyy,nzz
    integer :: ill,ii,lmnp0, lmnp1, lmnp2
    real(kind=DP), pointer, dimension(:,:) :: pa0,pb0,pb
    integer,       pointer, dimension(:,:) :: ka0,ka2
    integer,       pointer, dimension(:,:) :: ip2cub_wk
    integer,       pointer, dimension(:)   :: nstar2

    if(ipri_kp >= 2) write(nfout,*) '-- create k-table for tetrahedron method --'

    if(way_ksample == SKPS_DIRECT_IN) then
        call phase_error_with_msg(nfout,'way_ksample = SKPS_DIRECT_IN case is not'// &
             & 'supportted in tetrahedron method.',__LINE__,__FILE__)
    end if

    if(nbztyp.eq.WHOLE_BZ) then

        np0 = (k_sample_mesh(1,2)+1)*(k_sample_mesh(2,2)+1)*(k_sample_mesh(3,2)+1)
        np2 = k_sample_mesh(1,2)*k_sample_mesh(2,2)*k_sample_mesh(3,2)
        np1=np2

        if(ipri_kp>=1) write(nfout,'(" np0,np1,np2 = ",3i6)') np0,np1,np2

! ============================== modified by K. Tagami ============== 11.0
!        if(np2 /= kv3/nspin) then
!           if(ipri_kp>=1) then
!              write(nfout,*) ' np0,np1,np2 ',np0,np1,np2
!              write(nfout,*) ' kv3/nspin ',kv3/nspin
!              write(nfout,*) ' np2  /=  (kv3/nspin)'
!           end if
!           stop
!        endif
!
        if ( noncol ) then
          if(np2 /= kv3/ndim_spinor) then
             if(ipri_kp>=1) then
                write(nfout,*) ' np0,np1,np2 ',np0,np1,np2
                write(nfout,*) ' kv3/ndim_spinor ',kv3/ndim_spinor
                write(nfout,*) ' np2  /=  (kv3/ndim_spinor)'
             end if
             call phase_error_with_msg(nfout,'np2 /= kv3/ndim_spinor',__LINE__,__FILE__)
          endif
        else
          if(np2 /= kv3/nspin) then
             if(ipri_kp>=1) then
                write(nfout,*) ' np0,np1,np2 ',np0,np1,np2
                write(nfout,*) ' kv3/nspin ',kv3/nspin
                write(nfout,*) ' np2  /=  (kv3/nspin)'
             end if
             call phase_error_with_msg(nfout,'np2 /= kv3/nspin',__LINE__,__FILE__)
          endif
        endif
! ==================================================================== 11.0

        allocate(ip20(np0)) ; ip20 = 0
	if(ipri_kp>=2) write(nfout,'(" !kp ip20 is allocated")')
        allocate(iwt(np2)) ; iwt = 0
        allocate(ip2cub(np1)) ; ip2cub = 0
        allocate(nxyz_tetra(3)) ; nxyz_tetra = 0

        nx = k_sample_mesh(1,2)
        ny = k_sample_mesh(2,2)
        nz = k_sample_mesh(3,2)
        ill = 1

        lmnp0=np0
        lmnp1=np1

        allocate(ip01(np1)) ; ip01 = 0
        allocate(pa0(3,np0)) ; pa0 = 0.d0

        call nskma0(ill,nx,ny,nz,nxyz_tetra(1),nxyz_tetra(2) &
                    &,nxyz_tetra(3),nx1,ny1,nz1,nd)

        call nskp00(nxyz_tetra(1),nxyz_tetra(2),nxyz_tetra(3) &
                    &,nx1,ny1,nz1,nd,lmnp0,np0,pa0)

!! When nbztyp=1, ip10 in FLAPW program should be ip20 in this program !!

        call nskpbm(np0,lmnp0,lmnp1,pa0,np1,ip20,ip01)

        deallocate(pa0)
        deallocate(ip01)

    else

      if(nbztyp_spg == GENERAL .or.nbztyp_spg == GENERAL_LARGER) then

         nx=k_sample_mesh(1,1); ny=k_sample_mesh(2,1); nz=k_sample_mesh(3,1)

         call nskma0(il,nx,ny,nz,nxx,nyy,nzz,nx1,ny1,nz1,nd)
         lmnp0=(nxx+1)*(nyy+1)*(nzz+1)
         lmnp1=lmnp0
         lmnp2=lmnp0

         allocate(nxyz_tetra(3)) ; nxyz_tetra = 0
         allocate(ip10(lmnp0))  ; ip10 = 0
         allocate(ip20(lmnp0))  ; ip20 = 0
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
         allocate(ka2(4,lmnp2))  ; ka2 = 0

         ipri_kp_t = ipri_kp*(1-ipri_kp_count)
         if(ipri_kp_count == 0) ipri_kp_count = 1

! === KT_mod === 2015/02/18
!         call setkp0_n(il,ngen,inv,igen,jgen,a,b,c,ca,cb,cc &
!              & ,np2,np1,np0,lmnp0,lmnp1,lmnp2 &
!              & ,nx,ny,nz,nxyz_tetra(1),nxyz_tetra(2),nxyz_tetra(3) &
!              & ,ip10,ip20,ip01,ip02,ip21,ip12,iu21,iv21 &
!              & ,nstar2,pa0,pb0,pb,ka0,ka2 &
!              & ,ipri_kp_t)
         if ( gen_tetramesh_mode == 0 ) then
            call setkp0_n(il,ngen,inv,igen,jgen,a,b,c,ca,cb,cc &
                 & ,np2,np1,np0,lmnp0,lmnp1,lmnp2 &
                 & ,nx,ny,nz,nxyz_tetra(1),nxyz_tetra(2),nxyz_tetra(3) &
                 & ,ip10,ip20,ip01,ip02,ip21,ip12,iu21,iv21 &
                 & ,nstar2,pa0,pb0,pb,ka0,ka2 &
                 & ,ipri_kp_t, itrs )
         else if ( gen_tetramesh_mode == 1 ) then
            call setkp0_n_kt(il,ngen,inv,igen,jgen,a,b,c,ca,cb,cc &
                 & ,np2,np1,np0,lmnp0,lmnp1,lmnp2 &
                 & ,nx,ny,nz,nxyz_tetra(1),nxyz_tetra(2),nxyz_tetra(3) &
                 & ,ip10,ip20,ip01,ip02,ip21,ip12,iu21,iv21 &
                 & ,nstar2,pa0,pb0,pb,ka0,ka2 &
                 & ,ipri_kp_t, &
                 &  use_altv_rltv, altv, rltv, itrs, &
                 &  gen_name_in_carts )
         endif
! ================ 2015/02/18

      else
           nx = k_sample_mesh(1,1)
           ny = k_sample_mesh(2,1)
           nz = k_sample_mesh(3,1)

           if(ipri_kp>=2) write(nfout,'(" !kp nx,ny,nz = ",3i6)') nx,ny,nz
!!$           if(nbztyp_spg.eq.SIMPLE_CUBIC) il = 1
!!$           if(nbztyp_spg.eq.BCC) il = 3
!!$           if(nbztyp_spg.eq.FCC) il = 2
!!$           if(nbztyp_spg.eq.DIAMOND) il = 2
!!$           if(nbztyp_spg.eq.HEXAGONAL) il = 0

           call nskma0(il,nx,ny,nz,nxx,nyy,nzz,nx1,ny1,nz1,nd)
           if(ipri_kp >=2) write(nfout,'(" !kp nxx,nyy,nzz = ",3i6)') nxx,nyy,nzz

           lmnp0=(nxx+1)*(nyy+1)*(nzz+1)
           lmnp1=lmnp0
           lmnp2=lmnp0

           allocate(nxyz_tetra(3)) ; nxyz_tetra = 0
           allocate(ip10(lmnp0))  ; ip10 = 0
           allocate(ip20(lmnp0))  ; ip20 = 0
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
!!$           allocate(ka2(4,lmnp2))  ; ka2 = 0

!!$           call setkp0_default_n(il,ngen,inv,igen,jgen,a,b,c,ca,cb,cc &
!!$                &               ,nbztyp_spg,altv,nx,ny,nz &
!!$                &               ,np2,np1,np0,lmnp0,lmnp1,lmnp2 &
!!$                &               ,nxyz_tetra(1),nxyz_tetra(2),nxyz_tetra(3)&
!!$                &               ,ip10,ip20,ip01,ip02,ip21,ip12,iu21,iv21 &
!!$                &               ,nstar2,pa0,pb0,pb,ka0,ka2 &
!!$                &               ,ipri_kp)
           ipri_kp_t = ipri_kp*(1-ipri_kp_count)
           if(ipri_kp_count == 0) ipri_kp_count = 1

! ===================================== modified by K. Tagami ============== 12.0A
!           call setkp0_default_n(il,ngen,inv,igen,jgen,a,b,c,ca,cb,cc &
!                &               ,nx,ny,nz &
!                &               ,np2,np1,np0,lmnp0,lmnp1,lmnp2 &
!                &               ,nxyz_tetra(1),nxyz_tetra(2),nxyz_tetra(3)&
!                &               ,ip10,ip20,ip01,ip02,ip21,ip12,iu21,iv21 &
!                &               ,nstar2,pa0,pb0,pb,ka0 &
!                &               ,ipri_kp_t)
           if ( gen_tetramesh_mode == 0 ) then
              call setkp0_default_n(il,ngen,inv,igen,jgen,a,b,c,ca,cb,cc &
                   &               ,nx,ny,nz &
                   &               ,np2,np1,np0,lmnp0,lmnp1,lmnp2 &
                   &               ,nxyz_tetra(1),nxyz_tetra(2),nxyz_tetra(3)&
                   &               ,ip10,ip20,ip01,ip02,ip21,ip12,iu21,iv21 &
                   &               ,nstar2,pa0,pb0,pb,ka0 &
                   &               ,ipri_kp_t, itrs )
           else if ( gen_tetramesh_mode == 1 ) then
              call setkp0_default_n_kt(il,ngen,inv,igen,jgen,a,b,c,ca,cb,cc &
                   &               ,nx,ny,nz &
                   &               ,np2,np1,np0,lmnp0,lmnp1,lmnp2 &
                   &               ,nxyz_tetra(1),nxyz_tetra(2),nxyz_tetra(3)&
                   &               ,ip10,ip20,ip01,ip02,ip21,ip12,iu21,iv21 &
                   &               ,nstar2,pa0,pb0,pb,ka0 &
                   &               ,ipri_kp_t, &
                   &               use_altv_rltv, altv, rltv, itrs, &
                   &               gen_name_in_carts )
           else if ( gen_tetramesh_mode == 2 ) then
              call setkp0_default_n_kt2(il,ngen,inv,igen,jgen,a,b,c,ca,cb,cc &
                   &               ,nx,ny,nz &
                   &               ,np2,np1,np0,lmnp0,lmnp1,lmnp2 &
                   &               ,nxyz_tetra(1),nxyz_tetra(2),nxyz_tetra(3)&
                   &               ,ip10,ip20,ip01,ip02,ip21,ip12,iu21,iv21 &
                   &               ,nstar2,pa0,pb0,pb,ka0 &
                   &               ,ipri_kp_t, &
                   &               use_altv_rltv, altv, rltv, itrs, &
                   &               gen_name_in_carts )
           endif
! ========================================================================== 12.0A
        end if

        if(ipri_kp>=2) write(nfout,'(" np0,np1,np2 = ",3i6)') np0,np1,np2
        if(ipri_kp>=1) then
           write(nfout,'("!Kp nxyz_tetra(1:3) = ",3i8)') nxyz_tetra(1:3)
        end if

! ==================================== modified by K. Tagami =========== 11.0
!        if(np2 /= kv3/nspin) then
!           if(ipri_kp>=1) then
!              write(nfout,*) ' np0,np1,np2 ',np0,np1,np2
!              write(nfout,*) ' kv3/nspin ',kv3/nspin
!              write(nfout,*) ' np2  /=  (kv3/nspin)'
!           end if
!           stop
!        endif
!
        if ( noncol ) then
          if(np2 /= kv3/ndim_spinor) then
             if(ipri_kp>=1) then
                write(nfout,*) ' np0,np1,np2 ',np0,np1,np2
                write(nfout,*) ' kv3/ndim_spinor ',kv3/ndim_spinor
                write(nfout,*) ' np2  /=  (kv3/ndim_spinor)'
             end if
             call phase_error_with_msg(nfout,'np2 /= kv3/ndim_spinor',__LINE__,__FILE__)
          endif
        else
          if(np2 /= kv3/nspin) then
             if(ipri_kp>=1) then
                write(nfout,*) ' np0,np1,np2 ',np0,np1,np2
                write(nfout,*) ' kv3/nspin ',kv3/nspin
                write(nfout,*) ' np2  /=  (kv3/nspin)'
             end if
             call phase_error_with_msg(nfout,'np2 /= kv3/nspin',__LINE__,__FILE__)
          endif
        endif
! =================================================================== 11.0

        allocate(iwt(np2)) ; iwt = 0
        allocate(ip2cub(np1)) ; ip2cub = 0

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
         if(nbztyp_spg == GENERAL .or.nbztyp_spg == GENERAL_LARGER) &
              & deallocate(ka2)
    end if

    if(ipri_kp >= 2) then
       write(nfout,*) 'ip20'
       write(nfout,'(" np0 = ",i8)') np0
       write(nfout,*) (ip20(ii),ii=1,np0)
       write(nfout,'(" nx1,ny1,nz1,nd = ",4i6)') nx1,ny1,nz1,nd
    end if

    allocate(ip2cub_wk(9,nxyz_tetra(1)*nxyz_tetra(2)*nxyz_tetra(3)))
    ip2cub_wk = 0.d0
    call wtetra &
      &  (nxyz_tetra(1),nxyz_tetra(2),nxyz_tetra(3),np0,np2,ip20 &
      &  ,iwt,ip2cub &
      &  ,ip2cub_wk)
    deallocate(ip2cub_wk)
#endif
  end subroutine m_Kp_cr_kpoints_table

  subroutine change_kpoints_order
    integer :: nx,ny,nz,ik,ix,iy,iz &
          & , k0(k_sample_mesh(1,2),k_sample_mesh(2,2),k_sample_mesh(3,2)) &
          & , k1(kv3)
    real(kind=DP) :: wk(1:kv3,3)

    nx = k_sample_mesh(1,2)
    ny = k_sample_mesh(2,2)
    nz = k_sample_mesh(3,2)

! ======================================= modified by K. Tagami ======== 11.0
!    do ik = 1, kv3/nspin
!       wk(ik,1) = vkxyz(ik,1,CARTS)
!       wk(ik,2) = vkxyz(ik,2,CARTS)
!       wk(ik,3) = vkxyz(ik,3,CARTS)
!    end do

    if ( noncol ) then
      do ik = 1, kv3/ndim_spinor
         wk(ik,1) = vkxyz(ik,1,CARTS)
         wk(ik,2) = vkxyz(ik,2,CARTS)
         wk(ik,3) = vkxyz(ik,3,CARTS)
      end do
    else
      do ik = 1, kv3/nspin
         wk(ik,1) = vkxyz(ik,1,CARTS)
         wk(ik,2) = vkxyz(ik,2,CARTS)
         wk(ik,3) = vkxyz(ik,3,CARTS)
      end do
    endif
! ====================================================================== 11.0

    ik = 0
    do ix = 1,nx
    do iy = 1,ny
    do iz = 1,nz
        ik = ik+1
        k0(ix,iy,iz) = ik
        k1(ik) = 0
    end do
    end do
    end do

    ik = 0
    do iz = 1,nz
    do iy = 1,ny
    do ix = 1,nx
        ik = ik+1
        k1(ik) = k0(ix,iy,iz)
    end do
    end do
    end do

! ===================================== modified by K. Tagami ============= 11.0
!    do ik = 1, kv3/nspin
!      vkxyz(ik,1,CARTS) = wk(k1(ik),1)
!      vkxyz(ik,2,CARTS) = wk(k1(ik),2)
!      vkxyz(ik,3,CARTS) = wk(k1(ik),3)
!    end do

    if ( noncol ) then
      do ik = 1, kv3/ndim_spinor
        vkxyz(ik,1,CARTS) = wk(k1(ik),1)
        vkxyz(ik,2,CARTS) = wk(k1(ik),2)
        vkxyz(ik,3,CARTS) = wk(k1(ik),3)
      end do
    else
      do ik = 1, kv3/nspin
        vkxyz(ik,1,CARTS) = wk(k1(ik),1)
        vkxyz(ik,2,CARTS) = wk(k1(ik),2)
        vkxyz(ik,3,CARTS) = wk(k1(ik),3)
      end do
    endif
! ======================================================================== 11.0

  end subroutine change_kpoints_order


  !  <<< Monkhorst-Pack special kpoints generator >>>
  !  coded by Takenori YAMAMOTO (Univ. Tokyo)
  !  last update : Aug. 12, 2003
! ========================================== modified by K. Tagami ======== 11.0P
!  subroutine special_kpoints(printlevel,output,nsym,rot,b1,b2,b3,a1,a2,a3 &
!    & ,nkpoint,nkpmax,kpoint,weight,magnetic,use_trs )

  subroutine special_kpoints( printlevel, output, nsym, rot, b1, b2, b3, a1, a2, a3, &
       &                      nkpoint, nkpmax, kpoint, weight, magnetic, use_trs, &
       &                      nkpoint_full, kpoint_full, to_ibz_from_fbz_for_kpoint, &
       &                      nsym_before_sym_reduction, rot_before_sym_reduction )
! ========================================================================== 11.0P
    implicit none

    integer, intent(in) :: printlevel
    integer, intent(in) :: output
    integer, intent(in) :: nsym
    real(kind=DP), intent(in) :: rot(3,3,nsym)
    real(kind=DP), intent(in) :: b1(3),b2(3),b3(3)
    real(kind=DP), intent(in) :: a1(3),a2(3),a3(3)
    integer, intent(inout) :: nkpoint
    integer, intent(in) :: nkpmax
    real(kind=DP), intent(inout) :: kpoint(3,*),weight(*)
    logical, intent(in) :: magnetic
    logical, intent(in) :: use_trs

! === KT_add === 2014/09/30
    integer, optional, intent(inout) :: nkpoint_full
    real(kind=DP), optional, intent(inout) :: kpoint_full(3,*)
    integer, optional, intent(inout) :: to_ibz_from_fbz_for_kpoint(*)
! ============== 2014/09/30

! ============================ aded by K. Tagami ======================== 11.0P
    integer, optional, intent(in) :: nsym_before_sym_reduction
    real(kind=DP), optional, intent(in) :: rot_before_sym_reduction(3,3,*)
! ======================================================================= 11.0P

    ! local variables
    integer :: i,nkmesh,nmp_kmesh,nfbz_kmesh
    real(kind=DP), allocatable :: kmesh(:,:)
    real(kind=DP), allocatable :: mp_kmesh(:,:),fbz_kmesh(:,:)
    integer :: nsym_trs
    real(kind=DP), allocatable :: rot_trs(:,:,:),grot(:,:,:)
    real(kind=DP) :: kshf(3)
    integer :: face(3,nface_max),nface

! ============================= added by K. Tagami =============== 11.0P
    integer :: nsym_trs_before_sym_reduction
    real(kind=DP), allocatable :: rot_trs_before_sym_reduction(:,:,:)
    real(kind=DP), allocatable :: grot_before_sym_reduction(:,:,:)

    logical :: use_rot_before_sym_reduction
! ================================================================= 11.0P

    real(kind=DP) :: kshift_norm2
    if(printable .and. printlevel .ge. 0) then
      write(output,*) 'Special k-points generation starts.'
    end if

    ! debug
    !if(printable .and. printlevel .ge. 3) then
    !  write(output,'("Symmetry operation matirx")')
    !  do i=1,nsym
    !    write(output,*) 'n=',i
    !    write(output,'(3(1x,f10.5))') rot(1,1:3,i)
    !    write(output,'(3(1x,f10.5))') rot(2,1:3,i)
    !    write(output,'(3(1x,f10.5))') rot(3,1:3,i)
    !  end do
    !end if
    ! end debug

    call set_number_nsym_trs( nsym, rot, magnetic, use_trs, nsym_trs )
    allocate( rot_trs(3,3,nsym_trs), grot(3,3,nsym_trs) )
    call gen_rot_trs( nsym, rot, nsym_trs, rot_trs )
    call get_grot(a1,a2,a3,b1,b2,b3,nsym_trs,rot_trs,grot)

! ============================= added by K. Tagami =============== 11.0P
    if ( noncol .and. present( nsym_before_sym_reduction ) ) then
       use_rot_before_sym_reduction = .true.
    else
       use_rot_before_sym_reduction = .false.
    endif

    if ( use_rot_before_sym_reduction ) then
       call set_number_nsym_trs( nsym_before_sym_reduction, rot_before_sym_reduction, &
            &                    magnetic, use_trs, nsym_trs_before_sym_reduction )
       allocate( rot_trs_before_sym_reduction(3,3,nsym_trs_before_sym_reduction) )
       allocate( grot_before_sym_reduction(3,3,nsym_trs_before_sym_reduction) )
       call gen_rot_trs( nsym_before_sym_reduction, rot_before_sym_reduction, &
            &            nsym_trs_before_sym_reduction, rot_trs_before_sym_reduction )
       call get_grot( a1, a2, a3, b1, b2, b3, nsym_trs_before_sym_reduction, &
            &         rot_trs_before_sym_reduction, grot_before_sym_reduction )
    endif
! =================================================================== 11.0P

    ! debug
    !if(printable .and. printlevel .ge. 3) then
    !  write(output,'("Symmetry operation matirx for reciprocal lattice indeces")')
    !  do n=1,nsym
    !    write(output,*) 'n=',n
    !    write(output,'(3(1x,f10.5))') grot(1,1:3,n)
    !    write(output,'(3(1x,f10.5))') grot(2,1:3,n)
    !    write(output,'(3(1x,f10.5))') grot(3,1:3,n)
    !  end do
    !end if
    ! end debug

    if (nkpoint == 0) then

      nkmesh = (2*mp_index(1)+1)*(2*mp_index(2)+1)*(2*mp_index(3)+1)
      allocate(kmesh(3,nkmesh*nsym_trs))
      call gen_gmesh(mp_index,nkmesh,kmesh)

      kshf(1:3) = kshift(1:3)/mp_index(1:3)
      do i=1,nkmesh
        kmesh(1:3,i) = kmesh(1:3,i) + kshf(1:3)
      end do

      ! debug
      !do i=1,nkmesh
      !  write(output,'(3(1x,f10.5))') kmesh(1:3,i)
      !end do
      !stop 'debug kmesh'
      ! end debug

      call fbz_face(b1,b2,b3,nface,face)

      allocate(mp_kmesh(3,nkmesh))
      allocate(fbz_kmesh(3,nkmesh))

      if ( sw_kpt_segmentation == ON ) then
         call first_bz_segment(b1,b2,b3,nkmesh,kmesh,nmp_kmesh,mp_kmesh,nface,face)
      else
         call first_bz(b1,b2,b3,nkmesh,kmesh,nmp_kmesh,mp_kmesh,nface,face)
      endif

      if(printable .and. printlevel .ge. 0) &
      &  write(output,'(1x,"number of k-points in MP mesh        = ",i8)') nmp_kmesh

      ! debug
      !do i=1,nmp_kmesh
      !  write(output,'(3(1x,f10.5))') mp_kmesh(1:3,i)
      !end do
      !stop 'debug mp_kmesh'
      ! end debug

! == KT_add == 2014/09/30
#ifndef EXPERIMENTAL_REGEN_KMESH
      if ( present( kpoint_full ) ) then
         if ( sw_kpt_segmentation == ON ) then
            call regen_kmesh_using_grot_segment( nkmesh, kmesh, mp_kmesh, &
                 &                               1, grot(1:3,1:3,1) )
         else
            call regen_kmesh_using_grot( nkmesh, kmesh, mp_kmesh, 1, grot(1:3,1:3,1) )
         endif
         if ( sw_kpt_segmentation == ON ) then
            call first_bz_segment(b1,b2,b3,nkmesh,kmesh,nfbz_kmesh,fbz_kmesh,nface,face)
         else
            call first_bz(b1,b2,b3,nkmesh,kmesh,nfbz_kmesh,fbz_kmesh,nface,face)
         endif
         nkpoint_full = nfbz_kmesh
         kpoint_full(1:3,1:nkpoint_full) = fbz_kmesh(1:3,1:nkpoint_full)
      endif
#endif
! ============ 2014/09/30

      if ( use_rot_before_sym_reduction ) then
         call regen_kmesh_using_grot( nkmesh, kmesh, mp_kmesh, &
              &                       nsym_trs_before_sym_reduction, &
              &                       grot_before_sym_reduction )
      else
         if ( sw_kpt_segmentation == ON ) then
            call regen_kmesh_using_grot_segment( nkmesh, kmesh, mp_kmesh, &
                 &                               nsym_trs, grot )
         else
            call regen_kmesh_using_grot( nkmesh, kmesh, mp_kmesh, nsym_trs, grot )
         endif
      endif

      if ( sw_kpt_segmentation == ON ) then
         call first_bz_segment(b1,b2,b3,nkmesh,kmesh,nfbz_kmesh,fbz_kmesh,nface,face)
      else
         call first_bz(b1,b2,b3,nkmesh,kmesh,nfbz_kmesh,fbz_kmesh,nface,face)
      endif

      if(printable .and. printlevel .ge. 0) &
      &  write(output,'(1x,"number of k-points in full BZ        = ",i8)') nfbz_kmesh

      ! debug
      !do i=1,nfbz_kmesh
      !  write(output,'(3(1x,f10.5))') fbz_kmesh(1:3,i)
      !end do
      !stop 'debug fbz_kmesh'
      ! end debug

#ifndef EXPERIMENTAL_REGEN_KMESH
      if ( present( kpoint_full ) ) then
         if ( sw_kpt_segmentation == ON ) then
            call gen_spk2_segment( nsym_trs, grot, nkpoint_full, kpoint_full, nkpoint, &
                 &                 nkpmax, kpoint, weight, to_ibz_from_fbz_for_kpoint )
         else
            call gen_spk2( nsym_trs, grot, nkpoint_full, kpoint_full, nkpoint, nkpmax, &
                 &        kpoint, weight, to_ibz_from_fbz_for_kpoint )
         endif
      else
         call gen_spk( nsym_trs, grot, nfbz_kmesh, fbz_kmesh, nkpoint, nkpmax, &
              &        kpoint, weight )
      endif
#else
! == KT_add == 2014/09/30
      if ( present( kpoint_full ) ) then
         nkpoint_full = nfbz_kmesh
         kpoint_full(1:3,1:nkpoint_full) = fbz_kmesh(1:3,1:nkpoint_full)
      endif
      call gen_spk2( nsym_trs, grot, nfbz_kmesh, fbz_kmesh, nkpoint, nkpmax, &
           &         kpoint, weight, to_ibz_from_fbz_for_kpoint )
! ============ 2014/09/30
#endif

! === ASMS ====
!      call positive(nsym_trs,grot,nkpoint,kpoint)

      kshift_norm2 = kshift(1)**2 +kshift(2)**2 +kshift(3)**2
      if ( sw_hybrid_functional == ON .and. kshift_norm2 > 0 ) then
         sw_fix_ibz_on_fbz_mesh = ON
         write(nfout,*) "!** sw_fix_ibz_on_fbz_mesh is forced to ON"
      endif
      if ( sw_fix_ibz_on_fbz_mesh == OFF ) then
         if ( charge_symm_mode >= chg_symm_level1 ) then
            sw_fix_ibz_on_fbz_mesh = ON
            write(nfout,*) "!** sw_fix_ibz_on_fbz_mesh is forced to ON"
         endif
      endif

      if ( present( kpoint_full ) ) then
         if ( sw_fix_ibz_on_fbz_mesh == ON ) then
            call positive2( nsym_trs, grot, nkpoint, kpoint, nkpoint_full, kpoint_full )
         else
            call positive(nsym_trs,grot,nkpoint,kpoint)
         endif
      else
         call positive(nsym_trs,grot,nkpoint,kpoint)
      endif
! === ASMS ====

      deallocate(rot_trs,mp_kmesh,kmesh,fbz_kmesh)

    end if

    !!call accuracy(a1,a2,a3,nsym_trs,grot,nkpoint,kpoint,weight)

    deallocate(grot)

! ============================= added by K. Tagami ======================== 11.0P
    if ( allocated(grot_before_sym_reduction) ) deallocate( grot_before_sym_reduction )
    if ( allocated(rot_trs_before_sym_reduction) ) deallocate( rot_trs_before_sym_reduction )
! ========================================================================= 11.0P

    if(printable .and. printlevel .ge. 0) then
      write(output,*) 'Special k-points generation ends.'
    end if

  contains

! ================================ added by K. Tagami ================= 11.0P
    subroutine regen_kmesh_using_grot( nkmesh, kmesh, mp_kmesh, nsym_trs, grot )
      integer, intent(in) :: nsym_trs
      integer, intent(inout) :: nkmesh
      real(kind=DP), intent(inout) :: kmesh(3,nkmesh*nsym_trs)
      real(kind=DP), intent(in) :: mp_kmesh(3,nkmesh)
      real(kind=DP), intent(in) :: grot(3,3,nsym_trs)

      real(kind=DP), allocatable :: kt_mesh(:,:)

      integer :: i, j, n
      logical :: new_kp
      real(kind=DP) :: v(3)

      allocate( kt_mesh(3,nmp_kmesh) ); kt_mesh = 0

#ifndef EXPERIMENTAL_REGEN_KMESH
      Do i=1, nmp_kmesh
         kt_mesh(1,i) = mp_kmesh(1,i)
         kt_mesh(2,i) = mp_kmesh(2,i)
         kt_mesh(3,i) = mp_kmesh(3,i)
      End Do
#else
      Do i=1, nmp_kmesh
         kt_mesh(1,i) = mp_kmesh(1,i) -kshf(1)
         kt_mesh(2,i) = mp_kmesh(2,i) -kshf(2)
         kt_mesh(3,i) = mp_kmesh(3,i) -kshf(3)
      End Do
#endif

      nkmesh=1
      kmesh(1:3,nkmesh) = kt_mesh(1:3,1)

      do i=1,nmp_kmesh
         do n=1,nsym_trs
            v(1:3) = grot(1:3,1,n)*kt_mesh(1,i) +grot(1:3,2,n)*kt_mesh(2,i) &
                 & + grot(1:3,3,n)*kt_mesh(3,i)

            new_kp = .true.
            do j=1,nkmesh
               if (abs(v(1)-kmesh(1,j)) < 1.d-10 .and. &
                    &abs(v(2)-kmesh(2,j)) < 1.d-10 .and. &
                    &abs(v(3)-kmesh(3,j)) < 1.d-10 ) then
                  new_kp = .false.
                  exit
               end if
            end do
            if(new_kp) then
               nkmesh = nkmesh + 1
               kmesh(1:3,nkmesh) = v(1:3)
            end if
         end do
      end do
#ifdef EXPERIMENTAL_REGEN_KMESH
      Do i=1, nkmesh
         kmesh(1,i) = kmesh(1,i) +kshf(1)
         kmesh(2,i) = kmesh(2,i) +kshf(2)
         kmesh(3,i) = kmesh(3,i) +kshf(3)
      End Do
#endif
      deallocate( kt_mesh )

    end subroutine regen_kmesh_using_grot

    subroutine regen_kmesh_using_grot_segment( nkmesh, kmesh, mp_kmesh, nsym_trs, grot )
      integer, intent(in) :: nsym_trs
      integer, intent(inout) :: nkmesh
      real(kind=DP), intent(inout) :: kmesh(3,nkmesh*nsym_trs)
      real(kind=DP), intent(in) :: mp_kmesh(3,nkmesh)
      real(kind=DP), intent(in) :: grot(3,3,nsym_trs)

      real(kind=DP), allocatable :: kt_mesh(:,:)

      integer :: i, j, n
      logical :: new_kp
      real(kind=DP) :: v(3)
! ===
      integer :: ix, iy, iz, ndiv_x, ndiv_y, ndiv_z
      integer :: nx, ny, nz
      integer :: range_xs, range_ys, range_zs, range_xe, range_ye, range_ze
      integer :: m1, m2, num1, num2, num3, j1, j2, n1, n2
      real(kind=DP) :: w(3), resol_x, resol_y, resol_z, dx, dy, dz
      real(kind=DP) :: eps = 1.0D-10
      real(kind=DP) :: eps2 = 1.0D-3

      type t_kpt_segment
         integer :: nkpt
         integer, allocatable :: list(:,:)
         integer :: nkpt_reg
         integer, allocatable :: list_reg(:,:)
      end type t_kpt_segment
      type(t_kpt_segment), allocatable :: kpt_segment(:,:,:)
      if ( kpt_save_memory_mode == 1 ) then
         if ( npes > 1 .and. mype /= 0 ) goto 1000
      endif

      allocate( kt_mesh(3,nmp_kmesh) ); kt_mesh = 0

#ifndef EXPERIMENTAL_REGEN_KMESH
      Do i=1, nmp_kmesh
         kt_mesh(1,i) = mp_kmesh(1,i)
         kt_mesh(2,i) = mp_kmesh(2,i)
         kt_mesh(3,i) = mp_kmesh(3,i)
      End Do
#else
      Do i=1, nmp_kmesh
         kt_mesh(1,i) = mp_kmesh(1,i) -kshf(1)
         kt_mesh(2,i) = mp_kmesh(2,i) -kshf(2)
         kt_mesh(3,i) = mp_kmesh(3,i) -kshf(3)
      End Do
#endif

      resol_x = 1.0d0 /dble(ndiv_segment(1))
      resol_y = 1.0d0 /dble(ndiv_segment(2))
      resol_z = 1.0d0 /dble(ndiv_segment(3))
!
#if 1
      ndiv_x = nint( ( 1.0d0 +kshift(1) ) /resol_x )
      ndiv_y = nint( ( 1.0d0 +kshift(2) ) /resol_y )
      ndiv_z = nint( ( 1.0d0 +kshift(3) ) /resol_z )
#else
      ndiv_x = int( ( 1.0d0 +kshift(1) ) /resol_x ) +1
      ndiv_y = int( ( 1.0d0 +kshift(2) ) /resol_y ) +1
      ndiv_z = int( ( 1.0d0 +kshift(3) ) /resol_z ) +1
#endif

      allocate( kpt_segment( -ndiv_x:ndiv_x, -ndiv_y:ndiv_y, -ndiv_z:ndiv_z ) )
      kpt_segment(:,:,:)%nkpt = 0

      do i=1,nmp_kmesh
         do n=1,nsym_trs
#if 1
            if ( n==1 ) then
               v(1:3) = kt_mesh(1:3,i)
            else
               v(1:3) = grot(1:3,1,n)*kt_mesh(1,i) +grot(1:3,2,n)*kt_mesh(2,i) &
                    & + grot(1:3,3,n)*kt_mesh(3,i)
            endif
#else
            v(1:3) = grot(1:3,1,n)*kt_mesh(1,i) +grot(1:3,2,n)*kt_mesh(2,i) &
                 & + grot(1:3,3,n)*kt_mesh(3,i)
#endif
            ix = nint( v(1) /resol_x )
            iy = nint( v(2) /resol_y )
            iz = nint( v(3) /resol_z )
            kpt_segment(ix,iy,iz)%nkpt = kpt_segment(ix,iy,iz)%nkpt +1
         end do
      end do
      Do ix=-ndiv_x, ndiv_x
         Do iy=-ndiv_y, ndiv_y
            Do iz=-ndiv_z, ndiv_z
               num1 = kpt_segment(ix,iy,iz)%nkpt
               if ( num1 > 0 ) allocate( kpt_segment(ix,iy,iz)%list(num1,2) )
            End Do
         ENd Do
      end Do
      kpt_segment(:,:,:)%nkpt = 0

      do i=1,nmp_kmesh
         do n=1,nsym_trs
#if 1
            if ( n==1 ) then
               v(1:3) = kt_mesh(1:3,i)
            else
               v(1:3) = grot(1:3,1,n)*kt_mesh(1,i) +grot(1:3,2,n)*kt_mesh(2,i) &
                    & + grot(1:3,3,n)*kt_mesh(3,i)
            endif
#else
            v(1:3) = grot(1:3,1,n)*kt_mesh(1,i) +grot(1:3,2,n)*kt_mesh(2,i) &
                 & + grot(1:3,3,n)*kt_mesh(3,i)
#endif
            ix = nint( v(1) /resol_x )
            iy = nint( v(2) /resol_y )
            iz = nint( v(3) /resol_z )

            kpt_segment(ix,iy,iz)%nkpt = kpt_segment(ix,iy,iz)%nkpt +1
            num1 = kpt_segment(ix,iy,iz)%nkpt
            kpt_segment(ix,iy,iz)%list(num1,1) = i
            kpt_segment(ix,iy,iz)%list(num1,2) = n
         End do
      End do

      nkmesh = 0

      kpt_segment(:,:,:)%nkpt_reg = 0
      Do ix=-ndiv_x, ndiv_x
         Do iy=-ndiv_y, ndiv_y
            Do iz=-ndiv_z, ndiv_z
               num1 = kpt_segment(ix,iy,iz)%nkpt
               if ( num1 > 0 ) allocate( kpt_segment(ix,iy,iz)%list_reg(num1,2 ) )
            End Do
         ENd Do
      end Do

#if 1
      nkmesh=1;     kmesh(1:3,nkmesh) = kt_mesh(1:3,1)
      v(1:3) = kt_mesh(1:3,1)
      ix = nint( v(1) /resol_x )
      iy = nint( v(2) /resol_y )
      iz = nint( v(3) /resol_z )
      kpt_segment(ix,iy,iz)%nkpt_reg = kpt_segment(ix,iy,iz)%nkpt_reg +1
      num1 = kpt_segment(ix,iy,iz)%nkpt_reg
      kpt_segment(ix,iy,iz)%list_reg(num1,1) = 1
      kpt_segment(ix,iy,iz)%list_reg(num1,2) = 1       ! E
#endif

      do i=1,nmp_kmesh
         do n=1,nsym_trs
#if 1
            if ( n==1 ) then
               v(1:3) = kt_mesh(1:3,i)
            else
               v(1:3) = grot(1:3,1,n)*kt_mesh(1,i) +grot(1:3,2,n)*kt_mesh(2,i) &
                    & + grot(1:3,3,n)*kt_mesh(3,i)
            endif
#else
            v(1:3) = grot(1:3,1,n)*kt_mesh(1,i) +grot(1:3,2,n)*kt_mesh(2,i) &
                 & + grot(1:3,3,n)*kt_mesh(3,i)
#endif
            ix = nint( v(1) /resol_x )
            iy = nint( v(2) /resol_y )
            iz = nint( v(3) /resol_z )
!
            dx = v(1) -ix *resol_x
            dy = v(2) -iy *resol_y
            dz = v(3) -iz *resol_z
!
            range_xs = 0;  range_ys = 0;    range_zs = 0
            range_xe = 0;  range_ye = 0;    range_ze = 0
            if ( dx < -resol_x /2. +eps2 ) range_xs = -1
            if ( dx >  resol_x /2. -eps2 ) range_xe =  1
            if ( dy < -resol_y /2. +eps2 ) range_ys = -1
            if ( dy >  resol_y /2. -eps2 ) range_ye =  1
            if ( dz < -resol_z /2. +eps2 ) range_zs = -1
            if ( dz >  resol_z /2. -eps2 ) range_ze =  1

            new_kp = .true.

            nxloop: Do nx=range_xs, range_xe
               Do ny=range_ys, range_ye
                  Do nz=range_zs, range_ze
                     num1 = ix +nx
                     num2 = iy +ny
                     num3 = iz +nz
                     if ( num1 < -ndiv_x .or. num1 > ndiv_x ) cycle
                     if ( num2 < -ndiv_y .or. num2 > ndiv_y ) cycle
                     if ( num3 < -ndiv_z .or. num3 > ndiv_z ) cycle

                     do j1=1, kpt_segment(num1,num2,num3)%nkpt_reg
                        m1 = kpt_segment(num1,num2,num3)%list_reg(j1,1)
                        m2 = kpt_segment(num1,num2,num3)%list_reg(j1,2)
                        w(1:3) = grot(1:3,1,m2) *kt_mesh(1,m1) &
                             &  +grot(1:3,2,m2) *kt_mesh(2,m1) &
                             &  +grot(1:3,3,m2) *kt_mesh(3,m1)
                        if ( abs(v(1)-w(1)) < eps &
                             & .and. abs(v(2)-w(2)) < eps &
                             & .and. abs(v(3)-w(3)) < eps ) then
                           new_kp = .false.;  exit nxloop
                        end if
                     end do

                  end Do
               end Do
            end Do nxloop

            if (new_kp) then
               nkmesh = nkmesh +1;
               kmesh(1:3,nkmesh) = v(1:3)
               kpt_segment(ix,iy,iz)%nkpt_reg = kpt_segment(ix,iy,iz)%nkpt_reg +1

               num1 = kpt_segment(ix,iy,iz)%nkpt_reg
               kpt_segment(ix,iy,iz)%list_reg(num1,1) = i
               kpt_segment(ix,iy,iz)%list_reg(num1,2) = n
            end if
         End do
      End do
      deallocate( kpt_segment )

#ifdef EXPERIMENTAL_REGEN_KMESH
      Do i=1, nkmesh
         kmesh(1,i) = kmesh(1,i) +kshf(1)
         kmesh(2,i) = kmesh(2,i) +kshf(2)
         kmesh(3,i) = kmesh(3,i) +kshf(3)
      End Do
#endif
      deallocate( kt_mesh )

1000  continue
      if ( kpt_save_memory_mode == 1 ) then
         if ( npes > 1 ) then
            call mpi_bcast( nkmesh, 1, mpi_integer, 0, MPI_CommGroup, ierr )
            call mpi_bcast( kmesh, nkmesh, mpi_double_precision, &
                 &          0, MPI_CommGroup, ierr )
         endif
      endif

    end subroutine regen_kmesh_using_grot_segment

    subroutine set_number_nsym_trs( nsym, rot, magnetic, use_trs, nsym_trs )
      integer, intent(in) :: nsym
      logical, intent(in) :: magnetic, use_trs
      real(kind=DP), intent(in) :: rot(3,3,nsym)

      integer, intent(out) :: nsym_trs

      integer :: i
      logical :: with_inverse

      real(kind=DP), parameter :: eps = 1.d-8

! check if the inverse operation is in rot
      with_inverse = .false.

      do i=1,nsym
         if ( abs(rot(1,1,i)+1.d0) < eps .and. &
              & abs(rot(2,2,i)+1.d0) < eps .and. &
              & abs(rot(3,3,i)+1.d0) < eps .and. &
              & abs(rot(1,2,i)) < eps .and. &
              & abs(rot(2,1,i)) < eps .and. &
              & abs(rot(2,3,i)) < eps .and. &
              & abs(rot(3,2,i)) < eps .and. &
              & abs(rot(3,1,i)) < eps .and. &
              & abs(rot(1,3,i)) < eps ) then
            with_inverse = .true.
         end if
      end do

! ----------------------------------------------
!    if(magnetic) with_inverse = .true.
! ---------------------------------------------

#if 0
      if ( magnetic ) then
         if ( .not. noncol ) with_inverse = .true.
      endif
#endif

! ----------------------------------------------
!    if(with_inverse) then
!      nsym_trs = nsym
!    else
!      if(use_trs) then
!         nsym_trs = nsym*2
!      else
!         nsym_trs = nsym
!      end if
!    end if
!    allocate(rot_trs(3,3,nsym_trs),grot(3,3,nsym_trs))
!    do i=1,nsym
!      rot_trs(1:3,1:3,i) = rot(1:3,1:3,i)
!    end do
!    if(use_trs .and. .not.with_inverse) then
!      do i=nsym+1,nsym_trs
!        rot_trs(1:3,1:3,i) = -rot(1:3,1:3,i-nsym)
!      end do
!    end if
! ------------------------------------------------------

      if ( noncol ) then
         if(with_inverse) then
            nsym_trs = nsym
         else
            if(use_trs) then
               nsym_trs = nsym*2
            else
               nsym_trs = nsym
            end if
         end if
      else
         if(with_inverse) then
            nsym_trs = nsym
         else
            if(use_trs) then
               nsym_trs = nsym*2
            else
               nsym_trs = nsym
            end if
         end if
      endif

    end subroutine set_number_nsym_trs

    subroutine gen_rot_trs( nsym, rot, nsym_trs, rot_trs )
      integer, intent(in) :: nsym, nsym_trs
      real(kind=DP), intent(in) :: rot(3,3,nsym)
      real(kind=DP), intent(out) :: rot_trs(3,3,nsym_trs)

      integer :: i

      do i=1,nsym
         rot_trs(1:3,1:3,i) = rot(1:3,1:3,i)
      end do
      if ( nsym_trs > nsym ) then
         do i=nsym+1,nsym_trs
            rot_trs(1:3,1:3,i) = -rot(1:3,1:3,i-nsym)
         end do
      end if

    end subroutine gen_rot_trs
! ====================================================================== 11.0P

    subroutine gen_gmesh(mp_index,nmesh,gmesh)
      implicit none

      integer, intent(in) :: mp_index(3)
      integer, intent(in) :: nmesh
      real(kind=DP), intent(out) :: gmesh(3,nmesh)

      ! local variables
      integer :: i,j,k,n

      n=0
      do i=-mp_index(1),mp_index(1)
        do j=-mp_index(2),mp_index(2)
          do k=-mp_index(3),mp_index(3)
            n=n+1
            gmesh(1,n)=dble(i)/mp_index(1)
            gmesh(2,n)=dble(j)/mp_index(2)
            gmesh(3,n)=dble(k)/mp_index(3)
          end do
        end do
      end do

    end subroutine gen_gmesh

    subroutine first_bz(b1,b2,b3,nkmesh,kmesh,nfbz_kmesh,fbz_kmesh &
    & ,nface,face)
      implicit none

      real(kind=DP), intent(in) :: b1(3),b2(3),b3(3)
      integer, intent(in) :: nkmesh
      real(kind=DP), intent(in) :: kmesh(3,nkmesh)
      integer, intent(out) :: nfbz_kmesh
      real(kind=DP), intent(out) :: fbz_kmesh(3,nkmesh)
      integer, intent(in) :: nface,face(3,nface_max)

      ! local variables
      integer :: i,j,k
      real(kind=DP) :: g(3),q(3),gg,qg,g1(3),g2(3),g3(3),diff
      logical :: exists_in_fbz(nkmesh),trans_sym
      real(kind=DP), parameter :: eps = 1.d-8

      g1(1:3)=b1(1:3)*0.5d0
      g2(1:3)=b2(1:3)*0.5d0
      g3(1:3)=b3(1:3)*0.5d0

      ! restrict k-points within FBZ
      exists_in_fbz(1:nkmesh) = .true.
      do i=1,nface
        g(1:3) = face(1,i)*g1(1:3) + face(2,i)*g2(1:3) + face(3,i)*g3(1:3)
        gg = sum(g(1:3)*g(1:3))
        do j=1,nkmesh
          q(1:3) = kmesh(1,j)*b1(1:3) + kmesh(2,j)*b2(1:3) + kmesh(3,j)*b3(1:3)
          qg = sum(q(1:3)*g(1:3))
          if(qg > gg + eps) exists_in_fbz(j) = .false.
        end do
      end do

      ! traslational symmetry
      loop_1: do i=1,nkmesh
        if(.not.exists_in_fbz(i)) cycle loop_1
        loop_2: do j=i+1,nkmesh
          if(.not.exists_in_fbz(j)) cycle loop_2
          trans_sym = .true.
          do k=1,3
            diff =  kmesh(k,i) - kmesh(k,j)
            if(abs(dble(nint(diff))-diff) .gt. eps) trans_sym = .false.
          end do
          if(trans_sym) exists_in_fbz(j) = .false.
        end do loop_2
      end do loop_1

      nfbz_kmesh=0
      do i=1,nkmesh
        if(exists_in_fbz(i)) then
          nfbz_kmesh=nfbz_kmesh+1
          fbz_kmesh(1:3,nfbz_kmesh)=kmesh(1:3,i)
        end if
      end do

    end subroutine first_bz

    subroutine first_bz_segment(b1,b2,b3,nkmesh,kmesh,nfbz_kmesh,fbz_kmesh &
         & ,nface,face)
      implicit none

      real(kind=DP), intent(in) :: b1(3),b2(3),b3(3)
      integer, intent(in) :: nkmesh
      real(kind=DP), intent(in) :: kmesh(3,nkmesh)
      integer, intent(out) :: nfbz_kmesh
      real(kind=DP), intent(out) :: fbz_kmesh(3,nkmesh)
      integer, intent(in) :: nface,face(3,nface_max)

      ! local variables
      integer :: i,j,k
      real(kind=DP) :: g(3),q(3),gg,qg,g1(3),g2(3),g3(3),diff
      logical :: trans_sym
      real(kind=DP), parameter :: eps = 1.d-8
      logical, allocatable :: exists_in_fbz(:)
      real(kind=DP) :: eps2 = 1.0D-3

      integer :: ndiv_x, ndiv_y, ndiv_z
      integer :: ix, iy, iz, nx, ny, nz, j1, m1, mx, my, mz, jx, jy, jz
      integer :: num1, num2, num3, range_x2, range_y2, range_z2
      integer :: range_xs, range_xe, range_ys, range_ye, range_zs, range_ze
      real(kind=DP) :: resol_x, resol_y, resol_z, v(3)
      real(kind=DP) :: dx, dy, dz

      type t_kpt_segment
         integer :: nkpt
         integer, allocatable :: list(:,:)
         integer :: nkpt_reg
         integer, allocatable :: list_reg(:,:)
      end type t_kpt_segment
      type(t_kpt_segment), allocatable :: kpt_segment(:,:,:)

      if ( kpt_save_memory_mode == 1 ) then
         if ( npes > 1 .and. mype /= 0 ) goto 1000
      endif

      allocate( exists_in_fbz(nkmesh) )

      g1(1:3)=b1(1:3)*0.5d0
      g2(1:3)=b2(1:3)*0.5d0
      g3(1:3)=b3(1:3)*0.5d0

      ! restrict k-points within FBZ
      exists_in_fbz(1:nkmesh) = .true.
      do i=1,nface
        g(1:3) = face(1,i)*g1(1:3) + face(2,i)*g2(1:3) + face(3,i)*g3(1:3)
        gg = sum(g(1:3)*g(1:3))
        do j=1,nkmesh
          q(1:3) = kmesh(1,j)*b1(1:3) + kmesh(2,j)*b2(1:3) + kmesh(3,j)*b3(1:3)
          qg = sum(q(1:3)*g(1:3))
          if(qg > gg + eps) exists_in_fbz(j) = .false.
        end do
      end do

      resol_x = 1.0d0 /dble(ndiv_segment(1))
      resol_y = 1.0d0 /dble(ndiv_segment(2))
      resol_z = 1.0d0 /dble(ndiv_segment(3))

      ndiv_x = nint( ( 1.0d0 +kshift(1) ) /resol_x )
      ndiv_y = nint( ( 1.0d0 +kshift(2) ) /resol_y )
      ndiv_z = nint( ( 1.0d0 +kshift(3) ) /resol_z )
!      ndiv_x = nint( ( 2.0d0 +kshift(1) ) /resol_x )
!      ndiv_y = nint( ( 2.0d0 +kshift(2) ) /resol_y )
!      ndiv_z = nint( ( 2.0d0 +kshift(3) ) /resol_z )

      allocate( kpt_segment( -ndiv_x:ndiv_x, -ndiv_y:ndiv_y, -ndiv_z:ndiv_z ) )
      kpt_segment(:,:,:)%nkpt = 0

      do i=1,nkmesh
         v(1:3) = kmesh(1:3,i)
         ix = nint( v(1) /resol_x )
         iy = nint( v(2) /resol_y )
         iz = nint( v(3) /resol_z )
         kpt_segment(ix,iy,iz)%nkpt = kpt_segment(ix,iy,iz)%nkpt +1
      end do
      Do ix=-ndiv_x, ndiv_x
         Do iy=-ndiv_y, ndiv_y
            Do iz=-ndiv_z, ndiv_z
               num1 = kpt_segment(ix,iy,iz)%nkpt
               if ( num1 > 0 ) allocate( kpt_segment(ix,iy,iz)%list(num1,1) )
            End Do
         ENd Do
      end Do
      kpt_segment(:,:,:)%nkpt = 0

      do i=1, nkmesh
         v(1:3) = kmesh(1:3,i)
         ix = nint( v(1) /resol_x )
         iy = nint( v(2) /resol_y )
         iz = nint( v(3) /resol_z )
         kpt_segment(ix,iy,iz)%nkpt = kpt_segment(ix,iy,iz)%nkpt +1
         num1 = kpt_segment(ix,iy,iz)%nkpt
         kpt_segment(ix,iy,iz)%list(num1,1) = i
      End do

      range_x2 = 1;      range_y2 = 1;      range_z2 = 1
!      range_x2 = 2;      range_y2 = 2;      range_z2 = 2

      ! traslational symmetry
      loop_1: do i=1,nkmesh
         if (.not.exists_in_fbz(i)) cycle loop_1
         v(1:3) = kmesh(1:3,i)
         ix = nint( v(1) /resol_x )
         iy = nint( v(2) /resol_y )
         iz = nint( v(3) /resol_z )

         dx = v(1) -ix *resol_x
         dy = v(2) -iy *resol_y
         dz = v(3) -iz *resol_z
!
         range_xs = 0;  range_ys = 0;    range_zs = 0
         range_xe = 0;  range_ye = 0;    range_ze = 0
         if ( dx < -resol_x /2. +eps2 ) range_xs = -1
         if ( dx >  resol_x /2. -eps2 ) range_xe =  1
         if ( dy < -resol_y /2. +eps2 ) range_ys = -1
         if ( dy >  resol_y /2. -eps2 ) range_ye =  1
         if ( dz < -resol_z /2. +eps2 ) range_zs = -1
         if ( dz >  resol_z /2. -eps2 ) range_ze =  1
!
         Do mx=range_xs, range_xe
            Do my=range_ys, range_ye
               Do mz=range_zs, range_ze
                  jx = ix +mx;   jy = iy +my;   jz = iz +mz
                  do nx=-range_x2, range_x2
                     do ny=-range_y2, range_y2
                        do nz=-range_z2, range_z2
                           num1 = jx +ndiv_segment(1) *nx
                           num2 = jy +ndiv_segment(2) *ny
                           num3 = jz +ndiv_segment(3) *nz
                           if ( num1 < -ndiv_x .or. num1 > ndiv_x ) cycle
                           if ( num2 < -ndiv_y .or. num2 > ndiv_y ) cycle
                           if ( num3 < -ndiv_z .or. num3 > ndiv_z ) cycle

                           loop_2: do j1=1, kpt_segment(num1,num2,num3)%nkpt
                              m1 = kpt_segment(num1,num2,num3)%list(j1,1)
                              if ( m1 <= i ) cycle
                              if (.not.exists_in_fbz(m1)) cycle loop_2

                              trans_sym = .true.
                              do k=1,3
                                 diff = kmesh(k,i) - kmesh(k,m1)
                                 if ( abs( dble( nint(diff)) -diff ) .gt. eps) &
                                      &        trans_sym = .false.
                              end do
                              if(trans_sym) exists_in_fbz(m1) = .false.
                           end do loop_2
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do loop_1

      deallocate( kpt_segment )

      nfbz_kmesh=0
      do i=1,nkmesh
        if(exists_in_fbz(i)) then
          nfbz_kmesh=nfbz_kmesh+1
          fbz_kmesh(1:3,nfbz_kmesh)=kmesh(1:3,i)
        end if
      end do
      deallocate( exists_in_fbz )

1000  continue
      if ( kpt_save_memory_mode == 1 ) then
         if ( npes > 1 ) then
            call mpi_bcast( nfbz_kmesh, 1, mpi_integer, 0, MPI_CommGroup, ierr )
            call mpi_bcast( fbz_kmesh, 3*nfbz_kmesh, mpi_double_precision, &
                 &          0, MPI_CommGroup, ierr )
         endif
      endif

    end subroutine first_bz_segment

    subroutine fbz_face(b1,b2,b3,nface,face)
      implicit none

      real(kind=DP), intent(in) :: b1(3),b2(3),b3(3)
      integer, intent(out) :: nface
      integer, intent(out) :: face(3,nface_max)

      ! local variables
      integer :: i,j,k,n
      real(kind=DP) :: gg_ref,gg
      real(kind=DP) :: g1(3),g2(3),g3(3)
      real(kind=DP) :: gface(3,nface_max)
      integer :: iface(3,nface_max)
      logical :: exists_in_fbz(nface_max)
      real(kind=DP), parameter :: lambda = 1.d0+1.d-8

      g1(1:3)=b1(1:3)*0.5d0
      g2(1:3)=b2(1:3)*0.5d0
      g3(1:3)=b3(1:3)*0.5d0

      n=0
      do i=-2,2
        do j=-2,2
          do k=-2,2
            if(i == 0 .and. j == 0 .and. k == 0) cycle
            n=n+1
            gface(1:3,n) = i*g1(1:3) + j*g2(1:3) + k*g3(1:3)
            iface(1,n)=i
            iface(2,n)=j
            iface(3,n)=k
          end do
        end do
      end do

      exists_in_fbz(1:nface_max) = .true.
      do i=1,n
        if(exists_in_fbz(i)) then
          gg_ref = gface(1,i)**2 + gface(2,i)**2 + gface(3,i)**2
          do j=1,n
            if(i /= j) then
              gg = gface(1,i)*gface(1,j) + gface(2,i)*gface(2,j) + gface(3,i)*gface(3,j)
              if(gg*lambda > gg_ref) exists_in_fbz(j) = .false.
            end if
          end do
        end if
      end do

      nface=0
      do i=1,n
        if(exists_in_fbz(i)) then
          nface=nface+1
          face(1:3,nface)=iface(1:3,i)
        end if
      end do

      ! debug
      !if(printable .and. printlevel .ge. 3) then
      !  write(output,*) 'nface=',nface
      !  write(output,'(4(1x,a3))') ' no',' n1',' n2',' n3'
      !  do i=1,nface
      !    write(output,'(4(1x,i3))') i,face(1:3,i)
      !  end do
      !end if
      !stop 'debug fbz_face'
      ! end debug

    end subroutine fbz_face

    subroutine get_grot(a1,a2,a3,b1,b2,b3,nsym,rot,grot)
      implicit none

      real(kind=DP), intent(in) :: a1(3),a2(3),a3(3)
      real(kind=DP), intent(in) :: b1(3),b2(3),b3(3)
      integer, intent(in) :: nsym
      real(kind=DP), intent(in) :: rot(3,3,nsym)
      real(kind=DP), intent(out) :: grot(3,3,nsym)

      ! local variables
      integer :: n,i,j
      real(kind=DP) :: amat(3,3),bmat(3,3)

      amat(1:3,1) = a1(1:3)
      amat(1:3,2) = a2(1:3)
      amat(1:3,3) = a3(1:3)
      bmat(1:3,1) = b1(1:3)
      bmat(1:3,2) = b2(1:3)
      bmat(1:3,3) = b3(1:3)

      do n=1,nsym
        do j=1,3
          do i=1,3
            grot(i,j,n) = element(amat(1,i),bmat(1,j),rot(1,1,n))
          end do
        end do
      end do


! === KT_add === 2014/08/14 & 13.2S
      do n=1, nsym
         if ( noncol .or. sw_use_magnetic_symmetry == ON ) then
            if(magmom_dir_inversion_opr_flag(n) == -1 ) then
               grot(:,:,n) = -grot(:,:,n)
            endif
         endif
      end do
! ============== 2014/08/14 & 13.2S

    end subroutine get_grot

    function element(a,b,rot) result(value)
      implicit none

      real(kind=DP) :: value

      real(kind=DP), intent(in) :: a(3),b(3),rot(3,3)

      ! local variables
      integer :: i,j

      value = 0.d0

      do j=1,3
        do i=1,3
          value = value + a(i)*rot(i,j)*b(j)
        end do
      end do

    end function element

    subroutine gen_spk( nsym, grot, nfbz_kmesh, fbz_kmesh, nkpoint, nkpmax,  &
         &              kpoint, weight )
      implicit none

      integer, intent(in) :: nsym
      real(kind=DP), intent(in) :: grot(3,3,nsym)
      integer, intent(in) :: nfbz_kmesh
      real(kind=DP), intent(in) :: fbz_kmesh(3,nfbz_kmesh)
      integer, intent(inout) :: nkpoint
      integer, intent(in) :: nkpmax
      real(kind=DP), intent(inout) :: kpoint(3,nkpmax), weight(nkpmax)

      ! local variables
      integer :: i,j,n
      integer :: nequiv(nkpmax)
      real(kind=DP) :: g(3),q(3)
!      real(kind=DP), parameter :: eps = 1.d-3
      real(kind=DP), parameter :: eps = 1.d-5
      logical :: not_included(nfbz_kmesh)

      not_included(1:nfbz_kmesh) = .true.

      nkpoint = 0
      do i=1,nfbz_kmesh
        if(not_included(i)) then
          nkpoint = nkpoint+1
          if(nkpoint > nkpmax) call phase_error_with_msg(nfout,'gen_spk: nkpoint > nkpmax',__LINE__,__FILE__)
          g(1:3) = fbz_kmesh(1:3,i)
          kpoint(1:3,nkpoint) = g(1:3)
          nequiv(nkpoint) = 1

          do j=i+1,nfbz_kmesh
            do n=1,nsym
              q(1:3) = grot(1:3,1,n)*fbz_kmesh(1,j) + grot(1:3,2,n)*fbz_kmesh(2,j) + grot(1:3,3,n)*fbz_kmesh(3,j)

              if(sum((g(1:3)-q(1:3))**2) < eps .and. not_included(j)) then
                nequiv(nkpoint) = nequiv(nkpoint)+1
                not_included(j) = .false.
              end if

              ! debug
              !write(output,'(i2,3(1x,f10.5),2x,3(1x,f10.5),2x,l1)') n,g(1:3),q(1:3),not_included(j)
              ! end debug
            end do
          end do
        end if
      end do

      ! debug
      !do i=1,nfbz_kmesh
      !  write(output,'(3(1x,f10.5),1x,l1)') fbz_kmesh(1:3,i), not_included(i)
      !end do
      ! end debug

      n=sum(nequiv(1:nkpoint))
      do i=1,nkpoint
        weight(i) = dble(nequiv(i))/dble(n)
      end do

    end subroutine gen_spk

! === KT_add ==== 2014/09/30
    subroutine gen_spk2( nsym, grot, nfbz_kmesh, fbz_kmesh, nkpoint, nkpmax,  &
         &               kpoint, weight, kpoint_index_in_ibz )
      implicit none

      integer, intent(in) :: nsym
      real(kind=DP), intent(in) :: grot(3,3,nsym)
      integer, intent(in) :: nfbz_kmesh
      real(kind=DP), intent(in) :: fbz_kmesh(3,nfbz_kmesh)
      integer, intent(inout) :: nkpoint
      integer, intent(in) :: nkpmax
      real(kind=DP), intent(inout) :: kpoint(3,nkpmax), weight(nkpmax)
      integer, optional, intent(inout) :: kpoint_index_in_ibz(nfbz_kmesh)

      ! local variables
      integer :: i,j,n
      integer :: nequiv(nkpmax)
      real(kind=DP) :: g(3), q(3), vec_tmp(3)
!      real(kind=DP), parameter :: eps = 1.d-3
      real(kind=DP), parameter :: eps = 1.d-5
      logical :: not_included(nfbz_kmesh)

      not_included(1:nfbz_kmesh) = .true.

      nkpoint = 0
      do i=1,nfbz_kmesh
         if(not_included(i)) then
            nkpoint = nkpoint+1
            if(nkpoint > nkpmax) call phase_error_with_msg(nfout,'gen_spk: nkpoint > nkpmax',__LINE__,__FILE__)
            g(1:3) = fbz_kmesh(1:3,i)
            kpoint(1:3,nkpoint) = g(1:3)
            nequiv(nkpoint) = 1
! ============ KT_add === 2014/09/30
            if ( present(kpoint_index_in_ibz) ) kpoint_index_in_ibz(i) = nkpoint
! ======================= 2014/09/30

            do j=i+1,nfbz_kmesh
               do n=1,nsym
                  q(1:3) = grot(1:3,1,n)*fbz_kmesh(1,j) &
                       &    +grot(1:3,2,n)*fbz_kmesh(2,j) &
                       &    +grot(1:3,3,n)*fbz_kmesh(3,j)

                  vec_tmp(1:3) = g(1:3) -q(1:3)
                  vec_tmp(1:3) = vec_tmp(1:3) -nint( vec_tmp(1:3) )

                  if ( sum( vec_tmp(1:3)**2 ) < eps .and. not_included(j) ) then
                     nequiv(nkpoint) = nequiv(nkpoint)+1
                     not_included(j) = .false.
! ============ KT_add === 2014/09/30
                     if ( present(kpoint_index_in_ibz) ) &
                          &        kpoint_index_in_ibz( j ) = nkpoint
! ======================= 2014/09/30
                  end if
               end do
            end do
         end if
      end do

      ! debug
      !do i=1,nfbz_kmesh
      !  write(output,'(3(1x,f10.5),1x,l1)') fbz_kmesh(1:3,i), not_included(i)
      !end do
      ! end debug

      n=sum(nequiv(1:nkpoint))
      do i=1,nkpoint
        weight(i) = dble(nequiv(i))/dble(n)
      end do

    end subroutine gen_spk2
! =============== 2014/09/30

    subroutine gen_spk2_segment( nsym, grot, nfbz_kmesh, fbz_kmesh, nkpoint, nkpmax,  &
         &                       kpoint, weight, kpoint_index_in_ibz )
      implicit none

      integer, intent(in) :: nsym
      real(kind=DP), intent(in) :: grot(3,3,nsym)
      integer, intent(in) :: nfbz_kmesh
      real(kind=DP), intent(in) :: fbz_kmesh(3,nfbz_kmesh)
      integer, intent(inout) :: nkpoint
      integer, intent(in) :: nkpmax
      real(kind=DP), intent(inout) :: kpoint(3,nkpmax), weight(nkpmax)
      integer, optional, intent(inout) :: kpoint_index_in_ibz(nfbz_kmesh)

      ! local variables
      integer :: i,n
      real(kind=DP) :: g(3), q(3), vec_tmp(3)
!      real(kind=DP), parameter :: eps = 1.d-3
      real(kind=DP), parameter :: eps = 1.d-5
      real(kind=DP) :: eps2 = 1.0D-3

      integer, allocatable :: nequiv(:)
      logical, allocatable :: not_included(:)

      integer :: ndiv_x, ndiv_y, ndiv_z
      integer :: ix, iy, iz, nx, ny, nz, j1, m1, m2, mx, my, mz, jx, jy, jz
      integer :: num1, num2, num3, range_x2, range_y2, range_z2
      integer :: range_xs, range_xe, range_ys, range_ye, range_zs, range_ze
      real(kind=DP) :: resol_x, resol_y, resol_z
      real(kind=DP) :: dx, dy, dz

      type t_kpt_segment
         integer :: nkpt
         integer, allocatable :: list(:,:)
         integer :: nkpt_reg
         integer, allocatable :: list_reg(:,:)
      end type t_kpt_segment
      type(t_kpt_segment), allocatable :: kpt_segment(:,:,:)

      if ( kpt_save_memory_mode == 1 ) then
         if ( npes > 1 .and. mype /= 0 ) goto 1000
      endif

      allocate( nequiv(nkpmax) );  nequiv = 0
      allocate( not_included(nfbz_kmesh) );   not_included(1:nfbz_kmesh) = .true.

      resol_x = 1.0d0 /dble(ndiv_segment(1))
      resol_y = 1.0d0 /dble(ndiv_segment(2))
      resol_z = 1.0d0 /dble(ndiv_segment(3))

      ndiv_x = nint( ( 1.0d0 +kshift(1) ) /resol_x )
      ndiv_y = nint( ( 1.0d0 +kshift(2) ) /resol_y )
      ndiv_z = nint( ( 1.0d0 +kshift(3) ) /resol_z )

      allocate( kpt_segment( -ndiv_x:ndiv_x, -ndiv_y:ndiv_y, -ndiv_z:ndiv_z ) )
      kpt_segment(:,:,:)%nkpt = 0

      do i=1, nfbz_kmesh
         do n=1, nsym
            q(1:3) = grot(1:3,1,n)*fbz_kmesh(1,i) +grot(1:3,2,n)*fbz_kmesh(2,i) &
                 & + grot(1:3,3,n)*fbz_kmesh(3,i)
            ix = nint( q(1) /resol_x )
            iy = nint( q(2) /resol_y )
            iz = nint( q(3) /resol_z )
            kpt_segment(ix,iy,iz)%nkpt = kpt_segment(ix,iy,iz)%nkpt +1
         end do
      end do
      Do ix=-ndiv_x, ndiv_x
         Do iy=-ndiv_y, ndiv_y
            Do iz=-ndiv_z, ndiv_z
               num1 = kpt_segment(ix,iy,iz)%nkpt
               if ( num1 > 0 ) allocate( kpt_segment(ix,iy,iz)%list(num1,2) )
            End Do
         ENd Do
      end Do
      kpt_segment(:,:,:)%nkpt = 0

      do i=1, nfbz_kmesh
         do n=1, nsym
            q(1:3) = grot(1:3,1,n)*fbz_kmesh(1,i) +grot(1:3,2,n)*fbz_kmesh(2,i) &
                 & + grot(1:3,3,n)*fbz_kmesh(3,i)
            ix = nint( q(1) /resol_x )
            iy = nint( q(2) /resol_y )
            iz = nint( q(3) /resol_z )
            kpt_segment(ix,iy,iz)%nkpt = kpt_segment(ix,iy,iz)%nkpt +1
            num1 = kpt_segment(ix,iy,iz)%nkpt
            kpt_segment(ix,iy,iz)%list(num1,1) = i
            kpt_segment(ix,iy,iz)%list(num1,2) = n
         End do
      End do

      range_x2 = 1;  range_y2 = 1;  range_z2 = 1
!      range_x2 = 2;  range_y2 = 2;  range_z2 = 2

      nkpoint = 0
      do i=1,nfbz_kmesh
         if (not_included(i)) then
            nkpoint = nkpoint+1
            if(nkpoint > nkpmax) call phase_error_with_msg(nfout,'gen_spk: nkpoint > nkpmax',__LINE__,__FILE__)
            g(1:3) = fbz_kmesh(1:3,i)
            kpoint(1:3,nkpoint) = g(1:3)
            nequiv(nkpoint) = 1
            if ( present(kpoint_index_in_ibz) ) kpoint_index_in_ibz(i) = nkpoint

            ix = nint( g(1) /resol_x )
            iy = nint( g(2) /resol_y )
            iz = nint( g(3) /resol_z )

            dx = g(1) -ix *resol_x
            dy = g(2) -iy *resol_y
            dz = g(3) -iz *resol_z
!
            range_xs = 0;  range_ys = 0;    range_zs = 0
            range_xe = 0;  range_ye = 0;    range_ze = 0
            if ( dx < -resol_x /2. +eps2 ) range_xs = -1
            if ( dx >  resol_x /2. -eps2 ) range_xe =  1
            if ( dy < -resol_y /2. +eps2 ) range_ys = -1
            if ( dy >  resol_y /2. -eps2 ) range_ye =  1
            if ( dz < -resol_z /2. +eps2 ) range_zs = -1
            if ( dz >  resol_z /2. -eps2 ) range_ze =  1
!
            Do mx=range_xs, range_xe
               Do my=range_ys, range_ye
                  Do mz=range_zs, range_ze
                     jx = ix +mx;   jy = iy +my;   jz = iz +mz

                     Do nx=-range_x2, range_x2
                        Do ny=-range_y2, range_y2
                           Do nz=-range_z2, range_z2
                              num1 = jx +ndiv_segment(1) *nx
                              num2 = jy +ndiv_segment(2) *ny
                              num3 = jz +ndiv_segment(3) *nz
                              if ( num1 < -ndiv_x .or. num1 > ndiv_x ) cycle
                              if ( num2 < -ndiv_y .or. num2 > ndiv_y ) cycle
                              if ( num3 < -ndiv_z .or. num3 > ndiv_z ) cycle

                              do j1=1, kpt_segment(num1,num2,num3)%nkpt
                                 m1 = kpt_segment(num1,num2,num3)%list(j1,1)
                                 m2 = kpt_segment(num1,num2,num3)%list(j1,2)
                                 if ( m1 <= i ) cycle

                                 q(1:3) = grot(1:3,1,m2)*fbz_kmesh(1,m1) &
                                      &    +grot(1:3,2,m2)*fbz_kmesh(2,m1) &
                                      &    +grot(1:3,3,m2)*fbz_kmesh(3,m1)
                                 vec_tmp(1:3) = g(1:3) -q(1:3)
                                 vec_tmp(1:3) = vec_tmp(1:3) -nint( vec_tmp(1:3) )

                                 if ( sum( vec_tmp(1:3)**2 ) < eps &
                                      &   .and. not_included(m1) ) then
                                    nequiv(nkpoint) = nequiv(nkpoint)+1
                                    not_included(m1) = .false.
                                    if ( present(kpoint_index_in_ibz) ) &
                                         &        kpoint_index_in_ibz(m1) = nkpoint
                                 end if
                              end do
                           End Do
                        end do
                     end Do
                  end Do
               end Do
            end do
         end if
      end do
      deallocate( kpt_segment )

      n=sum(nequiv(1:nkpoint))
      do i=1,nkpoint
        weight(i) = dble(nequiv(i))/dble(n)
      end do
      deallocate( nequiv );   deallocate( not_included )

1000  continue

      if ( kpt_save_memory_mode == 1 ) then
         if ( npes > 1 ) then
            call mpi_bcast( nkpoint, 1, mpi_integer, 0, MPI_CommGroup, ierr )
            call mpi_bcast( kpoint, 3*nkpoint, mpi_double_precision, &
                 &          0, MPI_CommGroup, ierr )
            call mpi_bcast( weight, nkpoint, mpi_double_precision, &
                 &          0, MPI_CommGroup, ierr )
            if ( present(kpoint_index_in_ibz) ) then
               call mpi_bcast( kpoint_index_in_ibz, nfbz_kmesh, mpi_integer, &
                    &          0, MPI_CommGroup, ierr )
            endif
         endif
      endif

    end subroutine gen_spk2_segment

    subroutine positive(nsym,grot,nkpoint,kpoint)
      implicit none

      integer, intent(in) :: nsym
      real(kind=DP), intent(in) :: grot(3,3,nsym)
      integer, intent(in) :: nkpoint
      real(kind=DP), intent(inout) :: kpoint(3,nkpoint)

      ! local variables
      integer :: i,n
      real(kind=DP) :: g(3)

      do i=1,nkpoint
        if(kpoint(1,i) /= 0.d0 .or. kpoint(2,i) /= 0.d0 .or. kpoint(3,i) /= 0.d0) then
          do n=1,nsym
            g(1:3) = grot(1:3,1,n)*kpoint(1,i) + grot(1:3,2,n)*kpoint(2,i) + grot(1:3,3,n)*kpoint(3,i)
            if(g(1) >= 0.d0 .and. g(2) >= 0.d0 .and. g(3) >= 0.d0) then
              kpoint(1:3,i)=g(1:3)
              exit
            end if
          end do
        end if
      end do

    end subroutine positive

! === ASMS ====
    subroutine positive2( nsym, grot, nkpoint, kpoint, nkpoint_full, kpoint_full )
      implicit none

      integer, intent(in) :: nsym
      real(kind=DP), intent(in) :: grot(3,3,nsym)
      integer, intent(in) :: nkpoint, nkpoint_full
      real(kind=DP), intent(inout) :: kpoint(3,nkpoint)
      real(kind=DP), intent(in) :: kpoint_full(3,nkpoint_full)

      integer :: iopr, ik, ikbz
      real(kind=DP) :: v0(3), dk(3)

      Do ik=1, nkpoint
         Do iopr=1, nsym
            v0 = matmul( grot(:,:,iopr), kpoint(:,ik) )

!            write(nfout,*) "ik iopr ", ik, iopr
!            write(nfout,'(A,3F20.10)') "kp ", kpoint(1:3,ik)
!            write(nfout,'(A,3F20.10)') "v0 ", v0(1:3)

            if ( v0(1) >= 0.0d0 .and. v0(2) >= 0.0d0 .and. v0(3) >= 0.0d0 ) then
               LOOP_S: Do ikbz=1, nkpoint_full
                  dk(1:3) = v0(1:3) -kpoint_full(1:3,ikbz)

                  if( abs(dk(1)-nint(dk(1)))<DELTA07  &
                       &      .and. abs(dk(2)-nint(dk(2)))<DELTA07  &
                       &      .and. abs(dk(3)-nint(dk(3)))<DELTA07 ) then
                     kpoint(:,ik) = v0(:)
!                     write(nfout,'(A,3F20.10)') "kp new ", kpoint(1:3,ik)
                     exit LOOP_S
                  endif
               End do LOOP_S
            endif
         End Do
      End Do

    end subroutine positive2
! === ASMS ====

    subroutine accuracy(a1,a2,a3,nsym,grot,nkpoint,kpoint,weight)
      implicit none

      real(kind=DP), intent(in) :: a1(3),a2(3),a3(3)
      integer, intent(in) :: nsym
      real(kind=DP), intent(in) :: grot(3,3,nsym)
      integer, intent(in) :: nkpoint
      real(kind=DP), intent(in) :: kpoint(3,nkpoint),weight(nkpoint)

      ! local variables
      real(kind=DP) :: rgrid(0:3,ngrid)
      type(shell), allocatable :: rshell(:)
      type(shell) :: rshtmp
      integer, allocatable :: id_sub(:)
      integer :: i,j,k,nsub
      real(kind=DP) :: sm
      integer :: n,istart
      integer :: none_zero
      real(kind=DP) :: sm_nz,length_nz


      if(printable .and. printlevel .ge. 0) &
      &  write(output,'(1x,"number of k-points in irreducible BZ = ",i8)') nkpoint


      call gen_rgrid(a1,a2,a3,rgrid)

      ! debug
      !do i=1,ngrid
      !  write(output,'(3(1x,i4),1x,f10.5)') nint(rgrid(1:3,i)), rgrid(0,i)
      !end do
      !stop 'debug rgrid'
      ! end debug

      istart=2
      n=0
      shell_loop : do
        if( n+1 > nshell ) exit shell_loop
        call get_shell(rgrid,rshtmp,istart)
        allocate(id_sub(rshtmp%ndim))
        call get_subshell(nsym,grot,rshtmp,id_sub,nsub)

        ! debug
        !write(output,*) 'n=',n+1
        !do i=1,rshtmp%ndim
        !  write(output,'(3(1x,i5),1x,i2)') nint(rshtmp%rindex(1:3,i)),id_sub(i)
        !end do
        !stop 'debug rshtmp'
        ! end debug
        allocate(rshell(nsub))
        do i=1,nsub
          allocate(rshell(i)%rindex(3,48))
          rshell(i)%ndim = 0
          do j=1,rshtmp%ndim
            if(id_sub(j) == i) then
              rshell(i)%ndim = rshell(i)%ndim+1
              rshell(i)%rindex(1:3,rshell(i)%ndim) = rshtmp%rindex(1:3,j)
            end if
          end do
          rshell(i)%length=rshtmp%length
        end do
        do i=1,nsub
          if( n+1 > nshell ) exit shell_loop
          n=n+1
          ! debug
          !write(output,'(i5,1x,f10.5,4(1x,i5))') n,rshell(i)%length,nint(rshell(i)%rindex(1:3,1)),rshell(i)%ndim
          ! end debug
          do j=1,rshell(i)%ndim
            if(nint(rshell(i)%rindex(1,j)) == ng ) call phase_error_with_msg(nfout,'accuracy: rindex >= ng',__LINE__,__FILE__)
            if(nint(rshell(i)%rindex(2,j)) == ng ) call phase_error_with_msg(nfout,'accuracy: rindex >= ng',__LINE__,__FILE__)
            if(nint(rshell(i)%rindex(3,j)) == ng ) call phase_error_with_msg(nfout,'accuracy: rindex >= ng',__LINE__,__FILE__)
          end do
          sm=0.d0
          do k=1,nkpoint
             sm = sm + weight(k)*amk(rshell(i)%rindex,rshell(i)%ndim,kpoint(1,k))
          end do

          ! debug
          !write(output,'(i5,2(1x,f10.5))') n, rshell(i)%length, sm
          ! end debug

          if(abs(sm) > 1.d-8) then
            none_zero = n
            sm_nz = sm/dble(2*nkpoint)
            length_nz = rshell(i)%length
            do j=1,nsub
              deallocate(rshell(j)%rindex)
            end do
            deallocate(rshtmp%rindex,rshell,id_sub)
            exit shell_loop
          end if
        end do
        do i=1,nsub
          deallocate(rshell(i)%rindex)
        end do
        deallocate(rshtmp%rindex,rshell,id_sub)
      end do shell_loop

      if( n == nshell ) call phase_error_with_msg(nfout,'accuracy: n == nshell, accuracy check failure',__LINE__,__FILE__)

      if(printable .and. printlevel .ge. 0 ) then
        write(output,'(1x,"Index of the first none zero shell = ",i8)') none_zero
        write(output,'(1x,"|R| of the first none zero shell   = ",f10.5)') length_nz
        write(output,'(1x,"phase sum = ",f10.5)') sm_nz
        write(output,'(1x,"efficiency = ",f8.2)') dble(none_zero)/dble(nkpoint)
      end if

      if(printable .and. printlevel .ge. 1) then
        write(output,'(5x,4(1x,a10))') 'k1','k2','k3','weight'
        do i=1,nkpoint
          write(output,'(i8,4(1x,f10.5))') i,kpoint(1:3,i),weight(i)
        end do
      end if

    end subroutine accuracy

    subroutine gen_rgrid(a1,a2,a3,rgrid)
      implicit none

      real(kind=DP), intent(in) :: a1(3),a2(3),a3(3)
      real(kind=DP), intent(out) :: rgrid(0:3,ngrid)

      ! local variables
      integer :: n,i1,i2,i3
      real(kind=DP) :: x(3)

      n=0
      loop_1: do i1=-ng,ng
        loop_2: do i2=-ng,ng
          loop_3: do i3=-ng,ng
            n=n+1
            x(1:3) = i1*a1(1:3) + i2*a2(1:3) + i3*a3(1:3)
            rgrid(0,n) = sqrt(sum(x(1:3)**2))
            rgrid(1,n) = dble(i1)
            rgrid(2,n) = dble(i2)
            rgrid(3,n) = dble(i3)
          end do loop_3
        end do loop_2
      end do loop_1

      call heap_sort(4,ngrid,rgrid,1)

    end subroutine gen_rgrid

    subroutine get_shell(rgrid,rshell,istart)
      implicit none

      real(kind=DP), intent(in) :: rgrid(0:3,ngrid)
      type(shell), intent(out) :: rshell
      integer, intent(inout) :: istart

      ! local variables
      real(kind=DP) :: eps = 1.d-10
      integer :: i,ndim

      ndim=0
      do
        ndim=ndim+1
        if(istart+ndim>ngrid) exit
        if(abs(rgrid(0,istart+ndim)-rgrid(0,istart+ndim-1)) > eps) exit
      end do

      allocate(rshell%rindex(3,ndim))
      do i=1,ndim
        rshell%rindex(1:3,i) = rgrid(1:3,istart+i-1)
      end do
      rshell%length = rgrid(0,istart)
      rshell%ndim = ndim

      istart=istart+ndim

    end subroutine get_shell

    subroutine get_subshell(nsym,grot,rshell,id_sub,nsub)
      implicit none

      integer, intent(in) :: nsym
      real(kind=DP), intent(in) :: grot(3,3,nsym)
      type(shell), intent(in) :: rshell
      integer, intent(out) :: id_sub(rshell%ndim)
      integer, intent(out) :: nsub

      ! local variables
      integer :: i,j,n
      real(kind=DP) :: eps=1.d-8
      real(kind=DP) :: v(3)

      id_sub(1:rshell%ndim)=0
      nsub=0
      do i=1,rshell%ndim
        if(id_sub(i) /= 0) cycle
        nsub=nsub+1
        id_sub(i)=nsub
        do n=1,nsym
          v(1:3) =   rshell%rindex(1,i)*grot(1,1:3,n) &
                 & + rshell%rindex(2,i)*grot(2,1:3,n) &
                 & + rshell%rindex(3,i)*grot(3,1:3,n)
          do j=1,rshell%ndim
            if(abs(v(1)-rshell%rindex(1,j)) < eps .and. &
               abs(v(2)-rshell%rindex(2,j)) < eps .and. &
               abs(v(3)-rshell%rindex(3,j)) < eps ) then
              id_sub(j)=id_sub(i)
            end if
          end do
        end do
      end do

    end subroutine get_subshell

    function amk(ri,ndim,kp) result(value)
      use m_Const_Parameters, only : PAI
      implicit none

      real(kind=DP) :: value

      integer, intent(in) :: ndim
      real(kind=DP), intent(in) :: ri(3,ndim),kp(3)

      ! local variables
      integer :: i
      real(kind=DP), parameter :: pi = PAI, pi2=2*pi

      value=0.d0
      do i=1,ndim
        value = value + cos(pi2*sum(kp(1:3)*ri(1:3,i)))
      end do
      value = value/dble(ndim)

    end function amk

  end subroutine special_kpoints

  subroutine gen_special_kpoints( paramset, nfout, printlevel, nkpmax, &
       &                          nfmatbp, nkpoint, vkxyz, qwgt, imag, itrs, &
       &                          nkpoint_full, vkxyz_full, to_ibz_from_fbz_for_kpoint )

    use m_Const_Parameters, only : CARTS, CRDTYP, DP, PAI
    use m_Crystal_Structure, only : altv, rltv, nopr, op

! ========================= added by K. Tagami =========================== 11.0P
    use m_CS_Magnetic,   only : nopr_before_sym_reduction, op_before_sym_reduction
! ======================================================================== 11.0P

    implicit none

    logical, intent(in)        :: paramset
    integer, intent(in)        :: nfout, printlevel, nkpmax
    integer, intent(in)        :: nfmatbp
    integer, intent(inout)     :: nkpoint
    ! if nkpoint == 0 then generating special k-points
    ! else ckecking accuracy of special k-points
    real(kind=DP), intent(out) :: vkxyz(nkpmax,3,CRDTYP),qwgt(nkpmax)
    integer, intent(in)        :: imag
    integer, intent(in)        :: itrs

! === KT_add === 2014/09/30
    integer, intent(inout)     :: nkpoint_full
    real(kind=DP), intent(out) :: vkxyz_full(nkpoint_full,3,CRDTYP)
    integer, intent(out) :: to_ibz_from_fbz_for_kpoint(nkpoint_full)
! ============== 2014/09/30

    ! local variables
    real(kind=DP), parameter :: pi = PAI, pi2=2*pi
    integer       :: i,nk,pl
    real(kind=DP), allocatable, dimension(:,:) :: trmat,trbp,trpb,mat1,mat2
    real(kind=DP) :: v(3)
    real(kind=DP) :: b1(3),b2(3),b3(3)
    real(kind=DP), allocatable :: kpoint(:,:), weight(:)

! === KT_add ==== 2014/09/30
    real(kind=DP), allocatable :: kpoint_full(:,:)
    integer,       allocatable :: kpoint_num_in_ibz(:)
! =============== 2014/09/30

    !!integer, parameter :: nkpmax0 = 10000
    integer :: nkpmax0
    logical :: magnetic
    logical :: use_trs

    if(imag/=PARA) then
       magnetic = .true.
    else
       magnetic = .false.
    end if
    if(itrs == 1) then
       use_trs = .true.
    else
       use_trs = .false.
    end if

    b1(1:3) = rltv(1:3,1)/pi2
    b2(1:3) = rltv(1:3,2)/pi2
    b3(1:3) = rltv(1:3,3)/pi2

    if(paramset) then
      pl = -1
    else
      pl = printlevel
    end if
    !write(nfout,*) 'paramset=',paramset

    if(.not.(allocated(kpoint).or.allocated(weight))) then
      nkpoint = 0 ! for generating special k-points
      nkpmax0 = 2*mp_index(1)*mp_index(2)*mp_index(3)
      ng= 2 * maxval(mp_index)
      ngrid=(2*ng+1)**3
      nshell = nkpmax0
      allocate(kpoint(3,nkpmax0));   allocate(weight(nkpmax0))
! === KT_add === 2014/09/30
      allocate(kpoint_full(3,nkpmax0))
      allocate(kpoint_num_in_ibz(nkpmax0))
! ============== 2014/09/30
    end if

! ==================================== modified by K. Tagami ============== 11.0P
!    call special_kpoints( pl, nfout, nopr, op, &
!    &                     b1, b2, b3, altv(1,1), altv(1,2), altv(1,3), &
!    &                     nkpoint, nkpmax0, kpoint, weight, magnetic, use_trs )
!
    if ( noncol .and. use_op_before_sym_reduction ) then
       call special_kpoints( pl, nfout, nopr, op, &
            &                b1, b2, b3, altv(1,1), altv(1,2), altv(1,3), &
            &                nkpoint, nkpmax0, kpoint, weight, magnetic, use_trs, &
            &                nkpoint_full, kpoint_full, kpoint_num_in_ibz, &
            &                nopr_before_sym_reduction, op_before_sym_reduction )
    else
       call special_kpoints( pl, nfout, nopr, op, &
            &                b1, b2, b3, altv(1,1), altv(1,2), altv(1,3), &
            &                nkpoint, nkpmax0, kpoint, weight, magnetic, use_trs, &
            &                nkpoint_full, kpoint_full, kpoint_num_in_ibz )
    endif
! ========================================================================== 11.0P

    ! exists in m_SpecicalKpoints

    ! debug
    ! stop 'debug b_SpecailKpoints when first calling'
    ! end debug


    if(paramset) then
      deallocate(kpoint,weight);
      deallocate(kpoint_full); deallocate(kpoint_num_in_ibz)
      return
    end if

    allocate(trmat(3,3))
    allocate(trbp(3,3))
    allocate(trpb(3,3))
    allocate(mat1(3,3))
    allocate(mat2(3,3))

    call get_trmat  !-(contained here) ->(trmat)

    do nk = 1, nkpoint
      v(1:3) = kpoint(1:3,nk)
      vkxyz(nk,1:3,CARTS) = matmul(trmat,v)
      qwgt(nk) = weight(nk)
    enddo

! ==== KT_add === 2014/09/30
    do nk = 1, nkpoint_full
      v(1:3) = kpoint_full(1:3,nk)
      vkxyz_full(nk,1:3,CARTS) = matmul(trmat,v)
      to_ibz_from_fbz_for_kpoint(nk) = kpoint_num_in_ibz(nk)
    end do
    deallocate( kpoint_full ); deallocate( kpoint_num_in_ibz )
! =============== 2014/09/30

    deallocate(kpoint,weight)

    deallocate(trmat); deallocate(trbp ); deallocate(trpb )
    deallocate(mat1 ); deallocate(mat2 )

    if(printable .and. printlevel >= 2) then
       write(nfout,'(/,"  << gen_special_kpoints >>")')
       write(nfout,*) ' !Total Generated Kpoints = ',nkpoint
       do nk = 1, nkpoint
         write(nfout,'(i6," ",3f12.6," : ",f12.6)') &
         & nk,(vkxyz(nk,i,CARTS),i=1,3),qwgt(nk)
       end do
    end if

  contains
    subroutine get_trmat
    !    make translation matrix  trpb (P -> B)
      integer  :: i
      logical :: open

      trbp = 0.d0
      do i = 1,3
        trbp(i,i) = 1.d0
      end do

      inquire(unit=nfmatbp, opened = open)
      if(open) then
        rewind nfmatbp
        do i = 1, 3
          read(nfmatbp,*,end=1,err=1) trbp(i,1),trbp(i,2),trbp(i,3)
        enddo
        goto 2
1       trbp = 0.d0
        do i = 1,3
          trbp(i,i) = 1.d0
        end do
2       continue
      end if

      call inver3n(3,trbp,trpb)

      mat1 = transpose(trpb)
      call inver3n(3,mat1,mat2)
      call matpr3(rltv,mat2,trmat)

    end subroutine get_trmat
  end subroutine gen_special_kpoints

  subroutine read_mpindex(input)
    implicit none

    integer, intent(in) :: input

    read(input,*) mp_index(1:3)
    read(input,*) kshift(1:3)

  end subroutine read_mpindex

  subroutine m_Kp_cp_kxyz_to_vkxyz(nk,kxyz)
    integer, intent(in) :: nk
    real(kind=DP),intent(in) :: kxyz(3,nk)
    integer :: i,ik

    if(.not.allocated(vkxyz)) then
       if(ipri_kp >= 1) write(nfout,'(" vkxyz is not allocated")')
       call phase_error_with_msg(nfout,'vkxyz is not allocated',__LINE__,__FILE__)
    end if

    if(ipri_kp>=1) write(nfout,'(" nk = ",i6," <<m_Kp_cp_kxyz_to_vkxyz>>")') nk
    do ik = 1, nk
       vkxyz(ik,1:3,BUCS) = kxyz(1:3,ik)
       if(ipri_kp>=1) write(nfout,'(" kxyz = ",3f10.6)') kxyz(1:3,ik)
    end do

    do ik = 1,nk
       do i = 1, 3
          vkxyz(ik,i,CARTS) = rltv(i,1)*vkxyz(ik,1,BUCS) &
               &            + rltv(i,2)*vkxyz(ik,2,BUCS) &
               &            + rltv(i,3)*vkxyz(ik,3,BUCS)
       end do
    end do

  end subroutine m_Kp_cp_kxyz_to_vkxyz

  subroutine m_Kp_set_kv3(nk)
    integer, intent(in) ::nk
    kv3 = nk
    if(ipri_kp>=1) write(nfout,'(" !kp kv3 = ",i8," <<m_kp_set_kv3>>")') kv3
  end subroutine m_Kp_set_kv3

  subroutine m_Kp_set_mesh_super
    integer :: n
    if( way_ksample == MESH ) then
       n = k_sample_mesh(1,1)
       k_sample_mesh(1,1) = max(nint(dble(n)/n1_sc),1)
!       n = k_sample_mesh(1,1)
       n = k_sample_mesh(2,1)
       k_sample_mesh(2,1) = max(nint(dble(n)/n2_sc),1)
!       n = k_sample_mesh(1,1)
       n = k_sample_mesh(3,1)
       k_sample_mesh(3,1) = max(nint(dble(n)/n3_sc),1)
       k_sample_mesh(:,2) = k_sample_mesh(:,1)

       if(printable) then
          write(nfout,*) 'k-point mesh will be changed.'
          write(nfout,'("mesh=",6i3)') k_sample_mesh
       end if
    else if( way_ksample == MONKHORST_PACK ) then
       n = mp_index(1)
       mp_index(1) = max(nint(dble(n)/n1_sc),1)
       n = mp_index(2)
       mp_index(2) = max(nint(dble(n)/n2_sc),1)
       n = mp_index(3)
       mp_index(3) = max(nint(dble(n)/n3_sc),1)
       if(printable) then
          write(nfout,*) 'k-point mesh will be changed.'
          write(nfout,'("mesh=",3i3)') mp_index
       end if
    end if
  end subroutine m_Kp_set_mesh_super

  subroutine m_Kp_dealloc
    if(allocated(vkxyz)) deallocate(vkxyz)
    if(allocated(qwgt)) deallocate(qwgt)
    if(allocated(k_symmetry)) deallocate(k_symmetry)
    if(allocated(nxyz_tetra)) deallocate(nxyz_tetra)
    if(allocated(ip20)) deallocate(ip20)
    if(allocated(iwt)) deallocate(iwt)
    if(allocated(ip2cub)) deallocate(ip2cub)

! === KT_add ==== 2014/09/30
    if ( allocated( vkxyz_fbz ) ) deallocate( vkxyz_fbz )
    if ( allocated( to_ibz_from_fbz_for_kpoint ) ) then
       deallocate( to_ibz_from_fbz_for_kpoint )
    endif
    if ( allocated( iopr_k_fbz_to_ibz ) ) deallocate( iopr_k_fbz_to_ibz )
    if ( allocated( num_star_of_k ) ) deallocate( num_star_of_k )
    if ( allocated( star_of_k ) ) deallocate( star_of_k )
! ================ 2014/09/30

  end subroutine m_Kp_dealloc

! === KT_add ==== 2014/09/30
!        === modified by T. Y. 2015/01/16 ===
!!$  subroutine gen_kpoint_in_full_BZ( paramset, kv3_fbz, vkxyz_fbz )
  subroutine gen_kpoint_in_full_BZ( paramset)
    logical, intent(in) :: paramset
!!$    integer, intent(in) :: kv3_fbz
!!$    real(kind=DP), intent(inout) :: vkxyz_fbz(kv3_fbz,3,CARTS)
!        ====================================

    integer :: i, ik
    real(kind=DP) :: b1(3), b2(3), b3(3)

    integer :: nface, face(3,nface_max)

    if ( paramset ) return

    b1(1:3) = rltv(1:3,1) /PAI2
    b2(1:3) = rltv(1:3,2) /PAI2
    b3(1:3) = rltv(1:3,3) /PAI2

    do i = 1, 3
       do ik = 1, kv3_fbz
          vkxyz_fbz(ik,i,BUCS) = (altv(1,i) *vkxyz_fbz(ik,1,CARTS) &
               &                + altv(2,i) *vkxyz_fbz(ik,2,CARTS) &
               &                + altv(3,i) *vkxyz_fbz(ik,3,CARTS) ) /PAI2
       enddo
    enddo

    call fbz_face(b1,b2,b3,nface,face)
    call move_kpt_to_inside_fbz( nface, face )

    do i = 1, 3
       do ik = 1, kv3_fbz
          vkxyz_fbz(ik,i,CARTS) = (rltv(i,1) *vkxyz_fbz(ik,1,BUCS) &
               &                + rltv(i,2) *vkxyz_fbz(ik,2,BUCS) &
               &                + rltv(i,3) *vkxyz_fbz(ik,3,BUCS) )
       enddo
    enddo

  contains

    subroutine move_kpt_to_inside_fbz( nface, iface_bz )
      integer, intent(in) :: nface
      integer, intent(in) :: iface_bz(3,nface_max)

      integer :: i, j, ik
      real(kind=DP) :: mtr(3,3)
      real(kind=DP), parameter :: lambda = 1.d0+1.d-8
      real(kind=DP) :: gg, gg_ref, ko(3)
!        === modified by T. Y. 2015/01/16 ===
      integer,parameter :: max_icount = 10000
      integer :: icount
!        ====================================

      do i=1,3
         do j=1,3
            mtr(i,j) = dot_product(rltv(:,i),rltv(:,j))
         end do
      end do

      Do ik=1, kv3_fbz
         ko(1:3) = vkxyz_fbz(ik,1:3,BUCS)
         do i=1,nface
            gg_ref = dot_product(iface_bz(:,i),matmul(mtr,iface_bz(:,i))) * 0.5d0
            icount = 0
            do
               gg = dot_product(iface_bz(:,i),matmul(mtr,ko))
               if(gg > gg_ref*lambda) then
                  ko(1:3) = ko(1:3) - dble(iface_bz(1:3,i))
               else
                  exit
               end if
!        === modified by T. Y. 2015/01/16 ===
               icount = icount+1
               if(icount > max_icount) then
                  write(nfout,'(" icount > max_icount in <<move_kpt_to_inside_fbz>>")')
                  write(nfout,'(" icount = ",i8)') icount
                  call phase_error_with_msg(nfout,' stop at <<move_kpt_to_inside_fbz>>. Icount is too large',__LINE__,__FILE__)
               end if
!        ====================================
            end do
         end do
         vkxyz_fbz(ik,1:3,BUCS) = ko(1:3)
      End do
    end subroutine move_kpt_to_inside_fbz

    subroutine fbz_face(b1,b2,b3,nface,face)
      implicit none

      real(kind=DP), intent(in) :: b1(3),b2(3),b3(3)
      integer, intent(out) :: nface
      integer, intent(out) :: face(3,nface_max)

      ! local variables
      integer :: i,j,k,n
      real(kind=DP) :: gg_ref,gg
      real(kind=DP) :: g1(3),g2(3),g3(3)
      real(kind=DP) :: gface(3,nface_max)
      integer :: iface(3,nface_max)
      logical :: exists_in_fbz(nface_max)
      real(kind=DP), parameter :: lambda = 1.d0+1.d-8

      g1(1:3)=b1(1:3)*0.5d0
      g2(1:3)=b2(1:3)*0.5d0
      g3(1:3)=b3(1:3)*0.5d0

      n=0
      do i=-2,2
        do j=-2,2
          do k=-2,2
            if(i == 0 .and. j == 0 .and. k == 0) cycle
            n=n+1
            gface(1:3,n) = i*g1(1:3) + j*g2(1:3) + k*g3(1:3)
            iface(1,n)=i
            iface(2,n)=j
            iface(3,n)=k
          end do
        end do
      end do

      exists_in_fbz(1:nface_max) = .true.
      do i=1,n
        if(exists_in_fbz(i)) then
          gg_ref = gface(1,i)**2 + gface(2,i)**2 + gface(3,i)**2
          do j=1,n
            if(i /= j) then
              gg = gface(1,i)*gface(1,j) + gface(2,i)*gface(2,j) + gface(3,i)*gface(3,j)
              if(gg*lambda > gg_ref) exists_in_fbz(j) = .false.
            end if
          end do
        end if
      end do

      nface=0
      do i=1,n
        if(exists_in_fbz(i)) then
          nface=nface+1
          face(1:3,nface)=iface(1:3,i)
        end if
      end do

    end subroutine fbz_face

  end subroutine gen_kpoint_in_full_BZ

  subroutine m_Kp_force_kpoint_into_BZ
    integer :: i, ik
    real(kind=DP) :: b1(3), b2(3), b3(3)

    integer :: nface, face(3,nface_max)

    b1(1:3) = rltv(1:3,1) /PAI2
    b2(1:3) = rltv(1:3,2) /PAI2
    b3(1:3) = rltv(1:3,3) /PAI2

    call fbz_face(b1,b2,b3,nface,face)
    call move_kpt_to_inside_fbz( nface, face )

    do i = 1, 3
       do ik = 1, kv3
          vkxyz(ik,i,CARTS) = (rltv(i,1) *vkxyz(ik,1,BUCS) &
               &                + rltv(i,2) *vkxyz(ik,2,BUCS) &
               &                + rltv(i,3) *vkxyz(ik,3,BUCS) )
      enddo
    enddo

  contains

    subroutine move_kpt_to_inside_fbz( nface, iface_bz )
      integer, intent(in) :: nface
      integer, intent(in) :: iface_bz(3,nface_max)

      integer :: i, j, ik
      real(kind=DP) :: mtr(3,3)
      real(kind=DP), parameter :: lambda = 1.d0+1.d-8
      real(kind=DP) :: gg, gg_ref, ko(3)
!        === modified by T. Y. 2015/01/16 ===
      integer,parameter :: max_icount = 10000
      integer :: icount
!        ====================================

      do i=1,3
         do j=1,3
            mtr(i,j) = dot_product(rltv(:,i),rltv(:,j))
         end do
      end do

      Do ik=1, kv3
         ko(1:3) = vkxyz(ik,1:3,BUCS)
         do i=1,nface
            gg_ref = dot_product(iface_bz(:,i),matmul(mtr,iface_bz(:,i))) * 0.5d0
            icount = 0
            do
               gg = dot_product(iface_bz(:,i),matmul(mtr,ko))
               if(gg > gg_ref*lambda) then
                  ko(1:3) = ko(1:3) - dble(iface_bz(1:3,i))
               else
                  exit
               end if
!        === modified by T. Y. 2015/01/16 ===
               icount = icount+1
               if(icount > max_icount) then
                  write(nfout,'(" icount > max_icount in <<move_kpt_to_inside_fbz>>")')
                  write(nfout,'(" icount = ",i8)') icount
                  call phase_error_with_msg(nfout,' stop at <<move_kpt_to_inside_fbz>>. Icount is too large',__LINE__,__FILE__)
               end if
!        ====================================
            end do
         end do
         GvecTrans_kpt(ik,1:3) = nint( ko(1:3) -vkxyz(ik,1:3,BUCS) )
         vkxyz(ik,1:3,BUCS) = ko(1:3)
      ENd Do

    end subroutine move_kpt_to_inside_fbz

    subroutine fbz_face(b1,b2,b3,nface,face)
      implicit none

      real(kind=DP), intent(in) :: b1(3),b2(3),b3(3)
      integer, intent(out) :: nface
      integer, intent(out) :: face(3,nface_max)

      ! local variables
      integer :: i,j,k,n
      real(kind=DP) :: gg_ref,gg
      real(kind=DP) :: g1(3),g2(3),g3(3)
      real(kind=DP) :: gface(3,nface_max)
      integer :: iface(3,nface_max)
      logical :: exists_in_fbz(nface_max)
      real(kind=DP), parameter :: lambda = 1.d0+1.d-8

      g1(1:3)=b1(1:3)*0.5d0
      g2(1:3)=b2(1:3)*0.5d0
      g3(1:3)=b3(1:3)*0.5d0

      n=0
      do i=-2,2
        do j=-2,2
          do k=-2,2
            if(i == 0 .and. j == 0 .and. k == 0) cycle
            n=n+1
            gface(1:3,n) = i*g1(1:3) + j*g2(1:3) + k*g3(1:3)
            iface(1,n)=i
            iface(2,n)=j
            iface(3,n)=k
          end do
        end do
      end do

      exists_in_fbz(1:nface_max) = .true.
      do i=1,n
        if(exists_in_fbz(i)) then
          gg_ref = gface(1,i)**2 + gface(2,i)**2 + gface(3,i)**2
          do j=1,n
            if(i /= j) then
              gg = gface(1,i)*gface(1,j) + gface(2,i)*gface(2,j) + gface(3,i)*gface(3,j)
              if(gg*lambda > gg_ref) exists_in_fbz(j) = .false.
            end if
          end do
        end if
      end do

      nface=0
      do i=1,n
        if(exists_in_fbz(i)) then
          nface=nface+1
          face(1:3,nface)=iface(1:3,i)
        end if
      end do

    end subroutine fbz_face

  end subroutine m_Kp_force_kpoint_into_BZ

  subroutine m_Kp_set_star_of_k
    integer :: iopr, i, j
    integer :: ikbz, jtrs, ik, jk, is
    real(kind=DP) :: x(3), y(3), rkxyz(3), dk(3)
    real(kind=DP), allocatable :: oprec(:,:,:)
    real(kind=DP), allocatable :: wk_star_of_k(:,:)
    integer :: id_sname = -1

    call tstatc0_begin('m_Kp_set_star_of_k ',id_sname)
    allocate(oprec(3,3,nopr))
    do iopr=1,nopr
       do j=1,3
          x = matmul(op(:,:,iopr),rltv(:,j))
          do i=1,3
             y = altv(:,i)
             oprec(i,j,iopr) = dot_product(y,x)/PAI2
          end do
       end do
    end do

    if ( allocated( num_star_of_k ) ) deallocate( num_star_of_k )
    if ( allocated( star_of_k ) )     deallocate( star_of_k )
    if ( allocated( iopr_k_fbz_to_ibz ) )   deallocate( iopr_k_fbz_to_ibz )
    if ( allocated( trev_k_fbz_to_ibz ) )   deallocate( trev_k_fbz_to_ibz )

    allocate( num_star_of_k( kv3 ) );         num_star_of_k = 0;
    allocate( iopr_k_fbz_to_ibz( kv3_fbz ) );  iopr_k_fbz_to_ibz = 0
    allocate( trev_k_fbz_to_ibz( kv3_fbz ) );  trev_k_fbz_to_ibz = 0

! --------------
    Do is=1, nspin
       do ikbz=is, kv3_fbz, nspin
          LOOP_S1: do jtrs=0,1
             do iopr=1,nopr
                rkxyz = (1-2*jtrs) * matmul(oprec(:,:,iopr),vkxyz_fbz(ikbz,1:3,BUCS))
                if ( allocated(magmom_dir_inversion_opr_flag) ) then
                   if ( magmom_dir_inversion_opr_flag(iopr) == -1 ) rkxyz = -rkxyz
                endif

                do jk=is, kv3, nspin
                   dk(1:3) = vkxyz(jk,1:3,BUCS) -rkxyz(1:3)

                   if( abs(dk(1)-nint(dk(1)))<DELTA07  &
                        &      .and. abs(dk(2)-nint(dk(2)))<DELTA07  &
                        &      .and. abs(dk(3)-nint(dk(3)))<DELTA07 ) then

                      iopr_k_fbz_to_ibz(ikbz) = iopr
                      trev_k_fbz_to_ibz(ikbz) = jtrs
                      num_star_of_k(jk) = num_star_of_k(jk) +1
                      exit LOOP_S1
                   end if
                end do
             end do
          end do LOOP_S1
       End do
    End Do
    max_num_star_of_k = maxval( num_star_of_k )

    allocate( star_of_k( kv3, max_num_star_of_k ) );  star_of_k = 0

    num_star_of_k = 0
    Do is=1, nspin
       do ikbz=is, kv3_fbz, nspin
          LOOP_S2: do jtrs=0,1
             do iopr=1,nopr
                rkxyz = (1-2*jtrs) * matmul(oprec(:,:,iopr),vkxyz_fbz(ikbz,1:3,BUCS))
                if ( allocated(magmom_dir_inversion_opr_flag) ) then
                   if ( magmom_dir_inversion_opr_flag(iopr) == -1 ) rkxyz = -rkxyz
                endif

                do jk=is, kv3, nspin
                   dk(1:3) = vkxyz(jk,1:3,BUCS) -rkxyz(1:3)

                   if( abs(dk(1)-nint(dk(1)))<DELTA07  &
                        &      .and. abs(dk(2)-nint(dk(2)))<DELTA07  &
                        &      .and. abs(dk(3)-nint(dk(3)))<DELTA07 ) then

                      iopr_k_fbz_to_ibz(ikbz) = iopr
                      trev_k_fbz_to_ibz(ikbz) = jtrs
                      num_star_of_k(jk) = num_star_of_k(jk) +1
                      star_of_k( jk, num_star_of_k(jk) ) = ikbz
                      exit LOOP_S2
                   end if
                end do
             end do
          end do LOOP_S2
       End do
    End Do
    call tstatc0_end(id_sname)

  end subroutine m_Kp_set_star_of_k

! ==== ASMS === 2018/02/26
  subroutine m_Kp_chk_opr_from_fbz_to_ibz
    integer :: iopr, i, j
    integer :: ikbz, jtrs, ik, jk, is
    real(kind=DP) :: x(3), y(3), rkxyz(3), dk(3), g0(3)
    real(kind=DP), allocatable :: oprec(:,:,:)

    if ( allocated( flg_opr_from_fbz_to_ibz) ) deallocate( flg_opr_from_fbz_to_ibz )
    allocate( flg_opr_from_fbz_to_ibz( nopr ) );  flg_opr_from_fbz_to_ibz = 0

#ifdef USE_ALL_OPR_FOR_CHARGE_SYMMETRIZATION
    nopr_from_fbz_to_ibz = nopr
    flg_opr_from_fbz_to_ibz = 1
    return
#endif

    allocate(oprec(3,3,nopr))
    do iopr=1,nopr
       do j=1,3
          x = matmul(op(:,:,iopr),rltv(:,j))
          do i=1,3
             y = altv(:,i)
             oprec(i,j,iopr) = dot_product(y,x)/PAI2
          end do
       end do
    end do

    if ( way_ksample == SKPS_DIRECT_IN .or. way_ksample == FILE ) then
! --- E operation --
       nopr_from_fbz_to_ibz = 1
       flg_opr_from_fbz_to_ibz(1) = 1
! -----------------

       write(nfout,*) "** nopr_from_fbz_to_ibz is assumed to be 1 "
       return
    endif

    nopr_from_fbz_to_ibz = 0
    Do is=1, nspin
       do ikbz=is, kv3_fbz, nspin
          LOOP_S: do jtrs=0,1
             do iopr=1,nopr
                rkxyz = (1-2*jtrs) * matmul(oprec(:,:,iopr),vkxyz_fbz(ikbz,1:3,BUCS))
                if ( allocated(magmom_dir_inversion_opr_flag) ) then
                   if ( magmom_dir_inversion_opr_flag(iopr) == -1 ) rkxyz = -rkxyz
                endif

                do jk=is, kv3, nspin
                   dk(1:3) = vkxyz(jk,1:3,BUCS) -rkxyz(1:3)
                   if( abs(dk(1)-nint(dk(1)))<DELTA07  &
                        &      .and. abs(dk(2)-nint(dk(2)))<DELTA07  &
                        &      .and. abs(dk(3)-nint(dk(3)))<DELTA07 ) then
                      g0(1:3) = nint(dk(1:3))
                      flg_opr_from_fbz_to_ibz(iopr) = 1
                      exit LOOP_S
                   end if
                end do
             end do
          end do LOOP_S
       End do
    End Do
    Do iopr=1, nopr
       if ( flg_opr_from_fbz_to_ibz(iopr) == 1 ) then
          nopr_from_fbz_to_ibz = nopr_from_fbz_to_ibz +1
       endif
    End do

    if ( ipri_kp >= 2 ) then
       write(nfout,'(A,I5)') &
            &  "=== Number of sym. operations for reducing from fbz to ibz : ", &
            &  nopr_from_fbz_to_ibz
       write(nfout,'(A)') " no.     Flag "
       Do iopr=1, nopr
          write(nfout,'(I5,I5)') iopr, flg_opr_from_fbz_to_ibz(iopr)
       End do
       write(nfout,*)
    endif

    deallocate( oprec )

  end subroutine m_Kp_chk_opr_from_fbz_to_ibz
! ==== ASMS === 2018/02/26

  subroutine m_Kp_wd_BandSymInput(ipri,nfout,nf)
    integer, intent(in) ::ipri, nfout, nf
    character(len=4),dimension(3) :: spinstate
    data spinstate/"    ","  UP","DOWN"/
    integer :: ip, ik
    integer :: id_sname = -1
    call tstatc0_begin('m_Kp_wd_BandSymInput ',id_sname)
    if(mype == 0) then
       write(nf,'("##KPOINTS")')
       ip = 1
       if(allocated(vkxyz_ek) .and. kv3_ek >= kv3) then
          do ik = 1, kv3_ek
             if(nspin==2) ip = 3 - mod(ik,2)
             write(nf,'(i6,3f14.6,2x,a4)') ik,vkxyz_ek(ik,1:3,BUCS),spinstate(ip)
          end do
       else
          do ik = 1, kv3
             if(nspin==2) ip = 3 - mod(ik,2)
             write(nf,'(i6,3f14.6,2x,a4)') ik,vkxyz(ik,1:3,BUCS),spinstate(ip)
          end do
       end if
       write(nf,'("##"/)')
    end if
    call tstatc0_end(id_sname)
  end subroutine m_Kp_wd_BandSymInput

  subroutine m_Kp_get_nkmesh(nkmesh)
    integer, intent(out), dimension(3) :: nkmesh

    if ( way_ksample == MONKHORST_PACK ) then
       nkmesh(1:3) = mp_index(1:3)
    else
       nkmesh(1:3) = k_sample_mesh(1:3,1)
    endif
  end subroutine m_Kp_get_nkmesh

  subroutine m_Kp_get_kptable_bxsf( nkpt_in, kvtab, flg_add_boundary )
    integer, intent(in) :: nkpt_in
    integer, intent(out) :: kvtab(nkpt_in)
    integer, intent(in) :: flg_add_boundary

    integer :: iopr, i, j, nadd_k
    integer :: jtrs, ik, jk, is, num
    integer :: i1, i2, i3, nbz_mesh(3), nkpt_wk
    real(kind=DP) :: x(3), y(3), rkxyz(3), dk(3)
    real(kind=DP), allocatable :: oprec(:,:,:), vkxyz_wk(:,:)
    logical :: found

    allocate(oprec(3,3,nopr))
    do iopr=1,nopr
       do j=1,3
          x = matmul(op(:,:,iopr),rltv(:,j))
          do i=1,3
             y = altv(:,i)
             oprec(i,j,iopr) = dot_product(y,x)/PAI2
          end do
       end do
    end do

    call m_Kp_get_nkmesh( nbz_mesh )
    if ( flg_add_boundary == 1 ) then
       nadd_k = 1
    else
       nadd_k = 0
    endif

    nkpt_wk = ( nbz_mesh(1) +nadd_k ) &
         &   *( nbz_mesh(2) +nadd_k ) &
         &   *( nbz_mesh(3) +nadd_k ) *nspin

    if ( nkpt_wk /= nkpt_in ) call phase_error_with_msg(nfout,"nkpt_wk /= nkpt_in",__LINE__,__FILE__)

    allocate( vkxyz_wk(nkpt_wk,3) )

    num = 0
    Do i1=0, nbz_mesh(1) -1+nadd_k
       Do i2=0, nbz_mesh(2) -1 +nadd_k
          Do i3=0, nbz_mesh(3) -1 +nadd_k
             Do is=1, nspin
                num = num +1
                vkxyz_wk(num,1) = dble(i1)/dble(nbz_mesh(1))
                vkxyz_wk(num,2) = dble(i2)/dble(nbz_mesh(2))
                vkxyz_wk(num,3) = dble(i3)/dble(nbz_mesh(3))
             End do
          End do
       End do
    End Do

    Do is=1, nspin
       do ik=is, nkpt_wk, nspin
          found = .false.
          LOOP_S: do jtrs=0,1
             do iopr=1,nopr
                rkxyz = (1-2*jtrs) * matmul(oprec(:,:,iopr),vkxyz_wk(ik,1:3))
                if ( allocated(magmom_dir_inversion_opr_flag) ) then
                   if ( magmom_dir_inversion_opr_flag(iopr) == -1 ) rkxyz = -rkxyz
                endif

                if ( ekmode == ON ) then
                   do jk=is, kv3_ek, nspin
                      dk(1:3) = vkxyz_ek(jk,1:3,BUCS) -rkxyz(1:3)
                      if( abs(dk(1)-nint(dk(1)))<DELTA07  &
                           &      .and. abs(dk(2)-nint(dk(2)))<DELTA07  &
                           &      .and. abs(dk(3)-nint(dk(3)))<DELTA07 ) then
                         kvtab(ik) = jk;   found = .true.;  exit LOOP_S
                      end if
                   end do
                else
                   do jk=is, kv3, nspin
                      dk(1:3) = vkxyz(jk,1:3,BUCS) -rkxyz(1:3)
                      if( abs(dk(1)-nint(dk(1)))<DELTA07  &
                           &      .and. abs(dk(2)-nint(dk(2)))<DELTA07  &
                           &      .and. abs(dk(3)-nint(dk(3)))<DELTA07 ) then
                         kvtab(ik) = jk;   found = .true.;  exit LOOP_S
                      end if
                   end do
                endif
             End do
          End do LOOP_S
          if ( .not. found ) then
             write(nfout,*) "Not found k point", ik, vkxyz_wk(ik,1:3)
             call phase_error_with_msg(nfout,"Not found k point",__LINE__,__FILE__)
          endif
       End do
    End Do
    deallocate( oprec ); deallocate( vkxyz_wk )

  end subroutine m_Kp_get_kptable_bxsf

end module m_Kpoints
