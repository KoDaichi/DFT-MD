!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 621 $)
!
!  MODULE: m_UnitCell
!
!  AUTHOR(S): J. KOGA @ asms
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
! *************************************************************
!
! =========== Contributions ===================================
!
! Through the courtesy of contributors, the following functions are added.
!
! Company:  ASMS Co.,Ltd.
! Functions:  [Identifier: 13.0B]
!                 cell optimizations with angles fixed is available
!                 ( external_stress is included in the stress term )
!
! =============================================================

module m_UnitCell
! $Id: m_UnitCell.f90 621 2020-05-20 06:51:42Z jkoga $

use m_Const_Parameters, only : DP,QUENCHED_MD,STEEPEST_DESCENT,BFGS,OFF,PAI,CONST_kB, &
                             & PT_CONTROL,P_CONTROL,VELOCITY_SCALING,LANGEVIN, &
                             & VOLUME, METRIC_TENSOR, LATTICE_VECTOR, PAI2, NOSE_HOOVER, T_CONTROL
use m_Stress, only : m_Stress_get_stressmx,m_Stress_get_curr_stress,m_Stress_set_stressmx
use m_IterationNumbers, only : iteration_unit_cell,iteration_ionic
use m_Control_Parameters, only : nhistory_stress,printable,delta_stress,lattice_optimization_method &
                                ,max_stress,ipri,eigenvalue_threshold,sw_uniform,ipriunitcell,external_stress &
                                ,dtio,imdalg, max_mdstep
use m_Files, only : nfout,nfdynm, nflatconst, nfmetric, m_Files_flush_pcontrol_files
use m_Ionic_System, only : altv,pos,natm,cps,m_baro,target_pressure,cpd_l,qmass,tkb,ityp,amion &
                         & ,t_ctrl_method,scale_velocity, mobility_cell, control_pressure &
                         & ,sw_shift_velocities, shift_velocities, sw_average_diagonal, nrsv, imdtyp &
                         & ,p_ctrl_method
use m_Crystal_Structure, only : m_CS_altv_2_rltv,p2bmat,a,b,c,ca,cb,cc,il,kt_iuctype,rltv,univol,rvol, b2pmat
use m_Parallelization, only : mype,npes,MPI_CommGroup

! ====================================== KT_add ============= 13.0B &13.1AS
use m_Control_Parameters,  only : sw_fix_lattice_angles, fix_angle_alpha, &
     &                            fix_angle_beta, fix_angle_gamma, &
     &                            sw_fix_lattice_lengths, fix_length_a, &
     &                            fix_length_b, fix_length_c, &
     &                            sw_neglect_stress_offdiagonal, &
     &                            sw_fix_lattice_shapes, fix_shape_ab_plane, &
     &                            fix_shape_bc_plane, fix_shape_ac_plane, &
     &                            sw_optimize_coords_sametime, lattice_coords_opt_mode, &
     &                            icond, omega_for_cellopt, stress_force_mul_factor
  use m_Const_Parameters,   only : ON, INITIAL
! =========================================================== 13.0B &13.1AS

  use m_Total_Energy, only : etotal
  use m_Ionic_System, only : natm, imdtypxyz
  use m_Force,        only : forc_l, m_Forces_are_Converged_core
  use mpi

  implicit none
!  include 'mpif.h'

  real(kind=DP), private, dimension(3,3) :: stress_forc
  real(kind=DP), private, dimension(3,3) :: stress_forc_old
  real(kind=DP), private, dimension(3,3) :: stress_velocity
  real(kind=DP), private, allocatable, dimension(:,:,:) :: cellvec_history
  real(kind=DP), private, allocatable, dimension(:,:,:) :: stress_forc_history
  real(kind=DP), private, allocatable, dimension(:) :: volume_history
  real(kind=DP), private, allocatable, dimension(:) :: volume_forc_history

  logical,private :: first_call = .true.
  logical,private :: H0_has_been_set = .false.
  character(len('unit_cell_optimization')),private,parameter :: tag_unit_cell_optimization='unit_cell_optimization'
  character(len("lattice_vector")),private,parameter :: tag_lattice_vector = "lattice_vector"
  character(len("pressure_control")), private, parameter :: tag_pressure_control='pressure_control'
  real(kind=DP) :: maxoptforc_old=1000.d0
  real(kind=DP) :: univolo

! for NPT and NPH simulation based on the metric tensor
  real(kind=DP), allocatable, dimension(:,:) :: p_l
  real(kind=DP), allocatable, dimension(:,:) :: p_l_lattice
!real(kind=DP) :: s_thermo,p_thermo,m_thermo
  real(kind=DP), dimension(3,3) :: altv_inv,metric,metric_inv,dudmetric
!real(kind=DP), dimension(3,3) :: p_baro, p_baro_tmp
  real(kind=DP), dimension(3,3) :: factor
  real(kind=DP), allocatable, dimension(:,:) :: forc_l_lattice
  real(kind=DP), allocatable, dimension(:,:) :: p_l_metric

  integer :: nmax = 1000
  real(kind=DP) :: eps = 1.d-8
  real(kind=DP) :: H0=0.d0

  integer :: nrsv_t

  type thermostats_t
     integer :: natm
     integer, pointer, dimension(:) :: target_atoms
     real(kind=DP) :: m_thermo, target_temperature
     real(kind=DP) :: s_thermo, p_thermo
     real(kind=DP), dimension(3,3) :: p_baro,p_baro_tmp
     real(kind=DP) :: gfree
  end type thermostats_t

  type(thermostats_t), allocatable, dimension(:) :: thermostats

  real(kind=DP) :: vol=0.d0, vol_old,econst
  real(kind=DP) :: pv_baro,pv_baro_old
  real(kind=DP),dimension(3,3) :: latv

! ------------------
! simultaneous optimization of lattice and coords
!
  real(kind=DP), private, allocatable :: cps_history(:,:,:)
  real(kind=DP), private, allocatable :: atom_forc_history(:,:,:)

  real(kind=DP),dimension(3,3) :: altv0
  logical :: altv0_is_set = .false.
! ------------------
contains

function get_volume(avec,bvec,cvec) result(vol)
  real(kind=DP), dimension(3), intent(in) :: avec,bvec,cvec
  real(kind=DP), dimension(3) :: tmpvec
  real(kind=DP) :: vol
  tmpvec(1) = bvec(2)*cvec(3)-bvec(3)*cvec(2)
  tmpvec(2) = bvec(3)*cvec(1)-bvec(1)*cvec(3)
  tmpvec(3) = bvec(1)*cvec(2)-bvec(2)*cvec(1)
  vol = dot_product(avec,tmpvec)
end function get_volume

logical function m_UC_converged()
  integer :: i,j
  real(kind=DP) :: maxs,astress
  real(kind=DP), dimension(3,3) :: stress_tensor, stress_work
  real(kind=DP), dimension(3,3) :: altv_brav

  if(imdalg == PT_CONTROL .or. imdalg == P_CONTROL)then
    !m_UC_converged = .false.
    m_UC_converged = iteration_unit_cell.ge.max_mdstep
    return;
  endif
  stress_tensor = m_Stress_get_curr_stress()
  maxs = 0.d0

  if(sw_uniform==OFF)then
    do i=1,3
      do j=1,3
! ============ KT mod ============================= 13.0B
!         astress = abs(stress_tensor(i,j)-external_stress(i,j))
         astress = abs(stress_tensor(i,j))
! ================================================= 13.0B
         if(astress>maxs) maxs = astress
      enddo
    enddo
  else
! ============ KT mod ============================= 13.0B
!       maxs = maxs+abs(stress_tensor(i,i)-external_stress(i,i))
!       maxs = maxs+abs(stress_tensor(i,i))
    maxs = abs(stress_tensor(1,1)+stress_tensor(2,2)+stress_tensor(3,3))/3.0d0
! ================================================= 13.0B
  endif

! === KT_add === 13.1AS
  if ( sw_neglect_stress_offdiagonal == ON ) then
     stress_work = 0.0d0;
     stress_work(1,1) = stress_tensor(1,1)
     stress_work(2,2) = stress_tensor(2,2)
     stress_work(3,3) = stress_tensor(3,3)

     maxs = 0.0d0
     do i=1,3
        astress = sqrt( stress_work(i,1)**2 +stress_work(i,2)**2 +stress_work(i,3)**2 )
        if(astress>maxs) maxs = astress
     enddo
  endif

  if ( sw_fix_lattice_angles==ON .or. sw_fix_lattice_lengths == ON ) then
     stress_work = stress_tensor

     if ( sw_fix_lattice_shapes == ON ) then
        if ( fix_shape_ab_plane ) then
           stress_work(1,1) = ( stress_tensor(1,1)+stress_tensor(2,2) )/2.0d0
           stress_work(2,2) = ( stress_tensor(1,1)+stress_tensor(2,2) )/2.0d0
           stress_work(1,2) = 0.0d0
           stress_work(2,1) = 0.0d0
        endif
        if ( fix_shape_ac_plane ) then
           stress_work(1,1) = ( stress_tensor(1,1)+stress_tensor(3,3) )/2.0d0
           stress_work(3,3) = ( stress_tensor(1,1)+stress_tensor(3,3) )/2.0d0
           stress_work(1,3) = 0.0d0
           stress_work(3,1) = 0.0d0
        endif
        if ( fix_shape_bc_plane ) then
           stress_work(2,2) = ( stress_tensor(2,2)+stress_tensor(3,3) )/2.0d0
           stress_work(3,3) = ( stress_tensor(2,2)+stress_tensor(3,3) )/2.0d0
           stress_work(2,3) = 0.0d0
           stress_work(3,2) = 0.0d0
        endif
     endif

     call get_latvec_from_prim_to_brav( p2bmat, altv, altv_brav )

     stress_work = projected_stress_along_latvecs( altv_brav, stress_work )
     maxs = 0.0d0

     if ( sw_fix_lattice_angles == ON ) then
        stress_work = constrain_on_force_type1( altv_brav, stress_work )
     endif
     if ( sw_fix_lattice_lengths==ON ) then
        stress_work = constrain_on_force_type2( altv_brav, stress_work )
     endif

     do i=1,3
        astress = sqrt( stress_work(i,1)**2 +stress_work(i,2)**2 +stress_work(i,3)**2 )
        if(astress>maxs) maxs = astress
     enddo
  endif
! ========= 13.1AS
  m_UC_converged = max_stress>maxs
  if ( sw_optimize_coords_sametime == ON ) then
     m_UC_converged = m_UC_converged .and. m_Forces_are_Converged_core()
  endif

  if(printable) then
    write(nfout,'(a)')         '!**'
    write(nfout,'(a)')         '!** checking convergence of the stress tensor ...'
    write(nfout,'(a,2f20.10)') '!** max(stress-external_stress), tolerance : ',maxs,max_stress
    write(nfout,'(a,l2)')      '!** convergence : ',m_UC_converged
    write(nfout,'(a)')         '!**'
  endif
end function m_UC_converged

subroutine m_UC_init()
   stress_forc = 0.0d0
   stress_forc_old = 0.0d0
   stress_velocity = 0.0d0
   allocate(cellvec_history(3,3,nhistory_stress));cellvec_history=0.d0
   allocate(stress_forc_history(3,3,nhistory_stress));stress_forc_history=0.0d0
   allocate(volume_history(nhistory_stress));volume_history=0.d0
   allocate(volume_forc_history(nhistory_stress));volume_forc_history=0.d0

   if ( sw_optimize_coords_sametime == ON ) then
      allocate( cps_history( natm, 3, nhistory_stress ) )
      allocate( atom_forc_history( natm, 3, nhistory_stress ) )
      cps_history = 0.d0;  atom_forc_history = 0.0d0
   endif
end subroutine m_UC_init

subroutine m_UC_wd_cntn_data(nfcntn)
   integer,intent(in) :: nfcntn
   integer :: i, j, nbox
   if(iteration_unit_cell-1<=nhistory_stress)then
      nbox = iteration_unit_cell-1
   else
      nbox = nhistory_stress
   endif
   if(first_call.or.nbox==0) return
   if(mype==0)then
      write(nfcntn,*) tag_unit_cell_optimization
      write(nfcntn,*) 'history of the unit cell'
      do i=1,nbox
         if(sw_uniform==OFF)then
            write(nfcntn,'(3f25.15)') cellvec_history(1,1,i),cellvec_history(1,2,i),cellvec_history(1,3,i)
            write(nfcntn,'(3f25.15)') cellvec_history(2,1,i),cellvec_history(2,2,i),cellvec_history(2,3,i)
            write(nfcntn,'(3f25.15)') cellvec_history(3,1,i),cellvec_history(3,2,i),cellvec_history(3,3,i)
         else
            write(nfcntn,'(f25.15)') volume_history(i)
         endif
      enddo
      write(nfcntn,*) 'history of the force acting on the unit cell'
      do i=1,nbox
         if(sw_uniform==OFF)then
            write(nfcntn,'(3f25.15)') stress_forc_history(1,1,i),stress_forc_history(1,2,i),stress_forc_history(1,3,i)
            write(nfcntn,'(3f25.15)') stress_forc_history(2,1,i),stress_forc_history(2,2,i),stress_forc_history(2,3,i)
            write(nfcntn,'(3f25.15)') stress_forc_history(3,1,i),stress_forc_history(3,2,i),stress_forc_history(3,3,i)
         else
            write(nfcntn,'(f25.15)') volume_forc_history(i)
         endif
      enddo
      write(nfcntn,*) 'max. stress'
      write(nfcntn,'(f25.15)') m_Stress_get_stressmx()
      write(nfcntn,*) 'maxoptforc'
      write(nfcntn,'(f25.15)') maxoptforc_old

      if ( sw_optimize_coords_sametime == ON ) then
         write(nfcntn,*) 'initial lattice coords'
         write(nfcntn,'(3f25.15)') altv0(1,1:3)
         write(nfcntn,'(3f25.15)') altv0(2,1:3)
         write(nfcntn,'(3f25.15)') altv0(3,1:3)
         write(nfcntn,*) 'history of the atomic coords'
         do i=1,nbox
            Do j=1, natm
               write(nfcntn,'(3F25.15)') cps_history(j,1:3,i)
            End Do
         end do
         write(nfcntn,*) 'history of the atomic forces'
         do i=1,nbox
            Do j=1, natm
               write(nfcntn,'(3F25.15)') atom_forc_history(j,1:3,i)
            End Do
         end do
      endif
   endif
end subroutine m_UC_wd_cntn_data

subroutine m_UC_rd_cntn_data(nfcntn)
   integer, intent(in) :: nfcntn
   integer :: i, j, nbox, ierr
   logical             :: tag_is_found, EOF_reach
   integer,      parameter   :: len_str = 132
   character(len=len_str)       :: str
   real(kind=DP) :: strmx
   if(first_call) then
      call m_UC_init()
      first_call = .false.
   endif
   first_call = .false.
   if(mype==0)then
      if(iteration_unit_cell-1<=nhistory_stress)then
         nbox = iteration_unit_cell-1
      else
         nbox = nhistory_stress
      endif
      call rewind_to_tag0(nfcntn,len(tag_unit_cell_optimization), &
            &  tag_unit_cell_optimization &
            &, EOF_reach, tag_is_found, str,len_str)
      if(tag_is_found)then
         read(nfcntn,*)
         do i=1,nbox
            if(sw_uniform==OFF)then
               read(nfcntn,*) cellvec_history(1,1,i),cellvec_history(1,2,i),cellvec_history(1,3,i)
               read(nfcntn,*) cellvec_history(2,1,i),cellvec_history(2,2,i),cellvec_history(2,3,i)
               read(nfcntn,*) cellvec_history(3,1,i),cellvec_history(3,2,i),cellvec_history(3,3,i)
            else
               read(nfcntn,*) volume_history(i)
            endif
         enddo
         read(nfcntn,*)
         do i=1,nbox
            if(sw_uniform==OFF)then
               read(nfcntn,*) stress_forc_history(1,1,i),stress_forc_history(1,2,i),stress_forc_history(1,3,i)
               read(nfcntn,*) stress_forc_history(2,1,i),stress_forc_history(2,2,i),stress_forc_history(2,3,i)
               read(nfcntn,*) stress_forc_history(3,1,i),stress_forc_history(3,2,i),stress_forc_history(3,3,i)
            else
               read(nfcntn,*) volume_forc_history(i)
            endif
         enddo
         read(nfcntn,*)
         read(nfcntn,*) strmx
         read(nfcntn,*)
         read(nfcntn,*) maxoptforc_old

         if ( sw_optimize_coords_sametime == ON ) then
            read(nfcntn,*)
            read(nfcntn,*) altv0(1,1:3)
            read(nfcntn,*) altv0(2,1:3)
            read(nfcntn,*) altv0(3,1:3)
            read(nfcntn,*)
            altv0_is_set = .true.

            do i=1,nbox
               Do j=1, natm
                  read(nfcntn,*) cps_history(j,1:3,i)
               End Do
            end do
            read(nfcntn,*)
            do i=1,nbox
               Do j=1, natm
                  read(nfcntn,*) atom_forc_history(j,1:3,i)
               End Do
            end do
         endif
      endif
      call rewind_to_tag0(nfcntn,len(tag_lattice_vector),tag_lattice_vector &
           &,  EOF_reach, tag_is_found, str, len_str)
      if(tag_is_found)then
         read(nfcntn,*) altv(1,1),altv(2,1),altv(3,1)
         read(nfcntn,*) altv(1,2),altv(2,2),altv(3,2)
         read(nfcntn,*) altv(1,3),altv(2,3),altv(3,3)
      endif
   endif
   if(npes>1)then
      call mpi_bcast(tag_is_found,1,mpi_logical,0,MPI_CommGroup,ierr)
      if(tag_is_found)then
         if(sw_uniform==OFF)then
            call mpi_bcast(cellvec_history,9*nhistory_stress,mpi_double_precision,0,MPI_CommGroup,ierr)
            call mpi_bcast(stress_forc_history,9*nhistory_stress,mpi_double_precision,0,MPI_CommGroup,ierr)
         else
            call mpi_bcast(volume_history,nhistory_stress,mpi_double_precision,0,MPI_CommGroup,ierr)
            call mpi_bcast(volume_forc_history,nhistory_stress,mpi_double_precision,0,MPI_CommGroup,ierr)
         endif
         call mpi_bcast(strmx,1,mpi_double_precision,0,MPI_CommGroup,ierr)
         call mpi_bcast(maxoptforc_old,1,mpi_double_precision,0,MPI_CommGroup,ierr)
         call mpi_bcast(altv,9,mpi_double_precision,0,MPI_CommGroup,ierr)
         if ( sw_optimize_coords_sametime == ON ) then
            call mpi_bcast( cps_history, natm *3 *nhistory_stress, &
                 &          mpi_double_precision, 0, MPI_CommGroup, ierr )
            call mpi_bcast( atom_forc_history, natm *3 *nhistory_stress, &
                 &          mpi_double_precision, 0, MPI_CommGroup, ierr )
            call mpi_bcast( altv0, 9, mpi_double_precision, 0, MPI_CommGroup, ierr )
            call mpi_bcast( altv0_is_set, 1, mpi_logical, 0, MPI_CommGroup, ierr )
         endif
      endif
   endif
   call m_CS_altv_2_rltv()
   call m_Stress_set_stressmx(strmx)
end subroutine m_UC_rd_cntn_data

function get_averaged_cell_forc(stress_tensor) result(res)
   real(kind=DP), dimension(:,:), intent(in) :: stress_tensor
   real(kind=DP) :: stav,exstav
   real(kind=DP) :: res
   stav = stress_tensor(1,1)+stress_tensor(2,2)+stress_tensor(3,3)
! ================ KT_mod ================= 13.0B
!   exstav = external_stress(1,1)+external_stress(2,2)+external_stress(3,3)
!   res = (stav-exstav)/3.0d0
!
   res = stav /3.0d0
! ========================================= 13.0B
end function get_averaged_cell_forc

subroutine m_UC_doit()
   integer :: i
   real(kind=DP), dimension(3,3) :: stress_tensor,stress_tensor_p
   real(kind=DP), dimension(3,3) :: stress_tensor_tmp
   real(kind=DP), dimension(3,3) :: work
   real(kind=DP), dimension(3,3) :: altv_brav
   real(kind=DP), dimension(3,3) :: rltv_t
   real(kind=DP), dimension(3) :: avec,bvec,cvec
   real(kind=DP) :: old_length(3), new_length(3)
   real(kind=DP) :: vol,dv,factor,mstress
   logical :: update_lattice_only

   if ( first_call ) then
      call m_UC_init();    first_call = .false.
   endif

   if(printable) write(nfout,'(a)')      ' -- lattice optimization --'
   stress_tensor = m_Stress_get_curr_stress()
   stress_tensor_p = 0.0d0

   stress_forc = get_force( altv, stress_tensor )

! ================================= KT_add ================== 13.0B &13.1AS
   if ( sw_neglect_stress_offdiagonal == ON ) then
      stress_tensor_tmp = 0.0d0
      stress_tensor_tmp(1,1) = stress_tensor(1,1)
      stress_tensor_tmp(2,2) = stress_tensor(2,2)
      stress_tensor_tmp(3,3) = stress_tensor(3,3)

      stress_forc = get_force( altv, stress_tensor_tmp )
   endif
   if ( sw_fix_lattice_angles==ON .or. sw_fix_lattice_lengths == ON ) then
      stress_tensor_tmp = stress_tensor

      if ( sw_fix_lattice_shapes == ON ) then
         if ( fix_shape_ab_plane ) then
            stress_tensor_tmp(1,1) = ( stress_tensor(1,1)+stress_tensor(2,2) )/2.0d0
            stress_tensor_tmp(2,2) = ( stress_tensor(1,1)+stress_tensor(2,2) )/2.0d0
            stress_tensor_tmp(1,2) = 0.0d0
            stress_tensor_tmp(2,1) = 0.0d0
         endif
         if ( fix_shape_ac_plane ) then
            stress_tensor_tmp(1,1) = ( stress_tensor(1,1)+stress_tensor(3,3) )/2.0d0
            stress_tensor_tmp(3,3) = ( stress_tensor(1,1)+stress_tensor(3,3) )/2.0d0
            stress_tensor_tmp(1,3) = 0.0d0
            stress_tensor_tmp(3,1) = 0.0d0
         endif
         if ( fix_shape_bc_plane ) then
            stress_tensor_tmp(2,2) = ( stress_tensor(2,2)+stress_tensor(3,3) )/2.0d0
            stress_tensor_tmp(3,3) = ( stress_tensor(2,2)+stress_tensor(3,3) )/2.0d0
            stress_tensor_tmp(2,3) = 0.0d0
            stress_tensor_tmp(3,2) = 0.0d0
         endif
      endif

      call get_latvec_from_prim_to_brav( p2bmat, altv, altv_brav )
      stress_forc = get_force( altv_brav, stress_tensor_tmp )

     if ( sw_fix_lattice_angles==ON ) then
        stress_forc = constrain_on_force_type1( altv_brav, stress_forc )
     endif
     if ( sw_fix_lattice_lengths==ON ) then
        stress_forc = constrain_on_force_type2( altv_brav, stress_forc )
     endif
  endif
! =========================================================== 13.0B &13.1AS

   mstress = m_Stress_get_stressmx()
   if(printable)then
      write(nfout,'(a)')
      write(nfout,'(a)')         ' -- current lattice --'
      write(nfout,'(a,3f20.10)') '    a_vector ',altv(1,1),altv(2,1),altv(3,1)
      write(nfout,'(a,3f20.10)') '    b_vector ',altv(1,2),altv(2,2),altv(3,2)
      write(nfout,'(a,3f20.10)') '    c_vector ',altv(1,3),altv(2,3),altv(3,3)
      write(nfout,'(a)')

      write(nfout,'(a)')         ' -- current volume --'
      write(nfout,'(a,f22.12)')  '    univol ', univol
      write(nfout,'(a)')

      write(nfout,'(a,i8,a)')    ' -- stress tensor obtained from iteration_unit_cell ',iteration_unit_cell-1,' --'
      do i=1,3
         write(nfout,'(3f20.10)') stress_tensor(i,1),stress_tensor(i,2),stress_tensor(i,3)
      enddo
      write(nfout,'(a)')         ' -- current cps and pos --'
      do i=1,natm
         write(nfout,'(6f20.10)') cps(i,1),cps(i,2),cps(i,3),pos(i,1),pos(i,2),pos(i,3)
      enddo
      write(nfout,'(a,f20.10,a)') ' -- max. stress : ',mstress,' --'
      write(nfout,'(a)')
      write(nfout,'(a)')          ' -- force acting on the unit cell --'
      write(nfout,'(a,3f20.10)')  '    a_vector ',stress_forc(1,1),stress_forc(1,2),stress_forc(1,3)
      write(nfout,'(a,3f20.10)')  '    b_vector ',stress_forc(2,1),stress_forc(2,2),stress_forc(2,3)
      write(nfout,'(a,3f20.10)')  '    c_vector ',stress_forc(3,1),stress_forc(3,2),stress_forc(3,3)
      if(lattice_optimization_method==QUENCHED_MD)then
         write(nfout,'(a)')          ' -- velocity of the unit cell --'
         write(nfout,'(a,3f20.10)')  '    a_vector ',stress_velocity(1,1),stress_velocity(1,2),stress_velocity(1,3)
         write(nfout,'(a,3f20.10)')  '    b_vector ',stress_velocity(2,1),stress_velocity(2,2),stress_velocity(2,3)
         write(nfout,'(a,3f20.10)')  '    c_vector ',stress_velocity(3,1),stress_velocity(3,2),stress_velocity(3,3)
      endif
   endif

   if ( sw_fix_lattice_angles==ON .or. sw_fix_lattice_lengths == ON ) then
      altv = altv_brav
      old_length(1) = sqrt( altv(1,1)**2 +altv(2,1)**2 +altv(3,1)**2 )  ! a0
      old_length(2) = sqrt( altv(1,2)**2 +altv(2,2)**2 +altv(3,2)**2 )  ! b0
      old_length(3) = sqrt( altv(1,3)**2 +altv(2,3)**2 +altv(3,3)**2 )  ! c0
   endif

   if ( sw_optimize_coords_sametime == ON ) then
      if ( lattice_optimization_method == QUENCHED_MD &
           &  .or. lattice_optimization_method == STEEPEST_DESCENT ) then
         call phase_error_with_msg( nfout, &
              &  ' sw_optimize_coords_sametime works with BFGS' &
              &   ,__LINE__,__FILE__)
      endif
      if ( sw_uniform == ON ) then
         call phase_error_with_msg( nfout, &
              &  ' sw_optimize_coords_sametime works with sw_uniform=OFF' &
              &   ,__LINE__,__FILE__)
      endif
   endif

   update_lattice_only = .true.
   if ( sw_optimize_coords_sametime == ON ) then
      update_lattice_only = .false.
!      if ( lattice_coords_opt_mode == 0 .and. m_Forces_are_Converged_core() ) then
!         update_lattice_only = .true.
!      endif
   else
      update_lattice_only = .true.
   endif

   if (sw_uniform==OFF)then
      call store_cell_and_stress()
      if(lattice_optimization_method==QUENCHED_MD .or. &
      &  lattice_optimization_method==STEEPEST_DESCENT) then
         call update_lattice()
         call update_velocities()
      else if (lattice_optimization_method==BFGS) then
         if ( sw_optimize_coords_sametime == ON ) then
            if ( update_lattice_only ) then
               call update_lattice_by_bfgs()
            else
               call update_lat_and_coord_by_bfgs( lattice_coords_opt_mode )
            endif
         else
            call update_lattice_by_bfgs()
         endif
      endif
   else
      avec=altv(:,1);      bvec=altv(:,2);      cvec=altv(:,3)
      vol = get_volume(avec,bvec,cvec)
      dv = get_averaged_cell_forc(stress_tensor)
      call store_volume(vol,dv)
      call update_volume_by_bfgs()
   endif

   if ( sw_fix_lattice_angles==ON .or. sw_fix_lattice_lengths == ON ) then
      new_length(1) = sqrt( altv(1,1)**2 +altv(2,1)**2 +altv(3,1)**2 ) ! a0
      new_length(2) = sqrt( altv(1,2)**2 +altv(2,2)**2 +altv(3,2)**2 ) ! b0
      new_length(3) = sqrt( altv(1,3)**2 +altv(2,3)**2 +altv(3,3)**2 ) ! c0
      if ( fix_length_a ) altv(:,1) = altv(:,1) *old_length(1) /new_length(1)
      if ( fix_length_b ) altv(:,2) = altv(:,2) *old_length(2) /new_length(2)
      if ( fix_length_c ) altv(:,3) = altv(:,3) *old_length(3) /new_length(3)

      altv_brav = altv
      call get_latvec_from_brav_to_prim( b2pmat, altv_brav, altv )
   endif

   stress_forc_old = stress_forc
   univolo = univol
   call altv_2_rltv(altv,rltv,univol,rvol)  ! in b_CS

   if ( sw_optimize_coords_sametime == ON .and. (.not.update_lattice_only) ) then
      rltv_t = transpose(rltv)/PAI2
      call change_of_coordinate_system(rltv_t,cps,natm,natm,pos)
   else
      call change_of_coordinate_system(altv,pos,natm,natm,cps) !-(b_I.S.) pos -> cps
   endif
   call primitive2bravais(nfout,p2bmat,altv(:,1),altv(:,2),altv(:,3),a,b,c,ca,cb,cc,il)
                                                                ! in b_CS
   if(printable)then
      write(nfout,'(a)')         ' -- new lattice --'
      write(nfout,'(a,3f20.10)') '    a_vector ',altv(1,1),altv(2,1),altv(3,1)
      write(nfout,'(a,3f20.10)') '    b_vector ',altv(1,2),altv(2,2),altv(3,2)
      write(nfout,'(a,3f20.10)') '    c_vector ',altv(1,3),altv(2,3),altv(3,3)
      write(nfout,'(a)')
      write(nfout,'(a)')         ' -- new volume --'
      write(nfout,'(a,f22.12)')  '    univol ', univol
      write(nfout,'(a)')
      write(nfout,'(a)')         ' -- new cps and pos --'
      do i=1,natm
         write(nfout,'(6f20.10)') cps(i,1),cps(i,2),cps(i,3),pos(i,1),pos(i,2),pos(i,3)
      enddo
!
      if ( kt_iuctype == 0 .and. ipriunitcell >=2 ) then   ! bravias
         write(nfout,'(a)')
         write(nfout,'(a)') ' --- new lattice (bravais) ---- '
         write(nfout,'(3(a,f18.12))')  '  a     = ', a, '  b    = ', b, '  c     = ', c
         write(nfout,'(3(a,f18.12))')  '  alpha = ', acos(ca)/PAI*180.0d0, &
              &                           '  beta = ', acos(cb)/PAI*180.0d0, &
              &                           '  gamma = ', acos(cc)/PAI*180.0d0
         write(nfout,'(a)')

         work(1,:) = altv(:,1);    work(2,:) = altv(:,2);    work(3,:) = altv(:,3)
         work = matmul( p2bmat, work )
         write(nfout,'(a,3f20.10)') '    avec = ',work(1,1:3)
         write(nfout,'(a,3f20.10)') '    bvec = ',work(2,1:3)
         write(nfout,'(a,3f20.10)') '    cvec = ',work(3,1:3)
         write(nfout,'(a)') ' ----------------------------- '
      endif

   endif

end subroutine m_UC_doit

function m_UC_get_univol_old() result(res)
   real(kind=DP) :: res
   res = univolo
end function m_UC_get_univol_old

subroutine get_bfgs_force(nfree,nbox,coord,forc,forc_bfgs,forc_bfgs_max)
   integer, intent(in) :: nfree
   real(kind=DP),intent(in),dimension(nfree,nbox) :: coord,forc
   real(kind=DP),intent(out),dimension(nfree) :: forc_bfgs
   real(kind=DP),intent(out) :: forc_bfgs_max
   integer :: nbox,i,j,k,info
   real(kind=DP),dimension(nfree) :: xdelta,gdelta,gdotinvh
   real(kind=DP),dimension(nfree,nfree) :: ihess
   real(kind=DP),dimension(nfree*(nfree+1)/2) :: amat
   real(kind=DP),dimension(nfree) :: eigv
   real(kind=DP),dimension(nfree,nfree) :: eigvec
   real(kind=DP),dimension(nfree*3) :: workar
   real(kind=DP) :: xgi,gihg
   logical :: corrected_eig

   ihess=0.0d0
   do i=1,nfree
      ihess(i,i)=1.0d0
   enddo
   do i=2,nbox
      do j=1,nfree
         xdelta(j) = coord(j,i) - coord(j,i-1)
         gdelta(j) = -forc(j,i) + forc(j,i-1)
      enddo
      xgi = 1.0d0/dot_product(xdelta,gdelta)
      if(xgi<0)then
         if(printable)then
            write(nfout,'(a,i3)') '!** WARNING dx dot dg is negative for history ',i
            write(nfout,'(a)') 'skipping this update.'
         endif
         cycle
      endif
      do j=1,nfree
         gdotinvh(j)=dot_product(ihess(j,:),gdelta(:))
      enddo
      gihg = dot_product(gdelta,gdotinvh)
      do j=1,nfree
         do k=1,nfree
            ihess(j,k) = ihess(j,k)+xgi*xgi*(1.0d0/xgi+gihg)*xdelta(j)*xdelta(k) &
  &                    - (gdotinvh(j)*xdelta(k)+gdotinvh(k)*xdelta(j))*xgi
         enddo
      enddo
   enddo

!  correct bad eigenvalues present in the Hessian
!  we will always do this, since the matrix size is fixed & small
   corrected_eig=.false.
   amat=0.0d0;eigv=0.0d0;eigvec=0.0d0;workar=0.0d0
   do i=1,nfree
      do j=i,nfree
         amat(i + (j-1)*j/2) = ihess(i,j)
      enddo
   enddo
   call dspev('V','U',nfree,amat,eigv,eigvec,nfree,workar,info)
   if(printable .and. ipriunitcell>=2) write(nfout,'(a)') '--- eigenvalues for the approximate Hessian ---'
   do i=1,nfree
      if(printable.and.ipriunitcell>=2) write(nfout,'(i8,f20.10)') i,1.0d0/eigv(i)
      if (1.0d0/eigv(i)<eigenvalue_threshold)then
         eigv(i) = 1.0d0/eigenvalue_threshold
         if(printable.and.ipriunitcell>=2) write(nfout,'(a,i8,a,f20.10)') &
         &  'corrected the eigenvalue for the ',i,'-th element to : ',1.0d0/eigv(i)
         corrected_eig=.true.
      endif
   enddo
   if(corrected_eig)then
      ihess=0.d0
      do i=1,nfree
         do j=1,nfree
            do k=1,nfree
               ihess(i,j) = ihess(i,j)+eigvec(i,k)*eigvec(j,k)*eigv(k)
            enddo
         enddo
      enddo
   endif

   forc_bfgs=0.0d0
   forc_bfgs_max=0.0d0
   do i=1,nfree
      do j=1,nfree
         forc_bfgs(i) = forc_bfgs(i)+ihess(i,j)*forc(j,nbox)
      enddo
      if(abs(forc_bfgs(i))>forc_bfgs_max) forc_bfgs_max = abs(forc_bfgs(i))
   enddo

end subroutine get_bfgs_force

subroutine update_volume_by_bfgs()
   integer :: nbox
   real(kind=DP),dimension(1) :: forc_bfgs
   real(kind=DP) :: forc_bfgs_max,factor
   if(iteration_unit_cell-1<=nhistory_stress)then
      nbox = iteration_unit_cell-1
   else
      nbox = nhistory_stress
   endif
   call get_bfgs_force(1,nbox,volume_history,volume_forc_history,forc_bfgs,forc_bfgs_max)
   if(printable) write(nfout,*) 'max. bfgs force : ',forc_bfgs_max
   if(maxoptforc_old*10<forc_bfgs_max)then
      if (printable) write(nfout,'(a)') &
      & 'the estimated force acting on the cell seems to be very large; &
      & update will be done by the steepest-descent method'
      factor = (1+volume_forc_history(nbox)*delta_stress)**(1.d0/3.d0)
   else
      factor = (1+forc_bfgs(1))**(1.d0/3.d0)
   endif
   altv(:,:) = altv(:,:)*factor
   maxoptforc_old = forc_bfgs_max
end subroutine update_volume_by_bfgs

subroutine update_lattice_by_bfgs()
   integer :: nbox,i,j,k,icount
   real(kind=DP),dimension(3,3) :: tmpforc
   real(kind=DP),allocatable,dimension(:,:) :: coord,forc
   real(kind=DP),dimension(9) :: forc_bfgs
   real(kind=DP) :: forc_bfgs_max
   if(iteration_unit_cell-1<=nhistory_stress)then
      nbox = iteration_unit_cell-1
   else
      nbox = nhistory_stress
   endif
   allocate(coord(9,nbox))
   allocate(forc(9,nbox))
   do i=1,nbox
      icount=0
      do j=1,3
         do k=1,3
            icount = icount + 1
            coord(icount,i) = cellvec_history(j,k,i)
            forc(icount,i)  = stress_forc_history(k,j,i)
         enddo
      enddo
   enddo

   call get_bfgs_force(9,nbox,coord,forc,forc_bfgs,forc_bfgs_max)

   tmpforc=0.0d0
   icount=0
   do i=1,3
      do j=1,3
         icount=icount+1
         tmpforc(j,i) = forc_bfgs(icount)
      enddo
   enddo
   if(printable)then
      write(nfout,'(a)') ' -- BFGS force acting on the unit cell --'
      write(nfout,'(a,3f20.10)')  '    a_vector ',tmpforc(1,1),tmpforc(1,2),tmpforc(1,3)
      write(nfout,'(a,3f20.10)')  '    b_vector ',tmpforc(2,1),tmpforc(2,2),tmpforc(2,3)
      write(nfout,'(a,3f20.10)')  '    c_vector ',tmpforc(3,1),tmpforc(3,2),tmpforc(3,3)
      write(nfout,'(a,f20.10)')   '    max: ',forc_bfgs_max
   endif
   if(maxoptforc_old*10<forc_bfgs_max)then
      if (printable) write(nfout,'(a)') &
      & 'the estimated force acting on the cell seems to be very large; &
      & update will be done by the steepest-descent method'
      do i=1,3
         do j=1,3
            altv(i,j) = altv(i,j) + stress_forc_history(j,i,nbox)*delta_stress
         enddo
      enddo
   else
      do i=1,3
         do j=1,3
            altv(i,j) = altv(i,j) + tmpforc(j,i)
         enddo
      enddo
   endif
   maxoptforc_old = forc_bfgs_max
   deallocate(coord)
   deallocate(forc)
end subroutine update_lattice_by_bfgs

subroutine update_lat_and_coord_by_bfgs( opt_mode )
  integer, intent(in) :: opt_mode

  integer :: nbox,i,j,k,icount
  real(kind=DP), dimension(3,3) :: tmpforc
  real(kind=DP), allocatable,dimension(:,:) :: coord,forc
  real(kind=DP), allocatable :: forc_bfgs(:)
  real(kind=DP) :: forc_bfgs_max
  real(kind=DP), allocatable :: tmpforc_cps(:,:)
  integer :: nmobile, ndim
  real(kind=DP) :: c1, avec(3), bvec(3), cvec(3)
  real(kind=DP), allocatable :: wk_cps(:,:), wk_atom_forc(:,:)
  real(kind=DP), allocatable :: wk_cellvec(:,:), wk_stress_forc(:,:)

  if ( .not. altv0_is_set ) then
     altv0 = altv;  altv0_is_set = .true.
  endif

  if(iteration_unit_cell-1<=nhistory_stress)then
     nbox = iteration_unit_cell-1
  else
     nbox = nhistory_stress
  endif

  allocate( tmpforc_cps(natm,3) )
  nmobile=0
  do i=1, natm
     do j=1,3
        if (imdtypxyz(i,j)==0) cycle
        nmobile = nmobile+1
     enddo
  enddo

  ndim = 9 +nmobile
  allocate(forc_bfgs(ndim));   allocate(coord(ndim,nbox));  allocate(forc(ndim,nbox))

  allocate( wk_cps(natm,3) )
  allocate( wk_atom_forc(natm,3) )
  allocate( wk_cellvec(3,3) )
  allocate( wk_stress_forc(3,3) )

  do i=1,nbox
     icount=0
     if ( opt_mode == 0 ) then
        wk_cps(:,:)         = cps_history(:,:,i)
        wk_atom_forc(:,:)   = atom_forc_history(:,:,i)
        wk_cellvec(:,:)     = cellvec_history(:,:,i)

        avec = wk_cellvec(:,1);    bvec = wk_cellvec(:,2);   cvec = wk_cellvec(:,3)
        c1 = get_volume( avec, bvec, cvec )
        c1 = c1**(1./3.)
        wk_stress_forc(:,:) = stress_forc_history(:,:,i) *c1 *stress_force_mul_factor

     else if ( opt_mode == 1 ) then
        call convert_to_q_and_Atilde( natm, altv0, 1, &
             &           cps_history(:,:,i), atom_forc_history(:,:,i), &
             &           cellvec_history(:,:,i), stress_forc_history(:,:,i), &
             &           wk_cps, wk_atom_forc, wk_cellvec, wk_stress_forc )
     else if ( opt_mode == 2 ) then
        call convert_to_q_and_Atilde( natm, altv0, 2, &
             &           cps_history(:,:,i), atom_forc_history(:,:,i), &
             &           cellvec_history(:,:,i), stress_forc_history(:,:,i), &
             &           wk_cps, wk_atom_forc, wk_cellvec, wk_stress_forc )
     endif

     do j=1,3
        do k=1,3
           icount = icount + 1
           coord(icount,i) = wk_cellvec(j,k)
           forc(icount,i)  = wk_stress_forc(k,j)
        enddo
     enddo
     do j=1, natm
        Do k=1, 3
           if(imdtypxyz(j,k)==0)cycle
           icount = icount +1
           coord(icount,i) = wk_cps( j,k )
           forc(icount,i)  = wk_atom_forc( j,k )
        end Do
     end do
  enddo

  call get_bfgs_force( ndim, nbox, coord, forc, forc_bfgs, forc_bfgs_max )

  tmpforc=0.0d0
  icount=0
  do i=1,3
     do j=1,3
        icount=icount+1
        tmpforc(j,i) = forc_bfgs(icount)
     enddo
  enddo
  do i=1, natm
     do j=1, 3
        if (imdtypxyz(i,j)==0 ) cycle
        icount=icount+1
        tmpforc_cps(i,j) = forc_bfgs(icount)
     end do
  end do

  if(printable)then
     write(nfout,'(a)') ' -- BFGS force acting on the unit cell --'
     write(nfout,'(a,3f20.10)')  '    a_vector ',tmpforc(1,1),tmpforc(1,2),tmpforc(1,3)
     write(nfout,'(a,3f20.10)')  '    b_vector ',tmpforc(2,1),tmpforc(2,2),tmpforc(2,3)
     write(nfout,'(a,3f20.10)')  '    c_vector ',tmpforc(3,1),tmpforc(3,2),tmpforc(3,3)
     write(nfout,'(a,f20.10)')   '    max: ',forc_bfgs_max
  endif
  if (maxoptforc_old*10<forc_bfgs_max) then
     if (printable) write(nfout,'(a)') &
          & 'the estimated force acting on the cell seems to be very large; &
          & update will be done by the steepest-descent method'
     do i=1,3
        do j=1,3
           wk_cellvec(i,j) = wk_cellvec(i,j) &
                &        + wk_stress_forc(j,i) *delta_stress
        enddo
     enddo
     do i=1, natm
        do j=1, 3
           if (imdtypxyz(i,j)==0) cycle
           wk_cps(i,j) = wk_cps(i,j) +wk_atom_forc(i,j)
        end do
     end do
  else
     do i=1,3
        do j=1,3
           wk_cellvec(i,j) = wk_cellvec(i,j) + tmpforc(j,i)
        enddo
     enddo
     do i=1, natm
        do j=1, 3
           if (imdtypxyz(i,j)==0) cycle
           wk_cps(i,j) = wk_cps(i,j) +tmpforc_cps(i,j)
        end do
     end do
  endif

  if ( opt_mode == 0 ) then
     cps(:,:)  = wk_cps(:,:)
     altv(:,:) = wk_cellvec(:,:)
  else if ( opt_mode == 1 ) then
     call convert_from_q_and_Atilde( natm, altv0, 1, &
          &           wk_cps, wk_cellvec, cps, altv )
  else if ( opt_mode == 2 ) then
     call convert_from_q_and_Atilde( natm, altv0, 2, &
          &           wk_cps, wk_cellvec, cps, altv )
  endif

  deallocate( wk_cps );     deallocate( wk_atom_forc )
  deallocate( wk_cellvec ); deallocate( wk_stress_forc )
  deallocate( coord );      deallocate( forc );
  deallocate( forc_bfgs );  deallocate( tmpforc_cps )

  maxoptforc_old = forc_bfgs_max

end subroutine update_lat_and_coord_by_bfgs

  subroutine convert_from_q_and_Atilde( natm, altv0, mode, &
       &           coord_in,  cellvec_in, &
       &           coord_out, cellvec_out )
    integer, intent(in) :: natm, mode
    real(kind=DP), intent(in) :: altv0(3,3)
    real(kind=DP), intent(in) :: coord_in(natm,3)
    real(kind=DP), intent(in) :: cellvec_in(3,3)
    real(kind=DP), intent(out) :: coord_out(natm,3)
    real(kind=DP), intent(out) :: cellvec_out(3,3)
!
    real(kind=DP) :: rltv0_t(3,3), rltv0(3,3)
    real(kind=DP) :: tmp_univol, tmp_rvol, length0(3), fac
    real(kind=DP), allocatable :: wk1(:,:)
!
    call altv_2_rltv( altv0, rltv0, tmp_univol, tmp_rvol )
    rltv0_t = transpose(rltv0)/PAI2

    length0(1) = sqrt( altv0(1,1)**2 +altv0(2,1)**2 +altv0(3,1)**2 )  ! a0
    length0(2) = sqrt( altv0(1,2)**2 +altv0(2,2)**2 +altv0(3,2)**2 )  ! b0
    length0(3) = sqrt( altv0(1,3)**2 +altv0(2,3)**2 +altv0(3,3)**2 )  ! c0

    coord_out(:,:)       = coord_in(:,:)

    fac = omega_for_cellopt *sqrt( dble(natm) )

    if ( mode == 1 ) then
       cellvec_out = cellvec_in
    else if ( mode == 2 ) then
       cellvec_out = cellvec_in /fac
       cellvec_out(:,1) = cellvec_out(:,1) *length0(1)
       cellvec_out(:,2) = cellvec_out(:,2) *length0(2)
       cellvec_out(:,3) = cellvec_out(:,3) *length0(3)
    endif
!
    allocate( wk1(natm,3) )
    call change_of_coordinate_system( rltv0_t, coord_out, natm, natm, wk1 )
    call change_of_coordinate_system( cellvec_out, wk1, natm, natm, coord_out )
    deallocate(wk1)

  end subroutine convert_from_q_and_Atilde

  subroutine convert_to_q_and_Atilde( natm, altv0, mode, &
       &           coord_in,  atom_forc_in,  cellvec_in,  stress_forc_in, &
       &           coord_out, atom_forc_out, cellvec_out, stress_forc_out )
    integer, intent(in) :: natm, mode
    real(kind=DP), intent(in) :: altv0(3,3)
    real(kind=DP), intent(in) :: coord_in(natm,3)
    real(kind=DP), intent(in) :: atom_forc_in(natm,3)
    real(kind=DP), intent(in) :: cellvec_in(3,3)
    real(kind=DP), intent(in) :: stress_forc_in(3,3)
    real(kind=DP), intent(out) :: coord_out(natm,3)
    real(kind=DP), intent(out) :: atom_forc_out(natm,3)
    real(kind=DP), intent(out) :: cellvec_out(3,3)
    real(kind=DP), intent(out) :: stress_forc_out(3,3)
!
    real(kind=DP) :: rltv0_t(3,3), rltv0(3,3), rltv_t(3,3)
    real(kind=DP) :: tmp_altv(3,3), tmp_rltv(3,3)
    real(kind=DP) :: tmp_univol, tmp_rvol, length0(3), fac
    real(kind=DP), allocatable :: wk1(:,:)
!
    call altv_2_rltv( altv0, rltv0, tmp_univol, tmp_rvol )
    rltv0_t = transpose(rltv0)/PAI2

    length0(1) = sqrt( altv0(1,1)**2 +altv0(2,1)**2 +altv0(3,1)**2 )  ! a0
    length0(2) = sqrt( altv0(1,2)**2 +altv0(2,2)**2 +altv0(3,2)**2 )  ! b0
    length0(3) = sqrt( altv0(1,3)**2 +altv0(2,3)**2 +altv0(3,3)**2 )  ! c0

    allocate( wk1(natm,3) )

    coord_out(:,:)       = coord_in(:,:)
    atom_forc_out(:,:)   = atom_forc_in(:,:)
    tmp_altv(:,:)        = cellvec_in(:,:)         ! alat
    stress_forc_out(:,:) = stress_forc_in(:,:) *stress_force_mul_factor

    call altv_2_rltv( tmp_altv, tmp_rltv, tmp_univol, tmp_rvol )
    rltv_t = transpose(tmp_rltv) /PAI2

    call change_of_coordinate_system( rltv_t, coord_out, natm, natm, wk1 )
    call change_of_coordinate_system( altv0, wk1, natm, natm, coord_out )
                                                              ! coord_out = q_i
    call change_of_coordinate_system( rltv0_t, atom_forc_out, natm, natm, wk1 )
    call change_of_coordinate_system( tmp_altv, wk1, natm, natm, atom_forc_out )
!
    fac = omega_for_cellopt *sqrt( dble(natm) )

    if ( mode == 1 ) then
       cellvec_out = tmp_altv
       stress_forc_out = stress_forc_out *univol**(1./3.)
    else if ( mode == 2 ) then
       cellvec_out = fac *tmp_altv
       cellvec_out(:,1) = cellvec_out(:,1) /length0(1)
       cellvec_out(:,2) = cellvec_out(:,2) /length0(2)
       cellvec_out(:,3) = cellvec_out(:,3) /length0(3)
!
       stress_forc_out = stress_forc_out /fac *univol
       stress_forc_out(:,1) = stress_forc_out(:,1) /length0(1)
       stress_forc_out(:,2) = stress_forc_out(:,2) /length0(2)
       stress_forc_out(:,3) = stress_forc_out(:,3) /length0(3)
    endif
    deallocate(wk1)

  end subroutine convert_to_q_and_Atilde

subroutine update_velocities()
   integer :: i,j
   if(lattice_optimization_method==QUENCHED_MD)then
      stress_velocity(:,:) = stress_velocity(:,:) + 0.5*delta_stress * (stress_forc(:,:) + stress_forc_old(:,:))
      do i=1,3
         do j=1,3
            if(stress_velocity(j,i)*stress_forc(j,i)<0) then
               stress_velocity(j,i) = 0.0d0 ! quench
               if(printable) write(nfout,'(a,2i2)') 'quenched stress component : ',j,i
            endif
         enddo
      enddo
   endif
end subroutine update_velocities

subroutine store_volume(vol,dv)
   real(kind=DP),intent(in) :: vol
   real(kind=DP),intent(in) :: dv
   integer :: i
   if(iteration_unit_cell-1<=nhistory_stress)then
      volume_history(iteration_unit_cell-1) = vol
      volume_forc_history(iteration_unit_cell-1) = dv
   else
      do i=2,nhistory_stress
         volume_history(i-1) = volume_history(i)
         volume_forc_history(i-1) = volume_forc_history(i)
      enddo
      volume_history(nhistory_stress) = vol
      volume_forc_history(nhistory_stress) = dv
   endif
end subroutine store_volume

subroutine store_cell_and_stress()
   integer :: ihist,i,nn

   if (iteration_unit_cell-1<=nhistory_stress) then
      cellvec_history(:,:,iteration_unit_cell-1) = altv(:,:)
      stress_forc_history(:,:,iteration_unit_cell-1) = stress_forc(:,:)
      nn = iteration_unit_cell-1

      if ( sw_optimize_coords_sametime == ON ) then
         cps_history(:,:,iteration_unit_cell-1) = cps(:,:)
         atom_forc_history(:,:,iteration_unit_cell-1) = forc_l(:,:)
      endif

   else
      do i=2,nhistory_stress
         cellvec_history(:,:,i-1) = cellvec_history(:,:,i)
         stress_forc_history(:,:,i-1) = stress_forc_history(:,:,i)
      enddo
      cellvec_history(:,:,nhistory_stress) = altv(:,:)
      stress_forc_history(:,:,nhistory_stress) = stress_forc(:,:)
      nn = nhistory_stress

      if ( sw_optimize_coords_sametime == ON ) then
         do i=2,nhistory_stress
            cps_history(:,:,i-1) = cps_history(:,:,i)
            atom_forc_history(:,:,i-1) = atom_forc_history(:,:,i)
         end do
         cps_history(:,:,nhistory_stress) = cps(:,:)
         atom_forc_history(:,:,nhistory_stress) = forc_l(:,:)
      endif
   endif

   if (printable.and.ipriunitcell>=2) then
      write(nfout,'(a)') 'cellvec history'
      do i=1,nn
         write(nfout,'(a,i5)') 'step no ',i
         write(nfout,'(a,3f20.10)') 'avec ',cellvec_history(1,1,i),cellvec_history(2,1,i),cellvec_history(3,1,i)
         write(nfout,'(a,3f20.10)') 'bvec ',cellvec_history(1,2,i),cellvec_history(2,2,i),cellvec_history(3,2,i)
         write(nfout,'(a,3f20.10)') 'cvec ',cellvec_history(1,3,i),cellvec_history(2,3,i),cellvec_history(3,3,i)
      enddo
   endif
end subroutine store_cell_and_stress

subroutine update_lattice()
   integer :: i,j
   if(lattice_optimization_method==QUENCHED_MD)then
      do i=1,3
         do j=1,3
            altv(i,j) = altv(i,j) + delta_stress*stress_velocity(j,i) + 0.5*stress_forc(j,i)*delta_stress*delta_stress
         enddo
      enddo
   else if (lattice_optimization_method==STEEPEST_DESCENT)then
      do i=1,3
         do j=1,3
            altv(i,j) = altv(i,j) + stress_forc(j,i)*delta_stress
         enddo
      enddo
   endif
end subroutine update_lattice

function get_force(altv_in,stress) result(ret)
   real(kind=DP), dimension(3,3), intent(in) :: stress, altv_in
   real(kind=DP), dimension(3,3) :: ret
   integer :: i
   ret = 0.0d0
   do i=1,3
! =========================== KT_mod ================== 13.0B
!      ret(1,i) = dot_product(stress(i,:)-external_stress(i,:),altv(:,1)) ! a-vec
!      ret(2,i) = dot_product(stress(i,:)-external_stress(i,:),altv(:,2)) ! b-vec
!      ret(3,i) = dot_product(stress(i,:)-external_stress(i,:),altv(:,3)) ! c-vec
!
      ret(1,i) = dot_product(stress(i,:),altv_in(:,1)) ! a-vec
      ret(2,i) = dot_product(stress(i,:),altv_in(:,2)) ! b-vec
      ret(3,i) = dot_product(stress(i,:),altv_in(:,3)) ! c-vec
! ====================================================== 13.0B
   enddo
end function get_force

function get_max(stress) result(res)
   real(kind=DP), dimension(3,3), intent(in) :: stress
   integer :: i,j
   real(kind=DP) :: res
   res=dabs(stress(1,1))
   do i=1,3
      do j=1,3
         if(dabs(stress(j,i))>res) res = dabs(stress(j,i))
      enddo
   enddo
end function get_max

! === KT_add ==== 13.1AS
function projected_stress_along_latvecs(altv_in,stress) result(ret)
  real(kind=DP), dimension(3,3), intent(in) :: stress, altv_in
  real(kind=DP), dimension(3,3) :: ret

  integer :: i
  real(kind=DP) :: normalized_lat_vec(3,3), c1
!
  Do i=1, 3
     c1 = sqrt( altv_in(1,i)**2 +altv_in(2,i)**2 + altv_in(3,i)**2 )
     normalized_lat_vec(:,i) = altv_in(:,i) / c1
  End do

  ret = 0.0d0

  do i=1,3
     ret(1,i) = dot_product(stress(i,:),normalized_lat_vec(:,1)) ! a-vec
     ret(2,i) = dot_product(stress(i,:),normalized_lat_vec(:,2)) ! b-vec
     ret(3,i) = dot_product(stress(i,:),normalized_lat_vec(:,3)) ! c-vec
  enddo
end function projected_stress_along_latvecs

! =========================== KT_add ======================== 13.0B &13.1AS
function constrain_on_force_type1(altv_in,f_in) result(f_out)
  real(kind=DP), intent(in) :: f_in(3,3), altv_in(3,3)
  real(kind=DP) :: f_out(3,3)

  integer :: i
  real(kind=DP) :: c1, c2, c3
  real(kind=DP) :: normalized_lat_vec(3,3), vectmp(3), vectmp2(3)
!
  f_out = f_in
!
  Do i=1, 3
     c1 = sqrt( altv_in(1,i)**2 +altv_in(2,i)**2 + altv_in(3,i)**2 )
     normalized_lat_vec(:,i) = altv_in(:,i) / c1
  End do

! -----------------------------------
  if ( fix_angle_alpha .and. fix_angle_beta .and. fix_angle_gamma ) then
     Do i=1, 3
        c1 = dot_product( f_out(i,:), normalized_lat_vec(:,i) )
        f_out(i,:) = c1 *normalized_lat_vec(:,i)
     End do
     return
  endif

! -----------------------------------
  if ( fix_angle_alpha ) then
     vectmp(:) = normalized_lat_vec(:,2) + normalized_lat_vec(:,3)

     Do i=2, 3
        c1 = dot_product( vectmp(:), normalized_lat_vec(:,i) )
        vectmp2(:) = vectmp(:) - c1*normalized_lat_vec(:,i)

        c2 = sqrt( vectmp2(1)**2 + vectmp2(2)**2 + vectmp2(3)**2 )
        vectmp2 = vectmp2 /c2

        c3 = dot_product( f_out(i,:), vectmp2(:) )
        f_out(i,:) = f_out(i,:) - c3 *vectmp2(:)
     End do
  endif

! -----------------------------------
  if ( fix_angle_beta ) then
     vectmp(:) = normalized_lat_vec(:,1) + normalized_lat_vec(:,3)

     Do i=1, 3, 2
        c1 = dot_product( vectmp(:), normalized_lat_vec(:,i) )
        vectmp2(:) = vectmp(:) - c1*normalized_lat_vec(:,i)

        c2 = sqrt( vectmp2(1)**2 + vectmp2(2)**2 + vectmp2(3)**2 )
        vectmp2 = vectmp2 /c2

        c3 = dot_product( f_out(i,:), vectmp2(:) )
        f_out(i,:) = f_out(i,:) - c3 *vectmp2(:)
     End do
  endif

! -----------------------------------
  if ( fix_angle_gamma ) then
     vectmp(:) = normalized_lat_vec(:,1) + normalized_lat_vec(:,2)

     Do i=1, 2
        c1 = dot_product( vectmp(:), normalized_lat_vec(:,i) )
        vectmp2(:) = vectmp(:) - c1*normalized_lat_vec(:,i)

        c2 = sqrt( vectmp2(1)**2 + vectmp2(2)**2 + vectmp2(3)**2 )
        vectmp2 = vectmp2 /c2

        c3 = dot_product( f_out(i,:), vectmp2(:) )
        f_out(i,:) = f_out(i,:) - c3 *vectmp2(:)
     End do

  endif

end function constrain_on_force_type1
! ================================================================ 13.0B

! === KT_add === 13.1AS
function constrain_on_force_type2(altv_in,f_in) result(f_out)
  real(kind=DP), intent(in) :: f_in(3,3), altv_in(3,3)
  real(kind=DP) :: f_out(3,3)

  integer :: i
  real(kind=DP) :: c1, c2, c3
  real(kind=DP) :: normalized_lat_vec(3,3), vectmp(3), vectmp2(3)
!
  f_out = f_in
!
  Do i=1, 3
     c1 = sqrt( altv_in(1,i)**2 +altv_in(2,i)**2 + altv_in(3,i)**2 )
     normalized_lat_vec(:,i) = altv_in(:,i) / c1
  End do

! -----------------------------------
  if ( fix_length_a ) then
     c1 = dot_product( f_out(1,:), normalized_lat_vec(:,1) )
     f_out(1,:) = f_out(1,:) -c1 *normalized_lat_vec(:,1)
  endif
  if ( fix_length_b ) then
     c1 = dot_product( f_out(2,:), normalized_lat_vec(:,2) )
     f_out(2,:) = f_out(2,:) -c1 *normalized_lat_vec(:,2)
  endif
  if ( fix_length_c ) then
     c1 = dot_product( f_out(3,:), normalized_lat_vec(:,3) )
     f_out(3,:) = f_out(3,:) -c1 *normalized_lat_vec(:,3)
  endif
end function constrain_on_force_type2
! ====== 13.1AS

subroutine init_thermostats()
  integer :: i,j,ir,iat
  integer :: icnstrnt_typ, count_nat
  nrsv_t = 1
  allocate(thermostats(nrsv_t))
  do i=1,nrsv_t
    thermostats(i)%m_thermo = qmass(i)
    thermostats(i)%target_temperature = tkb(i)/CONST_kB
    thermostats(i)%s_thermo = 1.d0
    thermostats(i)%p_thermo = 0.d0
    count_nat=0
    do j=1,natm
       ir = icnstrnt_typ(imdtyp(j),T_CONTROL)
       if (ir == i) count_nat = count_nat+1
    enddo
    allocate(thermostats(i)%target_atoms(count_nat))
    thermostats(i)%natm = count_nat
    thermostats(i)%gfree = 3*count_nat
!    if(thermostats(i)%gfree<1) thermostats(i)%gfree = 3*natm
    count_nat=0
    do j=1,natm
       ir = icnstrnt_typ(imdtyp(j),T_CONTROL)
       if (ir == i) then
         count_nat = count_nat+1
         thermostats(i)%target_atoms(count_nat) = j
       endif
    enddo
    thermostats(i)%p_baro=0.d0
    thermostats(i)%p_baro_tmp=0.d0
  enddo
  if(printable .and. ipriunitcell>=1) then
    do i=1,nrsv_t
      write(nfout,'(a,i5)') '!** status of thermostat no ',i
      write(nfout,'(a,f20.3)') '!** temperature            ',thermostats(i)%target_temperature
      write(nfout,'(a,f20.3)') '!** mass                   ',thermostats(i)%m_thermo
      write(nfout,'(a,i8)')    '!** number of target atoms ',thermostats(i)%natm
      write(nfout,'(8i8)')   (thermostats(i)%target_atoms(j), j=1,thermostats(i)%natm)
    enddo
  endif
  return
end subroutine init_thermostats

subroutine m_UnitCell_init_md()
  integer :: i,j,ir,iat
  real(kind=DP), dimension(3) :: avec,bvec,cvec
  integer :: icnstrnt_typ, count_nat
  if(.not.first_call) return
  allocate(p_l(natm,3));p_l = 0.d0
  allocate(p_l_lattice(natm,3));p_l_lattice = 0.d0
  allocate(forc_l_lattice(natm,3));forc_l_lattice = 0.d0
  allocate(p_l_metric(natm,3));p_l_metric = 0.d0
  if(p_ctrl_method==METRIC_TENSOR) then
  metric = matmul(transpose(altv),altv)
  else
  metric = transpose(altv)
  endif
  call init_thermostats()
  call m_UnitCell_init_md_per_step()
  do i=1,nrsv_t
    do j=1,thermostats(i)%natm
      iat = thermostats(i)%target_atoms(j)
      p_l(iat,:) = amion(ityp(iat))*cpd_l(iat,:) * thermostats(i)%s_thermo
    enddo
  enddo
  !do i=1,natm
  !   p_l(i,:) = amion(ityp(i)) * cpd_l(i,:) * s_thermo
  !enddo
  if(p_ctrl_method==METRIC_TENSOR) then
  call convert_coords(natm,p_l,altv,p_l_lattice)
  else
  p_l_lattice = p_l
  endif
  call convert_coords(natm,p_l_lattice,metric_inv,p_l_metric)
  do i=1,3
     do j=1,3
        factor(j,i) = 1.d0
        if(i.eq.j) factor(j,i) = 0.5d0
     enddo
  enddo
  first_call = .false.
end subroutine m_UnitCell_init_md

subroutine build_dudmetric()
  real(kind=DP), dimension(3,3) :: stress_tensor
  real(kind=DP) :: avstress
  integer :: i
  stress_tensor = m_Stress_get_curr_stress()
  if(sw_average_diagonal == ON) then
    avstress = 0.d0
    do i=1, 3
      avstress = avstress+stress_tensor(i,i)
    enddo
    avstress = avstress/3.d0
    do i=1, 3
      stress_tensor(i,i) = avstress
    enddo
  endif
  if(p_ctrl_method==METRIC_TENSOR) then
  dudmetric = -univol * 0.5d0 * matmul(altv_inv,matmul(stress_tensor,transpose(altv_inv)))
  else
  dudmetric = -univol * 0.5d0 * get_force(altv,stress_tensor)
  endif
end subroutine build_dudmetric

subroutine convert_coords(natm,inar,cmat,outar)
  integer, intent(in) :: natm
  real(kind=DP), intent(in) :: inar(natm,3)
  real(kind=DP), intent(in) :: cmat(3,3)
  real(kind=DP), intent(out) :: outar(natm,3)
  integer :: i,j
  outar = 0.d0
  do i=1,3
     do j=1,natm
        outar(j,i)  = dot_product(cmat(i,:),inar(j,:))
     enddo
  enddo
end subroutine convert_coords

function get_kinetic_atom_thermo(ithermo,natm,p_atm,mat) result(ret)
  integer, intent(in) :: ithermo,natm
  real(kind=DP), intent(in), dimension(natm,3) :: p_atm
  real(kind=DP), intent(in), dimension(3,3) :: mat
  real(kind=DP) :: ret
  real(kind=DP),dimension(3) :: tmpp
  integer :: i,j,k,iat
  ret = 0.d0
  do i=1,thermostats(ithermo)%natm
    iat = thermostats(ithermo)%target_atoms(i)
    if(p_ctrl_method==METRIC_TENSOR)then
    do j=1,3
      tmpp(j) = dot_product(mat(j,:),p_atm(iat,:))
    enddo
    else
      tmpp(:) = p_atm(iat,:)
    endif
    ret = ret + dot_product(p_atm(iat,:),tmpp(:))/amion(ityp(iat))
  enddo
  ret = 0.5d0 * ret/(thermostats(ithermo)%s_thermo*thermostats(ithermo)%s_thermo)
end function get_kinetic_atom_thermo

!function get_kinetic_atom(natm,p_atm, mat) result(ret)
!  integer, intent(in) :: natm
!  real(kind=DP), intent(in), dimension(natm,3) :: p_atm
!  real(kind=DP), intent(in), dimension(3,3) :: mat
!  real(kind=DP) :: ret
!  real(kind=DP),dimension(3) :: tmpp
!  integer :: i,j
!  ret = 0.d0
!  do i=1,natm
!     do j=1,3
!        tmpp(j) = dot_product(mat(j,:),p_atm(i,:))
!     enddo
!     ret = ret + dot_product(p_atm(i,:),tmpp(:))/amion(ityp(i))
!  enddo
!  ret = 0.5d0 * ret/(s_thermo*s_thermo)
!end function get_kinetic_atom

function get_kinetic_thermo_thermo(ithermo) result(ret)
  integer, intent(in) :: ithermo
  real(kind=DP) :: ret
  ret = 0.5d0 * thermostats(ithermo)%p_thermo * thermostats(ithermo)%p_thermo / thermostats(ithermo)%m_thermo
end function get_kinetic_thermo_thermo

function get_kinetic_thermo(ithermo,p_thermo) result (ret)
  integer, intent(in) :: ithermo
  real(kind=DP), intent(in) :: p_thermo
  real(kind=DP) :: ret
  ret = 0.5d0 * p_thermo * p_thermo /thermostats(ithermo)%m_thermo
end function get_kinetic_thermo

function get_kinetic_baro(p_baro) result(ret)
  real(kind=DP), intent(in), dimension(3,3) :: p_baro
  real(kind=DP) :: ret
  integer :: i
  ret = 0.d0
  do i=1,3
     ret = ret + dot_product(p_baro(i,:),p_baro(:,i))
  enddo
  ret = 0.5d0 * ret/(m_baro * univol * univol)
end function get_kinetic_baro

function get_potential_baro() result(ret)
  real(kind=DP) :: ret
  ret = target_pressure * univol
end function get_potential_baro

subroutine m_UnitCell_init_md_per_step()
  real(kind=DP), external :: deter3
  if(p_ctrl_method==METRIC_TENSOR) then
  univol = sqrt(deter3(metric))
  else
  univol=dabs(deter3(altv))
  endif
  call inver3n(3,metric,metric_inv)
  call inver3n(3,altv,altv_inv)
  call build_dudmetric()
  H0 = get_H0()
end subroutine m_UnitCell_init_md_per_step

function tensor_rank2_equals(t1,t2) result(ret)
  real(kind=DP), dimension(3,3) :: t1,t2
  logical :: ret
  integer :: i,j
  do i=1,3
     do j=1,3
        if(abs(t1(j,i)-t2(j,i)).gt.eps) then
           ret = .false.
           return
        endif
     enddo
  enddo
  ret = .true.
end function tensor_rank2_equals

subroutine metric_to_altv()
  real(kind=DP) :: r1,r2,r3,cosa,cosb,cosc,sina,sinb,sinc

  if(p_ctrl_method==METRIC_TENSOR) then
  altv = 0.d0

  r1=dsqrt(metric(1,1))
  r2=dsqrt(metric(2,2))
  r3=dsqrt(metric(3,3))

  cosa=metric(2,1)/(r2*r1)
  sina=dsqrt(dabs(1.d0-cosa**2))
  cosb=metric(3,1)/(r3*r1)
  sinb=dsqrt(dabs(1d0-cosb**2))
  cosc=metric(3,2)/(r3*r2)
  sinc=dsqrt(dabs(1d0-cosc**2))

  altv(1,1)=r1
  altv(1,2)=r2*cosa
  altv(2,2)=r2*sina
  altv(1,3)=r3*cosb
  altv(2,3)=(metric(3,2)-r3*r2*cosa*cosb)/(r2*sina)
  altv(3,3)=dsqrt(dabs(metric(3,3)-altv(2,3)**2-altv(1,3)**2))

  else
  altv = transpose(metric)
  endif

end subroutine metric_to_altv

function get_econst_volume() result(res)
  real(kind=DP) :: res
  integer :: i
  res = etotal
  do i=1,natm
    res = res+0.5d0*amion(ityp(i))*dot_product(cpd_l(i,:),cpd_l(i,:))
  enddo
  res = res+0.5d0*pv_baro*pv_baro/m_baro+target_pressure*vol
  if(imdalg==PT_CONTROL .and. t_ctrl_method==NOSE_HOOVER) then
    do i=1,nrsv_t
      res = res + 0.5d0*thermostats(i)%m_thermo*(thermostats(i)%p_thermo*thermostats(i)%s_thermo)**2 &
                + thermostats(i)%gfree*log(thermostats(i)%s_thermo)*thermostats(i)%target_temperature*CONST_kB
    enddo
  endif
  return
end function get_econst_volume

function m_UC_get_hamiltonian() result(res)
  real(kind=DP) :: res
  real(kind=DP) :: k_atm,k_baro,k_thermo,pot_baro,pot_thermo
  logical, save :: first_call = .true.
  integer :: it
  if(p_ctrl_method==VOLUME) then
!    res = get_econst_volume()
    res = econst
    return
  endif
  call m_UnitCell_init_md()
  k_atm = 0.d0
  do it=1,nrsv_t
    k_atm = k_atm+get_kinetic_atom_thermo(it,natm,p_l_lattice,metric_inv)
  enddo
  k_baro = 0.d0
  do it=1,nrsv_t
    k_baro = k_baro + get_kinetic_baro(matmul(metric,thermostats(it)%p_baro))
  enddo
  pot_baro = get_potential_baro()
  if(control_pressure == OFF) then
    k_baro = 0.d0
    pot_baro = 0.d0
  endif
  if(imdalg == PT_CONTROL .and. t_ctrl_method /= VELOCITY_SCALING)then
     k_thermo = 0.d0
     pot_thermo = 0.d0
     do it=1,nrsv_t
       k_thermo = k_thermo + get_kinetic_thermo_thermo(it)
       pot_thermo = pot_thermo &
                  + thermostats(it)%gfree * thermostats(it)%target_temperature*CONST_kB*log(thermostats(it)%s_thermo)
     enddo
     !k_thermo = get_kinetic_thermo(p_thermo)
     !pot_thermo = gfree * target_temperature * CONST_kB * log(s_thermo)
  else if (imdalg == P_CONTROL .or. (imdalg == PT_CONTROL .and. t_ctrl_method == VELOCITY_SCALING)) then
     k_thermo = 0.d0
     pot_thermo = 0.d0
  endif
!  res = s_thermo * (k_atm + etotal + k_thermo + pot_thermo + k_baro + pot_baro)
  res = k_atm + etotal + k_thermo + pot_thermo + k_baro + pot_baro
  if(.not.H0_has_been_set)then
    H0 = res
    H0_has_been_set = .true.
  endif
  res = res-H0
end function m_UC_get_hamiltonian

function get_H0() result (res)
  real(kind=DP) :: res
  real(kind=DP) :: k_atm,k_baro,k_thermo,pot_baro,pot_thermo
  integer :: it
  k_atm = 0.d0
  do it=1,nrsv_t
    k_atm = k_atm + get_kinetic_atom_thermo(it,natm,p_l_lattice,metric_inv)
  enddo
  k_baro = 0.d0
  do it=1,nrsv_t
    k_baro = k_baro + get_kinetic_baro(matmul(metric,thermostats(it)%p_baro))
  enddo
  pot_baro = get_potential_baro()
  if(control_pressure == OFF) then
    k_baro = 0.d0
    pot_baro = 0.d0
  endif
  if(imdalg == PT_CONTROL .and. t_ctrl_method /= VELOCITY_SCALING)then
     k_thermo = 0.d0
     pot_thermo = 0.d0
     do it=1,nrsv_t
       k_thermo = k_thermo+get_kinetic_thermo_thermo(it)
       pot_thermo = pot_thermo &
                  + thermostats(it)%gfree * thermostats(it)%target_temperature*CONST_kB*log(thermostats(it)%s_thermo)
     enddo
     !k_thermo = get_kinetic_thermo(p_thermo)
     !pot_thermo = gfree * target_temperature * CONST_kB * log(s_thermo)
  else if (imdalg == P_CONTROL .or. (imdalg == PT_CONTROL .and. t_ctrl_method == VELOCITY_SCALING)) then
     k_thermo = 0.d0
     pot_thermo = 0.d0
  endif
!  res = s_thermo * (k_atm + etotal + k_thermo + pot_thermo + k_baro + pot_baro)
  res = k_atm + etotal + k_thermo + pot_thermo + k_baro + pot_baro
  return
end function get_H0

function get_curr_pressure_from_for(forc) result(res)
  real(kind=DP), dimension(natm,3), intent(in) :: forc
  integer :: i
  real(kind=DP) :: res, rf
  real(kind=DP), dimension(3,3) :: stress_tensor
  res = 0.d0
  rf = 0.d0
  do i=1,natm
    res = res+amion(ityp(i))*dot_product(cpd_l(i,:),cpd_l(i,:))+dot_product(cps(i,:),forc(i,:))
  enddo
  res = res/univol/3.d0
end function get_curr_pressure_from_for

function get_curr_pressure_from_stress() result(res)
  integer :: i
  real(kind=DP) :: res
  real(kind=DP), dimension(3,3) :: stress_tensor
  res = 0.d0
  do i=1,natm
    res = res+amion(ityp(i))*dot_product(cpd_l(i,:),cpd_l(i,:))
  enddo
  res = res/univol

  stress_tensor = m_Stress_get_curr_stress()
  do i=1,3
    res = res + stress_tensor(i,i)
  enddo
  res = res/3.d0
  return
end function get_curr_pressure_from_stress

function get_curr_temperature() result(res)
  integer :: i
  real(kind=DP) :: res
  res = 0.d0
  do i=1,natm
    res = res+amion(ityp(i))*dot_product(cpd_l(i,:),cpd_l(i,:))
  enddo
  res = res/dble(3*natm)/const_kB
  return
end function get_curr_temperature

function m_UC_get_curr_pressure() result(res)
  real(kind=DP) :: curr_pressure, res
  integer :: i,it
  real(kind=DP), dimension(3,3) :: stress
  if(p_ctrl_method==VOLUME .or. p_ctrl_method==LATTICE_VECTOR) then
    res = get_curr_pressure_from_stress()
    return
  endif
  curr_pressure = 0.d0
  do it=1,nrsv_t
    curr_pressure = curr_pressure+2.d0*get_kinetic_atom_thermo(it,natm,p_l_lattice,metric_inv)
  enddo
  !!curr_pressure = 2.d0 * get_kinetic_atom(natm,p_l_lattice,metric_inv)
!  if(p_ctrl_method==METRIC_TENSOR) then
  do i=1,3
    curr_pressure = curr_pressure - 2.d0 * dudmetric(i,i) * altv(i,i) * altv(i,i)
  enddo
!  else
!  write(nfout,'(a,f20.10)') 'currpressure ke ',curr_pressure/3.d0/univol
!  stress = m_Stress_get_curr_stress()
!  do i=1,3
!    curr_pressure = curr_pressure - 2.d0 * stress(i,i)
!  enddo
!  write(nfout,'(a,f20.10)') 'currpressure ke+stress ',curr_pressure/3.d0/univol
!  endif
  res = curr_pressure/3.d0/univol
end function m_UC_get_curr_pressure

function m_UC_get_curr_temperature() result(res)
  integer :: it
  real(kind=DP) :: ke,res
  if(p_ctrl_method==VOLUME) then
    res = get_curr_temperature()
    return
  endif
  ke = 0.d0
  do it=1,nrsv_t
    ke = ke+ 2.d0 * get_kinetic_atom_thermo(it,natm,p_l_lattice,metric_inv)/dble(thermostats(it)%gfree)
  enddo
!  res = ke/gfree/CONST_kB
  res = ke/CONST_kB
end function m_UC_get_curr_temperature

subroutine m_UC_rd_md_cntn_data(nfcntn)
   integer, intent(in) :: nfcntn
   integer :: i,ierr
   logical             :: tag_is_found, EOF_reach
   integer,      parameter   :: len_str = 132
   character(len=len_str)       :: str

   call m_UnitCell_init_md()

   if(mype==0)then
      call rewind_to_tag0(nfcntn,len(tag_pressure_control), &
            &  tag_pressure_control, &
            &  EOF_reach, tag_is_found, str,len_str)
      if(tag_is_found)then
         do i=1,nrsv_t
           read(nfcntn,*)
           read(nfcntn,*) thermostats(i)%s_thermo
           read(nfcntn,*)
           read(nfcntn,*) thermostats(i)%p_thermo
         enddo
         read(nfcntn,*)
         read(nfcntn,*) H0
         call read_33('altv',altv)
         call read_33('altv_inv',altv_inv)
         call read_33('metric',metric)
         call read_33('metric_inv',metric_inv)
         call read_33('dudmetric',dudmetric)
         do i=1,nrsv_t
           call read_33('p_baro',thermostats(i)%p_baro)
           call read_33('p_baro_tmp',thermostats(i)%p_baro_tmp)
         enddo
         call read_natm3('forc_l_lattice',forc_l_lattice,allocated(forc_l_lattice))
         call read_natm3('p_l',p_l,allocated(p_l))
         call read_natm3('p_l_lattice',p_l_lattice,allocated(p_l_lattice))
         call read_natm3('p_l_metric',p_l_metric,allocated(p_l_metric))
      endif
   endif

   if(npes>1)then
      call mpi_bcast(tag_is_found,1,mpi_logical,0,MPI_CommGroup,ierr)
      if(tag_is_found) then
         do i=1,nrsv_t
           call mpi_bcast(thermostats(i)%s_thermo,1,mpi_double_precision,0,MPI_CommGroup,ierr)
           call mpi_bcast(thermostats(i)%p_thermo,1,mpi_double_precision,0,MPI_CommGroup,ierr)
         enddo
         call mpi_bcast(H0,1,mpi_double_precision,0,MPI_CommGroup,ierr)
         call mpi_bcast(altv,9,mpi_double_precision,0,MPI_CommGroup,ierr)
         call mpi_bcast(altv_inv,9,mpi_double_precision,0,MPI_CommGroup,ierr)
         call mpi_bcast(metric,9,mpi_double_precision,0,MPI_CommGroup,ierr)
         call mpi_bcast(metric_inv,9,mpi_double_precision,0,MPI_CommGroup,ierr)
         do i=1,nrsv_t
           call mpi_bcast(thermostats(i)%p_baro,9,mpi_double_precision,0,MPI_CommGroup,ierr)
           call mpi_bcast(thermostats(i)%p_baro_tmp,9,mpi_double_precision,0,MPI_CommGroup,ierr)
         enddo
         if(allocated(forc_l_lattice)) &
         &  call mpi_bcast(forc_l_lattice,natm*3,mpi_double_precision,0,MPI_CommGroup,ierr)
         if(allocated(p_l)) &
         &  call mpi_bcast(p_l,natm*3,mpi_double_precision,0,MPI_CommGroup,ierr)
         if(allocated(p_l_lattice)) &
         &  call mpi_bcast(p_l_lattice,natm*3,mpi_double_precision,0,MPI_CommGroup,ierr)
         if(allocated(p_l_metric)) &
         &  call mpi_bcast(p_l_metric,natm*3,mpi_double_precision,0,MPI_CommGroup,ierr)
         call altv_2_rltv(altv,rltv,univol,rvol)  ! in b_CS
         call change_of_coordinate_system(altv,pos,natm,natm,cps) !-(b_I.S.) pos -> cps
         call primitive2bravais(nfout,p2bmat,altv(:,1),altv(:,2),altv(:,3),a,b,c,ca,cb,cc,il) ! in b_CS
         if(iteration_unit_cell>1) then
           H0_has_been_set = .true.
         endif
      endif
   endif

   contains

   subroutine read_33(key,tensor)
      character(len=*), intent(in) :: key
      real(kind=DP), intent(out), dimension(3,3) :: tensor
      read(nfcntn,*)
      read(nfcntn,*) tensor(1,1),tensor(2,1),tensor(3,1)
      read(nfcntn,*) tensor(1,2),tensor(2,2),tensor(3,2)
      read(nfcntn,*) tensor(1,3),tensor(2,3),tensor(3,3)
   end subroutine read_33

   subroutine read_natm3(key,ar,allcd)
      character(len=*), intent(in) :: key
      real(kind=DP), intent(out), dimension(natm,3) :: ar
      logical, intent(in) :: allcd
      integer :: i
      read(nfcntn,*)
      if(allcd)then
         do i=1,natm
            read(nfcntn,*) ar(i,1),ar(i,2),ar(i,3)
         enddo
      else
         do i=1,natm
            read(nfcntn,*)
         enddo
      endif
   end subroutine read_natm3

end subroutine m_UC_rd_md_cntn_data

subroutine m_UC_wd_md_cntn_data(nfcntn)
   integer, intent(in) :: nfcntn
   integer :: i
   if(first_call)then
      call m_UnitCell_init_md()
   endif
   if(mype==0)then
      write(nfcntn,*) tag_pressure_control
      do i=1,nrsv_t
        write(nfcntn,*) 's_thermo'
        write(nfcntn,'(3f25.15)') thermostats(i)%s_thermo
        write(nfcntn,*) 'p_thermo'
        write(nfcntn,'(3f25.15)') thermostats(i)%p_thermo
      enddo
      write(nfcntn,*) 'H0'
      write(nfcntn,'(3f25.15)') H0
      call dump_33('altv',altv)
      call dump_33('altv_inv',altv_inv)
      call dump_33('metric',metric)
      call dump_33('metric_inv',metric_inv)
      call dump_33('dudmetric',dudmetric)
      do i=1,nrsv_t
        call dump_33('p_baro',thermostats(i)%p_baro)
        call dump_33('p_baro_tmp',thermostats(i)%p_baro_tmp)
      enddo
      call dump_natm3('forc_l_lattice',forc_l_lattice,allocated(forc_l_lattice))
      call dump_natm3('p_l',p_l,allocated(p_l))
      call dump_natm3('p_l_lattice',p_l_lattice,allocated(p_l_lattice))
      call dump_natm3('p_l_metric',p_l_metric,allocated(p_l_metric))
   endif

   contains

   subroutine dump_33(key,tensor)
      character(len=*), intent(in) :: key
      real(kind=DP), intent(in), dimension(3,3) :: tensor
      write(nfcntn,*) key
      write(nfcntn,'(3f25.15)') tensor(1,1),tensor(2,1),tensor(3,1)
      write(nfcntn,'(3f25.15)') tensor(1,2),tensor(2,2),tensor(3,2)
      write(nfcntn,'(3f25.15)') tensor(1,3),tensor(2,3),tensor(3,3)
   end subroutine dump_33

   subroutine dump_natm3(key,ar,allcd)
      character(len=*), intent(in) :: key
      real(kind=DP), intent(in), dimension(natm,3) :: ar
      logical, intent(in) :: allcd
      integer :: i
      write(nfcntn,*) key
      if(allcd)then
         do i=1,natm
            write(nfcntn,'(3f25.15)') ar(i,1),ar(i,2),ar(i,3)
         enddo
      else
         do i=1,natm
            write(nfcntn,'(3f25.15)') 0.d0,0.d0,0.d0
         enddo
      endif
   end subroutine dump_natm3

end subroutine m_UC_wd_md_cntn_data

! perform npt-md based on the metric tensor and Nose-Poincare thermostat
! original theory from I. Souza and J. L. Martins PRB 55 8733
! formulation from E. Herna´ndez JCP 115 10282
subroutine m_UnitCell_md(forc_l)
  use m_Ionic_System, only : temperature_profile, set_new_temperature, set_new_pressure, sw_temperature_profile &
                           , m_IS_add_friction_force, m_IS_add_random_force, pressure_profile, sw_pressure_profile
  real(kind=DP), intent(in), dimension(natm,3) :: forc_l
  real(kind=DP), dimension(3,3) :: p33,p_baro33,p_baro33_2,p_baro_tmp_m,p_baro_tmp2
  real(kind=DP), dimension(3,3) :: metric_tmp,metric_tmp2
  integer :: i,j
  real(kind=DP) :: k_baro,k_atm,pot_s,pot_baro,a,b,c,k_baro2,uvol2t,s_thermo_tmp,k_thermo,pot_thermo
  real(kind=DP), external :: deter3

  if (p_ctrl_method == VOLUME) then
    call uniform_md(forc_l)
    return
  endif
  call m_UnitCell_init_md()

  call m_UnitCell_init_md_per_step()
  if(sw_temperature_profile==ON) then
    call set_new_temperature(nfout,nrsv_t,temperature_profile)
    do i=1,nrsv_t
      thermostats(i)%m_thermo = qmass(i)
      thermostats(i)%target_temperature = tkb(i)/CONST_kB
    enddo
  endif
  if(sw_pressure_profile==ON) then
    call set_new_pressure(nfout,pressure_profile)
  endif
  call second_half()
  call first_half()

  call post()

  contains

  subroutine post()
    integer :: i,j,iat
    real(kind=DP),dimension(3,3) :: stress
    if(p_ctrl_method==METRIC_TENSOR) then
    call convert_coords(natm,p_l_lattice,altv_inv,p_l)
    else
    p_l = p_l_lattice
    endif
    do i=1,nrsv_t
      do j=1,thermostats(i)%natm
        iat = thermostats(i)%target_atoms(j)
        cpd_l(iat,:) = p_l(iat,:)/amion(ityp(iat))/thermostats(i)%s_thermo
      enddo
    enddo
    !do i=1,natm
    !   cpd_l(i,:) = p_l(i,:)/amion(ityp(i))/s_thermo
    !enddo
  end subroutine post

  subroutine first_half()
    integer :: j, it, iat, iatm
    real(kind=DP) :: dotp

    if(imdalg == PT_CONTROL .and. t_ctrl_method == VELOCITY_SCALING) then
      if(p_ctrl_method==METRIC_TENSOR) then
      call convert_coords(natm,p_l_lattice,altv_inv,p_l)
      else
      p_l = p_l_lattice
      endif
      do it=1,nrsv_t ! thermo loop
        do j=1,thermostats(it)%natm
          iat = thermostats(it)%target_atoms(j)
          cpd_l(iat,:) = p_l(iat,:)/amion(ityp(iat))/thermostats(it)%s_thermo
        enddo
      enddo
      call scale_velocity()
      do it=1,nrsv_t ! thermo loop
        do j=1,thermostats(it)%natm
          iat = thermostats(it)%target_atoms(j)
          p_l(iat,:) = amion(ityp(iat)) * cpd_l(iat,:) * thermostats(it)%s_thermo
        enddo
      enddo
      if(p_ctrl_method==METRIC_TENSOR) then
      call convert_coords(natm,p_l,altv,p_l_lattice)
      else
      p_l_lattice = p_l
      endif
    endif

    do it=1,nrsv_t ! thermo loop
      do iatm=1,thermostats(it)%natm
        iat = thermostats(it)%target_atoms(iatm)
        ! eq(18a)
        p_l_lattice(iat,:) = p_l_lattice(iat,:) + 0.5d0 * dtio * thermostats(it)%s_thermo * forc_l_lattice(iat,:)
      enddo
      ! eq(18b)
      !if(p_ctrl_method==METRIC_TENSOR)then
      call convert_coords(natm,p_l_lattice,metric_inv,p_l_metric)
      !else
      !p_l_metric = p_l_lattice
      !endif
      if(control_pressure == ON) then
        p33 = 0.d0
        do iatm=1,thermostats(it)%natm
          iat = thermostats(it)%target_atoms(iatm)
          do i=1,3
             do j=1,3
                dotp = p_l_metric(iat,j)*p_l_metric(iat,i)/amion(ityp(iat))/(thermostats(it)%s_thermo*thermostats(it)%s_thermo)
                p33(j,i) = p33(j,i) + 0.5d0 * dotp
             enddo
          enddo
        enddo
        thermostats(it)%p_baro_tmp = thermostats(it)%p_baro
        do i=1,nmax
          p_baro_tmp_m = matmul(metric,thermostats(it)%p_baro_tmp)
          p_baro33 = matmul(thermostats(it)%p_baro_tmp,p_baro_tmp_m)/(m_baro * univol * univol)
          k_baro = get_kinetic_baro(p_baro_tmp_m)
          pot_baro = 0.5d0 * get_potential_baro()
          p_baro_tmp2(:,:) = thermostats(it)%p_baro(:,:) - 0.5d0 * dtio * thermostats(it)%s_thermo * &
          & (dudmetric(:,:) + 2.d0 * factor(:,:)*(-p33(:,:)+p_baro33(:,:) + &
          &  metric_inv(:,:)*(pot_baro - k_baro)))
          if(tensor_rank2_equals(thermostats(it)%p_baro_tmp,p_baro_tmp2)) then
             thermostats(it)%p_baro_tmp = p_baro_tmp2
             exit
          endif
          thermostats(it)%p_baro_tmp = p_baro_tmp2
          if(i .eq. nmax) then
             write(nfout,'(a)') 'failed to reach convergence for the iteration of the barostat momentum'
             call flush(nfout)
             call phase_error_with_msg(nfout,'failed to reach convergence for the iteration of the barostat momentum'&
                                      ,__LINE__,__FILE__)
          endif
        enddo
      endif

      !eq(18c,d)
      if(imdalg == PT_CONTROL .and. t_ctrl_method /= VELOCITY_SCALING) then
         k_atm = get_kinetic_atom_thermo(it,natm,p_l_lattice,metric_inv)
         p_baro_tmp_m = matmul(metric,thermostats(it)%p_baro_tmp)
         k_baro = get_kinetic_baro(p_baro_tmp_m)
         pot_baro = get_potential_baro()
         if(control_pressure == OFF) then
             k_baro = 0.d0
             pot_baro = 0.d0
         endif
         pot_thermo = (1.d0+log(thermostats(it)%s_thermo)) * thermostats(it)%gfree &
                    & * thermostats(it)%target_temperature * CONST_kB
         a = 0.25d0 * dtio/thermostats(it)%m_thermo
         b = 1.d0
         c = 0.5 * dtio * (-k_atm + pot_thermo + etotal + k_baro + pot_baro-H0) - thermostats(it)%p_thermo
         thermostats(it)%p_thermo = (-b+dsqrt(b*b-4.d0*a*c))/(2.d0*a)
         s_thermo_tmp = thermostats(it)%s_thermo * (1.d0+0.5d0*dtio*thermostats(it)%p_thermo/thermostats(it)%m_thermo)/ &
                                                   (1.d0-0.5d0*dtio*thermostats(it)%p_thermo/thermostats(it)%m_thermo)
      else if (imdalg == P_CONTROL .or. (imdalg == PT_CONTROL .and. t_ctrl_method == VELOCITY_SCALING)) then
         s_thermo_tmp = 1.d0
         thermostats(it)%p_thermo = 0.d0
      endif

      if(control_pressure ==  ON)then
        !eq(18e)
        metric_tmp = metric
        p_baro33 = matmul(metric,matmul(thermostats(it)%p_baro_tmp,metric))/(m_baro*univol * univol)
        do i=1,nmax
           if(p_ctrl_method==METRIC_TENSOR)then
           uvol2t = dabs(deter3(metric_tmp))
           else
           uvol2t = dabs(deter3(metric_tmp))*dabs(deter3(metric_tmp))
           endif
           p_baro33_2 = matmul(metric_tmp,matmul(thermostats(it)%p_baro_tmp,metric_tmp))/(m_baro*uvol2t)
           metric_tmp2(:,:) = metric(:,:) + mobility_cell(:,:)*factor(:,:)*dtio &
           &                *(thermostats(it)%s_thermo*p_baro33(:,:) +s_thermo_tmp * p_baro33_2(:,:))
           if(tensor_rank2_equals(metric_tmp2,metric_tmp)) then
              metric_tmp = metric_tmp2
              exit
           endif
           metric_tmp = metric_tmp2
           if(i .eq. nmax) then
              write(nfout,'(a)') 'failed to reach convergence for the metric iteration'
              call flush(nfout)
              call phase_error_with_msg(nfout,'failed to reach convergence for the metric iteration',__LINE__,__FILE__)
           endif
        enddo
        metric = metric_tmp
      endif

      !eq(18f)
      !do i=1,natm
      do iatm=1,thermostats(it)%natm
         iat = thermostats(it)%target_atoms(iatm)
         pos(iat,:) = pos(iat,:) + 0.5d0*dtio*(p_l_metric(iat,:)/(amion(ityp(iat))*thermostats(it)%s_thermo) +  &
                    &                          p_l_metric(iat,:)/(amion(ityp(iat))*s_thermo_tmp))
      enddo

      thermostats(it)%s_thermo = s_thermo_tmp
    enddo ! thermo loop

    if(control_pressure==ON)then
    call metric_to_altv()
    call inver3n(3,metric,metric_inv)
    call inver3n(3,altv,altv_inv)
    call altv_2_rltv(altv,rltv,univol,rvol)  ! in b_CS
    call change_of_coordinate_system(altv,pos,natm,natm,cps) !-(b_I.S.) pos -> cps
    call primitive2bravais(nfout,p2bmat,altv(:,1),altv(:,2),altv(:,3),a,b,c,ca,cb,cc,il) ! in b_CS
    endif
    if(mype==0)then
      write(nflatconst,'(i10,7f20.10)') iteration_unit_cell,a,b,c, &
      &     acos(ca)/PAI*180.d0,acos(cb)/PAI*180.d0,acos(cc)/PAI*180.d0,univol
      write(nfmetric,'(i10)')     iteration_unit_cell
      write(nfmetric,'(3f20.10)') metric(1,1),metric(1,2),metric(1,3)
      write(nfmetric,'(3f20.10)') metric(2,1),metric(2,2),metric(2,3)
      write(nfmetric,'(3f20.10)') metric(3,1),metric(3,2),metric(3,3)
      call m_Files_flush_pcontrol_files()
    endif
  end subroutine first_half

  subroutine second_half()
    integer :: it,iatm,iat
    real(kind=DP) :: ptherm,dotp
    real(kind=DP), allocatable, dimension(:,:) :: forc
    !eq(18g)
    if(imdalg == PT_CONTROL)then
       k_atm = 0.d0
       do it=1,nrsv_t
         k_atm = k_atm + get_kinetic_atom_thermo(it,natm,p_l_lattice,metric_inv)
       enddo
       k_baro = 0.d0
       do it=1,nrsv_t
         p_baro_tmp_m = matmul(metric,thermostats(it)%p_baro_tmp)
         k_baro = k_baro + get_kinetic_baro(p_baro_tmp_m)
       enddo
       pot_baro = get_potential_baro()
       if(control_pressure == OFF) then
          k_baro = 0.d0
          pot_baro = 0.d0
       endif
       pot_thermo = 0.d0
       k_thermo  = 0.d0
       do it=1,nrsv_t
         ptherm     = (1.d0+log(thermostats(it)%s_thermo)) * thermostats(it)%gfree * thermostats(it)%target_temperature * CONST_kB
         pot_thermo = pot_thermo+ptherm
         k_thermo   = k_thermo+get_kinetic_thermo(it,ptherm)
         thermostats(it)%p_thermo = thermostats(it)%p_thermo &
         + 0.5d0*dtio*(k_atm - etotal - k_thermo - pot_thermo - k_baro - pot_baro + H0)
       enddo
    endif

    if(control_pressure == ON)then
      do it=1,nrsv_t
        !eq(18h)
        p_baro_tmp_m = matmul(metric,thermostats(it)%p_baro_tmp)
        p_baro33 = matmul(thermostats(it)%p_baro_tmp,p_baro_tmp_m)/(m_baro*univol * univol)
        k_baro = get_kinetic_baro(p_baro_tmp_m)
        pot_baro = 0.5d0 * get_potential_baro()
        thermostats(it)%p_baro = thermostats(it)%p_baro_tmp
        p33 = 0.d0
        do iatm=1,thermostats(it)%natm
          iat = thermostats(it)%target_atoms(iatm)
          do i=1,3
             do j=1,3
                dotp = p_l_metric(iat,i)*p_l_metric(iat,j)/amion(ityp(iat))/(thermostats(it)%s_thermo*thermostats(it)%s_thermo)
                !p33(j,i) = p33(j,i)+0.5d0*dot_product(p_l_metric(:,i),p_l_metric(:,j)/amion(ityp(:))/(s_thermo*s_thermo))
                p33(j,i) = p33(j,i)+dotp
             enddo
          enddo
        enddo
        thermostats(it)%p_baro(:,:) = thermostats(it)%p_baro(:,:) - 0.5d0 * dtio * thermostats(it)%s_thermo * &
        & (dudmetric(:,:) + 2.d0 * factor(:,:)*(-p33(:,:)+p_baro33(:,:) +     &
        &  metric_inv(:,:)*(pot_baro - k_baro)))
      enddo
    endif

    !eq(18i)
    if(t_ctrl_method==LANGEVIN) then
      allocate(forc(natm,3));forc=forc_l
      call m_IS_add_friction_force(forc)
      call m_IS_add_random_force(forc)
      if(p_ctrl_method==METRIC_TENSOR)then
      call convert_coords(natm,forc,altv,forc_l_lattice)
      else
      forc_l_lattice = forc
      endif
      deallocate(forc)
    else
      if(p_ctrl_method==METRIC_TENSOR)then
      call convert_coords(natm,forc_l,altv,forc_l_lattice)
      else
      forc_l_lattice = forc_l
      endif
    endif
    do it=1,nrsv_t
      do iatm=1,thermostats(it)%natm
        iat = thermostats(it)%target_atoms(iatm)
        p_l_lattice(iat,:) = p_l_lattice(iat,:) + 0.5d0 * dtio * thermostats(it)%s_thermo * forc_l_lattice(iat,:)
      enddo
    enddo

  end subroutine second_half

end subroutine m_UnitCell_md

subroutine init_uniform_md()
  logical, save :: first_call = .true.
  integer :: i
  if(.not.first_call) return
  pv_baro = 0.d0
  pv_baro_old = 0.d0
  vol     = univol
  vol_old = univol
  latv    = altv
!  call init_thermostats()
  if(ipriunitcell>=2) write(nfout,'(a,3f20.10)') '!** initial vol_old vol ',vol_old,vol
  first_call = .false.
end subroutine init_uniform_md

subroutine uniform_md(forc_l)
  use m_Ionic_System, only : temperature_profile, set_new_temperature, set_new_pressure, sw_temperature_profile &
                           , m_IS_add_friction_force, m_IS_add_random_force, pressure_profile,sw_pressure_profile
  real(kind=DP), dimension(natm,3), intent(in) :: forc_l
  real(kind=DP), dimension(:,:), allocatable :: forc
  integer :: i
  call init_uniform_md()
  allocate(forc(natm,3));forc=0.d0
  if(t_ctrl_method==LANGEVIN) then
      call m_IS_add_friction_force(forc)
      call m_IS_add_random_force(forc)
  endif
  call second_half()
  if(sw_temperature_profile==ON) then
    call set_new_temperature(nfout,nrsv_t,temperature_profile)
    do i=1,nrsv_t
      thermostats(i)%m_thermo = qmass(i)
      thermostats(i)%target_temperature = tkb(i)/CONST_kB
    enddo
  endif
  if(sw_pressure_profile==ON) then
    call set_new_pressure(nfout,pressure_profile)
  endif
  if(t_ctrl_method==VELOCITY_SCALING) call scale_velocity()
  call first_half()
  deallocate(forc)

  contains

  subroutine first_half()
    integer :: i,j,ia,na,ii
    real(kind=DP) :: pressure, pthermot, pthermot2
    real(kind=DP), allocatable, dimension(:,:) :: cpdt
    real(kind=DP) :: sfactor
    sfactor = 1.d0
    if(imdalg==PT_CONTROL.and. t_ctrl_method==NOSE_HOOVER) sfactor = thermostats(1)%s_thermo
    pressure = get_curr_pressure_from_stress()
    pv_baro_old = pv_baro
    pv_baro = pv_baro+0.5d0*dtio*(pressure-target_pressure)*pv_baro
    vol_old = vol
    vol = vol+dtio*(pv_baro/m_baro)*sfactor
    do i=1,natm
!      cpd_l(i,:) = (cpd_l(i,:)+dtio*0.5*forc_l(i,:)/amion(ityp(i)))/(1.d0+0.5*dtio*pv_baro/(3.d0*m_baro*vol))
      cpd_l(i,:) = cpd_l(i,:)+dtio*0.5*forc_l(i,:)/amion(ityp(i))
    enddo
    if(t_ctrl_method==LANGEVIN)then
      do i=1,natm
        cpd_l(i,:) = cpd_l(i,:)+dtio*0.5*forc(i,:)/amion(ityp(i))
      enddo
    endif
!    do i=1,natm
!      cpd_l(i,:) = cpd_l(i,:)+dtio*0.5*forc(i,:)/amion(ityp(i))
!    enddo

    if(imdalg==PT_CONTROL .and. t_ctrl_method==NOSE_HOOVER) then
      allocate(cpdt(natm,3))
      do i=1,nrsv_t
        na = thermostats(i)%natm
        cpdt = cpd_l
        pthermot = thermostats(i)%p_thermo
        pthermot2 = 0.d0
        do ii=1,nmax
          pthermot = thermostats(i)%p_thermo+0.5d0*dtio*(2.d0/thermostats(i)%m_thermo)*(get_ke(i,cpdt)-get_nkbt(i))
          do j=1,na
            ia = thermostats(i)%target_atoms(j)
            cpdt(ia,:) = cpd_l(ia,:)-cpdt(ia,:)*dtio*0.5d0*pthermot
          enddo
          if(abs(pthermot-pthermot2)<eps)then
            cpd_l = cpdt
            thermostats(i)%p_thermo = pthermot
            exit
          endif
          pthermot2 = pthermot
        enddo
      enddo
      deallocate(cpdt)
    endif

    do i=1,natm
      cps(i,:) = (cps(i,:)+0.5d0*dtio*(2*cpd_l(i,:)+cps(i,:)*pv_baro/(3.d0*m_baro*vol))) &
               & /(1.d0-0.5d0*dtio*pv_baro/(3.d0*m_baro*vol))
    enddo
    call uniform_update_of_altv(vol/vol_old)
    if(mype==0)then
      write(nflatconst,'(i10,7f20.10)') iteration_unit_cell,a,b,c, &
      &     acos(ca)/PAI*180.d0,acos(cb)/PAI*180.d0,acos(cc)/PAI*180.d0,univol
      !write(nfmetric,'(i10)')     iteration_unit_cell
      !write(nfmetric,'(3f20.10)') latv(1,1),latv(1,2),latv(1,3)
      !write(nfmetric,'(3f20.10)') latv(2,1),latv(2,2),latv(2,3)
      !write(nfmetric,'(3f20.10)') latv(3,1),latv(3,2),latv(3,3)
      call m_Files_flush_pcontrol_files()
    endif
  end subroutine first_half

  subroutine second_half()
    integer :: i,j,na,ia,ii
    real(kind=DP) :: pressure,pthermot,pthermot2
    real(kind=DP), allocatable, dimension(:,:) :: cpdt
    real(kind=DP) :: sfactor
    sfactor = 1.d0
    if(imdalg==PT_CONTROL.and. t_ctrl_method==NOSE_HOOVER) sfactor = thermostats(1)%s_thermo
    pressure = get_curr_pressure_from_stress()
    pv_baro = pv_baro+0.5d0*dtio*(pressure-target_pressure)*sfactor
    do i=1,natm
      !cpd_l(i,:) = cpd_l(i,:)+0.5d0*dtio*(forc(i,:)/amion(ityp(i))-cpd_l(i,:)*pv_baro/(3.d0*vol*m_baro))
!      cpd_l(i,:) = (cpd_l(i,:)+dtio*0.5*forc_l(i,:)/amion(ityp(i)))/(1.d0+0.5*dtio*pv_baro/(3.d0*m_baro*vol))
      cpd_l(i,:) = cpd_l(i,:)+dtio*0.5*forc_l(i,:)/amion(ityp(i))
    enddo
    if(t_ctrl_method==LANGEVIN)then
      do i=1,natm
        cpd_l(i,:) = cpd_l(i,:)+dtio*0.5*forc(i,:)/amion(ityp(i))
      enddo
    endif
    if(imdalg==PT_CONTROL .and. t_ctrl_method==NOSE_HOOVER) then
      allocate(cpdt(natm,3))
      do i=1,nrsv_t
        na = thermostats(i)%natm
        cpdt = cpd_l
        pthermot = thermostats(i)%p_thermo
        pthermot2 = 0.d0
        do ii=1,nmax
          pthermot = thermostats(i)%p_thermo+0.5d0*dtio*(2.d0/thermostats(i)%m_thermo)*(get_ke(i,cpdt)-get_nkbt(i))
          do j=1,na
            ia = thermostats(i)%target_atoms(j)
            cpdt(ia,:) = cpd_l(ia,:)-cpdt(ia,:)*dtio*0.5d0*pthermot
          enddo
          if(abs(pthermot-pthermot2)<eps)then
            cpd_l = cpdt
            thermostats(i)%p_thermo = pthermot
            exit
          endif
          pthermot2 = pthermot
        enddo
        thermostats(i)%s_thermo = thermostats(i)%s_thermo/(1.d0-dtio*thermostats(i)%p_thermo)
      enddo
      deallocate(cpdt)
    endif
    if(.not.H0_has_been_set) then
      H0 = get_econst_volume()
      H0_has_been_set = .true.
    endif
    econst = H0-get_econst_volume()
  end subroutine second_half

  subroutine uniform_update_of_altv(volfac)
    real(kind=DP), intent(in) :: volfac
    real(kind=DP) :: vecfac,rltv_t(3,3)
    integer :: i
    vecfac = volfac**(1.d0/3.d0)
    latv = latv*vecfac
    altv = latv
    call altv_2_rltv(altv,rltv,univol,rvol)  ! in b_CS
    rltv_t = transpose(rltv)/PAI2
    call change_of_coordinate_system(rltv_t,cps,natm,natm,pos) !-(b_I.S.) cps -> pos
    call primitive2bravais(nfout,p2bmat,altv(:,1),altv(:,2),altv(:,3),a,b,c,ca,cb,cc,il) ! in b_CS
  end subroutine uniform_update_of_altv

  function get_ke(ithermo,cpd_l) result(ret)
    integer, intent(in) :: ithermo
    real(kind=DP), dimension(natm,3), intent(in) :: cpd_l
    integer :: i,nat,ia
    real(kind=DP) :: ret
    nat = thermostats(ithermo)%natm
    ret = 0.d0
    do i=1,nat
      ia = thermostats(ithermo)%target_atoms(i)
      ret = ret+amion(ityp(ia))*dot_product(cpd_l(ia,:),cpd_l(ia,:))
    enddo
    ret = 0.5d0*ret
  end function get_ke

  function get_nkbt(ithermo) result(res)
    integer, intent(in) :: ithermo
    real(kind=DP) :: res
    res = 0.5d0*thermostats(ithermo)%gfree*thermostats(ithermo)%target_temperature*CONST_kB
  end function get_nkbt

end subroutine uniform_md

end module m_UnitCell

logical function unitcell_can_change()
  use m_Control_Parameters, only : sw_optimize_lattice,imdalg,sw_stress_correction
  use m_Const_Parameters, only : ON,P_CONTROL, PT_CONTROL
  use m_Stress, only : m_Stress_in_correction
  implicit none
  if(imdalg == P_CONTROL .or. imdalg == PT_CONTROL) then
     unitcell_can_change = .true.
     return
  endif
!  if(m_Stress_in_correction()) then
!     unitcell_can_change = .true.
!     return
!  endif
  if(sw_stress_correction==ON)then
    unitcell_can_change = .true.
    return
  endif
  unitcell_can_change = sw_optimize_lattice==ON
end function unitcell_can_change

