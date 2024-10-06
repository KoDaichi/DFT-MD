!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 604 $)
!
!  SUBROUINE: constrained_dynamics
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
!!!!!BRANCH_P ORG_Parallel
#ifndef DISABLE_CONSTRAINTS
! ==============================================================================
subroutine constrained_dynamics()
  use m_velocity_verlet, only : m_vv_increment_md_step, m_vv_do_dynamics, &
   &  m_vv_get_curr_md_step
  use m_constraints, only : &
   &  m_cnstr_get_n_reac_coords, m_cnstr_set_reac_coord, &
   &  m_cnstr_get_reac_coords_s
  use m_Control_Parameters, only : printable, icond, icond_org, terminated_because
  use m_Files, only : nfout
  use m_Parallelization, only : nrank_conf, mype_conf, conf_para
  use m_Const_Parameters, only : STREVL_ITERATION,REAC_ITERATION, MAX_MDSTEPS_REACHED
  use mpi
  implicit none

!  include 'mpif.h'

  integer :: stat
  integer :: ireac,istart,iend,ierr
  logical :: initialization_required
  
  call Preparation(0)                          ! Basis set, symmetry check etc.
  call Preparation_for_mpi(1)     ! mpi
  call PseudoPotential_Construction
#ifdef ENABLE_ESM_PACK
  if(initialization_required())then
     call Preparation_for_ESM
  endif
#endif
  call Ewald_and_Structure_Factor
  call Initial_Electronic_Structure

  call Initial_MD_Condition()

  icond=icond_org

!!$  call constrained_dynamics_init(nfout)

  istart=m_cnstr_get_reac_coords_s()
  if(conf_para) istart=1
  iend=m_cnstr_get_n_reac_coords()
  reaction_coordinate: do ireac=istart,iend
    if(conf_para.and.mod(ireac-1,nrank_conf)/=mype_conf) cycle reaction_coordinate
    if(m_cnstr_set_reac_coord(ireac)) cycle reaction_coordinate
    md_step: do while(m_vv_increment_md_step())
      if(printable) write(nfout,'(a,i7)') 'entering MD step no. ',m_vv_get_curr_md_step()
      stat=0
      call m_vv_do_dynamics(stat)
      if (stat==1) then ! relaxation converged
        call postproc_md_step(printable,nfout)
        exit md_step
      else if(stat==2)then ! end of program
        call postproc_reac_coord(stat)
        exit reaction_coordinate
      endif
      call postproc_md_step(printable,nfout)
      if(printable) write(nfout,'(a,i7)') 'done MD step no.     ',m_vv_get_curr_md_step()
      call Checkpoint_File(STREVL_ITERATION)
    enddo md_step
    call Stress()
    call postproc_reac_coord(stat)
    call Checkpoint_File(REAC_ITERATION)
  enddo reaction_coordinate
  if (ireac > iend ) terminated_because =  MAX_MDSTEPS_REACHED  ! Calculation of all reaction coordinates has been completed.

  call mpi_barrier(mpi_comm_world, ierr)

  call finalize()

end subroutine constrained_dynamics

subroutine postproc_md_step(printable,nfout)
  use m_constraints, only : m_cnstr_print_status
  use m_velocity_verlet, only : m_vv_io_per_step,m_vv_get_curr_md_step
  use m_Control_Parameters, only : postproc_frequency
  implicit none
  logical,intent(in) :: printable
  integer, intent(in) :: nfout
  if(printable) write(nfout,'(a)') 'doing post-proc for each MD step'
  call m_cnstr_print_status()
  call m_vv_io_per_step()
  if(postproc_frequency<=0) return
  if(mod(m_vv_get_curr_md_step(),postproc_frequency)==0)then
     call Postprocessing(.true.)
  endif
end subroutine postproc_md_step

subroutine postproc_reac_coord(stat)
  use m_constraints, only : m_cnstr_pp_reac_coords
  implicit none
  integer,intent(in) :: stat
  call m_cnstr_pp_reac_coords(stat)
end subroutine postproc_reac_coord

subroutine finalize()
  use m_routines, only : close_all_opened_files
  implicit none
  call scf_finalize()
  call constrained_dynamics_finalize()
  call close_all_opened_files()
end subroutine finalize

subroutine constrained_dynamics_dump()
  use m_Files, only : F_CNTN, nfcntn
  use m_Parallelization, only : mype
  use m_velocity_verlet, only : m_vv_wd_vv_variables
  use m_constraints, only : m_cnstr_wd_constraints
  implicit none
  if (mype/=0) return
  open(nfcntn,&
& file=F_CNTN,status='unknown',position='append')
  call m_vv_wd_vv_variables(nfcntn)
  call m_cnstr_wd_constraints(nfcntn)
  close(nfcntn)
end subroutine constrained_dynamics_dump

subroutine constrained_dynamics_finalize()
  use m_velocity_verlet, only : m_vv_finalize
  use m_constraints, only : m_cnstr_finalize
  use m_variables_for_atoms, only : m_vfa_dealloc
  implicit none
  call constrained_dynamics_dump()
  call m_vv_finalize()
  call m_cnstr_finalize()
  call m_vfa_dealloc()
end subroutine constrained_dynamics_finalize

subroutine constrained_dynamics_init(nfout)
  use m_variables_for_atoms, only : m_vfa_alloc_atoms
  use m_variables_for_dynamics, only : m_vfd_init
  use m_Files, only :  F_ZAJ_in, F_CNTN_in, F_CNTN_BIN_in, F_CHGT_in &
 &   , F_ZAJ,F_CNTN,F_CNTN_BIN,F_CHGT, nfinp, F_INP, nfcntn, open0
  use m_Const_Parameters, only : QUENCHED_MD, T_CONTROL, VERLET, DAMPED_MD, CONTINUATION, INITIAL &
 &   , formatted, check_file_name_on, old, ON, VELOCITY_SCALING
  use m_constraints, only : m_cnstr_constraints_exist, m_cnstr_initialize, &
 &      m_cnstr_parse_input, m_cnstr_rd_constraints, &
 &      F_ZAJ_org, F_CNTN_BIN_org, F_CNTN_org, F_CHGT_org
  use m_velocity_verlet, only : m_vv_init, m_vv_parse_input, &
 &      m_vv_rd_vv_variables
  use m_Parallelization, only : conf_para, sw_wdir, workdir
  use m_Control_Parameters, only: icond,imdalg,printable
  use mpi
  implicit none
  integer, intent(in) :: nfout
  integer ::  f_selectTop, f_selectBlock, f_openInputFile, f_closeInputFile
  integer :: iiret,i
  logical :: ex,op
  logical, allocatable, dimension(:) :: ltmp
  logical, allocatable, dimension(:) :: ltmp0
  character(len=2) :: cid
  character(len=2) :: cid_tmp
  integer :: mpierr
  integer :: ifile
  integer :: nline
  logical :: statfile_exists
  integer :: ireac
  logical :: fin

!  include 'mpif.h'

  call m_vfa_alloc_atoms()
  call m_vfd_init()

  inquire(unit=nfinp,opened=op)
  if(op)close(nfinp,status='keep')
  if(sw_wdir == ON) then
    iiret = f_openInputFile(trim(workdir)//F_INP)
  else
    iiret = f_openInputFile(F_INP)
  end if
  if( iiret < 0) then
    !if(printable) write(nfout,'(" !!! Error in opening of inputfile: ",a32)') F_INP
    !stop
    call phase_error_with_msg(nfout,'!!! Error in opening of inputfile:'//trim(F_INP),__LINE__,__FILE__)
  else if(iiret > 0) then
    if(printable) write(nfout,'(" !!! There is something wrong in the input file: ",a32)') F_INP
  end if

  iiret = f_selectTop()
  iiret = f_selectBlock('structure')

  call m_cnstr_parse_input()
  if(m_cnstr_constraints_exist()) call m_cnstr_initialize()


  call m_vv_init()
  call m_vv_parse_input()

  iiret = f_selectTop()

  if(imdalg.ne.QUENCHED_MD.and. &
  &  imdalg.ne.T_CONTROL  .and. &
  &  imdalg.ne.VERLET     .and. &
  &  imdalg.ne.DAMPED_MD  .and. &
  &  imdalg.ne.VELOCITY_SCALING) then
    if(printable) write(nfout,'(a,i5,a)') '**WARN imdalg :      ',imdalg, &
                  ' not implemented (yet) for the constraints driver'
    if(printable) write(nfout,'(a,i5)')   '**WARN defaults to : ',QUENCHED_MD
    imdalg = QUENCHED_MD
    !!$stop
  endif

  F_ZAJ_org      = F_ZAJ_in
  F_CNTN_org     = F_CNTN_in
  F_CNTN_BIN_org = F_CNTN_BIN_in
  F_CHGT_org     = F_CHGT_in

!  iiret = f_closeInputFile()
  call open0(nfinp, F_INP, 'F_INP     ',    old,   formatted,check_file_name_on)

  if(icond==CONTINUATION .and. .not. conf_para)then
    call m_vv_rd_vv_variables(nfcntn,.true.)
    call m_cnstr_rd_constraints(nfcntn)
  endif

end subroutine constrained_dynamics_init

subroutine prepare_continuation_of_reac(nfcntn,sk)
  use m_Files, only : m_Files_reopen_nfcntn
  use m_Ionic_System, only : m_IS_rd_pos_and_v
  use m_velocity_verlet, only : m_vv_rd_vv_variables
  use m_constraints, only : m_cnstr_rd_constraints
  use m_IterationNumbers, only : m_Iter_rd_iteration_numbers
  use m_Total_Energy, only : m_TE_rd_total_energy
  use m_Control_Parameters, only : m_CtrlP_rd_isolver,icond
  implicit none
  integer, intent(in) :: nfcntn
  logical, intent(in) :: sk
  integer :: i
  call m_Files_reopen_nfcntn()
  if(.not.sk)then
  call m_TE_rd_total_energy(nfcntn)
  call m_CtrlP_rd_isolver(nfcntn)
  call m_Iter_rd_iteration_numbers(nfcntn,icond)
  endif
  call m_IS_rd_pos_and_v(nfcntn)
  call m_vv_rd_vv_variables(nfcntn,sk)
  call m_cnstr_rd_constraints(nfcntn)
end subroutine prepare_continuation_of_reac
! === "phase" repository merge: To make 3D_Parallel by tkato 2011/11/10 ========
#endif
!!!!!BRANCH_P_END ORG_Parallel
! ==============================================================================
