!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 622 $)
!
!  SUBROUINE: meta_dynamics, init_mtd
!
!  AUTHOR(S): J. Koga 2010
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
! ==============================================================================
#ifndef DISABLE_CONSTRAINTS
subroutine meta_dynamics()
  use m_velocity_verlet, only : m_vv_increment_md_step, m_vv_do_dynamics, m_vv_get_curr_md_step
  use m_meta_dynamics, only : m_mtd_update_bias, m_mtd_finalize
  use m_Const_Parameters, only : CONTINUATION
  use m_Control_Parameters, only : printable,icond
  use m_Files, only : nfout, F_CNTN, F_CNTN_BIN,F_CHGT,F_ZAJ,F_CNTN_bak,F_CNTN_BIN_bak,F_CHGT_bak,F_ZAJ_bak &
 &                  , F_CNTN_in, F_CNTN_BIN_in, F_CHGT_in, F_ZAJ_in &
 &                  , file_existence_contfiles, file_existence_3contfiles &
 &                  , m_Files_reopen_nfcntn, m_Files_reopen_nfcntn_bin, m_Files_reopen_nfzaj, m_Files_reopen_nfchgt
  use m_Parallelization, only : mype_conf,MPI_CommGroup
  use mpi

  implicit none

!  include 'mpif.h'

  integer :: stat, ierr
  character(len=256) :: to_string
  logical :: initialization_required

  call Preparation(0)                          ! Basis set, symmetry check etc.
  call Preparation_for_mpi(1)     ! mpi
  call init_mtd()
  call PseudoPotential_Construction
#ifdef ENABLE_ESM_PACK
  if(initialization_required())then
     call Preparation_for_ESM
  endif
#endif
  call Ewald_and_Structure_Factor
  call Initial_Electronic_Structure
  call Initial_MD_Condition()

  if(mype_conf>0 .and. icond==CONTINUATION .and. .not. &
 &  (file_existence_contfiles.and.file_existence_3contfiles))then
    F_CNTN = trim(F_CNTN_bak)//'_conf'//trim(to_string(mype_conf,2))
    F_CNTN_BIN = trim(F_CNTN_BIN_bak)//'_conf'//trim(to_string(mype_conf,2))
    F_CHGT = trim(F_CHGT_bak)//'_conf'//trim(to_string(mype_conf,2))
    F_ZAJ = trim(F_ZAJ_bak)//'_conf'//trim(to_string(mype_conf,2))

    F_CNTN_in = trim(F_CNTN_bak)//'_conf'//trim(to_string(mype_conf,2))
    F_CNTN_BIN_in = trim(F_CNTN_BIN_bak)//'_conf'//trim(to_string(mype_conf,2))
    F_CHGT_in = trim(F_CHGT_bak)//'_conf'//trim(to_string(mype_conf,2))
    F_ZAJ_in = trim(F_ZAJ_bak)//'_conf'//trim(to_string(mype_conf,2))

    if(printable)then
      write(nfout,'(a)')
      write(nfout,'(a)') 'restored continuation files'
      write(nfout,'(a)') '  F_CNTN     : '//trim(F_CNTN)
      write(nfout,'(a)') '  F_CNTN_BIN : '//trim(F_CNTN_BIN)
      write(nfout,'(a)') '  F_CHGT     : '//trim(F_CHGT)
      write(nfout,'(a)') '  F_ZAJ      : '//trim(F_ZAJ)
      write(nfout,'(a)')
    endif

    call m_Files_reopen_nfcntn()
    call m_Files_reopen_nfcntn_bin()
    call m_Files_reopen_nfzaj()
    call m_Files_reopen_nfchgt()
  endif

  bias:do
    if(.not.m_mtd_update_bias())then
      exit bias
    endif
    md_step: do while(m_vv_increment_md_step())
      if(printable) write(nfout,'(a,i7)') 'entering MD step no. ',m_vv_get_curr_md_step()
      stat=0
      call m_vv_do_dynamics(stat)
      if(stat==2)then
        exit bias
      endif
      call postproc_md_step_mtd(printable,nfout)
      if(printable) write(nfout,'(a,i7)') 'done MD step no.     ',m_vv_get_curr_md_step()
      if(printable) write(nfout,'(a)')    ''
    enddo md_step
    call postproc_bias()
  enddo bias

  call mpi_barrier(MPI_CommGroup, ierr)
  call mpi_barrier(mpi_comm_world, ierr)

  call finalize_mtd()

end subroutine meta_dynamics

subroutine postproc_md_step_mtd(printable,nfout)
  use m_velocity_verlet, only : m_vv_io_per_step
  implicit none
  logical, intent(in) :: printable
  integer, intent(in) :: nfout
  if(printable) write(nfout,'(a)') 'doing post-proc for each MD step'
  call m_vv_io_per_step()
end subroutine postproc_md_step_mtd

subroutine postproc_bias()
end subroutine postproc_bias

subroutine finalize_mtd()
  use m_routines, only : close_all_opened_files
  implicit none
  call close_all_opened_files()
  call scf_finalize()
  call dump_mtd()
end subroutine finalize_mtd

subroutine dump_mtd()
  use m_Files, only : F_CNTN, nfcntn
  use m_Parallelization, only : mype
  use m_velocity_verlet, only : m_vv_wd_vv_variables
  use m_meta_dynamics, only : m_mtd_wd_mtd
  implicit none
  if(mype/=0)return
  open(nfcntn,&
& file=F_CNTN,status='unknown',position='append')
  call m_mtd_wd_mtd(nfcntn)
  call m_vv_wd_vv_variables(nfcntn)
  close(nfcntn)
end subroutine dump_mtd

subroutine init_mtd()
  use m_Const_Parameters, only : ON,old,check_file_name_on,formatted, CONTINUATION,INITIAL
  use m_Control_Parameters, only : printable, icond
  use m_variables_for_atoms, only : m_vfa_alloc_atoms
  use m_variables_for_dynamics, only : m_vfd_init
  use m_velocity_verlet, only : m_vv_init, m_vv_parse_input, m_vv_rd_vv_variables
  use m_meta_dynamics, only : m_mtd_parse_input, m_mtd_print_status, m_mtd_rd_mtd, m_mtd_gen_bias_mode &
  &                         , m_mtd_set_bias_gen_only, m_mtd_open_mtd_files,m_mtd_wd_bias_potential &
  &                         , m_mtd_finalize
  use m_Files, only : nfinp, F_INP, nfout, open0,nfcntn, F_ZAJ_in, F_CNTN_in, F_CNTN_BIN_in, F_CHGT_in &
  &                 , F_ZAJ, F_CNTN, F_CNTN_BIN, F_CHGT, m_Files_check_file_existence &
  &                 , file_existence_contfiles, file_existence_3contfiles, F_CNTN_bak
  use m_Parallelization, only : workdir,sw_wdir, mype_conf

  implicit none
  integer ::  f_selectTop, f_selectBlock, f_openInputFile, f_closeInputFile
  integer :: iiret
  logical :: ex,op

  if(printable) write(nfout,'(a)')
  if(printable) write(nfout,'(a)') 'start initialization of meta-dynamics'
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
    if(printable) write(nfout,'(" !!! Error in opening of inputfile: ",a32)') F_INP
    call phase_error_with_msg(nfout,'Error in opening of inputfile',__LINE__,__FILE__)
  else if(iiret > 0) then
    if(printable) write(nfout,'(" !!! There is something wrong in the input file: ",a32)') F_INP
  end if

  iiret = f_selectTop()

  call m_mtd_parse_input()

  if(m_mtd_gen_bias_mode())then
    call m_mtd_print_status()
    write(nfout,'(a)') ' -- generate and output the bias potential --'
    call m_mtd_wd_bias_potential()
    call m_mtd_finalize()
    write(nfout,'(a)') '... exiting.'
    stop
  else
    call m_vv_init()
    call m_vv_parse_input()

    iiret = f_selectTop()

    call m_mtd_print_status()

!    iiret = f_closeInputFile()

    if(printable) write(nfout,'(a)') 'end initialization of meta-dynamics'
    if(printable) write(nfout,'(a)')

    if(icond==CONTINUATION)then
      if(printable) write(nfout,'(a)') 'loading continuation data for the meta-dynamics...'
      call m_vv_rd_vv_variables(nfcntn,.true.)
      call m_mtd_rd_mtd(nfcntn,file_existence_contfiles.and.file_existence_3contfiles)
      if(printable) write(nfout,'(a)') '...done'
    endif
  endif

end subroutine init_mtd
#endif
! === "phase" repository merge: To make 3D_Parallel by tkato 2011/11/10 ========
! ==============================================================================

