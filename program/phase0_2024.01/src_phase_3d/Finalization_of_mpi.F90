!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 593 $)
!
!  SUBROUINE: Finalization_of_mpi
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!      Further modification by T. Yamasaki and M. Saito,   Feb-May. 2004
!
!  Contact address :  Phase System Consortium
!                     E-mail: phase_system@nims.go.jp URL https://azuma.nims.go.jp
!
!
!
!=======================================================================
subroutine Finalization_of_mpi
! $Id: Finalization_of_mpi.F90 593 2019-06-20 03:47:31Z jkoga $
  use m_Parallelization, only :       m_Parallel_end_mpi
  use m_PlaneWaveBasisSet, only :     m_pwBS_dealloc_ngpt_igfp_gr &
       &                            , m_pwBS_dealloc_ylm_l
#ifdef __EDA__
  use m_Ionic_System, only :          m_IS_dealloc_zfm3_3D &
       &                            , m_IS_dealloc_eewald_per_atom, m_IS_dealloc_zfm3_EDA
  use m_PseudoPotential, only :       m_PP_dealloc_psc_qitg_rhpcg, m_PP_dealloc_PP_per_atom_etc &
         &                            , m_PP_dealloc_paw,flg_paw
#else
  use m_Ionic_System, only :          m_IS_dealloc_zfm3_3D
  use m_PseudoPotential, only :       m_PP_dealloc_psc_qitg_rhpcg &
         &                            , m_PP_dealloc_paw,flg_paw
#endif
#ifndef PARAMSET
  use m_Charge_Density,       only : m_CD_dealloc_chgq
  use m_Electronic_Structure, only : m_ES_dealloc_vlhxc, m_ES_dealloc_vlhxcQ &
                                   , m_ES_dealloc_Dhub
  use m_ES_WF_by_SDorCG,      only : m_ESsd_dealloc_dzajn2, m_ESsd_dealloc_zaj_old
  use m_PAW_ChargeDensity,    only : m_PAWCD_dealloc
#ifdef __EDA__
  use m_XC_Potential_2D,      only : m_XC_dealloc_exc_on_a_grid
#endif
  use m_XC_Potential,         only : m_XC_dealloc_vxc_3D
#endif
#ifdef _POSITRON_
  use m_epc_potential,        only : m_epc_dealloc, m_epc_dealloc_vlhxc_p
  use m_Positron_Wave_Functions,only:m_pWF_deallocate_pzaj_etc
#endif

! ======== KT_add =============== 2013/10/31
  use m_Control_Parameters,   only : noncol, SpinOrbit_Mode
  use m_Const_Parameters,     only : Neglected, BuiltIn, ByPawPot, ZeffApprox, &
       &                             ByProjector, ReadFromPP, ON
  use m_SpinOrbit_Potential,  only : m_SO_dealloc_Dsoc, m_SO_dealloc_Mat_SOC_Strenth
! =============================== 2013/10/31


  use m_Control_Parameters,  only : sw_local_approx_trans_moment, sw_corelevel_spectrum
  use m_CLS_dipquad,         only : m_CLS_dealloc_dipquad, m_CLS_dealloc_transmom_ek
#ifdef KMATH_FFT3D
  use m_Control_Parameters,  only : sw_kmath_fft3d
  use m_Parallelization,     only : m_Parallel_kmath3d_finalize
#endif

#ifdef __EDA__
  use m_Control_Parameters,  only : sw_eda
#endif

  implicit none

  call m_pwBS_dealloc_ngpt_igfp_gr()
  call m_IS_dealloc_zfm3_3D()
  call m_PP_dealloc_psc_qitg_rhpcg()
  if(flg_paw) then
     call m_PP_dealloc_paw()
     call m_PAWCD_dealloc()
  end if
#ifndef PARAMSET
  call m_pwBS_dealloc_ylm_l()
  call m_CD_dealloc_chgq()
  call m_ESsd_dealloc_dzajn2()
  call m_ESsd_dealloc_zaj_old()
  call m_ES_dealloc_vlhxc()
  call m_ES_dealloc_vlhxcQ()
  call m_XC_dealloc_vxc_3D()
#ifdef __EDA__
  if(sw_eda==ON) then
! -----  ascat starts modifying  -----
  call m_XC_dealloc_exc_on_a_grid()
  call m_IS_dealloc_zfm3_EDA()
  call m_IS_dealloc_eewald_per_atom()
  call m_PP_dealloc_PP_per_atom_etc()
!  call m_CD_dealloc_chgq_EDA()
  endif
! -----  ascat ceases modifying  -----
#endif
#endif
#ifdef _POSITRON_
  call m_epc_dealloc()
  call m_epc_dealloc_vlhxc_p
  call m_pWF_deallocate_pzaj_etc()
#endif
  call m_ES_dealloc_Dhub()

! ==== KT_add ============== 2013/10/31
  if ( noncol ) then
     if ( SpinOrbit_Mode == ByPawPot .or. SpinOrbit_Mode == ZeffApprox &
          &                          .or. SpinOrbit_Mode == ReadFromPP ) then
        call m_SO_dealloc_Mat_SOC_strenth
        call m_SO_dealloc_Dsoc
     endif
     if ( SpinOrbit_Mode == ByProjector ) then
        call m_SO_dealloc_Dsoc
     endif
  endif
! ========================== 2013/10/31


  if ( sw_corelevel_spectrum == ON .and. sw_local_approx_trans_moment == ON ) then
     call m_CLS_dealloc_dipquad
     call m_CLS_dealloc_transmom_ek
  endif

#ifdef KMATH_FFT3D
  if(sw_kmath_fft3d==ON)then
    call m_Parallel_kmath3d_finalize()
  endif
#endif

  call m_Parallel_end_mpi           ! MPI
end subroutine Finalization_of_mpi

subroutine mpi_stop(nf)
  use mpi
  integer, intent(in) :: nf
!  include 'mpif.h'             ! MPI
  integer :: ierror
  integer :: errorcode

  errorcode = 1000
  write(nf,*) ' ','',''
  call flush(nf)
  call mpi_abort(mpi_comm_world,errorcode,ierror)
  call mpi_finalize(ierror)
  stop
end subroutine mpi_stop

!!$#ifdef DEBUG_ERRORS
subroutine phase_error(ierrNO, nfout, nf, filename, line, modulefile)
!!$#else
!!$subroutine phase_error(ierrNO, nfout, nf, filename)
!!$#endif
  use m_Control_Parameters, only : ipri
  use m_ErrorMessages
  use m_Parallelization,    only : MPI_CommGroup
  implicit none
  integer, intent(in) :: ierrNO,nfout,nf
  character(len=*), intent(in)          :: filename
!!$#ifdef DEBUG_ERRORS
  integer, intent(in),optional            :: line
  character(len=*), intent(in),optional   :: modulefile
!!$#endif

!!$  include 'mpif.h'             ! MPI

  character(len=255) :: I_name
  integer :: ierror, len_filename, i

  ierror = 0
  if(ipri>=1) then
     len_filename = min(255,len(trim(filename))) !ASMS
     do i = 1, len_filename
        I_name(i:i) = filename(i:i)
     end do
     do i = len_filename + 1 , 255 !ASMS
        I_name(i:i) = " "          !ASMS
     end do                        !ASMS

     if(ierrNO == F_POT_FILE_NOT_EXIST) then
        write(nfout,'("#### ERROR(",i4,"): No F_POT(",i3,") (=",a,")")') &
             & ierrNO, nf, trim(filename) !ASMS
     else if(ierrNO == FILE_NOT_EXIST) then
        write(nfout,'("#### ERROR(", i4, "): File (",a &
             & ,") assigned to a file no",i4," does not exist.")') &
             & ierrNO, trim(filename),nf !ASMS
     else if(ierrNO == F_ZAJ_FILE_NOT_EXIST) then    ! 1204
        write(nfout,'("#### ERROR(", i4, "):",a)') &
             & ierrNO, msg_1204_error
     else if(ierrNO == F_CHGT_FILE_NOT_EXIST) then   ! 1205
        write(nfout,'("#### ERROR(", i4, "):",a)') &
             & ierrNO, msg_1205_error
     else if(ierrNO == ERROR_IN_INPUTFILE_OPENING) then ! 1211
        write(nfout,'("#### ERROR(", i4, "):",a,a)') &
             & ierrNO, msg_1211_error, trim(filename) !ASMS
     else if(ierrNO == F_CHGT_FILE_NOT_EXIST_EK) then   ! 1212
        write(nfout,'("#### ERROR(",i4,"): F_CHGT(file no",i3,") does not exist. ",a)') &
             & ierrNO, nf, msg_1212_error
     end if
  end if

#ifdef DEBUG_ERRORS
  if(present(line)) then
     if(present(modulefile)) then
        write(nfout,'(a,i8,a)') ' PHASE experienced an error at the line ' &
             & , line, ' of source file '//trim(adjustl(modulefile))
     end if
  end if
#endif
  call flush(nfout)
!  stop ' after flush'
  call mpi_abort(MPI_CommGroup,261,ierror)
  call mpi_finalize(ierror)
end subroutine phase_error

subroutine phase_error_wo_filename(ierrNO, nfout, nf, line, modulefile)

  use m_Control_Parameters, only : ipri
  use m_Parallelization,    only : MPI_CommGroup
  use m_ErrorMessages
  use mpi

  integer, intent(in) :: ierrNO,nfout
  integer, intent(in),optional :: nf
  integer, intent(in),optional :: line
  character(len=*), intent(in),optional ::modulefile

!  include 'mpif.h'             ! MPI

  logical :: I_opened
  character(len=255) :: I_name
  integer :: ierror

  ierror = 0

  if(ipri>=1) then
!!$     write(nfout,'("#### ERROR(",i4,")")') ierrNo
     if(ierrNO == EOF_REACHED) then
        if(present(nf)) then
           inquire(unit=nf, OPENED=I_opened, NAME=I_name)
           if(I_opened) then
              write(nfout,'("#### ERROR(",i4,"): EOF_reached. File No = " &
                   & ,i4,", File Name = ",a)') ierrNO, nf, trim(adjustl(I_name))
           else
              write(nfout,'("#### ERROR(",i4,"): EOF_reached. File No = " &
                   & ,i4)') ierrNO, nf
           end if
        else
           write(nfout,'("#### ERROR(",i4,"): EOF_reached. But, File No is not given.")') &
                & ierrNO
        end if
     else if(ierrNO == FORMAT_ERROR) then
        if(present(nf)) then
           inquire(unit=nf, OPENED=I_opened, NAME=I_name)
           if(I_opened) then
              write(nfout,'("#### ERROR(",i4,"): ",a," File No = ",i4," File Name = ",a)') &
                   &                     ierrNO, msg_format_error, nf, trim(adjustl(I_name))
           else
              write(nfout,'("#### ERROR(",i4,"): ",a," File No = ",i4)') &
                   &                     ierrNO, msg_format_error, nf
           end if
        else
           write(nfout,'("#### ERROR(",i4,"): ",a," But, File No is unknown.")') &
                &                     ierrNO, msg_format_error
        end if
     else if(ierrNO == CPP_DEFINE_ERROR_1) then
        write(nfout,'("#### ERROR(",i4,"): ",a)') ierrNO,msg_6100_error
        write(nfout,'("                      (",a,")")') msg_6101_error
     else
        write(nfout,'("#### ERROR(",i4,"):msg")') ierrNO
     end if
#ifdef DEBUG_ERRORS
     if(present(line)) then
        if(present(modulefile)) then
           write(nfout,'(a,i8,a)') ' PHASE aborted at the line ' &
                & , line, ' of source file '//trim(adjustl(modulefile))
        end if
     end if
#endif
  end if
  call flush(nfout)
!!$  if(mode>0) then
     call mpi_abort(MPI_CommGroup,261,ierror)
     write(nfout,'(" ierror = ",i8)') ierror
!!$  end if
  call flush(nfout)
  call mpi_finalize(ierror)
  stop
end subroutine phase_error_wo_filename

subroutine phase_error_with_msg(nfout, msg, line, modulefile)
!!$#else
!!$subroutine phase_error(ierrNO, nfout, nf, filename)
!!$#endif
  use m_Control_Parameters, only : ipri, printable, iprijobstatus, jobstatus_format, jobstatus_series
  use m_Files,              only : m_Files_open_nfstatus, m_Files_skiptoend, nfstatus &
  &                              , m_Files_close_nfstatus
  use m_ErrorMessages
  use m_Parallelization,    only : MPI_CommGroup, mype
  use m_Files, only : m_Files_close_all, m_Files_close_logfile
  use m_IterationNumbers,   only : iteration, iteration_ionic, iteration_electronic
  use m_Const_Parameters,   only : ON, ERROR
  use m_Timing,             only : m_Timing_wd_status

  implicit none

  integer,          intent(in)   :: nfout
  character(len=*), intent(in)   :: msg
  integer, intent(in)            :: line
  character(len=*), intent(in)   :: modulefile

!!$  include 'mpif.h'             ! MPI

  character(len=255) :: I_name
  integer :: ierror, len_filename, i
  integer :: status_wdmode

  ierror = 0
  call mpi_barrier(MPI_CommGroup,ierror)
  if(mype==0) then
    write(nfout,'(a,i8,a)') ' PHASE/0 terminated abnormally at line ', line, ' of file '//trim(adjustl(modulefile))
    write(nfout,'(a)')      ' reason : '//trim(adjustl(msg))

    write(nfout,'(a)') '----------------------------------------'
    write(nfout,'(a)') ' EEEEEE  RRRRR   RRRRR    OOOO   RRRRR  '
    write(nfout,'(a)') ' E       R    R  R    R  O    O  R    R '
    write(nfout,'(a)') ' EEEEE   R    R  R    R  O    O  R    R '
    write(nfout,'(a)') ' E       RRRRR   RRRRR   O    O  RRRRR  '
    write(nfout,'(a)') ' E       R   R   R   R   O    O  R   R  '
    write(nfout,'(a)') ' EEEEEE  R    R  R    R   OOOO   R    R '
    write(nfout,'(a)') '----------------------------------------'

    write(0,'(a,i8,a)') ' PHASE/0 terminated abnormally at line ', line, ' of file '//trim(adjustl(modulefile))
    write(0,'(a)')      ' reason : '//trim(adjustl(msg))
    call flush(nfout)
    call flush(0)
  endif

  if(iprijobstatus >= 1) then
!!$     write(nfout,'(" WDD onto Files")')
     call m_Files_open_nfstatus()
     if(jobstatus_series == ON) then
!!$        switch_header = OFF
        call m_Files_skiptoend(nfstatus)
     else
!!$        switch_header = ON
     end if
     status_wdmode = ERROR
     call m_Timing_wd_status(nfstatus,jobstatus_format,jobstatus_series,status_wdmode &
          &                , iteration,iteration_ionic,iteration_electronic)
  end if
  call m_Files_close_all()
  call m_Files_close_logfile()
  call mpi_abort(MPI_CommGroup,261,ierror)
  call mpi_finalize(ierror)
end subroutine phase_error_with_msg

subroutine phase_execution_error(ierrNO)
  use m_Control_Parameters, only : ipri
  use m_Parallelization,    only : MPI_CommGroup, mype, ierr
  use m_Files,              only : nfout
  use m_ErrorMessages
  use m_Const_Parameters,   only : FMAXTAGLEN

  implicit none
  integer f_check_input_duplicated_tag
  integer, intent(in) :: ierrNO
  character(len=FMAXTAGLEN) :: duplicated_block_name
  integer :: len_name,i
!!$  include 'mpif.h'             ! MPI

  if(ipri>=1 .and. mype == 0) then
     write(nfout,*)  ! BLANK LINE
     select case (ierrNO)
     case (CONT_FILES_NOT_EXIST)      ! 1209
        write(nfout,'("###ERROR(",i4,") file(s) for continuation do not exist.")') &
             & CONT_FILES_NOT_EXIST
     case (FILENAMES_FORMAT_ERROR)      ! 1209
        write(nfout,'("###ERROR(",i4,") format error in file_names.data.")') &
             & FILENAMES_FORMAT_ERROR
     case (FILENAMES_FORMAT_ERROR_NEB)  ! 1210
        write(nfout,'("###ERROR(",i4,") format error in file_names.data. (nebfiles)")') &
             & FILENAMES_FORMAT_ERROR_NEB
     case (ERROR_IN_INPUTFILE_OPENING)  ! 1211
        write(nfout,'("###ERROR(",i4,") syntax error in F_INP.")') &
             & ERROR_IN_INPUTFILE_OPENING
        if(f_check_input_duplicated_tag().ne.0) then              ! (input_interface.F90)
           call f_cp_duplicated_block_name(duplicated_block_name) ! (input_interface.F90)
           write(nfout,'(" #-- Syntax Error: Block ",a1,a,a1," is duplicately defined in F_INP.")') &
                & char(34),trim(duplicated_block_name),char(34)
        end if
     case (INVALID_CHARGE_MIXING)       ! 1302
        write(nfout,'("###ERROR(",i4,") charge-mixing in the wavefunction_solver block does not exist.")') &
             & INVALID_CHARGE_MIXING
     case (INVALID_ATOMIC_NUMBER)       ! 1401
        write(nfout,'("###ERROR(",i4,") The atomic-number is inconsistent with the pseudopotential file.")') &
             & INVALID_ATOMIC_NUMBER
     case (PARALLELIZATION_INVALID_2D)  ! 2101
        write(nfout,'("###ERROR(",i4,") Number of parallel process (npes) must be ne*nk.")')  &
             & PARALLELIZATION_INVALID_2D
     case (PARALLELIZATION_INVALID_3D)  ! 2102
        continue
     case (PARALLELIZATION_INVALID_NK)  ! 2103
        write(nfout,'("###ERROR(",i4,") Number of kpoint-parallelization is greater than number of kpoints(kv3).")') &
             & PARALLELIZATION_INVALID_NK
     case (PARALLELIZATION_INVALID_NE)  ! 2104
        write(nfout,'("###ERROR(",i4,") Number of band-parallelization is greater than number of bands(neg).")') &
             & PARALLELIZATION_INVALID_NE
     end select
  end if
  call flush(nfout)
  call mpi_barrier(MPI_CommGroup,ierr)
  call mpi_abort(MPI_CommGroup,275,ierrNO)
  call mpi_finalize(ierrNO)
end subroutine phase_execution_error
