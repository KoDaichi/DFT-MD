!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 590 $)
!
!  MODULE: m_IterationNumbers
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
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
! *************************************************************
!
module m_IterationNumbers
!     (m_Iter)
! $Id: m_IterationNumbers.F90 590 2018-11-06 04:20:58Z jkoga $
  use m_Const_Parameters, only : len_tag_iteration, tag_iteration &
       &                       , len_tag_iters_and_nk, tag_iters_and_nk &
       &                       , EK_CONVERGED, COORDINATE_CONTINUATION
  use m_Parallelization, only  : MPI_CommGroup,mype,npes,ierr,nrank_k
  use mpi
  implicit none

  integer :: iteration = 0            ! #total iteration
  integer :: first_iteration_of_this_job = 0
  integer :: iteration_ionic = 1      ! #timestep for MD of ions
  integer :: iteration_scdft = 1
  integer :: iteration_scdft_initial = 1
  logical :: iteration_scdft_initial_is_set = .false.
  integer :: iteration_electronic = 0 ! #iteration for electronic states convergence
  integer :: iteration_rmm_start  = 0 ! #iteration where rmm starts in iteration_electonic
  integer :: iteration_for_cmix   = 0 ! #iteration for charge mixing.
  !      This is almost same with iteration_electronic.

  integer :: iteration_unit_cell  = 1 ! #iteration for the unit cell, when unit cell optimization is enabled

  integer :: iteration_uramp = 1 ! #iteration for the U-ramping method, when the U-ramping method is enabled

  integer :: iteration_dimer = 1 ! #iteration for the dimer method

  !      Unique difference is that this value is reset when the job starts to run.
  integer :: nk_in_the_process    = 0 ! #kpoint calculating. This is effective only
  !                                     when ekmode (in m_Control_Parameters) is ON.
  integer :: nk_converged         = 0
  integer :: nkgroup              = 0
  integer :: first_kpoint_in_this_job = 0
  integer :: first_iteration_electronic = 0

  integer :: iteration_positron_wf = 0

  integer :: iteration_stress_correction = 1

  integer, parameter :: len_str = 132
  character(len=len_str) :: str

!  include 'mpif.h'
contains

  subroutine m_Iter_reset_iter_ionic()
     iteration_ionic = 1
  end subroutine m_Iter_reset_iter_ionic

  subroutine m_Iter_reset_iter_scdft()
    iteration_scdft = 1
  end subroutine m_Iter_reset_iter_scdft

  subroutine m_Iter_reset_iter_electronic(nfout)
    integer, intent(in) :: nfout
    iteration_electronic = 0
    first_iteration_electronic = iteration_electronic
    if(mype == 0) then
       write(nfout,'(" iteration_electronic is rest 0")')
    end if
  end subroutine m_Iter_reset_iter_electronic

  subroutine m_Iter_wd_iteration_numbers(nfcntn)
    integer, intent(in) :: nfcntn
    logical :: unitcell_can_change,Uramping,isSCDFT
    if(mype==0) then
       write(nfcntn,*) tag_iteration
       if(.not.unitcell_can_change())then
          if (Uramping()) then
            write(nfcntn,'(4i10)') iteration, iteration_uramp,iteration_ionic, iteration_electronic
          else if (isSCDFT()) then
            write(nfcntn,'(4i10)') iteration, iteration_scdft,iteration_ionic, iteration_electronic
          else
            write(nfcntn,'(3i10)') iteration, iteration_ionic, iteration_electronic
          endif
       else
          write(nfcntn,'(4i10)') iteration, iteration_unit_cell, iteration_ionic, iteration_electronic
       endif
    endif
  end subroutine m_Iter_wd_iteration_numbers

  subroutine m_Iter_wd_iteration_numbers_b(nf,F_partitioned)
    integer, intent(in) :: nf
    logical, intent(in) :: F_partitioned
    call mpi_barrier(MPI_CommGroup,ierr)
    if(F_partitioned) then
       if(mype == 0) then
          write(nf,*) iteration, iteration_ionic, iteration_electronic
       end if
    else
       write(nf,*) iteration, iteration_ionic, iteration_electronic
    end if
  end subroutine m_Iter_wd_iteration_numbers_b

  subroutine m_Iter_rd_iteration_numbers(nfcntn,icond)
    integer, intent(in) :: nfcntn
    integer, intent(in) :: icond
    logical             :: EOF_reach, tag_is_found
    logical :: unitcell_can_change,Uramping,isSCDFT

    if(mype==0) then
       call rewind_to_tag0(nfcntn,len_tag_iteration,tag_iteration &
            &, EOF_reach, tag_is_found, str,len_str)
       if(.not.tag_is_found) then
          call phase_error_with_msg(6,' tag_iteration is not found',__LINE__,__FILE__)
       else
          if(.not.unitcell_can_change())then
             if(Uramping()) then
               read(nfcntn,*) iteration, iteration_uramp,iteration_ionic, iteration_electronic
             else if (isSCDFT())then
               read(nfcntn,*) iteration, iteration_scdft,iteration_ionic, iteration_electronic
             else
               read(nfcntn,*) iteration, iteration_ionic, iteration_electronic
             endif
          else
             read(nfcntn,*) iteration, iteration_unit_cell, iteration_ionic, iteration_electronic
          endif
!          if(icond==COORDINATE_CONTINUATION.and.iteration_electronic>0) then
!             iteration=iteration-iteration_electronic
!             iteration_electronic=0
!          endif
       endif
    endif
    if(npes > 1) then
       call mpi_bcast(iteration,1 &
            & ,mpi_integer,0,MPI_CommGroup,ierr)
       call mpi_bcast(iteration_ionic,1 &
            & ,mpi_integer,0,MPI_CommGroup,ierr)
       call mpi_bcast(iteration_electronic,1 &
            & ,mpi_integer,0,MPI_CommGroup,ierr)
       if(unitcell_can_change()) call mpi_bcast(iteration_unit_cell,1 &
            & ,mpi_integer,0,MPI_CommGroup,ierr)
       if(Uramping()) call mpi_bcast(iteration_uramp,1 &
            & ,mpi_integer,0,MPI_CommGroup,ierr)
       if(isSCDFT()) call mpi_bcast(iteration_scdft,1 &
            & ,mpi_integer,0,MPI_CommGroup,ierr)
    end if
    first_iteration_of_this_job = iteration
!!$    iteration_ionic =  iteration_ionic + 1
  end subroutine m_Iter_rd_iteration_numbers

!!$  subroutine m_Iter_rd_iters_and_nk(icon,nfcntn,nspin)
  subroutine m_Iter_rd_iters_and_nk(nfcntn,nspin,nfout)
    integer, intent(in) :: nfcntn, nspin, nfout
    logical             :: EOF_reach, tag_is_found

!!$    write(nfout,'(" !! iconvergence = ",i6)') icon
!!$    if(icon < EK_CONVERGED) then
       if(mype==0) then
          call rewind_to_tag0(nfcntn,len_tag_iters_and_nk,tag_iters_and_nk &
               &, EOF_reach, tag_is_found, str,len_str)
          if(.not.tag_is_found) then
             call phase_error_with_msg(nfout,' tag_iters_and_nk is not found',__LINE__,__FILE__)
          else
             read(nfcntn,*) iteration, iteration_electronic, nk_in_the_process
          endif
       endif
       if(npes > 1) then
          call mpi_bcast(iteration,           1,mpi_integer,0,MPI_CommGroup,ierr)
          call mpi_bcast(iteration_electronic,1,mpi_integer,0,MPI_CommGroup,ierr)
          call mpi_bcast(nk_in_the_process,   1,mpi_integer,0,MPI_CommGroup,ierr)
       end if
!!$    else
!!$       iteration = 1; iteration_electronic = 1; nk_in_the_process = 1
!!$    end if
    if(mype==0) then
         write(nfout,'(" !! iteration, iteration_electronic, nk_in_the_process = ",3i6)')&
         & iteration, iteration_electronic, nk_in_the_process
    end if
    first_iteration_of_this_job = iteration
    first_iteration_electronic  = iteration_electronic
    first_kpoint_in_this_job    = nk_in_the_process
    nk_converged                = first_kpoint_in_this_job - 1
    iteration_electronic = iteration_electronic - 1
    nk_in_the_process    = nk_in_the_process - nspin*nrank_k
!!$    iteration_ionic =  iteration_ionic + 1
  end subroutine m_Iter_rd_iters_and_nk

  subroutine m_Iter_wd_iters_and_nk(nfcntn)
    integer, intent(in) :: nfcntn
    if(mype==0) then
       write(nfcntn,*) tag_iters_and_nk
       write(nfcntn,'(3i10)') iteration, iteration_electronic, nk_in_the_process
    end if
  end subroutine m_Iter_wd_iters_and_nk

  subroutine m_Iter_set_converged_nk(nspin)
    integer, intent(in) :: nspin
    nk_converged = nk_in_the_process + nspin-1
  end subroutine m_Iter_set_converged_nk

  subroutine m_Iter_nk_incre(incre)
    integer, intent(in) :: incre
    integer :: aincre
    if(incre == 2 .and. nk_in_the_process == 0) then
       aincre = 1
    else
       aincre = incre
    end if
    nk_in_the_process = nk_in_the_process + aincre
  end subroutine m_Iter_nk_incre

  subroutine m_Iter_nk_incre2(incre,kv3_ek)
    integer, intent(in) :: incre,kv3_ek
    integer :: aincre
    if(nk_in_the_process == 0) then
       aincre = 1
    else
       aincre = nrank_k*incre
    end if
    nk_in_the_process = nk_in_the_process + aincre
    if(nk_in_the_process > kv3_ek+1) nk_in_the_process = kv3_ek + 1
    nkgroup = nkgroup+1
  end subroutine m_Iter_nk_incre2

  subroutine m_Iter_nkgroup_incre()
    nkgroup=nkgroup+1
  end subroutine m_Iter_nkgroup_incre

  subroutine m_Iter_wd_nk(nfout)
    integer, intent(in) :: nfout
    if(mype==0) then
       write(nfout,*)
       write(nfout,'(" ------------------------")')
       write(nfout,'(" nk_in_the_process = ",i6)') nk_in_the_process
    end if
  end subroutine m_Iter_wd_nk

  subroutine m_Iter_wd_nk2(nfout,kv3)
    integer, intent(in) :: nfout,kv3
    if(mype==0) then
       write(nfout,*)
       write(nfout,'(" ------------------------")')
       write(nfout,'(" nk_in_the_process = ",i6, " nk in the k_point_group = ",i6)') &
            & nk_in_the_process, nk_in_the_process+kv3-1
       write(nfout,'(" nkgroup = ",i6)') nkgroup
    end if
  end subroutine m_Iter_wd_nk2

  subroutine m_Iter_mdIterN_increment(nfout)
    integer, intent(in) :: nfout
    iteration_ionic   = iteration_ionic + 1
    if(mype==0) then
       write(nfout,*) ' ! iteration_ionic = ', iteration_ionic
    end if
  end subroutine m_Iter_mdIterN_increment

  subroutine m_Iter_dimer_incre(nfout)
    integer, intent(in) :: nfout
    iteration_dimer  = iteration_dimer + 1
    if(mype==0) then
       write(nfout,*) ' ! iteration_dimer = ', iteration_dimer
    end if
  end subroutine m_Iter_dimer_incre

  subroutine m_Iter_scdftIterN_increment(nfout)
    integer, intent(in) :: nfout
    iteration_scdft  = iteration_scdft + 1
    if(mype==0) then
       write(nfout,*) ' ! iteration_scdft = ', iteration_scdft
    end if
  end subroutine m_Iter_scdftIterN_increment

  subroutine m_Iter_unitcell_increment()
    iteration_unit_cell = iteration_unit_cell + 1
  end subroutine m_Iter_unitcell_increment

  subroutine m_Iter_stress_correction_incre()
    iteration_stress_correction = iteration_stress_correction + 1
  end subroutine m_Iter_stress_correction_incre

  subroutine m_Iter_uramp_increment()
    iteration_uramp = iteration_uramp + 1
  end subroutine m_Iter_uramp_increment

  subroutine m_Iter_electronic_reset()
    iteration_electronic = 0
    iteration_for_cmix   = 0
  end subroutine m_Iter_electronic_reset

  subroutine m_Iter_cmix_reset()
    iteration_for_cmix  = 0
  end subroutine m_Iter_cmix_reset

  subroutine m_Iter_electronic_set()
    iteration_electronic = first_iteration_electronic
  end subroutine m_Iter_electronic_set

!!$  subroutine m_Iter_nkgroup_set(kv3)
!!$    integer, intent(in) :: kv3
!!$    nkgroup = (nk_in_the_process-1)/kv3
!!$  end subroutine m_Iter_nkgroup_set

  subroutine m_Iter_wd_electronic(nfout)
    integer, intent(in) :: nfout
    if(mype==0) then
       write(nfout,'(" iteration_electronic = ",i6)') iteration_electronic
    end if
  end subroutine m_Iter_wd_electronic

  subroutine m_Iter_initialize
    iteration = 0
  end subroutine m_Iter_initialize

  subroutine m_Iter_electronic_incre
    iteration_electronic = iteration_electronic + 1
    iteration_for_cmix   = iteration_for_cmix   + 1
  end subroutine m_Iter_electronic_incre

  subroutine m_Iter_total_increment
    iteration = iteration + 1
  end subroutine m_Iter_total_increment

  subroutine m_Iter_set_rmm_start
    iteration_rmm_start = iteration_electronic
  end subroutine m_Iter_set_rmm_start

  subroutine m_Iter_positron_incre()
    iteration_positron_wf = iteration_positron_wf + 1
  end subroutine m_Iter_positron_incre

  subroutine m_Iter_positron_set()
    iteration_positron_wf = 0
  end subroutine m_Iter_positron_set

  subroutine m_Iter_reset_nk_in_the_process()
    nk_in_the_process = 0
    nkgroup = 0
  end subroutine m_Iter_reset_nk_in_the_process

end module m_IterationNumbers
