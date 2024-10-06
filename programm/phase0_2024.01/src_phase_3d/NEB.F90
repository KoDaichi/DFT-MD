!
! =========== Contributions ===================================
!
! Through the courtesy of contributors, the following functions are added.
!
! Company:  ASMS Co.,Ltd.
! Functions:  [Identifier: 13.0D]
!                     NEB works with PAW
!
! =============================================================

subroutine do_neb()
  use m_Files,                only : nfneb,nfout
  use m_Const_Parameters, only : OFF,ON,PREPARATION_ONLY,RETURN_AFTER_SYMMCHECK,SKIP_SYMMCHECK
  use m_IterationNumbers, only : iteration,iteration_electronic,iteration_ionic
  use m_Replica
  use m_Parallelization, only : MPI_CommGroup, mype, mype_conf, m_Parallel_end_mpi &
                              , m_Parallel_init_mpi_paw_3D

! ======================== KT_add ================= 13.0D
  use m_Files,         only :  m_Files_open_nfcntn_bin_paw, nfout, nfcntn_bin_paw &
                            ,  m_Files_checkpoint_dir_neb
  use m_PseudoPotential,  only : flg_paw, m_PP_rd_PAW_parameters
  use m_Charge_Density,    only :  m_CD_alloc_chgsoft
! ================================================= 13.0D

  use m_Control_Parameters, only : sw_fcp, cpt_neb_iteration, icond, sw_hubbard &
                                 , ipriparallel, PAW_switch
  use m_Fcp, only :  m_Fcp_set_neb_replica_totch, &
     m_Fcp_set_neb_replica_force, m_Fcp_update_neb, m_Fcp_set_neb_parameter, &
     m_Fcp_allreduce_neb_force, m_Fcp_initialize_neb_totch
  use m_Ionic_System, only : natm
  use m_PseudoPotential, only : mmesh
  use mpi

  implicit none
  logical  :: ChargeDensity_is_Converged, TotalEnergy_is_Divergent
  logical  :: Already_Converged, Positron_bulk
  logical  :: Hubbard_model
  logical  :: Forces_are_Converged, Ending_Time, Ending_Time2, Force_errors_are_tolerable
  logical  :: replica_converged

  integer i, j, itr, ierr
  integer mpi_err
  logical all_conv_flag

  integer prepare_communicator
  logical :: pp_generated

  logical :: dump_checkpoint
  logical :: to_dimer = .false.

!  include 'mpif.h'

  call initialize_neb

  !call Initialization(OFF)
  !call InputData_Analysis(OFF) ! rd_cntn_data = ON

  call set_neb_condition
  call create_replica

  itr = 0
  all_conv_flag = all(neb%image(:)%scf_convergence)  ! for restart
  prepare_communicator = ON
  pp_generated=.false.

  do i=1,neb%number_of_images
     call set_neb_replica_parameter(i)
     call set_neb_replica_coordinate(i)
     if(.not.(i==1 .or. i==neb%number_of_images)) call InputData_Analysis_neb()
     call Preparation(RETURN_AFTER_SYMMCHECK)
     call set_sym(i)
     call deallocate_array_partial()
     call remove_annoying_files(nfout)
  enddo
  call m_Replica_resolve_lowest_symmetry(nfneb)
  if(icond == PREPARATION_ONLY)then
     call mpi_barrier(MPI_CommGroup,ierr)
     call m_Parallel_end_mpi()
     stop 'The preparation has been done.'
  endif
  NEBIteration: do

     itr = itr + 1
     if(all(neb%image(:)%scf_convergence)) then
        neb%step = neb%step + 1
     end if
     call set_neb_parameter(itr)
     if( sw_fcp == ON ) call m_Fcp_set_neb_parameter(itr)

     MultiReplica: do i = 1, neb%number_of_images
        if(mod(i-1,nrank_r) /= mype_r) then
           neb%image(i)%scf_convergence = .true.
           cycle
        end if
        call set_lowest_symmetry_op_tau()
        if(neb%cond%condition == 1) then
           if(itr == 1 .and. (.not. all_conv_flag)) then
              if(neb%image(i)%scf_convergence) then
                 !neb%image(i)%energy = neb%image(i)%energy0
                 !neb%image(i)%force = neb%image(i)%force0
                 !neb%image(i)%force_org = neb%image(i)%force0
                 cycle
              end if
           end if
        end if

        if(neb%cond%end0_energy_given .and. i==1)then
            neb%image(1)%energy = neb%image(1)%energy0
            cycle
        endif
        if(neb%cond%end1_energy_given .and. i==neb%number_of_images)then
            neb%image(neb%number_of_images)%energy = neb%image(neb%number_of_images)%energy0
            cycle
        endif
        if(itr > 1 .and. (i==1 .or. i==neb%number_of_images)) then
           neb%image(i)%energy = neb%image(i)%energy0
           neb%image(i)%force = neb%image(i)%force0
           neb%image(i)%force_org = neb%image(i)%force0
           cycle
        end if

!        write(nfneb,*) 'neb: ', neb%step, i, mype
        if (mype==0) write(nfneb,*) 'neb: ', neb%step, i

        call set_neb_replica_parameter(i)
        call set_neb_replica_coordinate(i)

        if(.not.neb%image(i)%scf_convergence) then
           call InputData_Analysis_neb()
        end if

        call Preparation(SKIP_SYMMCHECK)                   ! Basis set, symmetry check etc.
        call Preparation_for_mpi(prepare_communicator)           ! mpi

! ======================== KT_mod ================= 13.0D
!        if (.not.pp_generated) then
!          call PseudoPotential_Construction
!          pp_generated = .true.
!        endif

        if (.not.pp_generated) then
           call PseudoPotential_Construction
           pp_generated = .true.
#ifdef ENABLE_ESM_PACK
           call Preparation_for_ESM()
#endif
        else
           if ( flg_paw ) then
              if ( itr > 1 .or. neb%cond%condition == 1 ) then
                 call m_PP_rd_PAW_parameters(nfout,nfcntn_bin_paw)
              endif
              call m_CD_alloc_chgsoft
           endif
        endif
        if(PAW_switch==ON) call m_Parallel_init_mpi_paw_3D(nfout,ipriparallel,natm,mmesh)
! ================================================= 13.0D

        call Ewald_and_Structure_Factor
        call Initial_Electronic_Structure
        call Initial_MD_Condition

        if( sw_fcp == ON ) then
           call m_Fcp_initialize_neb_totch()
           call m_Fcp_set_neb_replica_totch(i)
        end if

        if(neb%cond%condition == 1 .and. neb%image(i)%scf_convergence) then
           iteration = 0
           iteration_electronic = 0
           iteration_ionic = 0
        end if
        neb%image(i)%scf_convergence = .true.

        if(.not.Already_Converged()) then
           AtomicConfiguration: do
              ChargeDensity:    do
                 call IterationNumber_Setting
                 call Renewal_of_WaveFunctions
                 call ChargeDensity_Construction(1)
                 call ChargeDensity_Mixing
                 call Renewal_of_Potential
                 if(Hubbard_model()) then
                    call Renewal_of_Hubbard_Potential
                 end if
!                 if(Ending_Time())                 exit ChargeDensity
                 if(Ending_Time()) then
                    neb%image(i)%scf_convergence = .false.
                    exit AtomicConfiguration
                 end if
                 if(TotalEnergy_is_Divergent())    exit ChargeDensity
                 if(ChargeDensity_is_Converged())  exit ChargeDensity
              enddo ChargeDensity
#ifdef LIBRARY_BUILD
              call Forces_phase0
#else
              call Forces
#endif
!              if(Forces_are_Converged()) exit AtomicConfiguration
              exit AtomicConfiguration
!!$              if(Force_errors_are_tolerable()) then
!!$                 call Move_Ions
!!$                 call Ewald_and_Structure_Factor
!!$                 call MDIterationNumber_Setting
!!$              end if

           enddo AtomicConfiguration

        end if

!        if(.not.Ending_Time2()) call set_neb_replica_energy_force(i)
        call set_neb_replica_energy_force(i)
        if( sw_fcp == ON ) call m_Fcp_set_neb_replica_force(i)

!!$        if(Positron_bulk()) then
!!$           call Renewal_of_pPotential()
!!$           call Solve_pWaveFunctions()
!!$        end if

!!$        call Postprocessing
        call WriteDownData_onto_Files(.true.)

        call deallocate_array

        if (mype==0) write(nfneb,*) 'scf convergence: ', i, neb%image(i)%scf_convergence

!!        prepare_communicator=OFF
        first_replica_done = .true.
        call flush(nfneb)
     end do MultiReplica

     call allreduce_scf_convergence
     call allreduce_neb_energy_force
     if( sw_fcp == ON ) call m_Fcp_allreduce_neb_force(nrank_r)

     if(.not. all(neb%image(:)%scf_convergence)) then
        if( mype==0) write(nfneb,*) 'scf not convergence...'
        exit NEBIteration
     end if

!     call allreduce_neb_energy_force

     call force_neb
     call atom_update_neb
     if( sw_fcp == ON ) call m_Fcp_update_neb(neb%dt, neb%cond%time_integral)
!     call atom_update_pos_neb

     if(replica_converged()) then
        exit NEBIteration
     end if
     if(neb%max_iteration>0 .and. neb%step >= neb%max_iteration) then
        exit NEBIteration
     end if

     if(dump_checkpoint(itr)) then
        call m_Files_checkpoint_dir_neb()
        call write_down_data_replica(.false.)
     endif

     if(curr_thres < neb%cond%dimer_thres) then
       to_dimer = .true.
       exit NEBIteration
     endif

  end do NEBIteration

  call m_Files_checkpoint_dir_neb()
  call write_down_data_replica(.not. to_dimer)

  if(.not. to_dimer) then
    call finalize_neb()
    call mpi_finalize(mpi_err)
  else
    call transition_to_dimer_method()
  endif
end subroutine do_neb

subroutine transition_to_dimer_method()
  use m_Ionic_System, only : cps_end0, cps_end1
  use m_Parallelization, only : mype
  use m_Files, only : nfneb
  use m_dimer, only : m_dm_set_pp_generated, m_dm_set_standard_out_opened, m_dm_do_dimer_method
  use m_Control_Parameters, only : icond
  use m_Replica
  implicit none
  integer :: k, targ, mpi_err
  real(kind=8) :: max_energy
  if(mype==0) write(nfneb, *) 'transitioning to the dimer method'
  max_energy = maxval( neb%image(:)%energy )
  targ = -1
  do k=2,neb%number_of_images-1
     if (neb%image(k)%energy==max_energy) then
        targ = k
        exit
     endif
  enddo
  if (targ<0) then
    if(mype==0) write(nfneb,*) 'failed to find the highest-energy replica'
    call finalize_neb()
    call mpi_finalize(mpi_err)
  else
    if(mype==0) then
      write(nfneb,*) 'create initial dimer from replica ',targ-1,targ+1
    endif
    cps_end0 = neb%image(targ-1)%cps
    cps_end1 = neb%image(targ+1)%cps
    icond = 0
!    call m_Files_open_standardout()
    call finalize_neb()
    call m_dm_set_pp_generated(.true.)
    call m_dm_set_standard_out_opened(.false.)
    call m_dm_do_dimer_method()
  endif

end subroutine transition_to_dimer_method

subroutine remove_annoying_files(nf)
  use m_Parallelization, only : mype
  use m_Control_Parameters,  only : ipriparadeb
  implicit none
  integer, intent(in) :: nf
#ifdef NEB_NEW_FILENAMES
  return
#endif
  if ( mype == 0 ) then
     close(nf)
  else if(ipriparadeb == 0) then
     close(nf, status='delete')
  else
     close(nf, status='keep')
  end if
end subroutine remove_annoying_files

subroutine set_sym(i)
  use m_Replica
  use m_Const_Parameters, only : DP
  use m_Crystal_Structure, only : nopr,op,tau
  implicit none
  integer, intent(in) :: i
  neb%image(i)%symmetry%nopr = nopr
  allocate(neb%image(i)%symmetry%op(3,3,nopr))
  allocate(neb%image(i)%symmetry%tau(3,nopr,2))
  neb%image(i)%symmetry%op = op
  neb%image(i)%symmetry%tau = tau
end subroutine set_sym

subroutine set_lowest_symmetry_op_tau()
  use m_Replica
  use m_Crystal_Structure, only : nopr,op,tau
  implicit none
  if(.not.allocated(op)) allocate(op(3,3,nopr_lowest))
  if(.not.allocated(tau)) allocate(tau(3,nopr_lowest,2))
  nopr = nopr_lowest
  op = op_lowest
  tau = tau_lowest
end subroutine set_lowest_symmetry_op_tau

subroutine set_sym_back(i)
  use m_Replica
  use m_Const_Parameters, only : DP
  use m_Crystal_Structure, only : nopr,op,tau
  implicit none
  integer, intent(in) :: i
  integer :: j
  nopr = neb%image(i)%symmetry%nopr
  do j=1,nopr
     op(1:3,1:3,j) = neb%image(i)%symmetry%op(1:3,1:3,j)
     tau(1:3,j,1) = neb%image(i)%symmetry%tau(1:3,j,1)
     tau(1:3,j,2) = neb%image(i)%symmetry%tau(1:3,j,2)
  enddo
end subroutine set_sym_back

logical function dump_checkpoint(itr)
  use m_Control_Parameters, only : cpt_neb_iteration, cpt_time, m_CtrlP_get_elpsd_time
  use m_Const_Parameters, only : DP
  use m_Parallelization, only : MPI_CommGroup,mype
  use m_Files, only : nfneb
  implicit none
  integer, intent(in) :: itr
  real(kind=DP), save :: last_time=0.d0
  real(kind=DP) :: curr_time
  integer :: ierr
  logical :: write_chkpnt, lneb, ltime
  write_chkpnt=.false.;lneb = .false.;ltime=.false.
  if(cpt_neb_iteration>0)then
     if(mod(itr,cpt_neb_iteration)==0) then
       write_chkpnt = .true.
       lneb = .true.
     endif
  endif
  curr_time = m_CtrlP_get_elpsd_time()
  if (cpt_time>0 .and. (curr_time-last_time)>cpt_time) then
     write_chkpnt = .true.
     ltime = .true.
     last_time = curr_time
     call mpi_barrier(MPI_CommGroup,ierr)
  endif
  if(write_chkpnt) then
     if(mype==0) then
        write(nfneb,'(a)')            '!** dumped checkpoint files because '
        if (lneb)  write(nfneb,'(a)') '!** neb iteration '
        if (ltime) write(nfneb,'(a)') '!** elapsed time '
        write(nfneb,'(a)')            '!** met the criterion'
     endif
  endif
  dump_checkpoint = write_chkpnt
end function dump_checkpoint

!   Nudged Elastic Band Method
!   Nobutaka Nishikawa

subroutine initialize_neb

  use m_Parallelization,      only : MPI_CommGroup, nrank_e, nrank_k, nrank_conf, mype_conf
  use m_Files,                only : nfneb, F_NEB_OUT
  use m_Replica
  use mpi

  implicit none
  integer i, j
  integer npes,mype,icolor,ikey,mpi_err
  integer, allocatable ::  new_comm_world(:)
!  integer nrank_r, mype_r
  integer mype_e, npes_e
  integer iargc
  character(100) arg
  character(4) image_name

!  include 'mpif.h'

  nrank_r = nrank_conf
  mype_r  = mype_conf
  neb%step = 0

!  write(6,*) 'program start...'
!  write(6,*)
!
!  call mpi_init(mpi_err)
!
!  call mpi_comm_size(mpi_comm_world,npes,mpi_err)
!!  call mpi_comm_rank(mpi_comm_world,mype,mpi_err)
!
!  nrank_r = 0
!  if(mype==0) then
!     do i = 1, iargc()
!        call getarg(i,arg)
!        if( arg(1:3) == 'nr=' ) then
!           read(arg(4:),*) nrank_r
!        end if
!     end do
!  end if
!  if(nrank_r == 0) nrank_r = npes
!  call mpi_bcast(nrank_r,1,mpi_integer,0,mpi_comm_world,mpi_err)
!
!  allocate(new_comm_world(0:nrank_r-1))
!  do i=0,nrank_r-1
!     icolor = 0
!     ikey = 0
!     do j=0, npes/nrank_r-1
!        if(mype == j + (npes/nrank_r)*i ) then
!           icolor = 1
!           ikey = i
!        end if
!     end do
!     call mpi_comm_split( mpi_comm_world, icolor, ikey, new_comm_world(i), mpi_err )
!  end do
!
!  mype_r = mype/(npes/nrank_r)
!  MPI_CommGroup = new_comm_world(mype_r)
!
!  call mpi_comm_size(MPI_CommGroup,npes_e,mpi_err)
!  call mpi_comm_rank(MPI_CommGroup,mype_e,mpi_err)
!
!  write(6,'("mpi: pe=",i4,",  npes=",i4,",  comm=",i20)') mype, npes, mpi_comm_world
!  write(6,'("   neb: ",i4,",  nr=",i4,"    energy: ",i4,",  ne*nk=",i4,",  comm=",i20)') &
!       mype_r, nrank_r, mype_e, npes_e, MPI_CommGroup
!  write(6,*)
!
!  neb%step = 0

end subroutine initialize_neb

subroutine initialize_neb_ep

  use m_Parallelization,      only : MPI_CommGroup
  use m_Files,                only : nfneb, F_NEB_OUT
  use m_Replica
  use mpi

  implicit none

  integer mype, npes, mpi_err

!  include 'mpif.h'

  call mpi_init(mpi_err)

  call mpi_comm_size(mpi_comm_world,npes,mpi_err)
  call mpi_comm_rank(mpi_comm_world,mype,mpi_err)

  MPI_CommGroup = mpi_comm_world

  neb%step = 0

  write(6,'(a6,4i5)') ' mpi: ', mype, npes, mpi_comm_world, MPI_CommGroup

end subroutine initialize_neb_ep

subroutine set_neb_condition

  use m_Const_Parameters, only : ON
  use m_Control_Parameters, only : icond, multiple_replica_max_iteration, sw_hubbard
  use m_Ionic_System,  only : neb_max_iteration, neb_dt, ci_neb &
                            , sp_k_init, sp_k_min, sp_k_max, sp_k_variable &
                            , penalty_function &
                            , neb_convergence_condition, neb_convergence_threshold &
                            , neb_time_integral,ci_index,ci_thres,sd_factor &
                            , end0_energy_given, end1_energy_given, to_dimer_thres
  use m_Files,         only : F_OUT, F_DYNM, F_ENF, F_CHGT, F_CNTN, F_CNTN_BIN &
                            , F_ZAJ, F_CHR, nfneb, F_NEB_OUT &
                            , nfnebenf, F_NEB_ENF, nfnebdynm, F_NEB_DYNM &
                            , nfdynm, nfenf, nfout
  use m_Parallelization,      only : MPI_CommGroup
  use m_Replica

! ======================== KT_add ================= 13.0D
  use m_PseudoPotential,  only : flg_paw
  use m_Files,         only : F_CNTN_BIN_PAW
! ================================================= 13.0D
  use m_Files,         only : F_OCCMAT

! ====== KT_add === 2014/07/20
  use m_Control_parameters, only : ipriparadeb, reuse_nfout_for_nfneb
  use m_Files,         only : F_OUT_BASE, F_CONF_EXTENSION, F_Para_EXTENSION
! ================= 2014/07/20

  use m_dimer,         only : set_dimer_condition
  use mpi

  implicit none

  integer mype, npes, mpi_err
  integer mype_e, npes_e
  character(4) image_name

  logical :: existence

!  include 'mpif.h'

  call mpi_comm_size(mpi_comm_world,npes,mpi_err)
  call mpi_comm_rank(mpi_comm_world,mype,mpi_err)
  call mpi_comm_size(MPI_CommGroup,npes_e,mpi_err)
  call mpi_comm_rank(MPI_CommGroup,mype_e,mpi_err)

! ====== KT_add ======== 2014/07/20
  if ( mype_e == 0 ) then
     close(nfout)
  else if(ipriparadeb == 0) then
     close(nfout, status='delete')
  else
     close(nfout, status='keep')
  end if
! ===================== 2014/07/20

! === KT_mod ========== 2014/07/20
#ifdef NEB_NEW_FILENAMES
  if ( reuse_nfout_for_nfneb == 1 ) then
     neb%file%F_NEB_OUT = trim(F_OUT_BASE) // trim(adjustl(F_CONF_EXTENSION)) &
          &                                // trim(adjustl(F_PARA_EXTENSION))
     if (mype_e==0) then
        open(nfneb,file=neb%file%F_NEB_OUT,status='unknown',position='append')
     endif
  else
!     neb%file%F_NEB_OUT = trim(F_NEB_OUT) // trim(adjustl(F_CONF_EXTENSION)) &
!          &                               // trim(adjustl(F_PARA_EXTENSION))
     neb%file%F_NEB_OUT = trim(F_OUT_BASE) // '_neb' &
          &                                // trim(adjustl(F_CONF_EXTENSION)) &
          &                                // trim(adjustl(F_PARA_EXTENSION))
     if (mype_e==0) open(nfneb,file=neb%file%F_NEB_OUT)
  endif
#else
  if(npes > 1) then
     write(image_name,'(a1,i3.3)') 'p', mype
     neb%file%F_NEB_OUT = trim(F_NEB_OUT)//'_'//image_name
  else
     neb%file%F_NEB_OUT = trim(F_NEB_OUT)
  end if
  open(nfneb,file=neb%file%F_NEB_OUT)
#endif
! ===================== 2014/07/20

  !write(nfneb,'(a6,4i5)') ' mpi: ', mype, npes, mpi_comm_world, MPI_CommGroup
  !write(nfneb,'(a6,8i8)') ' mpi: ', mype, npes, mpi_comm_world, &
  !     mype_r, nrank_r, &
  !     mype_e, npes_e, MPI_CommGroup

  neb%cond%condition = icond
  if (mype_e==0) write(nfneb,*) 'neb_condition: ', neb%cond%condition

  neb%cond%end0_energy_given = end0_energy_given
  neb%cond%end1_energy_given = end1_energy_given

  !neb%max_iteration = neb_max_iteration
  neb%max_iteration = multiple_replica_max_iteration
  if (mype_e==0) write(nfneb,*) 'max_iteration: ', neb%max_iteration

  neb%dt = neb_dt
  if (mype_e==0) write(nfneb,*) 'dt: ', neb%dt

  neb%cond%sd_factor = sd_factor

  neb%cond%k_spring_init = sp_k_init
  neb%cond%k_spring_min = sp_k_min
  neb%cond%k_spring_max = sp_k_max
  if (mype_e==0) write(nfneb,'(a11,3f10.5)') 'k_spring: ', &
       neb%cond%k_spring_init,neb%cond%k_spring_min,neb%cond%k_spring_max

  select case(sp_k_variable)
  case(0)
     neb%cond%k_variable = .false.
  case(1)
     neb%cond%k_variable = .true.
  end select
  if (mype_e==0) write(nfneb,*) 'sp_k_variable: ', neb%cond%k_variable

  select case(ci_neb)
  case(0)
     neb%cond%climbing_image = .false.
  case(1)
     neb%cond%climbing_image = .true.
  end select
  neb%cond%ci_index = ci_index
  neb%cond%ci_thres = ci_thres
  neb%cond%dimer_thres = to_dimer_thres
  if (mype_e==0) then
    write(nfneb,*) 'ci_neb: ', neb%cond%climbing_image
    if (neb%cond%ci_index>1 .and. neb%cond%ci_index<neb%number_of_images)then
      write(nfneb,*)'ci_index: ',neb%cond%ci_index
    endif
    write(nfneb,*)'ci_thres: ',neb%cond%ci_thres
    if(neb%cond%dimer_thres>0) then
      write(nfneb,'(a,f10.5)') 'switch to the dimer method if convergence is better than ' &
                      , neb%cond%dimer_thres
    endif
  endif

  if(neb%cond%dimer_thres>0) then
    call set_dimer_condition()
  endif

  neb%cond%ci_neb_start = 0
  neb%cond%ci_neb_end = 0

  select case(penalty_function)
  case(0)
     neb%cond%penalty_function = .false.
  case(1)
     neb%cond%penalty_function = .true.
  end select
  if (mype_e==0) write(nfneb,*) 'penalty_function: ', neb%cond%penalty_function

  neb%cond%k_damping = .false.
  neb%cond%k_damping_factor = 1.0d0

  neb%cond%convergence_condition = neb_convergence_condition
  if (mype_e==0) write(nfneb,*) 'convergence_condition: ', neb%cond%convergence_condition
  neb%cond%convergence_threshold = neb_convergence_threshold
  if (mype_e==0) write(nfneb,*) 'convergence_threshold: ', neb%cond%convergence_threshold

  !neb%cond%time_integral = 'verlet'
  neb%cond%time_integral = 'velocity_verlet'
  !neb%cond%time_integral = 'steepest_descent'
  select case(neb_time_integral)
  case(2)
     neb%cond%time_integral = 'velocity_verlet'
  case(12)
     neb%cond%time_integral = 'steepest_descent'
  case(20)
     neb%cond%time_integral = 'lbfgs'
  case(16)
     neb%cond%time_integral = 'bfgs'
  case(11)
     neb%cond%time_integral = 'cg'
  case(-7)
     neb%cond%time_integral = 'fire'
  end select
  if (mype_e==0) write(nfneb,*) 'time_integral: ', neb%cond%time_integral

  neb%file%F_OUT = F_OUT
  !neb%file%F_DYNM = F_DYNM
  !neb%file%F_ENF = F_ENF
  neb%file%F_CHGT = F_CHGT
  neb%file%F_CNTN = F_CNTN
  neb%file%F_CNTN_BIN = F_CNTN_BIN
  neb%file%F_ZAJ = F_ZAJ
  neb%file%F_CHR = F_CHR

! ======================== KT_add ================= 13.0D
  neb%file%F_CNTN_BIN_PAW = F_CNTN_BIN_PAW
! ================================================= 13.0D
  neb%file%F_OCCMAT = F_OCCMAT

  if (mype_e==0) then
     write(nfneb,*) 'F_OUT: ', trim(neb%file%F_OUT)
  !write(nfneb,*) 'F_DYNM: ', trim(neb%file%F_DYNM)
  !write(nfneb,*) 'F_ENF: ', trim(neb%file%F_ENF)
     write(nfneb,*) 'F_CHGT: ', trim(neb%file%F_CHGT)
     write(nfneb,*) 'F_CNTN: ', trim(neb%file%F_CNTN)
     write(nfneb,*) 'F_CNTN_BIN: ', trim(neb%file%F_CNTN_BIN)
     write(nfneb,*) 'F_ZAJ: ', trim(neb%file%F_ZAJ)
     write(nfneb,*) 'F_CHR: ', trim(neb%file%F_CHR)

! ======================== KT_add ================= 13.0D
     if ( flg_paw ) then
        write(nfneb,*) 'F_CNTN_BIN_PAW: ', trim(neb%file%F_CNTN_BIN_PAW)
     endif
! ================================================= 13.0D
     if (sw_hubbard==ON)then
        write(nfneb,*) 'F_OCCMAT: ', trim(neb%file%F_OCCMAT)
     endif
  end if

  if(mype == 0) then
     if(neb%cond%condition==1)then
       open(nfnebenf, file=F_NEB_ENF,access='append')
       open(nfnebdynm, file=F_NEB_DYNM,access='append')
     else
       open(nfnebenf, file=F_NEB_ENF)
       open(nfnebdynm, file=F_NEB_DYNM)

       write(nfnebenf,*) '#step  image  image_distance  energy   force_org  force_neb  force_normal'
       write(nfnebdynm,*) '#step  image  atom  cps'
     endif
  end if

end subroutine set_neb_condition

subroutine set_neb_parameter(itr)

  use m_Replica
  use m_Files, only : nfneb

  implicit none
  integer, intent(in) :: itr
  integer i

  do i = 1, neb%number_of_images
     if(itr==1.and.neb%image(i)%scf_convergence.and.neb%cond%condition==1) cycle
     neb%image(i)%energy0 = neb%image(i)%energy
     neb%image(i)%energy = 0.0d0
     neb%image(i)%force0 = neb%image(i)%force
     neb%image(i)%force = 0.0d0
     neb%image(i)%force_org = 0.0d0
!     neb%image(i)%scf_convergence = .true.
  end do

end subroutine set_neb_parameter

subroutine set_neb_replica_parameter(image_id)

  use m_Const_Parameters, only : ON
  use m_Control_parameters, only : icond, ipriparadeb, sw_hubbard
  use m_IterationNumbers, only : iteration,iteration_electronic,iteration_ionic
  use m_Replica
  use m_Files
  use m_Parallelization,  only : mype

! ======================== KT_add ================= 13.0D
  use m_PseudoPotential,  only : flg_paw
! ================================================= 13.0D

  use m_dimer, only : dimvars

  implicit none

  integer image_id
  character(4) image_name

  character(100) :: F_OUT_NEW

  if( neb%step > 1 .and. icond == 0) then
     icond = 1
  end if

  write(image_name,'(a1,i3.3)') 'r', image_id

  F_OUT = trim(neb%file%F_OUT)//'_'//image_name

! ===== KT_add === 2014/07/20
#ifdef NEB_NEW_FILENAMES
!  F_OUT = trim(F_OUT_BASE)//'_'//image_name &
!       &              //trim(F_CONF_EXTENSION)//trim(F_PARA_EXTENSION)
  F_OUT = trim(adjustl(F_OUT_BASE)) // '_' //image_name &
       &                             // trim(adjustl(F_PARA_EXTENSION))
! ================ 2014/07/20
#endif

  !F_DYNM = trim(neb%file%F_DYNM)//'_'//image_name
  !F_ENF = trim(neb%file%F_ENF)//'_'//image_name
  F_CHGT = trim(neb%file%F_CHGT)//'_'//image_name
  F_CNTN = trim(neb%file%F_CNTN)//'_'//image_name
  F_CNTN_BIN = trim(neb%file%F_CNTN_BIN)//'_'//image_name
  F_ZAJ = trim(neb%file%F_ZAJ)//'_'//image_name
  F_CHR = trim(neb%file%F_CHR)//'_'//image_name

! ======================== KT_add ================= 13.0D
  F_CNTN_BIN_PAW = trim(neb%file%F_CNTN_BIN_PAW)//'_'//image_name
! ================================================= 13.0D
  F_OCCMAT = trim(neb%file%F_OCCMAT)//'_'//image_name

  F_ZAJ_in      = F_ZAJ
  F_CHGT_in     = F_CHGT
  F_CNTN_BIN_in = F_CNTN_BIN
  F_CNTN_in     = F_CNTN

! ======================== KT_add ================= 13.0D
  F_CNTN_BIN_PAW_in = F_CNTN_BIN_PAW
! ================================================= 13.0D

  if ( mype == 0 ) then
     write(nfneb,*) 'F_OUT: ', trim(F_OUT)
  !write(nfneb,*) 'F_DYNM: ', trim(F_DYNM)
  !write(nfneb,*) 'F_ENF: ', trim(F_ENF)
     write(nfneb,*) 'F_CHGT: ', trim(F_CHGT)
     write(nfneb,*) 'F_CNTN: ', trim(F_CNTN)
     write(nfneb,*) 'F_CNTN_BIN: ', trim(F_CNTN_BIN)
     write(nfneb,*) 'F_ZAJ: ', trim(F_ZAJ)
     write(nfneb,*) 'F_CHR: ', trim(F_CHR)

! ======================== KT_add ================= 13.0D
     if ( flg_paw ) write(nfneb,*) 'F_CNTN_BIN_PAW: ', trim(F_CNTN_BIN_PAW)
! ================================================= 13.0D
     if (sw_hubbard==ON) write(nfneb,*) 'F_OCCMAT: ',trim(F_OCCMAT)
  endif

  if( neb%step == 0) then
     close(nfout)
     !close(nfdynm)
     !close(nfenf)
     open(nfout,file=trim(F_OUT))
     !open(nfdynm,file=trim(F_DYNM))
     !open(nfenf,file=trim(F_ENF))
  else
     close(nfout)
     !close(nfdynm)
     !close(nfenf)
     open(nfout,file=trim(F_OUT),position='append')
     !open(nfdynm,file=trim(F_DYNM),position='append')
     !open(nfenf,file=trim(F_ENF),position='append')
  end if

  close(nfcntn)
  close(nfcntn_bin)
  open(nfcntn,file=trim(F_CNTN))
  open(nfcntn_bin,file=trim(F_CNTN_BIN),form='unformatted')

! ======================== KT_add ================= 13.0D
  if ( flg_paw ) then
!     close(nfcntn_bin_paw)
!     open(nfcntn_bin_paw,file=trim(F_CNTN_BIN_PAW),form='unformatted')

     call m_Files_close_nfcntn_bin_paw()
     call m_Files_open_nfcntn_bin_paw()
  endif
! ================================================= 13.0D

  iteration = 0
  iteration_electronic = 0
  iteration_ionic = 0

end subroutine set_neb_replica_parameter

subroutine create_replica

  use m_Ionic_System,  only : number_of_replicas, natm &
                            , pos_end0, pos_end1, cps_end0, cps_end1 &
                            , replica_howtogive_coordinates, replica_endpoints &
                            , pos_image, cps_image &
                            , imdtyp, sw_path_from_dynm, end0_energy, end1_energy
  use m_Files,         only : nfneb, nfnebcntn, F_NEB_CNTN
  use m_Const_Parameters,  only : FILE, PROPORTIONAL, ON, PREPARATION_ONLY
  use m_Replica
  use m_Parallelization,  only : mype
  use m_Control_Parameters, only : sw_fcp, icond
  use m_Fcp,                only : m_Fcp_create_replica, m_Fcp_read_totch_replica

  implicit none

  integer i, j
  integer end0, end1
  character token

  neb%number_of_images = number_of_replicas + 2
  if (mype==0) write(nfneb,*) 'number_of_images: ', neb%number_of_images

  ! temporary
!  allocate(neb%energy(neb%max_iteration,neb%number_of_images))
!  neb%energy = 0.0d0

  allocate(neb%image(neb%number_of_images))
  do i = 1, neb%number_of_images
     neb%image(i)%num_atom = natm

     neb%image(i)%energy = 0.0d0
     neb%image(i)%energy0 = 0.0d0

     allocate(neb%image(i)%pos(neb%image(i)%num_atom,3))
     allocate(neb%image(i)%pos0(neb%image(i)%num_atom,3))
     allocate(neb%image(i)%cps(neb%image(i)%num_atom,3))
     allocate(neb%image(i)%cps0(neb%image(i)%num_atom,3))
     allocate(neb%image(i)%tau(neb%image(i)%num_atom,3))
     allocate(neb%image(i)%force(neb%image(i)%num_atom,3))
     allocate(neb%image(i)%force0(neb%image(i)%num_atom,3))
     allocate(neb%image(i)%force_org(neb%image(i)%num_atom,3))
     allocate(neb%image(i)%spring_force(neb%image(i)%num_atom,3))
     allocate(neb%image(i)%true_force(neb%image(i)%num_atom,3))
     allocate(neb%image(i)%fix_flag(neb%image(i)%num_atom))
     allocate(neb%image(i)%velocity(neb%image(i)%num_atom,3))
     allocate(neb%image(i)%velocity0(neb%image(i)%num_atom,3))
     neb%image(i)%pos = 0.0d0
     neb%image(i)%pos0 = 0.0d0
     neb%image(i)%cps = 0.0d0
     neb%image(i)%cps0 = 0.0d0
     neb%image(i)%tau = 0.0d0
     neb%image(i)%force = 0.0d0
     neb%image(i)%force0 = 0.0d0
     neb%image(i)%force_org = 0.0d0
     neb%image(i)%spring_force = 0.0d0
     neb%image(i)%true_force = 0.0d0
     neb%image(i)%fix_flag = 0
     neb%image(i)%velocity = 0.0d0
     neb%image(i)%velocity0 = 0.0d0

     neb%image(i)%scf_convergence = .true.
  end do

  end0 = 1
  end1 = neb%number_of_images

  neb%image(end0)%pos = pos_end0
  neb%image(end0)%cps = cps_end0
  neb%image(end1)%pos = pos_end1
  neb%image(end1)%cps = cps_end1

  do i = 2, neb%number_of_images-1
     do j = 1, neb%image(i)%num_atom
        neb%image(i)%cps(j,:) = neb%image(end0)%cps(j,:) &
             + real(i-1) * (neb%image(end1)%cps(j,:) - neb%image(end0)%cps(j,:)) &
             / real(end1-end0)
        neb%image(i)%pos(j,:) = neb%image(end0)%pos(j,:) &
             + real(i-1) * (neb%image(end1)%pos(j,:) - neb%image(end0)%pos(j,:)) &
             / real(end1-end0)
     end do
  end do

  do i = 2, neb%number_of_images-1
     if(replica_howtogive_coordinates(i-1) == FILE .or. sw_path_from_dynm==ON) then
        neb%image(i)%pos = pos_image(i-1,:,:)
        neb%image(i)%cps = cps_image(i-1,:,:)
     end if
  end do

  neb%image(1)%fix_flag(:) = 0
  neb%image(neb%number_of_images)%fix_flag(:) = 0
  do i = 2, neb%number_of_images-1
        neb%image(i)%fix_flag(:) = imdtyp(:)
  end do
!  imdtyp = 0

  if( sw_fcp == ON ) call m_Fcp_create_replica(neb%number_of_images)

  if(neb%cond%condition == 1) then
     if (mype==0) write(nfneb,*) 'read data replica...'

     open(nfnebcntn,file=F_NEB_CNTN)
     read(nfnebcntn,*)
     read(nfnebcntn,*) neb%step
     read(nfnebcntn,*)
     read(nfnebcntn,*)

     do i = 1, neb%number_of_images
        read(nfnebcntn,*)
        read(nfnebcntn,*)
        read(nfnebcntn,*)
        read(nfnebcntn,*)
        read(nfnebcntn,*)
        do j = 1, neb%image(i)%num_atom
           read(nfnebcntn,*) neb%image(i)%pos(j,1:3)
        end do
        read(nfnebcntn,*)
        do j = 1, neb%image(i)%num_atom
           read(nfnebcntn,*) neb%image(i)%cps(j,1:3)
        end do
        read(nfnebcntn,*)
        if(mod(i-1,nrank_r) == mype_r) then
          read(nfnebcntn,*) neb%image(i)%energy
          neb%image(i)%energy0 = neb%image(i)%energy
          read(nfnebcntn,*)
        else
          read(nfnebcntn,*)
          read(nfnebcntn,*)
        endif
        do j = 1, neb%image(i)%num_atom
           read(nfnebcntn,*) neb%image(i)%force_org(j,1:3)
        end do
        neb%image(i)%force = neb%image(i)%force_org
        read(nfnebcntn,*)
        read(nfnebcntn,*) token
        select case(token)
        case('0')
           if (mype==0) write(nfneb,*) 'image: ', i, '  scf not convegence'
           neb%image(i)%scf_convergence = .false.
        case('1')
           if (mype==0) write(nfneb,*) 'image: ', i, '  scf convegence'
           neb%image(i)%scf_convergence = .true.
        end select
        if( sw_fcp == ON ) call m_Fcp_read_totch_replica(i)
     end do
     read(nfnebcntn,*,err=90,end=90)
     read(nfnebcntn,*,err=90,end=90) changed_to_cineb
90   continue

     close(nfnebcntn)
  end if

  if (mype==0) then
     write(nfneb,*) 'create replica...'
     do i = 1, neb%number_of_images
        write(nfneb,*) 'image: ', i, ' pos, cps'
        do j = 1, neb%image(i)%num_atom
           write(nfneb,'(6f10.5,i3)') neb%image(i)%pos(j,1:3), neb%image(i)%cps(j,1:3), &
                neb%image(i)%fix_flag(j)
        end do
     end do
  endif

  if(neb%cond%end0_energy_given)then
    if(mod(0,nrank_r) == mype_r) then
      neb%image(1)%energy0 = end0_energy
      neb%image(1)%energy = end0_energy
      write(nfneb,*) 'energy for end0 given from the input file : ',end0_energy
    endif
  endif
  if(neb%cond%end1_energy_given)then
    if(mod(neb%number_of_images-1,nrank_r) == mype_r) then
      neb%image(neb%number_of_images)%energy0 = end1_energy
      neb%image(neb%number_of_images)%energy = end1_energy
      write(nfneb,*) 'energy for end1 given from the input file : ',end1_energy
    endif
  endif

  if(neb%cond%condition==0 .or. icond == PREPARATION_ONLY) call write_result(.false.)

end subroutine create_replica

subroutine set_neb_replica_coordinate(image_id)

  use m_Ionic_System, only : pos, cps
  use m_Replica

  implicit none
  integer image_id

  cps = neb%image(image_id)%cps
  pos = neb%image(image_id)%pos

end subroutine set_neb_replica_coordinate

subroutine allreduce_neb_energy_force

  use m_Replica
  use m_Parallelization,      only : MPI_CommGroup, mype_conf
  use mpi

  implicit none

  integer i
  integer mype, npes, mpi_err
  real(8),allocatable :: work1(:)
  real(8),allocatable :: work2(:,:)

!  include 'mpif.h'

  if (nrank_r.lt.1) return

!  call mpi_comm_rank(mpi_comm_world,mype,mpi_err)
  call mpi_comm_size(MPI_CommGroup,npes,mpi_err)

!  write(220+mype,'(i5,5e20.10)') mype, neb%image(:)%energy

  allocate(work1(neb%number_of_images));work1=0.d0
  call mpi_allreduce(neb%image(:)%energy,work1,neb%number_of_images, &
       mpi_double_precision,mpi_sum,mpi_comm_world,mpi_err)
  neb%image(:)%energy = work1
  deallocate(work1)
  do i=1, neb%number_of_images
    allocate(work2(neb%image(i)%num_atom,3))
     work2=0.d0
     call mpi_allreduce(neb%image(i)%force,work2,neb%image(i)%num_atom*3, &
          mpi_double_precision,mpi_sum,mpi_comm_world,mpi_err)
     neb%image(i)%force = work2
     work2=0.d0
     call mpi_allreduce(neb%image(i)%force_org,work2,neb%image(i)%num_atom*3, &
          mpi_double_precision,mpi_sum,mpi_comm_world,mpi_err)
     neb%image(i)%force_org = work2
     deallocate(work2)
  end do
  ! temporary    npes = ne*nk
  neb%image(:)%energy = neb%image(:)%energy/npes
  do i=1, neb%number_of_images
     neb%image(i)%force = neb%image(i)%force/npes
     neb%image(i)%force_org = neb%image(i)%force_org/npes
  end do

!  write(230+mype,'(i5,5e20.10)') mype, neb%image(:)%energy

end subroutine allreduce_neb_energy_force

subroutine set_neb_replica_energy_force(image_id)

  use m_Force, only : forc_l
  use m_Total_Energy, only : etotal
  use m_Replica
  use m_Files, only : nfneb

  implicit none

  integer image_id

  integer :: i

  neb%image(image_id)%energy = etotal
  neb%image(image_id)%force = forc_l
  neb%image(image_id)%force_org = forc_l

end subroutine set_neb_replica_energy_force

subroutine force_neb

  use m_Files,                only : nfneb
  use m_Replica
  use m_Parallelization,  only : mype

  implicit none

  integer i, j

  if ( mype==0 ) then
     write(nfneb,*) 'force neb...'

     do i = 1, neb%number_of_images
        write(nfneb,*) 'image: ', i, neb%image(i)%energy, &
             neb%image(i)%energy-neb%image(i)%energy0
        write(nfneb,*) ' cps, force '
        do j = 1, neb%image(i)%num_atom
           write(nfneb,'(3f10.5,3e11.3)') neb%image(i)%cps(j,1:3), neb%image(i)%force(j,1:3)
        end do
     end do
  endif

!  do i = 1, neb%number_of_images
!     neb%image(i)%force_org = neb%image(i)%force
!  end do

  call elastic_constant
  call local_tangent
  call neb_force

  if ( mype==0 ) then
     do i = 1, neb%number_of_images
        write(nfneb,*) 'image: ', i, neb%image(i)%energy, &
             neb%image(i)%energy-neb%image(i)%energy0
!     write(nfneb,*) ' cps, pos '
!     do j = 1, neb%image(i)%num_atom
!        write(nfneb,'(6f10.5)') neb%image(i)%cps(j,1:3), neb%image(i)%pos(j,1:3)
!     end do
        write(nfneb,*) ' force(neb), force(org), force(neb-org) '
        do j = 1, neb%image(i)%num_atom
           write(nfneb,'(9e11.3)') neb%image(i)%force(j,1:3), neb%image(i)%force_org(j,1:3), &
                neb%image(i)%force(j,1:3)-neb%image(i)%force_org(j,1:3)
        end do
     end do
  endif
end subroutine force_neb

subroutine elastic_constant

  ! elastic constant between images

  use m_Files,                only : nfneb
  use m_Replica
  use m_Parallelization,     only : mype

  implicit none

  integer  i

  real(8)  max_energy, ref_energy, energy

  real(8)  k_spring, k_spring0
  real(8)  k_spring_init, k_spring_max, k_spring_min, k_spring_delta

  if (mype==0) write(nfneb,*) 'elastic constant...'

  max_energy = maxval( neb%image(:)%energy )
  ref_energy = max( neb%image(1)%energy, neb%image(neb%number_of_images)%energy )

  do i = 2, neb%number_of_images

    if( .not. neb%cond%k_variable ) then
       k_spring = neb%cond%k_spring_init
       k_spring0= neb%cond%k_spring_init
    elseif( neb%cond%k_variable ) then
       k_spring_max = neb%cond%k_spring_max
       k_spring_min = neb%cond%k_spring_min
       energy = max( neb%image(i-1)%energy, neb%image(i)%energy )
       k_spring_delta = k_spring_max - k_spring_min
       if( energy > ref_energy ) then
          k_spring = k_spring_max - k_spring_delta &
               * ( max_energy - energy ) / ( max_energy - ref_energy )
          k_spring0 = k_spring
       elseif( energy <= ref_energy ) then
          k_spring = k_spring_max - k_spring_delta
          k_spring0 = k_spring
       end if
    end if

    if( neb%cond%k_damping ) then
       k_spring0 = k_spring
       k_spring = k_spring * neb%cond%k_damping_factor ** ( ( neb%step-1) )
    end if

    neb%image(i-1)%k_spring = k_spring
    neb%image(i-1)%k_spring0= k_spring0

    if (mype==0) write(nfneb,'(i5,a10,2e20.12)') i, '-th spring',  k_spring, k_spring0

 end do

 return
end subroutine elastic_constant

subroutine local_tangent

  ! local tangent of minimum energy path

  use m_Files,                only : nfneb
  use m_Replica
  use m_Parallelization,     only : mype

  implicit none

  integer i, j
  integer num_atom
  real(8) distance21
  real(8) distance01, distance20
  real(8) abs_tau2, abs_tau3
  real(8) angle12, angle23, angle31
  real(8) energy0, energy1, energy2
  real(8),allocatable :: tau(:,:)
  real(8),allocatable :: tau1(:,:), tau2(:,:), tau3(:,:)
  real(8) phi
  real(8),allocatable ::  rx0(:), ry0(:), rz0(:)
  real(8),allocatable ::  rx1(:), ry1(:), rz1(:)
  real(8),allocatable ::  rx2(:), ry2(:), rz2(:)
  real(8) :: eps

  if (mype==0) write(nfneb,*) 'local tangent...'

  eps = 1.d-10
  do j = 2, neb%number_of_images-1

     num_atom = neb%image(j)%num_atom

     allocate( tau(num_atom,3) )
     allocate( tau1(num_atom,3) )
     allocate( tau2(num_atom,3) )
     allocate( tau3(num_atom,3) )

     allocate( rx0(num_atom) )
     allocate( ry0(num_atom) )
     allocate( rz0(num_atom) )
     allocate( rx1(num_atom) )
     allocate( ry1(num_atom) )
     allocate( rz1(num_atom) )
     allocate( rx2(num_atom) )
     allocate( ry2(num_atom) )
     allocate( rz2(num_atom) )

     tau1(:,:) = 0.0d0
     tau2(:,:) = 0.0d0
     tau3(:,:) = 0.0d0

     rx0 = neb%image(j)%cps(:,1)
     ry0 = neb%image(j)%cps(:,2)
     rz0 = neb%image(j)%cps(:,3)
     rx1 = neb%image(j-1)%cps(:,1)
     ry1 = neb%image(j-1)%cps(:,2)
     rz1 = neb%image(j-1)%cps(:,3)
     rx2 = neb%image(j+1)%cps(:,1)
     ry2 = neb%image(j+1)%cps(:,2)
     rz2 = neb%image(j+1)%cps(:,3)

     energy0 = neb%image(j)%energy
     energy1 = neb%image(j-1)%energy
     energy2 = neb%image(j+1)%energy

     distance21 = sqrt ( dot_product( rx2 - rx1, rx2 - rx1 ) &
          + dot_product( ry2 - ry1, ry2 - ry1 ) &
          + dot_product( rz2 - rz1, rz2 - rz1 ) )
     distance01 = sqrt ( dot_product( rx0 - rx1, rx0 - rx1 ) &
          + dot_product( ry0 - ry1, ry0 - ry1 ) &
          + dot_product( rz0 - rz1, rz0 - rz1 ) )
     distance20 = sqrt ( dot_product( rx2 - rx0, rx2 - rx0 ) &
          + dot_product( ry2 - ry0, ry2 - ry0 ) &
          + dot_product( rz2 - rz0, rz2 - rz0 ) )

     if( distance21.lt.eps .or. distance01.lt.eps .or. distance20.lt.eps ) cycle

     phi = ( dot_product( rx2 - rx0, rx0 - rx1 ) &
          + dot_product( ry2 - ry0, ry0 - ry1 ) &
          + dot_product( rz2 - rz0, rz0 - rz1 ) ) / ( distance20 * distance01 )

     ! for regure NEB

     do i = 1, num_atom
        tau1(i, 1) = ( rx2(i) - rx1(i)) / distance21
        tau1(i, 2) = ( ry2(i) - ry1(i)) / distance21
        tau1(i, 3) = ( rz2(i) - rz1(i)) / distance21
        tau2(i, 1) = ( rx0(i) - rx1(i)) / distance01 &
             + ( rx2(i) - rx0(i) ) / distance20
        tau2(i, 2) = ( ry0(i) - ry1(i)) / distance01 &
             + ( ry2(i) - ry0(i) ) / distance20
        tau2(i, 3) = ( rz0(i) - rz1(i)) / distance01 &
             + ( rz2(i) - rz0(i) ) / distance20
     end do

     abs_tau2 = sqrt( dot_product(tau2(:,1), tau2(:,1)) &
          + dot_product(tau2(:,2), tau2(:,2)) &
          + dot_product(tau2(:,3), tau2(:,3)))
     if( abs_tau2 /= 0.0d0 ) then
        tau2(:, 1) = tau2(:, 1) /  abs_tau2
        tau2(:, 2) = tau2(:, 2) /  abs_tau2
        tau2(:, 3) = tau2(:, 3) /  abs_tau2
     end if

!       tau(:,:) = tau1(:,:)
!       tau(:,:) = tau2(:,:)

     ! for CI-NEB  ( method of python )

     if( energy2 > energy0 ) then
        if( energy0 > energy1 ) then
           do i = 1, num_atom
              tau3(i, 1) = ( rx2(i) - rx0(i) )
              tau3(i, 2) = ( ry2(i) - ry0(i) )
              tau3(i, 3) = ( rz2(i) - rz0(i) )
           end do
        elseif ( energy0 <= energy1 ) then
           if( energy2 > energy1 ) then
              do i = 1, num_atom
                 tau3(i, 1) = ( rx2(i) - rx0(i) ) &
                      * max( abs( energy2-energy0 ), abs( energy1 -energy0 ) ) &
                      + ( rx0(i) - rx1(i) ) &
                      * min( abs( energy2-energy0 ), abs( energy1 -energy0 ) )
                 tau3(i, 2) = ( ry2(i) - ry0(i) ) &
                      * max( abs( energy2-energy0 ), abs( energy1 -energy0 ) ) &
                      + ( ry0(i) - ry1(i) ) &
                      * min( abs( energy2-energy0 ), abs( energy1 -energy0 ) )
                 tau3(i, 3) = ( rz2(i) - rz0(i) ) &
                      * max( abs( energy2-energy0 ), abs( energy1 -energy0 ) ) &
                      + ( rz0(i) - rz1(i) ) &
                      * min( abs( energy2-energy0 ), abs( energy1 -energy0 ) )
              end do
           else if( energy2 <= energy1 ) then
              do i = 1, num_atom
                 tau3(i, 1) = ( rx2(i) - rx0(i) ) &
                      * min( abs( energy2-energy0 ), abs( energy1 -energy0 ) ) &
                      + ( rx0(i) - rx1(i) ) &
                      * max( abs( energy2-energy0 ), abs( energy1 -energy0 ) )
                 tau3(i, 2) = ( ry2(i) - ry0(i) ) &
                      * min( abs( energy2-energy0 ), abs( energy1 -energy0 ) ) &
                      + ( ry0(i) - ry1(i) ) &
                      * max( abs( energy2-energy0 ), abs( energy1 -energy0 ) )
                 tau3(i, 3) = ( rz2(i) - rz0(i) ) &
                      * min( abs( energy2-energy0 ), abs( energy1 -energy0 ) ) &
                      + ( rz0(i) - rz1(i) ) &
                      * max( abs( energy2-energy0 ), abs( energy1 -energy0 ) )
              end do
           end if
        end if
     elseif( energy2 <= energy0 ) then
        if( energy0 < energy1 ) then
           do i = 1, num_atom
              tau3(i, 1) = ( rx0(i) - rx1(i) )
              tau3(i, 2) = ( ry0(i) - ry1(i) )
              tau3(i, 3) = ( rz0(i) - rz1(i) )
           end do
        elseif( energy0 >= energy1 ) then
           if( energy2 > energy1 ) then
              do i = 1, num_atom
                 tau3(i, 1) = ( rx2(i) - rx0(i) ) &
                      * max( abs( energy2-energy0 ), abs( energy1 -energy0 ) ) &
                      + ( rx0(i) - rx1(i) ) &
                      * min( abs( energy2-energy0 ), abs( energy1 -energy0 ) )
                 tau3(i, 2) = ( ry2(i) - ry0(i) ) &
                      * max( abs( energy2-energy0 ), abs( energy1 -energy0 ) ) &
                      + ( ry0(i) - ry1(i) ) &
                      * min( abs( energy2-energy0 ), abs( energy1 -energy0 ) )
                 tau3(i, 3) = ( rz2(i) - rz0(i) ) &
                      * max( abs( energy2-energy0 ), abs( energy1 -energy0 ) ) &
                      + ( rz0(i) - rz1(i) ) &
                      * min( abs( energy2-energy0 ), abs( energy1 -energy0 ) )
              end do
           elseif( energy2 <= energy1 ) then
              do i = 1, num_atom
                 tau3(i, 1) = ( rx2(i) - rx0(i) ) &
                      * min( abs( energy2-energy0 ), abs( energy1 -energy0 ) ) &
                      + ( rx0(i) - rx1(i) ) &
                      * max( abs( energy2-energy0 ), abs( energy1 -energy0 ) )
                 tau3(i, 2) = ( ry2(i) - ry0(i) ) &
                      * min( abs( energy2-energy0 ), abs( energy1 -energy0 ) ) &
                      + ( ry0(i) - ry1(i) ) &
                      * max( abs( energy2-energy0 ), abs( energy1 -energy0 ) )
                 tau3(i, 3) = ( rz2(i) - rz0(i) ) &
                      * min( abs( energy2-energy0 ), abs( energy1 -energy0 ) ) &
                      + ( rz0(i) - rz1(i) ) &
                      * max( abs( energy2-energy0 ), abs( energy1 -energy0 ) )
              end do
           end if
        end if
     end if

     abs_tau3 = sqrt( dot_product(tau3(:,1), tau3(:,1)) &
          + dot_product(tau3(:,2), tau3(:,2)) &
          + dot_product(tau3(:,3), tau3(:,3)))
     if( abs_tau3 /= 0.0d0 ) then
        tau3(:, 1) = tau3(:, 1) /  abs_tau3
        tau3(:, 2) = tau3(:, 2) /  abs_tau3
        tau3(:, 3) = tau3(:, 3) /  abs_tau3
     end if

     tau(:,:) = tau3(:,:)

     angle12 = ( dot_product( tau1(:,1), tau2(:,1) ) &
          + dot_product( tau1(:,2), tau2(:,2) ) &
          + dot_product( tau1(:,3), tau2(:,3) ) )
     angle23 = ( dot_product( tau2(:,1), tau3(:,1) ) &
          + dot_product( tau2(:,2), tau3(:,2) ) &
          + dot_product( tau2(:,3), tau3(:,3) ) )

     angle31 = ( dot_product( tau3(:,1), tau1(:,1) ) &
          + dot_product( tau3(:,2), tau1(:,2) ) &
          + dot_product( tau3(:,3), tau1(:,3) ) )

     ! output coordinate and local tangent

     if (mype==0) then
        write(nfneb,*) 'image : ', j
        write(nfneb,*) 'angle between local tangent'
        write(nfneb,'(a11,e12.5)') ' angle12 : ', angle12
        write(nfneb,'(a11,e12.5)') ' angle23 : ', angle23
        write(nfneb,'(a11,e12.5)') ' angle31 : ', angle31
        write(nfneb,*) ' cps '
        do i = 1, num_atom
           write(nfneb,'(3f10.5)') rx0(i), ry0(i), rz0(i)
        end do
        write(nfneb,*) ' local tangent : tau1, tau2, tau3'
        do i = 1, num_atom
           write(nfneb,'(9e11.3)') tau1(i,1), tau1(i,2), tau1(i,3), &
                tau2(i,1), tau2(i,2), tau2(i,3), &
                tau3(i,1), tau3(i,2), tau3(i,3)
        end do
!     write(nfneb,'(a19,3f15.5)') ' phi, dis01, dis20 ', phi, distance01, distance20
        write(nfneb,*) ''
     endif

     do i = 1, num_atom
        neb%image(j)%tau(i,:) = tau(i,:)
     end do

     deallocate( tau, tau1, tau2, tau3 )
     deallocate( rx0, ry0, rz0, rx1, ry1, rz1, rx2, ry2, rz2 )

  end do

  return
end subroutine local_tangent

subroutine neb_force

  ! constrain force of atoms of CI-NEB method

  use m_Files,                only : nfneb
  use m_Replica
  use m_Const_Parameters, only : PAI
  use m_Parallelization,     only : mype

  implicit none

  integer  num_atom
!  integer, intent(in) :: num_image
  integer  image_id

  real(8) energy
  real(8) max_energy

  real(8),allocatable :: tau(:,:)
  real(8),allocatable :: k_spring(:), k_spring0(:)

  integer i, k, iat

!  integer :: mype
!  integer :: energy_index(num_image)
  integer,allocatable :: update_id(:)
  real(8) :: distance01, distance20
  real(8) :: dot_product_f0_tau
  real(8) :: f_spring_abs
  real(8) :: f_spring_para_sum = 0.0d0
  real(8) :: f_spring_perp_sum = 0.0d0
  real(8) :: f_true_perp_sum = 0.0d0
  real(8) :: f0_perp_sum = 0.0d0, f0_para_sum = 0.0d0, f_sum = 0.0d0
  real(8) :: cosphi, switch = 0.0d0
  real(8), allocatable :: fx_spring(:), fy_spring(:), fz_spring(:)
  real(8), allocatable :: fx0_para(:), fy0_para(:), fz0_para(:)
  real(8), allocatable :: fx0_perp(:), fy0_perp(:), fz0_perp(:)
  real(8), allocatable :: fx_spring_para(:), fy_spring_para(:), fz_spring_para(:)
  real(8), allocatable :: fx_spring_perp(:), fy_spring_perp(:), fz_spring_perp(:)
  real(8), allocatable :: rx0(:), ry0(:), rz0(:)
  real(8), allocatable :: rx1(:), ry1(:), rz1(:)
  real(8), allocatable :: rx2(:), ry2(:), rz2(:)
  real(8), allocatable :: fx0(:), fy0(:), fz0(:)
  real(8), allocatable :: fx1(:), fy1(:), fz1(:)

  integer :: ci_target

  if (mype==0) write(nfneb,*) 'neb force...'

  allocate( k_spring(neb%number_of_images-1) )
  allocate( k_spring0(neb%number_of_images-1) )

  do k = 1, neb%number_of_images-1
     k_spring(k) = neb%image(k)%k_spring
     k_spring0(k) = neb%image(k)%k_spring0
  end do

  max_energy = maxval( neb%image(:)%energy )
  ci_target=0
  if (neb%cond%climbing_image .and. .not. changed_to_cineb) then
    write(nfneb,*) 'curr thres : ',curr_thres
  endif
  if( ( neb%step >= neb%cond%ci_neb_start) .and. ( neb%cond%climbing_image )  &
  &   .and. (curr_thres<neb%cond%ci_thres) .and. .not. changed_to_cineb ) then
     if(.not.changed_to_cineb) write(nfneb,*) 'transitioning to CI-NEB'
     changed_to_cineb = .true.
  endif
  if(changed_to_cineb) then
     if (neb%cond%ci_index>1 .and. neb%cond%ci_index<neb%number_of_images)then
        ci_target = neb%cond%ci_index
        if(mype==0) write(nfneb,*) 'climbing image : ',ci_target
     else
        do k=2,neb%number_of_images-1
           if (neb%image(k)%energy==max_energy) then
              ci_target=k
              if(mype==0) write(nfneb,*) 'climbing image : ',ci_target
              exit
           endif
        enddo
     endif
  endif

  do k = 2, neb%number_of_images-1

     num_atom = neb%image(k)%num_atom
     image_id = k

     allocate( tau(num_atom,3) )
     allocate( update_id(num_atom) )

     allocate( rx0(num_atom) )
     allocate( ry0(num_atom) )
     allocate( rz0(num_atom) )
     allocate( rx1(num_atom) )
     allocate( ry1(num_atom) )
     allocate( rz1(num_atom) )
     allocate( rx2(num_atom) )
     allocate( ry2(num_atom) )
     allocate( rz2(num_atom) )
     allocate( fx0(num_atom) )
     allocate( fy0(num_atom) )
     allocate( fz0(num_atom) )
     allocate( fx1(num_atom) )
     allocate( fy1(num_atom) )
     allocate( fz1(num_atom) )
     allocate( fx_spring(num_atom) )
     allocate( fy_spring(num_atom) )
     allocate( fz_spring(num_atom) )
     allocate( fx0_para(num_atom) )
     allocate( fy0_para(num_atom) )
     allocate( fz0_para(num_atom) )

     allocate( fx0_perp(num_atom) )
     allocate( fy0_perp(num_atom) )
     allocate( fz0_perp(num_atom) )

     allocate( fx_spring_para(num_atom) )
     allocate( fy_spring_para(num_atom) )
     allocate( fz_spring_para(num_atom) )

     allocate( fx_spring_perp(num_atom) )
     allocate( fy_spring_perp(num_atom) )
     allocate( fz_spring_perp(num_atom) )

     rx0 = neb%image(image_id)%cps(:,1)
     ry0 = neb%image(image_id)%cps(:,2)
     rz0 = neb%image(image_id)%cps(:,3)
     rx1 = neb%image(image_id-1)%cps(:,1)
     ry1 = neb%image(image_id-1)%cps(:,2)
     rz1 = neb%image(image_id-1)%cps(:,3)
     rx2 = neb%image(image_id+1)%cps(:,1)
     ry2 = neb%image(image_id+1)%cps(:,2)
     rz2 = neb%image(image_id+1)%cps(:,3)

     fx0 = neb%image(image_id)%force(:,1)
     fy0 = neb%image(image_id)%force(:,2)
     fz0 = neb%image(image_id)%force(:,3)

     tau(:,1) = neb%image(image_id)%tau(:,1)
     tau(:,2) = neb%image(image_id)%tau(:,2)
     tau(:,3) = neb%image(image_id)%tau(:,3)

     update_id = neb%image(image_id)%fix_flag(:)

!    call get_energy_index( energy_index, num_image )

     ! spring force for penalty function

     fx_spring(:) = k_spring0(image_id) * ( rx2(:) - rx0(:) ) &
          + k_spring0(image_id - 1 ) * ( rx1(:) - rx0(:) )
     fy_spring(:) = k_spring0(image_id) * ( ry2(:) - ry0(:) ) &
          + k_spring0(image_id - 1 ) * ( ry1(:) - ry0(:) )
     fz_spring(:) = k_spring0(image_id) * ( rz2(:) - rz0(:) ) &
          + k_spring0(image_id - 1 ) * ( rz1(:) - rz0(:) )

     if(k/=ci_target)then

        ! get true force perpendicular to local tangent

        dot_product_f0_tau &
             = dot_product( fx0(:),  tau(:,1) ) &
             + dot_product( fy0(:),  tau(:,2) ) &
             + dot_product( fz0(:),  tau(:,3) )

        fx0_perp(:) = fx0(:) - dot_product_f0_tau * tau(:,1)
        fy0_perp(:) = fy0(:) - dot_product_f0_tau * tau(:,2)
        fz0_perp(:) = fz0(:) - dot_product_f0_tau * tau(:,3)

        ! get spring force along local tangent

        distance20 = sqrt( dot_product( rx2 - rx0, rx2 - rx0 ) &
             + dot_product( ry2 - ry0, ry2 - ry0 ) &
             + dot_product( rz2 - rz0, rz2 - rz0 ) )
        distance01 = sqrt( dot_product( rx0 - rx1, rx0 - rx1 ) &
             + dot_product( ry0 - ry1, ry0 - ry1 ) &
             + dot_product( rz0 - rz1, rz0 - rz1 ) )
        f_spring_abs = k_spring0(image_id) *  distance20 &
             - k_spring0(image_id - 1) * distance01

        fx_spring_para(:) = ( dot_product(fx_spring, tau(:,1)) &
             + dot_product(fy_spring, tau(:,2)) &
             + dot_product(fz_spring, tau(:,3)) ) * tau(:,1)
        fy_spring_para(:) = ( dot_product(fx_spring, tau(:,1)) &
             + dot_product(fy_spring, tau(:,2)) &
             + dot_product(fz_spring, tau(:,3)) ) * tau(:,2)
        fz_spring_para(:) = ( dot_product(fx_spring, tau(:,1)) &
             + dot_product(fy_spring, tau(:,2)) &
             + dot_product(fz_spring, tau(:,3)) ) * tau(:,3)

        ! penalty function

        switch=0.d0
        fx_spring_perp=0;fy_spring_perp=0;fz_spring_perp=0;cosphi=0
        if( neb%cond%penalty_function )  then
           fx_spring_perp(:) = ( fx_spring(:) - ( dot_product(fx_spring, tau(:,1)) &
                + dot_product(fy_spring, tau(:,2)) &
                + dot_product(fz_spring, tau(:,3)) ) * tau(:,1) ) &
                * k_spring(image_id-1) / k_spring0(image_id-1)
           fy_spring_perp(:) = ( fy_spring(:) - ( dot_product(fx_spring, tau(:,1)) &
                + dot_product(fy_spring, tau(:,2)) &
                + dot_product(fz_spring, tau(:,3)) ) * tau(:,2) ) &
                * k_spring(image_id-1) / k_spring0(image_id-1)
           fz_spring_perp(:) = ( fz_spring(:) - ( dot_product(fx_spring, tau(:,1)) &
                + dot_product(fy_spring, tau(:,2)) &
                + dot_product(fz_spring, tau(:,3)) ) * tau(:,3) ) &
                * k_spring(image_id-1) / k_spring0(image_id-1)
           cosphi = ( dot_product( rx2 - rx0, rx0 - rx1 ) &
                + dot_product( ry2 - ry0, ry0 - ry1 ) &
                + dot_product( rz2 - rz0, rz0 - rz1 ) ) &
                / ( distance20 * distance01 )

           if( cosphi > 0.0d0 ) then
              switch =  0.5d0 * ( 1.0d0 + cos(PAI*cosphi) )
           else
              switch = 1.0d0
           end if
        end if

        ! for debug
        f_spring_para_sum = 0.0d0
        f_spring_perp_sum = 0.0d0
        f_true_perp_sum = 0.0d0
        do i = 1, num_atom
           if( update_id(i) ==  1 ) then
              f_spring_para_sum = f_spring_para_sum + fx_spring_para(i) * fx_spring_para(i) &
                   + fy_spring_para(i) * fy_spring_para(i) + fz_spring_para(i) * fz_spring_para(i)
              f_spring_perp_sum = f_spring_perp_sum + fx_spring_perp(i) * fx_spring_perp(i) &
                   + fy_spring_perp(i) * fy_spring_perp(i) + fz_spring_perp(i) * fz_spring_perp(i)
              f_true_perp_sum = f_true_perp_sum + fx0_perp(i) * fx0_perp(i) &
                   + fy0_perp(i) * fy0_perp(i) + fz0_perp(i) * fz0_perp(i)
           end if
        end do

        ! get constained force

        fx1(:) = fx_spring_para(:) + fx0_perp(:) + switch * fx_spring_perp(:)
        fy1(:) = fy_spring_para(:) + fy0_perp(:) + switch * fy_spring_perp(:)
        fz1(:) = fz_spring_para(:) + fz0_perp(:) + switch * fz_spring_perp(:)

        do i = 1, num_atom
           neb%image(image_id)%spring_force(i,1) = fx_spring_para(i)
           neb%image(image_id)%spring_force(i,2) = fy_spring_para(i)
           neb%image(image_id)%spring_force(i,3) = fz_spring_para(i)

           neb%image(image_id)%true_force(i,1) = fx0_perp(i)
           neb%image(image_id)%true_force(i,2) = fy0_perp(i)
           neb%image(image_id)%true_force(i,3) = fz0_perp(i)
        end do

        f0_perp_sum = 0.0d0
        f_sum = 0.0d0
        do i = 1, num_atom
           if( update_id(i) == 1  ) then
              f0_perp_sum = f0_perp_sum &
                   +  ( fx0_perp(i)+switch*fx_spring_perp(i) ) ** 2 &
                   +  ( fy0_perp(i)+switch*fy_spring_perp(i) ) ** 2 &
                   +  ( fz0_perp(i)+switch*fz_spring_perp(i) ) ** 2
              f_sum = f_sum +  fx1(i)*fx1(i) +  fy1(i)*fy1(i) + fz1(i)*fz1(i)
           end if
        end do

        f_sum = sqrt( f_sum )
        f_spring_para_sum = sqrt( f_spring_para_sum )
        f_spring_perp_sum = sqrt( f_spring_perp_sum )
        f_true_perp_sum = sqrt( f_true_perp_sum )
        f0_perp_sum = sqrt( f0_perp_sum )


        energy = neb%image(image_id)%energy

        if ( mype==0 ) then
           write(nfneb,*) 'image: ', image_id
           write(nfneb,*) ' energy: ', neb%step, energy
           write(nfneb,'(a31,e20.12)') ' sum of force (neb)         : ', f_sum
           write(nfneb,'(a31,e20.12)') ' sum of spring force (para) : ', f_spring_para_sum
           write(nfneb,'(a31,e20.12)') ' sum of true force (perp)   : ', f0_perp_sum
           write(nfneb,'(a31,3e11.3)') ' max force                  : ', &
                maxval( abs( fx1 ) ), maxval( abs( fy1 ) ), maxval( abs( fz1 ) )
        endif
     else

        ! get true force along local tangent

        dot_product_f0_tau &
             = dot_product(fx0(:), tau(:,1)) + dot_product(fy0(:), tau(:,2))&
             + dot_product(fz0(:), tau(:,3))

        fx0_para(:) = dot_product_f0_tau * tau(:,1)
        fy0_para(:) = dot_product_f0_tau * tau(:,2)
        fz0_para(:) = dot_product_f0_tau * tau(:,3)

        ! get true force perpendicular to local tangent

        fx0_perp(:) = fx0(:) - dot_product_f0_tau * tau(:,1)
        fy0_perp(:) = fy0(:) - dot_product_f0_tau * tau(:,2)
        fz0_perp(:) = fz0(:) - dot_product_f0_tau * tau(:,3)

        ! get constained force

        fx1(:) = fx0_perp(:) - fx0_para(:)
        fy1(:) = fy0_perp(:) - fy0_para(:)
        fz1(:) = fz0_perp(:) - fz0_para(:)

        do i = 1, num_atom
           neb%image(image_id)%spring_force(i,1) = fx0_para(i)
           neb%image(image_id)%spring_force(i,2) = fy0_para(i)
           neb%image(image_id)%spring_force(i,3) = fz0_para(i)

           neb%image(image_id)%true_force(i,1) = fx0_perp(i)
           neb%image(image_id)%true_force(i,2) = fy0_perp(i)
           neb%image(image_id)%true_force(i,3) = fz0_perp(i)
        end do

        f0_para_sum = 0.0d0
        f0_perp_sum = 0.0d0
        f_sum = 0.0d0
        do i = 1, num_atom
           if( update_id(i) == 1 ) then
              f0_para_sum = f0_para_sum +  fx0_para(i) * fx0_para(i) &
                   + fy0_para(i) * fy0_para(i) + fz0_para(i) * fz0_para(i)
              f0_perp_sum = f0_perp_sum +  fx0_perp(i) * fx0_perp(i) &
                   + fy0_perp(i) * fy0_perp(i) + fz0_perp(i) * fz0_perp(i)
              f_sum = fx1(i) * fx1(i) + fy1(i) * fy1(i) + fz1(i) * fz1(i)
           end if
        end do

        f_sum = sqrt( f_sum )
        f0_para_sum = sqrt( f0_para_sum )
        f0_perp_sum = sqrt( f0_perp_sum )

        energy = neb%image(image_id)%energy

        if ( mype==0 ) then
           write(nfneb,*) 'image: ', image_id
           write(nfneb,*) ' energy: ', neb%step, energy
           write(nfneb,'(a31,e20.12)') ' sum of force (neb)         : ', f_sum
           write(nfneb,'(a31,e20.12)') ' sum of true force (para)   : ', f0_para_sum
           write(nfneb,'(a31,e20.12)') ' sum of true force (perp)   : ', f0_perp_sum
           write(nfneb,'(a31,3e11.3)') ' max force                  : ', &
                maxval( abs( fx1 ) ), maxval( abs( fy1 ) ), maxval( abs( fz1 ) )
        endif
     end if

     ! set force to module

     do i = 1, num_atom
        neb%image(image_id)%force(i,1) = fx1(i)
        neb%image(image_id)%force(i,2) = fy1(i)
        neb%image(image_id)%force(i,3) = fz1(i)
     end do

     deallocate( tau )
     deallocate( update_id )
     deallocate( rx0, ry0, rz0, rx1, ry1, rz1, rx2, ry2, rz2, &
          fx0, fy0, fz0, fx1, fy1, fz1, fx_spring, fy_spring, fz_spring, &
          fx0_para, fy0_para, fz0_para, fx0_perp, fy0_perp, fz0_perp, &
          fx_spring_para, fy_spring_para, fz_spring_para, &
          fx_spring_perp, fy_spring_perp, fz_spring_perp)

  end do

  deallocate( k_spring, k_spring0 )

  return

end subroutine neb_force

function get_energy_sum() result(energy)
  use m_Replica
  implicit none
  integer :: i
  real(kind=DP) :: energy
  energy = 0.d0
  do i=2,neb%number_of_images-1
     energy = energy + neb%image(i)%energy
  enddo
  return
end function get_energy_sum

subroutine atom_update_neb

  use m_Control_Parameters, only : kqnmditer_p
  use m_Crystal_Structure,  only : altv, rltv
  use m_Ionic_System,  only : ionic_mass,m_IS_cg2_core,m_IS_fire_core,m_IS_gdiis_alloc &
  &   , m_IS_rot_ncrspd, m_IS_stor_cps_forc, iter_gdiis, ncrspd, m_IS_do_bfgs
  use m_Files,                only : nfneb
  use m_Const_Parameters,      only : PAI2

  use m_Replica
  use m_Parallelization,  only : mype

  implicit none

  integer i, j, k
  real(8) mass
  real(8) dt
  real(8) energy,get_energy_sum
  real(8) r(3), r0(3), r1(3), v(3), v1(3), f(3), f0(3) , unitf(3)
  real(8), allocatable ::  rltv_t(:,:)
  real(kind=DP), save, allocatable, dimension(:) :: x,g,diag,w
  real(kind=DP), allocatable, dimension(:,:) :: cpst,forct,cpdt
  integer, allocatable, dimension(:,:) :: imdtypt
  real(kind=DP) :: ff,eps,xtol,gtol,stpmin,stpmax,sdfac
  integer :: ndim,msave,nwork,lp,mp
  logical :: diagco
  integer :: icount
  integer :: nsum,m,it
  integer :: iprint(2)
  integer, save :: if_allocated=0
  integer, save :: iflag=0
  EXTERNAL LB2
  COMMON /LB3/MP,LP,GTOL,STPMIN,STPMAX

  if (mype==0) write(nfneb,*) 'atom coordinate update...'

  dt = neb%dt
  allocate(rltv_t(3,3))
  rltv_t = transpose(rltv)/PAI2

  energy = get_energy_sum()
#ifdef _L_BFGS_
  if(neb%cond%time_integral == 'lbfgs') then
    msave = 5
    icount = 0
    do i=2,neb%number_of_images-1
       do j=1,neb%image(i)%num_atom
          if(neb%image(i)%fix_flag(j) == 1) then
             icount = icount+3
          endif
       enddo
    enddo
    ff = energy
    ndim = icount
    nwork = ndim*(2*msave+1)+2*msave
    if(.not.allocated(x)) then
      allocate(x(ndim));x=0.d0
      allocate(g(ndim));g=0.d0
      allocate(diag(ndim));diag=0.d0
      allocate(w(nwork));w=0.d0
    endif
    icount = 0
    do i=2,neb%number_of_images-1
       do j=1,neb%image(i)%num_atom
          if(neb%image(i)%fix_flag(j) == 1)then
             do k=1,3
                icount = icount+1
                x(icount) = neb%image(i)%cps(j,k)
                g(icount) = -neb%image(i)%force(j,k)
             enddo
          endif
       enddo
    enddo
    iprint(1) = 3
    iprint(2) = 3
    call lbfgs(ndim,msave,x,ff,g,.false.,diag,IPRINT,0.d0,1e-16,w,iflag)
    icount = 0
    do i=2,neb%number_of_images-1
       do j=1,neb%image(i)%num_atom
          if(neb%image(i)%fix_flag(j) == 1)then
             do k=1,3
                icount = icount+1
                neb%image(i)%cps(j,k) = x(icount)
             enddo
          endif
       enddo
    enddo
    if (iflag<0) then
      write(nfneb,*) 'encountered errors in lbfgs : ',iflag
    endif
  else if (neb%cond%time_integral == 'fire')then
#else
  if (neb%cond%time_integral == 'fire')then
#endif
    ndim = (neb%number_of_images-2)*neb%image(1)%num_atom
    allocate(cpst(ndim,3))
    allocate(cpdt(ndim,3))
    allocate(forct(ndim,3))
    allocate(imdtypt(ndim,3))
    icount = 0
    do i=2,neb%number_of_images-1
       do j=1,neb%image(i)%num_atom
          icount = icount + 1
          do k=1,3
            cpst(icount,k) = neb%image(i)%cps(j,k)
            cpdt(icount,k) = neb%image(i)%velocity(j,k)
            forct(icount,k) = neb%image(i)%force(j,k)
            imdtypt(icount,k) = neb%image(i)%fix_flag(j)
          enddo
       enddo
    enddo
    call m_IS_fire_core(nfneb,ndim,neb%step,cpst,cpdt,forct,imdtypt)
    icount = 0
    do i=2,neb%number_of_images-1
       do j=1,neb%image(i)%num_atom
          icount = icount+1
          neb%image(i)%cps(j,:) = cpst(icount,:)
          neb%image(i)%velocity(j,:) = cpdt(icount,:)
       enddo
    enddo
    deallocate(cpst)
    deallocate(cpdt)
    deallocate(forct)
    deallocate(imdtypt)
  else if (neb%cond%time_integral == 'bfgs')then
    ndim = (neb%number_of_images-2)*neb%image(1)%num_atom
    iter_gdiis = iter_gdiis + 1
    if(if_allocated==0) call m_IS_gdiis_alloc(1,if_allocated,ndim)
    m = mod(iter_gdiis,kqnmditer_p)
    if(m == 0) m = kqnmditer_p
    call m_IS_rot_ncrspd(nsum)  ! -> ncrspd, nsum
    it = ncrspd(nsum)
    if(it == 0) it = 1
    allocate(cpst(ndim,3))
    allocate(forct(ndim,3))
    allocate(imdtypt(ndim,3))
    icount = 0
    do i=2,neb%number_of_images-1
       do j=1,neb%image(i)%num_atom
          icount = icount + 1
          do k=1,3
            cpst(icount,k) = neb%image(i)%cps(j,k)
            forct(icount,k) = neb%image(i)%force(j,k)
            imdtypt(icount,k) = neb%image(i)%fix_flag(j)
          enddo
       enddo
    enddo
    call m_IS_stor_cps_forc(it,ndim,cpst,forct) ! cps,forc_g -> u_l, w_l
    call m_IS_do_bfgs(nfneb,nsum,ndim,get_energy_sum(),imdtypt,cpst,forct)
    icount = 0
    do i=2,neb%number_of_images-1
       do j=1,neb%image(i)%num_atom
          icount = icount+1
          neb%image(i)%cps(j,:) = cpst(icount,:)
       enddo
    enddo
    deallocate(cpst)
    deallocate(forct)
    deallocate(imdtypt)
  else if (neb%cond%time_integral == 'cg')then
    ndim = (neb%number_of_images-2)*neb%image(1)%num_atom
    allocate(cpst(ndim,3))
    allocate(forct(ndim,3))
    allocate(imdtypt(ndim,3))
    icount = 0
    do i=2,neb%number_of_images-1
       do j=1,neb%image(i)%num_atom
          icount = icount + 1
          do k=1,3
            cpst(icount,k) = neb%image(i)%cps(j,k)
            forct(icount,k) = neb%image(i)%force(j,k)
            imdtypt(icount,k) = neb%image(i)%fix_flag(j)
          enddo
       enddo
    enddo
    energy = get_energy_sum()
    call m_IS_cg2_core(nfneb,ndim,cpst,forct,imdtypt,energy)
    icount = 0
    do i=2,neb%number_of_images-1
       do j=1,neb%image(i)%num_atom
          icount = icount+1
          neb%image(i)%cps(j,:) = cpst(icount,:)
       enddo
    enddo
    deallocate(cpst)
    deallocate(forct)
    deallocate(imdtypt)
  else
    do i = 2, neb%number_of_images-1

       neb%image(i)%pos0 = neb%image(i)%pos
       neb%image(i)%cps0 = neb%image(i)%cps

       do j = 1, neb%image(i)%num_atom

          if(neb%image(i)%fix_flag(j) == 1) then

             r = neb%image(i)%cps(j,:)
             r0 = neb%image(i)%cps0(j,:)
             v = neb%image(i)%velocity(j,:)
             f = neb%image(i)%force(j,:)
             f0 = neb%image(i)%force0(j,:)
             energy = neb%image(i)%energy
             mass = ionic_mass(j)

             select case(neb%cond%time_integral)
             case('steepest_descent')
                r1 = r + dt*f/neb%cond%sd_factor
                v = 0.0d0
             case('verlet')
                r1 = 2.0d0*r - r0 + dt*dt*f/mass
                v1 = (r1-r)/(2.0d0*dt)
                v = 0.0d0
                if( dot_product(v1,f) >= 0.0d0) then
                   if( dot_product(f,f) /= 0.0d0 ) then
                      v = f * dot_product(v1,f)/dot_product(f,f)
                   end if
                end if
             case('velocity_verlet')
                v1 = v + 0.5d0*dt*(f0+f)/mass
                v = 0.0d0
                unitf=0.d0
                if(dot_product(f,f)/=0) unitf=f/dsqrt(dot_product(f,f))
                if( dot_product(v1,unitf) >= 0.0d0) then
                  v = unitf * dot_product(v1,unitf)
                else
                   if (mype==0) write(nfneb,'(a,i6)') 'quenched velocity of atom : ',j
                end if
                r1 = r + v*dt + f*dt**2/mass
             end select

             neb%image(i)%cps(j,:) = r1
             neb%image(i)%velocity(j,:) = v

          end if

        ! periodic tranformation
        !do k = 1, 3
        !   if( neb%image(i)%pos(j,k) < 0 ) then
        !      neb%image(i)%pos(j,k) = neb%image(i)%pos(j,k) + 1.d0
        !   else if( neb%image(i)%pos(j,k) > 1.d0 ) then
        !      neb%image(i)%pos(j,k) = neb%image(i)%pos(j,k) - 1.d0
        !   end if
        !enddo

        ! pos -> cps
        !do k = 1, 3
        !   neb%image(i)%cps(j,k) = altv(k,1)*neb%image(i)%pos(j,1) &
        !        + altv(k,2)*neb%image(i)%pos(j,2) &
        !        + altv(k,3)*neb%image(i)%pos(j,3)
        !enddo
       end do
    end do
  endif
  do i = 2, neb%number_of_images-1
    do j = 1, neb%image(i)%num_atom
      do k = 1, 3
         neb%image(i)%pos(j,k) = rltv_t(k,1)*neb%image(i)%cps(j,1) &
              + rltv_t(k,2)*neb%image(i)%cps(j,2) &
              + rltv_t(k,3)*neb%image(i)%cps(j,3)
      enddo
    enddo
    if ( mype==0 ) then
      write(nfneb,*) ' image: ', i
      write(nfneb,*) '  cps, cps0, cps-cps0 '
      do j = 1, neb%image(i)%num_atom
         write(nfneb,'(6f10.5,3e11.3)') neb%image(i)%cps(j,1:3), neb%image(i)%cps0(j,1:3), &
              neb%image(i)%cps(j,1:3)-neb%image(i)%cps0(j,1:3)
      end do
    endif
  enddo

  deallocate(rltv_t)

end subroutine atom_update_neb

subroutine atom_update_pos_neb

  use m_Crystal_Structure,  only : altv, rltv
  use m_Ionic_System,  only : ionic_mass
  use m_Files,                only : nfneb
  use m_Const_Parameters,      only : PAI2
  use m_Parallelization,     only : mype

  use m_Replica

  implicit none

  integer i, j, k
  real(8) mass
  real(8) dt
  real(8) r(3), r0(3), r1(3), v(3), v1(3), f(3), f0(3)

  if (mype==0) write(nfneb,*) 'atom coordinate update... (pos)'

  dt = neb%dt

  do i = 2, neb%number_of_images-1

     neb%image(i)%pos0 = neb%image(i)%pos
     neb%image(i)%cps0 = neb%image(i)%cps

     do j = 1, neb%image(i)%num_atom

        if(neb%image(i)%fix_flag(j) == 1) then

           r = neb%image(i)%pos(j,:)
           r0 = neb%image(i)%pos0(j,:)
           v = neb%image(i)%velocity(j,:)
           f = neb%image(i)%force(j,:)
           f0 = neb%image(i)%force0(j,:)
           mass = ionic_mass(j)

           select case(neb%cond%time_integral)
           case('steepest_descent')
              r1 = r + dt*f
              v = 0.0d0
           case('verlet')
              r1 = 2.0d0*r - r0 + dt*dt*f/mass
              v1 = (r1-r)/(2.0d0*dt)
              v = 0.0d0
              if( dot_product(v1,f) >= 0.0d0) then
                 if( dot_product(f,f) /= 0.0d0 ) then
                    v = f * dot_product(v1,f)/dot_product(f,f)
                 end if
              end if
           case('velocity_verlet')
              v1 = v + 0.5d0*dt*(f0+f)/mass
              v = 0.0d0
              if( dot_product(v1,f) >= 0.0d0) then
                 if( dot_product(f,f) /= 0.0d0 ) then
                    v = f * dot_product(v1,f)/dot_product(f,f)
                 end if
              end if
              r1 = r + v*dt + f*dt**2/mass
           end select

           neb%image(i)%pos(j,:) = r1
           neb%image(i)%velocity(j,:) = v

        end if

        ! periodic tranformation
        !do k = 1, 3
        !   if( neb%image(i)%pos(j,k) < 0 ) then
        !      neb%image(i)%pos(j,k) = neb%image(i)%pos(j,k) + 1.d0
        !   else if( neb%image(i)%pos(j,k) > 1.d0 ) then
        !      neb%image(i)%pos(j,k) = neb%image(i)%pos(j,k) - 1.d0
        !   end if
        !enddo

        ! pos -> cps
        do k = 1, 3
           neb%image(i)%cps(j,k) = altv(k,1)*neb%image(i)%pos(j,1) &
                + altv(k,2)*neb%image(i)%pos(j,2) &
                + altv(k,3)*neb%image(i)%pos(j,3)
        enddo

     end do

     if ( mype==0 ) then
        write(nfneb,*) ' image: ', i
        write(nfneb,*) '  cps, cps0, cps-cps0, velocity '
        do j = 1, neb%image(i)%num_atom
           write(nfneb,'(6f10.5,3e11.3,3f10.5)') neb%image(i)%cps(j,1:3), neb%image(i)%cps0(j,1:3), &
                neb%image(i)%cps(j,1:3)-neb%image(i)%cps0(j,1:3), &
                neb%image(i)%velocity(j,1:3)
        end do
     endif

  end do

end subroutine atom_update_pos_neb

subroutine write_down_data_replica(replace_nfstop)

  use m_Files,                only : nfneb, nfnebcntn, F_NEB_CNTN, nfstop, F_STOP
  use m_Replica
  use m_Parallelization,    only : mype
  use m_Control_Parameters, only : sw_fcp
  use m_Const_Parameters,   only : ON
  use m_Fcp,                only : m_Fcp_write_down_totch_replica

  implicit none
  logical, intent(in) :: replace_nfstop
  integer i,j
  if (mype==0) write(nfneb,*) 'write down data replica...'

!  if(mype_r == 0) then
  if(mype == 0) then

  open(nfnebcntn,file=F_NEB_CNTN)
  write(nfnebcntn,*) 'iteration_neb'
  write(nfnebcntn,*) neb%step
  write(nfnebcntn,*) 'number_of_images'
  write(nfnebcntn,*) neb%number_of_images

  do i = 1, neb%number_of_images
     write(nfnebcntn,*) 'image'
     write(nfnebcntn,*) i
     write(nfnebcntn,*) '(natm)'
     write(nfnebcntn,*) neb%image(i)%num_atom
     write(nfnebcntn,*) '(pos)'
     do j = 1, neb%image(i)%num_atom
        write(nfnebcntn,'(3d24.16)') neb%image(i)%pos(j,1:3)
     end do
     write(nfnebcntn,*) '(cps)'
     do j = 1, neb%image(i)%num_atom
        write(nfnebcntn,'(3d24.16)') neb%image(i)%cps(j,1:3)
     end do
     write(nfnebcntn,*) '(energy)'
     write(nfnebcntn,'(d24.16)') neb%image(i)%energy
     write(nfnebcntn,*) '(force_org)'
     do j = 1, neb%image(i)%num_atom
        write(nfnebcntn,'(3d24.16)') neb%image(i)%force_org(j,1:3)
     end do
     write(nfnebcntn,*) '(convergence)'
     if(neb%image(i)%scf_convergence) then
        write(nfnebcntn,'(i3)')  1
     else
        write(nfnebcntn,'(i3)')  0
     end if
     if( sw_fcp == ON ) call m_Fcp_write_down_totch_replica(i)
  end do
  write(nfnebcntn,*) '(changed_to_cineb)'
  write(nfnebcntn,'(l1)') changed_to_cineb

  close(nfnebcntn)

  if(replace_nfstop)then
    open(nfstop, file=F_STOP, status='replace', form='formatted')
    close(nfstop, status='keep')
  endif

  end if
end subroutine write_down_data_replica

subroutine allreduce_scf_convergence

  use m_Files,   only : nfneb
  use m_Replica
  use m_Parallelization,  only : mype
  use mpi

  implicit none

  integer i, mpi_err
  logical,allocatable :: work(:)

!  include 'mpif.h'

  if(nrank_r.ge.1)then
    allocate(work(neb%number_of_images))
    call mpi_allreduce(neb%image(:)%scf_convergence,work,neb%number_of_images, &
         mpi_logical,mpi_land,mpi_comm_world,mpi_err)
    neb%image(:)%scf_convergence = work
    deallocate(work)
  endif

  if ( mype==0 ) then
     write(nfneb,*) 'neb: ', neb%step
     do i=1, neb%number_of_images
        write(nfneb,*) ' image: ', i, ' scf convergence: ', neb%image(i)%scf_convergence
     end do
  endif
end subroutine allreduce_scf_convergence

logical function replica_converged()

  use m_Files,                only : nfneb, nfnebstop, F_NEB_STOP
  use m_Replica
  use m_Parallelization, only : mype
  use m_Control_Parameters, only : sw_fcp
  use m_Const_Parameters, only : ON
  use m_Fcp, only : m_Fcp_replica_Converged

  implicit none

  logical conv
  real(8), allocatable :: force(:)
  integer i, j, istop, id(1)
  real(8) :: fmax, ftmp

  if ( mype==0 ) write(nfneb,*) 'convergence check...', neb%step

  conv = .false.
  if (neb%cond%climbing_image .and. .not. changed_to_cineb) then
    if (mype==0) write(nfneb,*) &
    &  'we have not transitioned to CI-NEB yet, so convergence test will be ignored'
  endif
  select case(neb%cond%convergence_condition)

  case(1) ! energy
     curr_thres = maxval(abs(neb%image(:)%energy-neb%image(:)%energy0))
     if( curr_thres < neb%cond%convergence_threshold ) then
        conv = .true.
     end if
     if ( mype==0) then
        write(nfneb,*) ' energy : ', maxval(abs(neb%image(:)%energy-neb%image(:)%energy0)), &
             neb%cond%convergence_threshold, conv
     endif
  case(2) ! force(org)
     fmax = 0
     do i=2,neb%number_of_images-1
       ftmp = get_max_force(neb,i,'org')
       if(ftmp<fmax) fmax=ftmp
     enddo
     if(fmax<neb%cond%convergence_threshold) conv = .true.
     if ( mype==0 ) then
        write(nfneb,*) ' force (org): ', fmax, neb%cond%convergence_threshold, conv
     endif
     curr_thres = fmax
  case(3) ! force(neb)
     fmax = 0
     do i=2, neb%number_of_images-1
       ftmp = get_max_force(neb,i,'neb')
       if(fmax < ftmp) fmax=ftmp
     enddo
     if( fmax < neb%cond%convergence_threshold ) conv = .true.
     if ( mype==0 ) then
        write(nfneb,*) ' force (neb): ', fmax, neb%cond%convergence_threshold, conv
     endif
     curr_thres = fmax
  case(4) ! force(org) of max energy image
     id = maxloc(neb%image(:)%energy);
     fmax = 0
     do i=1,neb%image(id(1))%num_atom
       if(neb%image(id(1))%fix_flag(i)==0)cycle
       ftmp = dsqrt(dot_product(neb%image(id(1))%force(i,:),neb%image(id(1))%force(i,:)))
       if(fmax<ftmp)fmax=ftmp
     enddo
     if(fmax<neb%cond%convergence_threshold) conv=.true.
     if ( mype==0 ) then
        write(nfneb,*) ' force (transition state): ', fmax, &
             neb%cond%convergence_threshold, conv
     endif
     curr_thres = fmax
  case(5) ! true force (normal force)
     fmax=0
     do i=2,neb%number_of_images-1
       ftmp = get_max_force(neb,i,'true')
       if(fmax<ftmp)fmax=ftmp
     enddo
     if(fmax<neb%cond%convergence_threshold) conv=.true.
     if ( mype==0 ) then
        write(nfneb,*) ' force (normal): ', fmax, neb%cond%convergence_threshold, conv
     endif
     curr_thres = fmax
  end select

  if (neb%cond%climbing_image .and. .not. changed_to_cineb) then
    conv = .false.
  endif
  ! check nfnebstop.data
  istop = -1
  open(nfnebstop,file=F_NEB_STOP)
  rewind(nfnebstop)
  read(nfnebstop,*,end=1,err=1) istop
1 continue
  close(nfnebstop,status='keep')
  if(0 <= istop .and. neb%step >= istop) then
     conv = .true.
     if (mype==0) write(nfneb,*) ' nfnebstop.data: ', istop
  end if

  if(sw_fcp == ON) then
    if(m_Fcp_replica_Converged()) conv = .false.
  end if

  replica_converged = conv

  call write_result(.true.)

end function replica_converged

subroutine write_result(inc_energy)

  use m_Files,                only : nfnebenf, nfnebdynm, nfenf, nfdynm, F_DYNM, F_ENF
  use m_Replica
  use m_Ionic_System,         only : m_IS_wd_speciesname_etc
  use m_IterationNumbers,     only : iteration
  use mpi

  implicit none

  logical, intent(in) :: inc_energy
  integer i, j
  character(20) form
  real(8) d, sum_d
  integer npes,mype,mpi_err
  real(8) :: mforg,mf,mftrue

!  include 'mpif.h'

  call mpi_comm_rank(mpi_comm_world,mype,mpi_err)

  if(mype /= 0) return

  ! write energy, force, cps

  open(nfdynm,file=F_DYNM)
  open(nfenf,file=F_ENF)
  call m_IS_wd_speciesname_etc(nfdynm)
  call m_IS_wd_speciesname_etc(nfnebdynm)
  write(nfenf,*) '#step  image  image_distance  energy   force_org  force_neb  force_normal'
  sum_d = 0.0d0
  do i=1, neb%number_of_images
     if(inc_energy)then
     d = 0.0d0
     if(i>1) then
        d = 0.0d0
        do j=1, neb%image(i)%num_atom
           d = d + dot_product(neb%image(i-1)%cps(j,:)-neb%image(i)%cps(j,:), &
                neb%image(i-1)%cps(j,:)-neb%image(i)%cps(j,:))
        end do
        d = sqrt(d)
     end if
     sum_d = sum_d + d
     write(nfnebenf,'(2i5,5e20.10)')  neb%step, i, sum_d, neb%image(i)%energy, &
          get_max_force(neb,i,'org'),get_max_force(neb,i,'neb'),get_max_force(neb,i,'true')
     write(nfenf,'(2i5,5e20.10)')  neb%step, i, sum_d, neb%image(i)%energy, &
          get_max_force(neb,i,'org'),get_max_force(neb,i,'neb'),get_max_force(neb,i,'true')
     endif
     write(nfnebdynm,'(a,i5,i8,a)') ' cps and forc at (nebstep, iter_total=',neb%step,iteration,')'
     write(nfdynm,'(a,i5,i8,a)') ' cps and forc at (nebstep, iter_total=',neb%step,iteration,')'
     do j = 1, neb%image(i)%num_atom
!        write(nfnebdynm,'(3i5,3f20.10)') neb%step, i, j, neb%image(i)%cps(j,1:3)
        write(nfnebdynm,'(" ",i4,3f15.9,3f12.6)') j &
        &, neb%image(i)%cps(j,1),neb%image(i)%cps(j,2),neb%image(i)%cps(j,3) &
        &, neb%image(i)%force(j,1), neb%image(i)%force(j,2), neb%image(i)%force(j,3)
        write(nfdynm,'(" ",i4,3f15.9,3f12.6)') j &
        &, neb%image(i)%cps(j,1),neb%image(i)%cps(j,2),neb%image(i)%cps(j,3) &
        &, neb%image(i)%force(j,1), neb%image(i)%force(j,2), neb%image(i)%force(j,3)
     end do
  end do

  if(inc_energy) write(nfnebenf,*)

  call flush(nfnebdynm)
  call flush(nfnebenf)

  close(nfdynm)
  close(nfenf)
  ! temporary
  !!write(form,'(a4,i5,a7)') '(i5,',neb%number_of_images,'e20.10)'
  !!write(212,form) neb%step, neb%image(:)%energy

  !!neb%energy(neb%step,:) = neb%image(:)%energy
  !!open(213)
  !!write(form,'(a4,i5,a7)') '(i5,',neb%step,'e20.10)'
  !!do i=1, neb%number_of_images
  !!   write(213,form) i, neb%energy(1:neb%step,i)
  !!end do
  !!close(213)

end subroutine write_result

subroutine finalize_neb()
  use m_Replica, only : m_Replica_finalize, nrank_r
  use m_Files, only : nfnebstop, F_NEB_STOP, nfneb, m_Files_close_and_clear_nfstop
  use m_Control_Parameters,  only : ipriparadeb
  use m_Parallelization, only : mype
  use mpi

  implicit none
!  include 'mpif.h'

  integer :: mpi_err
  call m_Replica_finalize()

  if(nrank_r.gt.1) call mpi_barrier(mpi_comm_world, mpi_err)
  open(nfnebstop, file=F_NEB_STOP, status='replace', form='formatted')
  close(nfnebstop, status='keep')

! ==== KT_add ===== 2014/07/20
#ifdef NEB_NEW_FILENAMES
  if ( mype == 0 ) close(nfneb)
#else
  if ( mype == 0 ) then
     close(nfneb)
  else if(ipriparadeb == 0) then
     close(nfneb, status='delete')
  else
     close(nfneb, status='keep')
  end if
#endif
! ================= 2014/07/20

  call m_Files_close_and_clear_nfstop()
end subroutine finalize_neb

subroutine deallocate_array_partial
  use m_PlaneWaveBasisSet
  use m_Parallelization
  use m_Kpoints
  implicit none
  call m_pwBS_dealloc()
  call m_Parallel_dealloc(neb_mode=.true.)
  call m_Kp_dealloc
end subroutine deallocate_array_partial

subroutine deallocate_array

  use m_Control_Parameters
  use m_Crystal_Structure
  use m_Ionic_System
  use m_PlaneWaveBasisSet
  use m_Parallelization
  use m_Kpoints
  use m_Force
  use m_Charge_Density
  use m_XC_Potential
  use m_PseudoPotential
  use m_NonLocal_Potential
  use m_Electronic_Structure
  use m_ES_WF_by_SDorCG

!  call m_CtrlP_dealloc

  call m_CS_dealloc(neb_mode=.true.)
  call m_IS_dealloc(neb_mode=.true.)
  call m_pwBS_dealloc
  call m_Parallel_dealloc(neb_mode=.true.)
  call m_Parallel_dealloc_mpi_elec
  call m_Parallel_dealloc_mpi_fft_box
  call m_Parallel_dealloc_mpi_nval
  call m_Parallel_fft_onto_wf_dealloc_3D
  call m_Parallel_dealloc_mpi_kngp_B
  call m_Kp_dealloc
  call m_Force_dealloc
  call m_CD_dealloc
  call m_XC_dealloc_vxc_3D
  !call m_PP_dealloc
  !call m_NLP_dealloc
  call m_ES_dealloc
  call m_ESsd_dealloc
  call m_CtrlP_set_init_status(.true.)

end subroutine deallocate_array
