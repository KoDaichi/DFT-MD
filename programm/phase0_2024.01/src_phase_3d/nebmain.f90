!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 599 $)
!
!  MAIN PROGRAM: CPMD
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!
!  FURTHER MODIFICATION: T. Yamasaki, January/13/2004
!
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)
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
!
!  $Id: nebmain.f90 599 2019-12-03 05:25:28Z ktagami $
!
program PHASENEB

  use m_Files,                only : nfneb
  use m_Const_Parameters, only : OFF,ON
  use m_IterationNumbers, only : iteration,iteration_electronic,iteration_ionic
  use m_Replica
  use m_Parallelization, only : MPI_CommGroup
  use mpi

  implicit none
  logical  :: ChargeDensity_is_Converged, TotalEnergy_is_Divergent
  logical  :: Already_Converged, Positron_bulk
  logical  :: Hubbard_model
  logical  :: Forces_are_Converged, Ending_Time, Ending_Time2, Force_errors_are_tolerable
  logical  :: replica_converged

  integer i, j, itr
  integer mpi_err
  logical all_conv_flag

  integer prepare_communicator

!  include 'mpif.h'

  call initialize_neb

  call Initialization(OFF)
  call InputData_Analysis(OFF) ! rd_cntn_data = ON

  call set_neb_condition
  call create_replica

  itr = 0
  all_conv_flag = all(neb%image(:)%scf_convergence)  ! for restart
  prepare_communicator = ON

  NEBIteration: do

     itr = itr + 1
     if(all(neb%image(:)%scf_convergence)) then
        neb%step = neb%step + 1
     end if
     call set_neb_parameter

     MultiReplica: do i = 1, neb%number_of_images

        if(mod(i-1,nrank_r) /= mype_r) then
           neb%image(i)%scf_convergence = .true.
           cycle
        end if

        if(neb%cond%condition == 1) then
           if(itr == 1 .and. (.not. all_conv_flag)) then
              if(neb%image(i)%scf_convergence) then
                 neb%image(i)%energy = neb%image(i)%energy0
                 neb%image(i)%force = neb%image(i)%force0
                 neb%image(i)%force_org = neb%image(i)%force0
                 cycle
              end if
           end if
        end if

        if(itr > 1 .and. (i==1 .or. i==neb%number_of_images)) then
           neb%image(i)%energy = neb%image(i)%energy0
           neb%image(i)%force = neb%image(i)%force0
           neb%image(i)%force_org = neb%image(i)%force0
           cycle
        end if

!        write(nfneb,*) 'neb: ', neb%step, i, mype
        write(nfneb,*) 'neb: ', neb%step, i

        call set_neb_replica_parameter(i)
        call set_neb_replica_coordinate(i)

        if(.not.neb%image(i)%scf_convergence) then
           call InputData_Analysis_neb()
        end if

        call Preparation(0)                   ! Basis set, symmetry check etc.
        call Preparation_for_mpi(prepare_communicator)           ! mpi
        call PseudoPotential_Construction
        call Ewald_and_Structure_Factor
        call Initial_Electronic_Structure
        call Initial_MD_Condition

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
                 call ChargeDensity_Construction(0)
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

!!$        if(Positron_bulk()) then
!!$           call Renewal_of_pPotential()
!!$           call Solve_pWaveFunctions()
!!$        end if

!!$        call Postprocessing
        call WriteDownData_onto_Files

        call deallocate_array

        write(nfneb,*) 'scf convergence: ', i, neb%image(i)%scf_convergence

        prepare_communicator=OFF
     end do MultiReplica

     call allreduce_scf_convergence
     call allreduce_neb_energy_force

     if(.not. all(neb%image(:)%scf_convergence)) then
        write(nfneb,*) 'scf not convergence...'
        exit NEBIteration
     end if

!     call allreduce_neb_energy_force

     call force_neb
     call atom_update_neb
!     call atom_update_pos_neb

     if(replica_converged()) then
        exit NEBIteration
     end if
     if(neb%step >= neb%max_iteration) then
        exit NEBIteration
     end if

  end do NEBIteration

  call write_down_data_replica

  call finalize_neb()

  call mpi_finalize(mpi_err)

end program PHASENEB

!subroutine deallocate_array
!
!  use m_Control_Parameters
!  use m_Crystal_Structure
!  use m_Ionic_System
!  use m_PlaneWaveBasisSet
!  use m_Parallelization
!  use m_Kpoints
!  use m_Force
!  use m_Charge_Density
!  use m_XC_Potential
!  use m_PseudoPotential
!  use m_NonLocal_Potential
!  use m_Electronic_Structure
!  use m_ES_WF_by_SDorCG
!
!!  call m_CtrlP_dealloc
!
!  call m_CS_dealloc
!  call m_IS_dealloc
!  call m_pwBS_dealloc
!  call m_Parallel_dealloc
!  call m_Kp_dealloc
!  call m_Force_dealloc
!  call m_CD_dealloc
!  call m_XC_dealloc_vxc
!  call m_PP_dealloc
!  call m_NLP_dealloc
!  call m_ES_dealloc
!  call m_ESsd_dealloc
!
!
!end subroutine deallocate_array
