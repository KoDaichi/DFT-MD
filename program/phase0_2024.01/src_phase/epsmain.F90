!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 599 $)
!
!  PROGRAM: EPSMAIN
!
!  AUTHOR(S): T. Hamada   MAY/8/2007
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
program EPSMAIN
! $Id: epsmain.f90 599 2019-12-03 05:25:28Z ktagami $
!
  implicit none
! === DEBUG by tkato 2013/10/16 ================================================
! logical  :: EigenValues_are_Converged, AllKpoints_are_Calculated
  logical  :: EigenValues_are_Converged, AllKpoints_are_Calculated2
! ==============================================================================
  logical  :: Already_Converged
  logical  :: Ending_Time

  logical  :: all_conv = .true.
! === DEBUG by tkato 2013/10/16 ================================================
  integer :: nk = 0
! ==============================================================================
! --------------------T. Hamada 2021.9.28 --------------------
  logical  :: EPS_Prep_Check
!-------------------------------------------------------------

  call Initialization_set_ekmode_ON  ! set `ekmode' ON in m_Control_Parameters
  call Initialization(1)             ! Initialization of mpi and file-setting
  call InputData_Analysis
  call Preparation(0)                ! Basis set, symmetry check etc.
  call Preparation_for_mpi(1)        ! mpi
  call PseudoPotential_Construction
  call Ewald_and_Structure_Factor    ! Calculate Structure Factor
  call Initial_Electronic_Structure()! read Charge Density, (lclchh)
  call Initialization_Epsilon(0)     ! Epsilon
  call Shift_Kpoint                  ! for Effective mass calculation  ! Epsilon

  KPOINTS: do
! === DEBUG by tkato 2013/10/16 ================================================
!    call KpointNumber_Setting()
     call KpointNumber_Setting2()
! ==============================================================================
     call Preparation_ek             ! (basnum)
     call Preparation_for_mpi_ek     ! mpi  -> np_g1k, mp_g1k
     call PseudoPotential_ek_Epsilon ! (kbint)                         ! Epsilon
     call Initial_WaveFunctions_ek   ! (rndzaj|rdzaj),(fsrfsi),(lclchh)
     if(.not.Already_Converged()) then
        all_conv = .false.
        SolveWaveFunctions: do
!-------------------- T. Hamada 2021.9.28 --------------------
           if(EPS_Prep_Check())              exit KPOINTS
! ------------------------------------------------------------
           if(Ending_Time())                 exit KPOINTS
           call IterationNumber_Setting()
           call Renewal_of_WaveFunctions()
           if(EigenValues_are_Converged()) then                        ! Epsilon
              call Transition_moment_Epsilon                           ! Epsilon
              call Dealloc_Radr_and_Wos_Epsilon                        ! Epsilon
              exit SolveWaveFunctions
           end if
       enddo SolveWaveFunctions
       call Postprocessing_k()
! === DEBUG by tkato 2013/10/16 ================================================
!       if(AllKpoints_are_Calculated())  then
        if(AllKpoints_are_Calculated2(nk))  then
! ==============================================================================
           all_conv = .true.
           exit KPOINTS
        end if
     else
        exit KPOINTS
     end if
  enddo KPOINTS
!!$  else
!!$     write(6,'(" Already_Converged")')
!!$     call KpointNumber_Setting()
!!$     call Preparation_ek             ! (basnum)
!!$     call PseudoPotential_ek         ! (kbint)
!!$     call Initial_WaveFunctions_ek   ! (rndzaj|rdzaj),(fsrfsi),(lclchh)
!!$  end if

  call Postprocessing(.false.)
!-------------------- T. Hamada 2021.9.28 --------------------
!  if(all_conv.) then
  if(all_conv.or.EPS_Prep_Check()) then
!-------------------------------------------------------------
     call Reset_Kpoint                                                   ! Epsilon
     call Prep_for_Calc_Epsilon                                          ! Epsilon
     call Calc_Epsilon                                                   ! Epsilon
     call Calc_Nonlinear_optics                                          ! Epsilon
     call WriteDownData_onto_Files_Eps                                   ! Epsilon
  end if
  call WriteDownData_onto_Files_ek()
  call Finalization_of_mpi           ! mpi
end program EPSMAIN
