!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 606 $)
!
!  PROGRAM: EK
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!      Further modification by T. Yamasaki   Feb. 2004
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
! $Id: ekmain.f90 606 2020-04-15 06:45:49Z ktagami $
!
program EK
! This program was coded by T. Yamasaki(FUJITSU Laboratories Ltd.), 17th Feb. 2003.
!
  implicit none
  logical  :: EigenValues_are_Converged, AllKpoints_are_Calculated
  logical  :: Already_Converged
  logical  :: Ending_Time

  call Initialization_set_ekmode_ON  ! set `ekmode' ON in m_Control_Parameters
  call Initialization(1)             ! Initialization of mpi and file-setting
  call InputData_Analysis
  call Preparation(0)                ! Basis set, symmetry check etc.
  call Preparation_for_mpi(1)        ! mpi
  call PseudoPotential_Construction
#ifdef ENABLE_ESM_PACK
  call Preparation_for_ESM
#endif
  call Ewald_and_Structure_Factor    ! Calculate Structure Factor
  call Initial_Electronic_Structure()! read Charge Density, (lclchh)

  KPOINTS: do
!     call KpointNumber_Setting()
     call KpointNumber_Setting2()
     call Preparation_ek             ! (basnum)
     call Preparation_for_mpi_ek     ! mpi  -> np_g1k, mp_g1k
     call PseudoPotential_ek         ! (kbint)
     call Initial_WaveFunctions_ek   ! (rndzaj|rdzaj),(fsrfsi),(lclchh)
     if(.not.Already_Converged()) then
        SolveWaveFunctions: do
           if(Ending_Time())                 exit KPOINTS
           call IterationNumber_Setting()
           call Renewal_of_WaveFunctions()
           if(EigenValues_are_Converged())   exit SolveWaveFunctions
        enddo SolveWaveFunctions
        call Postprocessing_k()
        if(AllKpoints_are_Calculated())     exit KPOINTS
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
  call WriteDownData_onto_Files_ek()
  call Finalization_of_mpi           ! mpi
end program EK
