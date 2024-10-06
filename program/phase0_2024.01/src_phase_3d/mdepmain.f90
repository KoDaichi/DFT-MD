!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  MAIN PROGRAM: CPMD
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
program CPMD
! $Id: mdepmain.f90 570 2017-04-21 20:34:50Z yamasaki $
  implicit none
  logical  :: Forces_are_Converged, Ending_Time, MultiReplicaMode &
       &     ,AllForces_are_Converged

  call Initialization
  call InputData_Analysis
  MultiReplica: do
     call Preparation_ep
     call Initial_MD_Condition

     AtomicConfiguration: do
        call Forces
        if(Forces_are_Converged()) exit AtomicConfiguration
        call Move_Ions
        if(Ending_Time())          exit AtomicConfiguration
        call MDIterationNumber_Setting_ep
     enddo AtomicConfiguration

     call WriteDownData_onto_Files

     if(.not.MultiReplicaMode()) then
        exit MultiReplica
     else if( Ending_Time() .or. AllForces_are_Converged()) then
        exit MultiReplica
     end if

  end do MultiReplica

end program CPMD
