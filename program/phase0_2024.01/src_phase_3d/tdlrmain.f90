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

! ======================================================================
!  This is the main program for caclulating spectra based on the LR-TDDFT .
! ======================================================================

! ======================= history ======================================
!  ver 1.0 :  2010/3/31
!               applicable to the crystal such as Si.  
!  ver 2.0 :  2011/3/31
!               applicable to the isolated molecule such as C6H6.
!
! ======================================================================

program LinearResponseMain

  use m_LinearResponse_Control, only : sw_tddft, sw_LinearResponse, sw_LongWaveLimit,&
       &                               xc_kernel_type, ALDA_R, &
       &                               tddft_eqn_type, DYSON, BS
  use m_Control_Parameters,   only : Num_q_Points 

  implicit none

  logical  :: EigenValues_are_Converged, AllKpoints_are_Calculated
  logical  :: Already_Converged
  logical  :: Ending_Time

  logical  :: all_conv = .true.
  logical  :: Additional_EndFlag = .false.

  integer i

!!  write(*,*) 'First  ', Additional_EndFlag

  call Initialization_set_ekmode_ON  ! set `ekmode' ON in m_Control_Parameters
  call Initialization(1)                ! Initialization of mpi and file-setting
  call InputData_Analysis
  Call InputData_Analysis_LR
  if ( sw_LinearResponse==1 )  Call Gen_q_Points

!--
  call Preparation(0)                ! Basis set, symmetry check etc.
  call Preparation_for_mpi(1)           ! mpi

  call PseudoPotential_Construction
  call Ewald_and_Structure_Factor    ! Calculate Structure Factor

  call Initial_Electronic_Structure()! read Charge Density, (lclchh)

  call Initialization_LR

! -------------------------- kokoga memory kuu --
  call Initialization_Spectrum
  if ( sw_LinearResponse==1 ) then
     Call Alloc_WorkArray_LR1
     Call Alloc_WorkArray_LR2
  endif
  Additional_EndFlag = .false.
!
! ------------------------------------- main -----------
  KPOINTS: do
     call KpointNumber_Setting()
     call Preparation_ek             ! (basnum)
     call Preparation_for_mpi_ek     ! mpi  -> np_g1k, mp_g1k
     call PseudoPotential_ek_Epsilon ! (kbint)                  ! Epsilon
     call Initial_WaveFunctions_ek   ! (rndzaj|rdzaj),(fsrfsi),(lclchh)

     if(.not.Already_Converged()) then
        all_conv = .false.
        SolveWaveFunctions: do
           if(Ending_Time())                 exit KPOINTS
           call IterationNumber_Setting()
           call Renewal_of_WaveFunctions()
           if(EigenValues_are_Converged())   exit SolveWaveFunctions
        enddo SolveWaveFunctions
     else
        exit KPOINTS
     end if
     Call Copy_Wfn_etc_To_WorkArray( 0 ) 

     if ( sw_LongWaveLimit==1 ) then
        Call Calc_RhoTilde_LWLimit
        if ( tddft_eqn_type == DYSON .and. xc_kernel_type == ALDA_R ) then
           call Calc_MatrixG_LWLimit
        endif
     else
        call Calc_Band_Q_Nonzero( Additional_EndFlag )
        if ( Additional_EndFlag ) then
           exit KPOINTS
        endif
        Call Calc_RhoTilde_General
        if ( tddft_eqn_type == DYSON .and. xc_kernel_type == ALDA_R ) then
           call Calc_MatrixG_General
        endif
     endif
! ===============================================
100  continue 
!
     if(AllKpoints_are_Calculated())   then
        all_conv = .true.
        exit KPOINTS
     endif

  End do KPOINTS
!  stop
!! ----
  if ( tddft_eqn_type == BS ) then
     Call Calc_Spectrum_Mol
  else
     Call Calc_Spectrum
  endif
  call Free_WorkArray_LR1
! =============
  call Postprocessing(.false.)
  call WriteDownData_onto_Files_ek()
! ---
  call Free_WorkArray_LR2

999 continue

  if ( sw_LinearResponse==1 ) then
     call Finalization_LR
     call Finalization_Spectrum
     Call Del_q_points
  endif
  call Finalization_of_mpi           ! mpi

contains

  subroutine Calc_Band_Q_Nonzero( EndFlag )
    logical, intent(inout) :: EndFlag 
    integer  :: nq

    Do nq=1, Num_Q_Points
!============================================= k - q ======== 
       call KpointNumber_Setting()
       call Preparation_ek             ! (basnum)
       call Preparation_for_mpi_ek     ! mpi  -> np_g1k, mp_g1k
       call PseudoPotential_ek_Epsilon ! (kbint)                  ! Epsilon
       call Initial_WaveFunctions_ek   ! (rndzaj|rdzaj),(fsrfsi),(lclchh)
       
       if(.not.Already_Converged()) then
          all_conv = .false.
          SolveWaveFunctions2: do
             if(Ending_Time())   then
                EndFlag = .true.
                return
             endif
             call IterationNumber_Setting()
             call Renewal_of_WaveFunctions()
             if (EigenValues_are_Converged()) then
                exit SolveWaveFunctions2
             endif
          enddo SolveWaveFunctions2
       else
          ENDFlag = .true.
          return
       end if
       ! --------
       Call Copy_Wfn_etc_To_WorkArray( 1 )
!       write(*,*) 'nk _ in _ the pro d = ', nk_in_the_process, ENDFLAG
    End do
  end subroutine Calc_Band_Q_Nonzero

end program LinearResponseMain
