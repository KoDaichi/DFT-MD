!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  PROGRAM: TDLRMAIN
!
!  AUTHOR(S): K. Tagami et al   Aug. 1 2011
!
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)
!
!=======================================================================

! ======================================================================
!  This is a collection of subroutines for spectra called in the main program.
! ======================================================================

! ======================= history ======================================
!  ver 1.0 :  2010/3/31
!               applicable to the crystal such as Si.  
!  ver 2.0 :  2011/3/31
!               applicable to the isolated molecule such as C6H6.
!
! ======================================================================

!---------------------------------------------------
!!
!!!           Calc Density etc. 
!!
!---------------------------------------------------
subroutine Calc_RhoTilde_General
  use m_LinearResponse_Density,   only : Add_SoftPart_General, Add_HardPart_General
  use m_PseudoPotential,          only : modnrm
  use m_Const_Parameters,         only : EXECUT

  implicit none

  call Add_SoftPart_General
  if ( modnrm == EXECUT ) call Add_HardPart_General
end subroutine Calc_RhoTilde_General

subroutine Calc_RhoTilde_LWLimit
  use m_LinearResponse_Density,   only : Add_SoftPart_LWLimit, Add_HardPart_LWLimit, &
       &                                 Add_Correction_To_G0, Add_KageShira_To_G0
  use m_PseudoPotential,          only : modnrm
  use m_Const_Parameters,         only : EXECUT

  implicit none

  Call Add_SoftPart_LWLimit
  if ( modnrm == EXECUT )  Call Add_HardPart_LWLimit
  Call Add_Correction_To_G0
  Call Add_KageShira_To_G0
end subroutine Calc_RhoTilde_LWLimit

!------------------------------------------------------
!!
!!!      calc  Matrix G ( DYSON, xc_kernel_type == ALDA_R )
!!
!------------------------------------------------------
subroutine Calc_MatrixG_General
  use m_LinearResponse_ALDA,     only : Add_SoftPart_MatG_General

  implicit none

  Call Add_SoftPart_MatG_General
end subroutine Calc_MatrixG_General

subroutine Calc_MatrixG_LWLimit
  use m_LinearResponse_ALDA,     only : Add_SoftPart_MatG_LWLimit

  implicit none

  call Add_SoftPart_MatG_LWLimit

end subroutine Calc_MatrixG_LWLimit

!---------------------------------------------------
!!
!!!           Open/close file for Spectrum  
!!
!---------------------------------------------------
subroutine Open_File_Spectrum
  use m_Files,              only  : m_Files_open_nf_LR
  use m_Parallelization,    only :  mype
  implicit none

  if ( mype == 0 ) call m_Files_open_nf_LR()
end subroutine Open_File_Spectrum

subroutine Close_File_Spectrum
  use m_Files,              only  : m_Files_close_nf_LR
  use m_Parallelization,    only :  mype
  implicit none

  if ( mype == 0 ) call m_Files_close_nf_LR
end subroutine Close_File_Spectrum

!---------------------------------------------------
!!
!!!         Initialize / Finalize Spectrum  
!!
!---------------------------------------------------
subroutine Initialization_Spectrum
  use m_Const_Parameters,           only : Valence_plus_PC_Charge
  use m_LinearResponse_ALDA,        only : m_LR_alloc_Matrix_G, m_LR_alloc_Matrix_T,&
       &                                   m_LR_alloc_RKernels_ALDA, &
       &                                   m_LR_alloc_igfp, m_LR_set_igfp
  use m_LinearResponse_Control,     only : xc_kernel_type, ALDA_R, ALDA_G, &
       &                                   tddft_eqn_type, DYSON, BS
  use m_LinearResponse_Kernel,      only : m_LR_alloc_Coulomb_Kernels, &
       &                                   m_LR_alloc_XC_Kernels
!  use m_LinearResponse_Spectrum,    only : m_LR_alloc_Chi_Head
  use m_LinearResponse_NonInt,      only : m_LR_alloc_Mat_Chi0
  !!use m_LinearResponse_ALDA,        only : m_LR_alloc_igfp, m_LR_set_igfp
  !!use m_LinearResponse_Control,     only : tddft_eqn_type, DYSON, BS
  implicit none

! -------------------------------- start ----------------
  call m_LR_alloc_Coulomb_Kernels
  if ( tddft_eqn_type == DYSON ) call m_LR_alloc_XC_Kernels
  call m_LR_alloc_Mat_Chi0
  !call m_LR_alloc_Chi_Head
  call alloc_chi_head()

  if ( xc_kernel_type == ALDA_R .or. xc_kernel_type == ALDA_G ) then
     call m_LR_alloc_RKernels_ALDA
     call m_LR_alloc_igfp
     call m_LR_set_igfp
  endif
!
  if ( tddft_eqn_type == DYSON .and. xc_kernel_type == ALDA_R ) then
     Call m_LR_alloc_Matrix_G
     Call m_LR_alloc_Matrix_T
  endif
  call set_Coulomb_and_XC_Kernels
  call Open_File_Spectrum
! ------------------------------- end ----------------
end subroutine Initialization_Spectrum

subroutine Finalization_Spectrum
  use m_LinearResponse_ALDA,        only : m_LR_dealloc_Matrix_G, &
       &                                   m_LR_dealloc_Matrix_T, &
       &                                   m_LR_dealloc_RKernels_ALDA, &
       &                                   m_LR_dealloc_igfp
  use m_LinearResponse_Kernel,      only : m_LR_dealloc_Coulomb_Kernels, &
       &                                   m_LR_dealloc_XC_Kernels
  use m_LinearResponse_NonInt,      only : m_LR_dealloc_Mat_Chi0
!  use m_LinearResponse_Spectrum,    only : &
!       &                                   m_LR_dealloc_Chi_Head, &
!       &                                   m_LR_dealloc_MatL
  use m_LinearResponse_Control,     only : xc_kernel_type, ALDA_R, ALDA_G, &
       &                                   tddft_eqn_type, DYSON, BS
  !!use m_LinearResponse_ALDA,        only : m_LR_dealloc_igfp
  use m_LinearResponse_BS,        only : m_LR_dealloc_MatV, &
       &                                   m_LR_dealloc_MatK_ALDA, &
       &                                   m_LR_dealloc_MatXi, &
       &                                   m_LR_dealloc_MatL0 
  !!use m_LinearResponse_Control,     only : tddft_eqn_type, DYSON, BS

  implicit none

! ---------------------------------- start ---------
  Call m_LR_dealloc_Coulomb_Kernels
  if ( tddft_eqn_type == DYSON ) Call m_LR_dealloc_XC_Kernels
  Call m_LR_dealloc_Mat_Chi0
!  Call m_LR_dealloc_Chi_Head
  call dealloc_chi_head()

  if ( xc_kernel_type == ALDA_R .or. xc_kernel_type == ALDA_G ) then
     call m_LR_dealloc_RKernels_ALDA
     call m_LR_dealloc_igfp
  end if
  if ( tddft_eqn_type == DYSON .and. xc_kernel_type == ALDA_R ) then
     call m_LR_dealloc_Matrix_G
     call m_LR_dealloc_Matrix_T
  endif

  if ( tddft_eqn_type == BS ) then
     if ( xc_kernel_type == ALDA_R ) then
        call m_LR_dealloc_MatK_ALDA
     endif
     call m_LR_dealloc_MatV
     call m_LR_dealloc_MatXi
     call m_LR_dealloc_MatL0
!     call m_LR_dealloc_MatL
     call dealloc_matl()
  endif

  call Close_File_Spectrum
! --------------------------------- end -----------
end subroutine Finalization_Spectrum

!---------------------------------------------------
!!
!!!          Calc Spectrum  in the case of DYSON
!!
!---------------------------------------------------
subroutine Calc_Spectrum
!  use m_LinearResponse_Spectrum,    only :  Write_Spectrum_Header, &
!       &                                    Write_Spectrum_Data, &
!       &                                    Calc_Chi_Head, Mat_Chi0, &
!       &                                    set_nband_LR
  use m_LinearResponse_Density,     only : Calc_FermiLevel_ek
  use m_LinearResponse_Control,     only : nrd_efermi, e, nstep, sw_LongWaveLimit, &
       &                                   xc_kernel_type, ALDA_R
  use m_LinearResponse_NonInt,      only : Calc_Chi0_LWLimit, Calc_Chi0_General, &
       &                                   Check_Total_Charge
  use m_Const_Parameters,           only : Hartree, DP, Valence_plus_PC_Charge
  !!$use m_LinearResponse_Control,     only : xc_kernel_type, ALDA_R
  use m_LinearResponse_ALDA,        only : Calc_MatrixT_LWLimit, Calc_MatrixT_General
!
  !!use m_LinearResponse_Spectrum,  only : set_nband_LR
  use m_Files,                      only : nfout

  implicit none

  integer            :: istep
  real(kind=DP)      :: ene
! ------------------------------------ start -------------
  Call Calc_FermiLevel_ek( nrd_efermi )
! --debug
!  Call Check_Total_Charge
!  stop
! --
!  call Write_Spectrum_Header
!  call set_nband_LR
  call write_header()
  call set_nb()
! ------------------------------------- main ---------------
  Do istep=1, nstep
     ene = e(istep)
     if ( sw_LongWaveLimit==1 ) then
        Call Calc_Chi0_LWLimit( ene )
        if ( xc_kernel_type == ALDA_R ) call Calc_MatrixT_LWLimit( ene )
     else
        Call Calc_Chi0_General( ene )
        if ( xc_kernel_type == ALDA_R ) call Calc_MatrixT_General( ene )
     endif
!     Call Calc_Chi_Head( ene )
!     Call Write_Spectrum_Data( ene )
     Call chi_head( ene )
     call write_data(ene)
     write(nfout,*) '! spectrum at ', istep, '-th ( e=', ene ,' ) is finished.'
  End do
end subroutine Calc_Spectrum

!---------------------------------------------------
!!
!!!         Calc Spectrum in the case of BS && ALDA_R
!!
!---------------------------------------------------
subroutine Calc_Spectrum_Mol
!  use m_LinearResponse_Spectrum,    only :  Write_Spectrum_Header, &
!       &                                    Write_Spectrum_Data, &
!       &                                    Calc_Chi_Head, Mat_Chi0, &
!       &                                    set_nband_LR, m_LR_Alloc_MatL, &
!       &                                    Calc_MatL, Calc_Chi_Head_Mol
  use m_LinearResponse_Density,     only : Calc_FermiLevel_ek
  use m_LinearResponse_Control,     only : nrd_efermi, e, nstep, sw_LongWaveLimit, &
       &                                   xc_kernel_type, ALDA_R
  use m_LinearResponse_NonInt,      only : Calc_Chi0_LWLimit, Calc_Chi0_General, &
       &                                   Check_Total_Charge
  use m_Const_Parameters,           only : Hartree, DP, Valence_plus_PC_Charge
  !!$use m_LinearResponse_Control,     only : xc_kernel_type, ALDA_R
!
!  use m_LinearResponse_BS,  only : m_LR_Alloc_MatK_ALDA, m_LR_Alloc_MatV
!  use m_LinearResponse_BS,  only : m_LR_Alloc_MatXi
!  use m_LinearResponse_BS,  only : m_LR_Alloc_MatL0
!  use m_LinearResponse_BS,  only : Calc_MatV, Calc_MatK_ALDA, Calc_MatXi
!  use m_LinearResponse_BS,  only : Calc_MatL0
  use m_LinearResponse_BS,  only : m_LR_Alloc_MatK_ALDA, m_LR_Alloc_MatV &
   ,m_LR_Alloc_MatXi &
   ,m_LR_Alloc_MatL0 &
   ,Calc_MatV, Calc_MatK_ALDA, Calc_MatXi &
   ,Calc_MatL0 &
   ,Check_Charge_on_mesh, Check_exc_energy &
   ,Print_MatV, Print_MatK


  !!use m_LinearResponse_Spectrum,  only : set_nband_LR, m_LR_Alloc_MatL, &
  !!     &                                 Calc_MatL, Calc_Chi_Head_Mol
  use m_Files,                     only : nfout

  !use m_LinearResponse_BS,  only : Check_Charge_on_mesh, Check_exc_energy
  !use m_LinearResponse_BS,  only : Print_MatV, Print_MatK


! --
  implicit none
  integer :: istep
  Real(kind=DP) :: ene
! --------------------------
  Call Calc_FermiLevel_ek( nrd_efermi )
!  Call Check_Total_Charge
!  stop
! --
!  call Write_Spectrum_Header
!  call set_nband_LR
  call write_header()
  call set_nb()
! --
!  Call Check_Charge_on_mesh
!  Call Check_Exc_energy
!  stop

! ------------
  call m_LR_Alloc_MatV
  call Calc_MatV
  if ( xc_kernel_type == ALDA_R ) then
     call m_LR_Alloc_MatK_ALDA
     call Calc_MatK_ALDA
  endif
! --
!!!!!!!!!!!  Call print_MatV
!!!!!!!!!!!!  if ( xc_kernel_type == ALDA_R ) Call print_MatK
!  stop

  call m_LR_Alloc_MatXi
  call Calc_MatXi
!
  call m_LR_Alloc_MatL0
!  call m_LR_Alloc_MatL
  call alloc_matl()
! ------------------------------------- main ---------------
  Do istep=1, nstep
     ene = e(istep)
!     write(*,*) 'istep = ', istep, ene
     call Calc_MatL0( ene )
!     Call Calc_MatL
     call mat_l()
!     call Calc_Chi_Head_Mol( ene )
!     Call Write_Spectrum_Data( ene )
     call chi_head_mol( ene )
     call write_data(ene)
     write(nfout,*) '! spectrum at ', istep, '-th ( e=', ene ,' ) is finished.'
  End do
!
end subroutine Calc_Spectrum_Mol

! -----------------------------------------------------
!!
!!!               set Coulomb + XC Kernels
!!
!-------------------------------------------------------
subroutine set_Coulomb_and_XC_Kernels
  use m_LinearResponse_Control,   only : xc_kernel_type, RPA, ALDA_G, ALDA_R, &
       &                                LRC, LRC_alpha, &
       &                                sw_Coulomb_Screening
  use m_LinearResponse_Kernel,   only : set_Coulomb_Kernel, &
       &                                set_Screened_Coulomb_Kernel, &
       &                                set_XC_Kernel_RPA, &
       &                                set_XC_Kernel_ALDA_G, &
       &                                set_XC_Kernel_ALDA_R, &
       &                                set_XC_Kernel_LRC
  use m_Const_Parameters,         only : Valence_Plus_PC_Charge
  use m_Control_Parameters,       only : xctype

  implicit none

! -------------------------------- start ----------
!  write(*,*) 'xc_type = ',xctype, xc_kernel_type
  if ( sw_Coulomb_Screening==1 ) then
     Call set_Screened_Coulomb_Kernel
  else
     Call set_Coulomb_Kernel
  endif

  select case (xc_kernel_type)
  case (RPA)
     Call set_XC_Kernel_RPA
  case (LRC)
     Call set_XC_Kernel_LRC( LRC_alpha )
  case (ALDA_G)
     Call set_XC_kernel_ALDA_G( Valence_plus_PC_Charge )
  case (ALDA_R)
     Call set_XC_kernel_ALDA_R( Valence_plus_PC_Charge )
  end select

end subroutine set_Coulomb_and_XC_Kernels

subroutine alloc_chi_head()
  use m_LinearResponse_Spectrum,    only : m_LR_alloc_Chi_Head
  call m_LR_alloc_Chi_Head()
end subroutine alloc_chi_head

subroutine dealloc_chi_head()
  use m_LinearResponse_Spectrum, only : m_LR_dealloc_Chi_Head
  Call m_LR_dealloc_Chi_Head
end subroutine dealloc_chi_head

subroutine dealloc_matl()
  use m_LinearResponse_Spectrum, only : m_LR_dealloc_MatL
  call m_LR_dealloc_MatL
end subroutine dealloc_matl

subroutine write_header()
  use m_LinearResponse_Spectrum,    only :  Write_Spectrum_Header
  implicit none
  call Write_Spectrum_Header()
end subroutine write_header

subroutine chi_head(ene)
  use m_Const_Parameters, only : DP
  use m_LinearResponse_Spectrum, only : Calc_Chi_Head
  implicit none
  real(kind=DP), intent(in) :: ene   
  call Calc_Chi_Head(ene)
end subroutine chi_head

subroutine chi_head_mol(ene)
  use m_Const_Parameters, only : DP
  use m_LinearResponse_Spectrum, only : Calc_Chi_Head_Mol
  implicit none
  real(kind=DP), intent(in) :: ene   
  call Calc_Chi_Head_Mol(ene)
end subroutine chi_head_mol

subroutine write_data(ene)
  use m_LinearResponse_Spectrum,    only :  Write_Spectrum_Data
  use m_Const_Parameters, only : DP
  implicit none
  real(kind=DP), intent(in) :: ene
  call Write_Spectrum_Data(ene)
end subroutine write_data

subroutine set_nb()
  use m_LinearResponse_Spectrum, only : set_nband_LR
  implicit none
  call set_nband_LR
end subroutine set_nb

subroutine alloc_matl()
  use m_LinearResponse_Spectrum, only : m_LR_Alloc_MatL
  implicit none
  call m_LR_Alloc_MatL
end subroutine

subroutine mat_l()
  use m_LinearResponse_Spectrum, only : Calc_MatL
  implicit none
  call Calc_MatL
end subroutine mat_l

