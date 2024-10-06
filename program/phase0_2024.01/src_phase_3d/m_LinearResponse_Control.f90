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
!  This is a module for setting control parameters for the LR-TDDFT .
! ======================================================================

! ======================= history ======================================
!  ver 1.0 :  2010/3/31
!               applicable to the crystal such as Si.  
!  ver 2.0 :  2011/3/31
!               applicable to the isolated molecule such as C6H6.
!
! ======================================================================

module m_LinearResponse_Control

  use m_Electronic_Structure,       only : efermi
  use m_Control_Parameters,         only : ipriinputfile, printable, uvsormode
 
  use m_Const_parameters,           only : DP, SP, BUCS, ON, OFF,CARTS, & 
       &              PAI4, PAI2, PAI, FMAXVALLEN, LOWER, BOHR, &
       &              MESH

  use m_Files,                      only : nfout, m_Files_reopen_nfinp
  use m_Kpoints,                    only : way_ksample

  use m_Control_Parameters,          only : sw_LinearResponse

  Implicit None
  include 'mpif.h'
  
! -----------  General ------------------------
  integer                ::         iprilrspectr = 1         !     OFF=0, ON=1
  integer, parameter     ::         OPTICS=0,  PACS=1, EELS=2, IXSS=3
  integer                ::         spectrum_type = OPTICS
! ------------ momentum transfer -------
  real(kind=DP)          ::         vec_q(3), norm_q
  real(kind=DP)          ::         min_norm_q = 1.0D-3
! -----------  TDDFT -------------------------
  integer                ::         sw_tddft  = OFF
  integer                ::         sw_LongWaveLimit  = ON
! ------------ Solver --------
  integer, parameter     ::         DYSON = 0,  BS = 1     
  integer                ::         tddft_eqn_type = DYSON        ! default
! --
  integer                ::         sw_NODA = OFF       ! approximatin in BS eq.
! ------------ Coulomb Kernel ---------
  integer                ::         sw_NLF    = OFF
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< EXPERIMENTAL <<<<<<<<<<<<<<<<<
  integer                ::         sw_Coulomb_screening = OFF      ! default
  integer                ::         sw_Fxc_enhanced  = OFF      ! default
  Real(kind=DP)          ::         rcut_yukawa      = 1.0 / Bohr
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ----------- XC Kernel ---------
  integer, parameter     ::         RPA = 0, ALDA_G = 1,  LRC = 2, ALDA_R = 11
  integer                ::         xc_kernel_type = RPA     ! default
  real(kind=DP)          ::         LRC_alpha = 1.0
  real(kind=DP)          ::         LRC_beta  = 1.0
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< EXPERIMENTAL <<<<<<<<<<<<<<<<<
  Real(kind=DP)          ::         fxc_enhanced_facotr = 1.0
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------ Expansion ----------
  integer                ::         nmax_G_LR = 100        ! default
! ------------ BZ ---------------------
  integer, parameter     ::         LORENTZIAN_B=1,  GAUSSIAN_B=2, L_TETRAHEDRON=3
  integer                ::         way_BZintegral = LORENTZIAN_B
!
  real(kind=DP)          ::         eta = 1.0E-4    ! in hartree unit
!
  integer                ::         nistep, nstep
  real(kind=DP)          ::         width, deltaq, tetra_eps
  real(kind=DP)          ::         scissor = 0.0D0
! ------------- Energy range -------------
  real(kind=DP)          ::         e_low, e_high, e_step
! ------------- Fermi ----------------
  integer                ::         nrd_efermi = 0
  real(kind=DP)          ::         efermi_backup

! ----------------------- Arrays to be allocated -----------
  Real(kind=DP),    allocatable     :: vqxyz( :,:,: )
  real(kind=DP),    allocatable     :: e(:)
!
  integer                ::         ipriqpt = 0
! ---------------------------------------------------------
contains

  subroutine m_CtrlP_rd_LinearResponse( nfout )
    integer, intent(in) :: nfout
! ------------------------------ General ----------------------
    character(len("spectrum")), parameter  ::  tag_spectrum       = "spectrum"
    character(len("LongWaveLimit")),  parameter  ::  tag_longwave  = "LongWaveLimit"
! --
    character(len("type")),     parameter  ::  tag_spectrum_type  = "type"
    character(len("optics")),   parameter  ::  tag_optics         = "optics" 
    character(len("PACS")),     parameter  ::  tag_PACS           = "PACS" 
    character(len("EELS")),     parameter  ::  tag_EELS           = "EELS" 
    character(len("IXSS")),     parameter  ::  tag_IXSS           = "IXSS" 
! ------------------------------- Momentum Transfer ----------
    character(len("momentum_transfer")), parameter  ::  tag_momentum = "momentum_transfer"
    character(len("deltaq")),   parameter  ::  tag_deltaq   = "deltaq"
    character(len("nx")),       parameter  ::  tag_nx       = "nx"
    character(len("ny")),       parameter  ::  tag_ny       = "ny"
    character(len("nz")),       parameter  ::  tag_nz       = "nz"
! ---------------------------------- TDDFT --------------------
    character(len("tddft")),    parameter  ::  tag_tddft    = "tddft"
    character(len("sw_tddft")), parameter  ::  tag_sw_tddft    = "sw_tddft"
! ----------------------------------- Solver ----------------
    character(len("solver")),   parameter  ::  tag_solver  = "solver"
    character(len("equation")),   parameter  ::  tag_equation  = "equation"
    character(len("DYSON")),    parameter  ::  tag_DYSON   = "DYSON"
    character(len("BS")),       parameter  ::  tag_BS      = "BS"
    character(len("sw_NODA")),  parameter  ::  tag_sw_NODA   = "sw_NODA"
! ----------------------------------- Coulomb Kernel -----------
    character(len("Coulomb_Kernel")),parameter  ::  tag_CoulombKernel = "Coulomb_Kernel"
    character(len("sw_NLF")),   parameter  ::  tag_sw_NLF   = "sw_NLF"
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< experimenal <<<<<<<<<
    character(len("sw_Coulomb_screening")), parameter :: tag_sw_Coulomb_screening &
         &                                   = "sw_Coulomb_screening"
    character(len("screening_length")),  parameter  ::  tag_screening_length &
         &                                   = "screening_length"
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ---------------------------------- XC Kernel --------------------
    character(len("XC_Kernel")),parameter  ::  tag_XCKernel = "XC_Kernel"
    character(len("kernel_type")), parameter  ::  tag_kernel_type = "kernel_type"
    character(len("RPA")),      parameter  ::  tag_RPA      = "RPA"
    character(len("ALDA-G")),   parameter  ::  tag_ALDA_G   = "ALDA-G"
    character(len("ALDA-R")),   parameter  ::  tag_ALDA_R   = "ALDA-R"
    character(len("LRC")),      parameter  ::  tag_LRC      = "LRC"
    character(len("LRC_alpha")),parameter  ::  tag_LRC_alpha = "LRC_alpha"
    character(len("LRC_beta")), parameter  ::  tag_LRC_beta  = "LRC_beta"
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< experimenal <<<<<<<<<
    character(len("sw_Fxc_enhanced")), parameter :: tag_sw_Fxc_enhanced &
         &                                   = "sw_Fxc_enhanced"
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! -------------------------------- Expansion --------------
    character(len("Expansion")),parameter  ::  tag_expansion = "expansion"
    character(len("NumGVec")),  parameter  ::  tag_num_GVec  = "NumGvec"
! ----------------------------------- energy -------------------
    character(len("energy")),   parameter  ::  tag_energy   = "energy"
    character(len("low")),      parameter  ::  tag_low      = "low"
    character(len("high")),     parameter  ::  tag_high     = "high"
    character(len("step")),     parameter  ::  tag_step     = "step"
! --------------------------------- BZ integration ------------
    character(len("BZ_integration")), parameter :: tag_BZ_integration = "BZ_integration"
    character(len("method")), parameter         :: tag_method = "method"
!
    character(len("lorentzian")), parameter     :: tag_lorentzian = "lorentzian"
    character(len("l")), parameter              :: tag_l = "l"
!
    character(len("gaussian")), parameter       :: tag_gaussian = "gaussian"
    character(len("g")), parameter              :: tag_g = "g"
    character(len("width")), parameter          :: tag_width = "width"
! ----------------------------- --- not supported ----------------
    character(len("tetrahedron")), parameter    :: tag_tetrahedron = "tetrahedron"
    character(len("t"))                         :: tag_t = "t"
    character(len("nistep")), parameter         :: tag_nistep = "nistep"
    character(len("spin")),  parameter          :: tag_spin  = "spin"
    character(len("both")),  parameter          :: tag_both  = "both"
    character(len("tetra_eps")),parameter       :: tag_tetra_eps = "tetra_eps"
! ------------------------------ band gap correction ----------
    character(len("band_gap_correction")), parameter   ::  tag_band_gap_correction = "band_gap_correction"
    character(len("scissor_operator")), parameter      ::  tag_scissor  = "scissor_operator"
! ----------------------------------- Fermi --------------------
    character(len("on")), parameter              ::  tag_on  = "on"
    character(len("off")), parameter             ::  tag_off = "off"
    character(len("fermi_energy")), parameter    ::  tag_fermi_energy = "fermi_energy"
    character(len("read_efermi")), parameter     ::  tag_read_efermi = "read_efermi"
    character(len("efermi")),parameter           ::  tag_efermi = "efermi"
! --

    character( len=FMAXVALLEN ) :: rstr
    integer :: iret, f_selectBlock, f_selectParentBlock, f_selectTop
    integer :: f_getStringValue, f_getRealValue, f_getIntValue
    real(kind=DP) :: dret
    
! ---------------------------------------------------------------------
! ----------------------------------------- start ---------------------
! ---------------------------------------------------------------------
    iret = f_selectTop()
    if ( ipriinputfile >= 2 .and. printable ) then
       write(nfout,'(" !*  tag_spectrum")')
    endif
    if ( f_selectBlock( tag_spectrum ) == 0 ) then
       sw_LinearResponse = ON
       call Set_Spectrum_Type
       call Set_switch_LongWaveLimit
! ----- momentum transfer
       if ( f_selectBlock( tag_momentum ) == 0 ) then
          call Set_Momentum_Transfer(1)
          iret = f_selectParentBlock()
       else
          call Set_Momentum_Transfer(0)
       endif
! ---- tddft --
       if ( f_selectBlock( tag_tddft ) == 0 ) then
          call Set_switch_TDDFT
          if ( f_selectBlock( tag_solver ) == 0 ) then
             call set_Equation_type
             call Set_switch_NODA
             iret = f_selectParentBlock()
          endif
! ---------------------------- Coulomb Kernel --------
          if ( f_selectBlock( tag_CoulombKernel ) == 0 ) then
             call Set_switch_NLF
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< experimental <<<<<<<<<<<<<<<<<
             call set_switch_Coulomb_screening
             call Set_screening_length(1)
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
             iret = f_selectParentBlock()
          endif
! ---------------------------- XC Kernel --------------
          if ( f_selectBlock( tag_XCKernel ) == 0 ) then
             call Set_XCKernel_Type
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< experimental <<<<<<<<<<<<<<<<<
             call set_switch_Fxc_enhanced
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
             iret = f_selectParentBlock()
          endif
          call Limitation_Equation_Type
! ------------------------------ G-expansion ----------
          if ( f_selectBlock( tag_expansion ) == 0 ) then
             call Set_Num_GVectors(1)
             iret = f_selectParentBlock()
          else
             call Set_Num_GVectors(0)
          endif
          iret = f_selectParentBlock()
       endif
! ------------- Fermi -----------
       if ( f_selectBlock( tag_fermi_energy ) == 0 ) then
          if( f_getStringValue( tag_read_efermi, rstr, LOWER) == 0) then
             call set_efermi(rstr)
          endif
          iret = f_selectParentBlock()
       else
          call Set_FermiLevel_Later
       end if
! --- Energy --
       if ( f_selectBlock( tag_energy ) == 0 ) then
          call set_energy_range
          iret = f_selectParentBlock()
       end if
! --- BZ Integral
       if (f_selectBlock(tag_BZ_integration) == 0 ) then
          call Set_BZ_Integration_Type
          iret = f_selectParentBlock()
       end if
! --- read band gap correction option
       if ( f_selectBlock(tag_band_gap_correction) == 0 ) then
          call set_band_gap_correction
          iret = f_selectParentBlock()
       end if
       write(nfout,*) '***********************************************'
! --- 
    else
       stop ' tag_spectrum is not given in the input file <<CtrlP_rd_spectrum>> '
    end if
    
  contains

!----------------------------------------------------------------
!!
!!!               Spectrum type
!!
!----------------------------------------------------------------
    subroutine Set_Spectrum_Type
      logical :: tf

      write(nfout,*) '*************************** Spectrum ***********'

      if ( f_getStringValue( tag_spectrum_type, rstr, LOWER ) == 0 ) then
         call strncmp0( tag_optics, trim(rstr), tf )
         if(tf) then
            spectrum_type = OPTICS;   goto 1001
         end if
         call strncmp0( tag_PACS, trim(rstr), tf )
         if(tf) then
            spectrum_type = PACS;      goto 1001
         end if
         call strncmp0( tag_EELS, trim(rstr), tf )
         if(tf) then
            spectrum_type = EELS;      goto 1001
         end if
         call strncmp0( tag_IXSS, trim(rstr), tf )
         if(tf) then
            spectrum_type = IXSS;      goto 1001
         end if
         if ( printable ) call Comment_Spectrum_Type( 1 )

1001     continue
         if ( printable ) call Comment_Spectrum_Type( 0 )
      else
         call Comment_Spectrum_Type( -1 )
      endif
    end subroutine Set_Spectrum_Type

    subroutine Comment_Spectrum_Type( icode )
      integer, intent(in) :: icode

      select case( icode )
      case(0)
         write(nfout,'("!* Spectrum type = ",i3)') spectrum_type
         write(nfout,*) "!* ------ 0:Optics; 1:PACS, 2:EELS, 3:IXSS -- "
      case (-1)
         write(nfout,'(1x,"!* tag_spectrum_type is not found")')
         write(nfout,'(1x,"!* spectrum_type is set to OPTICS")')
      case (1)
         write(nfout,'(1x,"!* incorrect words for spectrum_type")')
         write(nfout,'(1x,"!* spectrum type is set to OPTICS")')
      end select

    end subroutine Comment_Spectrum_Type
    
!----------------------------------------------------------------
!!
!!!               LongWaveLimit
!!
!----------------------------------------------------------------
    subroutine Set_switch_LongWaveLimit
      logical :: tf

      if ( f_getStringValue( tag_longwave, rstr, LOWER ) == 0 ) then
         call strncmp0( tag_off, trim(rstr), tf )
         if (tf) then
            sw_LongWaveLimit = OFF; goto 1002
         endif
         call strncmp0( tag_on, trim(rstr), tf )
         if (tf) then
            sw_LongWaveLimit = ON; goto 1002
         endif
         if ( printable ) Call Comment_sw_LongWaveLimit(1)

1002     continue
         if ( printable ) Call Comment_sw_LongWaveLimit(0)
      else
         if ( printable ) Call Comment_sw_LongWaveLimit(-1)
      endif
    end subroutine Set_switch_LongWaveLimit

    subroutine Comment_sw_LongWaveLimit( icode )
      integer, intent(in) :: icode

      select case( icode )
      case(0)
         if ( sw_LongWaveLimit == ON ) then
            write(nfout,'(1x,"!*         LongWaveLimit set ON")')
         else
            write(nfout,'(1x,"!*         LongWaveLimit is set OFF")')
         endif
      case(-1)
         write(nfout,'(1x,"!* tag_longwave is not found")')
         write(nfout,'(1x,"!* LongWaveLimit is set ON as default")')
      case(1)
         write(nfout,'(1x,"!* incorrect words for LongWaveLimit")')
         write(nfout,'(1x,"!* LongWaveLimit is set ON as default")')
      end select

    end subroutine Comment_sw_LongWaveLimit

!----------------------------------------------------------------
!!
!!!               Momentum transfer
!!
!----------------------------------------------------------------
    subroutine Set_Momentum_Transfer( itype )
      integer :: itype
      real(kind=DP) c1
      integer :: f_getRealValue

      select case ( itype )
      case (0)
         vec_q = 0.0d0; vec_q(3) = min_norm_q;
         vec_q = vec_q *Bohr
         call Comment_MomentumTransfer( 0 )
      case (1)
         if ( f_getRealValue( tag_deltaq, dret," " ) == 0 ) norm_q = dret
         if ( f_getRealValue( tag_nx, dret," " ) == 0 ) vec_q(1) = dret
         if ( f_getRealValue( tag_ny, dret," " ) == 0 ) vec_q(2) = dret
         if ( f_getRealValue( tag_nz, dret," " ) == 0 ) vec_q(3) = dret
         
         if ( vec_q(1)**2 +vec_q(2)**2 +vec_q(3)**2 == 0.0d0 ) vec_q = 0.0d0

!      if ( qnorm < qnorm_min .or. sw_LongWaveLimit==ON ) qnorm = qnorm_min
         if ( sw_LongWaveLimit==ON ) norm_q = min_norm_q
! --
         call Comment_MomentumTransfer( 1 )
         norm_q = norm_q *Bohr
         c1 = sqrt( vec_q(1)**2 + vec_q(2)**2 + vec_q(3)**2 )
         vec_q = norm_q * vec_q / c1
      end select

    end subroutine Set_Momentum_Transfer

    subroutine Comment_MomentumTransfer( icode )
      integer, intent(in) :: icode

      select case( icode )
      case(0)
         write(nfout,'(1x,"!* momentum transfer is set to a default value  ",F12.5," Angstrom^(-1)")') min_norm_q
         write(nfout,'(1x,"!*         Direction = ",3F12.5)') 0, 0, 1.0
      case (1)
         write(nfout,'(1x,"!* momentum transfer = ",F12.5," Angstrom^(-1)")') norm_q
         write(nfout,'(1x,"!*         Direction = ",3F12.5)') vec_q(1), vec_q(2), vec_q(3)
      end select
      
    end subroutine Comment_MomentumTransfer
    
!----------------------------------------------------------------
!!
!!!               TDDFT
!!
!----------------------------------------------------------------
    subroutine set_switch_TDDFT
      logical :: tf

      if ( f_getStringValue( tag_sw_tddft, rstr, LOWER ) == 0 ) then
         call strncmp0( tag_off, trim(rstr), tf )
         if (tf) then
            sw_tddft = OFF; goto 1003
         endif
         call strncmp0( tag_on, trim(rstr), tf )
         if (tf) then
            sw_tddft = ON; goto 1003
         endif
         if ( printable ) Call Comment_sw_tddft(1)
1003     continue
         if ( printable ) Call Comment_sw_tddft(0)
      else
         if ( printable ) Call Comment_sw_tddft(-1)
      endif
    end subroutine set_switch_TDDFT

    subroutine Comment_sw_tddft( icode )
      integer, intent(in) :: icode

      select case(icode)
      case(0)
         if ( sw_tddft == ON ) then
            write(nfout,*) "!* sw_tddft is set to ON"
         else
            write(nfout,*) "!* sw_tddft is set to OFF"
         endif
      case(-1)
         write(nfout,'(1x,"!* tag_sw_tddft is not found")')
         write(nfout,'(1x,"!* sw_tddft is set OFF as default")')
      case(1)
         write(nfout,'(1x,"!* incorrect words for sw_tddft")')
         write(nfout,'(1x,"!* sw_tddft is set OFF as default")')
      end select
    end subroutine Comment_sw_tddft

!----------------------------------------------------------------
!!
!!!               Solver       ( Dyson or BS )
!!
!----------------------------------------------------------------
    subroutine Set_Equation_Type
      logical :: tf

      if ( f_getStringValue( tag_equation, rstr, LOWER ) == 0 ) then
         call strncmp0( tag_DYSON, trim(rstr), tf )
         if(tf) then
            tddft_eqn_type  = DYSON;   goto 1006
         end if
         call strncmp0( tag_BS, trim(rstr), tf )
         if(tf) then
            tddft_eqn_type = BS;   goto 1006
         end if

         call Comment_Equation_Type(1)
!
1006     continue
         call Comment_Equation_Type(0)
      else
         call Comment_Equation_Type(-1)
      endif
    end subroutine Set_Equation_Type

    subroutine Comment_Equation_Type( icode )
      integer, intent(in) :: icode

      select case(icode)
      case(0)
         write(nfout,'("!* TDDFT Solver = ",i3)') tddft_eqn_type
         write(nfout,'(1x,"!*         DYSON_like : 0,  BS_like : 1  ")')
      case(-1)
         write(nfout,'(1x,"!* tag_equation is not found")')
         write(nfout,'(1x,"!* tddft_eqn_type is set to DYSON_like")')
      case (1)
         write(nfout,'(1x,"!* incorrect words for tag_solver")')
         write(nfout,'(1x,"!* tddft_eqn_type is set to DYSON_like")')
      end select
    end subroutine Comment_Equation_Type

    subroutine Limitation_Equation_Type
      if ( tddft_eqn_type == BS ) then
         if ( xc_kernel_type /= ALDA_R .and. xc_kernel_type /= RPA ) then
            write(nfout,'(1x,"!** BS_Like equation is supported with xc_kernel = ALDA_R/RPA")')
            stop
         endif
      endif
    end subroutine Limitation_Equation_Type

    subroutine Set_switch_NODA           ! Neglect Off-Diagonal Approximation
      logical :: tf

      if ( tddft_eqn_type /= BS ) return

      if ( f_getStringValue( tag_sw_NODA, rstr, LOWER ) == 0 ) then
         call strncmp0( tag_off, trim(rstr), tf )
         if (tf) then
            sw_NODA = OFF; goto 1011
         endif
         call strncmp0( tag_on, trim(rstr), tf )
         if (tf) then
            sw_NODA = ON; goto 1011
         endif
         if ( printable ) Call Comment_sw_NODA(1)
1011     continue
         if ( printable ) Call Comment_sw_NODA(0)
      else
         if ( printable ) Call Comment_sw_NODA(-1)
      endif
    end subroutine Set_switch_NODA

    subroutine Comment_sw_NODA( icode )
      integer, intent(in) :: icode

      select case(icode)
      case(0)
         if ( sw_NODA == ON ) then
            write(nfout,*) "!* sw_NODA is set to ON"
         else
            write(nfout,*) "!* sw_NODA is set to OFF"
         endif
      case(-1)
         write(nfout,'(1x,"!* tag_sw_NODA is not found")')
         write(nfout,'(1x,"!* sw_NODA is set OFF as default")')
      case(1)
         write(nfout,'(1x,"!* incorrect words for sw_NODA")')
         write(nfout,'(1x,"!* sw_NODA is set OFF as default")')
      end select
    end subroutine Comment_sw_NODA

!----------------------------------------------------------------
!!
!!!              Coulomb Kernel
!!
!----------------------------------------------------------------
    subroutine Set_switch_NLF
      logical :: tf

      if ( f_getStringValue( tag_sw_NLF, rstr, LOWER ) == 0 ) then
         call strncmp0( tag_off, trim(rstr), tf )
         if (tf) then
            sw_NLF = OFF; goto 1004
         endif
         call strncmp0( tag_on, trim(rstr), tf )
         if (tf) then
            sw_NLF = ON; goto 1004
         endif
         if ( printable ) Call Comment_sw_NLF(1)
1004     continue
         if ( printable ) Call Comment_sw_NLF(0)
      else
         if ( printable ) Call Comment_sw_NLF(-1)
      endif
    end subroutine Set_switch_NLF

    subroutine Comment_sw_NLF( icode )
      integer, intent(in) :: icode

      select case(icode)
      case(0)
         if ( sw_NLF == ON ) then
            write(nfout,*) "!* sw_NLF is set to ON"
         else
            write(nfout,*) "!* sw_NLF is set to OFF"
         endif
      case(-1)
         write(nfout,'(1x,"!* tag_sw_NLF is not found")')
         write(nfout,'(1x,"!* sw_NLF is set OFF as default")')
      case(1)
         write(nfout,'(1x,"!* incorrect words for sw_NLF")')
         write(nfout,'(1x,"!* sw_NLF is set OFF as default")')
      end select
    end subroutine Comment_sw_NLF

! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< EXPERIMENTAL <<<<<<<<<<<<<<<<<<<<<<<
    subroutine Set_switch_Coulomb_Screening
      logical :: tf

      if ( f_getStringValue( tag_sw_Coulomb_screening, rstr, LOWER ) == 0 ) then
         call strncmp0( tag_off, trim(rstr), tf )
         if (tf) then
            sw_Coulomb_screening = OFF; goto 1012
         endif
         call strncmp0( tag_on, trim(rstr), tf )
         if (tf) then
            sw_Coulomb_screening = ON; goto 1012
         endif
         if ( printable ) Call Comment_sw_Coulomb_Screening(1)
1012     continue
         if ( printable ) Call Comment_sw_Coulomb_Screening(0)
      else
         if ( printable ) Call Comment_sw_Coulomb_Screening(-1)
      endif
    end subroutine Set_switch_Coulomb_Screening

    subroutine Comment_sw_Coulomb_Screening( icode )
      integer, intent(in) :: icode
      
      select case(icode)
      case(0)
         if ( sw_Coulomb_Screening == ON ) then
            write(nfout,*) "!* sw_Coulomb_Screening is set to ON"
         else
            write(nfout,*) "!* sw_Coulomb_Screening is set to OFF"
         endif
      case(-1)
         write(nfout,'(1x,"!* tag_sw_Coulomb_Screening is not found")')
         write(nfout,'(1x,"!* sw_Coulomb_Screening is set OFF as default")')
      case(1)
         write(nfout,'(1x,"!* incorrect words for sw_Coulomb_Screeing")')
         write(nfout,'(1x,"!* sw_Coulomb_Screening is set OFF as default")')
      end select
    end subroutine Comment_sw_Coulomb_Screening

    subroutine Set_screening_length( icode )
      integer, intent(in) :: icode
      integer :: f_getRealValue

      if ( sw_Coulomb_Screening == OFF ) return

      select case(icode)
      case(0)
         call Comment_screening_length(-2)
      case(1)
         if ( f_getRealValue( tag_screening_length, dret, "" ) == 0 ) then
            rcut_yukawa = dret
            call Comment_screening_length(0)
         else
            call Comment_screening_length(-1)
         endif
      end select

    end subroutine Set_screening_length

    subroutine Comment_screening_length( icode )
      integer, intent(in) :: icode
      
      select case(icode)
      case(0)
         write(nfout,*) "!* The screening lengh is set to ", rcut_yukawa, " bohr"
      case(-1)
         write(nfout,*) "!* tag_screening_length is not found"
         write(nfout,*) "!* The screening length is set to a default value ", rcut_yukawa
      case(-2)
         write(nfout,*) "!* tag_screening_length is not found"
         write(nfout,*) "!* The screening length is set to a default value ", rcut_yukawa
      end select
    end subroutine Comment_screening_length

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!----------------------------------------------------------------
!!
!!!              Exchange correlatio Kernel
!!
!----------------------------------------------------------------
    subroutine Set_XCKernel_Type
      logical :: tf

      if ( f_getStringValue( tag_kernel_type, rstr, LOWER ) == 0 ) then
         call strncmp0( tag_RPA, trim(rstr), tf )
         if(tf) then
            xc_kernel_type = RPA;   goto 1005
         end if
         call strncmp0( tag_ALDA_G, trim(rstr), tf )
         if(tf) then
            xc_kernel_type = ALDA_G;      goto 1005
         end if
         call strncmp0( tag_ALDA_R, trim(rstr), tf )
         if(tf) then
            xc_kernel_type = ALDA_R;      goto 1005
         end if
         call strncmp0( tag_LRC, trim(rstr), tf )
         if(tf) then
            xc_kernel_type = LRC;      goto 1005
         end if

         call Comment_XCKernel_Type(1)
!
1005     continue
         call Comment_XCKernel_Type(0)
         if ( xc_kernel_type == LRC ) call Set_XCKernel_Param_LRC
      else
         call Comment_XCKernel_Type(-1)
      endif
    end subroutine Set_XCKernel_Type

    subroutine Comment_XCKernel_Type( icode )
      integer, intent(in) :: icode

      select case(icode)
      case(0)
         write(nfout,'(" !* XC Kernel type = ",i3)') xc_kernel_type
         write(nfout,'(1x,"!*         RPA : 0,  ALDA-G : 1, LRC : 2  ")')
         write(nfout,'(1x,"!*                   ALAD-R :11  ")')
      case(-1)
         write(nfout,'(1x,"!* tag_xc_kernel is not found")')
         write(nfout,'(1x,"!* xc_kernel_type is set to RPA")')
      case (1)
         write(nfout,'(1x,"!* incorrect words for xc_kernel")')
         write(nfout,'(1x,"!* xc_kernel_type is set to RPA")')
      end select
    end subroutine Comment_XCKernel_Type

    subroutine Set_XCKernel_Param_LRC
      integer :: f_getRealValue

      if ( f_getRealValue( tag_LRC_alpha, dret,"" ) == 0 ) then
         LRC_alpha = dret
         call Comment_XCKernel_LRC(0)
      else
         call Comment_XCKernel_LRC(20)
      endif
    end subroutine Set_XCKernel_Param_LRC

    subroutine Comment_XCKernel_LRC( icode )
      integer, intent(in) :: icode
      
      select case(icode)
      case(0)
         write(nfout,*) "!* LRC_alpha is set to", LRC_alpha         
      case(20)
         write(nfout,*) "!* tag_LRC_alpha is not found"
         write(nfout,*) "!* The parameter LRC_alpha is set to a dafault value", LRC_alpha
      end select
    end subroutine Comment_XCKernel_LRC

! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< EXPERIMENTAL <<<<<<<<<<<<<<<<<<<<<<
    subroutine Set_switch_Fxc_enhanced
      logical :: tf

      if ( f_getStringValue( tag_sw_Fxc_enhanced, rstr, LOWER ) == 0 ) then
         call strncmp0( tag_off, trim(rstr), tf )
         if (tf) then
            sw_Fxc_enhanced = OFF; goto 1015
         endif
         call strncmp0( tag_on, trim(rstr), tf )
         if (tf) then
            sw_Fxc_enhanced = ON; goto 1015
         endif
         if ( printable ) Call Comment_sw_Fxc_enhanced(1)
1015     continue
         if ( printable ) Call Comment_sw_Fxc_enhanced(0)
      else
         if ( printable ) Call Comment_sw_Fxc_enhanced(-1)
      endif
    end subroutine Set_switch_Fxc_enhanced

    subroutine Comment_sw_Fxc_enhanced( icode )
      integer, intent(in) :: icode
      
      select case(icode)
      case(0)
         if ( sw_Fxc_enhanced == ON ) then
            write(nfout,*) "!* sw_Fxc_enhanced is set to ON"
         else
            write(nfout,*) "!* sw_Fxc_enhanced is set to OFF"
         endif
      case(-1)
         write(nfout,'(1x,"!* tag_sw_Fxc_enhanced is not found")')
         write(nfout,'(1x,"!* sw_Fxc_enhanced set OFF as default")')
      case(1)
         write(nfout,'(1x,"!* incorrect words for sw_Fxc_enhanced")')
         write(nfout,'(1x,"!* sw_Fxc_enhanced set OFF as default")')
      end select
    end subroutine Comment_sw_Fxc_enhanced
    
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!----------------------------------------------------------------
!!
!!!              G-expansion
!!
!----------------------------------------------------------------
    subroutine Set_Num_GVectors( icode )
      integer, intent(in) :: icode
      integer :: f_getIntValue

      select case(icode)
      case(0)
         call Comment_GVecTruncated(-2)
      case(1)
         if ( f_getIntValue( tag_Num_GVec, iret ) == 0 ) then
            nmax_G_LR = iret
            call Comment_GVecTruncated(0)
         else
            call Comment_GVecTruncated(-1)
         endif
      end select

    end subroutine Set_Num_GVectors

    subroutine Comment_GVecTruncated( icode )
      integer, intent(in) :: icode
      
      select case(icode)
      case(0)
         write(nfout,*) "!* The number of G-Vector is limited to ", nmax_G_LR
         write(nfout,*) "!* This value may be changed later, depending on kgp"
      case(-1)
         write(nfout,*) "!* tag_Num_GVec is not found"
         write(nfout,*) "!* The number of G-Vector is set to a default value ", nmax_G_LR
      case(-2)
         write(nfout,*) "!* tag_trunctated is not found"
         write(nfout,*) "!* The number of G-Vector is set to a default value ", nmax_G_LR
      end select
    end subroutine Comment_GVecTruncated

!----------------------------------------------------------------
!!
!!!               Fermi ????????????
!!
!----------------------------------------------------------------
    subroutine Set_FermiLevel_Later
      write(nfout,*) "!* Fermi energy is set later"
    end subroutine Set_FermiLevel_Later

! ------------------------------------------------------------------
!!
!!! ---------------  Fermi, Range, BZ, Scissor  ----------------
!!
! ------------------------------------------------------------------

    subroutine set_efermi(rstr)
      character(len=FMAXVALLEN),intent(in) :: rstr

      integer :: f_getRealValue
      logical :: tf

      call strncmp0(tag_off, trim(rstr), tf)
      if(tf) then
         nrd_efermi=0;  goto 1001
      end if
      call strncmp0(tag_on, trim(rstr), tf)
      if(tf) then
         nrd_efermi=1
         if( f_getRealValue( tag_efermi, dret,'hartree') == 0) efermi=dret
         goto 1001
      end if
1001  if ( printable ) then
         write(nfout,'(1x,"!* nrd_efermi = ",i3)') nrd_efermi
         write(nfout,'(1x,"!*     Fermi level is set to ",f15.8)') efermi
      endif
      
    end subroutine set_efermi
      
    subroutine set_energy_range
      integer :: f_getRealValue
      if( f_getRealValue( tag_low, dret, 'hartree') == 0)  e_low=dret
      if( f_getRealValue( tag_high, dret,'hartree') == 0) e_high=dret
      if( f_getRealValue( tag_step, dret,'hartree') == 0) e_step=dret
      !       write(*,*) 'e_low = ', e_low, e_high, e_step
    end subroutine set_energy_range

    subroutine Set_BZ_Integration_Type
      integer :: f_getRealValue, f_getIntValue

      if ( f_getStringValue( tag_method, rstr, LOWER) == 0 ) call set_integration_method(rstr)
      if ( way_BZintegral == LORENTZIAN_B ) then
         if( f_getRealValue( tag_width, dret, 'hartree' ) == 0) then
            if ( dret > 0.0 ) then 
               eta = dret
            endif
         endif
      end if
      if ( way_BZintegral == GAUSSIAN_B ) then
         if( f_getRealValue( tag_width, dret, 'hartree' ) == 0) width = dret
      end if
      if ( way_BZintegral==L_TETRAHEDRON ) then
         if( f_getIntValue( tag_nistep, iret) == 0) nistep=iret
         if( f_getRealValue( tag_tetra_eps, dret, 'hartree' ) == 0) tetra_eps=dret
      end if
    end subroutine Set_BZ_Integration_Type

    subroutine set_integration_method(rstr)
      character(len=FMAXVALLEN),intent(in) :: rstr
      logical :: tf
      way_BZintegral = LORENTZIAN_B
      
      call strncmp0(tag_lorentzian, trim(rstr), tf)
      if(tf) then
         way_BZintegral = LORENTZIAN_B
         goto 1001
      end if
      call strncmp0(tag_l, trim(rstr), tf)
      if(tf) then
         way_BZintegral = L_TETRAHEDRON
         goto 1001
      end if
!
      call strncmp0(tag_gaussian,trim(rstr),tf)
      if(tf) then
         way_BZintegral = GAUSSIAN_B
         goto 1001
      end if
      call strncmp0(tag_g,trim(rstr),tf)
      if(tf) then
         way_BZintegral = GAUSSIAN_B
         goto 1001
      end if

      call strncmp0(tag_tetrahedron, trim(rstr), tf)
      if(tf) then
         way_BZintegral = L_TETRAHEDRON
         goto 1001
      end if
      call strncmp0(tag_t, trim(rstr), tf)
      if(tf) then
         way_BZintegral = L_TETRAHEDRON
         goto 1001
      end if
      
      stop ' ! tag for BZ_integration method is invalid <<m_CtrlP_rd_LinearResponse>>'
1001  continue
      if ( ipriinputfile >= 2 .and. printable) then
         write(nfout,'(" !* Brillouin zone integration  method = ",a10)') trim(rstr)
      endif
      
    end subroutine set_integration_method

    subroutine set_band_gap_correction
      integer f_getRealValue
      if( f_getRealValue( tag_scissor, dret, 'hartree') == 0) then
         scissor=dret
      else
         scissor=0.0d0
      end if
    end subroutine set_band_gap_correction
    
   end subroutine m_CtrlP_rd_LinearResponse
   

 end module m_LinearResponse_Control
