!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  MODULE: m_FFT
!
!  AUTHOR(S): T. Yamasaki, K. Betsuyaku,   August/20/2003
!  
!  FURTHER MODIFICATION: T. Yamasaki, January/13/2004
!                                   , March/14/2005, April/10/2007
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

  type(DFTI_DESCRIPTOR), POINTER :: dh_WF, dh_CD 
  type(mklfft_rdh) :: rdh_WF, rdh_CD
  integer :: strides(4), Status

  ! ------- Positron start 
  type(DFTI_DESCRIPTOR), POINTER :: dh_pWF
  type(mklfft_rdh) :: rdh_pWF
  ! ------- Positron end

  integer :: sw_avoiding_odd_fftbox = OFF
  integer :: sw_zero_padding = OFF

!  include 'mpif.h'
contains

  subroutine m_FFT_alloc_WF_work
  end subroutine m_FFT_alloc_WF_work

  subroutine m_FFT_alloc_pWF_work()
  end subroutine m_FFT_alloc_pWF_work

  subroutine m_FFT_dealloc_WF_work
  end subroutine m_FFT_dealloc_WF_work

  subroutine m_FFT_alloc_CD_box
    integer :: istat
    if(sw_mpifft == OFF) then
       allocate(afft_CD(nfftp), stat=istat)
    end if
  end subroutine m_FFT_alloc_CD_box

  subroutine m_FFT_dealloc_CD_box
    integer :: istat
    if(sw_mpifft == OFF) then
       deallocate(afft_CD, stat=istat)
    end if
  end subroutine m_FFT_dealloc_CD_box

  subroutine m_FFT_setup(inversion_symmetry,paramset)
    integer, intent(in) :: inversion_symmetry
    logical, intent(in) :: paramset

    integer :: ipad
    integer :: id_sname = -1
    call tstatc0_begin('m_FFT_setup ',id_sname)

    if(inversion_symmetry == ON) then  ! kimg == 1
       ipad = 2
    else if(inversion_symmetry == OFF) then ! kimg == 2
       ipad = 1
    end if
! --- fft_box_size_WF, fft_box_size_pWF ---
    fft_box_size_WF(1,0) = fft_box_size_WF(1,1) + ipad
    fft_box_size_WF(2:3,0) = fft_box_size_WF(2:3,1)
    if(sw_positron /= OFF) then
       fft_box_size_pWF(1,0)   = fft_box_size_pWF(1,1) + ipad
       fft_box_size_pWF(2:3,0) = fft_box_size_pWF(2:3,1)
    end if

    nfft =   product(fft_box_size_WF(1:3,0)) * (2-inversion_symmetry)
! --- fft_box_size_CD ---
    fft_box_size_CD_nonpara(1,0)   = fft_box_size_CD(1,1) + ipad
    fft_box_size_CD_nonpara(2:3,0) = fft_box_size_CD(2:3,1)
    nfftp_nonpara  = product(fft_box_size_CD_nonpara(1:3,0)) * (2-inversion_symmetry)

    if(sw_mpifft == ON) then
       if(ipri >= 1) write(nfout,*) '!FFT_CD = MPIFFT <<m_FFT_setup>>'
       call set_mpifft_box_size_CD(inversion_symmetry) ! -> fft_box_size_CD, fft_box_size_CD_c, ny_d, nz_d, etc.
       nfftp  = product(fft_box_size_CD_c(1:3,0)) * (2-inversion_symmetry)
       nfftps = product(fft_box_size_CD(1:3,0)) * (2-inversion_symmetry)
    else
       fft_box_size_CD(1:3,0) = fft_box_size_CD_nonpara(1:3,0)
       fft_box_size_CD_c(1:3,0) = fft_box_size_CD(1:3,0)
       nfftp  = product(fft_box_size_CD(1:3,0)) * (2-inversion_symmetry)
       nfftps = nfftp
    end if

    if(ipri >= 1) call wd_FFTboxsizes(nfout)

    if(sw_positron /= OFF) &
         & nfft_pstrn = product(fft_box_size_pWF(1:3,0))*(2-inversion_symmetry)

    if(.not. paramset) then
! Initialization of the Wave-Function FFT
       call init_fft_coefficients_arrays_WF()

       if(sw_mpifft == OFF) then
! Initialization of the Charge-Density FFT
          call CDFFT_setup()
       end if
    endif

    call tstatc0_end(id_sname)
  contains
    subroutine init_fft_coefficients_arrays_WF()
      integer :: id, id_p

      if(kimg == 1) then
         call init_mkl_rfft(fft_box_size_WF(1,1) &
                         & ,fft_box_size_WF(2,1) &
                         & ,fft_box_size_WF(3,1) &
                         & ,fft_box_size_WF(1,0) &
                         & ,fft_box_size_WF(2,0) &
                         & ,fft_box_size_WF(3,0) &
                         & ,rdh_WF)
      else
         id = fft_box_size_WF(1,0)
         Status = DftiCreateDescriptor( dh_WF, DFTI_DOUBLE, &
                             & DFTI_COMPLEX, 3, fft_box_size_WF(1,1))
         strides(1) = 0
         strides(2) = 1
         strides(3) = id
         strides(4) = id*fft_box_size_WF(2,0)
         Status = DftiSetValue(dh_WF, DFTI_INPUT_STRIDES, strides)
         Status = DftiCommitDescriptor( dh_WF )
      end if
      if(sw_positron /= OFF) then
         if(kimg == 1) then
            call init_mkl_rfft(fft_box_size_pWF(1,1) &
                            & ,fft_box_size_pWF(2,1) &
                            & ,fft_box_size_pWF(3,1) &
                            & ,fft_box_size_pWF(1,0) &
                            & ,fft_box_size_pWF(2,0) &
                            & ,fft_box_size_pWF(3,0) &
                            & ,rdh_pWF)
         else
            id_p = fft_box_size_pWF(1,0)
            Status = DftiCreateDescriptor( dh_pWF, DFTI_DOUBLE, &
                     & DFTI_COMPLEX, 3, fft_box_size_pWF(1,1))
            strides(1) = 0
            strides(2) = 1
            strides(3) = id_p
            strides(4) = id_p*fft_box_size_pWF(2,0)
            Status = DftiSetValue(dh_pWF, DFTI_INPUT_STRIDES, strides)
            Status = DftiCommitDescriptor( dh_pWF )
         end if
      end if
    end subroutine init_fft_coefficients_arrays_WF

  end subroutine m_FFT_setup

  subroutine CDFFT_setup()

    call init_fft_coefficients_arrays_CD()
    CD_setup_is_done = YES
  contains
    subroutine init_fft_coefficients_arrays_CD
      integer :: id

   !  ---> FFT for Charge density
      if(kimg == 1) then
         call init_mkl_rfft(fft_box_size_CD(1,1) &
                         & ,fft_box_size_CD(2,1) &
                         & ,fft_box_size_CD(3,1) &
                         & ,fft_box_size_CD_nonpara(1,0) &
                         & ,fft_box_size_CD_nonpara(2,0) &
                         & ,fft_box_size_CD_nonpara(3,0) &
                         & ,rdh_CD)
      else
         id = fft_box_size_CD_nonpara(1,0)
         Status = DftiCreateDescriptor( dh_CD, DFTI_DOUBLE, &
                             & DFTI_COMPLEX, 3, fft_box_size_CD(1,1))
         strides(1) = 0
         strides(2) = 1
         strides(3) = id
         strides(4) = id*fft_box_size_CD_nonpara(2,0)
         Status = DftiSetValue(dh_CD, DFTI_INPUT_STRIDES, strides)
         Status = DftiCommitDescriptor( dh_CD )
      end if
    end subroutine init_fft_coefficients_arrays_CD
  end subroutine CDFFT_setup

  subroutine m_FFT_WF(electron_or_positron,nfout,afft,inverse_or_direct,switch)  ! G space --> R space
    integer, intent(in)          :: electron_or_positron
    integer, intent(in)          :: nfout
    real(kind=DP), intent(inout) :: afft(nfft)
    integer, intent(in)          :: inverse_or_direct
    integer, intent(in)          :: switch

    integer :: id_sname = -1

    if(electron_or_positron == ELECTRON) then
       call tstatc0_begin('m_FFT_WF ',id_sname)
    else if(electron_or_positron == POSITRON) then
       call tstatc0_begin('m_FFT_pWF ',id_sname)
    end if

    if(electron_or_positron == ELECTRON) then
       if(inverse_or_direct == DIRECT) then
          if(kimg == 1) then
             call mkl_rfft(1,rdh_WF,afft)
          else
             Status = DftiComputeForward( dh_WF, afft)
          end if
       else
          if(kimg == 1) then
             call mkl_rfft(-1,rdh_WF,afft)
          else
             Status = DftiComputeBackward( dh_WF, afft)
          end if
       end if
    else if(electron_or_positron == POSITRON) then
       if(inverse_or_direct == DIRECT) then
          if(kimg == 1) then
             call mkl_rfft(1,rdh_pWF,afft)
          else
             Status = DftiComputeForward( dh_pWF, afft)
          end if
       else
          if(kimg == 1) then
             call mkl_rfft(-1,rdh_pWF,afft)
          else
             Status = DftiComputeBackward( dh_pWF, afft)
          end if
       end if
    end if

    call tstatc0_end(id_sname)
  end subroutine m_FFT_WF

  subroutine fft_CD_inverse_core(afft_CD)
    real(kind=DP),intent(inout),dimension(nfftp_nonpara) :: afft_CD

    if(kimg == 1) then
       call mkl_rfft(-1,rdh_CD,afft_CD)
    else
       Status = DftiComputeBackward( dh_CD, afft_CD)
    end if
  end subroutine fft_CD_inverse_core

  subroutine fft_CD_direct_core(afft_CD)
    real(kind=DP),intent(inout),dimension(nfftp_nonpara) :: afft_CD
    if(kimg == 1) then
       call mkl_rfft(1,rdh_CD,afft_CD)
    else
       Status = DftiComputeForward( dh_CD, afft_CD)
    end if
  end subroutine fft_CD_direct_core
  ! --------------------------
