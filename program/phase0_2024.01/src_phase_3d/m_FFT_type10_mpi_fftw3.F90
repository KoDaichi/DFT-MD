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

  ! ---------------------
  integer(8) :: plan_WF(2),plan_CD(2)

  integer(8) :: plan_CD_exx(2)

  ! ------- Positron start 
  integer(8) :: plan_pWF(2)
  ! ------- Positron end

! ======= KT_add ============================== 13.0F
  integer(8) :: plan_exx(2)
! ============================================= 13.0F

  integer :: sw_avoiding_odd_fftbox = OFF
  integer :: sw_zero_padding = OFF

  include 'fftw3-mpi.f03'
  type(C_PTR) :: plan_distfftw_inverse_wf
  type(C_PTR) :: plan_distfftw_direct_wf
  logical :: distfftw_inverse_wf_flg = .false.
  logical :: distfftw_direct_wf_flg  = .false.

  complex(kind=DP), pointer :: afft_mpifftw(:,:,:)=>null()
  real(kind=DP),    pointer :: afft_mpifftw_vlocal(:,:,:)=>null()
  complex(kind=DP), pointer :: bfft_mpifftw(:,:,:)=>null()

  complex(kind=DP), pointer :: afft_mpifftw_kimg1(:,:,:)=>null()
  real(kind=DP),    pointer :: bfft_mpifftw_kimg1(:,:,:)=>null()

!  include 'mpif.h'

contains
  subroutine m_FFT_alloc_WF_work
  end subroutine m_FFT_alloc_WF_work

  subroutine m_FFT_alloc_pWF_work()
  end subroutine m_FFT_alloc_pWF_work

  subroutine m_FFT_dealloc_WF_work
  end subroutine m_FFT_dealloc_WF_work

! ======= KT_add ============================= 13.0F
  subroutine m_FFT_alloc_exx_work
  end subroutine m_FFT_alloc_exx_work

  subroutine m_FFT_dealloc_exx_work
  end subroutine m_FFT_dealloc_exx_work
! =========================================== 13.0F

  subroutine m_FFT_alloc_CD_box
    integer :: istat
    if(sw_mpifft == OFF) then
       allocate(afft_CD(nfftp), stat=istat)
    end if
  end subroutine m_FFT_alloc_CD_box

  subroutine m_FFT_dealloc_CD_box
    integer :: istat
    if(sw_mpifft == OFF) then
       deallocate(afft_CD)
    end if
  end subroutine m_FFT_dealloc_CD_box

  subroutine m_FFT_finalize()
    call dfftw_destroy_plan(plan_WF(1))    
    call dfftw_destroy_plan(plan_WF(2))    
!    call dfftw_destroy_plan(plan_CD(1))    
!    call dfftw_destroy_plan(plan_CD(2))    
    if(sw_positron/=OFF)then
       call dfftw_destroy_plan(plan_pWF(1))
       call dfftw_destroy_plan(plan_pWF(2))
    endif

!    call dfftw_cleanup()
  end subroutine m_FFT_finalize

  subroutine m_FFT_setup(inversion_symmetry,paramset)
    integer, intent(in) :: inversion_symmetry
    logical, intent(in) :: paramset

    real(kind=DP), allocatable, dimension(:)          :: ftw

    integer :: istat = 0
    integer :: nfft_t, ipad
    integer :: id_sname = -1
    call tstatc0_begin('m_FFT_setup ',id_sname)

    if(inversion_symmetry == ON) then  ! kimg == 1
       ipad = 2
    else if(inversion_symmetry == OFF) then ! kimg == 2
       ipad = 0
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

!$$!!#ifdef PARA3D
       fft_box_size_CD_3D(:,1) = fft_box_size_CD(:,1)
       fft_box_size_CD_3D(:,0) = fft_box_size_CD_nonpara(:,0)
       nfftps = product(fft_box_size_CD_3D(1:3,0)) * (2-inversion_symmetry)
!$$!!#endif

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
       nfft_t = nfft
       if(sw_positron /= OFF .and. nfft_t < nfft_pstrn) nfft_t = nfft_pstrn
       allocate(ftw(nfft_t), stat=istat) ! ftw is used only for initiallization
       if(istat /= 0) then
          if(ipri >= 1) then
             write(nfout,*) 'Allocation error for ftw in sub. m_FFT_setup'
             write(nfout,*) 'stat =', istat, 'nfft =', nfft
          end if
          stop
       end if

       call init_fft_coefficients_arrays_WF()

       deallocate(ftw,stat=istat)
       if(istat /= 0 ) then
          if(ipri >= 1) then
             write(nfout,*) 'Deallocation error for ftw in sub. m_FFT_setup'
             write(nfout,*) 'stat =', istat
          end if
          stop
       end if
!!$#ifndef _MPIFFT_
       if(sw_mpifft == OFF) then
! Initialization of the Charge-Density FFT
          call CDFFT_setup()
!!$#endif
       end if
    end if

    call tstatc0_end(id_sname)
  contains
    subroutine init_fft_coefficients_arrays_WF()
      integer :: nl, nm, nn
      integer :: nl_p, nm_p, nn_p
      nl = fft_box_size_WF(1,1)
      nm = fft_box_size_WF(2,1)
      nn = fft_box_size_WF(3,1)
      if(sw_positron /= OFF) then
         nl_p = fft_box_size_pWF(1,1)
         nm_p = fft_box_size_pWF(2,1)
         nn_p = fft_box_size_pWF(3,1)
      end if

      if(kimg == 1) then
         ! Forward FFT
         call dfftw_plan_dft_c2r_3d(plan_WF(1) &
       &                       ,nl,nm,nn &
       &                       ,ftw(1),ftw(1) &
       &                       ,FFTW_MEASURE)
         ! Inverse FFT
         call dfftw_plan_dft_r2c_3d(plan_WF(2) &
       &                       ,nl,nm,nn &
       &                       ,ftw(1),ftw(1) &
       &                       ,FFTW_MEASURE)
         if(sw_positron /= OFF) then
            ! Forward FFT
            call dfftw_plan_dft_c2r_3d(plan_pWF(1) &
       &                       ,nl_p,nm_p,nn_p &
       &                       ,ftw(1),ftw(1) &
       &                       ,FFTW_MEASURE)
            ! Inverse FFT
            call dfftw_plan_dft_r2c_3d(plan_pWF(2) &
       &                       ,nl_p,nm_p,nn_p &
       &                       ,ftw(1),ftw(1) &
       &                       ,FFTW_MEASURE)
         end if
      else
         ! Forward FFT
         call dfftw_plan_dft_3d(plan_WF(1),nl,nm,nn &
     &                         ,ftw(1),ftw(1) &
     &                         ,-1,FFTW_MEASURE)
         ! Inverse FFT
         call dfftw_plan_dft_3d(plan_WF(2),nl,nm,nn &
     &                         ,ftw(1),ftw(1) &
     &                         ,+1,FFTW_MEASURE)
         if(sw_positron /= OFF) then
            ! Forward FFT
            call dfftw_plan_dft_3d(plan_pWF(1),nl_p,nm_p,nn_p &
     &                         ,ftw(1),ftw(1) &
     &                         ,-1,FFTW_MEASURE)
            ! Inverse FFT
            call dfftw_plan_dft_3d(plan_pWF(2),nl_p,nm_p,nn_p &
     &                         ,ftw(1),ftw(1) &
     &                         ,+1,FFTW_MEASURE)
         end if
      end if
      if(ipri >= 1) then
         write(nfout,'(" !(init_fft_coef_WF) nl, nm, nn   = ",3i8)') nl, nm, nn
         write(nfout,'(" !(init_fft_coef_WF) plan_WF(1:2) = ",2i20)') plan_WF(1:2)
      end if
    end subroutine init_fft_coefficients_arrays_WF


  end subroutine m_FFT_setup

! ======= KT_add =========================================== 13.0F
  subroutine m_FFT_setup_exx(inversion_symmetry,paramset)
    integer, intent(in) :: inversion_symmetry
    logical, intent(in) :: paramset

    real(kind=DP), allocatable, dimension(:)          :: ftw

    integer :: istat = 0
    integer :: nfft_t, ipad

    integer :: id_sname = -1
    call tstatc0_begin('m_FFT_setup_exx ',id_sname)

    if(inversion_symmetry == ON) then  ! kimg == 1
       ipad = 2
    else if(inversion_symmetry == OFF) then ! kimg == 2
       ipad = 0
    end if

    fft_box_size_exx(1,0)   = fft_box_size_exx(1,1) + ipad
    fft_box_size_exx(2:3,0) = fft_box_size_exx(2:3,1)

    nfft_exx = product(fft_box_size_exx(1:3,0)) * (2-inversion_symmetry)

    if(ipri >= 1) call wd_FFTboxsizes_exx(nfout)

    if(.not. paramset) then
! Initialization of the Wave-Function FFT
       nfft_t = nfft_exx

       allocate(ftw(nfft_t), stat=istat) ! ftw is used only for initiallization
       if(istat /= 0) then
          if(ipri >= 1) then
             write(nfout,*) 'Allocation error for ftw in sub. m_FFT_setup_exx'
             write(nfout,*) 'stat =', istat, 'nfft_exx =', nfft_exx
          end if
          stop
       end if

       call init_fft_coefficients_arrays_WF()

       deallocate(ftw,stat=istat)
       if(istat /= 0 ) then
          if(ipri >= 1) then
             write(nfout,*) 'Deallocation error for ftw in sub. m_FFT_setup_exx'
             write(nfout,*) 'stat =', istat
          end if
          stop
       end if
    end if

    call tstatc0_end(id_sname)

  contains

    subroutine init_fft_coefficients_arrays_WF()
      integer :: nl_exx, nm_exx, nn_exx

      nl_exx = fft_box_size_exx(1,1)
      nm_exx = fft_box_size_exx(2,1)
      nn_exx = fft_box_size_exx(3,1)

      if(kimg == 1) then
         ! Forward FFT
         call dfftw_plan_dft_c2r_3d( plan_exx(1),nl_exx,nm_exx, nn_exx &
              &                     ,ftw(1),ftw(1) &
              &                     ,FFTW_MEASURE)
         ! Inverse FFT
         call dfftw_plan_dft_r2c_3d( plan_exx(2),nl_exx, nm_exx,nn_exx &
              &                     ,ftw(1),ftw(1) &
              &                     ,FFTW_MEASURE)
      else
         ! Forward FFT
         call dfftw_plan_dft_3d( plan_exx(1),nl_exx,nm_exx, nn_exx &
              &                         ,ftw(1),ftw(1) &
              &                         ,-1,FFTW_MEASURE)
         ! Inverse FFT
         call dfftw_plan_dft_3d( plan_exx(2),nl_exx, nm_exx,nn_exx &
              &                         ,ftw(1),ftw(1) &
              &                         ,+1,FFTW_MEASURE)
      end if

      if(ipri >= 1) then
         write(nfout,'(" !(init_fft_coef_exx) nl_exx, nm_exx, nn_exx   = ",3i8)') &
              &              nl_exx, nm_exx, nn_exx
         write(nfout,'(" !(init_fft_coef_exx) plan_exx(1:2) = ",2i20)') &
              &             plan_exx(1:2)
      end if
    end subroutine init_fft_coefficients_arrays_WF

  end subroutine m_FFT_setup_exx
! ======================================================== 13.0F

  subroutine CDFFT_setup_exx()
    allocate(afft_CD_exx(nfftp_exx_nonpara))
    call init_fft_coefficients_arrays_CDx()
    CD_setup_is_done_exx = YES
    deallocate(afft_CD_exx)
  contains
    subroutine init_fft_coefficients_arrays_CDx
      integer :: nl, nm, nn

   !  ---> FFT for Charge density
      nl = fft_box_size_CD_exx(1,1)
      nm = fft_box_size_CD_exx(2,1)
      nn = fft_box_size_CD_exx(3,1)
      if(kimg == 1) then
         ! Forward FFT
         call dfftw_plan_dft_c2r_3d(plan_CD_exx(1) &
       &                       ,nl,nm,nn &
       &                       ,afft_CD_exx(1),afft_CD_exx(1) &
       &                       ,FFTW_MEASURE)
         ! Inverse FFT
         call dfftw_plan_dft_r2c_3d(plan_CD_exx(2) &
       &                       ,nl,nm,nn &
       &                       ,afft_CD_exx(1),afft_CD_exx(1) &
       &                       ,FFTW_MEASURE)
      else
         ! Forward FFT
         call dfftw_plan_dft_3d(plan_CD_exx(1),nl,nm,nn &
     &                         ,afft_CD_exx(1),afft_CD_exx(1) &
     &                         ,-1,FFTW_MEASURE)
         ! Inverse FFT
         call dfftw_plan_dft_3d(plan_CD_exx(2),nl,nm,nn &
     &                         ,afft_CD_exx(1),afft_CD_exx(1) &
     &                         ,+1,FFTW_MEASURE)
      end if
      if(ipri >= 1) then
         write(nfout,'(" !(CDFFT_setup_exx) nl, nm, nn   = ",3i8)') nl, nm, nn
         write(nfout,'(" !(CDFFT_setup_exx) plan_CD_exx(1:2) = ",2i20)') plan_CD_exx(1:2)
      end if
    end subroutine init_fft_coefficients_arrays_CDx
  end subroutine CDFFT_setup_exx

  subroutine CDFFT_setup()
    allocate(afft_CD(nfftp_nonpara))
    call init_fft_coefficients_arrays_CD()
    CD_setup_is_done = YES
    deallocate(afft_CD)
  contains
    subroutine init_fft_coefficients_arrays_CD
      integer :: nl, nm, nn

   !  ---> FFT for Charge density
      nl = fft_box_size_CD(1,1)
      nm = fft_box_size_CD(2,1)
      nn = fft_box_size_CD(3,1)
      if(kimg == 1) then
         ! Forward FFT
         call dfftw_plan_dft_c2r_3d(plan_CD(1) &
       &                       ,nl,nm,nn &
       &                       ,afft_CD(1),afft_CD(1) &
       &                       ,FFTW_MEASURE)
         ! Inverse FFT
         call dfftw_plan_dft_r2c_3d(plan_CD(2) &
       &                       ,nl,nm,nn &
       &                       ,afft_CD(1),afft_CD(1) &
       &                       ,FFTW_MEASURE)
      else
         ! Forward FFT
         call dfftw_plan_dft_3d(plan_CD(1),nl,nm,nn &
     &                         ,afft_CD(1),afft_CD(1) &
     &                         ,-1,FFTW_MEASURE)
         ! Inverse FFT
         call dfftw_plan_dft_3d(plan_CD(2),nl,nm,nn &
     &                         ,afft_CD(1),afft_CD(1) &
     &                         ,+1,FFTW_MEASURE)
      end if
      if(ipri >= 1) then
         write(nfout,'(" !(CDFFT_setup) nl, nm, nn   = ",3i8)') nl, nm, nn
         write(nfout,'(" !(CDFFT_setup) plan_CD(1:2) = ",2i20)') plan_CD(1:2)
      end if
    end subroutine init_fft_coefficients_arrays_CD
  end subroutine CDFFT_setup

! === necessary to make 3D_Parallel, too!!! by tkato ===========================
!!BRANCH_P ORG_Parallel
! ==============================================================================
!$$#ifndef PARA3D
  subroutine m_FFT_WF(electron_or_positron,nfout,afft,inverse_or_direct,switch)  ! G space --> R space
    integer, intent(in)          :: electron_or_positron
    integer, intent(in)          :: nfout
    real(kind=DP), intent(inout) :: afft(nfft)
    integer, intent(in)          :: inverse_or_direct
    integer, intent(in)          :: switch

    integer(8) :: plan(2)
    integer :: id_sname = -1
    if(electron_or_positron == ELECTRON) then
       call tstatc0_begin('m_FFT_WF ',id_sname)
    else if(electron_or_positron == POSITRON) then
       call tstatc0_begin('m_FFT_pWF ',id_sname)
    end if

    if(electron_or_positron == ELECTRON) then
       plan(1:2) = plan_WF(1:2)
    else if(electron_or_positron == POSITRON) then
       plan(1:2) = plan_pWF(1:2)
    end if

    if(inverse_or_direct == DIRECT) then
       if(kimg==1) then
          call dfftw_execute_dft_c2r(plan(1),afft(1),afft(1))
       else
          call dfftw_execute_dft(plan(1),afft(1),afft(1))
       end if
    else ! INVERSE
       if(kimg==1) then
          call dfftw_execute_dft_r2c(plan(2),afft(1),afft(1))
       else
          call dfftw_execute_dft(plan(2),afft(1),afft(1))
       end if
    end if

    call tstatc0_end(id_sname)
  end subroutine m_FFT_WF
!$$#endif
! === necessary to make 3D_Parallel, too!!! by tkato ===========================
!!BRANCH_P_END ORG_Parallel
! ==============================================================================

! ========== KT_add ========================================= 13.0F
  subroutine m_FFT_exx( nfout,afft,inverse_or_direct,switch)  ! G space --> R space
    integer, intent(in)          :: nfout
    real(kind=DP), intent(inout) :: afft(nfft_exx)
    integer, intent(in)          :: inverse_or_direct
    integer, intent(in)          :: switch

    integer(8) :: plan(2)
    integer :: id_sname = -1

    call tstatc0_begin('m_FFT_exx ',id_sname)

    plan(1:2) = plan_exx(1:2)

    if(inverse_or_direct == DIRECT) then
       if(kimg==1) then
          call dfftw_execute_dft_c2r(plan(1),afft(1),afft(1))
       else
          call dfftw_execute_dft(plan(1),afft(1),afft(1))
       endif
    else ! INVERSE
       if(kimg==1) then
          call dfftw_execute_dft_r2c(plan(2),afft(1),afft(1))
       else
          call dfftw_execute_dft(plan(2),afft(1),afft(1))
       end if
    end if

    call tstatc0_end(id_sname)
  end subroutine m_FFT_exx
! ============================================================ 13.0F

  subroutine fft_CD_inverse_core(afft_CD)
    real(kind=DP),intent(inout),dimension(nfftp_nonpara) :: afft_CD

    if(kimg==1) then
       call dfftw_execute_dft_r2c(plan_CD(2),afft_CD(1),afft_CD(1))
    else
       call dfftw_execute_dft(plan_CD(2),afft_CD(1),afft_CD(1))
    endif
  end subroutine fft_CD_inverse_core

  subroutine fft_CD_direct_core(afft_CD)
    real(kind=DP),intent(inout),dimension(nfftp_nonpara) :: afft_CD
    if(kimg==1) then
       call dfftw_execute_dft_c2r(plan_CD(1),afft_CD(1),afft_CD(1))
    else
       call dfftw_execute_dft(plan_CD(1),afft_CD(1),afft_CD(1))
    end if
  end subroutine fft_CD_direct_core

  subroutine fft_CD_inverse_core_exx(afft_CD)
    real(kind=DP),intent(inout),dimension(nfftp_exx_nonpara) :: afft_CD

    if(kimg==1) then
       call dfftw_execute_dft_r2c(plan_CD_exx(2),afft_CD(1),afft_CD(1))
    else
       call dfftw_execute_dft(plan_CD_exx(2),afft_CD(1),afft_CD(1))
    endif
  end subroutine fft_CD_inverse_core_exx

  subroutine fft_CD_direct_core_exx(afft_CD)
    real(kind=DP),intent(inout),dimension(nfftp_exx_nonpara) :: afft_CD
    if(kimg==1) then
       call dfftw_execute_dft_c2r(plan_CD_exx(1),afft_CD(1),afft_CD(1))
    else
       call dfftw_execute_dft(plan_CD_exx(1),afft_CD(1),afft_CD(1))
    end if
  end subroutine fft_CD_direct_core_exx

! --------------
