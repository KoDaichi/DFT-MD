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

  integer                   :: idp, idp2, idp3, nlp, nmp, nnp

  complex(kind=CMPLDP),private,target,allocatable,dimension(:)    :: cw1,cw2,cw3

  complex(kind=CMPLDP),private,allocatable,dimension(:)    :: wlp,wmp,wnp
  ! ------- Positron start 
  complex(kind=CMPLDP),private,target,allocatable,dimension(:)    :: cw1_pstrn,cw2_pstrn,cw3_pstrn

  real(kind=DP),private,allocatable,dimension(:)           :: ftw

  integer, dimension(2) :: flag_gfft_CD = (/1, -1/)   ! = 1 for INVERSE, = -1 for DIRECT
  integer, dimension(2) :: flag_jrcat_r3ft = (/-2,  2/) ! = -2 for INVERSE, = 2 for DIRECT

  integer :: sw_avoiding_odd_fftbox = ON
  integer :: sw_zero_padding = OFF

  integer, private      :: istat = 0

!  include 'mpif.h'

contains
  subroutine fft_WFCD_work_alloc
    integer ::nnl,nnm,nnn
! ----> work arrays of fft for wave functions
    if(kimg == 1) then
       if(ipri >= 1) write(nfout,'(" !!allocation of cw1,cw2,cw3 (fft_WFCD_work_alloc)")')
       nnl  = fft_box_size_WF(1,1)
       nnm  = fft_box_size_WF(2,1);  nnn  = fft_box_size_WF(3,1)
       allocate(cw1(nnl), stat=istat)
       allocate(cw2(nnm), stat=istat)
       allocate(cw3(nnn), stat=istat)
       if(ipri >= 1) write(nfout,'(" !! WF_JRCATFFTg <<m_FFT.fft_WFCD_work_alloc>>")')
    else if(kimg == 2) then
    end if
! ------- Positron start 
! ----> and work arrays of fft for positron wave functions
    if(sw_positron /= OFF) then
       if(kimg == 1) then
          if(ipri >= 1) write(nfout,'(" !!allocation of cw1_pstrn,cw2_pstrn,cw3_pstrn (fft_WFCD_work_alloc)")')
          nnl  = fft_box_size_pWF(1,1)
          nnm  = fft_box_size_pWF(2,1);  nnn  = fft_box_size_pWF(3,1)
          allocate(cw1_pstrn(nnl), stat=istat)
          allocate(cw2_pstrn(nnm), stat=istat)
          allocate(cw3_pstrn(nnn), stat=istat)
       else if(kimg == 2) then
       end if
    end if
! ------- Positron end

! ----> work arrays of ifft for the charge density
    if(kimg == 1) then
       nnl = fft_box_size_CD(1,1)
       nnm = fft_box_size_CD(2,1);  nnn = fft_box_size_CD(3,1)
       allocate(wlp(nnl), stat=istat)
       allocate(wmp(nnm), stat=istat)
       allocate(wnp(nnn), stat=istat)
       if(ipri >= 1) write(nfout,'(" !! wlp,wmp,wnp are allocated, CD_JRCATFFT_WS  <<m_FFT.fft_WFCD_work_alloc>>")')
    else
    end if
  end subroutine fft_WFCD_work_alloc

  subroutine m_FFT_alloc_WF_work
    integer :: nid,nnl,nnm,nnn,nmax
    if(kimg == 1) then
       nid  = fft_box_size_WF(1,0)
       nnm  = fft_box_size_WF(2,1);  nnn  = fft_box_size_WF(3,1)
       nmax = max(nid,nnm,nnn)
       allocate(ftw(nmax*2*2), stat=istat)     ! T. Yamasaki, 28 Jun. 2004
    else if(kimg == 2) then
       allocate(ftw(nfft), stat=istat)
    end if
  end subroutine m_FFT_alloc_WF_work

  subroutine m_FFT_alloc_pWF_work()
    integer :: nid,nnm,nnn, nmax
    if(kimg == 1) then
       nid  = fft_box_size_pWF(1,0)
       nnm  = fft_box_size_pWF(2,1);  nnn  = fft_box_size_pWF(3,1)
       nmax = max(nid,nnm,nnn)
       allocate(ftw(nmax*2*2), stat=istat)
    else if(kimg == 2) then
       allocate(ftw(nfft_pstrn), stat=istat)
    end if
  end subroutine m_FFT_alloc_pWF_work

  subroutine m_FFT_dealloc_WF_work
    if(allocated(ftw)) then
       deallocate(ftw, stat=istat)
       if(istat /= 0 ) then
          if(ipri>=1) then
             write(nfout,*) 'Deallocation error for ftw in sub. m_FFT_dealloc_WF_work'
             write(nfout,*) 'stat=', istat
          end if
          stop
       end if
    end if
  end subroutine m_FFT_dealloc_WF_work

  subroutine m_FFT_alloc_CD_box
    integer :: nmax
    if(kimg == 1) then
       nmax = max(idp,nmp,nnp)
       allocate(ftw(nmax*2*2), stat=istat)
    else if(kimg == 2) then
       allocate(ftw(nfftp_nonpara), stat=istat)
    end if
    allocate(afft_CD(nfftp_nonpara), stat=istat)
  end subroutine m_FFT_alloc_CD_box

  subroutine m_FFT_dealloc_CD_box
    if(allocated(ftw)) then
       deallocate(ftw, stat=istat)
    end if
    deallocate(afft_CD, stat=istat)
  end subroutine m_FFT_dealloc_CD_box

  subroutine m_FFT_setup(inversion_symmetry,paramset)
    integer, intent(in) :: inversion_symmetry
    logical, intent(in) :: paramset

    real(kind=DP),allocatable,dimension(:)          :: cfft

    integer :: nfft_t
    integer :: id_sname = -1
    call tstatc0_begin('m_FFT_setup ',id_sname)

    if(ipri >= 1) write(nfout,'(" !fft ncache = ",i5)') ncache
! --- fft_box_size_WF, fft_box_size_pWF ---
    if(inversion_symmetry == ON) then  ! kimg == 1
       fft_box_size_WF(1,0) = fft_box_size_WF(1,1) + 2
       fft_box_size_WF(2:3,0) = fft_box_size_WF(2:3,1)
       if(sw_positron /= OFF) then
          fft_box_size_pWF(1,0)   = fft_box_size_pWF(1,1) + 2
          fft_box_size_pWF(2:3,0) = fft_box_size_pWF(2:3,1)
       end if
    else if(inversion_symmetry == OFF) then ! kimg == 2

       if(ipri >= 1) write(nfout,*) '!FFT_WF = GOEDECKER FFT <<m_FFT_setup>>'
       fft_box_size_WF(1:3,0) = fft_box_size_WF(1:3,1) + 1
       if(sw_positron /= OFF) then
          fft_box_size_pWF(1:3,0) = fft_box_size_pWF(1:3,1) + 1
       end if
    endif

    nfft =   product(fft_box_size_WF(1:3,0)) * (2-inversion_symmetry)
! --- fft_box_size_CD ---
    if(inversion_symmetry == ON) then  ! kimg == 1
       fft_box_size_CD_nonpara(1,0) = fft_box_size_CD(1,1) + 2
       fft_box_size_CD_nonpara(2:3,0) = fft_box_size_CD(2:3,1)
    else if(inversion_symmetry == OFF) then ! kimg == 2
       if(ipri >= 1) write(nfout,*) '!FFT_CD = GOEDECKER FFT <<m_FFT_setup>>'
       fft_box_size_CD_nonpara(1:3,0) = fft_box_size_CD(1:3,1) + 1
    endif
    nfftp_nonpara  = product(fft_box_size_CD_nonpara(1:3,0)) * (2-inversion_symmetry)

    fft_box_size_CD(1:3,0) = fft_box_size_CD_nonpara(1:3,0)
    fft_box_size_CD_c(1:3,0) = fft_box_size_CD(1:3,0)
    nfftp  = product(fft_box_size_CD(1:3,0)) * (2-inversion_symmetry)
    nfftps = nfftp

    if(ipri >= 1) call wd_FFTboxsizes(nfout)

    idp = fft_box_size_CD_nonpara(1,0)
    idp2 = fft_box_size_CD_nonpara(2,0)
    idp3 = fft_box_size_CD_nonpara(3,0)
    nlp = fft_box_size_CD(1,1)
    nmp = fft_box_size_CD(2,1)
    nnp = fft_box_size_CD(3,1)

    if(sw_positron /= OFF) &
         & nfft_pstrn = product(fft_box_size_pWF(1:3,0))*(2-inversion_symmetry)

    if(.not. paramset) then
       if(kimg == 1) then
          call fft_WFCD_work_alloc     ! <cw[123],cw[123]_pstrn,wlp,wmp,wnp> are allocated

          if(sw_positron /= OFF .and. nfft_pstrn > nfft) then
             call m_FFT_alloc_pWF_work()  ! <ftw> is allocated
             nfft_t = nfft_pstrn
          else
             call m_FFT_alloc_WF_work()   ! <ftw> is allocated
             nfft_t = nfft
          end if
          allocate(cfft(nfft_t), stat=istat) ! cfft is used only for initiallization
          if(istat /= 0) then
             if(ipri >= 1) then
                write(nfout,*) 'Allocation error for cfft in sub. m_FFT_setup'
                write(nfout,*) 'stat =', istat, 'nfft =', nfft
             end if
             stop
          end if

          call init_fft_coefficients_arrays_WF()

          call m_FFT_dealloc_WF_work()
          deallocate(cfft,stat=istat)
          if(istat /= 0 ) then
             if(ipri >= 1) then
                write(nfout,*) 'Deallocation error for cfft in sub. m_FFT_setup'
                write(nfout,*) 'stat =', istat
             end if
             stop
          end if

          call CDFFT_setup()

       else if(kimg == 2) then
       end if

    endif

    call tstatc0_end(id_sname)
  contains
    subroutine init_fft_coefficients_arrays_WF()
      integer :: id, nl, nm, nn, ierr
      integer :: id_p, nl_p, nm_p, nn_p
      id = fft_box_size_WF(1,0)
      nl = fft_box_size_WF(1,1)
      nm = fft_box_size_WF(2,1)
      nn = fft_box_size_WF(3,1)
      if(sw_positron /= OFF) then
         id_p = fft_box_size_pWF(1,0)
         nl_p = fft_box_size_pWF(1,1)
         nm_p = fft_box_size_pWF(2,1)
         nn_p = fft_box_size_pWF(3,1)
      end if
      if(kimg == 1) then
         call jrcat_r3ft(cfft,ftw,id,nl,nm,nn,cw1,cw2,cw3,0,0,0,0)
         if(sw_positron /= OFF) &
              & call jrcat_r3ft(cfft,ftw,id_p,nl_p,nm_p,nn_p &
              &                   ,cw1_pstrn,cw2_pstrn,cw3_pstrn,0,0,0,0)
      else  ! kimg == 2
      endif

    end subroutine init_fft_coefficients_arrays_WF
  end subroutine m_FFT_setup

  subroutine CDFFT_setup()
    if(kimg == 1) then
       call m_FFT_alloc_CD_box()       ! <ftw> is allocated
       call init_fft_coefficients_arrays_CD()
       call m_FFT_dealloc_CD_box()
    end if

    CD_setup_is_done = YES
  contains
    subroutine init_fft_coefficients_arrays_CD
      integer :: id, nl, nm, nn, ierr
   !  ---> FFT for Charge density
      call jrcat_r3ft(afft_CD,ftw,idp,nlp,nmp,nnp,wlp,wmp,wnp,0,0,0,0)

    end subroutine init_fft_coefficients_arrays_CD
  end subroutine CDFFT_setup

  subroutine m_FFT_WF(electron_or_positron,nfout,afft,inverse_or_direct,switch)  ! G space --> R space
    integer, intent(in)          :: electron_or_positron
    integer, intent(in)          :: nfout
    real(kind=DP), intent(inout) :: afft(nfft)
    integer, intent(in)          :: inverse_or_direct
    integer, intent(in)          :: switch

    complex(kind=CMPLDP),pointer,dimension(:) :: cw1_t, cw2_t, cw3_t
    integer :: id, nl, nm, nn
    integer :: nd2, nd3, inzee, isign

!!$    integer, dimension(2) :: flag_jrcat_r3ft = (/-2,  2/)
    integer kc1, kc2, kc3

    integer :: id_sname = -1
    if(electron_or_positron == ELECTRON) then
       call tstatc0_begin('m_FFT_WF ',id_sname)

       id = fft_box_size_WF(1,0)
       nl = fft_box_size_WF(1,1)
       nm = fft_box_size_WF(2,1)
       nn = fft_box_size_WF(3,1)
       if(kimg == 1) then
          cw1_t=>cw1; cw2_t=>cw2; cw3_t=>cw3 
       else if(kimg == 2) then
          nd2 = fft_box_size_WF(2,0)
          nd3 = fft_box_size_WF(3,0)
       end if
    else if(electron_or_positron == POSITRON) then
       call tstatc0_begin('m_FFT_pWF ',id_sname)

       id = fft_box_size_pWF(1,0)
       nl = fft_box_size_pWF(1,1)
       nm = fft_box_size_pWF(2,1)
       nn = fft_box_size_pWF(3,1)
       if(kimg==1) then
          cw1_t=>cw1_pstrn; cw2_t=>cw2_pstrn; cw3_t=>cw3_pstrn
       else if(kimg==2) then
          nd2 = fft_box_size_pWF(2,0)
          nd3 = fft_box_size_pWF(3,0)
       end if
    end if

    if(kimg == 1) then
       if(switch == ON) then
          kc1 = 0;       kc2 = nm/4 + 1;       kc3 = nn/4 + 1
       else
          kc1 = 0;       kc2 = 0       ;       kc3 = 0
       endif
       call jrcat_r3ft(afft,ftw,id,nl,nm,nn,cw1_t,cw2_t,cw3_t,kc1,kc2,kc3&
            & ,flag_jrcat_r3ft(inverse_or_direct))
    else if(kimg == 2) then
       inzee = 1
       if(inverse_or_direct == INVERSE) then     ! INVERSE = 1 
          isign = -1
       else if(inverse_or_direct == DIRECT) then ! DIRECT = 2
          isign = 1
       else
          if(ipri>=1) write(*,*) " isign error <<m_FFT_WF>>, GOEDECKER_FFT"
          stop
       end if
       call gfft(nl,nm,nn,id,nd2,nd3,afft,ftw,isign,inzee,ncache)
       if(inzee == 2) afft = ftw

    endif

    call tstatc0_end(id_sname)

  end subroutine m_FFT_WF

  subroutine fft_CD_inverse_core(afft_CD)
    real(kind=DP),intent(inout),dimension(nfftp_nonpara) :: afft_CD
    integer :: kc1 = 0, kc2 = 0, kc3 = 0
    integer :: nsize_ftw
    integer :: inzee

    nsize_ftw = 0
    if(.not.allocated(ftw)) then
       if(kimg==1) then
          nsize_ftw = max(idp,nmp,nnp)*2*2
       else
          nsize_ftw = nfftp
       end if
       allocate(ftw(nsize_ftw))
    end if

    if(kimg == 1) then
       call jrcat_r3ft(afft_CD,ftw,idp,nlp,nmp,nnp,wlp,wmp,wnp,kc1,kc2,kc3&
            & ,flag_jrcat_r3ft(INVERSE))
    else if(kimg == 2) then
       inzee = 1
       ftw = 0.d0
       call gfft(nlp,nmp,nnp,idp,idp2,idp3,afft_CD,ftw,flag_gfft_CD(INVERSE),inzee,ncache)
       if(inzee == 2) afft_CD = ftw
    endif
    if(nsize_ftw /= 0) deallocate(ftw)
  end subroutine fft_CD_inverse_core

  subroutine fft_CD_direct_core(afft_CD)
    real(kind=DP),intent(inout),dimension(nfftp_nonpara) :: afft_CD
    integer :: nsize_ftw
    integer :: kc1 = 0, kc2 = 0, kc3 = 0
    integer :: inzee

    nsize_ftw = 0
    if(.not.allocated(ftw)) then
       if(kimg==1) then
          nsize_ftw = max(idp,nmp,nnp)*2*2
       else
          nsize_ftw = nfftp
       end if
       allocate(ftw(nsize_ftw))
    end if

    if(kimg == 1) then
       call jrcat_r3ft(afft_CD,ftw,idp,nlp,nmp,nnp,wlp,wmp,wnp,kc1,kc2,kc3&
            & , flag_jrcat_r3ft(DIRECT))
    else if(kimg == 2) then
       inzee = 1
       ftw = 0.d0
       call gfft(nlp,nmp,nnp,idp,idp2,idp3,afft_CD,ftw,flag_gfft_CD(DIRECT),inzee,ncache)
       if(inzee == 2) afft_CD = ftw
    endif

    if(nsize_ftw /= 0) deallocate(ftw)
  end subroutine fft_CD_direct_core

! ------------------------
