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


  complex(kind=CMPLDP),private,target,allocatable,dimension(:)    :: cw1,cw2,cw3
  complex(kind=CMPLDP),private,allocatable,dimension(:)    :: wlp,wmp,wnp
  complex(kind=CMPLDP),private,allocatable,dimension(:)    :: wlp_exx,wmp_exx,wnp_exx

  ! ------- Positron start 
  complex(kind=CMPLDP),private,target,allocatable,dimension(:)    :: cw1_pstrn,cw2_pstrn,cw3_pstrn
  ! ------- Positron end
  real(kind=DP),private,allocatable,dimension(:)           :: ftw

  integer :: idp, nlp, nmp, nnp
  integer :: idp_exx, nlp_exx, nmp_exx, nnp_exx
  integer :: istat

! ======= KT_add ========================================= 13.0F
  complex(kind=CMPLDP),private,target,allocatable,dimension(:) :: &
       &                        cw1_exx, cw2_exx, cw3_exx
! ======================================================== 13.0F

  ! ----------------
  integer :: sw_avoiding_odd_fftbox = OFF
  integer :: sw_zero_padding = OFF

!  include 'mpif.h'
contains
  subroutine fft_WFCD_work_alloc
    integer :: nnl,nnm,nnn
! ----> work arrays of fft for wave functions
    nnl  = fft_box_size_WF(1,1)
    nnm  = fft_box_size_WF(2,1);  nnn  = fft_box_size_WF(3,1)

    if(ipri >= 1) write(nfout,'(" !!allocation of cw1,cw2,cw3 (fft_WFCD_work_alloc)")')

    allocate(cw1(nnl), stat=istat)
    allocate(cw2(nnm), stat=istat)
    allocate(cw3(nnn), stat=istat)
    if(ipri >= 1) write(nfout,'(" !! WF_JRCATFFTg <<m_FFT.fft_WFCD_work_alloc>>")')

! ------- Positron start 
! ----> and work arrays of fft for positron wave functions
    if(sw_positron /= OFF) then
       nnl  = fft_box_size_pWF(1,1)
       nnm  = fft_box_size_pWF(2,1);  nnn  = fft_box_size_pWF(3,1)

       if(ipri >= 1) write(nfout,'(" !!allocation of cw1_pstrn,cw2_pstrn,cw3_pstrn (fft_WFCD_work_alloc)")')
       allocate(cw1_pstrn(nnl), stat=istat)
       allocate(cw2_pstrn(nnm), stat=istat)
       allocate(cw3_pstrn(nnn), stat=istat)
    end if
! ------- Positron end

! ----> work arrays of ifft for the charge density
    nnl = fft_box_size_CD(1,1)
    nnm = fft_box_size_CD(2,1);  nnn = fft_box_size_CD(3,1)
    allocate(wlp(nnl), stat=istat)
    allocate(wmp(nnm), stat=istat)
    allocate(wnp(nnn), stat=istat)
    if(ipri >= 1) write(nfout,'(" !! wlp,wmp,wnp are allocated, CD_JRCATFFT_WS  <<m_FFT.fft_WFCD_work_alloc>>")')
  end subroutine fft_WFCD_work_alloc

! =========== KT_add ======================= 13.0F
  subroutine fft_exx_work_alloc
    integer :: nnl,nnm,nnn

! ----> work arrays of fft for wave functions
    nnl  = fft_box_size_exx(1,1)
    nnm  = fft_box_size_exx(2,1);  nnn  = fft_box_size_exx(3,1)

    if (ipri >= 1) then
       write(nfout,'(" !!allocation of cw1_exx,cw2_exx,cw3_exx (fft_exx_work_alloc)")')
    endif
    allocate(cw1_exx(nnl), stat=istat)
    allocate(cw2_exx(nnm), stat=istat)
    allocate(cw3_exx(nnn), stat=istat)
  end subroutine fft_exx_work_alloc
! ========================================== 13.0F

  subroutine m_FFT_alloc_WF_work
    integer :: nid,nnm,nnn
    integer :: nmax, nsize_ftw

    nid  = fft_box_size_WF(1,0)
    nnm  = fft_box_size_WF(2,1);  nnn  = fft_box_size_WF(3,1)
    if(kimg == 1) then
       nmax = max(nid,nnm,nnn)     ! T. Yamasaki, 28 Jun. 2004
       nsize_ftw = nmax*2*2
    else if(kimg == 2) then
       nsize_ftw = (nid+nnm+nnn)*kimg*2
    end if

    allocate(ftw(nsize_ftw), stat=istat)
  end subroutine m_FFT_alloc_WF_work

  subroutine m_FFT_alloc_pWF_work()
    integer :: nid,nnm,nnn, nmax, nsize_ftw
    nid  = fft_box_size_pWF(1,0)
    nnm  = fft_box_size_pWF(2,1);  nnn  = fft_box_size_pWF(3,1)

    if(kimg == 1) then
       nmax = max(nid,nnm,nnn)     ! T. Yamasaki, 28 Jun. 2004
       nsize_ftw = nmax*2*2
    else if(kimg == 2) then
       nsize_ftw = (nid+nnm+nnn)*kimg*2
    end if

    allocate(ftw(nsize_ftw), stat=istat)
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

! ======= KT_add ======================== 13.0F
  subroutine m_FFT_alloc_exx_work()
    integer :: nid,nnm,nnn, nmax, nsize_ftw

    nid  = fft_box_size_exx(1,0)
    nnm  = fft_box_size_exx(2,1);  nnn  = fft_box_size_exx(3,1)

    if(kimg == 1) then
       nmax = max(nid,nnm,nnn)
       nsize_ftw = nmax*2*2
    else if(kimg == 2) then
       nsize_ftw = (nid+nnm+nnn)*kimg*2
    end if
    allocate(ftw(nsize_ftw), stat=istat)

  end subroutine m_FFT_alloc_exx_work

  subroutine m_FFT_dealloc_exx_work
    if(allocated(ftw)) then
       deallocate(ftw, stat=istat)
       if(istat /= 0 ) then
          if(ipri>=1) then
             write(nfout,*) 'Deallocation error for ftw in sub. m_FFT_dealloc_exx_work'
             write(nfout,*) 'stat=', istat
          end if
          stop
       end if
    end if
  end subroutine m_FFT_dealloc_exx_work
! ======================================== 13.0F

  subroutine m_FFT_alloc_CD_box
    integer :: nid,nnl,nnm,nnn, nmax, nsize_ftw
    nid  = fft_box_size_CD(1,0);  nnl  = fft_box_size_CD(1,1)
    nnm  = fft_box_size_CD(2,1);  nnn  = fft_box_size_CD(3,1)
    if(kimg == 1) then
       nmax = max(nid,nnm,nnn)     ! T. Yamasaki, 28 Jun. 2004
       nsize_ftw = nmax*2*2
    else if(kimg == 2) then
       nsize_ftw = (nid+nnm+nnn)*kimg*2
    end if

    allocate(ftw(nsize_ftw), stat=istat)

    allocate(afft_CD(nfftp), stat=istat)

  end subroutine m_FFT_alloc_CD_box

  subroutine m_FFT_alloc_CDx_box
    integer :: nid,nnl,nnm,nnn, nmax, nsize_ftw
    nid  = fft_box_size_CD_exx(1,0);  nnl  = fft_box_size_CD_exx(1,1)
    nnm  = fft_box_size_CD_exx(2,1);  nnn  = fft_box_size_CD_exx(3,1)
    allocate(wlp_exx(nnl), stat=istat)
    allocate(wmp_exx(nnm), stat=istat)
    allocate(wnp_exx(nnn), stat=istat)

    if(kimg == 1) then
       nmax = max(nid,nnm,nnn)     ! T. Yamasaki, 28 Jun. 2004
       nsize_ftw = nmax*2*2
    else if(kimg == 2) then
       nsize_ftw = (nid+nnm+nnn)*kimg*2
    end if

    allocate(ftw(nsize_ftw), stat=istat)

    allocate(afft_CD_exx(nfftp_exx_nonpara), stat=istat)

  end subroutine m_FFT_alloc_CDx_box

  subroutine m_FFT_dealloc_CD_box
    if(allocated(ftw)) then
       deallocate(ftw, stat=istat)
    end if

    deallocate(afft_CD, stat=istat)
  end subroutine m_FFT_dealloc_CD_box

  subroutine m_FFT_setup(inversion_symmetry,paramset)
    integer, intent(in) :: inversion_symmetry
    logical, intent(in) :: paramset

    real(kind=DP),allocatable,dimension(:)    :: cfft

    integer :: nfft_t, ipad
    integer :: id_sname = -1
    call tstatc0_begin('m_FFT_setup ',id_sname)

! --- fft_box_size_WF, fft_box_size_pWF ---
    if(inversion_symmetry == ON) then  ! kimg == 1
       ipad = 2
    else if(inversion_symmetry == OFF) then ! kimg == 2
       ipad = 1
    end if

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

    idp = fft_box_size_CD(1,0)
    nlp = fft_box_size_CD(1,1)
    nmp = fft_box_size_CD(2,1)
    nnp = fft_box_size_CD(3,1)

    if(sw_positron /= OFF) &
         & nfft_pstrn = product(fft_box_size_pWF(1:3,0))*(2-inversion_symmetry)

    if(.not. paramset) then
       call fft_WFCD_work_alloc     ! <cw[123],cw[123]_pstrn,wlp,wmp,wnp> are allocated

! Initialization of the Wave-Function FFT
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
! Initialization of the Charge-Density FFT
       call CDFFT_setup()
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
         call jrcat_c3ft(cfft,ftw,id,nl,nm,nn,cw1,cw2,cw3,0,0,0,0)
         if(sw_positron /= OFF) call jrcat_c3ft(cfft,ftw,id_p,nl_p,nm_p,nn_p &
              &                                 ,cw1_pstrn,cw2_pstrn,cw3_pstrn,0,0,0,0)
      endif
    end subroutine init_fft_coefficients_arrays_WF

  end subroutine m_FFT_setup

! ========= KT_add ====================================== 13.0F
  subroutine m_FFT_setup_exx(inversion_symmetry,paramset)
    integer, intent(in) :: inversion_symmetry
    logical, intent(in) :: paramset

    real(kind=DP),allocatable,dimension(:)    :: cfft

    integer :: nfft_t, ipad
    integer :: id_sname = -1
    call tstatc0_begin('m_FFT_setup_exx ',id_sname)

! --- fft_box_size_WF, fft_box_size_pWF ---
    if(inversion_symmetry == ON) then  ! kimg == 1
       ipad = 2
    else if(inversion_symmetry == OFF) then ! kimg == 2
       ipad = 1
    end if

    fft_box_size_exx(1,0)   = fft_box_size_exx(1,1) + ipad
    fft_box_size_exx(2:3,0) = fft_box_size_exx(2:3,1)

    nfft_exx = product(fft_box_size_exx(1:3,0)) * (2-inversion_symmetry)

    if(ipri >= 1) call wd_FFTboxsizes_exx(nfout)

    if(.not. paramset) then
       call fft_exx_work_alloc           ! <cw[123]> are allocated

! Initialization of the Wave-Function FFT
       nfft_t = nfft_exx

      allocate(cfft(nfft_t), stat=istat) ! cfft is used only for initiallization
       if(istat /= 0) then
          if(ipri >= 1) then
             write(nfout,*) 'Allocation error for cfft in sub. m_FFT_setup_exx'
             write(nfout,*) 'stat =', istat, 'nfft_exx =', nfft_exx
          end if
          stop
       end if

       call init_fft_coefficients_arrays_WF()

       call m_FFT_dealloc_WF_work()
       deallocate(cfft,stat=istat)
       if(istat /= 0 ) then
          if(ipri >= 1) then
             write(nfout,*) 'Deallocation error for cfft in sub. m_FFT_setup_exx'
             write(nfout,*) 'stat =', istat
          end if
          stop
       end if
    endif

    call tstatc0_end(id_sname)
  contains
    subroutine init_fft_coefficients_arrays_WF()
      integer :: id_exx, nl_exx, nm_exx, nn_exx, ierr

      id_exx = fft_box_size_exx(1,0)
      nl_exx = fft_box_size_exx(1,1)
      nm_exx = fft_box_size_exx(2,1)
      nn_exx = fft_box_size_exx(3,1)

      if (kimg == 1) then
         call jrcat_r3ft( cfft, ftw, id_exx, nl_exx, nm_exx, nn_exx &
              &          ,cw1_exx, cw2_exx, cw3_exx, 0,0, 0,0 )
      else  ! kimg == 2
         call jrcat_c3ft( cfft, ftw,id_exx, nl_exx, nm_exx, nn_exx &
              &           ,cw1_exx, cw2_exx, cw3_exx, 0,0, 0,0 )
      endif
    end subroutine init_fft_coefficients_arrays_WF

  end subroutine m_FFT_setup_exx
! ====================================================== 13.0F

  subroutine CDFFT_setup()
    call m_FFT_alloc_CD_box()       ! <ftw> is allocated
    call init_fft_coefficients_arrays_CD()

    CD_setup_is_done = YES
    call m_FFT_dealloc_CD_box()
  contains
    subroutine init_fft_coefficients_arrays_CD
      if(kimg == 1) then
         call jrcat_r3ft(afft_CD,ftw,idp,nlp,nmp,nnp,wlp,wmp,wnp,0,0,0,0)
      else if(kimg == 2) then
         call jrcat_c3ft(afft_CD,ftw,idp,nlp,nmp,nnp,wlp,wmp,wnp,0,0,0,0)
      endif
    end subroutine init_fft_coefficients_arrays_CD
  end subroutine CDFFT_setup

  subroutine CDFFT_setup_exx()
    idp_exx = fft_box_size_CD_exx(1,0)
    nlp_exx = fft_box_size_CD_exx(1,1)
    nmp_exx = fft_box_size_CD_exx(2,1)
    nnp_exx = fft_box_size_CD_exx(3,1)
    call m_FFT_alloc_CDx_box()       ! <ftw> is allocated
    call init_fft_coefficients_arrays_CDx()

    CD_setup_is_done_exx = YES
  contains
    subroutine init_fft_coefficients_arrays_CDx
      if(kimg == 1) then
         call jrcat_r3ft(afft_CD_exx,ftw,idp_exx,nlp_exx,nmp_exx,nnp_exx,wlp_exx,wmp_exx,wnp_exx,0,0,0,0)
      else if(kimg == 2) then
         call jrcat_c3ft(afft_CD_exx,ftw,idp_exx,nlp_exx,nmp_exx,nnp_exx,wlp_exx,wmp_exx,wnp_exx,0,0,0,0)
      endif
    end subroutine init_fft_coefficients_arrays_CDx
  end subroutine CDFFT_setup_exx

  subroutine m_FFT_WF(electron_or_positron,nfout,afft,inverse_or_direct,switch)  ! G space --> R space
    integer, intent(in)          :: electron_or_positron
    integer, intent(in)          :: nfout
    real(kind=DP), intent(inout) :: afft(nfft)
    integer, intent(in)          :: inverse_or_direct
    integer, intent(in)          :: switch

    complex(kind=CMPLDP),pointer,dimension(:) :: cw1_t, cw2_t, cw3_t

    integer :: id, nl, nm, nn

    integer, dimension(2) :: flag_jrcat_r3ft = (/-2,  2/)
    integer, dimension(2) :: flag_jrcat_c3ft = (/ 1, -1/)
    integer kc1, kc2, kc3

    integer :: id_sname = -1

    if(electron_or_positron == ELECTRON) then
       call tstatc0_begin('m_FFT_WF ',id_sname)

       id = fft_box_size_WF(1,0)
       nl = fft_box_size_WF(1,1)
       nm = fft_box_size_WF(2,1)
       nn = fft_box_size_WF(3,1)

       cw1_t=>cw1; cw2_t=>cw2; cw3_t=>cw3 

    else if(electron_or_positron == POSITRON) then
       call tstatc0_begin('m_FFT_pWF ',id_sname)

       id = fft_box_size_pWF(1,0)
       nl = fft_box_size_pWF(1,1)
       nm = fft_box_size_pWF(2,1)
       nn = fft_box_size_pWF(3,1)

       cw1_t=>cw1_pstrn; cw2_t=>cw2_pstrn; cw3_t=>cw3_pstrn

    end if

    if(switch == ON) then
       kc1 = 0;       kc2 = nm/4 + 1;       kc3 = nn/4 + 1
    else
       kc1 = 0;       kc2 = 0       ;       kc3 = 0
    endif
    if(kimg == 1) then
       call jrcat_r3ft(afft,ftw,id,nl,nm,nn,cw1_t,cw2_t,cw3_t,kc1,kc2,kc3&
            & ,flag_jrcat_r3ft(inverse_or_direct))
    else if(kimg == 2) then
       call jrcat_c3ft(afft,ftw,id,nl,nm,nn,cw1_t,cw2_t,cw3_t,kc1,kc2,kc3&
            & ,flag_jrcat_c3ft(inverse_or_direct))
    endif

    call tstatc0_end(id_sname)
  end subroutine m_FFT_WF

! =========== KT_add ======================================= 13.0F
  subroutine m_FFT_exx(nfout,afft,inverse_or_direct,switch)  ! G space --> R space
    integer, intent(in)          :: nfout
    real(kind=DP), intent(inout) :: afft(nfft_exx)
    integer, intent(in)          :: inverse_or_direct
    integer, intent(in)          :: switch

    complex(kind=CMPLDP),pointer,dimension(:) :: cw1_t, cw2_t, cw3_t

    integer :: id, nl, nm, nn

    integer, dimension(2) :: flag_jrcat_r3ft = (/-2,  2/)
    integer, dimension(2) :: flag_jrcat_c3ft = (/ 1, -1/)
    integer kc1, kc2, kc3

    integer :: id_sname = -1

    call tstatc0_begin('m_FFT_exx ',id_sname)

    id = fft_box_size_exx(1,0)
    nl = fft_box_size_exx(1,1)
    nm = fft_box_size_exx(2,1)
    nn = fft_box_size_exx(3,1)

    cw1_t=>cw1_exx; cw2_t=>cw2_exx; cw3_t=>cw3_exx

    if(switch == ON) then
       kc1 = 0;       kc2 = nm/4 + 1;       kc3 = nn/4 + 1
    else
       kc1 = 0;       kc2 = 0       ;       kc3 = 0
    endif
    if(kimg == 1) then
       call jrcat_r3ft(afft,ftw,id,nl,nm,nn,cw1_t,cw2_t,cw3_t,kc1,kc2,kc3&
            & ,flag_jrcat_r3ft(inverse_or_direct))
    else if(kimg == 2) then
       call jrcat_c3ft(afft,ftw,id,nl,nm,nn,cw1_t,cw2_t,cw3_t,kc1,kc2,kc3&
            & ,flag_jrcat_c3ft(inverse_or_direct))
    endif

    call tstatc0_end(id_sname)
  end subroutine m_FFT_exx
! =========================================================== 13.0F

  subroutine fft_CD_inverse_core(afft_CD)
    real(kind=DP), intent(inout) :: afft_cd(nfftp)
    integer :: nsize_ftw
    integer :: kc1 = 0, kc2 = 0, kc3 = 0

    nsize_ftw = 0
    if(.not.allocated(ftw)) then
       if(kimg==1) then
          nsize_ftw = max(idp,nmp,nnp)*2*2
       else
          nsize_ftw = (idp+nmp+nnp)*kimg*2
       end if
       allocate(ftw(nsize_ftw))
    end if

    if(kimg == 1) then
       call jrcat_r3ft(afft_cd,ftw,idp,nlp,nmp,nnp,wlp,wmp,wnp,kc1,kc2,kc3&
            & ,-2)
    else if(kimg == 2) then
       call jrcat_c3ft(afft_cd,ftw,idp,nlp,nmp,nnp,wlp,wmp,wnp,kc1,kc2,kc3&
            & , 1)
    endif

    if(nsize_ftw /= 0) deallocate(ftw)
  end subroutine fft_CD_inverse_core

  subroutine fft_CD_direct_core(afft_CD)
    real(kind=DP), intent(inout) :: afft_cd(nfftp)
    integer :: nsize_ftw
    integer :: kc1 = 0, kc2 = 0, kc3 = 0

    nsize_ftw = 0
    if(.not.allocated(ftw)) then
       if(kimg==1) then
          nsize_ftw = max(idp,nmp,nnp)*2*2
       else
          nsize_ftw = (idp+nmp+nnp)*kimg*2
       end if
       allocate(ftw(nsize_ftw))
    end if

    if(kimg == 1) then
       call jrcat_r3ft(afft_CD,ftw,idp,nlp,nmp,nnp,wlp,wmp,wnp,kc1,kc2,kc3&
            & , 2)
    else if(kimg == 2) then
       call jrcat_c3ft(afft_CD,ftw,idp,nlp,nmp,nnp,wlp,wmp,wnp,kc1,kc2,kc3&
            & ,-1)
    endif

    if(nsize_ftw /= 0) deallocate(ftw)
  end subroutine fft_CD_direct_core

! --------------
  subroutine fft_CD_inverse_core_exx(afft_CD)
    real(kind=DP), intent(inout) :: afft_cd(nfftp_exx_nonpara)
    integer :: nsize_ftw
    integer :: kc1 = 0, kc2 = 0, kc3 = 0

    nsize_ftw = 0
    if(.not.allocated(ftw)) then
       if(kimg==1) then
          nsize_ftw = max(idp_exx,nmp_exx,nnp_exx)*2*2
       else
          nsize_ftw = (idp_exx+nmp_exx+nnp_exx)*kimg*2
       end if
       allocate(ftw(nsize_ftw))
    end if

    if(kimg == 1) then
       call jrcat_r3ft(afft_cd,ftw,idp_exx,nlp_exx,nmp_exx,nnp_exx,wlp_exx,wmp_exx,wnp_exx,kc1,kc2,kc3&
            & ,-2)
    else if(kimg == 2) then
       call jrcat_c3ft(afft_cd,ftw,idp_exx,nlp_exx,nmp_exx,nnp_exx,wlp_exx,wmp_exx,wnp_exx,kc1,kc2,kc3&
            & , 1)
    endif

    if(nsize_ftw /= 0) deallocate(ftw)
  end subroutine fft_CD_inverse_core_exx

  subroutine fft_CD_direct_core_exx(afft_CD)
    real(kind=DP), intent(inout) :: afft_cd(nfftp_exx_nonpara)
    integer :: nsize_ftw
    integer :: kc1 = 0, kc2 = 0, kc3 = 0
    nsize_ftw = 0
    if(.not.allocated(ftw)) then
       if(kimg==1) then
          nsize_ftw = max(idp_exx,nmp_exx,nnp_exx)*2*2
       else
          nsize_ftw = (idp_exx+nmp_exx+nnp_exx)*kimg*2
       end if
       allocate(ftw(nsize_ftw))
    end if

    if(kimg == 1) then
       call jrcat_r3ft(afft_CD,ftw,idp_exx,nlp_exx,nmp_exx,nnp_exx,wlp_exx,wmp_exx,wnp_exx,kc1,kc2,kc3&
            & , 2)
    else if(kimg == 2) then
       call jrcat_c3ft(afft_CD,ftw,idp_exx,nlp_exx,nmp_exx,nnp_exx,wlp_exx,wmp_exx,wnp_exx,kc1,kc2,kc3&
            & ,-1)
    endif

    if(nsize_ftw /= 0) deallocate(ftw)
  end subroutine fft_CD_direct_core_exx

