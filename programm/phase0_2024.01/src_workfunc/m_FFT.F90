!================================================
!  Software name : STM
!  Module : m_FFT
!  Subroutine(s) : m_FFT_Vlocal_W_fine, m_FFT_query_inversion_symmetry
!                  m_FFT_rd_fine_fft_box_size, m_FFT_rd_fft_box_size,
!                  m_FFT_set_fft_box_size_fine, Vlocal_W, fft_phase_work_alloc
!                  fft_WF_work_alloc, fft_WF_work_dealloc, m_FFT_setup,
!                  init_fft_coefficients_arrays, WD_FFTboxsizes, FFT_WF, FFT_WF_fine,
!                  FFT_fine, FFT_CD, check_of_negative_CD, W_Vlocal_W, m_FFT_dFFT_fine,
!                  WF_2dFFT_backward, WF_dFFT_forward, substitute_afft_with_data, pr_afft
!  Author(s)     : Takahiro Yamasaki and Koichi Kato (June 7, 2004)
!
!  Contact address :  IIS,The University of Tokyo RSS21 project
!  
!  "Multiscale Simulation System for Function Analysis of Nanomaterials"  
!
!================================================
!
!     The original version of this program "STM" was developed in the
!  period 1998-2001 at JRCAT by Koichi Kato (Toshiba Corporate
!  Research and Development Center) and Takahiro Yamasaki (Fujitsu 
!  Laboratories Ltd.) through a joint research between JRCAT and
!  Toshiba Corporate Research and Development Center and another joint
!  research between JRCAT and FUJITSU Laboratories Ltd.
!     Since 2003, this program has been tuned and revised as it
!  matches to the first-principles molecular dynamics program "PHASE"
!  as a part of the national project "Frontier Simulation Software for
!  Industrial Science (FSIS)", which is supported by the IT program of
!  the Ministry of Education, Culture, Sports, Science and Technology
!  (MEXT) of Japan.
!     Since 2006, this program set has been developed as a part of the
!  national project "Revolutionary Simulation Software (RSS21)", which
!  is supported by the next-generation IT program of MEXT of Japan.

module m_FFT
! $Id: m_FFT.F90,v 1.1.1.1 2004/06/26 11:21:37 yamasaki Exp $
  use m_ArraySize_Parameters, only : kimg
  use m_Const_Parameters, only : DP,CMPLDP, ON, OFF, INVERSE, DIRECT, PAI2
  use m_Control_Parameters,only: n_fc

  implicit none

  integer, dimension(0:3) :: fft_box_size_WF, fft_box_size_CD &
       &,   fft_box_size_fine
  integer                 :: nfft, nfftp, nfftpf
  integer                 :: nlpf

  complex(kind=CMPLDP), allocatable, dimension(:)    :: cw1,cw2,cw3
  complex(kind=CMPLDP), allocatable, dimension(:)    :: cw1p,cw2p,cw3p
  complex(kind=CMPLDP), allocatable, dimension(:)    :: cw1pf,cw2pf,cw3pf

  real(kind=DP), allocatable, dimension(:)           :: cfft

  real(kind=DP), allocatable, dimension(:)           :: ftw

contains
  subroutine m_FFT_Vlocal_W_fine(afft,bfft,cfft)
    real(kind=DP), intent(in),  dimension(nfftpf) :: afft,bfft
    real(kind=DP), intent(out), dimension(nfftpf) :: cfft

    integer :: i
    !  V_{local} := bfft   : real
    !  \Psi      := afft   : complex

!!$    print *, ' nfftpf = ', nfftpf
    do i = 1, nfftpf-1, 2
       cfft(i)   = bfft(i) * afft(i)
       cfft(i+1) = bfft(i) * afft(i+1)
    end do
  end subroutine M_FFT_Vlocal_W_fine

  subroutine m_FFT_query_inversion_symmetry(inversion_symmetry)
    integer,intent(out) :: inversion_symmetry
    if(kimg == 1) inversion_symmetry = ON
    if(kimg == 2) inversion_symmetry = OFF
  end subroutine m_FFT_query_inversion_symmetry

  subroutine m_FFT_rd_fine_fft_box_size(nfinp)
    integer,intent(in) :: nfinp
    read(nfinp,*) nlpf
  end subroutine m_FFT_rd_fine_fft_box_size

  subroutine m_FFT_rd_fft_box_size(nfcntn_bin)
    integer,intent(in) :: nfcntn_bin
    read(nfcntn_bin) fft_box_size_WF(1:3),fft_box_size_CD(1:3)
#ifdef DEBUG
    write(0,*) 'fft box size wf',fft_box_size_WF(1),fft_box_size_WF(2),fft_box_size_WF(3)
    write(0,*) 'fft box size cd',fft_box_size_CD(1),fft_box_size_CD(2),fft_box_size_CD(3)
#endif
  end subroutine m_FFT_rd_fft_box_size

  subroutine m_FFT_set_fft_box_size_fine
    if( kimg == 1 ) then
        nlpf = fft_box_size_CD(n_fc(1))* 0.5
    else if ( kimg == 2) then
        nlpf = fft_box_size_CD(n_fc(1))
    end if
    fft_box_size_fine(1) = fft_box_size_CD(n_fc(1))
    fft_box_size_fine(2) = fft_box_size_CD(n_fc(2))
    fft_box_size_fine(3) = fft_box_size_CD(n_fc(3))
  end subroutine m_FFT_set_fft_box_size_fine

  subroutine Vlocal_W(afft,bfft)
    real(kind=DP), intent(in),    dimension(nfft) :: afft
    real(kind=DP), intent(inout), dimension(nfft) :: bfft

    integer i
    do i = 1, nfft-1, 2
       bfft(i)   = afft(i) * bfft(i)
       bfft(i+1) = afft(i) * bfft(i+1)
    end do
    if(mod(nfft,2) == 1) bfft(nfft) = afft(nfft)*bfft(nfft)
  end subroutine Vlocal_W

  subroutine fft_phase_work_alloc
    integer :: nnl,nnm,nnn

    nnl  = fft_box_size_WF(1);  nnm  = fft_box_size_WF(2)
    nnn  = fft_box_size_WF(3)
    allocate(cw1(nnl)); allocate(cw2(nnm)); allocate(cw3(nnn))

    nnl = fft_box_size_CD(1);   nnm = fft_box_size_CD(2)
    nnn = fft_box_size_CD(3)
    allocate(cw1p(nnl)); allocate(cw2p(nnm)); allocate(cw3p(nnn))

    nnl = fft_box_size_fine(1); nnm = fft_box_size_fine(2)
    nnn = fft_box_size_fine(3)
    allocate(cw1pf(nnl));allocate(cw2pf(nnm));allocate(cw3pf(nnn))
  end subroutine fft_phase_work_alloc

  subroutine fft_WF_work_alloc
    integer nfftwk
    nfftwk = nfft
    if(nfftwk < nfftp)  nfftwk = nfftp
    if(nfftwk < nfftpf) nfftwk = nfftpf
    allocate(ftw(nfftwk))
  end subroutine fft_WF_work_alloc

  subroutine fft_WF_work_dealloc
    deallocate(ftw)
  end subroutine fft_WF_work_dealloc

  subroutine m_FFT_setup(nfout,inversion_symmetry)
    integer, intent(in)           :: nfout, inversion_symmetry

    integer :: id_sname = -1

    if(inversion_symmetry == ON) then
       fft_box_size_WF(0) = fft_box_size_WF(1) + 2
       fft_box_size_CD(0) = fft_box_size_CD(1) + 2
       fft_box_size_fine(0) = fft_box_size_fine(1) + 2
    else if(inversion_symmetry == OFF) then
       fft_box_size_WF(0) = fft_box_size_WF(1) + 1
       fft_box_size_CD(0) = fft_box_size_CD(1) + 1
       fft_box_size_fine(0) = fft_box_size_fine(1) + 1
    endif
    nfft = fft_box_size_WF(0)*fft_box_size_WF(2)*fft_box_size_WF(3)&
         &*(2-inversion_symmetry)
    nfftp= fft_box_size_CD(0)*fft_box_size_CD(2)*fft_box_size_CD(3)&
         &*(2-inversion_symmetry)
    nfftpf= fft_box_size_fine(0)*fft_box_size_fine(2)*fft_box_size_fine(3)&
         &*(2-inversion_symmetry)

    call fft_WF_work_alloc
    allocate(cfft(nfft))
    call fft_phase_work_alloc

    call init_fft_coefficients_arrays  ! -(contained here)

    deallocate(cfft)

  contains
    subroutine init_fft_coefficients_arrays
      integer id, nl, nm, nn
      id = fft_box_size_WF(0)
      nl = fft_box_size_WF(1)
      nm = fft_box_size_WF(2)
      nn = fft_box_size_WF(3)

      if(kimg == 1) then
         call jrcat_r3ft(cfft,ftw,id,nl,nm,nn,cw1,cw2,cw3,0,0,0,0)
      else
         call jrcat_c3ft(cfft,ftw,id,nl,nm,nn,cw1,cw2,cw3,0,0,0,0)
      endif

      id = fft_box_size_CD(0)
      nl = fft_box_size_CD(1)
      nm = fft_box_size_CD(2)
      nn = fft_box_size_CD(3)

      if(kimg == 1) then
         call jrcat_r3ft(cfft,ftw,id,nl,nm,nn,cw1p,cw2p,cw3p,0,0,0,0)
      else
         call jrcat_c3ft(cfft,ftw,id,nl,nm,nn,cw1p,cw2p,cw3p,0,0,0,0)
      endif

      id = fft_box_size_fine(0)
      nl = fft_box_size_fine(1)
      nm = fft_box_size_fine(2)
      nn = fft_box_size_fine(3)

      if(kimg == 1) then
         call jrcat_r3ft(cfft,ftw,id,nl,nm,nn,cw1pf,cw2pf,cw3pf,0,0,0,0)
      else
         call jrcat_c3ft(cfft,ftw,id,nl,nm,nn,cw1pf,cw2pf,cw3pf,0,0,0,0)
      endif

    end subroutine init_fft_coefficients_arrays

    subroutine WD_FFTboxsizes(nfout)
      integer, intent(in) :: nfout

      write(nfout,*) ' ---- (WF)FFT size ----'
      write(nfout,'(a44,i7)') &
           &     '  total FFT elements(including work area) = ', nfft
      write(nfout,'(a44,3i5)') &
           & ' actual FFT box size                      = ' &
           & , fft_box_size_WF(0), fft_box_size_WF(2), fft_box_size_WF(3)
      write(nfout,'(a44,3i5)') &
           & '   real FFT size                          = ' &
           & , fft_box_size_WF(1), fft_box_size_WF(2), fft_box_size_WF(3)

      write(nfout,*) ' ---- (CD)FFT size ----'
      write(nfout,'(a44,i7)') &
           & '  total FFT elements(including work area) = ', nfftp
      write(nfout,'(a44,3i5)') &
           & ' actual FFT box size                      = ' &
           & , fft_box_size_CD(0), fft_box_size_CD(2), fft_box_size_CD(3)
      write(nfout,'(a44,3i5)') &
           & '   real FFT size                          = ' &
           & , fft_box_size_CD(1), fft_box_size_CD(2), fft_box_size_CD(3)

      write(nfout,*) ' ---- (fine)FFT size ----'
      write(nfout,'(a44,i7)') &
           & '  total FFT elements(including work area) = ', nfftpf
      write(nfout,'(a44,3i5)') &
           & ' actual FFT box size                      = ' &
           & ,fft_box_size_fine(0),fft_box_size_fine(2),fft_box_size_fine(3)
      write(nfout,'(a44,3i5)') &
           & '   real FFT size                          = ' &
           & ,fft_box_size_fine(1),fft_box_size_fine(2),fft_box_size_fine(3)
    end subroutine WD_FFTboxsizes

  end subroutine m_FFT_setup

  subroutine FFT_WF(afft,inverse_or_direct,switch)  ! G space --> R space
    real(kind=DP), intent(inout) :: afft(nfft)
    integer, intent(in)          :: inverse_or_direct
    integer, intent(in)          :: switch
    integer :: id, nl, nm, nn, ier1 = 0

    integer, dimension(2) :: flag_jrcat_r3ft = (/-2,  2/)
    integer, dimension(2) :: flag_jrcat_c3ft = (/ 1, -1/)

    integer kc1, kc2, kc3

    integer :: id_sname = -1

    id = fft_box_size_WF(0)
    nl = fft_box_size_WF(1)
    nm = fft_box_size_WF(2)
    nn = fft_box_size_WF(3)

    if(switch == ON) then
       kc1 = 0;       kc2 = nm/4 + 1;       kc3 = nn/4 + 1
    else
       kc1 = 0;       kc2 = 0       ;       kc3 = 0
    endif
    if(kimg == 1) then
       call jrcat_r3ft(afft,ftw,id,nl,nm,nn,cw1,cw2,cw3,kc1,kc2,kc3&
            & ,flag_jrcat_r3ft(inverse_or_direct))
    else if(kimg == 2) then
       call jrcat_c3ft(afft,ftw,id,nl,nm,nn,cw1,cw2,cw3,kc1,kc2,kc3&
            & ,flag_jrcat_c3ft(inverse_or_direct))
    endif
    if (ier1 /= 0) then
       print *, ' !!(Error fft) ier1 = ',ier1
       stop
    end if
  end subroutine FFT_WF

  subroutine FFT_WF_fine(afft,inverse_or_direct)
    real(kind=DP), intent(inout) :: afft(nfftpf)
    integer, intent(in)          :: inverse_or_direct
!!$    call FFT_CD(afft,inverse_or_direct)
    call FFT_fine(afft,inverse_or_direct)
  end subroutine FFT_WF_fine

  subroutine FFT_fine(afft,inverse_or_direct)   ! R space --> G space
    real(kind=DP), intent(inout) :: afft(nfftpf)
    integer, intent(in)          :: inverse_or_direct

    integer :: idp,nlp, nmp, nnp, ier1 = 0, kc1 = 0, kc2 = 0, kc3 = 0
    integer :: id_sname = -1
    integer, dimension(2) :: flag_jrcat_r3ft = (/-2,  2/)
    integer, dimension(2) :: flag_jrcat_c3ft = (/ 1, -1/)

 
    idp = fft_box_size_fine(0)
    nlp = fft_box_size_fine(1)
    nmp = fft_box_size_fine(2)
    nnp = fft_box_size_fine(3)

    if(kimg == 1) then
       call jrcat_r3ft(afft,ftw,idp,nlp,nmp,nnp,cw1pf,cw2pf,cw3pf &
            & ,kc1,kc2,kc3, flag_jrcat_r3ft(inverse_or_direct))
    else if(kimg == 2) then
       call jrcat_c3ft(afft,ftw,idp,nlp,nmp,nnp,cw1pf,cw2pf,cw3pf &
            & ,kc1,kc2,kc3 ,flag_jrcat_c3ft(inverse_or_direct))
    endif

    if (ier1 /= 0) then
       print *, ' !!Error  fft(direct) ier1 = ',ier1
       stop
    end if
  end subroutine FFT_fine

  subroutine FFT_CD(afft,inverse_or_direct)   ! R space --> G space
    real(kind=DP), intent(inout) :: afft(nfftp)
    integer, intent(in)          :: inverse_or_direct

    integer :: idp,nlp, nmp, nnp, ier1 = 0, kc1 = 0, kc2 = 0, kc3 = 0
    integer :: id_sname = -1
    integer, dimension(2) :: flag_jrcat_r3ft = (/-2,  2/)
    integer, dimension(2) :: flag_jrcat_c3ft = (/ 1, -1/)


    idp = fft_box_size_CD(0)
    nlp = fft_box_size_CD(1)
    nmp = fft_box_size_CD(2)
    nnp = fft_box_size_CD(3)

    if(kimg == 1) then
       call jrcat_r3ft(afft,ftw,idp,nlp,nmp,nnp,cw1p,cw2p,cw3p,kc1,kc2,kc3&
            &, flag_jrcat_r3ft(inverse_or_direct))
    else if(kimg == 2) then
       call jrcat_c3ft(afft,ftw,idp,nlp,nmp,nnp,cw1p,cw2p,cw3p,kc1,kc2,kc3&
            & ,flag_jrcat_c3ft(inverse_or_direct))
    endif

    if (ier1 /= 0) then
       print *, ' !!Error  fft(direct) ier1 = ',ier1
       stop
    end if
  end subroutine FFT_CD

  subroutine check_of_negative_CD(afft,nfout,nspin,ispin)
    real(kind=DP), intent(inout) :: afft(nfftp)
    integer, intent(in)          :: nfout,nspin,ispin

    integer                   :: imodi, i
    integer                   :: icwarn_NEGA, icwarn_IMAG
    integer, parameter        :: MXWARN = 10
    real(kind=DP),parameter   :: chgdel = 8.d-5

    if(kimg.eq.1) then
       imodi = nfftp+10
    else if(kimg.eq.2) then
       imodi = fft_box_size_CD(0)
    end if

    icwarn_NEGA = 0
    icwarn_IMAG = 0
    do i = 1,nfftp-1,2
       if( (afft(i) <= -chgdel) .and. (mod(i+1,imodi) /= 0)) then
          icwarn_NEGA = icwarn_NEGA + 1
          if(icwarn_NEGA.le.MXWARN) then
             if(icwarn_NEGA.eq.1.and.nspin.eq.2) then
                write(nfout,*) ' #spin = ', ispin
             endif
             write(nfout,110) i,afft(i)
110          format(1h ,'**WARN CHG.DEN<0.0 AT ',I7,2D15.7,'***')
          endif
          afft(i) = 1.d-40
       else if (afft(i) <= 0.d0) then
          afft(i) = 1.d-40
       end if
       if( (abs(afft(i+1)) > chgdel).and. (mod(i+1,imodi) /= 0) ) then
          icwarn_IMAG = icwarn_IMAG + 1
          if(icwarn_IMAG.le.MXWARN) then
             if(icwarn_IMAG.eq.1.and.nspin.eq.2) then
                write(nfout,*) ' #spin = ', ispin
             endif
             write(nfout,120) i,afft(i+1)
120          format(1h ,'**WARN IMAG(CHG)>0.0 AT ',I7,2D15.7,'***')
          endif
          afft(i+1) = 0.d0
       else if (abs(afft(i+1)).gt.0.d0) then
          afft(i+1) = 0.d0
       end if
    end do
    if(icwarn_NEGA.gt.MXWARN) then
       if(nspin.eq.2) &
            &   write(nfout,*) ' #spin = ',ispin,' : 1 = UP, 2 = DOWN'
       write(nfout,*) &
            &  ' Warning of <<Negative Charge Density>> = ',icwarn_NEGA
    endif
    if(icwarn_IMAG.gt.MXWARN) then
       if(nspin.eq.2)  &
            &   write(nfout,*) ' #spin = ',ispin,' : 1 = UP, 2 = DOWN'
       write(nfout,*) &
            &  ' Warning of <<Imaginary Charge Density>> = ',icwarn_IMAG
    endif
  end subroutine check_of_negative_CD

  subroutine W_Vlocal_W(afft,bfft,eg)
    real(kind=DP), intent(in),    dimension(nfft) :: afft
    real(kind=DP), intent(inout), dimension(nfft) :: bfft
    real(kind=DP), intent(out)                 :: eg

    integer i, nl, nm, nn, id
    real(kind=DP) :: vcel, s

    id = fft_box_size_WF(0)
    nl = fft_box_size_WF(1)
    nm = fft_box_size_WF(2)
    nn = fft_box_size_WF(3)

    do i = 1, nfft-1, 2
       bfft(i) = afft(i)*(bfft(i)**2 + bfft(i+1)**2)
    end do
    vcel = 1.d0/(nl*nm*nn)
    s = 0.d0
    if(kimg == 1) then
       do i = 1, nfft-1, 2
          s = s + bfft(i)
       end do
       s = s + s
       do i = 1, nfft-1, id
          s = s - bfft(i)
       end do
       do i = id-1, nfft-1, id
          s = s - bfft(i)
       end do
    else if(kimg == 2) then
       do i = 1, nfft-1, 2
          s = s + bfft(i)
       end do
    end if
    eg = s * vcel
  end subroutine W_Vlocal_W

  subroutine m_FFT_WF_2dFFT_fine(afft,inverse_or_direct)
    real(kind=DP), intent(inout) :: afft(nfftpf)
    integer, intent(in)          :: inverse_or_direct
!	INVERSE: 1, DIRECT: 2

    integer, dimension(2) :: flag_jrcat_r2ft = (/-1, 1/)
    integer, dimension(2) :: flag_jrcat_c2ft = (/ 1, -1/)

    real(kind=DP) :: factor
    integer :: idp, nlp, nmp, nnp

    integer :: id_sname = -1


    idp = fft_box_size_fine(0)
    nlp = fft_box_size_fine(1)
    nmp = fft_box_size_fine(2)
    nnp = fft_box_size_fine(3)

    ftw = 0.d0
    if(kimg == 1) then
       call jrcat_r2ft_of_3D(afft,ftw,idp,nlp,nmp,nnp,cw2pf,cw3pf &
            & ,flag_jrcat_r2ft(inverse_or_direct))
    else if(kimg == 2) then
       call jrcat_c2ft_of_3D(afft,ftw,idp,nlp,nmp,nnp,cw2pf,cw3pf &
            & ,flag_jrcat_c2ft(inverse_or_direct))
    endif
    if(inverse_or_direct .eq. DIRECT) then
       factor = 1.d0/(nmp*nnp)
       afft = afft*factor
    end if

  end subroutine M_FFT_WF_2dFFT_fine

  subroutine WF_2dFFT_backward(afft)
    real(kind=DP), intent(inout) :: afft(nfftpf)

    real(kind=DP) :: factor
    integer :: idp, nlp, nmp, nnp
    integer :: id_sname = -1
    integer :: key

    key = +1

    idp = fft_box_size_fine(0)
    nlp = fft_box_size_fine(1)
    nmp = fft_box_size_fine(2)
    nnp = fft_box_size_fine(3)

    factor = 1.d0/(nmp*nnp)

    ftw = 0.d0
    if(kimg == 1) then
       call jrcat_r2ft_of_3D(afft,ftw,idp,nlp,nmp,nnp,cw2pf,cw3pf,key)
!!$    call jrcat_r2ft_of_3D_back(afft,ftw,idp,nlp,nmp,nnp,cw2pf,cw3pf)
    else if(kimg == 2) then
       call jrcat_c2ft_of_3D(afft,ftw,idp,nlp,nmp,nnp,cw2pf,cw3pf,key)
!!$       call jrcat_c2ft_of_3D_back(afft,ftw,idp,nlp,nmp,nnp,cw2pf,cw3pf)
    endif
    afft = afft*factor

  end subroutine WF_2dFFT_backward

  subroutine WF_2dFFT_forward(afft)
    real(kind=DP), intent(inout) :: afft(nfftpf)

    integer :: idp, nlp, nmp, nnp
    integer :: id_sname = -1
    integer :: key

    key = -1

    idp = fft_box_size_fine(0)
    nlp = fft_box_size_fine(1)
    nmp = fft_box_size_fine(2)
    nnp = fft_box_size_fine(3)

!!$    factor = 1.d0/(nlp*nmp*nnp)

    ftw = 0.d0
    if(kimg == 1) then
!!$       call jrcat_r1ft_of_3D_back(afft,ftw,idp,nlp,nmp,nnp,cw1pf)
!!$       call jrcat_r3ft(afft,ftw,idp,nlp,nmp,nnp,cw1pf,cw2pf,cw3pf,0,0,0,-1)
!!$    call jrcat_r2ft_of_3D_forward(afft,ftw,idp,nlp,nmp,nnp,cw2pf,cw3pf)
       call jrcat_r2ft_of_3D(afft,ftw,idp,nlp,nmp,nnp,cw2pf,cw3pf,key)
!!$       afft = afft*factor
    else if(kimg == 2) then
!!$       call jrcat_c1ft_of_3D_back(afft,ftw,idp,nlp,nmp,nnp,cw1p)
!!$       call jrcat_c3ft(afft,ftw,idp,nlp,nmp,nnp,cw1p,cw2p,cw3p,0,0,0,-1)
!!$    call jrcat_c2ft_of_3D_forward(afft,ftw,idp,nlp,nmp,nnp,cw2pf,cw3pf)
       call jrcat_c2ft_of_3D(afft,ftw,idp,nlp,nmp,nnp,cw2pf,cw3pf,key)
!!$       afft = afft*factor
    endif

  end subroutine WF_2dFFT_forward

  subroutine substitute_afft_with_data(p,q,afft)
    integer, intent(in)                           :: p,q
    real(kind=DP), intent(out), dimension(nfftpf) :: afft
    real(kind=DP), pointer, dimension(:)          :: d1,d2,d3

    integer       :: nlpfh, nlpf_t, nmpf, nnpf, i,j,k,ip,igfpf, imax
    real(kind=DP) :: phase

    nlpf_t = fft_box_size_fine(1)
    if(kimg == 1) then
       nlpfh = fft_box_size_fine(0)/2
    else
       nlpfh = fft_box_size_fine(0)
    end if
    nmpf = fft_box_size_fine(2); nnpf = fft_box_size_fine(3)

    if(kimg == 1) then
       imax = nlpfh
    else
       imax = nlpf_t
    end if

    allocate(d1(imax)); allocate(d2(nmpf)); allocate(d3(nnpf))

    do i = 1, imax
       d1(i) = (nlpf_t-i+1)/dble(nlpf_t)
    end do

    do i = 1, nmpf
       phase = PAI2*p*(i-1)/nmpf
       d2(i) = dcos(phase)
    end do

    do i = 1, nnpf
       phase = PAI2*q*(i-1)/nnpf
       d3(i) = dcos(phase)
    end do

    do k = 1, nnpf
       do j = 1, nmpf
          do i = 1, imax
             igfpf = i + (j-1)*nlpfh + (k-1)*nlpfh*nmpf
             ip = 2*igfpf-1
             afft(ip) = d1(i)*d2(j)*d3(k)
             afft(ip+1) = 0.d0
          end do
       end do
    end do

    deallocate(d1); deallocate(d2); deallocate(d3)
  end subroutine substitute_afft_with_data

  subroutine pr_afft(p,q,afft)
    integer, intent(in)                          :: p,q
    real(kind=DP), intent(in), dimension(nfftpf) :: afft

    integer :: i,j,k,imin,imax,jmin,jmax,kmin,kmax,nlpfh,nmpf
    integer, pointer, dimension(:) :: ip
    
    if(kimg == 1) then
       nlpfh = fft_box_size_fine(0)/2
    else
       nlpfh = fft_box_size_fine(0)
    end if
    nmpf  = fft_box_size_fine(2)

    imin = 1; imax = 3

    jmin = p-3; jmax = p+3
    if(jmin < 1) jmin = 1
    if(jmax > fft_box_size_fine(2)) jmax = fft_box_size_fine(2)
    allocate(ip(jmin:jmax))

    kmin = q-3; kmax = q+3
    if(kmin < 1) kmin = 1
    if(kmax > fft_box_size_fine(3)) kmax = fft_box_size_fine(3)


    print '(" p,q = ",i5,i5)',p,q
    do i = imin, imax
       print '(" i = ", i5)',i
       print '(" range : (",i5,":",i5,",",i5,":",i5,")")' &
            & , jmin,jmax,kmin,kmax
       do k = kmin, kmax
          do j = jmin, jmax
             ip(j) = 2*(i+(j-1)*nlpfh+(k-1)*nlpfh*nmpf) - 1
          end do
          print '(" (k = ",i5,")",8f10.4)',k,(afft(ip(j)),j=jmin,jmax)
       end do
    end do
  end subroutine pr_afft

end module m_FFT
