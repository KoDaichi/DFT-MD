!*********************************************************
!* 逆空間でのクーロン力を遮蔽するポテンシャルを作成するモジュールの
!* FFTサブモジュール
!*********************************************************

! PHASE の m_FFT.F90 からコードを拝借。
! kimg=1の空間反転対称の条件でコードを簡素化。
! 密度のFFTルーチンのみでコードを簡素化。

#ifdef CD_JRCATFFT
# define CD_JRCATFFTg
#endif

#ifdef CD_JRCATFFT_WS
# define CD_JRCATFFTg
#endif

#ifndef ASLFFT
# ifndef DECFFT
#  ifndef MKLFFT
#   ifndef FFTW3
#    ifndef CD_JRCATFFTg
#     define CD_IFFT
#    endif
#   endif
#  endif
# endif
#endif

!module m_Const_Parameters
!  integer, parameter :: DP = kind(1.d0)
!  integer, parameter :: CMPLDP = kind((1.d0, 0.d0))
!  integer, parameter :: OUTER = 1
!end module

module m_Screening_FFT
#ifdef MKLFFT
  use MKL_DFTI
#endif

  use m_Const_Parameters, only : DP, CMPLDP, OUTER

  implicit none


  integer, dimension(3,0:1) :: fft_box_size_CD
  integer                   :: nfftp
  integer, private          :: istat = 0

  real(kind=DP),allocatable,dimension(:) :: afft_CD
#ifndef DECFFT
# ifndef MKLFFT
#  ifndef FFTW3

#   ifdef ASLFFT
  real(kind=DP),private,allocatable,dimension(:)           :: trigs_CD
  integer,      private,allocatable,dimension(:)           ::ifax_CD
#   else
  complex(kind=CMPLDP),private,allocatable,dimension(:)    :: wlp,wmp,wnp
  integer,             private,allocatable,dimension(:,:)  :: iwork
  integer, private                                         :: iopt = 0
#   endif
! /* ASLFFT */

#  endif
! /* FFTW3 */
# endif
! /* MKLFFT */
#endif
! /* DECFFT */

#ifdef MKLFFT
  type(DFTI_DESCRIPTOR), POINTER :: dh_CD
  type(DFTI_DESCRIPTOR), POINTER :: dh_CD_back
  integer :: strides(4), Status
  complex(DP), allocatable :: cwork_mkl(:)
#endif

#ifdef FFTW3
  integer(8) :: plan_CD(2)
#endif

#ifdef CD_JRCATFFTg
  real(kind=DP),private,allocatable,dimension(:) :: ftw
#elif GOEDECKER_FFT
  real(kind=DP),private,allocatable,dimension(:) :: ftw
#elif SRFFT
  real(kind=DP),private,allocatable,dimension(:) :: ftw
#elif ASLFFT
  real(kind=DP),private,allocatable,dimension(:) :: ftw
#elif FFTW3
  real(kind=DP),private,allocatable,dimension(:) :: ftw
#endif

!  include 'mpif.h'
!#ifdef FFTE
!  include 'ffte.h'
!#endif

#ifdef FFTW3
  integer, parameter :: FFTW_ESTIMATE=64
#endif
contains

#ifndef DECFFT
# ifndef MKLFFT
#  ifndef FFTW3

  subroutine screening_fft_CD_work_alloc
    integer :: nid,nnl,nnm,nnn
    nid = fft_box_size_CD(1,0);  nnl = fft_box_size_CD(1,1)
    nnm = fft_box_size_CD(2,1);  nnn = fft_box_size_CD(3,1)
#   ifdef CD_IFFT
    allocate(wlp(3*nnl+7), stat=istat)
    allocate(wmp(2*nnm*(nid+1)+7), stat=istat)
    allocate(wnp(2*nnn+7), stat=istat)
    allocate(iwork(2,nnl+nnm+nnn), stat=istat)
#   elif CD_JRCATFFT
    allocate(wlp(nnl), stat=istat)
    allocate(wmp(nnm), stat=istat)
    allocate(wnp(nnn), stat=istat)
#   elif CD_JRCATFFT_WS
    allocate(wlp(nnl), stat=istat)
    allocate(wmp(nnm), stat=istat)
    allocate(wnp(nnn), stat=istat)
#   elif ASLFFT
    allocate(trigs_CD(nnl+2*(nnm+nnn)), stat=istat)
    allocate(ifax_CD(60), stat=istat)
#   endif
  end subroutine screening_fft_CD_work_alloc
#  endif
! /* <- FFTW3 */
# endif
! /* <- MKLFFT */
#endif
! /* <- DECFFT */

  subroutine m_Screening_FFT_alloc_CD_box
#ifndef DECFFT
# ifndef MKLFFT
#  ifndef FFTW3
    integer :: nid,nnl,nnm,nnn
    nid  = fft_box_size_CD(1,0);  nnl  = fft_box_size_CD(1,1)
    nnm  = fft_box_size_CD(2,1);  nnn  = fft_box_size_CD(3,1)
#   ifdef CD_JRCATFFT_WS
    allocate(ftw(nfftp*2), stat=istat)
#   elif  CD_JRCATFFT
    allocate(ftw(nfftp), stat=istat)
#   elif  ASLFFT
    allocate(ftw(nfftp), stat=istat)
#   endif

#  endif
# endif
#endif
    allocate(afft_CD(nfftp), stat=istat)
  end subroutine m_Screening_FFT_alloc_CD_box

  subroutine m_Screening_FFT_dealloc_CD_box
#ifndef DECFFT
# ifndef MKLFFT
#  ifndef FFTW3
    if(allocated(ftw)) then
       deallocate(ftw, stat=istat)
    end if
#  endif
# endif
#endif
    deallocate(afft_CD, stat=istat)
  end subroutine m_Screening_FFT_dealloc_CD_box

  subroutine m_Screening_FFT_set_box_sizes(n_rGpv,outer_or_inner)
    integer,intent(in),dimension(3) :: n_rGpv
    integer,intent(in)              :: outer_or_inner

    integer :: i, ip
    do i = 1, 3
!  -- FFT box size for Charge Density --
       ip = n_rGpv(i)
#ifdef GOEDECKER_FFT
       if(mod(ip,2) == 1) ip = ip + 1
#endif
#ifdef _PFFT_
       call decomp3p_2(ip,fft_box_size_CD(i,1))
#else
       call decomp3_2(outer_or_inner,ip,fft_box_size_CD(i,1))
#endif
    end do

  end subroutine m_Screening_FFT_set_box_sizes

  subroutine decomp3_2(oi,i,j)
    integer, intent(in)  :: oi, i
    integer, intent(out) :: j
    integer :: ipm, iab, ik, imini
    imini = i*189/200+1

    if(oi == OUTER) then
       ipm = +1
    else
       ipm = -1
    endif
100 iab=-1
110 iab=iab+1
    ik= i+ipm*iab
    if(ik < imini) then
       ipm=1
       goto 100
    end if
    j = ik*2
120 continue
    if(mod(j,2) == 0) then
       j=j/2
       goto 120
    else if(mod(j,3) == 0) then
       j=j/3
       goto 120
    else if(mod(j,5) == 0) then
       j=j/5
       goto 120
    else if(j == 1) then
       goto 130
    else
       goto 110
    end if
130 ik= i+ipm*iab
    j  = ik*2
  end subroutine decomp3_2



#ifdef _PFFT_
  subroutine decomp3p_2(i,j)
    integer, intent(in)  :: i
    integer, intent(out) :: j

    integer :: iab, ik

    iab=-1
110 iab=iab+1
    ik = i+iab

    j  = ik*2
120 continue
    if(mod(j,2) == 0) then
       j=j/2
       goto 120
    else if(mod(j,3) == 0) then
       j=j/3
       goto 120
    else if(mod(j,5) == 0) then
       j=j/5
       goto 120
    else if(j == 1) then
       goto 130
    else
       goto 110
    end if
130 ik= i+iab
    j  = ik*2
  end subroutine decomp3p_2
#endif



  subroutine m_Screening_FFT_setup(inversion_symmetry,paramset)
    integer, intent(in)           :: inversion_symmetry
    logical, intent(in), optional :: paramset

    integer :: nfft_t
    integer :: id_sname = -1

#ifdef ASLFFT
    if(mod(fft_box_size_CD(1,1),4) == 2) then
       fft_box_size_CD(1,0) = fft_box_size_CD(1,1) + 4
    else
       fft_box_size_CD(1,0) = fft_box_size_CD(1,1) + 2
    end if
    fft_box_size_CD(2:3,0) = fft_box_size_CD(2:3,1) + 1
#else
    fft_box_size_CD(1,0) = fft_box_size_CD(1,1) + 2
    fft_box_size_CD(2:3,0) = fft_box_size_CD(2:3,1)

    nfftp=   product(fft_box_size_CD(1:3,0)) * (2-inversion_symmetry)
#endif

#ifndef DECFFT
# ifndef MKLFFT
#  ifndef FFTW3
    call screening_fft_CD_work_alloc     ! <cw[123],cw[123]_pstrn,wlp,wmp,wnp> are allocated
#  endif
# endif
#endif

! Initialization of the Charge-Density FFT
#ifndef DECFFT
# ifndef MKLFFT
       call m_Screening_FFT_alloc_CD_box()       ! <ftw> is allocated
# endif
#endif
       call init_scr_fft_coef_arrays_CD()

#ifndef DECFFT
# ifndef MKLFFT
       call m_Screening_FFT_dealloc_CD_box()
# endif
#endif

  contains

    subroutine init_scr_fft_coef_arrays_CD
      integer :: id, nl, nm, nn, ierr
#ifdef ASLFFT
      integer :: id2, id3
#endif
   !  ---> FFT for Charge density
      id = fft_box_size_CD(1,0)
#ifdef ASLFFT
      id2 = fft_box_size_CD(2,0)
      id3 = fft_box_size_CD(3,0)
#endif
      nl = fft_box_size_CD(1,1)
      nm = fft_box_size_CD(2,1)
      nn = fft_box_size_CD(3,1)
#ifdef CD_IFFT
      call r3fft(afft_CD,id,nl,nm,nn,wlp,wmp,wnp,0,0,1,ftw,ierr)
#elif DECFFT
      call cd_init(id,nl,nm,nn)
#elif MKLFFT
      Status = DftiCreateDescriptor( dh_CD, DFTI_DOUBLE, &
           & DFTI_REAL, 3, fft_box_size_CD(1,1))
      Status = DftiCreateDescriptor( dh_CD_back, DFTI_DOUBLE, &
           & DFTI_REAL, 3, fft_box_size_CD(1,1))
      Status = DftiSetValue( dh_CD, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
      Status = DftiSetValue( dh_CD_back, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
      Status = DftiSetValue( dh_CD, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX)
      Status = DftiSetValue( dh_CD_back, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX)
   !!$Status = DftiSetValue( dh_CD, DFTI_FORWARD_DOMAIN, DFTI_COMPLEX)
   !!$Status = DftiCopyDescriptor( dh_CD, dh_CD_back)
      strides(1) = 0
      strides(2) = 1
      strides(3) = id/2
      strides(4) = id/2*fft_box_size_CD(2,0)
      Status = DftiSetValue(dh_CD, DFTI_INPUT_STRIDES, strides)
      Status = DftiSetValue(dh_CD_back, DFTI_OUTPUT_STRIDES, strides)
      strides(1) = 0
      strides(2) = 1
      strides(3) = id
      strides(4) = id*fft_box_size_CD(2,0)
      Status = DftiSetValue(dh_CD, DFTI_OUTPUT_STRIDES, strides)
      Status = DftiSetValue(dh_CD_back, DFTI_INPUT_STRIDES, strides)
      Status = DftiCommitDescriptor( dh_CD )
      Status = DftiCommitDescriptor( dh_CD_back )
#elif FFTW3
      ! Forward FFT
      call dfftw_plan_dft_c2r_3d(plan_CD(1) &
           &                       ,nl,nm,nn &
           &                       ,afft_CD(1),afft_CD(1) &
           &                       ,FFTW_ESTIMATE)
      ! Inverse FFT
      call dfftw_plan_dft_r2c_3d(plan_CD(2) &
           &                       ,nl,nm,nn &
           &                       ,afft_CD(1),afft_CD(1) &
           &                       ,FFTW_ESTIMATE)
#elif ASLFFT
      call dfr3fb(nl,nm,nn,afft_CD,id,id2,id3,0,ifax_CD,trigs_CD,ftw,ierr)
#else
# ifdef CD_JRCATFFTg
      call jrcat_r3ft(afft_CD,ftw,id,nl,nm,nn,wlp,wmp,wnp,0,0,0,0)
# endif
#endif
    end subroutine init_scr_fft_coef_arrays_CD

  end subroutine m_Screening_FFT_setup


  subroutine m_Screening_FFT_CD_inverse(nfout,afft)  ! G space --> R space
    integer, intent(in)          :: nfout
    real(kind=DP), intent(inout) :: afft(1:nfftp)
    integer :: idp,nlp, nmp, nnp, ier1 = 0
    integer :: id_sname = -1
#ifdef MKLFFT
    integer :: i
#endif
#ifdef ASLFFT
    integer :: idp2,idp3
#endif

#ifdef GOEDECKER_FFT
    integer :: nd2p, nd3p, inzee, isign
#endif

#ifdef CD_JRCATFFTg
    integer :: kc1 = 0, kc2 = 0, kc3 = 0
#endif

    idp = fft_box_size_CD(1,0)
#ifdef ASLFFT
    idp2 = fft_box_size_CD(2,0)
    idp3 = fft_box_size_CD(3,0)
#endif
    nlp = fft_box_size_CD(1,1)
    nmp = fft_box_size_CD(2,1)
    nnp = fft_box_size_CD(3,1)
#ifdef GOEDECKER_FFT
    nd2p = fft_box_size_CD(2,0)
    nd3p = fft_box_size_CD(3,0)
#endif

    afft_CD = afft
#ifdef CD_JRCATFFTg
    call jrcat_r3ft(afft_CD,ftw,idp,nlp,nmp,nnp,wlp,wmp,wnp,kc1,kc2,kc3&
            & ,-2)
#endif
#ifdef DECFFT
    call cd_fft(1,idp,nmp,afft_CD)
#elif MKLFFT
    allocate(cwork_mkl(nfftp/2))
    Status = DftiComputeForward( dh_CD_back, afft_CD, cwork_mkl)
    do i=1,nfftp,2
       afft_CD(i) = dble(cwork_mkl(i/2+1))
       afft_CD(i+1) = dimag(cwork_mkl(i/2+1))
    end do
    deallocate(cwork_mkl)
#elif FFTW3
    call dfftw_execute_dft_r2c(plan_CD(2),afft_CD(1),afft_CD(1))
#endif
#ifdef ASLFFT
    call dfr3bf(nlp,nmp,nnp,afft_CD,idp,idp2,idp3,+1,ifax_CD,trigs_CD,ftw,ier1)
!!$    call sxfft_zero_padding(afft_CD,nlp,nmp,nnp,idp,idp2,idp3)
    call sxfft_zero_padding2(afft_CD,nlp,nmp,nnp,idp,idp2,idp3)
#endif
#ifdef CD_IFFT
    call r3fft(afft_CD,idp,nlp,nmp,nnp,wlp,wmp,wnp,iopt,-1,1,iwork,ier1)
#endif
    if (ier1 /= 0) then
       stop
    end if

    afft(1:nfftp) = afft_CD(1:nfftp)

  end subroutine m_Screening_FFT_CD_inverse

#ifdef ASLFFT
subroutine sxfft_zero_padding2(afft,nl,nm,nn,nd1,nd2,nd3)
  implicit none
  real(kind=DP), intent(out) :: afft(*)
  integer, intent(in) :: nl,nm,nn
  integer, intent(in) :: nd1,nd2,nd3

  integer :: i,j,k,ip
  integer :: nlh,ndh

  integer :: imax
  imax = nd1*nd2*nd3

  nlh = nl/2
  ndh = nd1/2
  do j = nlh+2, ndh
     do i = 1, nd2*nn
        ip = ndh*(i-1) + j
        afft(ip*2-1) = 0.d0
        afft(ip*2) = 0.d0
     end do
  end do
  do j = nm+1, nd2
     do k = 1, nd3
        do i = 1, ndh
           ip = i + ndh*(j-1) + ndh*nd2*(k-1)
           afft(ip*2-1) = 0.d0
           afft(ip*2) = 0.d0
        end do
     end do
  end do
  do k = nn+1, nd3
     do i = 1, ndh*nd2
        ip = i + ndh*nd2*(k-1)
        afft(ip*2-1) = 0.d0
        afft(ip*2) = 0.d0
     end do
  end do

end subroutine sxfft_zero_padding2
#endif

end module m_Screening_FFT
