!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 635 $)
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

module m_FFT
! $Id: m_FFT.F90 635 2021-02-26 07:16:10Z jkoga $
!
  use m_Timing,            only : tstatc0_begin, tstatc0_end
  use m_Files,             only : nfout
  use m_Control_Parameters,only : kimg, ipri, iprinegativecharge, iprinegativecharge0 &
       &                        , max_warnings_negativecharge, sw_positron, ldx, ldy, ldz
  use m_Const_Parameters,  only : DP,CMPLDP, ON, OFF, INVERSE, DIRECT, OUTER, INNER &
       &                        , ELECTRON, POSITRON, YES, NO
  use m_Parallelization,   only : MPI_CommGroup,ista_sfftp,iend_sfftp,ista_sfftph,iend_sfftph &
       &                        , ista_fftp, iend_fftp &
       &                        , mype,npes,ierr,nel_sfftp,idisp_sfftp &
       &                        , myrank_cdfft,myrank_ggacmp,nrank_ggacmp,npes_cdfft &
       &                        , nel_fftp, idisp_fftp, mpi_cdfft_world
#ifdef MKLFFT
  use MKL_DFTI
  use mklfft
#endif

! ================================= added by K. Tagami ============== 11.0
  use m_Const_Parameters,       only : CMPLDP
! =================================================================== 11.0

! ======= KT_add ========================================== 13.0F
  use m_Control_Parameters, only : sw_hybrid_functional
! ========================================================= 13.0F
  use mpi

  implicit none
  integer, dimension(3,0:1) :: fft_box_size_WF, fft_box_size_CD
  integer, dimension(3,0:1) :: fft_box_size_WF_prev
  integer, dimension(3,0:1) :: fft_box_size_CD_prev
  integer, dimension(3,0:0) :: fft_box_size_CD_c, fft_box_size_CD_nonpara
  integer, dimension(3,0:0) :: fft_box_size_CD_nonpara_prev
  integer, dimension(3,0:1) :: fft_box_size_CD_exx,fft_box_size_CD_exx_nonpara
  integer               :: nfft, nfftp, nfftps, nfftp_nonpara, nfftp_nonpara_prev

  integer               :: nfftp_exx_nonpara

  integer               :: ny_d, nz_d, ly_d, lx, ly, lz_d, lz &
       &                 , ny_ds,nz_ds,ly_ds,lxs,lys,lz_ds,lzs

  ! ------- Positron start
  integer, dimension(3,0:1) :: fft_box_size_pWF
  integer ::                   nfft_pstrn
  ! ------- Positron end

! ======= KT_add =================================== 13.0F
  integer, dimension(3,0:1) :: fft_box_size_EXX
  integer                   :: nfft_exx
! ================================================== 13.0F

#ifdef _MPIFFT_
  integer :: sw_mpifft = ON
#else
  integer :: sw_mpifft = OFF
#endif

  integer, private      :: CD_setup_is_done = NO
  integer, private      :: CD_setup_is_done_EXX = NO

  real(kind=DP),private,allocatable,dimension(:)           :: afft_CD
  real(kind=DP),private,allocatable,dimension(:)           :: afft_CD_exx
#ifdef VPP
  integer, private, parameter :: ncache=0
#elif DEC
  integer, private, parameter :: ncache=6*1024
#elif HP
  integer, private, parameter :: ncache=1*1024
#elif SUN
  integer, private, parameter :: ncache=1*1024
#elif ONYX
  integer, private, parameter :: ncache=1*1024
#elif IRIX64
  integer, private, parameter :: ncache=1*1024
#elif CRAY
  integer, private, parameter :: ncache=1*1024
#elif HIUX
  integer, private, parameter :: ncache=1*1024
#elif SX
  integer, private, parameter :: ncache=0
#elif SP2
  integer, private, parameter :: ncache=1*1024
#elif Linux
  integer, private, parameter :: ncache=6*1024
#endif

! contains
#ifdef JRCATFFT_WS
include "m_FFT_type1_jrcatfft.F90"
#define _INCLUDE_EXX_
#elif FFTW3
include "m_FFT_type4_fftw3_rev.F90"
#define _INCLUDE_EXX_
#elif MKLFFT
include "m_FFT_type3_mklfft.F90"
#elif ASLFFT
#ifdef _MPIFFT_
#ifdef NEC_TUNE_FFT
include "m_FFT_type5_aslfft+mpifft.hybrid.f90"
#else
include "m_FFT_type5_aslfft+mpifft.mpi.f90"
#endif
#else
#ifdef NEC_TUNE_FFT
include "m_FFT_type5_aslfft.hybrid.f90"
#else
include "m_FFT_type5_aslfft.mpi.f90"
#endif
#endif
#elif WF_JRCATFFT
include "m_FFT_type6_jrcatfft+ifft.F90"
#elif GOEDECKER_FFT
include "m_FFT_type8_goedecker.F90"
#elif FFTE
include "m_FFT_type9_ffte.F90"
#endif
! ---------------------

  ! 6. (common) wd_FFTboxsizes(nfout)
  ! 7. (common) m_FFT_set_box_sizes(n_rGv,n_rGpv,n_rGv_pstrn,outer_or_inner)
  ! 8. (common) m_FFT_wd_box_sizes(nf_bin)
  ! 9. (common) decomp3(oi,i,j)
  !10. (common) decomp3p(i,j)
  !13. (common) m_FFT_Vlocal_W(afft,bfft)
  !21. (common) m_FFT_check_of_negative_CD(npes_in,ista_fftp,iend_fftp,ista_fftph,iend_fftph &
  !22. (common) m_FFT_W_Vlocal_W(electron_or_positron,nabfft,afft,bfft,eg)
  !28. (common) m_FFT_set_cdata(aa,a)
  !29. (common) m_FFT_set_scdata(aa,a)
  !30. (common) m_FFT_coef_CD_integration(f2or1)
  !31. (common) m_FFT_set_cdata
  !32. (common) m_FFT_set_scdata(aa,a)


  subroutine wd_FFTboxsizes(nfout)
    integer, intent(in) :: nfout


    write(nfout,*) ' ---- (WF)FFT size ----'
    write(nfout,'("     (nfft, fft_box_size_WF)")')
    write(nfout,100) nfft
    write(nfout,101) fft_box_size_WF(1:3,0)
    write(nfout,102) fft_box_size_WF(1:3,1)

    write(nfout,*) ' ---- (CD)FFT size ----'
    write(nfout,'("     (nfftps, fft_box_size_CD)")')
    write(nfout,100) nfftps
    write(nfout,101) fft_box_size_CD(1:3,0)
    write(nfout,102) fft_box_size_CD(1:3,1)


#ifdef _MPIFFT_
    write(nfout,*) ' ---- (CD)FFT size ----'
    write(nfout,'("     (nfftps, fft_box_size_CD_nonpara)")')
    write(nfout,100) nfftp_nonpara
    write(nfout,101) fft_box_size_CD_nonpara(1:3,0)
#endif

    if(sw_positron /= OFF) then
       write(nfout,*) ' ---- (pWF)FFT size ----'
       write(nfout,'("     (nfft_pstrn, fft_box_size_pWF)")')
       write(nfout,100) nfft_pstrn
       write(nfout,101) fft_box_size_pWF(1:3,0)
       write(nfout,102) fft_box_size_pWF(1:3,1)
    end if

    write(nfout,*) ' ---- (CD_c)FFT size ----'
    write(nfout,'("     (nfftp, fft_box_size_CD_c)")')
    write(nfout,100) nfftp
    write(nfout,101) fft_box_size_CD_c(1:3,0)

100 format(" FFT total elements(including work area) = ",i10)
101 format(" FFT box adjustable dimension            = ",3i5)
102 format(" FFT box real dimension                  = ",3i5)
  end subroutine wd_FFTboxsizes

! =============== KT_add ============================= 13.0F
  subroutine wd_FFTboxsizes_exx(nfout)
    integer, intent(in) :: nfout

    write(nfout,*) ' ---- (Exx)FFT size ----'
    write(nfout,'("     (nfft_exx, fft_box_size_Exx)")')
    write(nfout,100) nfft_exx
    write(nfout,101) fft_box_size_EXX(1:3,0)
    write(nfout,102) fft_box_size_EXX(1:3,1)

    write(nfout,*) ' ---- (Exx)FFT_CD size ----'
    write(nfout,'("     (nfftp_exx_nonpara, fft_box_size_CD_Exx)")')
    write(nfout,100) nfftp_exx_nonpara
    write(nfout,101) fft_box_size_CD_EXX(1:3,0)
    write(nfout,102) fft_box_size_CD_EXX(1:3,1)

100 format(" FFT total elements(including work area) = ",i10)
101 format(" FFT box adjustable dimension            = ",3i5)
102 format(" FFT box real dimension                  = ",3i5)

  end subroutine wd_FFTboxsizes_exx
! ==================================================== 13.0F

  subroutine m_FFT_set_box_sizes(n_rGv,n_rGpv,n_rGv_pstrn,outer_or_inner)
    integer,intent(in),dimension(3) :: n_rGv,n_rGpv,n_rGv_pstrn
    integer,intent(in)              :: outer_or_inner

    integer :: i, ip
    fft_box_size_WF_prev = fft_box_size_WF
    do i = 1, 3
!  -- FFT box size for Wave Functions --
       ip = n_rGv(i)
       if(sw_avoiding_odd_fftbox == ON) then
          if(mod(ip,2) == 1) ip = ip + 1
       endif
       call decomp3(outer_or_inner,ip,fft_box_size_WF(i,1))
!  -- FFT box size for Charge Density --
       ip = n_rGpv(i)
       if(sw_avoiding_odd_fftbox == ON) then
          if(mod(ip,2) == 1) ip = ip + 1
       end if

#ifdef _PFFT_
       call decomp3p(ip,fft_box_size_CD(i,1))
#else
       call decomp3(outer_or_inner,ip,fft_box_size_CD(i,1))
#endif
    end do
    if(ipri>=1) then
       write(nfout,*) ' --- fft_box_size <<m_FFT_set_box_sizes>> ---'
       write(nfout,'("  WF ",3i5)') fft_box_size_WF(1:3,1)
       write(nfout,'("  CD ",3i5)') fft_box_size_CD(1:3,1)
    end if
    ! ------- Positron start
    if(sw_positron /= OFF) then
       do i = 1, 3
          !  -- FFT box size for positron Wave Functions --
          ip = n_rGv_pstrn(i)
          if(sw_avoiding_odd_fftbox == ON) then
             if(mod(ip,2) == 1) ip = ip + 1
          end if
          call decomp3(outer_or_inner,ip,fft_box_size_pWF(i,1))
       end do
       if(ipri>=1) write(nfout,'(" pWF ",3i5)') fft_box_size_pWF(1:3,1)
    end if
    ! ------- Positron end

  end subroutine m_FFT_set_box_sizes

  subroutine m_FFT_set_box_size_prev( n_rGv_prev, outer_or_inner )
    integer,intent(in),dimension(3) :: n_rGv_prev
    integer,intent(in)              :: outer_or_inner

    integer :: i, ip

    do i = 1, 3
       ip = n_rGv_prev(i)
       if(sw_avoiding_odd_fftbox == ON) then
          if(mod(ip,2) == 1) ip = ip + 1
       end if
       call decomp3(outer_or_inner,ip,fft_box_size_WF_prev(i,1))
    end do
    if(ipri>=1) write(nfout,'(" prev ",3i5)') fft_box_size_WF_prev(1:3,1)

  end subroutine m_FFT_set_box_size_prev

  subroutine m_FFT_set_box_size_cd_prev( n_rGpv_prev, outer_or_inner )
    integer,intent(in),dimension(3) :: n_rGpv_prev
    integer,intent(in)              :: outer_or_inner

    integer :: i, ip, ipad

    do i = 1, 3
       !  -- FFT box size for Wave Functions (exx) --
       ip = n_rGpv_prev(i)

       if(sw_avoiding_odd_fftbox == ON) then
          if(mod(ip,2) == 1) ip = ip + 1
       end if
       call decomp3(outer_or_inner,ip,fft_box_size_CD_prev(i,1))
    end do
    if(ipri>=1) write(nfout,'(" exx_CD ",3i5)') fft_box_size_CD_prev(1:3,1)
    if(kimg == 1) then  ! kimg == 1
       ipad = 2
    else if(kimg == 2) then ! kimg == 2
       ipad = 0
    end if

    fft_box_size_CD_nonpara_prev(1,0)   = fft_box_size_CD_prev(1,1) + ipad
    fft_box_size_CD_nonpara_prev(2:3,0) = fft_box_size_CD_prev(2:3,1)

    fft_box_size_CD_prev(1,0)   = fft_box_size_CD_prev(1,1) + ipad
    fft_box_size_CD_prev(2:3,0) = fft_box_size_CD_prev(2:3,1)

    nfftp_nonpara_prev  = product(fft_box_size_CD_nonpara_prev(1:3,0)) * kimg

  end subroutine m_FFT_set_box_size_cd_prev

! ======= KT_add ======================================= 13.0F
  subroutine m_FFT_set_box_size_exx( n_rGv_exx, outer_or_inner )
    integer,intent(in),dimension(3) :: n_rGv_exx
    integer,intent(in)              :: outer_or_inner

    integer :: i, ip

    if(sw_hybrid_functional /= OFF) then
       do i = 1, 3
          !  -- FFT box size for Wave Functions (exx) --
          ip = n_rGv_exx(i)

          if(sw_avoiding_odd_fftbox == ON) then
             if(mod(ip,2) == 1) ip = ip + 1
          end if
          call decomp3(outer_or_inner,ip,fft_box_size_exx(i,1))
       end do
       if(ipri>=1) write(nfout,'(" exx ",3i5)') fft_box_size_exx(1:3,1)
    end if

  end subroutine m_FFT_set_box_size_exx

  subroutine m_FFT_set_box_size_cd_exx( n_rGpv_exx, outer_or_inner )
    integer,intent(in),dimension(3) :: n_rGpv_exx
    integer,intent(in)              :: outer_or_inner

    integer :: i, ip, ipad

    if(sw_hybrid_functional /= OFF) then
       do i = 1, 3
          !  -- FFT box size for Wave Functions (exx) --
          ip = n_rGpv_exx(i)

          if(sw_avoiding_odd_fftbox == ON) then
             if(mod(ip,2) == 1) ip = ip + 1
          end if
          call decomp3(outer_or_inner,ip,fft_box_size_CD_exx(i,1))
       end do
       if(ipri>=1) write(nfout,'(" exx_CD ",3i5)') fft_box_size_CD_exx(1:3,1)
       if(kimg == 1) then  ! kimg == 1
          ipad = 2
       else if(kimg == 2) then ! kimg == 2
          ipad = 0
       end if

       fft_box_size_CD_exx_nonpara(1,0)   = fft_box_size_CD_exx(1,1) + ipad
       fft_box_size_CD_exx_nonpara(2:3,0) = fft_box_size_CD_exx(2:3,1)

       fft_box_size_CD_exx(1,0)   = fft_box_size_CD_exx(1,1) + ipad
       fft_box_size_CD_exx(2:3,0) = fft_box_size_CD_exx(2:3,1)

       nfftp_exx_nonpara  = product(fft_box_size_CD_exx_nonpara(1:3,0)) * kimg
    end if

  end subroutine m_FFT_set_box_size_cd_exx
! ====================================================== 13.0F

  subroutine m_FFT_wd_box_sizes(nf_bin)
    integer, intent(in) :: nf_bin
    if(mype == 0) &
         & write(nf_bin) fft_box_size_WF(1:3,1),fft_box_size_CD(1:3,1)
  end subroutine m_FFT_wd_box_sizes

  subroutine decomp3(oi,i,j)
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
  end subroutine decomp3

#ifdef _PFFT_
  subroutine decomp3p(i,j)
    integer, intent(in)  :: i
    integer, intent(out) :: j

    integer :: iab, ik

    iab=-1
110 iab=iab+1
    ik = i+iab
    if(mod(ik,npes) /= 0) goto 110
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
  end subroutine decomp3p
#endif

  subroutine m_FFT_Vlocal_W(afft,bfft)
    real(kind=DP), intent(in),    dimension(nfft) :: afft
    real(kind=DP), intent(inout), dimension(nfft) :: bfft

    integer i

    do i = 1, nfft-1, 2
       bfft(i)   = afft(i) * bfft(i)
       bfft(i+1) = afft(i) * bfft(i+1)
    end do
    if(mod(nfft,2) == 1) bfft(nfft) = afft(nfft)*bfft(nfft)
  end subroutine m_FFT_Vlocal_W

! =============================== added by K. Tagami ============= 11.0
  subroutine m_FFT_Vlocal_W_noncl( afft_kt, bfft_kt, ndim_chgpot, ndim_spinor )
    integer,  intent(in) :: ndim_spinor, ndim_chgpot
    real(kind=DP), intent(inout) :: afft_kt( nfft,ndim_chgpot )
    real(kind=DP), intent(inout) :: bfft_kt( nfft,ndim_spinor )

    real(kind=DP) :: ctmp_r, ctmp_i

    real(kind=DP), allocatable :: ccfft( :,: )
    integer i

    allocate( ccfft( nfft, ndim_spinor ) ); ccfft = 0.0d0

    do i = 1, nfft-1, 2
       ctmp_r = afft_kt(i,1) *bfft_kt(i,1) - afft_kt(i+1,1) *bfft_kt(i+1,1) &
            & + afft_kt(i,2) *bfft_kt(i,2) - afft_kt(i+1,2) *bfft_kt(i+1,2)
       ctmp_i = afft_kt(i+1,1) *bfft_kt(i,1) + afft_kt(i,1) *bfft_kt(i+1,1) &
            & + afft_kt(i+1,2) *bfft_kt(i,2) + afft_kt(i,2) *bfft_kt(i+1,2)
!       afft_kt(i,  1) = ctmp_r
!       afft_kt(i+1,1) = ctmp_i

       ccfft( i,1 ) = ctmp_r;  ccfft( i+1, 1 ) = ctmp_i
    end do
!
    do i = 1, nfft-1, 2
       ctmp_r = afft_kt(i,3) *bfft_kt(i,1) - afft_kt(i+1,3) *bfft_kt(i+1,1) &
            & + afft_kt(i,4) *bfft_kt(i,2) - afft_kt(i+1,4) *bfft_kt(i+1,2)
       ctmp_i = afft_kt(i+1,3) *bfft_kt(i,1) + afft_kt(i,3) *bfft_kt(i+1,1) &
            & + afft_kt(i+1,4) *bfft_kt(i,2) + afft_kt(i,4) *bfft_kt(i+1,2)
!       afft_kt(i,  2) = ctmp_r
!       afft_kt(i+1,2) = ctmp_i

       ccfft( i,2 ) = ctmp_r;  ccfft( i+1, 2 ) = ctmp_i
    end do

!    bfft_kt(:,1) = afft_kt(:,1)
!    bfft_kt(:,2) = afft_kt(:,2)

    bfft_kt(:,1) = ccfft(:,1)
    bfft_kt(:,2) = ccfft(:,2)

    deallocate( ccfft )

  end subroutine m_FFT_Vlocal_W_noncl
! ========================================================================= 11.0

  subroutine m_FFT_check_of_negative_CD(npes_in,ista_fftp,iend_fftp,ista_fftph,iend_fftph &
       &                               ,afft,nfout,nspin,ispin)
    integer, intent(in)          :: npes_in, ista_fftp, iend_fftp, ista_fftph, iend_fftph
    real(kind=DP), intent(inout) :: afft(ista_fftp:iend_fftp)
    integer, intent(in)          :: nfout,nspin,ispin

    integer                   :: ip, i,j,k, nlp,nmp,nnp,idp,nd2p,nd3p, nlph, nd1h
    integer                   :: icwarn_NEGA, icwarn_IMAG, icwarn_NEGA2,icwarn_IMAG2 &
         &                     , icwarn_NEGA_total,icwarn_IMAG_total, icwarn_NEGA_0, icwarn_IMAG_0, isw
    integer, allocatable, dimension(:)       :: istatus_afft_nega, istatus_afft_imag
!!$    integer, parameter        :: MXWARN = 10
    real(kind=DP),parameter   :: chgdel = 8.d-5
    real(kind=DP),parameter   :: D_min  = 1.d-40
!     Revised by T. Yamasaki, 14th Mar. 2005
!         Write statements in if-blocks in a do-loop have been excluded for vector processors.
    integer :: ipad1, ipad2, ipad3

    integer ::        kk, nstart, nend, itmp

    if(kimg ==  1) then
       idp = nfftp+10
    end if

    nlp = fft_box_size_CD(1,1)
    nmp = fft_box_size_CD(2,1)
    nnp = fft_box_size_CD(3,1)
    if(npes_in == npes) then
       nd2p = fft_box_size_CD(2,0)
       nd3p = fft_box_size_CD(3,0)
       if(kimg == 2) idp = fft_box_size_CD(1,0)
    else if(npes_in == npes_cdfft) then
       nd2p = fft_box_size_CD_c(2,0)
       nd3p = fft_box_size_CD_c(3,0)
       if(kimg == 2) idp = fft_box_size_CD_c(1,0)
    end if

! -- Zero padding
    if(kimg == 1) then
       ipad1 = fft_box_size_CD(1,0) - fft_box_size_CD(1,1)
       ipad2 = fft_box_size_CD(2,0) - fft_box_size_CD(2,1)
       ipad3 = fft_box_size_CD(3,0) - fft_box_size_CD(3,1)

! === DEBUG by tkato 2011/09/12 ================================================
       nlph = nlp/2
       nd1h = fft_box_size_CD(1,0)/2
! ==============================================================================
       if(ipad1 >= 4) then
! === DEBUG by tkato 2011/09/12 ================================================
!         nlph = nlp/2
!         nd1h = fft_box_size_CD(1,0)/2
! ==============================================================================
          do j = nlph+2, nd1h
             do i = 1, nd2p*nnp
                itmp = nd1h*(i-1)+j
                if(ipri >= 2) write(nfout,*) 'ZERO1: i=',itmp
                ip = 2*(nd1h*(i-1)+j)-1
                if(ip >= ista_fftp .and. ip <= iend_fftp) then
                   afft(ip)   = 0.d0
                   afft(ip+1) = 0.d0
                end if
             end do
          end do
       end if
       if(ipad2 >= 2) then
          do j = nmp+1, nd2p
             do k = 1, nnp
                do i = 1, nlph
                   itmp = i + nd1h*(j-1) + nd1h*nd2p*(k-1)
                   if(ipri >= 2) write(nfout,*) 'ZERO2: i=',itmp
                   ip = 2*(i + nd1h*(j-1) + nd1h*nd2p*(k-1))-1
                   if(ip >= ista_fftp .and. ip <= iend_fftp) then
                      afft(ip)   = 0.d0
                      afft(ip+1) = 0.d0
                   end if
                end do
             end do
          end do
       end if
       if(ipad3 >= 2) then
          do k = nnp+1, nd3p
             do i = 1, nd1h*nd2p
                itmp = nd1h*nd2p*(k-1) + i
                if(ipri >= 2) write(nfout,*) 'ZERO3: i=',itmp
                ip = 2*(nd1h*nd2p*(k-1) + i) - 1
                if(ip >= ista_fftp .and. ip <= iend_fftp) then
                   afft(ip) = 0.d0
                   afft(ip+1) = 0.d0
                end if
             end do
          end do
       end if
    end if

    if(sw_mpifft == ON) then
       if(kimg == 2) then
          do j = nlp+1, idp       ! x
             do i = 1, nd2p*nd3p  ! y,z
                ip = 2*(idp*(i-1)+j)-1
                if(ip >= ista_fftp .and. ip+1 <= iend_fftp) then
                   afft(ip)   = 0.d0
                   afft(ip+1) = 0.d0
                end if
             end do
          end do
!!$       lz_d   = fft_box_size_CD_c(3,0)/npes_cdfft
          do j = nmp+1, nd2p      ! y
             do k = 1, lz_d       ! z
                do i = 1, nlp     ! x
                   ip = 2*(i+idp*(j-1)+idp*nd2p*(k-1)+idp*nd2p*lz_d*myrank_cdfft)-1
                   if(ip >= ista_fftp .and. ip+1 <= iend_fftp) then
                      afft(ip)   = 0.d0
                      afft(ip+1) = 0.d0
                   end if
                end do
             end do
          end do
          nstart = min(lz_d*myrank_cdfft+nz_d+1,nnp+myrank_cdfft*(lz_d-nz_d)+1) ! = nnp + (isize-1)*ld + 1
          nend   = (myrank_cdfft+1)*lz_d
          do k = nstart, nend     ! z
             kk = k - lz_d*myrank_cdfft + 1
             do j = 1, nmp        ! y
                do i = 1, nlp     ! x
                   ip = 2*(i+idp*(j-1)+idp*nd2p*(kk-1)+idp*nd2p*lz_d*myrank_cdfft)-1
                   if(ip >= ista_fftp .and. ip+1 <= iend_fftp) then
                      afft(ip)   = 0.d0
                      afft(ip+1) = 0.d0
                   end if
                end do
             end do
          end do
       end if
    else
!!$#ifndef ASLFFT
       if(sw_zero_padding == ON) then
          if(kimg == 2) then
             do j = nlp+1, idp        ! x
                do i = 1, nmp*nnp     ! y,z
                   ip = 2*(idp*(i-1)+j)-1
                   if(ip >= ista_fftp .and. ip <= iend_fftp) then
                      afft(ip)   = 0.d0
                      afft(ip+1) = 0.d0
                   end if
                end do
             end do
             do j = nmp+1, nd2p       ! y
                do k = 1, nnp         ! z
                   do i = 1, nlp      ! x
                      ip = 2*(i + idp*(j-1) + idp*nd2p*(k-1))-1
                      if(ip >= ista_fftp .and. ip <= iend_fftp) then
                         afft(ip)   = 0.d0
                         afft(ip+1) = 0.d0
                      end if
                   end do
                end do
             end do
          end if
       end if
!!$#endif
    end if

    icwarn_NEGA = 0
    icwarn_IMAG = 0
    if(iprinegativecharge0 < 1 .or. max_warnings_negativecharge <= 0) then
#ifdef NEC_TUNE_MXCP
!CDIR INNER
#endif
       do i = ista_fftph, iend_fftph
          ip = i*2-1
          if(afft(ip) <= -chgdel .and. mod(ip+1,idp) /= 0) then
             icwarn_NEGA = icwarn_NEGA + 1
             afft(ip) = D_min
          else if (afft(ip) <= 0.d0) then
             afft(ip) = D_min
          end if
          if( (abs(afft(ip+1)) > chgdel).and. (mod(ip+1,idp) /= 0) ) then
             icwarn_IMAG = icwarn_IMAG + 1
             afft(ip+1) = 0.d0
          else if (abs(afft(ip+1)) > 0.d0) then
             afft(ip+1) = 0.d0
          end if
       end do

    else if(iprinegativecharge0 >= 1 .and. max_warnings_negativecharge > 0) then
#ifdef _VECTOR_TUNING_
       allocate(istatus_afft_nega(ista_fftph:iend_fftph)); istatus_afft_nega = 0;
       allocate(istatus_afft_imag(ista_fftph:iend_fftph)); istatus_afft_imag = 0;

       icwarn_NEGA2 = 0
       icwarn_IMAG2 = 0
#ifdef NEC_TUNE_MXCP
!CDIR INNER
#endif
       do i = ista_fftph, iend_fftph
          ip = i*2-1
          if(afft(ip) <= -chgdel .and. mod(ip+1,idp) /= 0) then
             istatus_afft_nega(i) = 1
             icwarn_NEGA = icwarn_NEGA + 1
          else if (afft(ip) <= 0.d0) then
             icwarn_NEGA2 = icwarn_NEGA2 + 1
             istatus_afft_nega(i) = 2
          end if
          if( (abs(afft(ip+1)) > chgdel).and. (mod(ip+1,idp) /= 0) ) then
             istatus_afft_imag(i) = 1
             icwarn_IMAG = icwarn_IMAG + 1
          else if (abs(afft(ip+1)) > 0.d0) then
             icwarn_IMAG2 = icwarn_IMAG2 + 1
             istatus_afft_imag(i) = 2
          end if
       end do

       if(iprinegativecharge >= 1) then
          if(icwarn_NEGA >= 1) then
             do i = ista_fftph, iend_fftph
                if(istatus_afft_nega(i) == 1) exit
             end do
             icwarn_NEGA_0 = i
             if(iprinegativecharge >= 2) &
                  & write(nfout,'(" ista_fftph, iend_fftph, icwarn_NEGA_0 = ",3i8)') ista_fftph, iend_fftph, icwarn_NEGA_0

! ============================ modified by K. Tagami ================ 11.0
!             if(nspin == 2)  write(nfout,'(" #spin = ",i3)') ispin
!
             if ( nspin == 100 ) then
                write(nfout,'(" #loop = ",i3)') ispin
             else if ( nspin == 2 )  then
                write(nfout,'(" #spin = ",i3)') ispin
             endif
! ============================ modified by K. Tagami ================ 11.0

             j = 0
             do i = icwarn_NEGA_0, iend_fftph
                ip = i*2-1
                if(istatus_afft_nega(i) == 1) then
                   j = j+1
                   if(j > max_warnings_negativecharge) exit
                   write(nfout,'(" *** WARN CHG.DEN = ",d15.7, " < 0.0 AT ",i8," ***")') afft(ip),i
                end if
             end do
          end if
          if(icwarn_IMAG >= 1) then
             do i = ista_fftph, iend_fftph
                if(istatus_afft_imag(i) == 1) exit
             end do
             icwarn_IMAG_0 = i
             if(iprinegativecharge >= 2) &
                  & write(nfout,'(" ista_fftph, iend_fftph, icwarn_NEGA_0 = ",3i8)') ista_fftph, iend_fftph, icwarn_IMAG_0

! =================================== modified by K. Tagami ================ 11.0
!             if(nspin == 2) write(nfout,'(" #spin = ",i3)') ispin
!
             if ( nspin == 100 ) then
                write(nfout,'(" #loop = ",i3)') ispin
             else if ( nspin == 2 ) then
                write(nfout,'(" #spin = ",i3)') ispin
             endif
! ========================================================================== 11.0

             j = 0
             do i = icwarn_IMAG_0, iend_fftph
                ip = i*2-1
                if(istatus_afft_imag(i) == 1) then
                   j = j+1
                   if(j > max_warnings_negativecharge) exit
                   write(nfout,'(" *** WARN IMAG(CHG) = ",d15.7, " > 0.0 AT ",i8," ***")') afft(ip+1),i
                end if
             end do
          end if
       end if

       if(icwarn_NEGA+icwarn_NEGA2 >=1 .or. icwarn_IMAG+icwarn_IMAG2 >=1) then
          do i = ista_fftph, iend_fftph
             ip = i*2-1
             if(istatus_afft_nega(i) >= 1) afft(ip) = D_min
             if(istatus_afft_imag(i) >= 1) afft(ip+1) = 0.d0
          end do
       end if

       deallocate(istatus_afft_imag, istatus_afft_nega)

#else
       isw = 0
       do i = ista_fftph, iend_fftph
          ip = i*2-1
          if(afft(ip) <= -chgdel .and. mod(ip+1,idp) /= 0) then
             icwarn_NEGA = icwarn_NEGA + 1
             if(icwarn_NEGA <= max_warnings_negativecharge) then

! ==================================== modified by K. Tagami ============= 11.0
!                if(icwarn_NEGA == 1 .and. nspin == 2) then
!                   if(iprinegativecharge>=1) then
!                      write(nfout,'(" #spin = ",i3)') ispin
!                      isw = 1
!                   end if
!                endif
!
                if (icwarn_NEGA == 1) then
                   if ( nspin == 100 ) then
                      if(iprinegativecharge>=1) then
                         write(nfout,'(" #loop = ",i3)') ispin
                         isw = 1
                      end if
                   else if ( nspin == 2 ) then
                      if(iprinegativecharge>=1) then
                         write(nfout,'(" #spin = ",i3)') ispin
                         isw = 1
                      end if
                   end if
                end if
! ====================================================================== 11.0

                if(iprinegativecharge >=1) &
                     & write(nfout,'(" *** WARN CHG.DEN = ",d15.7, " < 0.0 AT ",i8," ***")') afft(ip),i
             endif
             afft(ip) = D_min
          else if (afft(ip) <= 0.d0) then
             afft(ip) = D_min
          end if
          if( (abs(afft(ip+1)) > chgdel).and. (mod(ip+1,idp) /= 0) ) then
             icwarn_IMAG = icwarn_IMAG + 1
             if(icwarn_IMAG <= max_warnings_negativecharge) then
! ==================================== modified by K. Tagami ============ 11.0
!                if(icwarn_IMAG == 1 .and. nspin == 2) then
!                   if(iprinegativecharge >=1 .and. isw == 0) write(nfout,'("#spin = ",i3)') ispin
!                endif
!
                if (icwarn_IMAG == 1 ) then
                   if ( nspin == 100 ) then
                      if(iprinegativecharge >=1 .and. isw == 0) &
                           &         write(nfout,'(" #loop = ",i3)') ispin
                   else if ( nspin == 2 ) then
                      if(iprinegativecharge >=1 .and. isw == 0) &
                           &         write(nfout,'(" #spin = ",i3)') ispin
                   endif
                endif
! ======================================================================== 11.0

                if(iprinegativecharge>=1) &
                     & write(nfout,'(" *** WARN IMAG(CHG) = ",d15.7, " > 0.0 AT ",i8," ***")') afft(ip+1),i
             endif
             afft(ip+1) = 0.d0
          else if (abs(afft(ip+1)) > 0.d0) then
             afft(ip+1) = 0.d0
          end if
       end do
#endif
    end if

    if(iprinegativecharge0 >= 1) then
       if(npes_in == npes .and. npes > 1) then
          call mpi_allreduce(icwarn_NEGA,icwarn_NEGA_total,1,mpi_integer,mpi_sum,MPI_CommGroup,ierr)
          call mpi_allreduce(icwarn_IMAG,icwarn_IMAG_total,1,mpi_integer,mpi_sum,MPI_CommGroup,ierr)
       else if(npes_in == npes_cdfft .and. npes_cdfft > 1) then
          call mpi_allreduce(icwarn_NEGA,icwarn_NEGA_total,1,mpi_integer,mpi_sum &
               & ,mpi_cdfft_world(myrank_ggacmp),ierr)
          call mpi_allreduce(icwarn_IMAG,icwarn_IMAG_total,1,mpi_integer,mpi_sum &
               & ,mpi_cdfft_world(myrank_ggacmp),ierr)
       else
          icwarn_NEGA_total = icwarn_NEGA
          icwarn_IMAG_total = icwarn_IMAG
       end if
    end if

    if(iprinegativecharge >=1) then
       if(icwarn_NEGA_total >= 1) then
! =============================== modified by K. Tagami =============== 11.0
!          if(nspin == 2) write(nfout,'(" #spin = ",i3," : 1 = UP, 2 = DOWN")') ispin
!
          if ( nspin == 100 ) then
             write(nfout,'(" #loop = ",i3)' ) ispin
          else if (nspin == 2) then
             write(nfout,'(" #spin = ",i3," : 1 = UP, 2 = DOWN")') ispin
          endif
! ===================================================================== 11.0

          if(npes == 1) then
             !write(nfout,'(" *** WARN  # of <<Negative Charge Density>>  = ",i9)') icwarn_NEGA_total
             write(nfout,'(" ### Warning(4202): Number of <<Negative Charge Density>> =",i9)') icwarn_NEGA_total

          else
             !write(nfout,'(" *** WARN  # of <<Negative Charge Density>>  = ",i9,", (node ",i3,") = ",i9)') &
             !     & icwarn_NEGA_total, mype,icwarn_NEGA
             write(nfout,'(" ### Warning(4202): Number of <<Negative Charge Density>> =",i9,", (node ",i3,") = ",i9)') &
                  & icwarn_NEGA_total, mype,icwarn_NEGA
          end if
       endif

       if(icwarn_IMAG_total >= 1) then
! =============================== modified by K. Tagami =============== 11.0
!          if(nspin == 2) write(nfout,'(" #spin = ",i3," : 1 = UP, 2 = DOWN")') ispin
!
          if ( nspin == 100 ) then
             write(nfout,'(" #loop = ",i3)' ) ispin
          else if (nspin == 2) then
             write(nfout,'(" #spin = ",i3," : 1 = UP, 2 = DOWN")') ispin
          endif
! ===================================================================== 11.0

          if(npes == 1) then
             write(nfout,'(" *** WARN  # of <<Imaginary Charge Density>> = ",i9)') icwarn_IMAG_total
          else
             write(nfout,'(" *** WARN  # of <<Imaginary Charge Density>> = ",i9,", (node ",i3,") = ",i9)') &
                  & icwarn_IMAG_total, mype,icwarn_IMAG
          end if
       end if
    end if
  end subroutine m_FFT_check_of_negative_CD

  subroutine m_FFT_W_Vlocal_W(electron_or_positron,nabfft,afft,bfft,eg)
    integer, intent(in) ::                             electron_or_positron
    integer, intent(in) ::                             nabfft
    real(kind=DP), intent(in),    dimension(nabfft) :: afft
    real(kind=DP), intent(inout), dimension(nabfft) :: bfft
    real(kind=DP), intent(out)                 :: eg

    integer i, id, j, k, ip, nl, nm, nn, nd2, nd3, nlh, idh
    real(kind=DP) :: s



    if(electron_or_positron == ELECTRON) then
       id  = fft_box_size_WF(1,0)
       nl  = fft_box_size_WF(1,1)
       nm  = fft_box_size_WF(2,1)
       nn  = fft_box_size_WF(3,1)
       nd2 = fft_box_size_WF(2,0)
       nd3 = fft_box_size_WF(3,0)
    else if(electron_or_positron == POSITRON) then
       id  = fft_box_size_pWF(1,0)
       nl  = fft_box_size_pWF(1,1)
       nm  = fft_box_size_pWF(2,1)
       nn  = fft_box_size_pWF(3,1)
       nd2 = fft_box_size_pWF(2,0)
       nd3 = fft_box_size_pWF(3,0)
    end if

    s = 0.d0
    if(kimg == 1) then
       do i = 1, id*nd2*nn-1, 2
          bfft(i) = afft(i)*(bfft(i)**2 + bfft(i+1)**2)
       end do
       nlh = nl/2
       do i = 1, id*nd2*nn-1, 2
          s = s + bfft(i)
       end do
       s = s + s
       do i = 1, id*nd2*nn-1, id
          s = s - bfft(i)
       end do
       do i = nlh*2+1, id*nd2*nn-1, id
          s = s - bfft(i)
       end do
!!$       write(6,'(" s = ", d20.8)') s
    else if(kimg == 2) then
       do i = 1, id*nd2*nn*kimg-1, 2
          s = s + afft(i)*(bfft(i)**2 + bfft(i+1)**2)
       end do
!!$       if(sw_zero_padding == ON) then
!!$          do j = nl+1, id
!!$             do i = 1, nd2*nn
!!$                ip = id*(i-1) + j
!!$                bfft(ip*2-1) = 0.d0
!!$                bfft(ip*2) = 0.d0
!!$             end do
!!$          end do
!!$          do j = nm+1, nd2
!!$             do k = 1, nn
!!$                do i = 1, nl
!!$                   ip = i + id*(j-1) + id*nd2*(k-1)
!!$                   bfft(ip*2-1) = 0.d0
!!$                   bfft(ip*2) = 0.d0
!!$                end do
!!$             end do
!!$          end do
    end if
    eg = s / product(fft_box_size_WF(1:3,1))
  end subroutine m_FFT_W_Vlocal_W

! =========================== added by K. Tagami =================== 11.0
  subroutine m_FFT_W_Vlocal_W_noncl(electron_or_positron,nabfft,afft_kt,bfft_kt,eg, &
       &     ndim_spinor, ndim_chgpot )
    integer, intent(in) ::     electron_or_positron
    integer, intent(in) ::     nabfft
    integer, intent(in) ::     ndim_chgpot, ndim_spinor

    real(kind=DP), intent(in)    :: afft_kt( nabfft,ndim_chgpot )
    real(kind=DP), intent(inout) :: bfft_kt( nabfft,ndim_spinor  )
    real(kind=DP), intent(out)                 :: eg
!
    real(kind=DP), allocatable ::  afft_tmp(:)
    integer i, id, j, k, ip, nl, nm, nn, nd2, nd3, nlh, idh
    real(kind=DP) :: s, ctmp_r, ctmp_i

!
    integer :: is1, is2, is_tmp
    complex(kind=CMPLDP) :: z1, z2, z3, zsum

    if(electron_or_positron == ELECTRON) then
       id  = fft_box_size_WF(1,0)
       nl  = fft_box_size_WF(1,1)
       nm  = fft_box_size_WF(2,1)
       nn  = fft_box_size_WF(3,1)
       nd2 = fft_box_size_WF(2,0)
       nd3 = fft_box_size_WF(3,0)
    else if(electron_or_positron == POSITRON) then
       id  = fft_box_size_pWF(1,0)
       nl  = fft_box_size_pWF(1,1)
       nm  = fft_box_size_pWF(2,1)
       nn  = fft_box_size_pWF(3,1)
       nd2 = fft_box_size_pWF(2,0)
       nd3 = fft_box_size_pWF(3,0)
    end if

    allocate( afft_tmp( nabfft ) ); afft_tmp = 0.0d0

    s = 0.d0
    if(kimg == 1) then
       do i = 1, id*nd2*nn-1, 2
          ctmp_r = bfft_kt( i,  1 ) *bfft_kt( i,2 ) + bfft_kt( i+1,1 )*bfft_kt( i+1,2 )
          ctmp_i = bfft_kt( i+1,1 ) *bfft_kt( i,2 ) - bfft_kt( i  ,1 )*bfft_kt( i+1,2 )

          afft_tmp(i) =  afft_kt(i,1)*( bfft_kt(i,1)**2 + bfft_kt(i+1,1)**2 ) &
               &        +afft_kt(i,4)*( bfft_kt(i,2)**2 + bfft_kt(i+1,2)**2 ) &
               &        +( afft_kt(i,2)*ctmp_r -afft_kt(i+1,2)*ctmp_i )*2.0d0

       end do

       nlh = nl/2
       do i = 1, id*nd2*nn-1, 2
          s = s + afft_tmp(i)
       end do
       s = s + s
       do i = 1, id*nd2*nn-1, id
          s = s - afft_tmp(i)
       end do
       do i = nlh*2+1, id*nd2*nn-1, id
          s = s - afft_tmp(i)
       end do
!!$       write(6,'(" s = ", d20.8)') s
    else if(kimg == 2) then

       do i = 1, id*nd2*nn*kimg -1, 2
          zsum = 0.0d0
          Do is1=1, ndim_spinor
             Do is2=1, ndim_spinor

                is_tmp = ( is1- 1 )*ndim_spinor + is2
                z1 = dcmplx( bfft_kt(i,is1), bfft_kt(i+1,is1) )
                z2 = dcmplx( bfft_kt(i,is2), bfft_kt(i+1,is2) )
                z3 = dcmplx( afft_kt(i,is_tmp), afft_kt(i+1,is_tmp) )
                zsum = zsum + conjg( z1 )*z2 *z3
             End do
          End do

          afft_tmp(i) = real(zsum)

!          ctmp_r = bfft_kt( i,  1 ) *bfft_kt( i,2 ) + bfft_kt( i+1,1 )*bfft_kt( i+1,2 )
!          ctmp_i = bfft_kt( i+1,1 ) *bfft_kt( i,2 ) - bfft_kt( i  ,1 )*bfft_kt( i+1,2 )

!          afft_tmp(i) =  afft_kt(i,1)*( bfft_kt(i,1)**2 + bfft_kt(i+1,1)**2 ) &
!               &        +afft_kt(i,4)*( bfft_kt(i,2)**2 + bfft_kt(i+1,2)**2 ) &
!               &        +( afft_kt(i,2)*ctmp_r -afft_kt(i+1,2)*ctmp_i )*2.0d0

       end do
       do i = 1, id*nd2*nn*kimg -1, 2
          s = s + afft_tmp(i)
       end do

    end if
    eg = s / product(fft_box_size_WF(1:3,1))

    if ( allocated( afft_tmp ) ) deallocate( afft_tmp )

  end subroutine m_FFT_W_Vlocal_W_noncl
!! ============================================================================ 11.0

  subroutine m_FFT_afft_allgatherv(afft_mpi0,afft)
    real(kind=DP),intent(in),dimension(:) :: afft_mpi0(ista_sfftp:iend_sfftp)
    real(kind=DP),intent(out),dimension(:):: afft(nfftp)
    integer  :: id_sname = -1
    if(npes >= 2) then
       call tstatc0_begin('m_FFT_afft_allgatherv(in m_FFT.) ',id_sname)
       call mpi_allgatherv(afft_mpi0,nel_sfftp(mype),mpi_double_precision &  ! MPI
            &       ,afft,nel_sfftp,idisp_sfftp,mpi_double_precision,MPI_CommGroup,ierr)
    else     ! npes == 1
       afft = afft_mpi0
    end if
    if(npes >= 2) call tstatc0_end(id_sname)
  end subroutine m_FFT_afft_allgatherv

  subroutine m_FFT_afft_allgatherv_c(afft_mpi0,afft)
    real(kind=DP),intent(in),dimension(:) :: afft_mpi0(ista_fftp:iend_fftp)
    real(kind=DP),intent(out),dimension(:):: afft(nfftp)
    integer  :: id_sname = -1
    if(npes_cdfft >= 2) then
       call tstatc0_begin('m_FFT_afft_allgatherv_c(in m_FFT.) ',id_sname)
       if(myrank_ggacmp <nrank_ggacmp) then
          call mpi_allgatherv(afft_mpi0,nel_fftp(myrank_cdfft),mpi_double_precision &  ! MPI
               &       ,afft,nel_fftp,idisp_fftp,mpi_double_precision,mpi_cdfft_world(myrank_ggacmp),ierr)
       end if
       call tstatc0_end(id_sname)
    else     ! npes == 1
       afft = afft_mpi0
    end if
  end subroutine m_FFT_afft_allgatherv_c

  subroutine m_FFT_Vlocal_pW(afft,bfft)
    real(kind=DP), intent(in),    dimension(nfft_pstrn) :: afft
    real(kind=DP), intent(inout), dimension(nfft_pstrn) :: bfft

    integer i
    do i = 1, nfft_pstrn-1, 2
       bfft(i)   = afft(i) * bfft(i)
       bfft(i+1) = afft(i) * bfft(i+1)
    end do
    if(mod(nfft_pstrn,2) == 1) bfft(nfft_pstrn) = afft(nfft_pstrn)*bfft(nfft_pstrn)
  end subroutine m_FFT_Vlocal_pW

! === EXP_CELLOPT ==== 2015/09/24
  subroutine m_FFT_coef_CD_integration_kt(ista,iend,f2or1)
                                            ! from m_XC_Potential.F90,  for electron ?
    integer, intent(in) :: ista, iend
    real(kind=DP),intent(out) :: f2or1( ista:iend )

    integer                  :: idp,nlp,nmp,nnp,nd2p,nd3p,ip, idph, nlph

    nlp  = fft_box_size_CD(1,1)
    nmp  = fft_box_size_CD(2,1)
    nnp  = fft_box_size_CD(3,1)
#ifdef _MPIFFT_
    idp  = fft_box_size_CD_c(1,0)
    nd2p = fft_box_size_CD_c(2,0)
    nd3p = fft_box_size_CD_c(3,0)
#else
    idp  = fft_box_size_CD(1,0)
    nd2p = fft_box_size_CD(2,0)
    nd3p = fft_box_size_CD(3,0)
#endif

    call set_f2or1( npes, ista, iend, f2or1 )

  contains

    subroutine set_f2or1(npes,ista,iend,f2or1)
      integer, intent(in) :: npes,ista, iend
      real(kind=DP), intent(out), dimension(ista:iend) :: f2or1
      integer :: idph,nlph,ip,i,j,k

      if(kimg == 1) then
         idph = idp/2
         nlph = nlp/2
#ifdef _MPIFFT_

         f2or1 = 0.d0
!!$         do j = 1, min(nz_d,nnp-nz_d*myrank_cdfft)*nmp
!!$         do j = 1, min(nz_d,nnp-nz_d*myrank_cdfft)*nd2p
         do k = 1, min(nz_d, nnp-nz_d*myrank_cdfft)
            do j = 1, nmp
               do i = 1, nlph
!!$                  ip = i + idph*(j-1) + idph*ly*lz_d*myrank_cdfft
                  ip = i + idph*(j-1) + idph*ly*(k-1) + idph*ly*lz_d*myrank_cdfft
                  f2or1(ip) = 2.d0
               end do
               ip = 1 + idph*(j-1) + idph*ly*(k-1) + idph*ly*lz_d*myrank_cdfft
               f2or1(ip) = 1.d0
               ip = nlph+1 + idph*(j-1) + idph*ly*(k-1) + idph*ly*lz_d*myrank_cdfft
               f2or1(ip) = 1.d0
            end do
         end do
!!$            ip = idph*(j-1) + 1 + idph*ly*lz_d*myrank_cdfft
!!$            f2or1(ip) = 1.d0
!!$            ip = idph*(j-1)+ nlph + 1 + idph*ly*lz_d*myrank_cdfft
!!$            f2or1(ip) = 1.d0
!!$         end do
#else
         f2or1 = 2.d0
         if(npes >= 2) then
            do i = 1, nd2p*nnp
               ip = idph*(i-1) + 1
               if(ip>= ista .and. ip <= iend) f2or1(ip) = 1.d0
            end do
            do i = 1, nd2p*nnp
               ip = idph*(i-1) + nlph + 1
               if(ip>= ista .and. ip <= iend) f2or1(ip) = 1.d0
            end do
            do j = nlph+2, idph
               do i = 1,nd2p*nnp
                  ip = idph*(i-1)+j
                  if(ip>= ista .and. ip <= iend) f2or1(ip) = 0.d0
               end do
            end do
            do j = nmp+1, nd2p
               do k = 1, nnp
                  do i = 1, nlph
                     ip = i + idph*(j-1) + idph*nd2p*(k-1)
                     if(ip >= ista .and. ip <= iend) f2or1(ip) = 0.d0
                  end do
               end do
            end do
! for SX6 ASL 20040817
            do k = nnp+1, nd3p
               do i = 1, idph*nd2p
                  ip = i + idph*nd2p*(k-1)
                  if(ip >= ista .and. ip <= iend) f2or1(ip) = 0.d0
               end do
            end do
         else
            do i = 1, nd2p*nnp
               ip = idph*(i-1) + 1
               f2or1(ip) = 1.d0
            end do
            do i = 1, nd2p*nnp
               ip = idph*(i-1) + nlph + 1
               f2or1(ip) = 1.d0
            end do
            do j = nlph+2, idph
               do i = 1,nd2p*nnp
                  ip = idph*(i-1)+j
                  f2or1(ip) = 0.d0
               end do
            end do
            do j = nmp+1, nd2p
               do k = 1, nnp
                  do i = 1, nlph
                     ip = i + idph*(j-1) + idph*nd2p*(k-1)
                     f2or1(ip) = 0.d0
                  end do
               end do
            end do
! for SX6 ASL 20040817
            do k = nnp+1, nd3p
               do i = 1, idph*nd2p
                  ip = i + idph*nd2p*(k-1)
                  f2or1(ip) = 0.d0
               end do
            end do
         end if
!!$       do i = ista, iend
!!$          if(mod(i*2,idp) == 2 .or. mod(i*2,idp) == 0) f2or1(i) = 1.d0
!!$       end do
#endif
      else
#ifdef _MPIFFT_
!!$         f2or1 = 0.d0                 ! f2or1 works to the fft data in Rspace.
!!$         do k = 1, min(nz_d,nnp-nz_d*myrank_cdfft)
!!$            do j = 1, nmp     ! nmp = fft_box_size_CD(2,1)
!!$               do i = 1, nlp  ! nlp = fft_box_size_CD(1,1)
!!$                  ip = i+(j-1)*idp+(k-1)*idp*nd2p+idp*nd2p*lz_d*myrank_cdfft
!!$                  f2or1(ip) = 1.d0
!!$               end do
!!$            end do
!!$         end do
!         if(iprixc >= 2 ) write(nfout,'(" ix kimg = 2 <<set_f2or1>>")')
         f2or1 = 1.d0
         do j = nlp+1, idp      ! x
            do i = 1, ly*nz_d
               ip = idp*(i-1)+j+ista-1
               f2or1(ip) = 0.d0
            end do
         end do
!         if(iprixc >= 2 ) write(nfout,'(" iy kimg = 2 <<set_f2or1>>")')
         do j = nmp+1,ly         ! y
            do k = 1, nz_d
               do i = 1, nlp
                  ip = i + idp*(j-1) + idp*ly*(k-1) + ista-1
                  f2or1(ip) = 0.d0
               end do
            end do
         end do
!         if(iprixc >= 2 ) write(nfout,'(" iz kimg = 2 <<set_f2or1>>")')
         do  k = nz_d+1, lz_d   ! z
            do i = 1, idp*ly
               ip = i + idp*ly*(k-1) + ista-1
               f2or1(ip) = 0.d0
            end do
         end do
#else
         f2or1 = 1.d0
         if(npes >= 2) then
            do j = nlp+1, idp    ! x
               do i = 1, nd2p*nnp
                  ip = idp*(i-1)+j
                  if(ip>= ista .and. ip <= iend) f2or1(ip) = 0.d0
               end do
            end do
            do j = nmp+1, nd2p   ! y
               do k = 1, nnp
                  do i = 1, nlp
                     ip = i + idp*(j-1) + idp*nd2p*(k-1)
                     if(ip >= ista .and. ip <= iend) f2or1(ip) = 0.d0
                  end do
               end do
            end do
! for SX6 ASL 20040817
            do k = nnp+1, nd3p   ! z
               do i = 1, idp*nd2p
                  ip = i + idp*nd2p*(k-1)
                  if(ip >= ista .and. ip <= iend) f2or1(ip) = 0.d0
               end do
            end do
         else
            do j = nlp+1, idp    ! x
               do i = 1, nd2p*nnp
                  ip = idp*(i-1)+j
! ================================ modifed by K. Tagami ====( uncertain )== 11.0
!                  f2or1(ip) = 0.d0
                  if(ip >= ista .and. ip <= iend) f2or1(ip) = 0.d0
! ===================================================================== 11.0
               end do
            end do
            do j = nmp+1, nd2p   ! y
               do k = 1, nnp
                  do i = 1, nlp
                     ip = i + idp*(j-1) + idp*nd2p*(k-1)
! ================================ modifed by K. Tagami ====( uncertain )== 11.0
!                     f2or1(ip) = 0.d0
                     if(ip >= ista .and. ip <= iend) f2or1(ip) = 0.d0
! ===================================================================== 11.0
                  end do
               end do
            end do
! for SX6 ASL 20040817
            do k = nnp+1, nd3p   ! z
               do i = 1, idp*nd2p
                  ip = i + idp*nd2p*(k-1)
! ================================ modifed by K. Tagami ====( uncertain )== 11.0
!                  f2or1(ip) = 0.d0
                  if(ip >= ista .and. ip <= iend) f2or1(ip) = 0.d0
! ===================================================================== 11.0
               end do
            end do
         end if
#endif
      end if
    end subroutine set_f2or1
  end subroutine m_FFT_coef_CD_integration_kt
! ==================== 2015/09/24

  subroutine m_FFT_coef_CD_integration(f2or1)
    real(kind=DP),intent(out), dimension(ista_sfftph:iend_sfftph) :: f2or1
    integer :: idp, nlp, nmp, nnp, nd2p, idph, nlph, ip, i, j, k

    idp  = fft_box_size_CD(1,0)
    nlp  = fft_box_size_CD(1,1)
    nmp  = fft_box_size_CD(2,1)
    nnp  = fft_box_size_CD(3,1)
    nd2p = fft_box_size_CD(2,0)


    if(kimg == 1) then
       f2or1 = 2.d0
       idph = idp/2
       nlph = nlp/2
       if(npes >= 2) then
          do i = 1, nd2p*nnp
             ip = idph*(i-1) + 1
             if(ip>= ista_sfftph .and. ip <= iend_sfftph) f2or1(ip) = 1.d0
          end do
          do i = 1, nd2p*nnp
             ip = idph*(i-1) + nlph + 1
             if(ip>= ista_sfftph .and. ip <= iend_sfftph) f2or1(ip) = 1.d0
          end do
          do j = nlph+2, idph
             do i = 1,nd2p*nnp
                ip = idph*(i-1)+j
                if(ip>= ista_sfftph .and. ip <= iend_sfftph) f2or1(ip) = 0.d0
             end do
          end do
          do j = nmp+1, nd2p
             do i = 1, nlph
                ip = i + idph*(j-1) + idph*nd2p*(k-1)
                if(ip >= ista_sfftph .and. ip <= iend_sfftph) f2or1(ip) = 0.d0
             end do
          end do
       else
          do i = 1, nd2p*nnp
             ip = idph*(i-1) + 1
             f2or1(ip) = 1.d0
          end do
          do i = 1, nd2p*nnp
             ip = idph*(i-1) + nlph + 1
             f2or1(ip) = 1.d0
          end do
          do j = nlph+2, idph
             do i = 1,nd2p*nnp
                ip = idph*(i-1)+j
                f2or1(ip) = 0.d0
             end do
          end do
          do j = nmp+1, nd2p
             do i = 1, nlph
                ip = i + idph*(j-1) + idph*nd2p*(k-1)
                f2or1(ip) = 0.d0
             end do
          end do
       end if
!!$       do i = ista_sfftph, iend_sfftph
!!$          if(mod(i*2,idp) == 2 .or. mod(i*2,idp) == 0) f2or1(i) = 1.d0
!!$       end do
    else
       f2or1 = 1.d0
       if(npes >= 2) then
          do j = nlp+1, idp
             do i = 1, nd2p*nnp
                ip = idp*(i-1)+j
                if(ip>= ista_sfftph .and. ip <= iend_sfftph) f2or1(ip) = 0.d0
             end do
          end do
          do j = nmp+1, nd2p
             do k = 1, nnp
                do i = 1, nlp
                   ip = i + idp*(j-1) + idp*nd2p*(k-1)
                   if(ip >= ista_sfftph .and. ip <= iend_sfftph) f2or1(ip) = 0.d0
                end do
             end do
          end do
       else
          do j = nlp+1, idp
             do i = 1, nd2p*nnp
                ip = idp*(i-1)+j
                f2or1(ip) = 0.d0
             end do
          end do
          do j = nmp+1, nd2p
             do k = 1, nnp
                do i = 1, nlp
                   ip = i + idp*(j-1) + idp*nd2p*(k-1)
                   f2or1(ip) = 0.d0
                end do
             end do
          end do
       end if
    end if
  end subroutine m_FFT_coef_CD_integration
#ifdef _MPIFFT_
  subroutine m_FFT_set_cdata(aa,a)
    real(kind=DP), dimension(nfftp), intent(in)                :: aa
    real(kind=DP), dimension(ista_fftp:iend_fftp), intent(out) :: a
    integer :: i,j,k,ii,jj,kk,nn(3),nstart(3),nend(3)
    integer :: lx, llx,lly,lxy,llxy, indexa, ipout, ipin
    indexa = 2
    if(kimg == 1) then
       a = 0.d0
       nn(1:3) = fft_box_size_CD(1:3,1)
       nn(2)   = ny_d
       nstart(1:3) = 1
       nend(1:3) = nn(1:3)
       nstart(indexa)=myrank_cdfft*nn(indexa)+1
       nend(indexa)  =min(fft_box_size_CD(indexa,1),(myrank_cdfft+1)*nn(indexa))
       lx = fft_box_size_CD_c(1,0)
!!$       ly = ly_d   !          =fft_box_size_CD_c(2,0)/npes_cdfft
       llx = fft_box_size_CD_c(1,0)
       lly = fft_box_size_CD_c(2,0)
       lxy = lx*ly_d
       llxy = llx*lly
       do i=nstart(1),nend(1)
          ii=i-nstart(1)+1
          do j=nstart(2),nend(2)
             jj=j-nstart(2)+1
             do k=nstart(3),nend(3)
                kk=k-nstart(3)+1
                ipout = ista_fftp-1+ (ii +(jj-1)* lx+(kk-1)* lxy)
                ipin  =              ( i +( j-1)*llx+( k-1)*llxy)
                A(ipout)   = AA(ipin)
!!$                A(ipout+1) = AA(ipin+1)
             enddo
          enddo
       enddo
    else if(kimg == 2) then
       A = 0.d0
       nn(1:3) = fft_box_size_CD(1:3,1)
       nn(2) = ny_d
       nstart(1:3) = 1
       nend(1:3) = nn(1:3)
       nstart(indexa)=myrank_cdfft*nn(indexa)+1
       nend(indexa)  =min(fft_box_size_CD(indexa,1),(myrank_cdfft+1)*nn(indexa))
       lx = fft_box_size_CD_c(1,0)
!!$       ly = ly_d   !          =fft_box_size_CD_c(2,0)/npes_cdfft
       llx = fft_box_size_CD_c(1,0)
       lly = fft_box_size_CD_c(2,0)
       lxy = lx*ly_d
       llxy = llx*lly
       if(ipri >= 2) then
          write(nfout,'(" nstart(1), nend(1) = ",2i8)') nstart(1), nend(1)
          write(nfout,'(" nstart(2), nend(2) = ",2i8)') nstart(2), nend(2)
          write(nfout,'(" nstart(3), nend(3) = ",2i8)') nstart(3), nend(3)
          write(nfout,'(" llx,lx,lly,ly_d    = ",4i8)') llx,lx,lly,ly_d
       end if
!!$       if(npes_cdfft == 1) nend(2) = nend(2)/2
!!$       if(myrank_cdfft == 0) then
       do i=nstart(1),nend(1)
          ii=i-nstart(1)+1
          do j=nstart(2),nend(2)
             jj=j-nstart(2)+1
             do k=nstart(3),nend(3)
                kk=k-nstart(3)+1
                ipout = ista_fftp-1 + 2*(ii +(jj-1)* lx+(kk-1)* lxy)-1
                ipin  =               2*( i +( j-1)*llx+( k-1)*llxy)-1
                A(ipout)   = AA(ipin)
                A(ipout+1) = AA(ipin+1)
             end do
          end do
       end do
    end if

  end subroutine m_FFT_set_cdata

  subroutine m_FFT_set_scdata(aa,a)
    real(kind=DP), dimension(*), intent(in)                    :: aa
    real(kind=DP), dimension(ista_fftp:iend_fftp), intent(out) :: a
    integer :: i,j,k,ii,jj,kk,nn(3),nstart(3),nend(3)
    integer :: lx,llx,lly,lxy,llxy,indexa
    indexa = 2
    if(kimg == 2) then
       nn(1:3) = fft_box_size_CD(1:3,1)
       nn(2) = ny_ds
       nstart(1:3) = 1
       nend(1:3) = nn(1:3)
       nstart(indexa)=mype*nn(indexa)+1
       nend(indexa)  =min(fft_box_size_CD(indexa,1),(mype+1)*nn(indexa))
       lx = fft_box_size_CD_c(1,0)
!!$       ly = fft_box_size_CD(2,0)/npes ! = ly_ds
       llx = fft_box_size_CD_c(1,0)
       lly = fft_box_size_CD_c(2,0)
       lxy = lx*ly_ds
       llxy = llx*lly
       do i=nstart(1),nend(1)
          ii=i-nstart(1)+1
          do j=nstart(2),nend(2)
             jj=j-nstart(2)+1
             do k=nstart(3),nend(3)
                kk=k-nstart(3)+1
                A(ista_fftp-1+ii+(jj-1)*lx+(kk-1)*lxy)=AA(i+(j-1)*llx+(k-1)*llxy-ista_fftp)
             enddo
          enddo
       enddo
    end if

  end subroutine m_FFT_set_scdata
#endif

#ifdef _MPIFFT_
  subroutine check_ierr(ierr0)
    integer, intent(in),dimension(2) :: ierr0
    IF( IERR0(1) .GE. 3000 ) THEN
       WRITE(6,*) '*** 3D FFT ERROR ***'
       WRITE(6,*) mype,':IERR =',IERR0(1)
       IF( IERR0(1) .EQ. -1 ) THEN
          WRITE(6,*) mype,':MPI ERROR =',IERR0(2)
       ENDIF
       call phase_error_with_msg(nfout,'error hyb fft',__LINE__,__FILE__)
    ENDIF
  end subroutine check_ierr

  subroutine wd_time(time)
    real(kind=DP), intent(in), dimension(5) :: time
    WRITE(6,6000)
    WRITE(6,6010) mype,TIME(1:5)
6000 FORMAT(' RANK     X-FFT       Y-FFT       Z-FFT       ', &
          &       'TRANS       TOTAL  (Direct)')
6010 FORMAT(1X,I4,':',6(1PD12.4),0PF9.1)
  end subroutine wd_time

  subroutine get_nxnynz(nx,ny,nz)
    integer, intent(out) :: nx,ny,nz
    NX = fft_box_size_CD(1,1)
    NY = fft_box_size_CD(2,1)
    NZ = fft_box_size_CD(3,1)
  end subroutine get_nxnynz

#endif

  subroutine m_FFT_CD_direct(nfout,afft)   ! R space --> G space
#ifdef NEC_TUNE_FFT
      use m_Parallelization,   only : itask
#endif
    integer, intent(in)       :: nfout
    real(kind=DP), intent(inout) :: afft(ista_sfftp:iend_sfftp)
#ifdef _MPIFFT_
    INTEGER :: NX,NY,NZ
    INTEGER :: NT,ISW,IERR0(2)
    integer, parameter :: INDEXA = 2, INDEXB = 3
    real(kind=DP),allocatable :: B(:)
    real(kind=DP) :: TIME(5)
    integer :: id_sname = -1
    call tstatc0_begin('m_FFT_CD_direct ',id_sname)

    call get_nxnynz(nx,ny,nz)

    allocate(B(ista_sfftp:iend_sfftp))
    ISW = 1
#ifdef NEC_TUNE_FFT
    NT = itask
#else
    NT = 1
#endif

    if(kimg==1) then
       call fft_r3d(nx,ny,nz,ny_ds,nz_ds,afft &
            &      ,lxs,lys,lzs,ly_ds,lz_ds,indexa,isw,time,b &
            &      ,nt, MPI_CommGroup,ierr0)
    else
       call fft_c3d(nx,ny,nz,ny_ds,nz_ds,afft &
            &      ,lxs,lys,lzs,ly_ds,lz_ds,indexa,isw,time,b &
            &      ,nt, MPI_CommGroup,ierr0)
    end if
    call check_ierr(ierr0)

    CALL MPI_BARRIER(MPI_CommGroup,IERR)

    if(ipri >= 2) call wd_time(time)

    deallocate(B)
#else
    integer :: id_sname = -1
    call tstatc0_begin('m_FFT_CD_direct ',id_sname)

    if(npes >= 2) then
       call m_FFT_afft_allgatherv(afft,afft_CD)
    else
       afft_CD = afft
    end if

    call fft_CD_direct_core(afft_CD)

    afft(ista_sfftp:iend_sfftp) = afft_CD(ista_sfftp:iend_sfftp)
#endif
    call tstatc0_end(id_sname)
  end subroutine m_FFT_CD_direct

  subroutine m_FFT_reset_CD_setup_stat()
    CD_setup_is_done = NO
    fft_box_size_CD_prev = fft_box_size_CD
    fft_box_size_CD_nonpara_prev = fft_box_size_CD_nonpara
  end subroutine m_FFT_reset_CD_setup_stat

  subroutine m_FFT_CD0(nfout,afftp,inverse_or_direct)  ! R space --> G space
    integer, intent(in)          :: nfout, inverse_or_direct
    real(kind=DP), intent(inout) :: afftp(nfftp_nonpara)
    integer :: id_sname = -1
    integer, dimension(2) :: flag_mklfft = (/-1, 1/)

    if(ipri >= 2) write(nfout,'(" <<m_FFT_CD0>>")')
    call tstatc0_begin('m_FFT_CD0 ', id_sname, 1)
    if(CD_setup_is_done == NO) then
       if(ipri >= 1) write(nfout,'(" <<CDFFT_setup>>")')
       call CDFFT_setup()
    end if

    if(inverse_or_direct == DIRECT) then
       call fft_CD_direct_core(afftp)
    else if(inverse_or_direct == INVERSE) then
       call fft_CD_inverse_core(afftp)
    end if

    call tstatc0_end(id_sname)
  end subroutine m_FFT_CD0

  subroutine m_FFT_CD0_exx(nfout,afftp,inverse_or_direct)  ! R space --> G space
    integer, intent(in)          :: nfout, inverse_or_direct
    real(kind=DP), intent(inout) :: afftp(nfftp_exx_nonpara)
    integer :: id_sname = -1
    integer, dimension(2) :: flag_mklfft = (/-1, 1/)
#if defined(JRCATFFT_WS) || defined(FFTW3)
    if(ipri >= 2) write(nfout,'(" <<m_FFT_CD0_exx>>")')
    call tstatc0_begin('m_FFT_CD0_exx ', id_sname)
    if(CD_setup_is_done_exx == NO) then
#ifdef _INCLUDE_EXX_
       if(ipri >= 1) write(nfout,'(" <<CDFFT_setup_exx>>")')
       call CDFFT_setup_exx()
#endif
    end if
    if(inverse_or_direct == DIRECT) then
       call fft_CD_direct_core_exx(afftp)
    else if(inverse_or_direct == INVERSE) then
       call fft_CD_inverse_core_exx(afftp)
    end if

    call tstatc0_end(id_sname)
#else
    call phase_error_with_msg(nfout,'EXX in rspace can only be used with JRCAT FFT or FFTW3',__LINE__,__FILE__)
#endif
  end subroutine m_FFT_CD0_exx

  subroutine m_FFT_CD_direct_c(nfout,afft,mode)   ! R space --> G space
#ifdef NEC_TUNE_FFT
      use m_Parallelization,   only : itask
#endif
    integer, intent(in)       :: nfout
    real(kind=DP), intent(inout) :: afft(ista_fftp:iend_fftp)
    integer, intent(in), optional :: mode
#ifdef _MPIFFT_

    INTEGER :: NX,NY,NZ
    INTEGER :: NT,ISW,IERR0(2)

    integer, parameter :: INDEXA = 2, INDEXB = 3
    real(kind=DP),allocatable :: B(:)
    real(kind=DP) :: TIME(5)
#endif

    integer :: id_sname = -1
    call tstatc0_begin('m_FFT_CD_direct_c ',id_sname)
    if(present(mode)) then
      if(mode == 0) goto 9999
    endif

#ifdef _MPIFFT_
    call get_nxnynz(nx,ny,nz)

    allocate(B(ista_fftp:iend_fftp))
    if(kimg == 1) then
       isw = -1
    else if(kimg == 2) then
       ISW = 1
    end if

#ifdef NEC_TUNE_FFT
    NT = itask
#else
    NT = 1
#endif

    if(kimg==1) then
       call fft_r3d(nx,ny,nz,ny_d,nz_d,afft &
            &      ,lx,ly,lz,ly_d,lz_d,indexb,isw,time,b &
            &      ,nt,mpi_cdfft_world(myrank_ggacmp),ierr0)
    else
       call fft_c3d(nx,ny,nz,ny_d,nz_d,afft &
            &      ,lx,ly,lz,ly_d,lz_d,indexb,isw,time,b &
            &      ,nt,mpi_cdfft_world(myrank_ggacmp),ierr0)
    end if
    IF( IERR0(1) .GE. 3000 ) THEN
       write(6,*) ' nx,ny,nz,ny_d,nz_d = ',nx,ny,nz,ny_d,nz_d
       write(6,*) ' lx,ly,lz,ly_d,lz_d = ',lx,ly,lz,ly_d,lz_d
    ENDIF
    call check_ierr(ierr0)

    CALL MPI_BARRIER(mpi_cdfft_world(myrank_ggacmp),IERR)
    if(ipri >= 2) call wd_time(time)

    deallocate(B)

#else
    call m_FFT_afft_allgatherv_c(afft,afft_CD)

    call fft_CD_direct_core(afft_CD)

    afft(ista_fftp:iend_fftp) = afft_CD(ista_fftp:iend_fftp)
#endif
9999 call tstatc0_end(id_sname)
  end subroutine m_FFT_CD_direct_c

  subroutine m_FFT_CD_inverse(nfout,afft)  ! G space --> R space
#ifdef NEC_TUNE_FFT
      use m_Parallelization,   only : itask
#endif
    integer, intent(in)          :: nfout
    real(kind=DP), intent(inout) :: afft(ista_sfftp:iend_sfftp)
#ifdef _MPIFFT_
    INTEGER :: NX,NY,NZ
    INTEGER :: NT,ISW,IERR0(2)
    integer, parameter :: INDEXA = 2, INDEXB = 3
    real(kind=DP),allocatable :: B(:)
    real(kind=DP) :: TIME(5)

    integer :: id_sname = -1
    call tstatc0_begin('m_FFT_CD_inverse ', id_sname)

    call get_nxnynz(nx,ny,nz)

    allocate(B(ista_sfftp:iend_sfftp))
    CALL MPI_BARRIER(MPI_CommGroup,IERR)
    ISW = -1

#ifdef NEC_TUNE_FFT
    NT = itask
#else
    NT = 1
#endif

    if(kimg==1) then
       call fft_r3d(nx,ny,nz,ny_ds,nz_ds,afft &
            &      ,lxs,lys,lzs,ly_ds,lz_ds,indexa,isw,time,b &
            &      ,nt, MPI_CommGroup,ierr0)
    else
       call fft_c3d(nx,ny,nz,ny_ds,nz_ds,afft &
            &      ,lxs,lys,lzs,ly_ds,lz_ds,indexa,isw,time,b &
            &      ,nt,MPI_CommGroup,ierr0)
    end if
    call check_ierr(ierr0)
    CALL MPI_BARRIER(MPI_CommGroup,IERR)
    if(ipri >= 2) call wd_time(time)
    deallocate(B)
#else
    integer :: id_sname = -1
    call tstatc0_begin('m_FFT_CD_inverse ', id_sname)

    if(npes >= 2) then
       call m_FFT_afft_allgatherv(afft,afft_CD)
    else
       afft_CD = afft
    end if

    call fft_CD_inverse_core(afft_CD)

    afft(ista_sfftp:iend_sfftp) = afft_CD(ista_sfftp:iend_sfftp)
#endif
    call tstatc0_end(id_sname)
  end subroutine m_FFT_CD_inverse

  subroutine m_FFT_CD_inverse_c(nfout,afft,mode)  ! G space --> R space
#ifdef NEC_TUNE_FFT
      use m_Parallelization,   only : itask
#endif
    integer, intent(in)          :: nfout
    real(kind=DP), intent(inout) :: afft(ista_fftp:iend_fftp)
    integer, intent(in), optional :: mode

#ifdef _MPIFFT_
    INTEGER :: NX,NY,NZ,lx,ly,lz
    INTEGER :: NT,ISW,IERR0(2)
    integer, parameter :: INDEXA = 2, INDEXB = 3
    real(kind=DP),allocatable :: B(:)
    integer :: i, ip, j, k
    real(kind=DP) :: TIME(5)
#endif

    integer :: id_sname = -1
    call tstatc0_begin('m_FFT_CD_inverse_c ', id_sname)
    if(present(mode)) then
      if(mode == 0) goto 9999
    endif

#ifdef _MPIFFT_
    call get_nxnynz(nx,ny,nz)

    lx = fft_box_size_CD_c(1,0)
    ly = fft_box_size_CD_c(2,0)
    lz = fft_box_size_CD_c(3,0)

    allocate(B(ista_fftp:iend_fftp))

    if(kimg == 1) then
       isw = 1
    else if(kimg == 2) then
       ISW = -1
    end if

#ifdef NEC_TUNE_FFT
    NT = itask
#else
    NT = 1
#endif

    if(kimg==1) then
       if(ipri >= 2) then
          write(nfout,'(" before fft_r3d <<m_FFT_CD_inverse_c>>")')
          write(nfout,'(" nx,ny,nz,ny_d,nz_d = ",5i5)') nx,ny,nz,ny_d,nz_d
          write(nfout,'(" lx,ly,lz,ly_d,lz_d = ",5i5)') lx,ly,lz,ly_d,lz_d
       end if
       call fft_r3d(nx,ny,nz,ny_d,nz_d,afft &
            &      ,lx,ly,lz,ly_d,lz_d,indexa,isw,time,b &
            &      ,nt,mpi_cdfft_world(myrank_ggacmp),ierr0)
    else if(kimg == 2) then
       if(ipri >= 2) then
          write(nfout,'(" before fft_c3d <<m_FFT_CD_inverse_c>>")')
          do k = 1, nz
             if(k >= 4 .and. k<= nz-3) cycle
             do j = 1, ny_d
                write(nfout,'(" j, k = ", 2i8)') j+ny_d*myrank_cdfft, k
                ip = 2*((j-1)*lx + (k-1)*lx*ly_d) + ista_fftp - 1
                write(nfout,'(8d12.4)') (afft( ip+i),i=1,fft_box_size_CD(1,1)*2)
             end do
          end do
       end if
       call mpi_barrier(mpi_cdfft_world(myrank_ggacmp),ierr)
       call fft_c3d(nx,ny,nz,ny_d,nz_d,afft &
            &      ,lx,ly,lz,ly_d,lz_d,indexa,isw,time,B &
            &      ,nt,mpi_cdfft_world(myrank_ggacmp),ierr0)
       if(ipri >= 2) then
          write(nfout,'(" after fft_c3d <<m_FFT_CD_inverse_c>>")')
          do k = 1, min(nz_d,8)
             do j = 1, ny
                write(nfout,'(" j, k = ", 2i8)') j, k+nz_d*myrank_cdfft
                ip = 2*((j-1)*lx + (k-1)*lx*ly) + ista_fftp - 1
                write(nfout,'(8d12.4)') (afft( ip+i),i=1,16)
             end do
          end do
       end if
    end if
    call check_ierr(ierr0)
    CALL MPI_BARRIER(MPI_cdfft_world(myrank_ggacmp),IERR)
    if(ipri >= 2) call wd_time(time)

    deallocate(B)
#else
    call m_FFT_afft_allgatherv_c(afft,afft_CD)

    call fft_CD_inverse_core(afft_CD)

    afft(ista_fftp:iend_fftp) = afft_CD(ista_fftp:iend_fftp)
#endif
9999 call tstatc0_end(id_sname)
  end subroutine m_FFT_CD_inverse_c

  subroutine m_FFT_CD_inverse0(nfout,afftp)  ! G space --> R space
    integer, intent(in)          :: nfout
    real(kind=DP), intent(inout) :: afftp(nfftp_nonpara)
    integer :: id_sname = -1
    call tstatc0_begin('m_FFT_CD_inverse0 ', id_sname)
    if(CD_setup_is_done == NO) call CDFFT_setup()

    call fft_CD_inverse_core(afftp)

    call tstatc0_end(id_sname)
  end subroutine m_FFT_CD_inverse0

  subroutine m_FFT_CD_direct0(nfout,afftp)  ! G space --> R space
    integer, intent(in)          :: nfout
    real(kind=DP), intent(inout) :: afftp(nfftp_nonpara)
    integer :: id_sname = -1
    call tstatc0_begin('m_FFT_CD_direct0 ', id_sname)
    if(CD_setup_is_done == NO) call CDFFT_setup()

    call fft_CD_direct_core(afftp)

    call tstatc0_end(id_sname)
  end subroutine m_FFT_CD_direct0

!!$#ifdef _MPIFFT_
  subroutine set_mpifft_box_size_CD(inversion_symmetry)
    integer, intent(in) :: inversion_symmetry
    integer :: nx, ny, nz, ic, isize
    nx = fft_box_size_CD(1,1); ny = fft_box_size_CD(2,1); nz = fft_box_size_CD(3,1)
    do ic = 1, 2
       if(ic == 1) isize = npes
       if(ic == 2) isize = npes_cdfft
       ny_d =  ny/isize
       nz_d =  nz/isize
       if(mod( ny,isize) .ne. 0) ny_d = ny_d + 1
       if(mod( nz,isize) .ne. 0) nz_d = nz_d + 1
       ly_d = ny_d;   lz_d = nz_d
       if(inversion_symmetry == ON) then  ! kimg == 1
!!$          if(mod(nx,4) == 2) then
!!$             lx = nx + 4
!!$          else
          lx = nx + 2
!!$          end if
       else if(inversion_symmetry == OFF) then ! kimg == 2
          lx = nx
          if(mod(lx,2) == 0)   lx = nx + ldx
       end if

       if(mod(ly_d,2) == 0) ly_d = ny_d + ldy
       if(mod(lz_d,2) == 0) lz_d = nz_d + ldz

       ly   = ly_d*isize
       lz   = lz_d*isize
       if(ic == 1) then
          fft_box_size_CD(1,0) = lx
          fft_box_size_CD(2,0) = ly
          fft_box_size_CD(3,0) = lz
          ny_ds = ny_d;          nz_ds = nz_d;          ly_ds = ly_d
          lys   = ly;            lz_ds = lz_d;          lzs   = lz
          lxs   = lx
       else
          fft_box_size_CD_c(1,0) = lx
          fft_box_size_CD_c(2,0) = ly
          fft_box_size_CD_c(3,0) = lz
       end if
    end do
    if(ipri >= 1) then
       write(nfout,'(" nx,  ny,  nz,  ny_d, nz_d  = ",5i8)') nx,ny,nz,ny_d,nz_d
       write(nfout,'("      ly,  lz,  ly_d, lz_d  = ",8x,4i8)') ly,lz,ly_d,lz_d
!!$       write(nfout,'(" ny,nz,ny_d, nz_d, ly, lz, ly_d, lz_d  = ",/8x,8i8)')  &
!!$            &          ny,nz,ny_d, nz_d, ly, lz, ly_d, lz_d
       write(nfout,'("                ny_ds,nz_ds = ",24x,2i8)') ny_ds, nz_ds
       write(nfout,'("      lys, lzs, ly_ds,lz_ds = ",8x,4i8)') lys,lzs,ly_ds,lz_ds
!!$       write(nfout,'("       ny_ds,nz_ds,lys,lzs,ly_ds,lz_ds = ",/24x,6i8)')&
!!$            &                ny_ds,nz_ds,lys,lzs,ly_ds,lz_ds
    end if
  end subroutine set_mpifft_box_size_CD
!!$#endif


  subroutine m_FFT_cp_afft_CD(afft)
    real(kind=DP),intent(out),dimension(:):: afft(nfftp)
    !!$afft(1:nfftp)=afft_CD(1:nfftp)
    if(allocated(afft_CD)) afft=afft_CD
  end subroutine m_FFT_cp_afft_CD


end module m_FFT
