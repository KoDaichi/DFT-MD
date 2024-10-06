! === DEBUG by tkato 2013/09/25 ================================================
#define FFTW_STRIDE
! ==============================================================================
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
#ifdef MPI_FFTW
  use, intrinsic :: iso_c_binding
#endif
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
  use m_Parallelization,   only : mpi_kg_world   &
       &                        , mpi_ke_world   &
       &                        , myrank_g, nrank_g   &
       &                        , myrank_e, nrank_e, np_e
  use m_Parallelization,   only : map_fft_x, map_fft_y, map_fft_z &
       &                        , mp_fft_x, mp_fft_y, mp_fft_z    &
       &                        , nel_fft_x, nel_fft_y, nel_fft_z &
       &                        , np_fft_x, np_fft_y, np_fft_z    &
       &                        , xyz_fft_x, xyz_fft_y, xyz_fft_z &
       &                        , map_fftcd_x, map_fftcd_y, map_fftcd_z &
       &                        , mp_fftcd_x, mp_fftcd_y, mp_fftcd_z    &
       &                        , nel_fftcd_x, nel_fftcd_y, nel_fftcd_z &
       &                        , np_fftcd_x, np_fftcd_y                &
       &                        , xyz_fftcd_x, xyz_fftcd_y, xyz_fftcd_z


  use m_Parallelization,   only : fft_X_z_dim, fft_X_y_dim &
       &                        , fft_Y_x_dim, fft_Y_z_dim &
       &                        , fft_Z_x_dim, fft_Z_y_dim &
       &                        , nis_fft_X_z, nis_fft_X_y, nie_fft_X_z, nie_fft_X_y &
       &                        , nis_fft_Y_x, nis_fft_Y_z, nie_fft_Y_x, nie_fft_Y_z &
       &                        , nis_fft_Z_x, nis_fft_Z_y, nie_fft_Z_x, nie_fft_Z_y &
!fj.2012s
       &                        , fft_X_x_dim, fft_X_x_nel, fft_X_y_nel, fft_X_z_nel &
!fj.2012e
       &                        , fftcd_X_z_dim, fftcd_X_y_dim &
       &                        , fftcd_Y_x_dim, fftcd_Y_z_dim &
       &                        , fftcd_Z_x_dim, fftcd_Z_y_dim &
!fj.2012s
       &                        , fftcd_X_x_dim, fftcd_X_x_nel, fftcd_X_y_nel, fftcd_X_z_nel &
!fj.2012e
       &                        , nis_fftcd_X_z, nis_fftcd_X_y, nie_fftcd_X_z, nie_fftcd_X_y &
       &                        , nis_fftcd_Y_x, nis_fftcd_Y_z, nie_fftcd_Y_x, nie_fftcd_Y_z &
       &                        , nis_fftcd_Z_x, nis_fftcd_Z_y, nie_fftcd_Z_x, nie_fftcd_Z_y

  use m_Parallelization,   only : mpi_fft_xz_world,   myrank_fft_xz,   nrank_fft_xz     &
       &                        , mpi_fft_zy_world,   myrank_fft_zy,   nrank_fft_zy     &
!fj.2012s
       &                        , mpi_fft_xy_world,   myrank_fft_xy,   nrank_fft_xy     &
       &                        , mpi_fft_yz_world,   myrank_fft_yz,   nrank_fft_yz     &
       &                        , mpi_fftcd_xy_world, myrank_fftcd_xy, nrank_fftcd_xy     &
       &                        , mpi_fftcd_yz_world, myrank_fftcd_yz, nrank_fftcd_yz     &
!fj.2012e
       &                        , mpi_fftcd_xz_world, myrank_fftcd_xz, nrank_fftcd_xz     &
       &                        , mpi_fftcd_zy_world, myrank_fftcd_zy, nrank_fftcd_zy
  use m_IterationNumbers,   only : iteration

! ================================= added by K. Tagami ============== 11.0
  use m_Const_Parameters,       only : CMPLDP
! =================================================================== 11.0

! ======= KT_add ========================================== 13.0F
  use m_Control_Parameters, only : sw_hybrid_functional
! ========================================================= 13.0F

#ifdef KMATH_FFT3D
  use m_Control_Parameters, only : nstage_fft3d,  sw_kmath_fft3d
  use m_Parallelization, only : kmath3d_handle_wf_inverse, kmath3d_handle_wf_direct, &
  &   nproc_fft3d_inverse, nproc_fft3d_direct, kmath3d_box_size
#endif
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
  integer, dimension(3,0:1) :: fft_box_size_CD_3D

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

  logical :: firstcall_direct_3d = .true.,firstcall_inverse_3d = .true.,firstcall_cd_direct_3d = .true. &
   &          ,firstcall_cd_inverse_3d = .true.,firstcall_direct_xyz_3d = .true.,firstcall_inverse_xyz_3d = .true. &
   &          ,firstcall_cd_direct_xyz_3d = .true.,firstcall_cd_inverse_xyz_3d = .true. &
   &          ,firstcall_direct_3ddiv = .true.,firstcall_inverse_3ddiv = .true.,firstcall_cd_direct_3ddiv = .true. &
   &          ,firstcall_cd_inverse_3ddiv = .true.,firstcall_inverse_xyz_3d_oop = .true.

! contains
#ifdef JRCATFFT_WS
include "m_FFT_type1_jrcatfft.F90"
#define _INCLUDE_EXX_
#elif FFTW3 && !MPI_FFTW
include "m_FFT_type4_fftw3_rev.F90"
#define _INCLUDE_EXX_
#elif MPI_FFTW
include "m_FFT_type10_mpi_fftw3.F90"
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

  subroutine m_FFT_reset_firstcall()
    firstcall_direct_3d = .true.;firstcall_inverse_3d = .true.;firstcall_cd_direct_3d = .true.
    firstcall_cd_inverse_3d = .true.;firstcall_direct_xyz_3d = .true.;firstcall_inverse_xyz_3d = .true.
    firstcall_cd_direct_xyz_3d = .true.;firstcall_cd_inverse_xyz_3d = .true.
    firstcall_direct_3ddiv = .true.;firstcall_inverse_3ddiv = .true.;firstcall_cd_direct_3ddiv = .true.
    firstcall_cd_inverse_3ddiv = .true.
#ifdef MPI_FFTW
    distfftw_inverse_wf_flg = .false.
    distfftw_direct_wf_flg = .false.
#endif
  end subroutine m_FFT_reset_firstcall

  subroutine wd_FFTboxsizes(nfout)
    integer, intent(in) :: nfout

#ifndef FFT_USE_SSL2
    write(nfout,'(" !FFT(D) FFT_USE_SSL2 is not defined ")')
#ifdef FFT_FFTW_OLD
    write(nfout,'(" !FFT(D) FFT_FFTW_OLD is defined ")')
#else
    write(nfout,'(" !FFT(D) FFT_FFTW_OLD is not defined ")')
#endif
#else
    write(nfout,'(" !FFT(D) FFT_USE_SSL2 is defined ")')
#ifdef FFT_FFTW_OLD
    write(nfout,'(" !FFT(D) FFT_FFTW_OLD is defined ")')
#else
    write(nfout,'(" !FFT(D) FFT_FFTW_OLD is not defined ")')
#endif
#endif

    write(nfout,*) ' ---- (WF)FFT size ----'
    write(nfout,'("     (nfft, fft_box_size_WF)")')
    write(nfout,100) nfft
    write(nfout,101) fft_box_size_WF(1:3,0)
    write(nfout,102) fft_box_size_WF(1:3,1)


    write(nfout,*) ' ---- (CD)FFT size ----'
    write(nfout,'("     (nfftps, fft_box_size_CD_3D)")')
    write(nfout,100) nfftps
    write(nfout,101) fft_box_size_CD_3D(1:3,0)
    write(nfout,102) fft_box_size_CD_3D(1:3,1)

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
       call phase_error_with_msg(nfout, 'error hyb fft',__LINE__,__FILE__)
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
       if(ipri >= 2) write(nfout,'(" <<CDFFT_setup>>")')
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
       if(ipri >= 2) write(nfout,'(" <<CDFFT_setup_exx>>")')
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
    if(ipri >= 2) then
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

#ifdef FFTW3
! --------------
!#ifdef PARA3D
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  subroutine m_FFT_WF_3D(electron_or_positron,nfout,afft_l,lsize,iesize,inverse_or_direct)
    integer, intent(in) :: electron_or_positron, nfout, lsize, iesize, inverse_or_direct
    real(kind=DP), dimension(lsize*kimg,iesize), intent(inout) :: afft_l
    integer               :: nsize, ierr

#ifdef __TIMER_SUB__
    call timer_sta(104)
#endif
    if(inverse_or_direct == DIRECT) then
       call m_FFT_Direct_3D (nfout, afft_l, lsize, iesize)
    else ! INVERSE
       call m_FFT_Inverse_3D(nfout, afft_l, lsize, iesize)
    end if
#ifdef __TIMER_SUB__
    call timer_end(104)
#endif

  end subroutine m_FFT_WF_3D

!------------------------------------------------------------------------------

#ifndef FFT_USE_SSL2
#ifdef FFT_FFTW_OLD
!------------------------------------------------------------------------------
  subroutine m_FFT_Direct_3D(nfout,wk_afft_l, wk_size, iesize)
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
use mod_timer
#endif
! === TIMERTIMERTIMER ==========================================================
    integer, intent(in) :: nfout, wk_size, iesize
    real(kind=DP), dimension(wk_size*kimg,iesize), intent(inout) :: wk_afft_l

    real(kind=DP), allocatable, dimension(:) :: wk_allfft
    integer,      allocatable, dimension(:)   :: wk_recvdsp
    integer, parameter :: FFTW_MEASURE=0
    integer, parameter :: FFTW_ESTIMATE=64
    integer, parameter :: FFTW_RANK=1
    integer            :: FFTW_FLAG
    integer,dimension(FFTW_RANK) :: NFFTW3, NEMBED, NEREAL
    integer :: nx, ny, nz, i, j, k, l, ri, ib, iaddr, lrank
    integer :: irank, itag, icnt_send, icnt_recv
    integer,save, allocatable, dimension(:)   :: req_r, req_s
    integer,save, allocatable, dimension(:,:) :: sta_r, sta_s
    integer,save, allocatable, dimension(:)   :: wk_recvcnt, wk_sendcnt
!   integer,save, allocatable, dimension(:)   :: wk_recvrank, wk_sendrank
!F  real(kind=DP),allocatable, dimension(:,:) :: wk_recv, wk_send
    real(kind=DP),allocatable, dimension(:,:),save :: wk_recv1, wk_send1
    real(kind=DP),allocatable, dimension(:,:),save :: wk_recv2, wk_send2
    integer,      allocatable, dimension(:)   :: wk_mp_fft_y
    real(kind=DP),allocatable, dimension(:,:),save :: wk_afft_y
    real(kind=DP),allocatable, dimension(:,:,:)   :: wk_gather
    integer :: mpi_comm, myrank, nmrank
    integer(kind=8),save :: planx1 = 0, planx2 = 0
    integer(kind=8),save :: plany1 = 0, plany2 = 0
    integer(kind=8),save :: planz1 = 0, planz2 = 0
    integer(kind=8) :: plan
    integer :: max_x, max_y, max_z
    integer,save :: max_elm = 0, savesize = 0

    integer,save, allocatable, dimension(:,:) :: y2z_recv, y2z_send, z2x_recv, z2x_send
    integer,save :: y2z_rrank, y2z_srank, z2x_rrank, z2x_srank
    integer,save :: y2z_rmax, y2z_smax, y2z_srmax, z2x_rmax, z2x_smax, z2x_srmax

    integer ::  nnx, nny, nnz, iadd0, iadd1, nis, nie, jadd, kadd
#ifdef FFT_ALLTOALL
    integer,dimension(2) :: max_send,max_recv
#endif

#ifdef __TIMER_SUB__
    call mpi_barrier(mpi_ke_world,ierr)
    call timer_sta(105)
#endif

    if (kimg == 1) then
       FFTW_FLAG = +1      ! FFTW_BACKWARD
    else
       FFTW_FLAG = -1      ! FFTW_FORWARD
    endif
    plan = 0
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
call start_timer('FFT-F reg0')
#endif
! === TIMERTIMERTIMER ==========================================================
!!  mpi_comm = mpi_kg_world
!!  myrank = myrank_e
!!  nmrank = nrank_e
    mpi_comm = mpi_ke_world
    myrank = myrank_g
    nmrank = nrank_g
    itag = 10

    if (np_fft_x==0)then
       nx = 0
       nx = 0
       nz = 0
       go to 1000
    endif

    if (firstcall_direct_3d) then
       savesize = 0
#ifdef __TIMER_DO__
  call timer_sta(118)
#endif
       max_x = maxval(nel_fft_x(:))
       max_y = maxval(nel_fft_y(:))
       max_z = maxval(nel_fft_z(:))
       max_elm = max(max_x,max_y,max_z)

       if(allocated(req_r)) deallocate(req_r)
       if(allocated(req_s)) deallocate(req_s)
       if(allocated(sta_r)) deallocate(sta_r)
       if(allocated(sta_s)) deallocate(sta_s)
       if(allocated(wk_recvcnt)) deallocate(wk_recvcnt)
       if(allocated(wk_sendcnt)) deallocate(wk_sendcnt)
       allocate(req_r(0:nmrank-1), stat=ierr)
       allocate(req_s(0:nmrank-1), stat=ierr)
       allocate(sta_r(MPI_STATUS_SIZE,0:nmrank-1), stat=ierr)
       allocate(sta_s(MPI_STATUS_SIZE,0:nmrank-1), stat=ierr)
       allocate(wk_recvcnt(0:nmrank-1), stat=ierr)
       allocate(wk_sendcnt(0:nmrank-1), stat=ierr)
!      allocate(wk_recvrank(0:nmrank-1), stat=ierr)
!      allocate(wk_sendrank(0:nmrank-1), stat=ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_FFT_Direct_3D :  Not allocate '
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 51, ierr)
        endif

       if(allocated(y2z_recv)) deallocate(y2z_recv)
       if(allocated(y2z_send)) deallocate(y2z_send)
       allocate(y2z_recv(2,0:nmrank-1))
       allocate(y2z_send(2,0:nmrank-1))
       y2z_send = 0
       y2z_recv = 0
       wk_recvcnt = 0
       wk_sendcnt = 0
       do i = 1, nel_fft_z(myrank)
          irank = map_fft_y(mp_fft_z(i)) - 1
          wk_recvcnt(irank) = wk_recvcnt(irank) + 1
       enddo
       do i = 1, nel_fft_y(myrank)
          irank = map_fft_z(mp_fft_y(i)) - 1
          wk_sendcnt(irank) = wk_sendcnt(irank) + 1
       enddo
       k = 0
       do i = 0, nmrank - 1
          if(wk_recvcnt(i) /= 0) then
             k = k + 1
             y2z_recv(1,i) = wk_recvcnt(i)
             y2z_recv(2,i) = k
          endif
       enddo
       y2z_rrank = k
       k = 0
       do i = 0, nmrank - 1
          if(wk_sendcnt(i) /= 0) then
             k = k + 1
             y2z_send(1,i) = wk_sendcnt(i)
             y2z_send(2,i) = k
          endif
       enddo
       y2z_srank = k
       y2z_rmax = maxval(wk_recvcnt)
       y2z_smax = maxval(wk_sendcnt)
       y2z_srmax = max(y2z_rmax,y2z_smax)

       if(allocated(z2x_recv)) deallocate(z2x_recv)
       if(allocated(z2x_send)) deallocate(z2x_send)
       allocate(z2x_recv(2,0:nmrank-1))
       allocate(z2x_send(2,0:nmrank-1))
       z2x_send = 0
       z2x_recv = 0
       wk_recvcnt = 0
       wk_sendcnt = 0
       do i = 1, nel_fft_x(myrank)
          irank = map_fft_z(mp_fft_x(i)) - 1
          wk_recvcnt(irank) = wk_recvcnt(irank) + 1
       enddo
       do i = 1, nel_fft_z(myrank)
          irank = map_fft_x(mp_fft_z(i)) - 1
          wk_sendcnt(irank) = wk_sendcnt(irank) + 1
       enddo
       k = 0
       do i = 0, nmrank - 1
          if(wk_recvcnt(i) /= 0) then
             k = k + 1
             z2x_recv(1,i) = wk_recvcnt(i)
             z2x_recv(2,i) = k
          endif
       enddo
       z2x_rrank = k
       k = 0
      do i = 0, nmrank - 1
          if(wk_sendcnt(i) /= 0) then
             k = k + 1
             z2x_send(1,i) = wk_sendcnt(i)
             z2x_send(2,i) = k
          endif
       enddo
       z2x_srank = k
       z2x_rmax = maxval(wk_recvcnt)
       z2x_smax = maxval(wk_sendcnt)
       z2x_srmax = max(z2x_rmax,z2x_smax)

#ifdef FFT_ALLTOALL
       max_send(1) = y2z_rmax
       max_send(2) = y2z_smax
       call mpi_allreduce(max_send,max_recv,2,mpi_integer,mpi_max,mpi_fft_zy_world,ierr)
       y2z_rmax = max_recv(1)
       y2z_smax = max_recv(2)
       max_send(1) = z2x_rmax
       max_send(2) = z2x_smax
       call mpi_allreduce(max_send,max_recv,2,mpi_integer,mpi_max,mpi_fft_xz_world,ierr)
       z2x_rmax = max_recv(1)
       z2x_smax = max_recv(2)
#endif

       if (ipri > 1) then
         write(nfout,'("m_FFT_Direct_3D   --   myrank_g=",i4)') myrank_g
         write(nfout,'("y2z_send")')
         write(nfout,'(10(i8,", "))') (y2z_send(1,i),i=0,nmrank-1)
         write(nfout,'(10(i8,", "))') (y2z_send(2,i),i=0,nmrank-1)
         write(nfout,'("y2z_recv")')
         write(nfout,'(10(i8,", "))') (y2z_recv(1,i),i=0,nmrank-1)
         write(nfout,'(10(i8,", "))') (y2z_recv(2,i),i=0,nmrank-1)
         write(nfout,'("z2x_send")')
         write(nfout,'(10(i8,", "))') (z2x_send(1,i),i=0,nmrank-1)
         write(nfout,'(10(i8,", "))') (z2x_send(2,i),i=0,nmrank-1)
         write(nfout,'("z2x_recv")')
         write(nfout,'(10(i8,", "))') (z2x_recv(1,i),i=0,nmrank-1)
         write(nfout,'(10(i8,", "))') (z2x_recv(2,i),i=0,nmrank-1)
       endif

#ifdef __TIMER_DO__
  call timer_sta(119)
#endif
       nx = xyz_fft_y(2,1)-xyz_fft_y(1,1)+1
       ny = xyz_fft_y(2,2)-xyz_fft_y(1,2)+1
       nz = xyz_fft_y(2,3)-xyz_fft_y(1,3)+1
!      if(kimg==1) then
          NFFTW3(1) = ny
          NEMBED(1) = ny
          call dfftw_plan_many_dft    (plany1, FFTW_RANK, NFFTW3, nz*nx/2, &
         &                             wk_afft_l, NEMBED, nz*nx/2, 1,    &
         &                             wk_afft_l, NEMBED, nz*nx/2, 1,    &
         &                             FFTW_FLAG, FFTW_MEASURE )
!      else
          NFFTW3(1) = ny
          NEMBED(1) = ny
          call dfftw_plan_many_dft    (plany2, FFTW_RANK, NFFTW3, nz*nx,   &
         &                             wk_afft_l, NEMBED, nz*nx, 1,      &
         &                             wk_afft_l, NEMBED, nz*nx, 1,      &
         &                             FFTW_FLAG, FFTW_MEASURE )
!      end if

       nx = xyz_fft_z(2,1)-xyz_fft_z(1,1)+1
       ny = xyz_fft_z(2,2)-xyz_fft_z(1,2)+1
       nz = xyz_fft_z(2,3)-xyz_fft_z(1,3)+1
!      if(kimg==1) then
          NFFTW3(1) = nz
          NEMBED(1) = nz
          call dfftw_plan_many_dft    (planz1, FFTW_RANK, NFFTW3, nx*ny/2, &
         &                             wk_afft_l, NEMBED, nx*ny/2, 1,     &
         &                             wk_afft_l, NEMBED, nx*ny/2, 1,     &
         &                             FFTW_FLAG, FFTW_MEASURE )
!      else
          NFFTW3(1) = nz
          NEMBED(1) = nz
          call dfftw_plan_many_dft    (planz2, FFTW_RANK, NFFTW3, nx*ny,   &
         &                             wk_afft_l, NEMBED, nx*ny, 1,       &
         &                             wk_afft_l, NEMBED, nx*ny, 1,       &
         &                             FFTW_FLAG, FFTW_MEASURE )
!      endif

       nx = xyz_fft_x(2,1)-xyz_fft_x(1,1)+1
       ny = xyz_fft_x(2,2)-xyz_fft_x(1,2)+1
       nz = xyz_fft_x(2,3)-xyz_fft_x(1,3)+1
!      if(kimg==1) then
          NFFTW3(1) = nx - 2
          NEMBED(1) = nx
          NEREAL(1) = nx
          call dfftw_plan_many_dft_c2r(planx1, FFTW_RANK, NFFTW3, ny*nz,   &
         &                             wk_afft_l, NEMBED, 1, nx/2,        &
         &                             wk_afft_l, NEREAL, 1, nx,          &
         &                             FFTW_MEASURE )
!      else
          NFFTW3(1) = nx
          NEMBED(1) = nx
          call dfftw_plan_many_dft    (planx2, FFTW_RANK, NFFTW3, ny*nz,   &
         &                             wk_afft_l, NEMBED, 1, nx,          &
         &                             wk_afft_l, NEMBED, 1, nx,          &
         &                             FFTW_FLAG, FFTW_MEASURE )
!      endif
#ifdef __TIMER_DO__
  call timer_end(119)
#endif

       firstcall_direct_3d = .false.
       go to 1000
#ifdef __TIMER_DO__
  call timer_end(118)
#endif
    endif

#ifdef FFT_ALLTOALL
    if (iesize /= savesize) then
#ifdef __TIMER_DO__
  call timer_sta(120)
#endif
       if (allocated(wk_recv1)) deallocate(wk_recv1)
       if (allocated(wk_send1)) deallocate(wk_send1)
       if (allocated(wk_recv2)) deallocate(wk_recv2)
       if (allocated(wk_send2)) deallocate(wk_send2)
       if (allocated(wk_afft_y)) deallocate(wk_afft_y)
       allocate(wk_recv1(y2z_rmax*kimg*iesize,y2z_rrank), stat=ierr)
       allocate(wk_send1(y2z_smax*kimg*iesize,y2z_srank), stat=ierr)
       allocate(wk_recv2(z2x_rmax*kimg*iesize,z2x_rrank), stat=ierr)
       allocate(wk_send2(z2x_smax*kimg*iesize,z2x_srank), stat=ierr)
       allocate(wk_afft_y(max_elm*kimg,iesize) ,stat=ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_FFT_Direct_3D :  Not allocate '
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 52, ierr)
        endif
       savesize = iesize
#ifdef __TIMER_DO__
  call timer_end(120)
#endif
    end if
#else
    if (iesize > savesize) then
#ifdef __TIMER_DO__
  call timer_sta(120)
#endif
       if (allocated(wk_recv1)) deallocate(wk_recv1)
       if (allocated(wk_send1)) deallocate(wk_send1)
       if (allocated(wk_recv2)) deallocate(wk_recv2)
       if (allocated(wk_send2)) deallocate(wk_send2)
       if (allocated(wk_afft_y)) deallocate(wk_afft_y)
       allocate(wk_recv1(y2z_rmax*kimg*iesize,y2z_rrank), stat=ierr)
       allocate(wk_send1(y2z_smax*kimg*iesize,y2z_srank), stat=ierr)
       allocate(wk_recv2(z2x_rmax*kimg*iesize,z2x_rrank), stat=ierr)
       allocate(wk_send2(z2x_smax*kimg*iesize,z2x_srank), stat=ierr)
       allocate(wk_afft_y(max_elm*kimg,iesize) ,stat=ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_FFT_Direct_3D :  Not allocate '
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 52, ierr)
        endif
       savesize = iesize
#ifdef __TIMER_DO__
  call timer_end(120)
#endif
    end if
#endif

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP          SHARED(     nfout,wk_afft_l,iesize,req_r,req_s,sta_r,sta_s,    &
!$OMP                 wk_recv1,wk_send1,wk_recv2,wk_send2,wk_afft_y,          &
!$OMP                 planx1,planx2,planz1,planz2,plany1,plany2,              &
!$OMP                 z2x_recv, z2x_send, y2z_recv, y2z_send,itag,            &
!$OMP                 xyz_fft_y,xyz_fft_z,xyz_fft_x,fft_Z_x_dim,nis_fft_Z_x,  &
!$OMP                 nie_fft_Z_x,fft_X_z_dim,nis_fft_X_z,nie_fft_X_z,        &
!$OMP                 fft_Y_z_dim,nis_fft_Y_z,nie_fft_Y_z,fft_Z_y_dim,        &
!$OMP                 nie_fft_Z_y,nis_fft_Z_y,                                &
#ifdef FFT_ALLTOALL
!$OMP                 mpi_fft_xz_world, z2x_rmax, z2x_smax,                   &
!$OMP                 mpi_fft_zy_world, y2z_rmax, y2z_smax,                   &
#endif
!$OMP                 mpi_comm,myrank,nmrank,kimg,ierr                    )   &
!$OMP          PRIVATE(nx,ny,nz,i,j,k,l,ri,ib,iaddr,lrank,icnt_send,icnt_recv, &
!$OMP                  irank,nnx,nny,nnz,iadd0,iadd1,nis,nie,jadd,kadd,plan   )

    nx = xyz_fft_y(2,1)-xyz_fft_y(1,1)+1
    ny = xyz_fft_y(2,2)-xyz_fft_y(1,2)+1
    nz = xyz_fft_y(2,3)-xyz_fft_y(1,3)+1
!
! Y-axis (z-x div)
!
    if(kimg==1) then
       plan = plany1
    else
       plan = plany2
    end if
#ifdef __TIMER_DO__
  call timer_sta(251)
#endif
!$OMP DO
    do ib = 1, iesize
       call dfftw_execute_dft(plan,wk_afft_l(1,ib),wk_afft_l(1,ib))
    end do
#ifdef __TIMER_DO__
  call timer_end(251)
#endif

#ifdef __TIMER_DO__
  call timer_sta(121)
#endif
!fj --------------------
#ifdef FFT_INTERCHANGE
    if (kimg == 1) then
       do ib = 1, iesize
!$OMP DO
          do k = 0, nz-1
             do j = 0, ny-1
                do i = 0, nx-1
                   do ri = 1, kimg
                      wk_afft_y((i+j*nx+k*nx*ny+1)*kimg+(ri-kimg),ib) = wk_afft_l((i+k*nx+j*nx*nz+1)*kimg+(ri-kimg),ib)
                   enddo
                enddo
             enddo
          enddo
       enddo
    else
       do ib = 1, iesize
!$OMP DO
          do k = 0, nz-1
             do j = 0, ny-1
                do i = 0, nx-1
                   wk_afft_y((i+j*nx+k*nx*ny+1)*2-1,ib) = wk_afft_l((i+k*nx+j*nx*nz+1)*2-1,ib)
                   wk_afft_y((i+j*nx+k*nx*ny+1)*2  ,ib) = wk_afft_l((i+k*nx+j*nx*nz+1)*2  ,ib)
                enddo
             enddo
          enddo
       enddo
    endif
#else
    if (kimg == 1) then
!$OMP DO
       do ib = 1, iesize
          do k = 0, nz-1
             do j = 0, ny-1
                do i = 0, nx-1
                   do ri = 1, kimg
                      wk_afft_y((i+j*nx+k*nx*ny+1)*kimg+(ri-kimg),ib) = wk_afft_l((i+k*nx+j*nx*nz+1)*kimg+(ri-kimg),ib)
                   enddo
                enddo
             enddo
          enddo
       enddo
    else
!$OMP DO
       do ib = 1, iesize
          do k = 0, nz-1
             do j = 0, ny-1
                do i = 0, nx-1
                   wk_afft_y((i+j*nx+k*nx*ny+1)*2-1,ib) = wk_afft_l((i+k*nx+j*nx*nz+1)*2-1,ib)
                   wk_afft_y((i+j*nx+k*nx*ny+1)*2  ,ib) = wk_afft_l((i+k*nx+j*nx*nz+1)*2  ,ib)
                enddo
             enddo
          enddo
       enddo
    endif
#endif
!fj --------------------
#ifdef __TIMER_DO__
  call timer_end(121)
#endif

#ifdef __TIMER_DO__
  call timer_sta(122)
#endif

    nx = xyz_fft_y(2,1)-xyz_fft_y(1,1)+1
    ny = xyz_fft_y(2,2)-xyz_fft_y(1,2)+1
    nz = xyz_fft_y(2,3)-xyz_fft_y(1,3)+1
#ifdef FFT_INTERCHANGE
    if (kimg == 1) then
!$OMP DO
       do ib = 1, iesize
          do k = 1, nz
             do l = 1, fft_Z_y_dim
                nis = nis_fft_Z_y(l)
                nie = nie_fft_Z_y(l)
                nny = nie - nis + 1
                iadd0 = nx*nny*nz*(ib-1)
                do j = nis, nie
                   do i = 1, nx
                      wk_send1(iadd0+(i+(j-nis)*nx+(k-1)*nx*nny),l) = wk_afft_y(i+(j-1)*nx+(k-1)*nx*ny,ib)
                   end do
                end do
             end do
          end do
       end do
    else
!$OMP DO
       do ib = 1, iesize
          do k = 1, nz
             do l = 1, fft_Z_y_dim
                nis = nis_fft_Z_y(l)
                nie = nie_fft_Z_y(l)
                nny = nie - nis + 1
                do j = nis, nie
                   do i = 1, nx
                      iadd1 = nx*nny*nz*(ib-1)*2+(i+((j-nis)*nx+(k-1)*nx*nny))*2
                      iaddr = (i+(j-1)*nx+(k-1)*nx*ny)*2
                      wk_send1(iadd1-1,l) = wk_afft_y(iaddr-1,ib)
                      wk_send1(iadd1  ,l) = wk_afft_y(iaddr  ,ib)
                   end do
                end do
             end do
          end do
       end do
    end if
#else
!$OMP DO
    do l = 1, fft_Z_y_dim
       nis = nis_fft_Z_y(l)
       nie = nie_fft_Z_y(l)
       nny = nie - nis + 1
       if (kimg == 1) then
          do ib = 1, iesize
             iadd0 = nx*nny*nz*(ib-1)
             do k = 1, nz
                do j = nis, nie
                   do i = 1, nx
                      wk_send1(iadd0+(i+(j-nis)*nx+(k-1)*nx*nny),l) = wk_afft_y(i+(j-1)*nx+(k-1)*nx*ny,ib)
                   end do
                end do
             end do
          end do
       else
          do ib = 1, iesize
             iadd0 = nx*nny*nz*(ib-1)*2
             do k = 1, nz
                kadd = (k-1)*nx*nny
                do j = nis, nie
                   jadd = (j-nis)*nx+kadd
                   do i = 1, nx
                      iadd1 = iadd0+(i+jadd)*2
                      iaddr = (i+(j-1)*nx+(k-1)*nx*ny)*2
!fj --------------------
!                     do ri = 1, kimg
!                       wk_send1(iadd1+(ri-kimg),l) = wk_afft_y(iaddr+(ri-kimg),ib)
!                     end do
                      wk_send1(iadd1-1,l) = wk_afft_y(iaddr-1,ib)
                      wk_send1(iadd1  ,l) = wk_afft_y(iaddr  ,ib)
!fj --------------------
                   end do
                end do
             end do
          end do
       end if
    end do
#endif

#ifdef __TIMER_DO__
  call timer_end(122)
#endif

#ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_fft_zy_world)
  call timer_sta(123)
#endif

! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call start_timer('*FFT-F y2z a2a')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================
!$OMP MASTER

    call MPI_ALLTOALL(wk_send1, y2z_smax*kimg*iesize, mpi_double_precision,   &
   &                  wk_recv1, y2z_rmax*kimg*iesize, mpi_double_precision,   &
   &                                                  mpi_fft_zy_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Direct_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 55, ierr)
     endif

!$OMP END MASTER
!$OMP BARRIER
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call stop_timer('*FFT-F y2z a2a')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================

#ifdef __TIMER_DO__
  call timer_sta(124)
#endif
#ifdef __TIMER_DO__
  call timer_end(124)
#endif

#else
!else ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_comm)
  call timer_sta(123)
#endif

! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call start_timer('*FFT-F y2z sr')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================
!$OMP MASTER
    icnt_recv = 0
    lrank = myrank + 1
    if (lrank > (nmrank-1)) lrank = 0
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank-1)) lrank = 0
       if ((lrank /= myrank) .and. (y2z_recv(1,lrank) /= 0)) then
          call mpi_irecv(wk_recv1(1,y2z_recv(2,lrank)), y2z_recv(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_r(icnt_recv), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Direct_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 55, ierr)
           endif
          icnt_recv = icnt_recv + 1
       endif
    enddo

    icnt_send = 0
    lrank = myrank
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank-1)) lrank = 0
       if ((lrank /= myrank) .and. (y2z_send(1,lrank) /= 0)) then
          call mpi_isend(wk_send1(1,y2z_send(2,lrank)), y2z_send(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_s(icnt_send), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Direct_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 56, ierr)
           endif
          icnt_send = icnt_send + 1
       endif
    enddo
!$OMP END MASTER
!$OMP BARRIER

#ifdef __TIMER_DO__
  call timer_sta(124)
#endif
!$OMP DO
    do i = 1, y2z_recv(1,myrank)*kimg*iesize
       wk_recv1(i,y2z_recv(2,myrank)) = wk_send1(i,y2z_send(2,myrank))
    enddo
#ifdef __TIMER_DO__
  call timer_end(124)
#endif

!$OMP MASTER
    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Direct_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 57, ierr)
     endif

    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Direct_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 58, ierr)
     endif
!$OMP END MASTER
!$OMP BARRIER
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call stop_timer('*FFT-F y2z sr')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================

#endif
!endif ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_end(123)
#endif

#ifdef __TIMER_DO__
  call timer_sta(125)
#endif
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call start_timer('*FFT-F y2z rbuf')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================
    nx = xyz_fft_z(2,1)-xyz_fft_z(1,1)+1
    ny = xyz_fft_z(2,2)-xyz_fft_z(1,2)+1
    nz = xyz_fft_z(2,3)-xyz_fft_z(1,3)+1
!$OMP DO
    do l = 1, fft_Y_z_dim
       nis = nis_fft_Y_z(l)
       nie = nie_fft_Y_z(l)
       nnz = nie - nis + 1
       if (kimg == 1) then
          do ib = 1, iesize
             iadd0 = nx*ny*nnz*(ib-1)
             do k = nis, nie
                do j = 1, ny
                   do i = 1, nx
                      wk_afft_l(i+(j-1)*nx+(k-1)*nx*ny,ib) = wk_recv1(iadd0+i+(j-1)*nx+(k-nis)*nx*ny,l)
                   end do
                end do
             end do
          end do
       else
          do ib = 1, iesize
             iadd0 = nx*ny*nnz*(ib-1)*2
             do k = nis, nie
                do j = 1, ny
                   do i = 1, nx
                      iadd1 = iadd0+(i+(j-1)*nx+(k-nis)*nx*ny)*2
                      iaddr = (i+(j-1)*nx+(k-1)*nx*ny)*2
                      wk_afft_l(iaddr-1,ib) = wk_recv1(iadd1-1,l)
                      wk_afft_l(iaddr  ,ib) = wk_recv1(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end if
    end do
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call stop_timer('*FFT-F y2z rbuf')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_DO__
  call timer_end(125)
#endif
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call stop_timer('FFT-F y2z')
call start_timer('FFT-F zfft')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================

!
! Z-axis (x-y div)
!
    nx = xyz_fft_z(2,1)-xyz_fft_z(1,1)+1
    ny = xyz_fft_z(2,2)-xyz_fft_z(1,2)+1
    nz = xyz_fft_z(2,3)-xyz_fft_z(1,3)+1
    if(kimg==1) then
       plan = planz1
    else
       plan = planz2
    end if
#ifdef __TIMER_DO__
  call timer_sta(252)
#endif
!$OMP DO
    do ib = 1, iesize
       call dfftw_execute_dft(plan,wk_afft_l(1,ib),wk_afft_l(1,ib))
    end do
#ifdef __TIMER_DO__
  call timer_end(252)
#endif
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call stop_timer('FFT-F zfft')
call start_timer('FFT-F z2x')
call start_timer('*FFT-F z2x sbuf')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================

#ifdef __TIMER_DO__
  call timer_sta(126)
#endif
    nx = xyz_fft_z(2,1)-xyz_fft_z(1,1)+1
    ny = xyz_fft_z(2,2)-xyz_fft_z(1,2)+1
    nz = xyz_fft_z(2,3)-xyz_fft_z(1,3)+1
!$OMP DO
    do l = 1, fft_X_z_dim
       nis = nis_fft_X_z(l)
       nie = nie_fft_X_z(l)
       nnz = nie - nis + 1
       if (kimg == 1) then
          do ib = 1, iesize
             iadd0 = nx*ny*nnz*(ib-1)
             do k = nis, nie
                do j = 1, ny
                   do i = 1, nx
                      wk_send2(iadd0+(i+(j-1)*nx+(k-nis)*nx*ny),l) = wk_afft_l(i+(j-1)*nx+(k-1)*nx*ny,ib)
                   end do
                end do
             end do
          end do
       else
          do ib = 1, iesize
             iadd0 = nx*ny*nnz*(ib-1)*2
             do k = nis, nie
                jadd = (k-nis)*nx*ny
                do j = 1, ny
                   kadd = (j-1)*nx+jadd
                   do i = 1, nx
                      iadd1 = iadd0+(i+kadd)*2
                      iaddr  = (i+(j-1)*nx+(k-1)*nx*ny)*2
!fj --------------------
!                     do ri = 1, kimg
!                        wk_send2(iadd1+(ri-kimg),l) = wk_afft_l(iaddr+(ri-kimg),ib)
!                     end do
                      wk_send2(iadd1-1,l) = wk_afft_l(iaddr-1,ib)
                      wk_send2(iadd1  ,l) = wk_afft_l(iaddr  ,ib)
!fj --------------------
                   end do
                end do
             end do
          end do
       end if
    end do

! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call stop_timer('*FFT-F z2x sbuf')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_DO__
  call timer_end(126)
#endif

#ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_fft_xz_world)
  call timer_sta(127)
#endif

! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call start_timer('*FFT-F z2x a2a')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================
!$OMP MASTER

    call MPI_ALLTOALL(wk_send2, z2x_smax*kimg*iesize, mpi_double_precision,   &
   &                  wk_recv2, z2x_rmax*kimg*iesize, mpi_double_precision,   &
   &                                                  mpi_fft_xz_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Direct_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 60, ierr)
     endif

!$OMP END MASTER
!$OMP BARRIER
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call stop_timer('*FFT-F z2x a2a')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================

#ifdef __TIMER_DO__
  call timer_sta(128)
#endif
#ifdef __TIMER_DO__
  call timer_end(128)
#endif

#else
!else ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_comm)
  call timer_sta(127)
#endif

! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call start_timer('*FFT-F z2x sr')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================
!$OMP MASTER
    icnt_recv = 0
    lrank = myrank + 1
    if (lrank > (nmrank-1)) lrank = 0
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank-1)) lrank = 0
       if ((lrank /= myrank) .and. (z2x_recv(1,lrank) /= 0)) then
          call mpi_irecv(wk_recv2(1,z2x_recv(2,lrank)), z2x_recv(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_r(icnt_recv), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Direct_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 60, ierr)
           endif
          icnt_recv = icnt_recv + 1
       endif
    enddo

    icnt_send = 0
    lrank = myrank
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank-1)) lrank = 0
       if ((lrank /= myrank) .and. (z2x_send(1,lrank) /= 0)) then
          call mpi_isend(wk_send2(1,z2x_send(2,lrank)), z2x_send(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_s(icnt_send), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Direct_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 61, ierr)
           endif
          icnt_send = icnt_send + 1
       endif
    enddo
!$OMP END MASTER
!$OMP BARRIER

#ifdef __TIMER_DO__
  call timer_sta(128)
#endif
!$OMP DO
    do i = 1, z2x_recv(1,myrank)*kimg*iesize
       wk_recv2(i,z2x_recv(2,myrank)) = wk_send2(i,z2x_recv(2,myrank))
    enddo
#ifdef __TIMER_DO__
  call timer_end(128)
#endif

!$OMP MASTER
    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Direct_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 62, ierr)
     endif

    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Direct_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 63, ierr)
     endif
!$OMP END MASTER
!$OMP BARRIER
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call stop_timer('*FFT-F z2x sr')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================
#endif
!endif ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_end(127)
#endif

#ifdef __TIMER_DO__
  call timer_sta(129)
#endif
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call start_timer('*FFT-F z2x rbuf')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================
    nx = xyz_fft_x(2,1)-xyz_fft_x(1,1)+1
    ny = xyz_fft_x(2,2)-xyz_fft_x(1,2)+1
    nz = xyz_fft_x(2,3)-xyz_fft_x(1,3)+1
#ifdef FFT_INTERCHANGE
    if (kimg == 1) then
!$OMP DO
       do ib = 1, iesize
          do k = 1, nz
             do j = 1, ny
                do l = 1, fft_Z_x_dim
                   nis = nis_fft_Z_x(l)
                   nie = nie_fft_Z_x(l)
                   nnx = nie - nis + 1
                   iadd0 = nnx*ny*nz*(ib-1)
                   do i = nis, nie
                      wk_afft_l(i+(j-1)*nx+(k-1)*nx*ny,ib) = wk_recv2(iadd0+i-nis+1+(j-1)*nnx+(k-1)*nnx*ny,l)
                   end do
                end do
             end do
          end do
       end do
    else
!$OMP DO
       do ib = 1, iesize
          do k = 1, nz
             do j = 1, ny
                do l = 1, fft_Z_x_dim
                   nis = nis_fft_Z_x(l)
                   nie = nie_fft_Z_x(l)
                   nnx = nie - nis + 1
                   do i = nis, nie
                      iadd1 = nnx*ny*nz*(ib-1)*2+(i-nis+1+(j-1)*nnx+(k-1)*nnx*ny)*2
                      iaddr = (i+(j-1)*nx+(k-1)*nx*ny)*2
                      wk_afft_l(iaddr-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_l(iaddr  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end do
    end if
#else
!$OMP DO
    do l = 1, fft_Z_x_dim
       nis = nis_fft_Z_x(l)
       nie = nie_fft_Z_x(l)
       nnx = nie - nis + 1
       if (kimg == 1) then
          do ib = 1, iesize
             iadd0 = nnx*ny*nz*(ib-1)
             do k = 1, nz
                do j = 1, ny
                   do i = nis, nie
                      wk_afft_l(i+(j-1)*nx+(k-1)*nx*ny,ib) = wk_recv2(iadd0+i-nis+1+(j-1)*nnx+(k-1)*nnx*ny,l)
                   end do
                end do
             end do
          end do
       else
          do ib = 1, iesize
             iadd0 = nnx*ny*nz*(ib-1)*2
             do k = 1, nz
                do j = 1, ny
                   do i = nis, nie
                      iadd1 = iadd0+(i-nis+1+(j-1)*nnx+(k-1)*nnx*ny)*2
                      iaddr = (i+(j-1)*nx+(k-1)*nx*ny)*2
                      wk_afft_l(iaddr-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_l(iaddr  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end if
    end do
#endif
#ifdef __TIMER_DO__
  call timer_end(129)
#endif
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call stop_timer('*FFT-F z2x rbuf')
call stop_timer('FFT-F z2x')
call start_timer('FFT-F xfft')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================

!
! X-axis (y-z div)
!
    nx = xyz_fft_x(2,1)-xyz_fft_x(1,1)+1
    ny = xyz_fft_x(2,2)-xyz_fft_x(1,2)+1
    nz = xyz_fft_x(2,3)-xyz_fft_x(1,3)+1
    if(kimg==1) then
       plan = planx1
    else
       plan = planx2
    end if
#ifdef __TIMER_DO__
  call timer_sta(253)
#endif
!$OMP DO
    do ib = 1, iesize
       if (kimg==1) then
          call dfftw_execute_dft_c2r(plan,wk_afft_l(1,ib),wk_afft_l(1,ib))
       else
          call dfftw_execute_dft    (plan,wk_afft_l(1,ib),wk_afft_l(1,ib))
       endif
    end do
#ifdef __TIMER_DO__
  call timer_end(253)
#endif
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call stop_timer('FFT-F xfft')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================

!$OMP ENDPARALLEL

1000 continue

!F  if (allocated(req_r))       deallocate(req_r)
!F  if (allocated(req_s))       deallocate(req_s)
!F  if (allocated(sta_r))       deallocate(sta_r)
!F  if (allocated(sta_s))       deallocate(sta_s)
!F  if (allocated(wk_recvcnt))    deallocate(wk_recvcnt)
!F  if (allocated(wk_sendcnt))    deallocate(wk_sendcnt)
!F  if (allocated(wk_recvrank))   deallocate(wk_recvrank)
!F  if (allocated(wk_sendrank))   deallocate(wk_sendrank)
!F  if (allocated(wk_afft_y))   deallocate(wk_afft_y)
    if (allocated(wk_mp_fft_y)) deallocate(wk_mp_fft_y)

#ifdef __TIMER_SUB__
    call timer_end(105)
#endif
  end subroutine m_FFT_Direct_3D

#else
!else ifdef FFT_FFTW_OLD

  subroutine m_FFT_Direct_3D(nfout,wk_afft_l, wk_size, iesize)
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
use mod_timer
#endif
! === TIMERTIMERTIMER ==========================================================
    integer, intent(in) :: nfout, wk_size, iesize
    real(kind=DP), dimension(wk_size*kimg,iesize), intent(inout) :: wk_afft_l

    real(kind=DP), allocatable, dimension(:) :: wk_allfft
    integer,      allocatable, dimension(:)   :: wk_recvdsp
    integer, parameter :: FFTW_MEASURE=0
    integer, parameter :: FFTW_ESTIMATE=64
    integer, parameter :: FFTW_RANK=1
    integer            :: FFTW_FLAG
    integer,dimension(FFTW_RANK) :: NFFTW3, NEMBED, NEREAL
    integer :: nx, ny, nz, i, j, k, l, ri, ib, iaddr, lrank
    integer :: irank, itag, icnt_send, icnt_recv
    integer,save, allocatable, dimension(:)   :: req_r, req_s
    integer,save, allocatable, dimension(:,:) :: sta_r, sta_s
    integer,save, allocatable, dimension(:)   :: wk_recvcnt, wk_sendcnt
!   integer,save, allocatable, dimension(:)   :: wk_recvrank, wk_sendrank
!F  real(kind=DP),allocatable, dimension(:,:) :: wk_recv, wk_send
    real(kind=DP),allocatable, dimension(:,:),save :: wk_recv1, wk_send1
    real(kind=DP),allocatable, dimension(:,:),save :: wk_recv2, wk_send2
    integer,      allocatable, dimension(:)   :: wk_mp_fft_y
    real(kind=DP),allocatable, dimension(:,:),save :: wk_afft_y
    real(kind=DP),allocatable, dimension(:,:,:)   :: wk_gather
    integer :: mpi_comm, myrank, nmrank
    integer(kind=8),save :: planx1 = 0, planx2 = 0
    integer(kind=8),save :: plany1 = 0, plany2 = 0
    integer(kind=8),save :: planz1 = 0, planz2 = 0
    integer(kind=8) :: plan
    integer :: max_x, max_y, max_z
    integer,save :: max_elm = 0, savesize = 0

    integer,save, allocatable, dimension(:,:) :: y2z_recv, y2z_send, z2x_recv, z2x_send
    integer,save :: y2z_rrank, y2z_srank, z2x_rrank, z2x_srank
    integer,save :: y2z_rmax, y2z_smax, y2z_srmax, z2x_rmax, z2x_smax, z2x_srmax

    integer ::  nnx, nny, nnz, iadd0, iadd1, nis, nie, jadd, kadd
#ifdef FFT_ALLTOALL
    integer,dimension(2) :: max_send,max_recv
#endif

#ifdef __TIMER_SUB__
    call mpi_barrier(mpi_ke_world,ierr)
    call timer_sta(105)
#endif

    if (kimg == 1) then
       FFTW_FLAG = +1      ! FFTW_BACKWARD
    else
       FFTW_FLAG = -1      ! FFTW_FORWARD
    endif
    plan = 0

! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
call start_timer('FFT-F reg0')
#endif
! === TIMERTIMERTIMER ==========================================================
!!  mpi_comm = mpi_kg_world
!!  myrank = myrank_e
!!  nmrank = nrank_e
    mpi_comm = mpi_ke_world
    myrank = myrank_g
    nmrank = nrank_g
    itag = 10

    if (np_fft_x==0)then
       nx = 0
       nx = 0
       nz = 0
       go to 1000
    endif

    if (firstcall_direct_3d) then
       savesize = 0
#ifdef __TIMER_DO__
  call timer_sta(118)
#endif
       max_x = maxval(nel_fft_x(:))
       max_y = maxval(nel_fft_y(:))
       max_z = maxval(nel_fft_z(:))
       max_elm = max(max_x,max_y,max_z)

       if(allocated(req_r)) deallocate(req_r)
       if(allocated(req_s)) deallocate(req_s)
       if(allocated(sta_r)) deallocate(sta_r)
       if(allocated(sta_s)) deallocate(sta_s)
       if(allocated(wk_recvcnt)) deallocate(wk_recvcnt)
       if(allocated(wk_sendcnt)) deallocate(wk_sendcnt)
       allocate(req_r(0:nmrank-1), stat=ierr)
       allocate(req_s(0:nmrank-1), stat=ierr)
       allocate(sta_r(MPI_STATUS_SIZE,0:nmrank-1), stat=ierr)
       allocate(sta_s(MPI_STATUS_SIZE,0:nmrank-1), stat=ierr)
       allocate(wk_recvcnt(0:nmrank-1), stat=ierr)
       allocate(wk_sendcnt(0:nmrank-1), stat=ierr)
!      allocate(wk_recvrank(0:nmrank-1), stat=ierr)
!      allocate(wk_sendrank(0:nmrank-1), stat=ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_FFT_Direct_3D :  Not allocate '
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 51, ierr)
        endif

       if(allocated(y2z_recv)) deallocate(y2z_recv)
       if(allocated(y2z_send)) deallocate(y2z_send)
       allocate(y2z_recv(2,0:nmrank-1))
       allocate(y2z_send(2,0:nmrank-1))
       y2z_send = 0
       y2z_recv = 0
       wk_recvcnt = 0
       wk_sendcnt = 0
       do i = 1, nel_fft_z(myrank)
          irank = map_fft_y(mp_fft_z(i)) - 1
          wk_recvcnt(irank) = wk_recvcnt(irank) + 1
       enddo
       do i = 1, nel_fft_y(myrank)
          irank = map_fft_z(mp_fft_y(i)) - 1
          wk_sendcnt(irank) = wk_sendcnt(irank) + 1
       enddo
       k = 0
       do i = 0, nmrank - 1
          if(wk_recvcnt(i) /= 0) then
             k = k + 1
             y2z_recv(1,i) = wk_recvcnt(i)
             y2z_recv(2,i) = k
          endif
       enddo
       y2z_rrank = k
       k = 0
       do i = 0, nmrank - 1
          if(wk_sendcnt(i) /= 0) then
             k = k + 1
             y2z_send(1,i) = wk_sendcnt(i)
             y2z_send(2,i) = k
          endif
       enddo
       y2z_srank = k
       y2z_rmax = maxval(wk_recvcnt)
       y2z_smax = maxval(wk_sendcnt)
       y2z_srmax = max(y2z_rmax,y2z_smax)

       if(allocated(z2x_recv)) deallocate(z2x_recv)
       if(allocated(z2x_send)) deallocate(z2x_send)
       allocate(z2x_recv(2,0:nmrank-1))
       allocate(z2x_send(2,0:nmrank-1))
       z2x_send = 0
       z2x_recv = 0
       wk_recvcnt = 0
       wk_sendcnt = 0
       do i = 1, nel_fft_x(myrank)
          irank = map_fft_z(mp_fft_x(i)) - 1
          wk_recvcnt(irank) = wk_recvcnt(irank) + 1
       enddo
       do i = 1, nel_fft_z(myrank)
          irank = map_fft_x(mp_fft_z(i)) - 1
          wk_sendcnt(irank) = wk_sendcnt(irank) + 1
       enddo
       k = 0
       do i = 0, nmrank - 1
          if(wk_recvcnt(i) /= 0) then
             k = k + 1
             z2x_recv(1,i) = wk_recvcnt(i)
             z2x_recv(2,i) = k
          endif
       enddo
       z2x_rrank = k
       k = 0
      do i = 0, nmrank - 1
          if(wk_sendcnt(i) /= 0) then
             k = k + 1
             z2x_send(1,i) = wk_sendcnt(i)
             z2x_send(2,i) = k
          endif
       enddo
       z2x_srank = k
       z2x_rmax = maxval(wk_recvcnt)
       z2x_smax = maxval(wk_sendcnt)
       z2x_srmax = max(z2x_rmax,z2x_smax)

#ifdef FFT_ALLTOALL
       max_send(1) = y2z_rmax
       max_send(2) = y2z_smax
       call mpi_allreduce(max_send,max_recv,2,mpi_integer,mpi_max,mpi_fft_zy_world,ierr)
       y2z_rmax = max_recv(1)
       y2z_smax = max_recv(2)
       max_send(1) = z2x_rmax
       max_send(2) = z2x_smax
       call mpi_allreduce(max_send,max_recv,2,mpi_integer,mpi_max,mpi_fft_xz_world,ierr)
       z2x_rmax = max_recv(1)
       z2x_smax = max_recv(2)
#endif

       if (ipri > 1) then
         write(nfout,'("m_FFT_Direct_3D   --   myrank_g=",i4)') myrank_g
         write(nfout,'("y2z_send")')
         write(nfout,'(10(i8,", "))') (y2z_send(1,i),i=0,nmrank-1)
         write(nfout,'(10(i8,", "))') (y2z_send(2,i),i=0,nmrank-1)
         write(nfout,'("y2z_recv")')
         write(nfout,'(10(i8,", "))') (y2z_recv(1,i),i=0,nmrank-1)
         write(nfout,'(10(i8,", "))') (y2z_recv(2,i),i=0,nmrank-1)
         write(nfout,'("z2x_send")')
         write(nfout,'(10(i8,", "))') (z2x_send(1,i),i=0,nmrank-1)
         write(nfout,'(10(i8,", "))') (z2x_send(2,i),i=0,nmrank-1)
         write(nfout,'("z2x_recv")')
         write(nfout,'(10(i8,", "))') (z2x_recv(1,i),i=0,nmrank-1)
         write(nfout,'(10(i8,", "))') (z2x_recv(2,i),i=0,nmrank-1)
       endif

#ifdef __TIMER_DO__
  call timer_sta(119)
#endif
       nx = xyz_fft_y(2,1)-xyz_fft_y(1,1)+1
       ny = xyz_fft_y(2,2)-xyz_fft_y(1,2)+1
       nz = xyz_fft_y(2,3)-xyz_fft_y(1,3)+1
!      if(kimg==1) then
          call dfftw_plan_many_dft    (plany1,         1,      ny, nx/2*nz, &
         &                             wk_afft_l,     ny,       1,      ny, &
         &                             wk_afft_l,     ny,       1,      ny, &
         &                             FFTW_FLAG, FFTW_MEASURE )
!      else
          call dfftw_plan_many_dft    (plany2,         1,      ny,   nz*nx, &
         &                             wk_afft_l,     ny,       1,      ny, &
         &                             wk_afft_l,     ny,       1,      ny, &
         &                             FFTW_FLAG, FFTW_MEASURE )
!      end if

       nx = xyz_fft_z(2,1)-xyz_fft_z(1,1)+1
       ny = xyz_fft_z(2,2)-xyz_fft_z(1,2)+1
       nz = xyz_fft_z(2,3)-xyz_fft_z(1,3)+1
!      if(kimg==1) then
          call dfftw_plan_many_dft    (planz1,         1,      nz, nx*ny/2, &
         &                             wk_afft_l,     nz, nx*ny/2,       1, &
         &                             wk_afft_l,     nz, nx*ny/2,       1, &
         &                             FFTW_FLAG, FFTW_MEASURE )
!      else
          call dfftw_plan_many_dft    (planz2,         1,      nz,   nx*ny, &
         &                             wk_afft_l,     nz,   nx*ny,       1, &
         &                             wk_afft_l,     nz,   nx*ny,       1, &
         &                             FFTW_FLAG, FFTW_MEASURE )
!      endif

       nx = xyz_fft_x(2,1)-xyz_fft_x(1,1)+1
       ny = xyz_fft_x(2,2)-xyz_fft_x(1,2)+1
       nz = xyz_fft_x(2,3)-xyz_fft_x(1,3)+1
!      if(kimg==1) then
          call dfftw_plan_many_dft_c2r(planx1,         1,    nx-2,   ny*nz, &
         &                             wk_afft_l,     nx,       1,    nx/2, &
         &                             wk_afft_l,     nx,       1,      nx, &
         &                             FFTW_MEASURE )
!      else
          call dfftw_plan_many_dft    (planx2,         1,      nx,   ny*nz, &
         &                             wk_afft_l,     nx,       1,      nx, &
         &                             wk_afft_l,     nx,       1,      nx, &
         &                             FFTW_FLAG, FFTW_MEASURE )
!      endif
#ifdef __TIMER_DO__
  call timer_end(119)
#endif

       firstcall_direct_3d = .false.
       go to 1000
#ifdef __TIMER_DO__
  call timer_end(118)
#endif
    endif

#ifdef FFT_ALLTOALL
    if (iesize /= savesize) then
#ifdef __TIMER_DO__
  call timer_sta(120)
#endif
       if (allocated(wk_recv1)) deallocate(wk_recv1)
       if (allocated(wk_send1)) deallocate(wk_send1)
       if (allocated(wk_recv2)) deallocate(wk_recv2)
       if (allocated(wk_send2)) deallocate(wk_send2)
       if (allocated(wk_afft_y)) deallocate(wk_afft_y)
       allocate(wk_recv1(y2z_rmax*kimg*iesize,y2z_rrank), stat=ierr)
       allocate(wk_send1(y2z_smax*kimg*iesize,y2z_srank), stat=ierr)
       allocate(wk_recv2(z2x_rmax*kimg*iesize,z2x_rrank), stat=ierr)
       allocate(wk_send2(z2x_smax*kimg*iesize,z2x_srank), stat=ierr)
       allocate(wk_afft_y(max_elm*kimg,iesize) ,stat=ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_FFT_Direct_3D :  Not allocate '
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 52, ierr)
        endif
       savesize = iesize
#ifdef __TIMER_DO__
  call timer_end(120)
#endif
    end if
#else
    if (iesize > savesize) then
#ifdef __TIMER_DO__
  call timer_sta(120)
#endif
       if (allocated(wk_recv1)) deallocate(wk_recv1)
       if (allocated(wk_send1)) deallocate(wk_send1)
       if (allocated(wk_recv2)) deallocate(wk_recv2)
       if (allocated(wk_send2)) deallocate(wk_send2)
       if (allocated(wk_afft_y)) deallocate(wk_afft_y)
       allocate(wk_recv1(y2z_rmax*kimg*iesize,y2z_rrank), stat=ierr)
       allocate(wk_send1(y2z_smax*kimg*iesize,y2z_srank), stat=ierr)
       allocate(wk_recv2(z2x_rmax*kimg*iesize,z2x_rrank), stat=ierr)
       allocate(wk_send2(z2x_smax*kimg*iesize,z2x_srank), stat=ierr)
       allocate(wk_afft_y(max_elm*kimg,iesize) ,stat=ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_FFT_Direct_3D :  Not allocate '
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 52, ierr)
        endif
       savesize = iesize
#ifdef __TIMER_DO__
  call timer_end(120)
#endif
    end if
#endif
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
call stop_timer('FFT-F reg0')
#endif
! === TIMERTIMERTIMER ==========================================================
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP          SHARED(     nfout,wk_afft_l,iesize,req_r,req_s,sta_r,sta_s,    &
!$OMP                 wk_recv1,wk_send1,wk_recv2,wk_send2,wk_afft_y,          &
!$OMP                 planx1,planx2,planz1,planz2,plany1,plany2,              &
!$OMP                 z2x_recv, z2x_send, y2z_recv, y2z_send,itag,            &
!$OMP                 xyz_fft_y,xyz_fft_z,xyz_fft_x,fft_Z_x_dim,nis_fft_Z_x,  &
!$OMP                 nie_fft_Z_x,fft_X_z_dim,nis_fft_X_z,nie_fft_X_z,        &
!$OMP                 fft_Y_z_dim,nis_fft_Y_z,nie_fft_Y_z,fft_Z_y_dim,        &
!$OMP                 nie_fft_Z_y,nis_fft_Z_y,                                &
#ifdef FFT_ALLTOALL
!$OMP                 mpi_fft_xz_world, z2x_rmax, z2x_smax,                   &
!$OMP                 mpi_fft_zy_world, y2z_rmax, y2z_smax,                   &
#endif
!$OMP                 mpi_comm,myrank,nmrank,kimg,ierr                    )   &
!$OMP          PRIVATE(nx,ny,nz,i,j,k,l,ri,ib,iaddr,lrank,icnt_send,icnt_recv, &
!$OMP                  irank,nnx,nny,nnz,iadd0,iadd1,nis,nie,jadd,kadd,plan   )

! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call start_timer('FFT-F yfft')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================
    nx = xyz_fft_y(2,1)-xyz_fft_y(1,1)+1
    ny = xyz_fft_y(2,2)-xyz_fft_y(1,2)+1
    nz = xyz_fft_y(2,3)-xyz_fft_y(1,3)+1
!
! Y-axis (z-x div)
!
    if(kimg==1) then
       plan = plany1
    else
       plan = plany2
    end if
#ifdef __TIMER_DO__
  call timer_sta(251)
#endif
!$OMP DO
    do ib = 1, iesize
       call dfftw_execute_dft(plan,wk_afft_l(1,ib),wk_afft_l(1,ib))
    end do
#ifdef __TIMER_DO__
  call timer_end(251)
#endif
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call stop_timer('FFT-F yfft')
call start_timer('FFT-F y2z')
call start_timer('*FFT-F y2z sbuf')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================

#ifdef __TIMER_DO__
  call timer_sta(121)
#endif
!fj --------------------
#ifdef FFT_INTERCHANGE
    if (kimg == 1) then
       do ib = 1, iesize
!$OMP DO
          do k = 1, nz
             do j = 1, ny
                do i = 1, nx/2
                   wk_afft_y((i+(j-1)*nx/2+(k-1)*nx/2*ny)*2-1,ib) = wk_afft_l((j+(i-1)*ny+(k-1)*ny*nx/2)*2-1,ib)
                   wk_afft_y((i+(j-1)*nx/2+(k-1)*nx/2*ny)*2  ,ib) = wk_afft_l((j+(i-1)*ny+(k-1)*ny*nx/2)*2  ,ib)
                enddo
             enddo
          enddo
       enddo
    else
       do ib = 1, iesize
!$OMP DO
          do k = 0, nz-1
             do j = 0, ny-1
                do i = 0, nx-1
                   wk_afft_y((i+j*nx+k*nx*ny+1)*2-1,ib) = wk_afft_l((j+i*ny+k*ny*nx+1)*2-1,ib)
                   wk_afft_y((i+j*nx+k*nx*ny+1)*2  ,ib) = wk_afft_l((j+i*ny+k*ny*nx+1)*2  ,ib)
                enddo
             enddo
          enddo
       enddo
    endif
#else
    if (kimg == 1) then
!$OMP DO
       do ib = 1, iesize
          do k = 1, nz
             do j = 1, ny
                do i = 1, nx/2
                   wk_afft_y((i+(j-1)*nx/2+(k-1)*nx/2*ny)*2-1,ib) = wk_afft_l((j+(i-1)*ny+(k-1)*ny*nx/2)*2-1,ib)
                   wk_afft_y((i+(j-1)*nx/2+(k-1)*nx/2*ny)*2  ,ib) = wk_afft_l((j+(i-1)*ny+(k-1)*ny*nx/2)*2  ,ib)
                enddo
             enddo
          enddo
       enddo
    else
!$OMP DO
       do ib = 1, iesize
          do k = 0, nz-1
             do j = 0, ny-1
                do i = 0, nx-1
                   wk_afft_y((i+j*nx+k*nx*ny+1)*2-1,ib) = wk_afft_l((j+i*ny+k*ny*nx+1)*2-1,ib)
                   wk_afft_y((i+j*nx+k*nx*ny+1)*2  ,ib) = wk_afft_l((j+i*ny+k*ny*nx+1)*2  ,ib)
                enddo
             enddo
          enddo
       enddo
    endif
#endif
!fj --------------------
#ifdef __TIMER_DO__
  call timer_end(121)
#endif

#ifdef __TIMER_DO__
  call timer_sta(122)
#endif

    nx = xyz_fft_y(2,1)-xyz_fft_y(1,1)+1
    ny = xyz_fft_y(2,2)-xyz_fft_y(1,2)+1
    nz = xyz_fft_y(2,3)-xyz_fft_y(1,3)+1
#ifdef FFT_INTERCHANGE
    if (kimg == 1) then
!$OMP DO
       do ib = 1, iesize
          do k = 1, nz
             do l = 1, fft_Z_y_dim
                nis = nis_fft_Z_y(l)
                nie = nie_fft_Z_y(l)
                nny = nie - nis + 1
                iadd0 = nx*nny*nz*(ib-1)
                do j = nis, nie
                   do i = 1, nx
                      wk_send1(iadd0+(i+(j-nis)*nx+(k-1)*nx*nny),l) = wk_afft_y(i+(j-1)*nx+(k-1)*nx*ny,ib)
                   end do
                end do
             end do
          end do
       end do
    else
!$OMP DO
       do ib = 1, iesize
          do k = 1, nz
             do l = 1, fft_Z_y_dim
                nis = nis_fft_Z_y(l)
                nie = nie_fft_Z_y(l)
                nny = nie - nis + 1
                do j = nis, nie
                   do i = 1, nx
                      iadd1 = nx*nny*nz*(ib-1)*2+(i+((j-nis)*nx+(k-1)*nx*nny))*2
                      iaddr = (i+(j-1)*nx+(k-1)*nx*ny)*2
                      wk_send1(iadd1-1,l) = wk_afft_y(iaddr-1,ib)
                      wk_send1(iadd1  ,l) = wk_afft_y(iaddr  ,ib)
                   end do
                end do
             end do
          end do
       end do
    end if
#else
!$OMP DO
    do l = 1, fft_Z_y_dim
       nis = nis_fft_Z_y(l)
       nie = nie_fft_Z_y(l)
       nny = nie - nis + 1
       if (kimg == 1) then
          do ib = 1, iesize
             iadd0 = nx*nny*nz*(ib-1)
             do k = 1, nz
                do j = nis, nie
                   do i = 1, nx
                      wk_send1(iadd0+(i+(j-nis)*nx+(k-1)*nx*nny),l) = wk_afft_y(i+(j-1)*nx+(k-1)*nx*ny,ib)
                   end do
                end do
             end do
          end do
       else
          do ib = 1, iesize
             iadd0 = nx*nny*nz*(ib-1)*2
             do k = 1, nz
                kadd = (k-1)*nx*nny
                do j = nis, nie
                   jadd = (j-nis)*nx+kadd
                   do i = 1, nx
                      iadd1 = iadd0+(i+jadd)*2
                      iaddr = (i+(j-1)*nx+(k-1)*nx*ny)*2
!fj --------------------
!                     do ri = 1, kimg
!                       wk_send1(iadd1+(ri-kimg),l) = wk_afft_y(iaddr+(ri-kimg),ib)
!                     end do
                      wk_send1(iadd1-1,l) = wk_afft_y(iaddr-1,ib)
                      wk_send1(iadd1  ,l) = wk_afft_y(iaddr  ,ib)
!fj --------------------
                   end do
                end do
             end do
          end do
       end if
    end do
#endif

#ifdef __TIMER_DO__
  call timer_end(122)
#endif
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call stop_timer('*FFT-F y2z sbuf')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================

#ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_fft_zy_world)
  call timer_sta(123)
#endif

! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call start_timer('*FFT-F y2z a2a')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================
!$OMP MASTER

    call MPI_ALLTOALL(wk_send1, y2z_smax*kimg*iesize, mpi_double_precision,   &
   &                  wk_recv1, y2z_rmax*kimg*iesize, mpi_double_precision,   &
   &                                                  mpi_fft_zy_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Direct_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 55, ierr)
     endif

!$OMP END MASTER
!$OMP BARRIER
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call stop_timer('*FFT-F y2z a2a')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================

#ifdef __TIMER_DO__
  call timer_sta(124)
#endif
#ifdef __TIMER_DO__
  call timer_end(124)
#endif

#else
!else ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_comm)
  call timer_sta(123)
#endif

! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call start_timer('*FFT-F y2z sr')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================
!$OMP MASTER
    icnt_recv = 0
    lrank = myrank + 1
    if (lrank > (nmrank-1)) lrank = 0
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank-1)) lrank = 0
       if ((lrank /= myrank) .and. (y2z_recv(1,lrank) /= 0)) then
          call mpi_irecv(wk_recv1(1,y2z_recv(2,lrank)), y2z_recv(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_r(icnt_recv), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Direct_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 55, ierr)
           endif
          icnt_recv = icnt_recv + 1
       endif
    enddo

    icnt_send = 0
    lrank = myrank
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank-1)) lrank = 0
       if ((lrank /= myrank) .and. (y2z_send(1,lrank) /= 0)) then
          call mpi_isend(wk_send1(1,y2z_send(2,lrank)), y2z_send(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_s(icnt_send), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Direct_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 56, ierr)
           endif
          icnt_send = icnt_send + 1
       endif
    enddo
!$OMP END MASTER
!$OMP BARRIER

#ifdef __TIMER_DO__
  call timer_sta(124)
#endif
!$OMP DO
    do i = 1, y2z_recv(1,myrank)*kimg*iesize
       wk_recv1(i,y2z_recv(2,myrank)) = wk_send1(i,y2z_send(2,myrank))
    enddo
#ifdef __TIMER_DO__
  call timer_end(124)
#endif

!$OMP MASTER
    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Direct_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 57, ierr)
     endif

    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Direct_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 58, ierr)
     endif
!$OMP END MASTER
!$OMP BARRIER
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call stop_timer('*FFT-F y2z sr')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================

#endif
!endif ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_end(123)
#endif

#ifdef __TIMER_DO__
  call timer_sta(125)
#endif
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call start_timer('*FFT-F y2z rbuf')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================
    nx = xyz_fft_z(2,1)-xyz_fft_z(1,1)+1
    ny = xyz_fft_z(2,2)-xyz_fft_z(1,2)+1
    nz = xyz_fft_z(2,3)-xyz_fft_z(1,3)+1
!$OMP DO
    do l = 1, fft_Y_z_dim
       nis = nis_fft_Y_z(l)
       nie = nie_fft_Y_z(l)
       nnz = nie - nis + 1
       if (kimg == 1) then
          do ib = 1, iesize
             iadd0 = nx*ny*nnz*(ib-1)
             do k = nis, nie
                do j = 1, ny
                   do i = 1, nx
                      wk_afft_l(i+(j-1)*nx+(k-1)*nx*ny,ib) = wk_recv1(iadd0+i+(j-1)*nx+(k-nis)*nx*ny,l)
                   end do
                end do
             end do
          end do
       else
          do ib = 1, iesize
             iadd0 = nx*ny*nnz*(ib-1)*2
             do k = nis, nie
                do j = 1, ny
                   do i = 1, nx
                      iadd1 = iadd0+(i+(j-1)*nx+(k-nis)*nx*ny)*2
                      iaddr = (i+(j-1)*nx+(k-1)*nx*ny)*2
                      wk_afft_l(iaddr-1,ib) = wk_recv1(iadd1-1,l)
                      wk_afft_l(iaddr  ,ib) = wk_recv1(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end if
    end do
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call stop_timer('*FFT-F y2z rbuf')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_DO__
  call timer_end(125)
#endif
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call stop_timer('FFT-F y2z')
call start_timer('FFT-F zfft')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================

!
! Z-axis (x-y div)
!
    nx = xyz_fft_z(2,1)-xyz_fft_z(1,1)+1
    ny = xyz_fft_z(2,2)-xyz_fft_z(1,2)+1
    nz = xyz_fft_z(2,3)-xyz_fft_z(1,3)+1
    if(kimg==1) then
       plan = planz1
    else
       plan = planz2
    end if
#ifdef __TIMER_DO__
  call timer_sta(252)
#endif
!$OMP DO
    do ib = 1, iesize
       call dfftw_execute_dft(plan,wk_afft_l(1,ib),wk_afft_l(1,ib))
    end do
#ifdef __TIMER_DO__
  call timer_end(252)
#endif
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call stop_timer('FFT-F zfft')
call start_timer('FFT-F z2x')
call start_timer('*FFT-F z2x sbuf')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================

#ifdef __TIMER_DO__
  call timer_sta(126)
#endif
    nx = xyz_fft_z(2,1)-xyz_fft_z(1,1)+1
    ny = xyz_fft_z(2,2)-xyz_fft_z(1,2)+1
    nz = xyz_fft_z(2,3)-xyz_fft_z(1,3)+1
!$OMP DO
    do l = 1, fft_X_z_dim
       nis = nis_fft_X_z(l)
       nie = nie_fft_X_z(l)
       nnz = nie - nis + 1
       if (kimg == 1) then
          do ib = 1, iesize
             iadd0 = nx*ny*nnz*(ib-1)
             do k = nis, nie
                do j = 1, ny
                   do i = 1, nx
                      wk_send2(iadd0+(i+(j-1)*nx+(k-nis)*nx*ny),l) = wk_afft_l(i+(j-1)*nx+(k-1)*nx*ny,ib)
                   end do
                end do
             end do
          end do
       else
          do ib = 1, iesize
             iadd0 = nx*ny*nnz*(ib-1)*2
             do k = nis, nie
                jadd = (k-nis)*nx*ny
                do j = 1, ny
                   kadd = (j-1)*nx+jadd
                   do i = 1, nx
                      iadd1 = iadd0+(i+kadd)*2
                      iaddr  = (i+(j-1)*nx+(k-1)*nx*ny)*2
!fj --------------------
!                     do ri = 1, kimg
!                        wk_send2(iadd1+(ri-kimg),l) = wk_afft_l(iaddr+(ri-kimg),ib)
!                     end do
                      wk_send2(iadd1-1,l) = wk_afft_l(iaddr-1,ib)
                      wk_send2(iadd1  ,l) = wk_afft_l(iaddr  ,ib)
!fj --------------------
                   end do
                end do
             end do
          end do
       end if
    end do

! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call stop_timer('*FFT-F z2x sbuf')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_DO__
  call timer_end(126)
#endif

#ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_fft_xz_world)
  call timer_sta(127)
#endif

! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call start_timer('*FFT-F z2x a2a')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================
!$OMP MASTER

    call MPI_ALLTOALL(wk_send2, z2x_smax*kimg*iesize, mpi_double_precision,   &
   &                  wk_recv2, z2x_rmax*kimg*iesize, mpi_double_precision,   &
   &                                                  mpi_fft_xz_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Direct_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 60, ierr)
     endif

!$OMP END MASTER
!$OMP BARRIER
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call stop_timer('*FFT-F z2x a2a')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================

#ifdef __TIMER_DO__
  call timer_sta(128)
#endif
#ifdef __TIMER_DO__
  call timer_end(128)
#endif

#else
!else ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_comm)
  call timer_sta(127)
#endif

! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call start_timer('*FFT-F z2x sr')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================
!$OMP MASTER
    icnt_recv = 0
    lrank = myrank + 1
    if (lrank > (nmrank-1)) lrank = 0
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank-1)) lrank = 0
       if ((lrank /= myrank) .and. (z2x_recv(1,lrank) /= 0)) then
          call mpi_irecv(wk_recv2(1,z2x_recv(2,lrank)), z2x_recv(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_r(icnt_recv), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Direct_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 60, ierr)
           endif
          icnt_recv = icnt_recv + 1
       endif
    enddo

    icnt_send = 0
    lrank = myrank
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank-1)) lrank = 0
       if ((lrank /= myrank) .and. (z2x_send(1,lrank) /= 0)) then
          call mpi_isend(wk_send2(1,z2x_send(2,lrank)), z2x_send(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_s(icnt_send), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Direct_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 61, ierr)
           endif
          icnt_send = icnt_send + 1
       endif
    enddo
!$OMP END MASTER
!$OMP BARRIER

#ifdef __TIMER_DO__
  call timer_sta(128)
#endif
!$OMP DO
    do i = 1, z2x_recv(1,myrank)*kimg*iesize
       wk_recv2(i,z2x_recv(2,myrank)) = wk_send2(i,z2x_recv(2,myrank))
    enddo
#ifdef __TIMER_DO__
  call timer_end(128)
#endif

!$OMP MASTER
    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Direct_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 62, ierr)
     endif

    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Direct_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 63, ierr)
     endif
!$OMP END MASTER
!$OMP BARRIER
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call stop_timer('*FFT-F z2x sr')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================
#endif
!endif ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_end(127)
#endif

#ifdef __TIMER_DO__
  call timer_sta(129)
#endif
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call start_timer('*FFT-F z2x rbuf')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================
    nx = xyz_fft_x(2,1)-xyz_fft_x(1,1)+1
    ny = xyz_fft_x(2,2)-xyz_fft_x(1,2)+1
    nz = xyz_fft_x(2,3)-xyz_fft_x(1,3)+1
#ifdef FFT_INTERCHANGE
    if (kimg == 1) then
!$OMP DO
       do ib = 1, iesize
          do k = 1, nz
             do j = 1, ny
                do l = 1, fft_Z_x_dim
                   nis = nis_fft_Z_x(l)
                   nie = nie_fft_Z_x(l)
                   nnx = nie - nis + 1
                   iadd0 = nnx*ny*nz*(ib-1)
                   do i = nis, nie
                      wk_afft_l(i+(j-1)*nx+(k-1)*nx*ny,ib) = wk_recv2(iadd0+i-nis+1+(j-1)*nnx+(k-1)*nnx*ny,l)
                   end do
                end do
             end do
          end do
       end do
    else
!$OMP DO
       do ib = 1, iesize
          do k = 1, nz
             do j = 1, ny
                do l = 1, fft_Z_x_dim
                   nis = nis_fft_Z_x(l)
                   nie = nie_fft_Z_x(l)
                   nnx = nie - nis + 1
                   do i = nis, nie
                      iadd1 = nnx*ny*nz*(ib-1)*2+(i-nis+1+(j-1)*nnx+(k-1)*nnx*ny)*2
                      iaddr = (i+(j-1)*nx+(k-1)*nx*ny)*2
                      wk_afft_l(iaddr-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_l(iaddr  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end do
    end if
#else
!$OMP DO
    do l = 1, fft_Z_x_dim
       nis = nis_fft_Z_x(l)
       nie = nie_fft_Z_x(l)
       nnx = nie - nis + 1
       if (kimg == 1) then
          do ib = 1, iesize
             iadd0 = nnx*ny*nz*(ib-1)
             do k = 1, nz
                do j = 1, ny
                   do i = nis, nie
                      wk_afft_l(i+(j-1)*nx+(k-1)*nx*ny,ib) = wk_recv2(iadd0+i-nis+1+(j-1)*nnx+(k-1)*nnx*ny,l)
                   end do
                end do
             end do
          end do
       else
          do ib = 1, iesize
             iadd0 = nnx*ny*nz*(ib-1)*2
             do k = 1, nz
                do j = 1, ny
                   do i = nis, nie
                      iadd1 = iadd0+(i-nis+1+(j-1)*nnx+(k-1)*nnx*ny)*2
                      iaddr = (i+(j-1)*nx+(k-1)*nx*ny)*2
                      wk_afft_l(iaddr-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_l(iaddr  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end if
    end do
#endif
#ifdef __TIMER_DO__
  call timer_end(129)
#endif
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call stop_timer('*FFT-F z2x rbuf')
call stop_timer('FFT-F z2x')
call start_timer('FFT-F xfft')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================

!
! X-axis (y-z div)
!
    nx = xyz_fft_x(2,1)-xyz_fft_x(1,1)+1
    ny = xyz_fft_x(2,2)-xyz_fft_x(1,2)+1
    nz = xyz_fft_x(2,3)-xyz_fft_x(1,3)+1
    if(kimg==1) then
       plan = planx1
    else
       plan = planx2
    end if
#ifdef __TIMER_DO__
  call timer_sta(253)
#endif
!$OMP DO
    do ib = 1, iesize
       if (kimg==1) then
          call dfftw_execute_dft_c2r(plan,wk_afft_l(1,ib),wk_afft_l(1,ib))
       else
          call dfftw_execute_dft    (plan,wk_afft_l(1,ib),wk_afft_l(1,ib))
       endif
    end do
#ifdef __TIMER_DO__
  call timer_end(253)
#endif
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call stop_timer('FFT-F xfft')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================

!$OMP ENDPARALLEL

1000 continue

!F  if (allocated(req_r))       deallocate(req_r)
!F  if (allocated(req_s))       deallocate(req_s)
!F  if (allocated(sta_r))       deallocate(sta_r)
!F  if (allocated(sta_s))       deallocate(sta_s)
!F  if (allocated(wk_recvcnt))    deallocate(wk_recvcnt)
!F  if (allocated(wk_sendcnt))    deallocate(wk_sendcnt)
!F  if (allocated(wk_recvrank))   deallocate(wk_recvrank)
!F  if (allocated(wk_sendrank))   deallocate(wk_sendrank)
!F  if (allocated(wk_afft_y))   deallocate(wk_afft_y)
    if (allocated(wk_mp_fft_y)) deallocate(wk_mp_fft_y)

#ifdef __TIMER_SUB__
    call timer_end(105)
#endif
  end subroutine m_FFT_Direct_3D

#endif
!endif ifdef FFT_FFTW_OLD

!------------------------------------------------------------------------------
#else
!ifndef FFT_USE_SSL2

!------------------------------------------------------------------------------
! For Fujitsu SSL2
  subroutine m_FFT_Direct_3D(nfout,wk_afft_l, wk_size, iesize)
    integer, intent(in) :: nfout, wk_size, iesize
    real(kind=DP), dimension(wk_size*kimg,iesize), intent(inout) :: wk_afft_l

    real(kind=DP), allocatable, dimension(:) :: wk_allfft
    integer,      allocatable, dimension(:)   :: wk_recvdsp
    integer, parameter :: FFTW_MEASURE=0
    integer, parameter :: FFTW_ESTIMATE=64
    integer, parameter :: FFTW_RANK=1
    integer            :: FFTW_FLAG
    integer,dimension(FFTW_RANK) :: NFFTW3, NEMBED, NEREAL
    integer :: nx, ny, nz, i, j, k, l, ri, ib, iadd, lrank
    integer :: irank, itag, icnt_send, icnt_recv
    integer,save, allocatable, dimension(:)   :: req_r, req_s
    integer,save, allocatable, dimension(:,:) :: sta_r, sta_s
    integer,save, allocatable, dimension(:)   :: wk_recvcnt, wk_sendcnt
!   integer,save, allocatable, dimension(:)   :: wk_recvrank, wk_sendrank
!F  real(kind=DP),allocatable, dimension(:,:) :: wk_recv, wk_send
    real(kind=DP),allocatable, dimension(:,:),save :: wk_recv1, wk_send1
    real(kind=DP),allocatable, dimension(:,:),save :: wk_recv2, wk_send2
    integer,      allocatable, dimension(:)   :: wk_mp_fft_y
    real(kind=DP),allocatable, dimension(:,:),save :: wk_afft_y
    real(kind=DP),allocatable, dimension(:,:,:)   :: wk_gather
    integer :: mpi_comm, myrank, nmrank
    integer(kind=8),save :: planx1 = 0, planx2 = 0
    integer(kind=8),save :: plany1 = 0, plany2 = 0
    integer(kind=8),save :: planz1 = 0, planz2 = 0
    integer(kind=8) :: plan
    integer :: max_x, max_y, max_z
    integer,save :: max_elm = 0, savesize = 0

    integer,save, allocatable, dimension(:,:) :: y2z_recv, y2z_send, z2x_recv, z2x_send
    integer,save :: y2z_rrank, y2z_srank, z2x_rrank, z2x_srank
    integer,save :: y2z_rmax, y2z_smax, y2z_srmax, z2x_rmax, z2x_smax, z2x_srmax

    integer ::  nnx, nny, nnz, iadd0, iadd1, nis, nie, jadd, kadd, nisx, niex, nxp, nyp, nzp

    integer,dimension(3) :: nsize, isin
#ifdef FFT_ALLTOALL
    integer,dimension(2) :: max_send,max_recv
#endif

#ifdef __TIMER_SUB__
    call mpi_barrier(mpi_ke_world,ierr)
    call timer_sta(105)
#endif

    if (kimg == 1) then
       FFTW_FLAG = +1      ! FFTW_BACKWARD
    else
       FFTW_FLAG = -1      ! FFTW_FORWARD
    endif
    plan = 0

!!  mpi_comm = mpi_kg_world
!!  myrank = myrank_e
!!  nmrank = nrank_e
    mpi_comm = mpi_ke_world
    myrank = myrank_g
    nmrank = nrank_g
    itag = 10

    if (np_fft_x==0)then
       nx = 0
       nx = 0
       nz = 0
       go to 1000
    endif

    if (firstcall_direct_3d) then
       savesize = 0
#ifdef __TIMER_DO__
  call timer_sta(118)
#endif
       max_x = maxval(nel_fft_x(:))
       max_y = maxval(nel_fft_y(:))
       max_z = maxval(nel_fft_z(:))
       max_elm = max(max_x,max_y,max_z)

       if(allocated(req_r)) deallocate(req_r)
       if(allocated(req_s)) deallocate(req_s)
       if(allocated(sta_r)) deallocate(sta_r)
       if(allocated(sta_s)) deallocate(sta_s)
       if(allocated(wk_recvcnt)) deallocate(wk_recvcnt)
       if(allocated(wk_sendcnt)) deallocate(wk_sendcnt)

       allocate(req_r(0:nmrank-1), stat=ierr)
       allocate(req_s(0:nmrank-1), stat=ierr)
       allocate(sta_r(MPI_STATUS_SIZE,0:nmrank-1), stat=ierr)
       allocate(sta_s(MPI_STATUS_SIZE,0:nmrank-1), stat=ierr)
       allocate(wk_recvcnt(0:nmrank-1), stat=ierr)
       allocate(wk_sendcnt(0:nmrank-1), stat=ierr)
!      allocate(wk_recvrank(0:nmrank-1), stat=ierr)
!      allocate(wk_sendrank(0:nmrank-1), stat=ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_FFT_Direct_3D :  Not allocate '
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 51, ierr)
        endif

       if(allocated(y2z_recv)) deallocate(y2z_recv)
       if(allocated(y2z_send)) deallocate(y2z_send)
       allocate(y2z_recv(2,0:nmrank-1))
       allocate(y2z_send(2,0:nmrank-1))
       y2z_send = 0
       y2z_recv = 0
       wk_recvcnt = 0
       wk_sendcnt = 0
       do i = 1, nel_fft_z(myrank)
          irank = map_fft_y(mp_fft_z(i)) - 1
          wk_recvcnt(irank) = wk_recvcnt(irank) + 1
       enddo
       do i = 1, nel_fft_y(myrank)
          irank = map_fft_z(mp_fft_y(i)) - 1
          wk_sendcnt(irank) = wk_sendcnt(irank) + 1
       enddo
       k = 0
       do i = 0, nmrank - 1
          if(wk_recvcnt(i) /= 0) then
             k = k + 1
             y2z_recv(1,i) = wk_recvcnt(i)
             y2z_recv(2,i) = k
          endif
       enddo
       y2z_rrank = k
       k = 0
       do i = 0, nmrank - 1
          if(wk_sendcnt(i) /= 0) then
             k = k + 1
             y2z_send(1,i) = wk_sendcnt(i)
             y2z_send(2,i) = k
          endif
       enddo
       y2z_srank = k
       y2z_rmax = maxval(wk_recvcnt)
       y2z_smax = maxval(wk_sendcnt)
       y2z_srmax = max(y2z_rmax,y2z_smax)

       if(allocated(z2x_recv)) deallocate(z2x_recv)
       if(allocated(z2x_send)) deallocate(z2x_send)
       allocate(z2x_recv(2,0:nmrank-1))
       allocate(z2x_send(2,0:nmrank-1))
       z2x_send = 0
       z2x_recv = 0
       wk_recvcnt = 0
       wk_sendcnt = 0
       do i = 1, nel_fft_x(myrank)
          irank = map_fft_z(mp_fft_x(i)) - 1
          wk_recvcnt(irank) = wk_recvcnt(irank) + 1
       enddo
       do i = 1, nel_fft_z(myrank)
          irank = map_fft_x(mp_fft_z(i)) - 1
          wk_sendcnt(irank) = wk_sendcnt(irank) + 1
       enddo
       k = 0
       do i = 0, nmrank - 1
          if(wk_recvcnt(i) /= 0) then
             k = k + 1
             z2x_recv(1,i) = wk_recvcnt(i)
             z2x_recv(2,i) = k
          endif
       enddo
       z2x_rrank = k
       k = 0
      do i = 0, nmrank - 1
          if(wk_sendcnt(i) /= 0) then
             k = k + 1
             z2x_send(1,i) = wk_sendcnt(i)
             z2x_send(2,i) = k
          endif
       enddo
       z2x_srank = k
       z2x_rmax = maxval(wk_recvcnt)
       z2x_smax = maxval(wk_sendcnt)
       z2x_srmax = max(z2x_rmax,z2x_smax)

#ifdef FFT_ALLTOALL
       max_send(1) = y2z_rmax
       max_send(2) = y2z_smax
       call mpi_allreduce(max_send,max_recv,2,mpi_integer,mpi_max,mpi_fft_zy_world,ierr)
       y2z_rmax = max_recv(1)
       y2z_smax = max_recv(2)
       max_send(1) = z2x_rmax
       max_send(2) = z2x_smax
       call mpi_allreduce(max_send,max_recv,2,mpi_integer,mpi_max,mpi_fft_xz_world,ierr)
       z2x_rmax = max_recv(1)
       z2x_smax = max_recv(2)
#endif

       if (ipri > 1) then
          write(nfout,'("m_FFT_Direct_3D   --   myrank_g=",i4)') myrank_g
          write(nfout,'("y2z_send")')
          write(nfout,'(10(i8,", "))') (y2z_send(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (y2z_send(2,i),i=0,nmrank-1)
          write(nfout,'("y2z_recv")')
          write(nfout,'(10(i8,", "))') (y2z_recv(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (y2z_recv(2,i),i=0,nmrank-1)
          write(nfout,'("z2x_send")')
          write(nfout,'(10(i8,", "))') (z2x_send(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (z2x_send(2,i),i=0,nmrank-1)
          write(nfout,'("z2x_recv")')
          write(nfout,'(10(i8,", "))') (z2x_recv(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (z2x_recv(2,i),i=0,nmrank-1)
       endif

       nx = xyz_fft_y(2,1)-xyz_fft_y(1,1)+1
       ny = xyz_fft_y(2,2)-xyz_fft_y(1,2)+1
       nz = xyz_fft_y(2,3)-xyz_fft_y(1,3)+1
!!     if(kimg==1) then
!X        NFFTW3(1) = ny
!X        NEMBED(1) = ny
!X        call dfftw_plan_many_dft    (plany1, FFTW_RANK, NFFTW3, nz*nx/2, &
!X       &                             wk_afft_l, NEMBED, nz*nx/2, 1,    &
!X       &                             wk_afft_l, NEMBED, nz*nx/2, 1,    &
!X       &                             FFTW_FLAG, FFTW_ESTIMATE)
!!     else
!X        NFFTW3(1) = ny
!X        NEMBED(1) = ny
!X        call dfftw_plan_many_dft    (plany2, FFTW_RANK, NFFTW3, nz*nx,   &
!X       &                             wk_afft_l, NEMBED, nz*nx, 1,      &
!X       &                             wk_afft_l, NEMBED, nz*nx, 1,      &
!X       &                             FFTW_FLAG, FFTW_ESTIMATE)
!!     end if

!X     nx = xyz_fft_z(2,1)-xyz_fft_z(1,1)+1
!X     ny = xyz_fft_z(2,2)-xyz_fft_z(1,2)+1
!X     nz = xyz_fft_z(2,3)-xyz_fft_z(1,3)+1
!!     if(kimg==1) then
!X        NFFTW3(1) = nz
!X        NEMBED(1) = nz
!X        call dfftw_plan_many_dft    (planz1, FFTW_RANK, NFFTW3, nx*ny/2, &
!X       &                             wk_afft_l, NEMBED, nx*ny/2, 1,     &
!X       &                             wk_afft_l, NEMBED, nx*ny/2, 1,     &
!X       &                             FFTW_FLAG, FFTW_ESTIMATE)
!!     else
!X        NFFTW3(1) = nz
!X        NEMBED(1) = nz
!X        call dfftw_plan_many_dft    (planz2, FFTW_RANK, NFFTW3, nx*ny,   &
!X       &                             wk_afft_l, NEMBED, nx*ny, 1,       &
!X       &                             wk_afft_l, NEMBED, nx*ny, 1,       &
!X       &                             FFTW_FLAG, FFTW_ESTIMATE)
!!     endif

!X     nx = xyz_fft_x(2,1)-xyz_fft_x(1,1)+1
!X     ny = xyz_fft_x(2,2)-xyz_fft_x(1,2)+1
!X     nz = xyz_fft_x(2,3)-xyz_fft_x(1,3)+1
!!     if(kimg==1) then
!X        NFFTW3(1) = nx - 2
!X        NEMBED(1) = nx
!X        NEREAL(1) = nx
!X        call dfftw_plan_many_dft_c2r(planx1, FFTW_RANK, NFFTW3, ny*nz,   &
!X       &                             wk_afft_l, NEMBED, 1, nx/2,        &
!X       &                             wk_afft_l, NEREAL, 1, nx,          &
!X       &                             FFTW_ESTIMATE)
!!     else
!X        NFFTW3(1) = nx
!X        NEMBED(1) = nx
!X        call dfftw_plan_many_dft    (planx2, FFTW_RANK, NFFTW3, ny*nz,   &
!X       &                             wk_afft_l, NEMBED, 1, nx,          &
!X       &                             wk_afft_l, NEMBED, 1, nx,          &
!X       &                             FFTW_FLAG, FFTW_ESTIMATE)
!!     endif

       firstcall_direct_3d = .false.
#ifdef __TIMER_DO__
  call timer_end(118)
#endif
    endif

#ifdef FFT_ALLTOALL
    if (iesize /= savesize) then
#ifdef __TIMER_DO__
  call timer_sta(120)
#endif
       if (allocated(wk_recv1)) deallocate(wk_recv1)
       if (allocated(wk_send1)) deallocate(wk_send1)
       if (allocated(wk_recv2)) deallocate(wk_recv2)
       if (allocated(wk_send2)) deallocate(wk_send2)
       if (allocated(wk_afft_y)) deallocate(wk_afft_y)
       allocate(wk_recv1(y2z_rmax*kimg*iesize,y2z_rrank), stat=ierr)
       allocate(wk_send1(y2z_smax*kimg*iesize,y2z_srank), stat=ierr)
       allocate(wk_recv2(z2x_rmax*kimg*iesize,z2x_rrank), stat=ierr)
       allocate(wk_send2(z2x_smax*kimg*iesize,z2x_srank), stat=ierr)
       allocate(wk_afft_y(max_elm*kimg,iesize) ,stat=ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_FFT_Direct_3D :  Not allocate '
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 52, ierr)
        endif
       savesize = iesize
#ifdef __TIMER_DO__
  call timer_end(120)
#endif
    end if
#else
    if (iesize > savesize) then
#ifdef __TIMER_DO__
  call timer_sta(120)
#endif
       if (allocated(wk_recv1)) deallocate(wk_recv1)
       if (allocated(wk_send1)) deallocate(wk_send1)
       if (allocated(wk_recv2)) deallocate(wk_recv2)
       if (allocated(wk_send2)) deallocate(wk_send2)
       if (allocated(wk_afft_y)) deallocate(wk_afft_y)
       allocate(wk_recv1(y2z_rmax*kimg*iesize,y2z_rrank), stat=ierr)
       allocate(wk_send1(y2z_smax*kimg*iesize,y2z_srank), stat=ierr)
       allocate(wk_recv2(z2x_rmax*kimg*iesize,z2x_rrank), stat=ierr)
       allocate(wk_send2(z2x_smax*kimg*iesize,z2x_srank), stat=ierr)
       allocate(wk_afft_y(max_elm*kimg,iesize) ,stat=ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_FFT_Direct_3D :  Not allocate '
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 52, ierr)
        endif
       savesize = iesize
#ifdef __TIMER_DO__
  call timer_end(120)
#endif
    end if
#endif

!
! Y-axis (z-x div)
!
#ifdef __TIMER_DO__
  call timer_sta(120)
#endif

    nx = xyz_fft_y(2,1)-xyz_fft_y(1,1)+1
    ny = xyz_fft_y(2,2)-xyz_fft_y(1,2)+1
    nz = xyz_fft_y(2,3)-xyz_fft_y(1,3)+1
    nyp = ny

    wk_afft_y(1:wk_size*kimg,1:iesize) = wk_afft_l(1:wk_size*kimg,1:iesize)
! (x,z,y) -> (y,x,z)
    if (kimg == 1) then
!OCL PARALLEL_STRONG
       do ib = 1, iesize
          do k = 1, nz
             do i = 1, nx/2
                do j = 1, ny
                   iadd0 = (i+(k-1)*nx/2+(j-1)*nx/2*nz)*2
                   iadd1 = (j+(i-1)*nyp+(k-1)*nyp*nx/2)*2
                   wk_afft_l(iadd1-1,ib) = wk_afft_y(iadd0-1,ib)
                   wk_afft_l(iadd1  ,ib) = wk_afft_y(iadd0  ,ib)
                enddo
             enddo
          enddo
       enddo
    else
!OCL PARALLEL_STRONG
       do ib = 1, iesize
!OCL SERIAL
          do k = 0, nz-1
!OCL SERIAL
             do i = 0, nx-1
!OCL SERIAL
                do j = 0, ny-1
                   wk_afft_l((j+i*nyp+k*nyp*nx+1)*2-1,ib) = wk_afft_y((i+k*nx+j*nx*nz+1)*2-1,ib)
                   wk_afft_l((j+i*nyp+k*nyp*nx+1)*2  ,ib) = wk_afft_y((i+k*nx+j*nx*nz+1)*2  ,ib)
                enddo
             enddo
          enddo
       enddo
    end if

#ifdef __TIMER_DO__
  call timer_end(120)
#endif
!
! Y-axis (x-z div)
!
#ifdef __TIMER_DO__
  call timer_sta(251)
#endif
    if(kimg==1) then
       nsize(1:3) = (/ny,(nx/2),nz/)
       isin(1:3)  = (/-1,0,0/)
!OCL INDEPENDENT[ DVMCF2 ]
!OCL PARALLEL_STRONG
       do ib = 1, iesize
          call DVMCF2(wk_afft_l(1,ib),nsize,3,isin,ierr)
       end do
    else
       nsize(1:3) = (/ny,nx,nz/)
       isin(1:3)  = (/1,0,0/)
!OCL INDEPENDENT[ DVMCF2 ]
!OCL PARALLEL_STRONG
       do ib = 1, iesize
          call DVMCF2(wk_afft_l(1,ib),nsize,3,isin,ierr)
       end do
    end if
#ifdef __TIMER_DO__
  call timer_end(251)
#endif

#ifdef __TIMER_DO__
  call timer_sta(121)
#endif
#ifdef __TIMER_DO__
  call timer_end(121)
#endif

#ifdef __TIMER_DO__
  call timer_sta(122)
#endif
    if (kimg == 1) then
!OCL NORECURRENCE
!OCL PARALLEL_STRONG
       do ib = 1, iesize
!OCL SERIAL
          do k = 1, nz
!OCL SERIAL
             do i = 1, nx/2
!OCL SERIAL
                do l = 1, fft_Z_y_dim
                   nis = nis_fft_Z_y(l)
                   nie = nie_fft_Z_y(l)
                   nny = nie - nis + 1
                   iadd0 = nx/2*nny*nz*(ib-1)*2
!OCL SERIAL
                   do j = nis, nie
                      iadd  = (j+(i-1)*nyp+(k-1)*nyp*nx/2)*2
                      iadd1 = iadd0+((j-nis+1)+(i-1)*nny+(k-1)*nny*nx/2)*2
                      wk_send1(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send1(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
    else
!OCL NORECURRENCE
!OCL PARALLEL_STRONG
       do ib = 1, iesize
!OCL SERIAL
          do k = 1, nz
!OCL SERIAL
             do i = 1, nx
!OCL SERIAL
                do l = 1, fft_Z_y_dim
                   nis = nis_fft_Z_y(l)
                   nie = nie_fft_Z_y(l)
                   nny = nie - nis + 1
                   iadd0 = nx*nny*nz*(ib-1)*2
!OCL SERIAL
                   do j = nis, nie
                      iadd  = (j+(i-1)*nyp+(k-1)*nyp*nx)*2
                      iadd1 = iadd0+((j-nis+1)+(i-1)*nny+(k-1)*nny*nx)*2
                      wk_send1(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send1(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
    end if
#ifdef __TIMER_DO__
  call timer_end(122)
#endif

#ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_fft_zy_world)
  call timer_sta(123)
#endif

    call MPI_ALLTOALL(wk_send1, y2z_smax*kimg*iesize, mpi_double_precision,   &
   &                  wk_recv1, y2z_rmax*kimg*iesize, mpi_double_precision,   &
   &                                                  mpi_fft_zy_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Direct_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 55, ierr)
     endif

#ifdef __TIMER_DO__
  call timer_sta(124)
#endif
#ifdef __TIMER_DO__
  call timer_end(124)
#endif

#else
!else ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_comm)
  call timer_sta(123)
#endif

    icnt_recv = 0
    lrank = myrank + 1
    if (lrank > (nmrank-1)) lrank = 0
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank-1)) lrank = 0
       if ((y2z_recv(1,lrank) /= 0)) then
          call mpi_irecv(wk_recv1(1,y2z_recv(2,lrank)), y2z_recv(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_r(icnt_recv), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Direct_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 55, ierr)
           endif
          icnt_recv = icnt_recv + 1
       endif
    enddo

    icnt_send = 0
    lrank = myrank
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank-1)) lrank = 0
       if ((y2z_send(1,lrank) /= 0)) then
          call mpi_isend(wk_send1(1,y2z_send(2,lrank)), y2z_send(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_s(icnt_send), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Direct_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 56, ierr)
           endif
          icnt_send = icnt_send + 1
       endif
    enddo

#ifdef __TIMER_DO__
  call timer_sta(124)
#endif
#ifdef __TIMER_DO__
  call timer_end(124)
#endif

    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Direct_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 57, ierr)
     endif

    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Direct_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 58, ierr)
     endif

#endif
!endif ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_end(123)
#endif

#ifdef __TIMER_DO__
  call timer_sta(125)
#endif
! (y,x,z) -> (z,x,y)
    nx = xyz_fft_z(2,1)-xyz_fft_z(1,1)+1
    ny = xyz_fft_z(2,2)-xyz_fft_z(1,2)+1
    nz = xyz_fft_z(2,3)-xyz_fft_z(1,3)+1
    nzp = nz
    if (kimg == 1) then
       do ib = 1, iesize
          do j = 1, ny
             do i = 1, nx/2
                do l = 1, fft_Y_z_dim
                   nis = nis_fft_Y_z(l)
                   nie = nie_fft_Y_z(l)
                   nnz = nie - nis + 1
                   iadd0 = nx/2*ny*nnz*(ib-1)*2
                   do k = nis, nie
                      iadd1 = iadd0 + (j+(i-1)*ny+(k-nis)*ny*nx/2)*2
                      iadd = (k+(i-1)*nzp+(j-1)*nzp*nx/2)*2
                      wk_afft_l(iadd-1,ib) = wk_recv1(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv1(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end do
    else
       do ib = 1, iesize
          do i = 1, nx
             do j = 1, ny
                do l = 1, fft_Y_z_dim
                   nis = nis_fft_Y_z(l)
                   nie = nie_fft_Y_z(l)
                   nnz = nie - nis + 1
                   iadd0 = nx*ny*nnz*(ib-1)*2
                   do k = nis, nie
                      iadd1 = iadd0+(j+(i-1)*ny+(k-nis)*ny*nx)*2
                      iadd  = (k+(i-1)*nzp+(j-1)*nzp*nx)*2
                      wk_afft_l(iadd-1,ib) = wk_recv1(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv1(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end do
    end if
#ifdef __TIMER_DO__
  call timer_end(125)
#endif
!
! Z-axis (x-y div)
!
#ifdef __TIMER_DO__
  call timer_sta(252)
#endif
    if(kimg==1) then
       nsize(1:3) = (/nz,(nx/2),ny/)
       isin(1:3)  = (/-1,0,0/)
!OCL INDEPENDENT[ DVMCF2 ]
!OCL PARALLEL_STRONG
       do ib = 1, iesize
          call DVMCF2(wk_afft_l(1,ib),nsize,3,isin,ierr)
       end do
    else
       nsize(1:3) = (/nz,nx,ny/)
       isin(1:3)  = (/1,0,0/)
!OCL INDEPENDENT[ DVMCF2 ]
!OCL PARALLEL_STRONG
       do ib = 1, iesize
          call DVMCF2(wk_afft_l(1,ib),nsize,3,isin,ierr)
       end do
    end if
#ifdef __TIMER_DO__
  call timer_end(252)
#endif

#ifdef __TIMER_DO__
  call timer_sta(126)
#endif
    if (kimg == 1) then
!OCL NORECURRENCE
!OCL PARALLEL_STRONG
       do ib = 1, iesize
!OCL SERIAL
          do j = 1, ny
!OCL SERIAL
             do i = 1, nx/2
!OCL SERIAL
                do l = 1, fft_X_z_dim
                   nis = nis_fft_X_z(l)
                   nie = nie_fft_X_z(l)
                   nnz = nie - nis + 1
                   iadd0 = nx/2*ny*nnz*(ib-1)*2
!OCL SERIAL
                   do k = nis, nie
                      iadd  = (k+(i-1)*nzp+(j-1)*nzp*nx/2)*2
                      iadd1 = iadd0 + (k-nis+1+(i-1)*nnz+(j-1)*nnz*nx/2)*2
                      wk_send2(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send2(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
    else
!OCL NORECURRENCE
!OCL PARALLEL_STRONG
       do ib = 1, iesize
!OCL SERIAL
          do j = 1, ny
!OCL SERIAL
             do i = 1, nx
!OCL SERIAL
                do l = 1, fft_X_z_dim
                   nis = nis_fft_X_z(l)
                   nie = nie_fft_X_z(l)
                   nnz = nie - nis + 1
                   iadd0 = nx*ny*nnz*(ib-1)*2
!OCL SERIAL
                   do k = nis, nie
                      iadd  = (k+(i-1)*nzp+(j-1)*nzp*nx)*2
                      iadd1 = iadd0+(k-nis+1+(i-1)*nnz+(j-1)*nnz*nx)*2
                      wk_send2(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send2(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
    end if

#ifdef __TIMER_DO__
  call timer_end(126)
#endif


#ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_fft_xz_world)
  call timer_sta(127)
#endif

    call MPI_ALLTOALL(wk_send2, z2x_smax*kimg*iesize, mpi_double_precision,   &
   &                  wk_recv2, z2x_rmax*kimg*iesize, mpi_double_precision,   &
   &                                                  mpi_fft_xz_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Direct_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 60, ierr)
     endif

#ifdef __TIMER_DO__
  call timer_sta(128)
#endif
#ifdef __TIMER_DO__
  call timer_end(128)
#endif

#else
!else ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_comm)
  call timer_sta(127)
#endif

    icnt_recv = 0
    lrank = myrank + 1
    if (lrank > (nmrank-1)) lrank = 0
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank-1)) lrank = 0
       if ((z2x_recv(1,lrank) /= 0)) then
          call mpi_irecv(wk_recv2(1,z2x_recv(2,lrank)), z2x_recv(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_r(icnt_recv), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Direct_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 60, ierr)
           endif
          icnt_recv = icnt_recv + 1
       endif
    enddo

    icnt_send = 0
    lrank = myrank
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank-1)) lrank = 0
       if ((z2x_send(1,lrank) /= 0)) then
          call mpi_isend(wk_send2(1,z2x_send(2,lrank)), z2x_send(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_s(icnt_send), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Direct_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 61, ierr)
           endif
          icnt_send = icnt_send + 1
       endif
    enddo

#ifdef __TIMER_DO__
  call timer_sta(128)
#endif
#ifdef __TIMER_DO__
  call timer_end(128)
#endif

    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Direct_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 62, ierr)
     endif

    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Direct_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 63, ierr)
     endif
#endif
!endif ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_end(127)
#endif

#ifdef __TIMER_DO__
  call timer_sta(129)
#endif
    wk_afft_l = 0.0d0

    nx = xyz_fft_x(2,1)-xyz_fft_x(1,1)+1
    ny = xyz_fft_x(2,2)-xyz_fft_x(1,2)+1
    nz = xyz_fft_x(2,3)-xyz_fft_x(1,3)+1
!.. nyp = wk_size/(nx*nz)
!.. nxp = wk_size/(ny*nz)
!  (z,x,y) -> (x,y,z)
    if (kimg == 1) then
       do ib = 1, iesize
          do k = 1, nz
             do j = 1, ny
                do l = 1, fft_Z_x_dim
                   nis = nis_fft_Z_x(l)
                   nie = nie_fft_Z_x(l)
                   nnx = nie - nis + 1
                   if(mod(nis,2)>0) then
                     nisx = nis/2 + 1
                   else
                     nisx = nis/2
                   endif
                   niex = nie/2
                   iadd0 = nnx/2*ny*nz*(ib-1)*2
                   do i = nisx, niex
                      iadd1 = iadd0 + (k+(i-nisx)*nz+(j-1)*nz*nnx/2)*2
!XX                   iadd  = (k+(i-1)*nz+(j-1)*nz*nx/2)*2
                      iadd  = (i+(j-1)*nx/2+(k-1)*nx/2*ny)*2
                      wk_afft_l(iadd-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end do
    else
       do ib = 1, iesize
          do k = 1, nz
             do j = 1, ny
                do l = 1, fft_Z_x_dim
                   nis = nis_fft_Z_x(l)
                   nie = nie_fft_Z_x(l)
                   nnx = nie - nis + 1
                   iadd0 = nnx*ny*nz*(ib-1)*2
                   do i = nis, nie
                      iadd1 = iadd0+(k+(i-nis)*nz+(j-1)*nz*nnx)*2
                      iadd  = (i+(j-1)*nx+(k-1)*nx*ny)*2
                      wk_afft_l(iadd-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end do
    end if
#ifdef __TIMER_DO__
  call timer_end(129)
#endif

!
! X-axis (y-z div)
!
#ifdef __TIMER_DO__
  call timer_sta(253)
#endif
    if(kimg==1) then
! (x,y,z) <- (z,x,y)
!XX    wk_afft_y = wk_afft_l
!XX    do ib = 1, iesize
!XX       do k = 1, nz
!XX          do j = 1, ny
!XX             do i = 1, nx/2
!XX                iadd0 = (k+(i-1)*nz+(j-1)*nz*nx/2)*2
!XX                iadd1 = (i+(j-1)*nx/2+(k-1)*nx/2*ny)*2
!XX                wk_afft_l(iadd1-1,ib) = wk_afft_y(iadd0-1,ib)
!XX                wk_afft_l(iadd1  ,ib) = wk_afft_y(iadd0  ,ib)
!XX             enddo
!XX          enddo
!XX       enddo
!XX    enddo
       nsize(1:3) = (/nx-2,ny,nz/)
       isin(1:3)  = (/-1,0,0/)
!OCL INDEPENDENT[ DVMRF2 ]
!OCL PARALLEL_STRONG
       do ib = 1, iesize
          call DVMRF2(wk_afft_l(1,ib),nsize,3,isin,-1,ierr)
       end do
    else
       nsize(1:3) = (/nx,ny,nz/)
       isin(1:3)  = (/1,0,0/)
!OCL INDEPENDENT[ DVMCF2 ]
!OCL PARALLEL_STRONG
       do ib = 1, iesize
          call DVMCF2(wk_afft_l(1,ib),nsize,3,isin,ierr)
       end do
    end if
#ifdef __TIMER_DO__
  call timer_end(253)
#endif

1000 continue

!F  if (allocated(req_r))       deallocate(req_r)
!F  if (allocated(req_s))       deallocate(req_s)
!F  if (allocated(sta_r))       deallocate(sta_r)
!F  if (allocated(sta_s))       deallocate(sta_s)
!F  if (allocated(wk_recvcnt))    deallocate(wk_recvcnt)
!F  if (allocated(wk_sendcnt))    deallocate(wk_sendcnt)
!F  if (allocated(wk_recvrank))   deallocate(wk_recvrank)
!F  if (allocated(wk_sendrank))   deallocate(wk_sendrank)
!F  if (allocated(wk_afft_y))   deallocate(wk_afft_y)
    if (allocated(wk_mp_fft_y)) deallocate(wk_mp_fft_y)

#ifdef __TIMER_SUB__
    call timer_end(105)
#endif
  end subroutine m_FFT_Direct_3D
!------------------------------------------------------------------------------
#endif
!ifndef FFT_USE_SSL2

#ifndef FFT_USE_SSL2

#ifdef FFT_FFTW_OLD

!------------------------------------------------------------------------------
  subroutine m_FFT_Inverse_3D(nfout,wk_afft_l, wk_size, iesize)
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
use mod_timer
#endif
! === TIMERTIMERTIMER ==========================================================
    integer, intent(in) :: nfout, wk_size, iesize
    real(kind=DP), dimension(wk_size*kimg,iesize), intent(inout) :: wk_afft_l

    real(kind=DP), allocatable, dimension(:) :: wk_allfft
    integer,      allocatable, dimension(:)   :: wk_recvdsp
    integer, parameter :: FFTW_MEASURE=0
    integer, parameter :: FFTW_ESTIMATE=64
    integer, parameter :: FFTW_RANK=1
    integer            :: FFTW_FLAG
    integer,dimension(1) :: NFFTW3, NEMBED, NEREAL
    integer :: nx, ny, nz, i, j, k, l, ri, ib, iadd, lrank
    integer :: irank, itag, icnt_send, icnt_recv
    integer,save, allocatable, dimension(:)   :: req_r, req_s
    integer,save, allocatable, dimension(:,:) :: sta_r, sta_s
    integer,save, allocatable, dimension(:)   :: wk_recvcnt, wk_sendcnt
!   real(kind=DP),allocatable, dimension(:,:) :: wk_recv, wk_send
    real(kind=DP),allocatable, dimension(:,:),save :: wk_recv1, wk_send1
    real(kind=DP),allocatable, dimension(:,:),save :: wk_recv2, wk_send2
    integer,      allocatable, dimension(:)   :: wk_mp_fft_y
    real(kind=DP),allocatable, dimension(:,:),save  :: wk_afft_y
    real(kind=DP),allocatable, dimension(:,:,:) :: wk_gather
    integer :: mpi_comm, myrank, nmrank
    integer(kind=8),save :: planx1 = 0, planx2 = 0
    integer(kind=8),save :: plany1 = 0, plany2 = 0
    integer(kind=8),save :: planz1 = 0, planz2 = 0
    integer(kind=8) :: plan
    integer :: max_x, max_y, max_z
    integer,save :: max_elm = 0, savesize = 0

    integer,save, allocatable, dimension(:,:) :: x2z_recv, x2z_send, z2y_recv, z2y_send
    integer,save :: x2z_rrank, x2z_srank, z2y_rrank, z2y_srank
    integer,save :: x2z_rmax, x2z_smax, x2z_srmax, z2y_rmax, z2y_smax, z2y_srmax

    integer ::  nnx, nny, nnz, iadd0, iadd1, nis, nie, jadd, kadd
#ifdef FFT_ALLTOALL
    integer,dimension(2) :: max_send,max_recv
#endif

#ifdef __TIMER_SUB__
    call timer_barrier(mpi_ke_world)
    call timer_sta(106)
#endif

    if (kimg == 1) then
       FFTW_FLAG = -1      ! FFTW_FORWARD
    else
       FFTW_FLAG = +1      ! FFTW_BACKWARD
    endif

! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
call start_timer('FFT-I reg0')
#endif
! === TIMERTIMERTIMER ==========================================================
!!  mpi_comm = mpi_kg_world
!!  myrank = myrank_e
!!  nmrank = nrank_e
    mpi_comm = mpi_ke_world
    myrank = myrank_g
    nmrank = nrank_g
    itag = 10

    if (np_fft_x==0)then
       nx = 0
       nx = 0
       nz = 0
!      go to 1000
       go to 2000
    endif

    if (firstcall_inverse_3d) then
       savesize = 0
#ifdef __TIMER_DO__
  call timer_sta(134)
#endif
       max_x = maxval(nel_fft_x(:))
       max_y = maxval(nel_fft_y(:))
       max_z = maxval(nel_fft_z(:))
       max_elm = max(max_x,max_y,max_z)
       if(allocated(req_r)) deallocate(req_r)
       if(allocated(req_s)) deallocate(req_s)
       if(allocated(sta_r)) deallocate(sta_r)
       if(allocated(sta_s)) deallocate(sta_s)
       if(allocated(wk_recvcnt)) deallocate(wk_recvcnt)
       if(allocated(wk_sendcnt)) deallocate(wk_sendcnt)
       allocate(req_r(0:nmrank-1), stat=ierr)
       allocate(req_s(0:nmrank-1), stat=ierr)
       allocate(sta_r(MPI_STATUS_SIZE,0:nmrank-1), stat=ierr)
       allocate(sta_s(MPI_STATUS_SIZE,0:nmrank-1), stat=ierr)
       allocate(wk_recvcnt(0:nmrank-1), stat=ierr)
       allocate(wk_sendcnt(0:nmrank-1), stat=ierr)
       if (ierr /= 0) then
          write(nfout,*)' m_FFT_Inverse_3D :  Not allocate '
          call flush(nfout)
          call mpi_abort(mpi_comm_world, 68, ierr)
       endif

       if(allocated(x2z_recv)) deallocate(x2z_recv)
       if(allocated(x2z_send)) deallocate(x2z_send)
       allocate(x2z_recv(2,0:nmrank-1))
       allocate(x2z_send(2,0:nmrank-1))
       x2z_send = 0
       x2z_recv = 0
       wk_recvcnt = 0
       wk_sendcnt = 0
       do i = 1, nel_fft_z(myrank)
          irank = map_fft_x(mp_fft_z(i)) - 1
          wk_recvcnt(irank) = wk_recvcnt(irank) + 1
       enddo
       do i = 1, nel_fft_x(myrank)
          irank = map_fft_z(mp_fft_x(i)) - 1
          wk_sendcnt(irank) = wk_sendcnt(irank) + 1
       enddo
       k = 0
       do i = 0, nmrank - 1
          if(wk_recvcnt(i) /= 0) then
             k = k + 1
             x2z_recv(1,i) = wk_recvcnt(i)
             x2z_recv(2,i) = k
          endif
       enddo
       x2z_rrank = k
       k = 0
       do i = 0, nmrank - 1
          if(wk_sendcnt(i) /= 0) then
             k = k + 1
             x2z_send(1,i) = wk_sendcnt(i)
             x2z_send(2,i) = k
          endif
       enddo
       x2z_srank = k
       x2z_rmax = maxval(wk_recvcnt)
       x2z_smax = maxval(wk_sendcnt)
       x2z_srmax = max(x2z_rmax,x2z_smax)

       if(allocated(z2y_recv)) deallocate(z2y_recv)
       if(allocated(z2y_send)) deallocate(z2y_send)
       allocate(z2y_recv(2,0:nmrank-1))
       allocate(z2y_send(2,0:nmrank-1))
       z2y_send = 0
       z2y_recv = 0
       wk_recvcnt = 0
       wk_sendcnt = 0

       nx = xyz_fft_y(2,1)-xyz_fft_y(1,1)+1
       ny = xyz_fft_y(2,2)-xyz_fft_y(1,2)+1
       nz = xyz_fft_y(2,3)-xyz_fft_y(1,3)+1
       allocate(wk_mp_fft_y(nx*ny*nz) ,stat=ierr)
       do k = 0, nz-1
          do j = 0, ny-1
             do i = 0, nx-1
                wk_mp_fft_y(i+j*nx+k*nx*ny+1) = mp_fft_y(i+k*nx+j*nx*nz+1)
             enddo
          enddo
       enddo

       do i = 1, nel_fft_y(myrank)
          irank = map_fft_z(wk_mp_fft_y(i)) - 1
          wk_recvcnt(irank) = wk_recvcnt(irank) + 1
       enddo
       do i = 1, nel_fft_z(myrank)
          irank = map_fft_y(mp_fft_z(i)) - 1
          wk_sendcnt(irank) = wk_sendcnt(irank) + 1
       enddo
       k = 0
       do i = 0, nmrank - 1
          if(wk_recvcnt(i) /= 0) then
             k = k + 1
             z2y_recv(1,i) = wk_recvcnt(i)
             z2y_recv(2,i) = k
          endif
       enddo
       z2y_rrank = k
       k = 0
      do i = 0, nmrank - 1
          if(wk_sendcnt(i) /= 0) then
             k = k + 1
             z2y_send(1,i) = wk_sendcnt(i)
             z2y_send(2,i) = k
          endif
       enddo
       z2y_srank = k
       z2y_rmax = maxval(wk_recvcnt)
       z2y_smax = maxval(wk_sendcnt)
       z2y_srmax = max(z2y_rmax,z2y_smax)

       deallocate(wk_mp_fft_y)

#ifdef FFT_ALLTOALL
       max_send(1) = x2z_rmax
       max_send(2) = x2z_smax
       call mpi_allreduce(max_send,max_recv,2,mpi_integer,mpi_max,mpi_fft_xz_world,ierr)
       x2z_rmax = max_recv(1)
       x2z_smax = max_recv(2)
       max_send(1) = z2y_rmax
       max_send(2) = z2y_smax
       call mpi_allreduce(max_send,max_recv,2,mpi_integer,mpi_max,mpi_fft_zy_world,ierr)
       z2y_rmax = max_recv(1)
       z2y_smax = max_recv(2)
#endif

       if (ipri > 1) then
          write(nfout,'("m_FFT_Inverse_3D   --   myrank_g=",i4)') myrank_g
          write(nfout,'("x2z_send")')
          write(nfout,'(10(i8,", "))') (x2z_send(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (x2z_send(2,i),i=0,nmrank-1)
          write(nfout,'("x2z_recv")')
          write(nfout,'(10(i8,", "))') (x2z_recv(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (x2z_recv(2,i),i=0,nmrank-1)
          write(nfout,'("z2y_send")')
          write(nfout,'(10(i8,", "))') (z2y_send(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (z2y_send(2,i),i=0,nmrank-1)
          write(nfout,'("z2y_recv")')
          write(nfout,'(10(i8,", "))') (z2y_recv(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (z2y_recv(2,i),i=0,nmrank-1)
          call flush(nfout)
       endif

#ifdef __TIMER_DO__
  call timer_sta(135)
#endif
       nx = xyz_fft_x(2,1)-xyz_fft_x(1,1)+1
       ny = xyz_fft_x(2,2)-xyz_fft_x(1,2)+1
       nz = xyz_fft_x(2,3)-xyz_fft_x(1,3)+1
       if(kimg==1) then
          NFFTW3(1) = nx - 2
          NEREAL(1) = nx - 2
          NEMBED(1) = (nx - 2)/2
          call dfftw_plan_many_dft_r2c(planx1, FFTW_RANK, NFFTW3, ny*nz,   &
         &                             wk_afft_l, NEREAL, 1, nx,          &
         &                             wk_afft_l, NEMBED, 1, nx/2,        &
         &                             FFTW_MEASURE )
       else
          NFFTW3(1) = nx
          NEMBED(1) = nx
          call dfftw_plan_many_dft    (planx2, FFTW_RANK, NFFTW3, ny*nz,   &
         &                             wk_afft_l, NEMBED, 1, nx,          &
         &                             wk_afft_l, NEMBED, 1, nx,          &
         &                             FFTW_FLAG, FFTW_MEASURE )
       endif

       nx = xyz_fft_z(2,1)-xyz_fft_z(1,1)+1
       ny = xyz_fft_z(2,2)-xyz_fft_z(1,2)+1
       nz = xyz_fft_z(2,3)-xyz_fft_z(1,3)+1
       if(kimg==1) then
          NFFTW3(1) = nz
          NEMBED(1) = nz
          call dfftw_plan_many_dft    (planz1, FFTW_RANK, NFFTW3, nx*ny/2, &
         &                             wk_afft_l, NEMBED, nx*ny/2, 1,     &
         &                             wk_afft_l, NEMBED, nx*ny/2, 1,     &
         &                             FFTW_FLAG, FFTW_MEASURE )
       else
          NFFTW3(1) = nz
          NEMBED(1) = nz
          call dfftw_plan_many_dft    (planz2, FFTW_RANK, NFFTW3, nx*ny,   &
         &                             wk_afft_l, NEMBED, nx*ny, 1,       &
         &                             wk_afft_l, NEMBED, nx*ny, 1,       &
         &                             FFTW_FLAG, FFTW_MEASURE )
       endif

       nx = xyz_fft_y(2,1)-xyz_fft_y(1,1)+1
       ny = xyz_fft_y(2,2)-xyz_fft_y(1,2)+1
       nz = xyz_fft_y(2,3)-xyz_fft_y(1,3)+1
       if(kimg==1) then
          NFFTW3(1) = ny
          NEMBED(1) = ny
          call dfftw_plan_many_dft    (plany1, FFTW_RANK, NFFTW3, nz*nx/2, &
         &                             wk_afft_l, NEMBED, nz*nx/2, 1,    &
         &                             wk_afft_l, NEMBED, nz*nx/2, 1,    &
         &                             FFTW_FLAG, FFTW_MEASURE )
       else
          NFFTW3(1) = ny
          NEMBED(1) = ny
          call dfftw_plan_many_dft    (plany2, FFTW_RANK, NFFTW3, nz*nx,   &
         &                             wk_afft_l, NEMBED, nz*nx, 1,      &
         &                             wk_afft_l, NEMBED, nz*nx, 1,      &
         &                             FFTW_FLAG, FFTW_MEASURE )
       end if
#ifdef __TIMER_DO__
  call timer_end(135)
#endif

       firstcall_inverse_3d = .false.
       go to 1000
#ifdef __TIMER_DO__
  call timer_end(134)
#endif
    endif

#ifdef FFT_ALLTOALL

    if (iesize /= savesize) then
#ifdef __TIMER_ETC__
  call timer_sta(264)
#endif
       if (allocated(wk_recv1)) deallocate(wk_recv1)
       if (allocated(wk_send1)) deallocate(wk_send1)
       if (allocated(wk_recv2)) deallocate(wk_recv2)
       if (allocated(wk_send2)) deallocate(wk_send2)
       if (allocated(wk_afft_y)) deallocate(wk_afft_y)
       allocate(wk_recv1(x2z_rmax*kimg*iesize,x2z_rrank), stat=ierr)
       allocate(wk_send1(x2z_smax*kimg*iesize,x2z_srank), stat=ierr)
       allocate(wk_recv2(z2y_rmax*kimg*iesize,z2y_rrank), stat=ierr)
       allocate(wk_send2(z2y_smax*kimg*iesize,z2y_srank), stat=ierr)
       allocate(wk_afft_y(max_elm*kimg,iesize) ,stat=ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_FFT_Inverse_3D :  Not allocate '
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 69, ierr)
        endif
       savesize = iesize
#ifdef __TIMER_ETC__
  call timer_end(264)
#endif
    end if

#else
!else ifdef FFT_ALLTOALL

    if (iesize > savesize) then
#ifdef __TIMER_ETC__
  call timer_sta(264)
#endif
       if (allocated(wk_recv1)) deallocate(wk_recv1)
       if (allocated(wk_send1)) deallocate(wk_send1)
       if (allocated(wk_recv2)) deallocate(wk_recv2)
       if (allocated(wk_send2)) deallocate(wk_send2)
       if (allocated(wk_afft_y)) deallocate(wk_afft_y)
       allocate(wk_recv1(x2z_rmax*kimg*iesize,x2z_rrank), stat=ierr)
       allocate(wk_send1(x2z_smax*kimg*iesize,x2z_srank), stat=ierr)
       allocate(wk_recv2(z2y_rmax*kimg*iesize,z2y_rrank), stat=ierr)
       allocate(wk_send2(z2y_smax*kimg*iesize,z2y_srank), stat=ierr)
       allocate(wk_afft_y(max_elm*kimg,iesize) ,stat=ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_FFT_Inverse_3D :  Not allocate '
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 69, ierr)
        endif
       savesize = iesize
#ifdef __TIMER_ETC__
  call timer_end(264)
#endif
    end if

#endif
!endif ifdef FFT_ALLTOALL

 2000 continue
#ifdef __TIMER_SUB__
    call mpi_barrier(mpi_ke_world, ierr)
    call timer_sta(1498)
#endif
    if (np_fft_x==0)then
       nx = 0
       nx = 0
       nz = 0
       go to 1000
    endif

! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
call stop_timer('FFT-I reg0')
#endif
! === TIMERTIMERTIMER ==========================================================
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP          SHARED(     nfout,wk_afft_l,iesize,req_r,req_s,sta_r,sta_s,    &
!$OMP                 wk_recv1,wk_send1,wk_recv2,wk_send2,wk_afft_y,          &
!$OMP                 planx1,planx2,planz1,planz2,plany1,plany2,              &
!$OMP                 x2z_recv, x2z_send, z2y_recv, z2y_send,itag,            &
!$OMP                 xyz_fft_y,xyz_fft_z,xyz_fft_x,fft_Z_x_dim,nis_fft_Z_x,  &
!$OMP                 nie_fft_Z_x,fft_X_z_dim,nis_fft_X_z,nie_fft_X_z,        &
!$OMP                 fft_Y_z_dim,nis_fft_Y_z,nie_fft_Y_z,fft_Z_y_dim,        &
!$OMP                 nie_fft_Z_y,nis_fft_Z_y,                                &
#ifdef FFT_ALLTOALL
!$OMP                 mpi_fft_xz_world, x2z_rmax, x2z_smax,                   &
!$OMP                 mpi_fft_zy_world, z2y_rmax, z2y_smax,                   &
#endif
!$OMP                 mpi_comm,myrank,nmrank,kimg,ierr                    )   &
!$OMP          PRIVATE(nx,ny,nz,i,j,k,l,ri,ib,iadd,lrank,icnt_send,icnt_recv, &
!$OMP                  irank,nnx,nny,nnz,iadd0,iadd1,nis,nie,jadd,kadd,plan   )

#ifdef __TIMER_SUB__
    call timer_sta(1499)
#endif
!
! X-axis (y-z div)
!
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call start_timer('FFT-I xfft')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================
    nx = xyz_fft_x(2,1)-xyz_fft_x(1,1)+1
    ny = xyz_fft_x(2,2)-xyz_fft_x(1,2)+1
    nz = xyz_fft_x(2,3)-xyz_fft_x(1,3)+1

    if (kimg == 1) then
       plan = planx1
    else
       plan = planx2
    endif

#ifdef __TIMER_DO__
  call timer_sta(254)
#endif
!$OMP DO
    do ib = 1, iesize
       if (kimg==1) then
          call dfftw_execute_dft_r2c(plan,wk_afft_l(1,ib),wk_afft_l(1,ib))
       else
          call dfftw_execute_dft    (plan,wk_afft_l(1,ib),wk_afft_l(1,ib))
       endif
    enddo

#ifdef __TIMER_DO__
  call timer_end(254)
#endif
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call stop_timer('FFT-I xfft')
call start_timer('FFT-I x2z')
call start_timer('*FFT-I x2z sbuf')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================

#ifdef __TIMER_DO__
  call timer_sta(136)
#endif
#ifdef FFT_INTERCHANGE
    if (kimg == 1) then
!$OMP DO
       do ib = 1, iesize
          do k = 1, nz
             do j = 1, ny
                do l = 1, fft_Z_x_dim
                   nis = nis_fft_Z_x(l)
                   nie = nie_fft_Z_x(l)
                   nnx = nie - nis + 1
                   iadd0 = nnx*ny*nz*(ib-1)
                   do i = nis, nie
                      wk_send1(iadd0+(i-nis+1+(j-1)*nnx+(k-1)*nnx*ny),l) = wk_afft_l(i+(j-1)*nx+(k-1)*nx*ny,ib)
                   end do
                end do
             end do
          end do
       end do
    else
!$OMP DO
       do ib = 1, iesize
          do k = 1, nz
             do j = 1, ny
                do l = 1, fft_Z_x_dim
                   nis = nis_fft_Z_x(l)
                   nie = nie_fft_Z_x(l)
                   nnx = nie - nis + 1
                   do i = nis, nie
                      iadd1 = nnx*ny*nz*(ib-1)*2+(i-nis+1+(j-1)*nnx+(k-1)*nnx*ny)*2
                      iadd  = (i+(j-1)*nx+(k-1)*nx*ny)*2
                      wk_send1(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send1(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
    end if
#else
!$OMP DO
    do l = 1, fft_Z_x_dim
       nis = nis_fft_Z_x(l)
       nie = nie_fft_Z_x(l)
       nnx = nie - nis + 1
       if (kimg == 1) then
          do ib = 1, iesize
             iadd0 = nnx*ny*nz*(ib-1)
             do k = 1, nz
                do j = 1, ny
                   do i = nis, nie
                      wk_send1(iadd0+(i-nis+1+(j-1)*nnx+(k-1)*nnx*ny),l) = wk_afft_l(i+(j-1)*nx+(k-1)*nx*ny,ib)
                   end do
                end do
             end do
          end do
       else
          do ib = 1, iesize
             iadd0 = nnx*ny*nz*(ib-1)*2
             do k = 1, nz
                jadd = (k-1)*nnx*ny
                do j = 1, ny
                   kadd = (j-1)*nnx+jadd
                   do i = nis, nie
                      iadd1 = iadd0+(i-nis+1+kadd)*2
                      iadd  = (i+(j-1)*nx+(k-1)*nx*ny)*2
!fj --------------------
!                     do ri = 1, kimg
!                       wk_send1(iadd1+(ri-kimg),l) = wk_afft_l(iadd+(ri-kimg),ib)
!                     end do
                      wk_send1(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send1(iadd1  ,l) = wk_afft_l(iadd  ,ib)
!fj --------------------
                   end do
                end do
             end do
          end do
       end if
    end do
#endif
#ifdef __TIMER_DO__
  call timer_end(136)
#endif
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call stop_timer('*FFT-I x2z sbuf')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================

#ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_fft_xz_world)
  call timer_sta(137)
#endif

! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call start_timer('*FFT-I x2z a2a')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================
!$OMP MASTER

    call MPI_ALLTOALL(wk_send1, x2z_smax*kimg*iesize, mpi_double_precision,   &
   &                  wk_recv1, x2z_rmax*kimg*iesize, mpi_double_precision,   &
   &                                                  mpi_fft_xz_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Inverse_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 70, ierr)
     endif

!$OMP END MASTER
!$OMP BARRIER
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call stop_timer('*FFT-I x2z a2a')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================

#ifdef __TIMER_DO__
  call timer_sta(138)
#endif
#ifdef __TIMER_DO__
  call timer_end(138)
#endif

#else
!else ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_comm)
  call timer_sta(137)
#endif

! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call start_timer('*FFT-I x2z sr')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================
!$OMP MASTER
    icnt_recv = 0
    lrank = myrank + 1
    if (lrank > (nmrank -1)) lrank = 0
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank -1)) lrank = 0
       if ((lrank /= myrank) .and. (x2z_recv(1,lrank) /= 0)) then
          call mpi_irecv(wk_recv1(1,x2z_recv(2,lrank)), x2z_recv(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_r(icnt_recv), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Inverse_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 70, ierr)
           endif
          icnt_recv = icnt_recv + 1
       endif
    enddo

    icnt_send = 0
    lrank = myrank
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank -1)) lrank = 0
       if ((lrank /= myrank) .and. (x2z_send(1,lrank) /= 0)) then
          call mpi_isend(wk_send1(1,x2z_send(2,lrank)), x2z_send(1,lrank)*kimg*iesize, &
                         mpi_double_precision, lrank, itag, mpi_comm, req_s(icnt_send), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Inverse_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 71, ierr)
           endif
          icnt_send = icnt_send + 1
       endif
    enddo
!$OMP END MASTER
!$OMP BARRIER

#ifdef __TIMER_DO__
  call timer_sta(138)
#endif
!$OMP DO
    do i = 1, x2z_recv(1,myrank)*kimg*iesize
       wk_recv1(i,x2z_recv(2,myrank)) = wk_send1(i,x2z_send(2,myrank))
    enddo
#ifdef __TIMER_DO__
  call timer_end(138)
#endif

!$OMP MASTER
    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Inverse_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 72, ierr)
     endif

    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Inverse_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 73, ierr)
     endif
!$OMP END MASTER
!$OMP BARRIER
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call stop_timer('*FFT-I x2z sr')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================

#endif
!endif ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_end(137)
#endif

#ifdef __TIMER_DO__
  call timer_sta(139)
#endif
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call start_timer('*FFT-I x2z rbuf')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================
    nx = xyz_fft_z(2,1)-xyz_fft_z(1,1)+1
    ny = xyz_fft_z(2,2)-xyz_fft_z(1,2)+1
    nz = xyz_fft_z(2,3)-xyz_fft_z(1,3)+1
!$OMP DO
    do l = 1, fft_X_z_dim
       nis = nis_fft_X_z(l)
       nie = nie_fft_X_z(l)
       nnz = nie - nis + 1
       if (kimg == 1) then
          do ib = 1, iesize
             iadd0 = nx*ny*nnz*(ib-1)
             do k = nis, nie
                do j = 1, ny
                   do i = 1, nx
                      wk_afft_l(i+(j-1)*nx+(k-1)*nx*ny,ib) = wk_recv1(iadd0+i+(j-1)*nx+(k-nis)*nx*ny,l)
                   end do
                end do
             end do
          end do
       else
          do ib = 1, iesize
             iadd0 = nx*ny*nnz*(ib-1)*2
             do k = nis, nie
                do j = 1, ny
                   do i = 1, nx
                      iadd1 = iadd0+(i+(j-1)*nx+(k-nis)*nx*ny)*2
                      iadd  = (i+(j-1)*nx+(k-1)*nx*ny)*2
                      wk_afft_l(iadd-1,ib) = wk_recv1(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv1(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end if
    end do
#ifdef __TIMER_DO__
  call timer_end(139)
#endif
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call stop_timer('*FFT-I x2z rbuf')
call stop_timer('FFT-I x2z')
call start_timer('FFT-I zfft')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================

!
! Z-axis (x-y div)
!
    nx = xyz_fft_z(2,1)-xyz_fft_z(1,1)+1
    ny = xyz_fft_z(2,2)-xyz_fft_z(1,2)+1
    nz = xyz_fft_z(2,3)-xyz_fft_z(1,3)+1

    if (kimg == 1) then
       plan = planz1
    else
       plan = planz2
    endif
#ifdef __TIMER_DO__
  call timer_sta(255)
#endif
!$OMP DO
    do ib = 1, iesize
       call dfftw_execute_dft(plan,wk_afft_l(1,ib),wk_afft_l(1,ib))
    end do
#ifdef __TIMER_DO__
  call timer_end(255)
#endif

! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call stop_timer('FFT-I zfft')
call start_timer('FFT-I z2y')
call start_timer('*FFT-I z2y sbuf')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_DO__
  call timer_sta(140)
#endif
#ifdef __TIMER_DO__
  call timer_end(140)
#endif

#ifdef __TIMER_DO__
  call timer_sta(141)
#endif
    nx = xyz_fft_z(2,1)-xyz_fft_z(1,1)+1
    ny = xyz_fft_z(2,2)-xyz_fft_z(1,2)+1
    nz = xyz_fft_z(2,3)-xyz_fft_z(1,3)+1
!$OMP DO
    do l = 1, fft_Y_z_dim
       nis = nis_fft_Y_z(l)
       nie = nie_fft_Y_z(l)
       nnz = nie - nis + 1
       if (kimg == 1) then
          do ib = 1, iesize
             iadd0 = nx*ny*nnz*(ib-1)
             do k = nis, nie
                do j = 1, ny
                   do i = 1, nx
                      wk_send2(iadd0+(i+(j-1)*nx+(k-nis)*nx*ny),l) = wk_afft_l(i+(j-1)*nx+(k-1)*nx*ny,ib)
                   end do
                end do
             end do
          end do
       else
          do ib = 1, iesize
             iadd0 = nx*ny*nnz*(ib-1)*2
             do k = nis, nie
                jadd = (k-nis)*nx*ny
                do j = 1, ny
                   kadd = (j-1)*nx+jadd
                   do i = 1, nx
                      iadd1 = iadd0+(i+kadd)*2
                      iadd  = (i+(j-1)*nx+(k-1)*nx*ny)*2
!fj --------------------
!                     do ri = 1, kimg
!                        wk_send2(iadd1+(ri-kimg),l) = wk_afft_l(iadd+(ri-kimg),ib)
!                     end do
                      wk_send2(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send2(iadd1  ,l) = wk_afft_l(iadd  ,ib)
!fj --------------------
                   end do
                end do
             end do
          end do
       end if
    end do
#ifdef __TIMER_DO__
  call timer_end(141)
#endif

! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call stop_timer('*FFT-I z2y sbuf')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================
#ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_fft_zy_world)
  call timer_sta(142)
#endif

! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call start_timer('*FFT-I z2y a2a')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================
!$OMP MASTER

    call MPI_ALLTOALL(wk_send2, z2y_smax*kimg*iesize, mpi_double_precision,   &
   &                  wk_recv2, z2y_rmax*kimg*iesize, mpi_double_precision,   &
   &                                                  mpi_fft_zy_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Inverse_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 76, ierr)
     endif

!$OMP END MASTER
!$OMP BARRIER
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call stop_timer('*FFT-I z2y a2a')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================

#ifdef __TIMER_DO__
  call timer_sta(143)
#endif
#ifdef __TIMER_DO__
  call timer_end(143)
#endif

#else
!else ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_comm)
  call timer_sta(142)
#endif

! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call start_timer('*FFT-I z2y sr')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================
!$OMP MASTER
    icnt_recv = 0
    lrank = myrank + 1
    if (lrank > (nmrank-1)) lrank = 0
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank-1)) lrank = 0
       if ((lrank /= myrank) .and. (z2y_recv(1,lrank) /= 0)) then
          call mpi_irecv(wk_recv2(1,z2y_recv(2,lrank)), z2y_recv(1,lrank)*kimg*iesize, &
                         mpi_double_precision, lrank, itag, mpi_comm, req_r(icnt_recv), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Inverse_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 76, ierr)
           endif
          icnt_recv = icnt_recv + 1
       endif
    enddo

    icnt_send = 0
    lrank = myrank
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank -1)) lrank = 0
       if ((lrank /= myrank) .and. (z2y_send(1,lrank) /= 0)) then
          call mpi_isend(wk_send2(1,z2y_send(2,lrank)), z2y_send(1,lrank)*kimg*iesize, &
                         mpi_double_precision, lrank, itag, mpi_comm, req_s(icnt_send), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Inverse_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 77, ierr)
           endif
          icnt_send = icnt_send + 1
       endif
    enddo
!$OMP END MASTER
!$OMP BARRIER
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call stop_timer('*FFT-I z2y sr')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================

#ifdef __TIMER_DO__
  call timer_sta(143)
#endif
!$OMP DO
    do i = 1, z2y_recv(1,myrank)*kimg*iesize
       wk_recv2(i,z2y_recv(2,myrank)) = wk_send2(i,z2y_send(2,myrank))
    enddo
#ifdef __TIMER_DO__
  call timer_end(143)
#endif

!$OMP MASTER
    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Inverse_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 78, ierr)
     endif

    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Inverse_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 79, ierr)
     endif
!$OMP END MASTER
!$OMP BARRIER

#endif
!endif ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_end(142)
#endif

#ifdef __TIMER_DO__
  call timer_sta(144)
#endif
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call start_timer('*FFT-I z2y rbuf')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================
    nx = xyz_fft_y(2,1)-xyz_fft_y(1,1)+1
    ny = xyz_fft_y(2,2)-xyz_fft_y(1,2)+1
    nz = xyz_fft_y(2,3)-xyz_fft_y(1,3)+1
#ifdef FFT_INTERCHANGE
    if (kimg == 1) then
!$OMP DO
       do ib = 1, iesize
          do k = 1, nz
             do l = 1, fft_Z_y_dim
                nis = nis_fft_Z_y(l)
                nie = nie_fft_Z_y(l)
                nny = nie - nis + 1
                iadd0 = nx*nny*nz*(ib-1)
                do j = nis, nie
                   do i = 1, nx
                      wk_afft_y(i+(j-1)*nx+(k-1)*nx*ny,ib) = wk_recv2(iadd0+i+(j-nis)*nx+(k-1)*nx*nny,l)
                   end do
                end do
             end do
          end do
       end do
    else
!$OMP DO
       do ib = 1, iesize
          do k = 1, nz
             do l = 1, fft_Z_y_dim
                nis = nis_fft_Z_y(l)
                nie = nie_fft_Z_y(l)
                nny = nie - nis + 1
                do j = nis, nie
                   do i = 1, nx
                      iadd1 = nx*nny*nz*(ib-1)*2+(i+(j-nis)*nx+(k-1)*nx*nny)*2
                      iadd  = (i+(j-1)*nx+(k-1)*nx*ny)*2
                      wk_afft_y(iadd-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_y(iadd  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end do
    end if
#else
!$OMP DO
    do l = 1, fft_Z_y_dim
       nis = nis_fft_Z_y(l)
       nie = nie_fft_Z_y(l)
       nny = nie - nis + 1
       if (kimg == 1) then
          do ib = 1, iesize
             iadd0 = nx*nny*nz*(ib-1)
             do k = 1, nz
                do j = nis, nie
                   do i = 1, nx
                      wk_afft_y(i+(j-1)*nx+(k-1)*nx*ny,ib) = wk_recv2(iadd0+i+(j-nis)*nx+(k-1)*nx*nny,l)
                   end do
                end do
             end do
          end do
       else
          do ib = 1, iesize
             iadd0 = nx*nny*nz*(ib-1)*2
             do k = 1, nz
                do j = nis, nie
                   do i = 1, nx
                      iadd1 = iadd0+(i+(j-nis)*nx+(k-1)*nx*nny)*2
                      iadd  = (i+(j-1)*nx+(k-1)*nx*ny)*2
                      wk_afft_y(iadd-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_y(iadd  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end if
    end do
#endif
#ifdef __TIMER_DO__
  call timer_end(144)
#endif

#ifdef __TIMER_DO__
  call timer_sta(145)
#endif
!fj --------------------
!   do ib = 1, iesize
!      do k = 0, nz-1
!         do j = 0, ny-1
!            do i = 0, nx-1
!   do ri = 1, kimg
!               wk_afft_l((i+k*nx+j*nx*nz+1)*kimg+(ri-kimg),ib) = wk_afft_y((i+j*nx+k*nx*ny+1)*kimg+(ri-kimg),ib)
!   enddo
!            enddo
!         enddo
!      enddo
!   enddo
#ifdef FFT_INTERCHANGE
    if (kimg == 1) then
       do ib = 1, iesize
!$OMP DO
          do k = 0, nz-1
             do j = 0, ny-1
                do i = 0, nx-1
                   wk_afft_l(i+k*nx+j*nx*nz+1,ib) = wk_afft_y(i+j*nx+k*nx*ny+1,ib)
                enddo
             enddo
          enddo
       enddo
    else
       do ib = 1, iesize
!$OMP DO
          do k = 0, nz-1
             do j = 0, ny-1
                do i = 0, nx-1
                   wk_afft_l((i+k*nx+j*nx*nz+1)*2-1,ib) = wk_afft_y((i+j*nx+k*nx*ny+1)*2-1,ib)
                   wk_afft_l((i+k*nx+j*nx*nz+1)*2  ,ib) = wk_afft_y((i+j*nx+k*nx*ny+1)*2  ,ib)
                enddo
             enddo
          enddo
       enddo
    end if
#else
    if (kimg == 1) then
!$OMP DO
       do ib = 1, iesize
          do k = 0, nz-1
             do j = 0, ny-1
                do i = 0, nx-1
                   wk_afft_l(i+k*nx+j*nx*nz+1,ib) = wk_afft_y(i+j*nx+k*nx*ny+1,ib)
                enddo
             enddo
          enddo
       enddo
    else
!$OMP DO
       do ib = 1, iesize
          do k = 0, nz-1
             do j = 0, ny-1
                do i = 0, nx-1
                   wk_afft_l((i+k*nx+j*nx*nz+1)*2-1,ib) = wk_afft_y((i+j*nx+k*nx*ny+1)*2-1,ib)
                   wk_afft_l((i+k*nx+j*nx*nz+1)*2  ,ib) = wk_afft_y((i+j*nx+k*nx*ny+1)*2  ,ib)
                enddo
             enddo
          enddo
       enddo
    end if
#endif
!fj --------------------
#ifdef __TIMER_DO__
  call timer_end(145)
#endif
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call stop_timer('*FFT-I z2y rbuf')
call stop_timer('FFT-I z2y')
call start_timer('FFT-I yfft')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================

!
! Y-axis (z-x div)
!
    if (kimg == 1) then
       plan = plany1
    else
       plan = plany2
    endif
#ifdef __TIMER_DO__
  call timer_sta(256)
#endif
!$OMP DO
    do ib = 1,iesize
       call dfftw_execute_dft(plan,wk_afft_l(1,ib),wk_afft_l(1,ib))
    end do
#ifdef __TIMER_DO__
  call timer_end(256)
#endif
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call stop_timer('FFT-I yfft')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================

#ifdef __TIMER_SUB__
  call timer_end(1498)
#endif

!$OMP END PARALLEL

#ifdef __TIMER_SUB__
  call timer_end(1499)
#endif

1000 continue

!F  if (allocated(req_r))       deallocate(req_r)
!F  if (allocated(req_s))       deallocate(req_s)
!F  if (allocated(sta_r))       deallocate(sta_r)
!F  if (allocated(sta_s))       deallocate(sta_s)
!F  if (allocated(wk_recvcnt))  deallocate(wk_recvcnt)
!F  if (allocated(wk_sendcnt))  deallocate(wk_sendcnt)
    if (allocated(wk_mp_fft_y)) deallocate(wk_mp_fft_y)
!F  if (allocated(wk_afft_y))   deallocate(wk_afft_y)

#ifdef __TIMER_SUB__
    call timer_end(106)
#endif
  end subroutine m_FFT_Inverse_3D

#else
!else ifdef FFT_FFTW_OLD

  subroutine m_FFT_Inverse_3D(nfout,wk_afft_l, wk_size, iesize)
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
use mod_timer
#endif
! === TIMERTIMERTIMER ==========================================================
    integer, intent(in) :: nfout, wk_size, iesize
    real(kind=DP), dimension(wk_size*kimg,iesize), intent(inout) :: wk_afft_l

    real(kind=DP), allocatable, dimension(:) :: wk_allfft
    integer,      allocatable, dimension(:)   :: wk_recvdsp
    integer, parameter :: FFTW_MEASURE=0
    integer, parameter :: FFTW_ESTIMATE=64
    integer, parameter :: FFTW_RANK=1
    integer            :: FFTW_FLAG
    integer,dimension(1) :: NFFTW3, NEMBED, NEREAL
    integer :: nx, ny, nz, i, j, k, l, ri, ib, iadd, lrank
    integer :: irank, itag, icnt_send, icnt_recv
    integer,save, allocatable, dimension(:)   :: req_r, req_s
    integer,save, allocatable, dimension(:,:) :: sta_r, sta_s
    integer,save, allocatable, dimension(:)   :: wk_recvcnt, wk_sendcnt
!   real(kind=DP),allocatable, dimension(:,:) :: wk_recv, wk_send
    real(kind=DP),allocatable, dimension(:,:),save :: wk_recv1, wk_send1
    real(kind=DP),allocatable, dimension(:,:),save :: wk_recv2, wk_send2
    integer,      allocatable, dimension(:)   :: wk_mp_fft_y
    real(kind=DP),allocatable, dimension(:,:),save  :: wk_afft_y
    real(kind=DP),allocatable, dimension(:,:,:) :: wk_gather
    integer :: mpi_comm, myrank, nmrank
    integer(kind=8),save :: planx1 = 0, planx2 = 0
    integer(kind=8),save :: plany1 = 0, plany2 = 0
    integer(kind=8),save :: planz1 = 0, planz2 = 0
    integer(kind=8) :: plan
    integer :: max_x, max_y, max_z
    integer,save :: max_elm = 0, savesize = 0

    integer,save, allocatable, dimension(:,:) :: x2z_recv, x2z_send, z2y_recv, z2y_send
    integer,save :: x2z_rrank, x2z_srank, z2y_rrank, z2y_srank
    integer,save :: x2z_rmax, x2z_smax, x2z_srmax, z2y_rmax, z2y_smax, z2y_srmax

    integer ::  nnx, nny, nnz, iadd0, iadd1, nis, nie, jadd, kadd
#ifdef FFT_ALLTOALL
    integer,dimension(2) :: max_send,max_recv
#endif

#ifdef __TIMER_SUB__
    call timer_barrier(mpi_ke_world)
    call timer_sta(106)
#endif

    if (kimg == 1) then
       FFTW_FLAG = -1      ! FFTW_FORWARD
    else
       FFTW_FLAG = +1      ! FFTW_BACKWARD
    endif

! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
call start_timer('FFT-I reg0')
#endif
! === TIMERTIMERTIMER ==========================================================
!!  mpi_comm = mpi_kg_world
!!  myrank = myrank_e
!!  nmrank = nrank_e
    mpi_comm = mpi_ke_world
    myrank = myrank_g
    nmrank = nrank_g
    itag = 10

    if (np_fft_x==0)then
       nx = 0
       nx = 0
       nz = 0
!      go to 1000
       go to 2000
    endif

    if (firstcall_inverse_3d) then
#ifdef __FAPP__
    call fapp_start('fftwf_inverse_firstcall',1,1)
#endif
       savesize = 0
#ifdef __TIMER_DO__
  call timer_sta(134)
#endif
       max_x = maxval(nel_fft_x(:))
       max_y = maxval(nel_fft_y(:))
       max_z = maxval(nel_fft_z(:))
       max_elm = max(max_x,max_y,max_z)
       if(allocated(req_r)) deallocate(req_r)
       if(allocated(req_s)) deallocate(req_s)
       if(allocated(sta_r)) deallocate(sta_r)
       if(allocated(sta_s)) deallocate(sta_s)
       if(allocated(wk_recvcnt)) deallocate(wk_recvcnt)
       if(allocated(wk_sendcnt)) deallocate(wk_sendcnt)
       allocate(req_r(0:nmrank-1), stat=ierr)
       allocate(req_s(0:nmrank-1), stat=ierr)
       allocate(sta_r(MPI_STATUS_SIZE,0:nmrank-1), stat=ierr)
       allocate(sta_s(MPI_STATUS_SIZE,0:nmrank-1), stat=ierr)
       allocate(wk_recvcnt(0:nmrank-1), stat=ierr)
       allocate(wk_sendcnt(0:nmrank-1), stat=ierr)
       if (ierr /= 0) then
          write(nfout,*)' m_FFT_Inverse_3D :  Not allocate '
          call flush(nfout)
          call mpi_abort(mpi_comm_world, 68, ierr)
       endif

       if(allocated(x2z_recv)) deallocate(x2z_recv)
       if(allocated(x2z_send)) deallocate(x2z_send)
       allocate(x2z_recv(2,0:nmrank-1))
       allocate(x2z_send(2,0:nmrank-1))
       x2z_send = 0
       x2z_recv = 0
       wk_recvcnt = 0
       wk_sendcnt = 0
       do i = 1, nel_fft_z(myrank)
          irank = map_fft_x(mp_fft_z(i)) - 1
          wk_recvcnt(irank) = wk_recvcnt(irank) + 1
       enddo
       do i = 1, nel_fft_x(myrank)
          irank = map_fft_z(mp_fft_x(i)) - 1
          wk_sendcnt(irank) = wk_sendcnt(irank) + 1
       enddo
       k = 0
       do i = 0, nmrank - 1
          if(wk_recvcnt(i) /= 0) then
             k = k + 1
             x2z_recv(1,i) = wk_recvcnt(i)
             x2z_recv(2,i) = k
          endif
       enddo
       x2z_rrank = k
       k = 0
       do i = 0, nmrank - 1
          if(wk_sendcnt(i) /= 0) then
             k = k + 1
             x2z_send(1,i) = wk_sendcnt(i)
             x2z_send(2,i) = k
          endif
       enddo
       x2z_srank = k
       x2z_rmax = maxval(wk_recvcnt)
       x2z_smax = maxval(wk_sendcnt)
       x2z_srmax = max(x2z_rmax,x2z_smax)

       if(allocated(z2y_recv)) deallocate(z2y_recv)
       if(allocated(z2y_send)) deallocate(z2y_send)
       allocate(z2y_recv(2,0:nmrank-1))
       allocate(z2y_send(2,0:nmrank-1))
       z2y_send = 0
       z2y_recv = 0
       wk_recvcnt = 0
       wk_sendcnt = 0

       nx = xyz_fft_y(2,1)-xyz_fft_y(1,1)+1
       ny = xyz_fft_y(2,2)-xyz_fft_y(1,2)+1
       nz = xyz_fft_y(2,3)-xyz_fft_y(1,3)+1
       allocate(wk_mp_fft_y(nx*ny*nz) ,stat=ierr)
       do k = 0, nz-1
          do j = 0, ny-1
             do i = 0, nx-1
                wk_mp_fft_y(i+j*nx+k*nx*ny+1) = mp_fft_y(i+k*nx+j*nx*nz+1)
             enddo
          enddo
       enddo

       do i = 1, nel_fft_y(myrank)
          irank = map_fft_z(wk_mp_fft_y(i)) - 1
          wk_recvcnt(irank) = wk_recvcnt(irank) + 1
       enddo
       do i = 1, nel_fft_z(myrank)
          irank = map_fft_y(mp_fft_z(i)) - 1
          wk_sendcnt(irank) = wk_sendcnt(irank) + 1
       enddo
       k = 0
       do i = 0, nmrank - 1
          if(wk_recvcnt(i) /= 0) then
             k = k + 1
             z2y_recv(1,i) = wk_recvcnt(i)
             z2y_recv(2,i) = k
          endif
       enddo
       z2y_rrank = k
       k = 0
      do i = 0, nmrank - 1
          if(wk_sendcnt(i) /= 0) then
             k = k + 1
             z2y_send(1,i) = wk_sendcnt(i)
             z2y_send(2,i) = k
          endif
       enddo
       z2y_srank = k
       z2y_rmax = maxval(wk_recvcnt)
       z2y_smax = maxval(wk_sendcnt)
       z2y_srmax = max(z2y_rmax,z2y_smax)

       deallocate(wk_mp_fft_y)

#ifdef FFT_ALLTOALL
       max_send(1) = x2z_rmax
       max_send(2) = x2z_smax
       call mpi_allreduce(max_send,max_recv,2,mpi_integer,mpi_max,mpi_fft_xz_world,ierr)
       x2z_rmax = max_recv(1)
       x2z_smax = max_recv(2)
       max_send(1) = z2y_rmax
       max_send(2) = z2y_smax
       call mpi_allreduce(max_send,max_recv,2,mpi_integer,mpi_max,mpi_fft_zy_world,ierr)
       z2y_rmax = max_recv(1)
       z2y_smax = max_recv(2)
#endif

       if (ipri > 1) then
          write(nfout,'("m_FFT_Inverse_3D   --   myrank_g=",i4)') myrank_g
          write(nfout,'("x2z_send")')
          write(nfout,'(10(i8,", "))') (x2z_send(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (x2z_send(2,i),i=0,nmrank-1)
          write(nfout,'("x2z_recv")')
          write(nfout,'(10(i8,", "))') (x2z_recv(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (x2z_recv(2,i),i=0,nmrank-1)
          write(nfout,'("z2y_send")')
          write(nfout,'(10(i8,", "))') (z2y_send(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (z2y_send(2,i),i=0,nmrank-1)
          write(nfout,'("z2y_recv")')
          write(nfout,'(10(i8,", "))') (z2y_recv(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (z2y_recv(2,i),i=0,nmrank-1)
          call flush(nfout)
       endif

#ifdef __TIMER_DO__
  call timer_sta(135)
#endif
       nx = xyz_fft_x(2,1)-xyz_fft_x(1,1)+1
       ny = xyz_fft_x(2,2)-xyz_fft_x(1,2)+1
       nz = xyz_fft_x(2,3)-xyz_fft_x(1,3)+1
       if(kimg==1) then
          call dfftw_plan_many_dft_r2c(planx1,           1,   nx-2,  ny*nz, &
         &                             wk_afft_l,     nx-2,      1,     nx, &
         &                             wk_afft_l, (nx-2)/2,      1,   nx/2, &
         &                             FFTW_MEASURE )
       else
          call dfftw_plan_many_dft    (planx2,           1,     nx,  ny*nz, &
         &                             wk_afft_l,       nx,      1,     nx, &
         &                             wk_afft_l,       nx,      1,     nx, &
         &                             FFTW_FLAG, FFTW_MEASURE )
       endif

       nx = xyz_fft_z(2,1)-xyz_fft_z(1,1)+1
       ny = xyz_fft_z(2,2)-xyz_fft_z(1,2)+1
       nz = xyz_fft_z(2,3)-xyz_fft_z(1,3)+1
       if(kimg==1) then
          call dfftw_plan_many_dft    (planz1,          1,      nz, nx*ny/2, &
         &                             wk_afft_l,      nz, nx*ny/2,       1, &
         &                             wk_afft_l,      nz, nx*ny/2,       1, &
         &                             FFTW_FLAG, FFTW_MEASURE )
       else
          call dfftw_plan_many_dft    (planz2,          1,      nz,   nx*ny, &
         &                             wk_afft_l,      nz,   nx*ny,       1, &
         &                             wk_afft_l,      nz,   nx*ny,       1, &
         &                             FFTW_FLAG, FFTW_MEASURE )
       endif

       nx = xyz_fft_y(2,1)-xyz_fft_y(1,1)+1
       ny = xyz_fft_y(2,2)-xyz_fft_y(1,2)+1
       nz = xyz_fft_y(2,3)-xyz_fft_y(1,3)+1
       if(kimg==1) then
          call dfftw_plan_many_dft    (plany1,          1,      ny, nx*nz/2, &
         &                             wk_afft_l,      ny,       1,      ny, &
         &                             wk_afft_l,      ny,       1,      ny, &
         &                             FFTW_FLAG, FFTW_MEASURE )
       else
          call dfftw_plan_many_dft    (plany2,          1,      ny,   nz*nx, &
         &                             wk_afft_l,      ny,       1,      ny, &
         &                             wk_afft_l,      ny,       1,      ny, &
         &                             FFTW_FLAG, FFTW_MEASURE )
       end if
#ifdef __TIMER_DO__
  call timer_end(135)
#endif

       firstcall_inverse_3d = .false.
#ifdef __FAPP__
    call fapp_stop('fftwf_inverse_firstcall',1,1)
#endif
       go to 1000
#ifdef __TIMER_DO__
  call timer_end(134)
#endif
    endif

#ifdef FFT_ALLTOALL

    if (iesize /= savesize) then
#ifdef __TIMER_ETC__
  call timer_sta(264)
#endif
       if (allocated(wk_recv1)) deallocate(wk_recv1)
       if (allocated(wk_send1)) deallocate(wk_send1)
       if (allocated(wk_recv2)) deallocate(wk_recv2)
       if (allocated(wk_send2)) deallocate(wk_send2)
       if (allocated(wk_afft_y)) deallocate(wk_afft_y)
       allocate(wk_recv1(x2z_rmax*kimg*iesize,x2z_rrank), stat=ierr)
       allocate(wk_send1(x2z_smax*kimg*iesize,x2z_srank), stat=ierr)
       allocate(wk_recv2(z2y_rmax*kimg*iesize,z2y_rrank), stat=ierr)
       allocate(wk_send2(z2y_smax*kimg*iesize,z2y_srank), stat=ierr)
       allocate(wk_afft_y(max_elm*kimg,iesize) ,stat=ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_FFT_Inverse_3D :  Not allocate '
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 69, ierr)
        endif
       savesize = iesize
#ifdef __TIMER_ETC__
  call timer_end(264)
#endif
    end if

#else
!else ifdef FFT_ALLTOALL

    if (iesize > savesize) then
#ifdef __TIMER_ETC__
  call timer_sta(264)
#endif
       if (allocated(wk_recv1)) deallocate(wk_recv1)
       if (allocated(wk_send1)) deallocate(wk_send1)
       if (allocated(wk_recv2)) deallocate(wk_recv2)
       if (allocated(wk_send2)) deallocate(wk_send2)
       if (allocated(wk_afft_y)) deallocate(wk_afft_y)
       allocate(wk_recv1(x2z_rmax*kimg*iesize,x2z_rrank), stat=ierr)
       allocate(wk_send1(x2z_smax*kimg*iesize,x2z_srank), stat=ierr)
       allocate(wk_recv2(z2y_rmax*kimg*iesize,z2y_rrank), stat=ierr)
       allocate(wk_send2(z2y_smax*kimg*iesize,z2y_srank), stat=ierr)
       allocate(wk_afft_y(max_elm*kimg,iesize) ,stat=ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_FFT_Inverse_3D :  Not allocate '
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 69, ierr)
        endif
       savesize = iesize
#ifdef __TIMER_ETC__
  call timer_end(264)
#endif
    end if

#endif
!endif ifdef FFT_ALLTOALL

 2000 continue
#ifdef __TIMER_SUB__
    call mpi_barrier(mpi_ke_world, ierr)
    call timer_sta(1498)
#endif
    if (np_fft_x==0)then
       nx = 0
       nx = 0
       nz = 0
       go to 1000
    endif

! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
call stop_timer('FFT-I reg0')
#endif
! === TIMERTIMERTIMER ==========================================================
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP          SHARED(     nfout,wk_afft_l,iesize,req_r,req_s,sta_r,sta_s,    &
!$OMP                 wk_recv1,wk_send1,wk_recv2,wk_send2,wk_afft_y,          &
!$OMP                 planx1,planx2,planz1,planz2,plany1,plany2,              &
!$OMP                 x2z_recv, x2z_send, z2y_recv, z2y_send,itag,            &
!$OMP                 xyz_fft_y,xyz_fft_z,xyz_fft_x,fft_Z_x_dim,nis_fft_Z_x,  &
!$OMP                 nie_fft_Z_x,fft_X_z_dim,nis_fft_X_z,nie_fft_X_z,        &
!$OMP                 fft_Y_z_dim,nis_fft_Y_z,nie_fft_Y_z,fft_Z_y_dim,        &
!$OMP                 nie_fft_Z_y,nis_fft_Z_y,                                &
#ifdef FFT_ALLTOALL
!$OMP                 mpi_fft_xz_world, x2z_rmax, x2z_smax,                   &
!$OMP                 mpi_fft_zy_world, z2y_rmax, z2y_smax,                   &
#endif
!$OMP                 mpi_comm,myrank,nmrank,kimg,ierr                    )   &
!$OMP          PRIVATE(nx,ny,nz,i,j,k,l,ri,ib,iadd,lrank,icnt_send,icnt_recv, &
!$OMP                  irank,nnx,nny,nnz,iadd0,iadd1,nis,nie,jadd,kadd,plan   )

#ifdef __TIMER_SUB__
    call timer_sta(1499)
#endif
!
! X-axis (y-z div)
!
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call start_timer('FFT-I xfft')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================
    nx = xyz_fft_x(2,1)-xyz_fft_x(1,1)+1
    ny = xyz_fft_x(2,2)-xyz_fft_x(1,2)+1
    nz = xyz_fft_x(2,3)-xyz_fft_x(1,3)+1

    if (kimg == 1) then
       plan = planx1
    else
       plan = planx2
    endif

#ifdef __TIMER_DO__
  call timer_sta(254)
#endif
#ifdef __FAPP__
    call fapp_start('fftwf_inverse_fftw1',1,1)
#endif
!$OMP DO
    do ib = 1, iesize
       if (kimg==1) then
          call dfftw_execute_dft_r2c(plan,wk_afft_l(1,ib),wk_afft_l(1,ib))
       else
          call dfftw_execute_dft    (plan,wk_afft_l(1,ib),wk_afft_l(1,ib))
       endif
    enddo
#ifdef __FAPP__
    call fapp_stop('fftwf_inverse_fftw1',1,1)
#endif

#ifdef __TIMER_DO__
  call timer_end(254)
#endif
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call stop_timer('FFT-I xfft')
call start_timer('FFT-I x2z')
call start_timer('*FFT-I x2z sbuf')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================

#ifdef __TIMER_DO__
  call timer_sta(136)
#endif
#ifdef FFT_INTERCHANGE
    if (kimg == 1) then
!$OMP DO
       do ib = 1, iesize
          do k = 1, nz
             do j = 1, ny
                do l = 1, fft_Z_x_dim
                   nis = nis_fft_Z_x(l)
                   nie = nie_fft_Z_x(l)
                   nnx = nie - nis + 1
                   iadd0 = nnx*ny*nz*(ib-1)
                   do i = nis, nie
                      wk_send1(iadd0+(i-nis+1+(j-1)*nnx+(k-1)*nnx*ny),l) = wk_afft_l(i+(j-1)*nx+(k-1)*nx*ny,ib)
                   end do
                end do
             end do
          end do
       end do
    else
!$OMP DO
       do ib = 1, iesize
          do k = 1, nz
             do j = 1, ny
                do l = 1, fft_Z_x_dim
                   nis = nis_fft_Z_x(l)
                   nie = nie_fft_Z_x(l)
                   nnx = nie - nis + 1
                   do i = nis, nie
                      iadd1 = nnx*ny*nz*(ib-1)*2+(i-nis+1+(j-1)*nnx+(k-1)*nnx*ny)*2
                      iadd  = (i+(j-1)*nx+(k-1)*nx*ny)*2
                      wk_send1(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send1(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
    end if
#else
!$OMP DO
    do l = 1, fft_Z_x_dim
       nis = nis_fft_Z_x(l)
       nie = nie_fft_Z_x(l)
       nnx = nie - nis + 1
       if (kimg == 1) then
          do ib = 1, iesize
             iadd0 = nnx*ny*nz*(ib-1)
             do k = 1, nz
                do j = 1, ny
                   do i = nis, nie
                      wk_send1(iadd0+(i-nis+1+(j-1)*nnx+(k-1)*nnx*ny),l) = wk_afft_l(i+(j-1)*nx+(k-1)*nx*ny,ib)
                   end do
                end do
             end do
          end do
       else
          do ib = 1, iesize
             iadd0 = nnx*ny*nz*(ib-1)*2
             do k = 1, nz
                jadd = (k-1)*nnx*ny
                do j = 1, ny
                   kadd = (j-1)*nnx+jadd
                   do i = nis, nie
                      iadd1 = iadd0+(i-nis+1+kadd)*2
                      iadd  = (i+(j-1)*nx+(k-1)*nx*ny)*2
!fj --------------------
!                     do ri = 1, kimg
!                       wk_send1(iadd1+(ri-kimg),l) = wk_afft_l(iadd+(ri-kimg),ib)
!                     end do
                      wk_send1(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send1(iadd1  ,l) = wk_afft_l(iadd  ,ib)
!fj --------------------
                   end do
                end do
             end do
          end do
       end if
    end do
#endif
#ifdef __TIMER_DO__
  call timer_end(136)
#endif
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call stop_timer('*FFT-I x2z sbuf')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================

#ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_fft_xz_world)
  call timer_sta(137)
#endif

! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call start_timer('*FFT-I x2z a2a')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================
!$OMP MASTER

    call MPI_ALLTOALL(wk_send1, x2z_smax*kimg*iesize, mpi_double_precision,   &
   &                  wk_recv1, x2z_rmax*kimg*iesize, mpi_double_precision,   &
   &                                                  mpi_fft_xz_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Inverse_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 70, ierr)
     endif

!$OMP END MASTER
!$OMP BARRIER
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call stop_timer('*FFT-I x2z a2a')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================

#ifdef __TIMER_DO__
  call timer_sta(138)
#endif
#ifdef __TIMER_DO__
  call timer_end(138)
#endif

#else
!else ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_comm)
  call timer_sta(137)
#endif

! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call start_timer('*FFT-I x2z sr')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================
!$OMP MASTER
    icnt_recv = 0
    lrank = myrank + 1
    if (lrank > (nmrank -1)) lrank = 0
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank -1)) lrank = 0
       if ((lrank /= myrank) .and. (x2z_recv(1,lrank) /= 0)) then
          call mpi_irecv(wk_recv1(1,x2z_recv(2,lrank)), x2z_recv(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_r(icnt_recv), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Inverse_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 70, ierr)
           endif
          icnt_recv = icnt_recv + 1
       endif
    enddo

    icnt_send = 0
    lrank = myrank
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank -1)) lrank = 0
       if ((lrank /= myrank) .and. (x2z_send(1,lrank) /= 0)) then
          call mpi_isend(wk_send1(1,x2z_send(2,lrank)), x2z_send(1,lrank)*kimg*iesize, &
                         mpi_double_precision, lrank, itag, mpi_comm, req_s(icnt_send), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Inverse_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 71, ierr)
           endif
          icnt_send = icnt_send + 1
       endif
    enddo
!$OMP END MASTER
!$OMP BARRIER

#ifdef __TIMER_DO__
  call timer_sta(138)
#endif
!$OMP DO
    do i = 1, x2z_recv(1,myrank)*kimg*iesize
       wk_recv1(i,x2z_recv(2,myrank)) = wk_send1(i,x2z_send(2,myrank))
    enddo
#ifdef __TIMER_DO__
  call timer_end(138)
#endif

!$OMP MASTER
    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Inverse_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 72, ierr)
     endif

    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Inverse_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 73, ierr)
     endif
!$OMP END MASTER
!$OMP BARRIER
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call stop_timer('*FFT-I x2z sr')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================

#endif
!endif ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_end(137)
#endif

#ifdef __TIMER_DO__
  call timer_sta(139)
#endif
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call start_timer('*FFT-I x2z rbuf')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================
    nx = xyz_fft_z(2,1)-xyz_fft_z(1,1)+1
    ny = xyz_fft_z(2,2)-xyz_fft_z(1,2)+1
    nz = xyz_fft_z(2,3)-xyz_fft_z(1,3)+1
!$OMP DO
    do l = 1, fft_X_z_dim
       nis = nis_fft_X_z(l)
       nie = nie_fft_X_z(l)
       nnz = nie - nis + 1
       if (kimg == 1) then
          do ib = 1, iesize
             iadd0 = nx*ny*nnz*(ib-1)
             do k = nis, nie
                do j = 1, ny
                   do i = 1, nx
                      wk_afft_l(i+(j-1)*nx+(k-1)*nx*ny,ib) = wk_recv1(iadd0+i+(j-1)*nx+(k-nis)*nx*ny,l)
                   end do
                end do
             end do
          end do
       else
          do ib = 1, iesize
             iadd0 = nx*ny*nnz*(ib-1)*2
             do k = nis, nie
                do j = 1, ny
                   do i = 1, nx
                      iadd1 = iadd0+(i+(j-1)*nx+(k-nis)*nx*ny)*2
                      iadd  = (i+(j-1)*nx+(k-1)*nx*ny)*2
                      wk_afft_l(iadd-1,ib) = wk_recv1(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv1(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end if
    end do
#ifdef __TIMER_DO__
  call timer_end(139)
#endif
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call stop_timer('*FFT-I x2z rbuf')
call stop_timer('FFT-I x2z')
call start_timer('FFT-I zfft')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================

!
! Z-axis (x-y div)
!
    nx = xyz_fft_z(2,1)-xyz_fft_z(1,1)+1
    ny = xyz_fft_z(2,2)-xyz_fft_z(1,2)+1
    nz = xyz_fft_z(2,3)-xyz_fft_z(1,3)+1

    if (kimg == 1) then
       plan = planz1
    else
       plan = planz2
    endif
#ifdef __TIMER_DO__
  call timer_sta(255)
#endif

#ifdef __FAPP__
    call fapp_start('fftwf_inverse_fftw2',1,1)
#endif
!$OMP DO
    do ib = 1, iesize
       call dfftw_execute_dft(plan,wk_afft_l(1,ib),wk_afft_l(1,ib))
    end do
#ifdef __FAPP__
    call fapp_stop('fftwf_inverse_fftw2',1,1)
#endif
#ifdef __TIMER_DO__
  call timer_end(255)
#endif

! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call stop_timer('FFT-I zfft')
call start_timer('FFT-I z2y')
call start_timer('*FFT-I z2y sbuf')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_DO__
  call timer_sta(140)
#endif
#ifdef __TIMER_DO__
  call timer_end(140)
#endif

#ifdef __TIMER_DO__
  call timer_sta(141)
#endif
    nx = xyz_fft_z(2,1)-xyz_fft_z(1,1)+1
    ny = xyz_fft_z(2,2)-xyz_fft_z(1,2)+1
    nz = xyz_fft_z(2,3)-xyz_fft_z(1,3)+1
!$OMP DO
    do l = 1, fft_Y_z_dim
       nis = nis_fft_Y_z(l)
       nie = nie_fft_Y_z(l)
       nnz = nie - nis + 1
       if (kimg == 1) then
          do ib = 1, iesize
             iadd0 = nx*ny*nnz*(ib-1)
             do k = nis, nie
                do j = 1, ny
                   do i = 1, nx
                      wk_send2(iadd0+(i+(j-1)*nx+(k-nis)*nx*ny),l) = wk_afft_l(i+(j-1)*nx+(k-1)*nx*ny,ib)
                   end do
                end do
             end do
          end do
       else
          do ib = 1, iesize
             iadd0 = nx*ny*nnz*(ib-1)*2
             do k = nis, nie
                jadd = (k-nis)*nx*ny
                do j = 1, ny
                   kadd = (j-1)*nx+jadd
                   do i = 1, nx
                      iadd1 = iadd0+(i+kadd)*2
                      iadd  = (i+(j-1)*nx+(k-1)*nx*ny)*2
!fj --------------------
!                     do ri = 1, kimg
!                        wk_send2(iadd1+(ri-kimg),l) = wk_afft_l(iadd+(ri-kimg),ib)
!                     end do
                      wk_send2(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send2(iadd1  ,l) = wk_afft_l(iadd  ,ib)
!fj --------------------
                   end do
                end do
             end do
          end do
       end if
    end do
#ifdef __TIMER_DO__
  call timer_end(141)
#endif

! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call stop_timer('*FFT-I z2y sbuf')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================
#ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_fft_zy_world)
  call timer_sta(142)
#endif

! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call start_timer('*FFT-I z2y a2a')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================
!$OMP MASTER

    call MPI_ALLTOALL(wk_send2, z2y_smax*kimg*iesize, mpi_double_precision,   &
   &                  wk_recv2, z2y_rmax*kimg*iesize, mpi_double_precision,   &
   &                                                  mpi_fft_zy_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Inverse_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 76, ierr)
     endif

!$OMP END MASTER
!$OMP BARRIER
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call stop_timer('*FFT-I z2y a2a')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================

#ifdef __TIMER_DO__
  call timer_sta(143)
#endif
#ifdef __TIMER_DO__
  call timer_end(143)
#endif

#else
!else ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_comm)
  call timer_sta(142)
#endif

! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call start_timer('*FFT-I z2y sr')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================
!$OMP MASTER
    icnt_recv = 0
    lrank = myrank + 1
    if (lrank > (nmrank-1)) lrank = 0
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank-1)) lrank = 0
       if ((lrank /= myrank) .and. (z2y_recv(1,lrank) /= 0)) then
          call mpi_irecv(wk_recv2(1,z2y_recv(2,lrank)), z2y_recv(1,lrank)*kimg*iesize, &
                         mpi_double_precision, lrank, itag, mpi_comm, req_r(icnt_recv), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Inverse_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 76, ierr)
           endif
          icnt_recv = icnt_recv + 1
       endif
    enddo

    icnt_send = 0
    lrank = myrank
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank -1)) lrank = 0
       if ((lrank /= myrank) .and. (z2y_send(1,lrank) /= 0)) then
          call mpi_isend(wk_send2(1,z2y_send(2,lrank)), z2y_send(1,lrank)*kimg*iesize, &
                         mpi_double_precision, lrank, itag, mpi_comm, req_s(icnt_send), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Inverse_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 77, ierr)
           endif
          icnt_send = icnt_send + 1
       endif
    enddo
!$OMP END MASTER
!$OMP BARRIER

#ifdef __TIMER_DO__
  call timer_sta(143)
#endif
!$OMP DO
    do i = 1, z2y_recv(1,myrank)*kimg*iesize
       wk_recv2(i,z2y_recv(2,myrank)) = wk_send2(i,z2y_send(2,myrank))
    enddo
#ifdef __TIMER_DO__
  call timer_end(143)
#endif

!$OMP MASTER
    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Inverse_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 78, ierr)
     endif

    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Inverse_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 79, ierr)
     endif
!$OMP END MASTER
!$OMP BARRIER
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call stop_timer('*FFT-I z2y sr')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================

#endif
!endif ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_end(142)
#endif

#ifdef __TIMER_DO__
  call timer_sta(144)
#endif
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call start_timer('*FFT-I z2y rbuf')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================
    nx = xyz_fft_y(2,1)-xyz_fft_y(1,1)+1
    ny = xyz_fft_y(2,2)-xyz_fft_y(1,2)+1
    nz = xyz_fft_y(2,3)-xyz_fft_y(1,3)+1
#ifdef FFT_INTERCHANGE
    if (kimg == 1) then
!$OMP DO
       do ib = 1, iesize
          do k = 1, nz
             do l = 1, fft_Z_y_dim
                nis = nis_fft_Z_y(l)
                nie = nie_fft_Z_y(l)
                nny = nie - nis + 1
                iadd0 = nx*nny*nz*(ib-1)
                do j = nis, nie
                   do i = 1, nx
                      wk_afft_y(i+(j-1)*nx+(k-1)*nx*ny,ib) = wk_recv2(iadd0+i+(j-nis)*nx+(k-1)*nx*nny,l)
                   end do
                end do
             end do
          end do
       end do
    else
!$OMP DO
       do ib = 1, iesize
          do k = 1, nz
             do l = 1, fft_Z_y_dim
                nis = nis_fft_Z_y(l)
                nie = nie_fft_Z_y(l)
                nny = nie - nis + 1
                do j = nis, nie
                   do i = 1, nx
                      iadd1 = nx*nny*nz*(ib-1)*2+(i+(j-nis)*nx+(k-1)*nx*nny)*2
                      iadd  = (i+(j-1)*nx+(k-1)*nx*ny)*2
                      wk_afft_y(iadd-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_y(iadd  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end do
    end if
#else
!$OMP DO
    do l = 1, fft_Z_y_dim
       nis = nis_fft_Z_y(l)
       nie = nie_fft_Z_y(l)
       nny = nie - nis + 1
       if (kimg == 1) then
          do ib = 1, iesize
             iadd0 = nx*nny*nz*(ib-1)
             do k = 1, nz
                do j = nis, nie
                   do i = 1, nx
                      wk_afft_y(i+(j-1)*nx+(k-1)*nx*ny,ib) = wk_recv2(iadd0+i+(j-nis)*nx+(k-1)*nx*nny,l)
                   end do
                end do
             end do
          end do
       else
          do ib = 1, iesize
             iadd0 = nx*nny*nz*(ib-1)*2
             do k = 1, nz
                do j = nis, nie
                   do i = 1, nx
                      iadd1 = iadd0+(i+(j-nis)*nx+(k-1)*nx*nny)*2
                      iadd  = (i+(j-1)*nx+(k-1)*nx*ny)*2
                      wk_afft_y(iadd-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_y(iadd  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end if
    end do
#endif
#ifdef __TIMER_DO__
  call timer_end(144)
#endif

#ifdef __TIMER_DO__
  call timer_sta(145)
#endif
!   do ib = 1, iesize
!      do k = 0, nz-1
!         do j = 0, ny-1
!            do i = 0, nx-1
!   do ri = 1, kimg
!               wk_afft_l((i+k*nx+j*nx*nz+1)*kimg+(ri-kimg),ib) = wk_afft_y((i+j*nx+k*nx*ny+1)*kimg+(ri-kimg),ib)
!   enddo
!            enddo
!         enddo
!      enddo
!   enddo
#ifdef FFT_INTERCHANGE
    if (kimg == 1) then
       do ib = 1, iesize
!$OMP DO
          do k = 1, nz
             do j = 1, ny
                do i = 1, nx/2
                   wk_afft_l((j+(i-1)*ny+(k-1)*ny*nx/2)*2-1,ib) = wk_afft_y((i+(j-1)*nx/2+(k-1)*nx/2*ny)*2-1,ib)
                   wk_afft_l((j+(i-1)*ny+(k-1)*ny*nx/2)*2  ,ib) = wk_afft_y((i+(j-1)*nx/2+(k-1)*nx/2*ny)*2  ,ib)
                enddo
             enddo
          enddo
       enddo
    else
       do ib = 1, iesize
!$OMP DO
          do k = 0, nz-1
             do j = 0, ny-1
                do i = 0, nx-1
                   wk_afft_l((j+i*ny+k*ny*nx+1)*2-1,ib) = wk_afft_y((i+j*nx+k*nx*ny+1)*2-1,ib)
                   wk_afft_l((j+i*ny+k*ny*nx+1)*2  ,ib) = wk_afft_y((i+j*nx+k*nx*ny+1)*2  ,ib)
                enddo
             enddo
          enddo
       enddo
    end if
#else
    if (kimg == 1) then
!$OMP DO
       do ib = 1, iesize
          do k = 1, nz
             do j = 1, ny
                do i = 1, nx/2
                   wk_afft_l((j+(i-1)*ny+(k-1)*ny*nx/2)*2-1,ib) = wk_afft_y((i+(j-1)*nx/2+(k-1)*nx/2*ny)*2-1,ib)
                   wk_afft_l((j+(i-1)*ny+(k-1)*ny*nx/2)*2  ,ib) = wk_afft_y((i+(j-1)*nx/2+(k-1)*nx/2*ny)*2  ,ib)
                enddo
             enddo
          enddo
       enddo
    else
!$OMP DO
       do ib = 1, iesize
          do k = 0, nz-1
             do j = 0, ny-1
                do i = 0, nx-1
                   wk_afft_l((j+i*ny+k*ny*nx+1)*2-1,ib) = wk_afft_y((i+j*nx+k*nx*ny+1)*2-1,ib)
                   wk_afft_l((j+i*ny+k*ny*nx+1)*2  ,ib) = wk_afft_y((i+j*nx+k*nx*ny+1)*2  ,ib)
                enddo
             enddo
          enddo
       enddo
    end if
#endif
!fj --------------------
#ifdef __TIMER_DO__
  call timer_end(145)
#endif
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call stop_timer('*FFT-I z2y rbuf')
call stop_timer('FFT-I z2y')
call start_timer('FFT-I yfft')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================

!
! Y-axis (z-x div)
!
    if (kimg == 1) then
       plan = plany1
    else
       plan = plany2
    endif
#ifdef __TIMER_DO__
  call timer_sta(256)
#endif
#ifdef __FAPP__
    call fapp_start('fftwf_inverse_fftw3',1,1)
#endif
!$OMP DO
    do ib = 1,iesize
       call dfftw_execute_dft(plan,wk_afft_l(1,ib),wk_afft_l(1,ib))
    end do
#ifdef __FAPP__
    call fapp_stop('fftwf_inverse_fftw3',1,1)
#endif
#ifdef __TIMER_DO__
  call timer_end(256)
#endif
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
!$OMP single
call stop_timer('FFT-I yfft')
!$OMP end single
#endif
! === TIMERTIMERTIMER ==========================================================

#ifdef __TIMER_SUB__
  call timer_end(1498)
#endif

!$OMP END PARALLEL

#ifdef __TIMER_SUB__
  call timer_end(1499)
#endif

1000 continue

!F  if (allocated(req_r))       deallocate(req_r)
!F  if (allocated(req_s))       deallocate(req_s)
!F  if (allocated(sta_r))       deallocate(sta_r)
!F  if (allocated(sta_s))       deallocate(sta_s)
!F  if (allocated(wk_recvcnt))  deallocate(wk_recvcnt)
!F  if (allocated(wk_sendcnt))  deallocate(wk_sendcnt)
    if (allocated(wk_mp_fft_y)) deallocate(wk_mp_fft_y)
!F  if (allocated(wk_afft_y))   deallocate(wk_afft_y)

#ifdef __TIMER_SUB__
    call timer_end(106)
#endif
  end subroutine m_FFT_Inverse_3D

#endif
!endif ifdef FFT_FFTW_OLD

!------------------------------------------------------------------------------
#else
!else FFT_USE_SSL2

!------------------------------------------------------------------------------
! For Fujitsu SSL2
  subroutine m_FFT_Inverse_3D(nfout,wk_afft_l, wk_size, iesize)
    integer, intent(in) :: nfout, wk_size, iesize
    real(kind=DP), dimension(wk_size*kimg,iesize), intent(inout) :: wk_afft_l

    real(kind=DP), allocatable, dimension(:) :: wk_allfft
    integer,      allocatable, dimension(:)   :: wk_recvdsp
    integer, parameter :: FFTW_MEASURE=0
    integer, parameter :: FFTW_ESTIMATE=64
    integer, parameter :: FFTW_RANK=1
    integer            :: FFTW_FLAG
    integer,dimension(1) :: NFFTW3, NEMBED, NEREAL
    integer :: nx, ny, nz, i, j, k, l, ri, ib, iadd, lrank
    integer :: irank, itag, icnt_send, icnt_recv
    integer,save, allocatable, dimension(:)   :: req_r, req_s
    integer,save, allocatable, dimension(:,:) :: sta_r, sta_s
    integer,save, allocatable, dimension(:)   :: wk_recvcnt, wk_sendcnt
!   real(kind=DP),allocatable, dimension(:,:) :: wk_recv, wk_send
    real(kind=DP),allocatable, dimension(:,:),save :: wk_recv1, wk_send1
    real(kind=DP),allocatable, dimension(:,:),save :: wk_recv2, wk_send2
    integer,      allocatable, dimension(:)   :: wk_mp_fft_y
    real(kind=DP),allocatable, dimension(:,:),save  :: wk_afft_y
    real(kind=DP),allocatable, dimension(:,:,:) :: wk_gather
    integer :: mpi_comm, myrank, nmrank
    integer(kind=8),save :: planx1 = 0, planx2 = 0
    integer(kind=8),save :: plany1 = 0, plany2 = 0
    integer(kind=8),save :: planz1 = 0, planz2 = 0
    integer(kind=8) :: plan
    integer :: max_x, max_y, max_z
    integer,save :: max_elm = 0, savesize = 0

    integer,save, allocatable, dimension(:,:) :: x2z_recv, x2z_send, z2y_recv, z2y_send
    integer,save :: x2z_rrank, x2z_srank, z2y_rrank, z2y_srank
    integer,save :: x2z_rmax, x2z_smax, x2z_srmax, z2y_rmax, z2y_smax, z2y_srmax

    integer ::  nnx, nny, nnz, iadd0, iadd1, nis, nie, jadd, kadd, nxp, nyp, nzp
    integer                  , dimension(3) :: isin,nsize
#ifdef FFT_ALLTOALL
    integer,dimension(2) :: max_send,max_recv
#endif

#ifdef __TIMER_SUB__
    call timer_barrier(mpi_ke_world)
    call timer_sta(106)
#endif

    if (kimg == 1) then
       FFTW_FLAG = -1      ! FFTW_FORWARD
    else
       FFTW_FLAG = +1      ! FFTW_BACKWARD
    endif

!!  mpi_comm = mpi_kg_world
!!  myrank = myrank_e
!!  nmrank = nrank_e
    mpi_comm = mpi_ke_world
    myrank = myrank_g
    nmrank = nrank_g
    itag = 10

    if (np_fft_x==0)then
       nx = 0
       nx = 0
       nz = 0
!      go to 1000
       go to 2000
    endif

    if (firstcall_inverse_3d) then
       savesize = 0
#ifdef __TIMER_DO__
  call timer_sta(134)
#endif
       max_x = maxval(nel_fft_x(:))
       max_y = maxval(nel_fft_y(:))
       max_z = maxval(nel_fft_z(:))
       max_elm = max(max_x,max_y,max_z)
       if(allocated(req_r)) deallocate(req_r)
       if(allocated(req_s)) deallocate(req_s)
       if(allocated(sta_r)) deallocate(sta_r)
       if(allocated(sta_s)) deallocate(sta_s)
       if(allocated(wk_recvcnt)) deallocate(wk_recvcnt)
       if(allocated(wk_sendcnt)) deallocate(wk_sendcnt)
       allocate(req_r(0:nmrank-1), stat=ierr)
       allocate(req_s(0:nmrank-1), stat=ierr)
       allocate(sta_r(MPI_STATUS_SIZE,0:nmrank-1), stat=ierr)
       allocate(sta_s(MPI_STATUS_SIZE,0:nmrank-1), stat=ierr)
       allocate(wk_recvcnt(0:nmrank-1), stat=ierr)
       allocate(wk_sendcnt(0:nmrank-1), stat=ierr)
       if (ierr /= 0) then
          write(nfout,*)' m_FFT_Inverse_3D :  Not allocate '
          call flush(nfout)
          call mpi_abort(mpi_comm_world, 68, ierr)
       endif

       if(allocated(x2z_recv)) deallocate(x2z_recv)
       if(allocated(x2z_send)) deallocate(x2z_send)
       allocate(x2z_recv(2,0:nmrank-1))
       allocate(x2z_send(2,0:nmrank-1))
       x2z_send = 0
       x2z_recv = 0
       wk_recvcnt = 0
       wk_sendcnt = 0
       do i = 1, nel_fft_z(myrank)
          irank = map_fft_x(mp_fft_z(i)) - 1
          wk_recvcnt(irank) = wk_recvcnt(irank) + 1
       enddo
       do i = 1, nel_fft_x(myrank)
          irank = map_fft_z(mp_fft_x(i)) - 1
          wk_sendcnt(irank) = wk_sendcnt(irank) + 1
       enddo
       k = 0
       do i = 0, nmrank - 1
          if(wk_recvcnt(i) /= 0) then
             k = k + 1
             x2z_recv(1,i) = wk_recvcnt(i)
             x2z_recv(2,i) = k
          endif
       enddo
       x2z_rrank = k
       k = 0
       do i = 0, nmrank - 1
          if(wk_sendcnt(i) /= 0) then
             k = k + 1
             x2z_send(1,i) = wk_sendcnt(i)
             x2z_send(2,i) = k
          endif
       enddo
       x2z_srank = k
       x2z_rmax = maxval(wk_recvcnt)
       x2z_smax = maxval(wk_sendcnt)
       x2z_srmax = max(x2z_rmax,x2z_smax)

       if(allocated(z2y_recv)) deallocate(z2y_recv)
       if(allocated(z2y_send)) deallocate(z2y_send)
       allocate(z2y_recv(2,0:nmrank-1))
       allocate(z2y_send(2,0:nmrank-1))
       z2y_send = 0
       z2y_recv = 0
       wk_recvcnt = 0
       wk_sendcnt = 0

       nx = xyz_fft_y(2,1)-xyz_fft_y(1,1)+1
       ny = xyz_fft_y(2,2)-xyz_fft_y(1,2)+1
       nz = xyz_fft_y(2,3)-xyz_fft_y(1,3)+1
       allocate(wk_mp_fft_y(nx*ny*nz) ,stat=ierr)
       do k = 0, nz-1
          do j = 0, ny-1
             do i = 0, nx-1
                wk_mp_fft_y(i+j*nx+k*nx*ny+1) = mp_fft_y(i+k*nx+j*nx*nz+1)
             enddo
          enddo
       enddo

       do i = 1, nel_fft_y(myrank)
          irank = map_fft_z(wk_mp_fft_y(i)) - 1
          wk_recvcnt(irank) = wk_recvcnt(irank) + 1
       enddo
       do i = 1, nel_fft_z(myrank)
          irank = map_fft_y(mp_fft_z(i)) - 1
          wk_sendcnt(irank) = wk_sendcnt(irank) + 1
       enddo
       k = 0
       do i = 0, nmrank - 1
          if(wk_recvcnt(i) /= 0) then
             k = k + 1
             z2y_recv(1,i) = wk_recvcnt(i)
             z2y_recv(2,i) = k
          endif
       enddo
       z2y_rrank = k
       k = 0
      do i = 0, nmrank - 1
          if(wk_sendcnt(i) /= 0) then
             k = k + 1
             z2y_send(1,i) = wk_sendcnt(i)
             z2y_send(2,i) = k
          endif
       enddo
       z2y_srank = k
       z2y_rmax = maxval(wk_recvcnt)
       z2y_smax = maxval(wk_sendcnt)
       z2y_srmax = max(z2y_rmax,z2y_smax)

       deallocate(wk_mp_fft_y)

#ifdef FFT_ALLTOALL
       max_send(1) = x2z_rmax
       max_send(2) = x2z_smax
       call mpi_allreduce(max_send,max_recv,2,mpi_integer,mpi_max,mpi_fft_xz_world,ierr)
       x2z_rmax = max_recv(1)
       x2z_smax = max_recv(2)
       max_send(1) = z2y_rmax
       max_send(2) = z2y_smax
       call mpi_allreduce(max_send,max_recv,2,mpi_integer,mpi_max,mpi_fft_zy_world,ierr)
       z2y_rmax = max_recv(1)
       z2y_smax = max_recv(2)
#endif

       if (ipri > 1) then
          write(nfout,'("m_FFT_Inverse_3D   --   myrank_g=",i4)') myrank_g
          write(nfout,'("x2z_send")')
          write(nfout,'(10(i8,", "))') (x2z_send(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (x2z_send(2,i),i=0,nmrank-1)
          write(nfout,'("x2z_recv")')
          write(nfout,'(10(i8,", "))') (x2z_recv(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (x2z_recv(2,i),i=0,nmrank-1)
          write(nfout,'("z2y_send")')
          write(nfout,'(10(i8,", "))') (z2y_send(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (z2y_send(2,i),i=0,nmrank-1)
          write(nfout,'("z2y_recv")')
          write(nfout,'(10(i8,", "))') (z2y_recv(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (z2y_recv(2,i),i=0,nmrank-1)
          call flush(nfout)
       endif

!X     nx = xyz_fft_x(2,1)-xyz_fft_x(1,1)+1
!X     ny = xyz_fft_x(2,2)-xyz_fft_x(1,2)+1
!X     nz = xyz_fft_x(2,3)-xyz_fft_x(1,3)+1
!X     if(kimg==1) then
!X        NFFTW3(1) = nx - 2
!X        NEREAL(1) = nx - 2
!X        NEMBED(1) = (nx - 2)/2
!X        call dfftw_plan_many_dft_r2c(planx1, FFTW_RANK, NFFTW3, ny*nz,   &
!X       &                             wk_afft_l, NEREAL, 1, nx,          &
!X       &                             wk_afft_l, NEMBED, 1, nx/2,        &
!X       &                             FFTW_ESTIMATE)
!X     else
!X        NFFTW3(1) = nx
!X        NEMBED(1) = nx
!X        call dfftw_plan_many_dft    (planx2, FFTW_RANK, NFFTW3, ny*nz,   &
!X       &                             wk_afft_l, NEMBED, 1, nx,          &
!X       &                             wk_afft_l, NEMBED, 1, nx,          &
!X       &                             FFTW_FLAG, FFTW_ESTIMATE)
!X     endif

!X     nx = xyz_fft_z(2,1)-xyz_fft_z(1,1)+1
!X     ny = xyz_fft_z(2,2)-xyz_fft_z(1,2)+1
!X     nz = xyz_fft_z(2,3)-xyz_fft_z(1,3)+1
!X     if(kimg==1) then
!X        NFFTW3(1) = nz
!X        NEMBED(1) = nz
!X        call dfftw_plan_many_dft    (planz1, FFTW_RANK, NFFTW3, nx*ny/2, &
!X       &                             wk_afft_l, NEMBED, nx*ny/2, 1,     &
!X       &                             wk_afft_l, NEMBED, nx*ny/2, 1,     &
!X       &                             FFTW_FLAG, FFTW_ESTIMATE)
!X     else
!X        NFFTW3(1) = nz
!X        NEMBED(1) = nz
!X        call dfftw_plan_many_dft    (planz2, FFTW_RANK, NFFTW3, nx*ny,   &
!X       &                             wk_afft_l, NEMBED, nx*ny, 1,       &
!X       &                             wk_afft_l, NEMBED, nx*ny, 1,       &
!X       &                             FFTW_FLAG, FFTW_ESTIMATE)
!X     endif
!
!X     nx = xyz_fft_y(2,1)-xyz_fft_y(1,1)+1
!X     ny = xyz_fft_y(2,2)-xyz_fft_y(1,2)+1
!X     nz = xyz_fft_y(2,3)-xyz_fft_y(1,3)+1
!X     if(kimg==1) then
!X        NFFTW3(1) = ny
!X        NEMBED(1) = ny
!X        call dfftw_plan_many_dft    (plany1, FFTW_RANK, NFFTW3, nz*nx, &
!X       &                             wk_afft_l, ny,     1,      ny/2 , &
!X       &                             wk_afft_l, ny,     1,      ny/2 , &
!X       &                             FFTW_FLAG, FFTW_ESTIMATE)
!X     else
!X        NFFTW3(1) = ny
!X        NEMBED(1) = ny
!X        call dfftw_plan_many_dft    (plany2, FFTW_RANK, NFFTW3, nz*nx,   &
!X       &                             wk_afft_l, NEMBED, nz*nx, 1,      &
!X       &                             wk_afft_l, NEMBED, nz*nx, 1,      &
!X       &                             FFTW_FLAG, FFTW_ESTIMATE)
!X     end if

       firstcall_inverse_3d = .false.

#ifdef __TIMER_DO__
  call timer_end(134)
#endif
    endif

#ifdef FFT_ALLTOALL

    if (iesize /= savesize) then
#ifdef __TIMER_ETC__
  call timer_sta(264)
#endif
       if (allocated(wk_recv1)) deallocate(wk_recv1)
       if (allocated(wk_send1)) deallocate(wk_send1)
       if (allocated(wk_recv2)) deallocate(wk_recv2)
       if (allocated(wk_send2)) deallocate(wk_send2)
       if (allocated(wk_afft_y)) deallocate(wk_afft_y)
       allocate(wk_recv1(x2z_rmax*kimg*iesize,x2z_rrank), stat=ierr)
       allocate(wk_send1(x2z_smax*kimg*iesize,x2z_srank), stat=ierr)
       allocate(wk_recv2(z2y_rmax*kimg*iesize,z2y_rrank), stat=ierr)
       allocate(wk_send2(z2y_smax*kimg*iesize,z2y_srank), stat=ierr)
       allocate(wk_afft_y(wk_size*kimg,iesize) ,stat=ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_FFT_Inverse_3D :  Not allocate '
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 69, ierr)
        endif
       savesize = iesize
#ifdef __TIMER_ETC__
  call timer_end(264)
#endif
    end if

#else
!else ifdef FFT_ALLTOALL

    if (iesize > savesize) then
#ifdef __TIMER_ETC__
  call timer_sta(264)
#endif
       if (allocated(wk_recv1)) deallocate(wk_recv1)
       if (allocated(wk_send1)) deallocate(wk_send1)
       if (allocated(wk_recv2)) deallocate(wk_recv2)
       if (allocated(wk_send2)) deallocate(wk_send2)
       if (allocated(wk_afft_y)) deallocate(wk_afft_y)
       allocate(wk_recv1(x2z_rmax*kimg*iesize,x2z_rrank), stat=ierr)
       allocate(wk_send1(x2z_smax*kimg*iesize,x2z_srank), stat=ierr)
       allocate(wk_recv2(z2y_rmax*kimg*iesize,z2y_rrank), stat=ierr)
       allocate(wk_send2(z2y_smax*kimg*iesize,z2y_srank), stat=ierr)
!      allocate(wk_afft_y(max_elm*kimg,iesize) ,stat=ierr)
       allocate(wk_afft_y(wk_size*kimg,iesize) ,stat=ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_FFT_Inverse_3D :  Not allocate '
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 69, ierr)
        endif
       savesize = iesize
#ifdef __TIMER_ETC__
  call timer_end(264)
#endif
    end if

#endif
!endif ifdef FFT_ALLTOALL

 2000 continue
#ifdef __TIMER_SUB__
    call mpi_barrier(mpi_ke_world, ierr)
    call timer_sta(1498)
#endif
    if (np_fft_x==0)then
       nx = 0
       nx = 0
       nz = 0
       go to 1000
    endif

#ifdef __TIMER_SUB__
    call timer_sta(1499)
#endif
!
! X-axis (y-z div)
!
    nx = xyz_fft_x(2,1)-xyz_fft_x(1,1)+1
    ny = xyz_fft_x(2,2)-xyz_fft_x(1,2)+1
    nz = xyz_fft_x(2,3)-xyz_fft_x(1,3)+1
    nxp = nx

#ifdef __TIMER_DO__
  call timer_sta(254)
#endif
    if (kimg==1) then
       nsize(1:3) = (/nx-2,ny,nz/)
       isin(1:3)  = (/1,0,0/)
!OCL INDEPENDENT[ DVMRF2 ]
!OCL PARALLEL_STRONG
       do ib = 1, iesize
          call DVMRF2(wk_afft_l(1,ib),nsize,3,isin,1,ierr)
       end do
    else
       nsize(1:3) = (/nx,ny,nz/)
       isin(1:3)  = (/-1,0,0/)
!OCL INDEPENDENT[ DVMCF2 ]
!OCL PARALLEL_STRONG
       do ib = 1, iesize
          call DVMCF2(wk_afft_l(1,ib),nsize,3,isin,ierr)
       end do
    end if
#ifdef __TIMER_DO__
  call timer_end(254)
#endif

#ifdef __TIMER_DO__
  call timer_sta(136)
#endif

    do l = 1, fft_Z_x_dim
       nis = nis_fft_Z_x(l)
       nie = nie_fft_Z_x(l)
       nnx = nie - nis + 1
       if (kimg == 1) then
          do ib = 1, iesize
             iadd0 = nnx*ny*nz*(ib-1)
             do k = 1, nz
                do j = 1, ny
                   do i = nis, nie
                      wk_send1(iadd0+(i-nis+1+(j-1)*nnx+(k-1)*nnx*ny),l) = wk_afft_l(i+(j-1)*nx+(k-1)*nx*ny,ib)
                   end do
                end do
             end do
          end do
       else
          do ib = 1, iesize
             iadd0 = nnx*ny*nz*(ib-1)*2
             do k = 1, nz
                do j = 1, ny
                   do i = nis, nie
                      iadd  = (i+(j-1)*nxp+(k-1)*nxp*ny)*2
                      iadd1 = iadd0+(i-nis+1+(j-1)*nnx+(k-1)*nnx*ny)*2
!fj --------------------
!                     do ri = 1, kimg
!                       wk_send1(iadd1+(ri-kimg),l) = wk_afft_l(iadd+(ri-kimg),ib)
!                     end do
                      wk_send1(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send1(iadd1  ,l) = wk_afft_l(iadd  ,ib)
!fj --------------------
                   end do
                end do
             end do
          end do
       end if
    end do

#ifdef __TIMER_DO__
  call timer_end(136)
#endif

#ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_fft_xz_world)
  call timer_sta(137)
#endif

    call MPI_ALLTOALL(wk_send1, x2z_smax*kimg*iesize, mpi_double_precision,   &
   &                  wk_recv1, x2z_rmax*kimg*iesize, mpi_double_precision,   &
   &                                                  mpi_fft_xz_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Inverse_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 70, ierr)
     endif

#ifdef __TIMER_DO__
  call timer_sta(138)
#endif
#ifdef __TIMER_DO__
  call timer_end(138)
#endif

#else
!else ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_comm)
  call timer_sta(137)
#endif

    icnt_recv = 0
    lrank = myrank + 1
    if (lrank > (nmrank -1)) lrank = 0
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank -1)) lrank = 0
       if ((x2z_recv(1,lrank) /= 0)) then
          call mpi_irecv(wk_recv1(1,x2z_recv(2,lrank)), x2z_recv(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_r(icnt_recv), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Inverse_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 70, ierr)
           endif
          icnt_recv = icnt_recv + 1
       endif
    enddo

    icnt_send = 0
    lrank = myrank
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank -1)) lrank = 0
       if ((x2z_send(1,lrank) /= 0)) then
          call mpi_isend(wk_send1(1,x2z_send(2,lrank)), x2z_send(1,lrank)*kimg*iesize, &
                         mpi_double_precision, lrank, itag, mpi_comm, req_s(icnt_send), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Inverse_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 71, ierr)
           endif
          icnt_send = icnt_send + 1
       endif
    enddo

#ifdef __TIMER_DO__
  call timer_sta(138)
#endif

#ifdef __TIMER_DO__
  call timer_end(138)
#endif

    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Inverse_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 72, ierr)
     endif

    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Inverse_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 73, ierr)
     endif

#endif
!endif ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_end(137)
#endif

#ifdef __TIMER_DO__
  call timer_sta(139)
#endif
    nx = xyz_fft_z(2,1)-xyz_fft_z(1,1)+1
    ny = xyz_fft_z(2,2)-xyz_fft_z(1,2)+1
    nz = xyz_fft_z(2,3)-xyz_fft_z(1,3)+1
    nzp = nz
    if (kimg == 1) then
       do ib = 1, iesize
          do j = 1, ny
             do i = 1, nx/2
                do l = 1, fft_X_z_dim
                   nis = nis_fft_X_z(l)
                   nie = nie_fft_X_z(l)
                   nnz = nie - nis + 1
                   iadd0 = nx/2*ny*nnz*(ib-1)*2
                   do k = nis, nie
                      iadd1 = iadd0 + (i+(j-1)*nx/2+(k-nis)*nx/2*ny)*2
                      iadd  = (k+(i-1)*nzp+(j-1)*nzp*nx/2)*2
#if 0
                      iadd  = (i+(j-1)*nx/2+(k-1)*nx/2*ny)*2
#endif
                      wk_afft_l(iadd-1,ib) = wk_recv1(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv1(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end do
#if 0
! (x,y,z) -> (z,x,y)
       wk_afft_y(1:wk_size,1:iesize) = wk_afft_l(1:wk_size,1:iesize)
       do ib = 1, iesize
          do k = 1, nz
             do j = 1, ny
                do i = 1, nx/2
                   iadd0 = (i+(j-1)*nx/2+(k-1)*nx/2*ny)*2
                   iadd1 = (k+(i-1)*nzp+(j-1)*nzp*nx/2)*2
                   wk_afft_l(iadd1-1,ib) = wk_afft_y(iadd0-1,ib)
                   wk_afft_l(iadd1  ,ib) = wk_afft_y(iadd0  ,ib)
                enddo
             enddo
          enddo
       enddo
#endif
    else
       do ib = 1, iesize
          do j = 1, ny
             do i = 1, nx
                do l = 1, fft_X_z_dim
                   nis = nis_fft_X_z(l)
                   nie = nie_fft_X_z(l)
                   nnz = nie - nis + 1
                   iadd0 = nx*ny*nnz*(ib-1)*2
                   do k = nis, nie
                      iadd1 = iadd0+(i+(j-1)*nx+(k-nis)*nx*ny)*2
                      iadd  = (k+(i-1)*nzp+(j-1)*nzp*nx)*2
#if 0
                      iadd  = (i+(j-1)*nx+(k-1)*nx*ny)*2
#endif
                      wk_afft_l(iadd-1,ib) = wk_recv1(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv1(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end do
! (x,y,z) -> (z,x,y)
#if 0
       wk_afft_y(1:wk_size*kimg,1:iesize) = wk_afft_l(1:wk_size*kimg,1:iesize)
       do ib = 1, iesize
          do k = 0, nz-1
             do i = 0, nx-1
                do j = 0, ny-1
                   wk_afft_l((k+i*nzp+j*nzp*nx+1)*2-1,ib) = wk_afft_y((i+j*nx+k*nx*ny+1)*2-1,ib)
                   wk_afft_l((k+i*nzp+j*nzp*nx+1)*2  ,ib) = wk_afft_y((i+j*nx+k*nx*ny+1)*2  ,ib)
                enddo
             enddo
          enddo
       enddo
#endif
    end if
#ifdef __TIMER_DO__
  call timer_end(139)
#endif

!
! Z-axis (x-y div)
!
#ifdef __TIMER_DO__
  call timer_sta(255)
#endif
    if (kimg == 1) then
       nsize(1:3) = (/nz,nx,(ny/2)/)
       isin(1:3)  = (/1,0,0/)
!OCL INDEPENDENT[ DVMCF2 ]
!OCL PARALLEL_STRONG
       do ib = 1, iesize
          call DVMCF2(wk_afft_l(1,ib),nsize,3,isin,ierr)
       end do
    else
       nsize(1:3) = (/nz,nx,ny/)
       isin(1:3)  = (/-1,0,0/)
!OCL INDEPENDENT[ DVMCF2 ]
!OCL PARALLEL_STRONG
       do ib = 1, iesize
          call DVMCF2(wk_afft_l(1,ib),nsize,3,isin,ierr)
       end do
    end if
#ifdef __TIMER_DO__
  call timer_end(255)
#endif

#ifdef __TIMER_DO__
  call timer_sta(140)
#endif
#ifdef __TIMER_DO__
  call timer_end(140)
#endif

#ifdef __TIMER_DO__
  call timer_sta(141)
#endif
    if (kimg == 1) then
       do ib = 1, iesize
          do j = 1, ny
             do i = 1, nx/2
                do l = 1, fft_Y_z_dim
                   nis = nis_fft_Y_z(l)
                   nie = nie_fft_Y_z(l)
                   nnz = nie - nis + 1
                   iadd0 = nx/2*ny*nnz*(ib-1)*2
                   do k = nis, nie
                      iadd  = (k+(i-1)*nzp+(j-1)*nzp*nx/2)*2
                      iadd1 = iadd0 + (k-nis+1+(i-1)*nnz+(j-1)*nnz*nx/2)*2
                      wk_send2(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send2(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
    else
!OCL PARALLEL_STRONG
!OCL NORECURRENCE
       do ib = 1, iesize
!OCL SERIAL
          do j = 1, ny
!OCL SERIAL
             do i = 1, nx
!OCL SERIAL
                do l = 1, fft_Y_z_dim
                   nis = nis_fft_Y_z(l)
                   nie = nie_fft_Y_z(l)
                   nnz = nie - nis + 1
                   iadd0 = nx*ny*nnz*(ib-1)*2
!OCL SERIAL
                   do k = nis, nie
                      iadd  = (k+(i-1)*nzp+(j-1)*nzp*nx)*2
                      iadd1 = iadd0+(k-nis+1+(i-1)*nnz+(j-1)*nnz*nx)*2
!fj --------------------
!                     do ri = 1, kimg
!                        wk_send2(iadd1+(ri-kimg),l) = wk_afft_l(iadd+(ri-kimg),ib)
!                     end do
                      wk_send2(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send2(iadd1  ,l) = wk_afft_l(iadd  ,ib)
!fj --------------------
                   end do
                end do
             end do
          end do
       end do
    end if
#ifdef __TIMER_DO__
  call timer_end(141)
#endif

#ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_fft_zy_world)
  call timer_sta(142)
#endif

    call MPI_ALLTOALL(wk_send2, z2y_smax*kimg*iesize, mpi_double_precision,   &
   &                  wk_recv2, z2y_rmax*kimg*iesize, mpi_double_precision,   &
   &                                                  mpi_fft_zy_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Inverse_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 76, ierr)
     endif

#ifdef __TIMER_DO__
  call timer_sta(143)
#endif
#ifdef __TIMER_DO__
  call timer_end(143)
#endif

#else
!else ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_comm)
  call timer_sta(142)
#endif

    icnt_recv = 0
    lrank = myrank + 1
    if (lrank > (nmrank-1)) lrank = 0
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank-1)) lrank = 0
       if ((z2y_recv(1,lrank) /= 0)) then
          call mpi_irecv(wk_recv2(1,z2y_recv(2,lrank)), z2y_recv(1,lrank)*kimg*iesize, &
                         mpi_double_precision, lrank, itag, mpi_comm, req_r(icnt_recv), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Inverse_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 76, ierr)
           endif
          icnt_recv = icnt_recv + 1
       endif
    enddo

    icnt_send = 0
    lrank = myrank
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank -1)) lrank = 0
       if ((z2y_send(1,lrank) /= 0)) then
          call mpi_isend(wk_send2(1,z2y_send(2,lrank)), z2y_send(1,lrank)*kimg*iesize, &
                         mpi_double_precision, lrank, itag, mpi_comm, req_s(icnt_send), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Inverse_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 77, ierr)
           endif
          icnt_send = icnt_send + 1
       endif
    enddo

#ifdef __TIMER_DO__
  call timer_sta(143)
#endif
#ifdef __TIMER_DO__
  call timer_end(143)
#endif

    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Inverse_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 78, ierr)
     endif

    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Inverse_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 79, ierr)
     endif

#endif
!endif ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_end(142)
#endif

#ifdef __TIMER_DO__
  call timer_sta(144)
#endif
    nx = xyz_fft_y(2,1)-xyz_fft_y(1,1)+1
    ny = xyz_fft_y(2,2)-xyz_fft_y(1,2)+1
    nz = xyz_fft_y(2,3)-xyz_fft_y(1,3)+1
    nyp = ny

    if (kimg == 1) then
       do ib = 1, iesize
          do k = 1, nz
             do i = 1, nx/2
                do l = 1, fft_Z_y_dim
                   nis = nis_fft_Z_y(l)
                   nie = nie_fft_Z_y(l)
                   nny = nie - nis + 1
                   iadd0 = nx/2*nny*nz*(ib-1)*2
                   do j = nis, nie
                      iadd1 = iadd0 + (k+(i-1)*nz+(j-nis)*nz*nx/2)*2
                      iadd  = (j+(i-1)*nyp+(k-1)*nyp*nx/2)*2
                      wk_afft_y(iadd-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_y(iadd  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end do
    else
       do ib = 1, iesize
          do k = 1, nz
             do i = 1, nx
                do l = 1, fft_Z_y_dim
                   nis = nis_fft_Z_y(l)
                   nie = nie_fft_Z_y(l)
                   nny = nie - nis + 1
                   iadd0 = nx*nny*nz*(ib-1)*2
                   do j = nis, nie
                      iadd1 = iadd0+(k+(i-1)*nz+(j-nis)*nz*nx)*2
                      iadd  = (j+(i-1)*nyp+(k-1)*nyp*nx)*2
                      wk_afft_y(iadd-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_y(iadd  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end do
    end if
#ifdef __TIMER_DO__
  call timer_end(144)
#endif

!
! Y-axis (z-x div)
!
#ifdef __TIMER_DO__
  call timer_sta(256)
#endif
    if (kimg == 1) then
       nsize(1:3) = (/ny,nx/2,nz/)
       isin(1:3)  = (/1,0,0/)
!OCL INDEPENDENT[ DVMCF2 ]
!OCL PARALLEL_STRONG
       do ib = 1, iesize
          call DVMCF2(wk_afft_y(1,ib),nsize,3,isin,ierr)
       end do
    else
       nsize(1:3) = (/ny,nx,nz/)
       isin(1:3)  = (/-1,0,0/)
!OCL INDEPENDENT[ DVMCF2 ]
!OCL PARALLEL_STRONG
       do ib = 1, iesize
          call DVMCF2(wk_afft_y(1,ib),nsize,3,isin,ierr)
       end do
    end if
#ifdef __TIMER_DO__
  call timer_end(256)
#endif

#ifdef __TIMER_DO__
  call timer_sta(145)
#endif
    if (kimg == 1) then
! (x,z,y) <- (y,x,z)
       do ib = 1, iesize
          do j = 1, ny
             do k = 1, nz
                do i = 1, nx/2
                   iadd0 = (j+(i-1)*nyp+(k-1)*nyp*nx/2)*2
                   iadd1 = (i+(k-1)*nx/2+(j-1)*nx/2*nz)*2
                   wk_afft_l(iadd1-1,ib) = wk_afft_y(iadd0-1,ib)
                   wk_afft_l(iadd1  ,ib) = wk_afft_y(iadd0  ,ib)
                enddo
             enddo
          enddo
       enddo
    else
!OCL PARALLEL_STRONG
       do ib = 1, iesize
!OCL SERIAL
          do j = 0, ny-1
!OCL SERIAL
             do k = 0, nz-1
!OCL SERIAL
                do i = 0, nx-1
                   wk_afft_l((i+k*nx+j*nx*nz+1)*2-1,ib) = wk_afft_y((j+i*nyp+k*nyp*nx+1)*2-1,ib)
                   wk_afft_l((i+k*nx+j*nx*nz+1)*2  ,ib) = wk_afft_y((j+i*nyp+k*nyp*nx+1)*2  ,ib)
                enddo
             enddo
          enddo
       enddo
    end if
#ifdef __TIMER_DO__
  call timer_end(145)
#endif

#ifdef __TIMER_SUB__
  call timer_end(1498)
#endif

#ifdef __TIMER_SUB__
  call timer_end(1499)
#endif

1000 continue

!F  if (allocated(req_r))       deallocate(req_r)
!F  if (allocated(req_s))       deallocate(req_s)
!F  if (allocated(sta_r))       deallocate(sta_r)
!F  if (allocated(sta_s))       deallocate(sta_s)
!F  if (allocated(wk_recvcnt))  deallocate(wk_recvcnt)
!F  if (allocated(wk_sendcnt))  deallocate(wk_sendcnt)
    if (allocated(wk_mp_fft_y)) deallocate(wk_mp_fft_y)
!F  if (allocated(wk_afft_y))   deallocate(wk_afft_y)

#ifdef __TIMER_SUB__
    call timer_end(106)
#endif
  end subroutine m_FFT_Inverse_3D
!------------------------------------------------------------------------------
#endif
!endif ifdef FFT_USE_SSL2

!------------------------------------------------------------------------------
  subroutine m_FFT_Vlocal_W_3D(afft_l,bfft_l,lsize,ibsize,nfft_l)
    integer, intent(in) :: lsize, ibsize, nfft_l
    real(kind=DP), intent(in)   , dimension(lsize*kimg)        :: afft_l
    real(kind=DP), intent(inout), dimension(lsize*kimg,ibsize) :: bfft_l
    integer i,j
#ifdef __FAPP__
    call fapp_start('fft_vlocal_w',1,1)
#endif
#ifdef __TIMER_SUB__
    call timer_sta(107)
#endif

#ifdef __TIMER_DO__
  call timer_sta(151)
#endif
    if (kimg == 1) then
!OCL NOFLTLD
       do j = 1, ibsize
          do i = 1, nfft_l, 2
             bfft_l(i  ,j) = afft_l(i) * bfft_l(i,  j)
             bfft_l(i+1,j) = afft_l(i) * bfft_l(i+1,j)
          end do
       end do
    else
     if (ibsize > 1) then
!$OMP PARALLEL DO PRIVATE(i,j)
       do j = 1, ibsize
          do i = 1, nfft_l
             bfft_l(i*2-1,j) = afft_l(i*2-1) * bfft_l(i*2-1,j)
             bfft_l(i*2  ,j) = afft_l(i*2-1) * bfft_l(i*2  ,j)
          end do
       end do
!$OMP END PARALLEL DO
     else
       do j = 1, ibsize
!$OMP PARALLEL DO PRIVATE(i)
          do i = 1, nfft_l
             bfft_l(i*2-1,j) = afft_l(i*2-1) * bfft_l(i*2-1,j)
             bfft_l(i*2  ,j) = afft_l(i*2-1) * bfft_l(i*2  ,j)
          end do
!$OMP END PARALLEL DO
       end do
     end if
    end if
#ifdef __TIMER_DO__
  call timer_end(151)
#endif

#ifdef __TIMER_SUB__
    call timer_end(107)
#endif
#ifdef __FAPP__
    call fapp_stop('fft_vlocal_w',1,1)
#endif
  end subroutine m_FFT_Vlocal_W_3D
!------------------------------------------------------------------------------

#ifdef MPI_FFTW

!------------------------------------------------------------------------------
  subroutine m_FFT_Vlocal_W_mpifftw(afft_l,lsize,nfft_l)
    integer, intent(in) :: lsize, nfft_l
    real(kind=DP), intent(in), dimension(lsize*kimg)  :: afft_l
!    real(kind=DP), intent(inout), dimension(lsize*kimg,ibsize) :: bfft_l
    integer i,j,ix,iy,iz
    integer :: i1,mm,mx,my,mz
    integer(C_INTPTR_T) :: alloc_local, lx, ly, lz, local_n, local_n_offset, lxh
    integer :: id_sname=-1
    call tstatc0_begin('m_FFT_Vlocal_W_mpifftw ',id_sname)
    lx = fft_box_size_WF(1,0)
    ly = fft_box_size_WF(2,0)
    lz = fft_box_size_WF(3,0)
    if(kimg==2) then
      alloc_local = fftw_mpi_local_size_3d(ly,lz,lx,mpi_ke_world,local_n,local_n_offset)
      do iy=1,local_n
        do iz=1,lz
          do ix=1,lx
            i1 = (iy-1)*lx*lz+(iz-1)*lx+ix
            afft_mpifftw(ix,iz,iy) = afft_l(i1*2-1)*afft_mpifftw(ix,iz,iy)
          enddo
        enddo
      enddo
    else
      lxh = lx/2
      alloc_local = fftw_mpi_local_size_3d(ly,lz,lxh,mpi_ke_world,local_n,local_n_offset)
      do iy=1,local_n
        do iz=1,lz
          do ix=1,lxh
            i1 = (iy-1)*lxh*lz+(iz-1)*lxh+ix
            afft_mpifftw_kimg1(ix,iz,iy) = afft_l(i1*2-1)*afft_mpifftw_kimg1(ix,iz,iy)
          enddo
        enddo
      enddo
    endif
    call tstatc0_end(id_sname)
  end subroutine m_FFT_Vlocal_W_mpifftw
!------------------------------------------------------------------------------

  subroutine m_FFT_Vlocal_W_mpifftw3d(afft_l_3d,lx,local_n,lz)
    integer(C_INTPTR_T), intent(in) :: lx,local_n,lz
    real(kind=DP), intent(in), dimension(lx,lz,local_n)  :: afft_l_3d
    integer :: ix, iy, iz
    integer :: id_sname = -1
    call tstatc0_begin('m_FFT_Vlocal_W_mpifftw ',id_sname)
    if(kimg==2) then
      afft_mpifftw = afft_l_3d * afft_mpifftw
    else
      afft_mpifftw_kimg1 = afft_l_3d * afft_mpifftw_kimg1
    endif
    call tstatc0_end(id_sname)
  end subroutine m_FFT_Vlocal_W_mpifftw3d
#endif

!------------------------------------------------------------------------------
  subroutine m_FFT_W_Vlocal_W_3D(electron_or_positron,afft,bfft,lsize,ibesize,s)
    integer, intent(in) ::                             electron_or_positron
    integer, intent(in) ::                             lsize, ibesize
    real(kind=DP), intent(in),    dimension(lsize*kimg) :: afft
    real(kind=DP), intent(inout), dimension(lsize*kimg,ibesize) :: bfft
    real(kind=DP), intent(out),   dimension(ibesize)    :: s

    integer i, id, j, k, ip, nl, nm, nn, nd2, nd3, nlh, idh, ib, ierr, nx, ny, nz
    real(kind=DP) :: prd
#ifdef __TIMER_SUB__
    call timer_sta(108)
#endif

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

    if (np_fft_y == 0) then
       nx = 0
       ny = 0
       nz = 0
    else
       nx = xyz_fft_y(2,1)-xyz_fft_y(1,1)+1
       ny = xyz_fft_y(2,2)-xyz_fft_y(1,2)+1
       nz = xyz_fft_y(2,3)-xyz_fft_y(1,3)+1
    endif

    s = 0.d0
    if(kimg == 1) then

       do ib = 1, ibesize
#ifdef __TIMER_DO__
  call timer_sta(152)
#endif
          do i = 1, nx*nz*ny-1, 2
             bfft(i,ib) = afft(i)*(bfft(i,ib)**2 + bfft(i+1,ib)**2)
          end do
#ifdef __TIMER_DO__
  call timer_end(152)
#endif
!         nlh = nl/2
#ifdef __TIMER_DO__
  call timer_sta(153)
#endif
          do i = 1, nx*nz*ny-1, 2
             s(ib) = s(ib) + bfft(i,ib)
          end do
#ifdef __TIMER_DO__
  call timer_end(153)
#endif
          s(ib) = s(ib) + s(ib)
#ifdef __TIMER_DO__
  call timer_sta(154)
#endif
          if(xyz_fft_y(1,1) == 1) then

#ifdef FFT_USE_SSL2
             do i = 1, nx*nz*ny-1, nx
                s(ib) = s(ib) - bfft(i,ib)
             end do
#else
#ifdef FFT_FFTW_OLD
             do i = 1, nx*nz*ny-1, nx
                s(ib) = s(ib) - bfft(i,ib)
             end do
#else
             do j = 1, nz
                do i = 1, ny*2-1, 2
                   s(ib) = s(ib) - bfft(i+(j-1)*ny*2*nx/2,ib)
                end do
             end do
#endif
#endif

          endif
#ifdef __TIMER_DO__
  call timer_end(154)
#endif
#ifdef __TIMER_DO__
  call timer_sta(155)
#endif
          if(xyz_fft_y(2,1) == id) then

#ifdef FFT_USE_SSL2
             do i = nx-1, nx*nz*ny-1, nx
                s(ib) = s(ib) - bfft(i,ib)
             enddo
#else
#ifdef FFT_FFTW_OLD
             do i = nx-1, nx*nz*ny-1, nx
                s(ib) = s(ib) - bfft(i,ib)
             enddo
#else
             do j = 1, nz
                do i = 1, ny*2-1, 2
                   s(ib) = s(ib) - bfft(i+j*ny*2*nx/2-ny*2,ib)
                enddo
             enddo
#endif
#endif

          endif
#ifdef __TIMER_DO__
  call timer_end(155)
#endif
       end do

    else if(kimg == 2) then
#ifdef __TIMER_DO__
  call timer_sta(156)
#endif
       do ib = 1, ibesize
          do i = 1, nx*nz*ny*kimg-1, 2
             s(ib) = s(ib) + afft(i)*(bfft(i,ib)**2 + bfft(i+1,ib)**2)
          end do
       end do
#ifdef __TIMER_DO__
  call timer_end(156)
#endif
    end if

!   call mpi_allreduce(s, s_mpi, ibesize, mpi_double_precision, mpi_sum, mpi_ke_world, ierr)

!   prd = product(fft_box_size_WF(1:3,1))
!   do ib = 1, ibesize
!      eg(ib) = s_mpi(ib) / prd
!   enddo
#ifdef __TIMER_SUB__
    call timer_end(108)
#endif
  end subroutine m_FFT_W_Vlocal_W_3D
!------------------------------------------------------------------------------

#ifndef FFT_USE_SSL2
!------------------------------------------------------------------------------
  subroutine m_FFT_CD_Direct_3D(nfout,wk_afft_l, wk_size, iesize)
    integer, intent(in) :: nfout, wk_size, iesize
    real(kind=DP), dimension(wk_size*kimg,iesize), intent(inout) :: wk_afft_l

    real(kind=DP), allocatable, dimension(:) :: wk_allfft
    integer,      allocatable, dimension(:)   :: wk_recvdsp
    integer, parameter :: FFTW_MEASURE=0
    integer, parameter :: FFTW_ESTIMATE=64
    integer, parameter :: FFTW_RANK=1
    integer            :: FFTW_FLAG
    integer,dimension(FFTW_RANK) :: NFFTW3, NEMBED, NEREAL
    integer :: nx, ny, nz, i, j, k, l, ri, ib, iaddr, lrank
    integer :: irank, itag, icnt_send, icnt_recv
    integer,save, allocatable, dimension(:)   :: req_r, req_s
    integer,save, allocatable, dimension(:,:) :: sta_r, sta_s
    integer,save, allocatable, dimension(:)   :: wk_recvcnt, wk_sendcnt
!   integer,      allocatable, dimension(:)   :: wk_recvrank, wk_sendrank
!   real(kind=DP),allocatable, dimension(:,:) :: wk_recv, wk_send
    real(kind=DP),allocatable, dimension(:,:),save :: wk_recv1, wk_send1
    real(kind=DP),allocatable, dimension(:,:),save :: wk_recv2, wk_send2
    integer,      allocatable, dimension(:)   :: wk_mp_fft_y
    real(kind=DP),allocatable, dimension(:,:),save :: wk_afft_y
    real(kind=DP),allocatable, dimension(:,:,:)   :: wk_gather
    integer :: mpi_comm, myrank, nmrank
    integer(kind=8) :: plan
    integer :: max_x, max_y, max_z
    integer,save :: max_elm = 0, savesize = 0

    integer,save, allocatable, dimension(:,:) :: y2z_recv, y2z_send, z2x_recv, z2x_send
    integer,save :: y2z_rrank, y2z_srank, z2x_rrank, z2x_srank
    integer,save :: y2z_rmax, y2z_smax, y2z_srmax, z2x_rmax, z2x_smax, z2x_srmax

    integer ::  nnx, nny, nnz, iadd0, iadd1, nis, nie, jadd, kadd
#ifdef FFT_ALLTOALL
    integer,dimension(2) :: max_send,max_recv
#endif
#ifdef __TIMER_SUB__
    call mpi_barrier(mpi_ke_world,ierr)
    call timer_sta(109)
#endif

    if (kimg == 1) then
       FFTW_FLAG = +1      ! FFTW_BACKWARD
    else
       FFTW_FLAG = -1      ! FFTW_FORWARD
    endif
    plan = 0

!fj!$$#ifdef CD_FFT_ALL
!fj!$$    mpi_comm = MPI_CommGroup
!fj!$$    myrank = mype
!fj!$$    nmrank = npes
!fj!$$#else
    mpi_comm = mpi_ke_world
    myrank = myrank_g
    nmrank = nrank_g
!fj!$$#endif
    itag = 10

    if (np_fftcd_x==0)then
       nx = 0
       nx = 0
       nz = 0
       go to 1000
    endif

    if (firstcall_cd_direct_3d) then
       savesize = 0
#ifdef __TIMER_DO__
  call timer_sta(157)
#endif
       max_x = maxval(nel_fftcd_x(:))
       max_y = maxval(nel_fftcd_y(:))
       max_z = maxval(nel_fftcd_z(:))
       max_elm = max(max_x,max_y,max_z)

       if(allocated(req_r)) deallocate(req_r)
       if(allocated(req_s)) deallocate(req_s)
       if(allocated(sta_r)) deallocate(sta_r)
       if(allocated(sta_s)) deallocate(sta_s)
       if(allocated(wk_recvcnt)) deallocate(wk_recvcnt)
       if(allocated(wk_sendcnt)) deallocate(wk_sendcnt)
       allocate(req_r(0:nmrank-1), stat=ierr)
       allocate(req_s(0:nmrank-1), stat=ierr)
       allocate(sta_r(MPI_STATUS_SIZE,0:nmrank-1), stat=ierr)
       allocate(sta_s(MPI_STATUS_SIZE,0:nmrank-1), stat=ierr)
       allocate(wk_recvcnt(0:nmrank-1), stat=ierr)
       allocate(wk_sendcnt(0:nmrank-1), stat=ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_FFT_CD_Direct_3D :  Not allocate '
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 51, ierr)
        endif

       if(allocated(y2z_recv)) deallocate(y2z_recv)
       if(allocated(y2z_send)) deallocate(y2z_send)
       allocate(y2z_recv(2,0:nmrank-1))
       allocate(y2z_send(2,0:nmrank-1))
       y2z_send = 0
       y2z_recv = 0
       wk_recvcnt = 0
       wk_sendcnt = 0
       do i = 1, nel_fftcd_z(myrank)
          irank = map_fftcd_y(mp_fftcd_z(i)) - 1
          wk_recvcnt(irank) = wk_recvcnt(irank) + 1
       enddo
       do i = 1, nel_fftcd_y(myrank)
          irank = map_fftcd_z(mp_fftcd_y(i)) - 1
          wk_sendcnt(irank) = wk_sendcnt(irank) + 1
       enddo
       k = 0
       do i = 0, nmrank - 1
          if(wk_recvcnt(i) /= 0) then
             k = k + 1
             y2z_recv(1,i) = wk_recvcnt(i)
             y2z_recv(2,i) = k
          endif
       enddo
       y2z_rrank = k
       k = 0
       do i = 0, nmrank - 1
          if(wk_sendcnt(i) /= 0) then
             k = k + 1
             y2z_send(1,i) = wk_sendcnt(i)
             y2z_send(2,i) = k
          endif
       enddo
       y2z_srank = k
       y2z_rmax = maxval(wk_recvcnt)
       y2z_smax = maxval(wk_sendcnt)
       y2z_srmax = max(y2z_rmax,y2z_smax)

       if(allocated(z2x_recv)) deallocate(z2x_recv)
       if(allocated(z2x_send)) deallocate(z2x_send)
       allocate(z2x_recv(2,0:nmrank-1))
       allocate(z2x_send(2,0:nmrank-1))
       z2x_send = 0
       z2x_recv = 0
       wk_recvcnt = 0
       wk_sendcnt = 0
       do i = 1, nel_fftcd_x(myrank)
          irank = map_fftcd_z(mp_fftcd_x(i)) - 1
          wk_recvcnt(irank) = wk_recvcnt(irank) + 1
       enddo
       do i = 1, nel_fftcd_z(myrank)
          irank = map_fftcd_x(mp_fftcd_z(i)) - 1
          wk_sendcnt(irank) = wk_sendcnt(irank) + 1
       enddo
       k = 0
       do i = 0, nmrank - 1
          if(wk_recvcnt(i) /= 0) then
             k = k + 1
             z2x_recv(1,i) = wk_recvcnt(i)
             z2x_recv(2,i) = k
          endif
       enddo
       z2x_rrank = k
       k = 0
      do i = 0, nmrank - 1
          if(wk_sendcnt(i) /= 0) then
             k = k + 1
             z2x_send(1,i) = wk_sendcnt(i)
             z2x_send(2,i) = k
          endif
       enddo
       z2x_srank = k
       z2x_rmax = maxval(wk_recvcnt)
       z2x_smax = maxval(wk_sendcnt)
       z2x_srmax = max(z2x_rmax,z2x_smax)

#ifdef FFT_ALLTOALL
       max_send(1) = y2z_rmax
       max_send(2) = y2z_smax
       call mpi_allreduce(max_send,max_recv,2,mpi_integer,mpi_max,mpi_fftcd_zy_world,ierr)
       y2z_rmax = max_recv(1)
       y2z_smax = max_recv(2)
       max_send(1) = z2x_rmax
       max_send(2) = z2x_smax
       call mpi_allreduce(max_send,max_recv,2,mpi_integer,mpi_max,mpi_fftcd_xz_world,ierr)
       z2x_rmax = max_recv(1)
       z2x_smax = max_recv(2)
#endif

       firstcall_cd_direct_3d = .false.
#ifdef __TIMER_DO__
  call timer_end(157)
#endif
    endif

#ifdef FFT_ALLTOALL
    if (iesize /= savesize) then
#else
    if (iesize > savesize) then
#endif
#ifdef __TIMER_DO__
  call timer_sta(158)
#endif
       if (allocated(wk_recv1)) deallocate(wk_recv1)
       if (allocated(wk_send1)) deallocate(wk_send1)
       if (allocated(wk_recv2)) deallocate(wk_recv2)
       if (allocated(wk_send2)) deallocate(wk_send2)
       if (allocated(wk_afft_y)) deallocate(wk_afft_y)
       allocate(wk_recv1(y2z_rmax*kimg*iesize,y2z_rrank), stat=ierr)
       allocate(wk_send1(y2z_smax*kimg*iesize,y2z_srank), stat=ierr)
       allocate(wk_recv2(z2x_rmax*kimg*iesize,z2x_rrank), stat=ierr)
       allocate(wk_send2(z2x_smax*kimg*iesize,z2x_srank), stat=ierr)
       allocate(wk_afft_y(max_elm*kimg,iesize) ,stat=ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_FFT_CD_Direct_3D :  Not allocate '
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 52, ierr)
        endif
       savesize = iesize
#ifdef __TIMER_DO__
  call timer_end(158)
#endif
    end if
!
! Y-axis (z-x div)
!
    nx = xyz_fftcd_y(2,1)-xyz_fftcd_y(1,1)+1
    ny = xyz_fftcd_y(2,2)-xyz_fftcd_y(1,2)+1
    nz = xyz_fftcd_y(2,3)-xyz_fftcd_y(1,3)+1

!XX!   allocate(wk_afft_y(nx*ny*nz*kimg,iesize) ,stat=ierr)
!XX!    if (ierr /= 0) then
!XX!       write(nfout,*)' m_FFT_CD_Direct_3D :  Not allocate '
!XX!       call flush(nfout)
!XX!       call mpi_abort(mpi_comm_world, 53, ierr)
!XX!    endif
!XX!   do ib = 1, iesize
!XX!      do k = 0, nz-1
!XX!         do j = 0, ny-1
!XX!            do i = 0, nx-1
!XX!   do ri = 1, kimg
!XX!               wk_afft_y((i+k*nx+j*nx*nz+1)*kimg+(ri-kimg),ib) = wk_afft_l((i+j*nx+k*nx*ny+1)*kimg+(ri-kimg),ib)
!XX!   enddo
!XX!            enddo
!XX!         enddo
!XX!      enddo
!XX!   enddo
!XX!   do ib = 1, iesize
!XX!      do k = 0, nz-1
!XX!         do j = 0, ny-1
!XX!            do i = 0, nx-1
!XX!   do ri = 1, kimg
!XX!               wk_afft_l((i+k*nx+j*nx*nz+1)*kimg+(ri-kimg),ib) = wk_afft_y((i+k*nx+j*nx*nz+1)*kimg+(ri-kimg),ib)
!XX!   enddo
!XX!            enddo
!XX!         enddo
!XX!      enddo
!XX!   enddo
!XX!   deallocate(wk_afft_y)

    if(kimg==1) then
       NFFTW3(1) = fft_box_size_CD_3D(2,1)
       NEMBED(1) = fft_box_size_CD_3D(2,1)
       call dfftw_plan_many_dft    (plan, FFTW_RANK, NFFTW3, nz*nx/2, &
      &                             wk_afft_l, NEMBED, nz*nx/2, 1,    &
      &                             wk_afft_l, NEMBED, nz*nx/2, 1,    &
      &                             FFTW_FLAG, FFTW_ESTIMATE)
    else
       NFFTW3(1) = ny
       NFFTW3(1) = fft_box_size_CD_3D(2,1)
       NEMBED(1) = ny
       NEMBED(1) = fft_box_size_CD_3D(2,0)
       call dfftw_plan_many_dft    (plan, FFTW_RANK, NFFTW3, nz*nx,   &
      &                             wk_afft_l, NEMBED, nz*nx, 1,      &
      &                             wk_afft_l, NEMBED, nz*nx, 1,      &
      &                             FFTW_FLAG, FFTW_ESTIMATE)
    end if
#ifdef __TIMER_DO__
  call timer_sta(257)
#endif
!OCL INDEPENDENT[dfftw_execute_dft]
!OCL PARALLEL_STRONG
    do ib = 1, iesize
       call dfftw_execute_dft(plan,wk_afft_l(1,ib),wk_afft_l(1,ib))
    end do
#ifdef __TIMER_DO__
  call timer_end(257)
#endif
    call dfftw_destroy_plan(plan)

#ifdef __TIMER_DO__
  call timer_sta(159)
#endif
#ifdef __TIMER_DO__
  call timer_end(159)
#endif

#ifdef __TIMER_DO__
  call timer_sta(160)
#endif
    if (kimg == 1) then
       do ib = 1, iesize
          do k = 0, nz-1
             do j = 0, ny-1
                do i = 0, nx-1
                   do ri = 1, kimg
                      wk_afft_y((i+j*nx+k*nx*ny+1)*kimg+(ri-kimg),ib) = wk_afft_l((i+k*nx+j*nx*nz+1)*kimg+(ri-kimg),ib)
                   enddo
                enddo
             enddo
          enddo
       enddo
    else
       do ib = 1, iesize
          do k = 0, nz-1
             do j = 0, ny-1
                do i = 0, nx-1
                   wk_afft_y((i+j*nx+k*nx*ny+1)*2-1,ib) = wk_afft_l((i+k*nx+j*nx*nz+1)*2-1,ib)
                   wk_afft_y((i+j*nx+k*nx*ny+1)*2  ,ib) = wk_afft_l((i+k*nx+j*nx*nz+1)*2  ,ib)
                enddo
             enddo
          enddo
       enddo
    end if
#ifdef __TIMER_DO__
  call timer_end(160)
#endif

#ifdef __TIMER_DO__
  call timer_sta(161)
#endif

    nx = xyz_fftcd_y(2,1)-xyz_fftcd_y(1,1)+1
    ny = xyz_fftcd_y(2,2)-xyz_fftcd_y(1,2)+1
    nz = xyz_fftcd_y(2,3)-xyz_fftcd_y(1,3)+1
    do l = 1, fftcd_Z_y_dim
       nis = nis_fftcd_Z_y(l)
       nie = nie_fftcd_Z_y(l)
       nny = nie - nis + 1
       if (kimg == 1) then
          do ib = 1, iesize
             iadd0 = nx*nny*nz*(ib-1)
             do k = 1, nz
                do j = nis, nie
                   do i = 1, nx
                      wk_send1(iadd0+(i+(j-nis)*nx+(k-1)*nx*nny),l) = wk_afft_y(i+(j-1)*nx+(k-1)*nx*ny,ib)
                   end do
                end do
             end do
          end do
       else
          do ib = 1, iesize
             iadd0 = nx*nny*nz*(ib-1)*2
             do k = 1, nz
                kadd = (k-1)*nx*nny
                do j = nis, nie
                   jadd = (j-nis)*nx+kadd
                   do i = 1, nx
                      iadd1 = iadd0+(i+jadd)*2
                      iaddr = (i+(j-1)*nx+(k-1)*nx*ny)*2
                      wk_send1(iadd1-1,l) = wk_afft_y(iaddr-1,ib)
                      wk_send1(iadd1  ,l) = wk_afft_y(iaddr  ,ib)
                   end do
                end do
             end do
          end do
       end if
    end do

#ifdef __TIMER_DO__
  call timer_end(161)
#endif

#ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_fftcd_zy_world)
  call timer_sta(162)
#endif

    call MPI_ALLTOALL(wk_send1, y2z_smax*kimg*iesize, mpi_double_precision,   &
   &                  wk_recv1, y2z_rmax*kimg*iesize, mpi_double_precision,   &
   &                                                  mpi_fftcd_zy_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Inverse_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 55, ierr)
     endif

#ifdef __TIMER_DO__
  call timer_sta(163)
#endif
#ifdef __TIMER_DO__
  call timer_end(163)
#endif

#else
!else ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_comm)
  call timer_sta(162)
#endif

    icnt_recv = 0
    do lrank = 0, nmrank - 1
       if ((lrank /= myrank) .and. (y2z_recv(1,lrank) /= 0)) then
          call mpi_irecv(wk_recv1(1,y2z_recv(2,lrank)), y2z_recv(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_r(icnt_recv), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_CD_Direct_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 55, ierr)
           endif
          icnt_recv = icnt_recv + 1
       endif
    enddo

    icnt_send = 0
    do lrank = 0, nmrank - 1
       if ((lrank /= myrank) .and. (y2z_send(1,lrank) /= 0)) then
          call mpi_isend(wk_send1(1,y2z_send(2,lrank)), y2z_send(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_s(icnt_send), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_CD_Direct_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 56, ierr)
           endif
          icnt_send = icnt_send + 1
       endif
    enddo

#ifdef __TIMER_DO__
  call timer_sta(163)
#endif
    do i = 1, y2z_recv(1,myrank)*kimg*iesize
       wk_recv1(i,y2z_recv(2,myrank)) = wk_send1(i,y2z_send(2,myrank))
    enddo
#ifdef __TIMER_DO__
  call timer_end(163)
#endif

    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_CD_Direct_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 57, ierr)
     endif

    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_CD_Direct_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 58, ierr)
     endif

#endif
!endif ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_end(162)
#endif

#ifdef __TIMER_DO__
  call timer_sta(164)
#endif
    nx = xyz_fftcd_z(2,1)-xyz_fftcd_z(1,1)+1
    ny = xyz_fftcd_z(2,2)-xyz_fftcd_z(1,2)+1
    nz = xyz_fftcd_z(2,3)-xyz_fftcd_z(1,3)+1
    do l = 1, fftcd_Y_z_dim
       nis = nis_fftcd_Y_z(l)
       nie = nie_fftcd_Y_z(l)
       nnz = nie - nis + 1
       if (kimg == 1) then
          do ib = 1, iesize
             iadd0 = nx*ny*nnz*(ib-1)
             do k = nis, nie
                do j = 1, ny
                   do i = 1, nx
                      wk_afft_l(i+(j-1)*nx+(k-1)*nx*ny,ib) = wk_recv1(iadd0+i+(j-1)*nx+(k-nis)*nx*ny,l)
                   end do
                end do
             end do
          end do
       else
          do ib = 1, iesize
             iadd0 = nx*ny*nnz*(ib-1)*2
             do k = nis, nie
                do j = 1, ny
                   do i = 1, nx
                      iadd1 = iadd0+(i+(j-1)*nx+(k-nis)*nx*ny)*2
                      iaddr = (i+(j-1)*nx+(k-1)*nx*ny)*2
                      wk_afft_l(iaddr-1,ib) = wk_recv1(iadd1-1,l)
                      wk_afft_l(iaddr  ,ib) = wk_recv1(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end if
    end do
#ifdef __TIMER_DO__
  call timer_end(164)
#endif

!
! Z-axis (x-y div)
!
    nx = xyz_fftcd_z(2,1)-xyz_fftcd_z(1,1)+1
    ny = xyz_fftcd_z(2,2)-xyz_fftcd_z(1,2)+1
    nz = xyz_fftcd_z(2,3)-xyz_fftcd_z(1,3)+1
    if(kimg==1) then
       NFFTW3(1) = fft_box_size_CD_3D(3,1)
       NEMBED(1) = fft_box_size_CD_3D(3,1)
       call dfftw_plan_many_dft    (plan, FFTW_RANK, NFFTW3, nx*ny/2, &
      &                             wk_afft_l, NEMBED, nx*ny/2, 1,     &
      &                             wk_afft_l, NEMBED, nx*ny/2, 1,     &
      &                             FFTW_FLAG, FFTW_ESTIMATE)
    else
       NFFTW3(1) = nz
       NFFTW3(1) = fft_box_size_CD_3D(3,1)
       NEMBED(1) = nz
       NEMBED(1) = fft_box_size_CD_3D(3,0)
       call dfftw_plan_many_dft    (plan, FFTW_RANK, NFFTW3, nx*ny,   &
      &                             wk_afft_l, NEMBED, nx*ny, 1,       &
      &                             wk_afft_l, NEMBED, nx*ny, 1,       &
      &                             FFTW_FLAG, FFTW_ESTIMATE)
    endif
#ifdef __TIMER_DO__
  call timer_sta(258)
#endif
!OCL INDEPENDENT[dfftw_execute_dft]
!OCL PARALLEL_STRONG
    do ib = 1, iesize
       call dfftw_execute_dft(plan,wk_afft_l(1,ib),wk_afft_l(1,ib))
    end do
#ifdef __TIMER_DO__
  call timer_end(258)
#endif
    call dfftw_destroy_plan(plan)

#ifdef __TIMER_DO__
  call timer_sta(165)
#endif
    nx = xyz_fftcd_z(2,1)-xyz_fftcd_z(1,1)+1
    ny = xyz_fftcd_z(2,2)-xyz_fftcd_z(1,2)+1
    nz = xyz_fftcd_z(2,3)-xyz_fftcd_z(1,3)+1
    do l = 1, fftcd_X_z_dim
       nis = nis_fftcd_X_z(l)
       nie = nie_fftcd_X_z(l)
       nnz = nie - nis + 1
       if (kimg == 1) then
          do ib = 1, iesize
             iadd0 = nx*ny*nnz*(ib-1)
             do k = nis, nie
                do j = 1, ny
                   do i = 1, nx
                      wk_send2(iadd0+(i+(j-1)*nx+(k-nis)*nx*ny),l) = wk_afft_l(i+(j-1)*nx+(k-1)*nx*ny,ib)
                   end do
                end do
             end do
          end do
       else
          do ib = 1, iesize
             iadd0 = nx*ny*nnz*(ib-1)*2
             do k = nis, nie
                jadd = (k-nis)*nx*ny
                do j = 1, ny
                   kadd = (j-1)*nx+jadd
                   do i = 1, nx
                      iadd1 = iadd0+(i+kadd)*2
                      iaddr  = (i+(j-1)*nx+(k-1)*nx*ny)*2
                      wk_send2(iadd1-1,l) = wk_afft_l(iaddr-1,ib)
                      wk_send2(iadd1  ,l) = wk_afft_l(iaddr  ,ib)
                   end do
                end do
             end do
          end do
       end if
    end do

#ifdef __TIMER_DO__
  call timer_end(165)
#endif

#ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_fftcd_xz_world)
  call timer_sta(166)
#endif

    call MPI_ALLTOALL(wk_send2, z2x_smax*kimg*iesize, mpi_double_precision,   &
   &                  wk_recv2, z2x_rmax*kimg*iesize, mpi_double_precision,   &
   &                                                  mpi_fftcd_xz_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_CD_Direct_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 60, ierr)
     endif

#ifdef __TIMER_DO__
  call timer_sta(167)
#endif
#ifdef __TIMER_DO__
  call timer_end(167)
#endif

#else
!else ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_comm)
  call timer_sta(166)
#endif

    icnt_recv = 0
    do lrank = 0, nmrank - 1
       if ((lrank /= myrank) .and. (z2x_recv(1,lrank) /= 0)) then
          call mpi_irecv(wk_recv2(1,z2x_recv(2,lrank)), z2x_recv(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_r(icnt_recv), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_CD_Direct_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 60, ierr)
           endif
          icnt_recv = icnt_recv + 1
       endif
    enddo

    icnt_send = 0
    do lrank = 0, nmrank - 1
       if ((lrank /= myrank) .and. (z2x_send(1,lrank) /= 0)) then
          call mpi_isend(wk_send2(1,z2x_send(2,lrank)), z2x_send(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_s(icnt_send), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_CD_Direct_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 61, ierr)
           endif
          icnt_send = icnt_send + 1
       endif
    enddo

#ifdef __TIMER_DO__
  call timer_sta(167)
#endif
    do i = 1, z2x_recv(1,myrank)*kimg*iesize
       wk_recv2(i,z2x_recv(2,myrank)) = wk_send2(i,z2x_recv(2,myrank))
    enddo
#ifdef __TIMER_DO__
  call timer_end(167)
#endif

    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_CD_Direct_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 62, ierr)
     endif

    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_CD_Direct_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 63, ierr)
     endif

#endif
!endif ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_end(166)
#endif

#ifdef __TIMER_DO__
  call timer_sta(168)
#endif
    nx = xyz_fftcd_x(2,1)-xyz_fftcd_x(1,1)+1
    ny = xyz_fftcd_x(2,2)-xyz_fftcd_x(1,2)+1
    nz = xyz_fftcd_x(2,3)-xyz_fftcd_x(1,3)+1
    do l = 1, fftcd_Z_x_dim
       nis = nis_fftcd_Z_x(l)
       nie = nie_fftcd_Z_x(l)
       nnx = nie - nis + 1
       if (kimg == 1) then
          do ib = 1, iesize
             iadd0 = nnx*ny*nz*(ib-1)
             do k = 1, nz
                do j = 1, ny
                   do i = nis, nie
                      wk_afft_l(i+(j-1)*nx+(k-1)*nx*ny,ib) = wk_recv2(iadd0+i-nis+1+(j-1)*nnx+(k-1)*nnx*ny,l)
                   end do
                end do
             end do
          end do
       else
          do ib = 1, iesize
             iadd0 = nnx*ny*nz*(ib-1)*2
             do k = 1, nz
                do j = 1, ny
                   do i = nis, nie
                      iadd1 = iadd0+(i-nis+1+(j-1)*nnx+(k-1)*nnx*ny)*2
                      iaddr = (i+(j-1)*nx+(k-1)*nx*ny)*2
                      wk_afft_l(iaddr-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_l(iaddr  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end if
    end do
#ifdef __TIMER_DO__
  call timer_end(168)
#endif

!
! X-axis (y-z div)
!
    nx = xyz_fftcd_x(2,1)-xyz_fftcd_x(1,1)+1
    ny = xyz_fftcd_x(2,2)-xyz_fftcd_x(1,2)+1
    nz = xyz_fftcd_x(2,3)-xyz_fftcd_x(1,3)+1
    if(kimg==1) then
       NFFTW3(1) = nx - 2
       NEMBED(1) = nx
       NEREAL(1) = nx
       call dfftw_plan_many_dft_c2r(plan, FFTW_RANK, NFFTW3, ny*nz,   &
      &                             wk_afft_l, NEMBED, 1, nx/2,        &
      &                             wk_afft_l, NEREAL, 1, nx,          &
      &                             FFTW_ESTIMATE)
    else
       NFFTW3(1) = fft_box_size_CD_3D(1,1)
       NEMBED(1) = nx
       NEMBED(1) = fft_box_size_CD_3D(1,1)
       call dfftw_plan_many_dft    (plan, FFTW_RANK, NFFTW3, ny*nz,   &
      &                             wk_afft_l, NEMBED, 1, nx,          &
      &                             wk_afft_l, NEMBED, 1, nx,          &
      &                             FFTW_FLAG, FFTW_ESTIMATE)
    endif
#ifdef __TIMER_DO__
  call timer_sta(259)
#endif
!OCL INDEPENDENT[dfftw_execute_dft_c2r,dfftw_execute_dft]
!OCL PARALLEL_STRONG
    do ib = 1, iesize
       if (kimg==1) then
          call dfftw_execute_dft_c2r(plan,wk_afft_l(1,ib),wk_afft_l(1,ib))
       else
          call dfftw_execute_dft    (plan,wk_afft_l(1,ib),wk_afft_l(1,ib))
       endif
    end do
#ifdef __TIMER_DO__
  call timer_end(259)
#endif
    call dfftw_destroy_plan(plan)

1000 continue

!F  if (allocated(req_r))       deallocate(req_r)
!F  if (allocated(req_s))       deallocate(req_s)
!F  if (allocated(sta_r))       deallocate(sta_r)
!F  if (allocated(sta_s))       deallocate(sta_s)
!F  if (allocated(wk_recvcnt))    deallocate(wk_recvcnt)
!F  if (allocated(wk_sendcnt))    deallocate(wk_sendcnt)
!   if (allocated(wk_recvrank))   deallocate(wk_recvrank)
!   if (allocated(wk_sendrank))   deallocate(wk_sendrank)
!F  if (allocated(wk_afft_y))   deallocate(wk_afft_y)
    if (allocated(wk_mp_fft_y)) deallocate(wk_mp_fft_y)

#ifdef __TIMER_SUB__
    call timer_end(109)
#endif
  end subroutine m_FFT_CD_Direct_3D

!------------------------------------------------------------------------------
#else
!else ifndef FFT_USE_SSL2
!------------------------------------------------------------------------------

  subroutine m_FFT_CD_Direct_3D(nfout,wk_afft_l, wk_size, iesize)
    integer, intent(in) :: nfout, wk_size, iesize
    real(kind=DP), dimension(wk_size*kimg,iesize), intent(inout) :: wk_afft_l

    real(kind=DP), allocatable, dimension(:) :: wk_allfft
    integer,      allocatable, dimension(:)   :: wk_recvdsp
    integer, parameter :: FFTW_MEASURE=0
    integer, parameter :: FFTW_ESTIMATE=64
    integer, parameter :: FFTW_RANK=1
    integer            :: FFTW_FLAG
    integer,dimension(FFTW_RANK) :: NFFTW3, NEMBED, NEREAL
    integer :: nx, ny, nz, i, j, k, l, ri, ib, iaddr, lrank
    integer :: irank, itag, icnt_send, icnt_recv
    integer,save, allocatable, dimension(:)   :: req_r, req_s
    integer,save, allocatable, dimension(:,:) :: sta_r, sta_s
    integer,save, allocatable, dimension(:)   :: wk_recvcnt, wk_sendcnt
!   integer,      allocatable, dimension(:)   :: wk_recvrank, wk_sendrank
!   real(kind=DP),allocatable, dimension(:,:) :: wk_recv, wk_send
    real(kind=DP),allocatable, dimension(:,:),save :: wk_recv1, wk_send1
    real(kind=DP),allocatable, dimension(:,:),save :: wk_recv2, wk_send2
    integer,      allocatable, dimension(:)   :: wk_mp_fft_y
    real(kind=DP),allocatable, dimension(:,:),save :: wk_afft_y
    real(kind=DP),allocatable, dimension(:,:,:)   :: wk_gather
    integer :: mpi_comm, myrank, nmrank
    integer(kind=8) :: plan
    integer :: max_x, max_y, max_z
    integer,save :: max_elm = 0, savesize = 0

    integer,save, allocatable, dimension(:,:) :: y2z_recv, y2z_send, z2x_recv, z2x_send
    integer,save :: y2z_rrank, y2z_srank, z2x_rrank, z2x_srank
    integer,save :: y2z_rmax, y2z_smax, y2z_srmax, z2x_rmax, z2x_smax, z2x_srmax

    integer ::  nnx, nny, nnz, iadd0, iadd1, nis, nie, jadd, kadd, nisx, niex, iadd
    integer,dimension(3) :: nsize, isin
#ifdef FFT_ALLTOALL
    integer,dimension(2) :: max_send,max_recv
#endif
#ifdef __TIMER_SUB__
    call mpi_barrier(mpi_ke_world,ierr)
    call timer_sta(109)
#endif

    if (kimg == 1) then
       FFTW_FLAG = +1      ! FFTW_BACKWARD
    else
       FFTW_FLAG = -1      ! FFTW_FORWARD
    endif
    plan = 0

!fj!$$#ifdef CD_FFT_ALL
!fj!$$    mpi_comm = MPI_CommGroup
!fj!$$    myrank = mype
!fj!$$    nmrank = npes
!fj!$$#else
    mpi_comm = mpi_ke_world
    myrank = myrank_g
    nmrank = nrank_g
!fj!$$#endif
    itag = 10

    if (np_fftcd_x==0)then
       nx = 0
       nx = 0
       nz = 0
       go to 1000
    endif

    if (firstcall_cd_direct_3d) then
       savesize = 0
#ifdef __TIMER_DO__
  call timer_sta(157)
#endif
       max_x = maxval(nel_fftcd_x(:))
       max_y = maxval(nel_fftcd_y(:))
       max_z = maxval(nel_fftcd_z(:))
       max_elm = max(max_x,max_y,max_z)

       if(allocated(req_r)) deallocate(req_r)
       if(allocated(req_s)) deallocate(req_s)
       if(allocated(sta_r)) deallocate(sta_r)
       if(allocated(sta_s)) deallocate(sta_s)
       if(allocated(wk_recvcnt)) deallocate(wk_recvcnt)
       if(allocated(wk_sendcnt)) deallocate(wk_sendcnt)
       allocate(req_r(0:nmrank-1), stat=ierr)
       allocate(req_s(0:nmrank-1), stat=ierr)
       allocate(sta_r(MPI_STATUS_SIZE,0:nmrank-1), stat=ierr)
       allocate(sta_s(MPI_STATUS_SIZE,0:nmrank-1), stat=ierr)
       allocate(wk_recvcnt(0:nmrank-1), stat=ierr)
       allocate(wk_sendcnt(0:nmrank-1), stat=ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_FFT_CD_Direct_3D :  Not allocate '
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 51, ierr)
        endif

       if(allocated(y2z_recv)) deallocate(y2z_recv)
       if(allocated(y2z_send)) deallocate(y2z_send)

       allocate(y2z_recv(2,0:nmrank-1))
       allocate(y2z_send(2,0:nmrank-1))
       y2z_send = 0
       y2z_recv = 0
       wk_recvcnt = 0
       wk_sendcnt = 0
       do i = 1, nel_fftcd_z(myrank)
          irank = map_fftcd_y(mp_fftcd_z(i)) - 1
          wk_recvcnt(irank) = wk_recvcnt(irank) + 1
       enddo
       do i = 1, nel_fftcd_y(myrank)
          irank = map_fftcd_z(mp_fftcd_y(i)) - 1
          wk_sendcnt(irank) = wk_sendcnt(irank) + 1
       enddo
       k = 0
       do i = 0, nmrank - 1
          if(wk_recvcnt(i) /= 0) then
             k = k + 1
             y2z_recv(1,i) = wk_recvcnt(i)
             y2z_recv(2,i) = k
          endif
       enddo
       y2z_rrank = k
       k = 0
       do i = 0, nmrank - 1
          if(wk_sendcnt(i) /= 0) then
             k = k + 1
             y2z_send(1,i) = wk_sendcnt(i)
             y2z_send(2,i) = k
          endif
       enddo
       y2z_srank = k
       y2z_rmax = maxval(wk_recvcnt)
       y2z_smax = maxval(wk_sendcnt)
       y2z_srmax = max(y2z_rmax,y2z_smax)

       if(allocated(z2x_recv)) deallocate(z2x_recv)
       if(allocated(z2x_send)) deallocate(z2x_send)
       allocate(z2x_recv(2,0:nmrank-1))
       allocate(z2x_send(2,0:nmrank-1))
       z2x_send = 0
       z2x_recv = 0
       wk_recvcnt = 0
       wk_sendcnt = 0
       do i = 1, nel_fftcd_x(myrank)
          irank = map_fftcd_z(mp_fftcd_x(i)) - 1
          wk_recvcnt(irank) = wk_recvcnt(irank) + 1
       enddo
       do i = 1, nel_fftcd_z(myrank)
          irank = map_fftcd_x(mp_fftcd_z(i)) - 1
          wk_sendcnt(irank) = wk_sendcnt(irank) + 1
       enddo
       k = 0
       do i = 0, nmrank - 1
          if(wk_recvcnt(i) /= 0) then
             k = k + 1
             z2x_recv(1,i) = wk_recvcnt(i)
             z2x_recv(2,i) = k
          endif
       enddo
       z2x_rrank = k
       k = 0
      do i = 0, nmrank - 1
          if(wk_sendcnt(i) /= 0) then
             k = k + 1
             z2x_send(1,i) = wk_sendcnt(i)
             z2x_send(2,i) = k
          endif
       enddo
       z2x_srank = k
       z2x_rmax = maxval(wk_recvcnt)
       z2x_smax = maxval(wk_sendcnt)
       z2x_srmax = max(z2x_rmax,z2x_smax)

#ifdef FFT_ALLTOALL
       max_send(1) = y2z_rmax
       max_send(2) = y2z_smax
       call mpi_allreduce(max_send,max_recv,2,mpi_integer,mpi_max,mpi_fftcd_zy_world,ierr)
       y2z_rmax = max_recv(1)
       y2z_smax = max_recv(2)
       max_send(1) = z2x_rmax
       max_send(2) = z2x_smax
       call mpi_allreduce(max_send,max_recv,2,mpi_integer,mpi_max,mpi_fftcd_xz_world,ierr)
       z2x_rmax = max_recv(1)
       z2x_smax = max_recv(2)
#endif

       firstcall_cd_direct_3d = .false.
#ifdef __TIMER_DO__
  call timer_end(157)
#endif
    endif

#ifdef FFT_ALLTOALL
    if (iesize /= savesize) then
#else
    if (iesize > savesize) then
#endif
#ifdef __TIMER_DO__
  call timer_sta(158)
#endif
       if (allocated(wk_recv1)) deallocate(wk_recv1)
       if (allocated(wk_send1)) deallocate(wk_send1)
       if (allocated(wk_recv2)) deallocate(wk_recv2)
       if (allocated(wk_send2)) deallocate(wk_send2)
       if (allocated(wk_afft_y)) deallocate(wk_afft_y)
       allocate(wk_recv1(y2z_rmax*kimg*iesize,y2z_rrank), stat=ierr)
       allocate(wk_send1(y2z_smax*kimg*iesize,y2z_srank), stat=ierr)
       allocate(wk_recv2(z2x_rmax*kimg*iesize,z2x_rrank), stat=ierr)
       allocate(wk_send2(z2x_smax*kimg*iesize,z2x_srank), stat=ierr)
       allocate(wk_afft_y(max_elm*kimg,iesize) ,stat=ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_FFT_CD_Direct_3D :  Not allocate '
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 52, ierr)
        endif
       savesize = iesize
#ifdef __TIMER_DO__
  call timer_end(158)
#endif
    end if
!
! Y-axis (z-x div)
!
    nx = xyz_fftcd_y(2,1)-xyz_fftcd_y(1,1)+1
    ny = xyz_fftcd_y(2,2)-xyz_fftcd_y(1,2)+1
    nz = xyz_fftcd_y(2,3)-xyz_fftcd_y(1,3)+1

#ifdef __TIMER_DO__
  call timer_sta(159)
#endif
! (x,z,y) -> (y,x,z)
    if (kimg == 1) then
       do ib = 1, iesize
          do k = 1, nz
             do i = 1, nx/2
                do j = 1, ny
                   iadd0 = (i+(k-1)*nx/2+(j-1)*nx/2*nz)*2
                   iadd1 = (j+(i-1)*ny+(k-1)*ny*nx/2)*2
                   wk_afft_y(iadd1-1,ib) = wk_afft_l(iadd0-1,ib)
                   wk_afft_y(iadd1  ,ib) = wk_afft_l(iadd0  ,ib)
                enddo
             enddo
          enddo
       enddo
    else
       do ib = 1, iesize
          do k = 0, nz-1
             do i = 0, nx-1
                do j = 0, ny-1
                   wk_afft_y((j+i*ny+k*ny*nx+1)*2-1,ib) = wk_afft_l((i+k*nx+j*nx*nz+1)*2-1,ib)
                   wk_afft_y((j+i*ny+k*ny*nx+1)*2  ,ib) = wk_afft_l((i+k*nx+j*nx*nz+1)*2  ,ib)
                enddo
             enddo
          enddo
       enddo
    end if
#ifdef __TIMER_DO__
  call timer_end(159)
#endif

#ifdef __TIMER_DO__
  call timer_sta(257)
#endif
    if(kimg==1) then
       do ib = 1, iesize
          call DM_V1DMCFT(wk_afft_y(1,ib),ny,ny,nx*nz/2,-1,ierr)
       end do
!!     call DM_V1DMCFT(wk_afft_y(1,ib),ny,ny,nx*nz/2*iesize,-1,ierr)
    else
       do ib = 1, iesize
          call DM_V1DMCFT(wk_afft_y(1,ib),ny,ny,nx*nz,1,ierr)
       end do
!!     call DM_V1DMCFT(wk_afft_y(1,1),ny,ny,nx*nz*iesize,1,ierr)
    end if
#ifdef __TIMER_DO__
  call timer_end(257)
#endif

#ifdef __TIMER_DO__
  call timer_sta(160)
#endif
#ifdef __TIMER_DO__
  call timer_end(160)
#endif

#ifdef __TIMER_DO__
  call timer_sta(161)
#endif
    if (kimg == 1) then
       do ib = 1, iesize
          do k = 1, nz
             do i = 1, nx/2
                do l = 1, fftcd_Z_y_dim
                   nis = nis_fftcd_Z_y(l)
                   nie = nie_fftcd_Z_y(l)
                   nny = nie - nis + 1
                   iadd0 = nx/2*nny*nz*(ib-1)*2
                   do j = nis, nie
                      iadd  = (j+(i-1)*ny+(k-1)*ny*nx/2)*2
                      iadd1 = iadd0+((j-nis+1)+(i-1)*nny+(k-1)*nny*nx/2)*2
                      wk_send1(iadd1-1,l) = wk_afft_y(iadd-1,ib)
                      wk_send1(iadd1  ,l) = wk_afft_y(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
    else
       do ib = 1, iesize
          do k = 1, nz
             do i = 1, nx
                do l = 1, fftcd_Z_y_dim
                   nis = nis_fftcd_Z_y(l)
                   nie = nie_fftcd_Z_y(l)
                   nny = nie - nis + 1
                   iadd0 = nx*nny*nz*(ib-1)*2
                   do j = nis, nie
                      iadd  = (j+(i-1)*ny+(k-1)*ny*nx)*2
                      iadd1 = iadd0+((j-nis+1)+(i-1)*nny+(k-1)*nny*nx)*2
                      wk_send1(iadd1-1,l) = wk_afft_y(iadd-1,ib)
                      wk_send1(iadd1  ,l) = wk_afft_y(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
    end if
#ifdef __TIMER_DO__
  call timer_end(161)
#endif

#ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_fftcd_zy_world)
  call timer_sta(162)
#endif

    call MPI_ALLTOALL(wk_send1, y2z_smax*kimg*iesize, mpi_double_precision,   &
   &                  wk_recv1, y2z_rmax*kimg*iesize, mpi_double_precision,   &
   &                                                  mpi_fftcd_zy_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Inverse_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 55, ierr)
     endif

#ifdef __TIMER_DO__
  call timer_sta(163)
#endif
#ifdef __TIMER_DO__
  call timer_end(163)
#endif

#else
!else ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_comm)
  call timer_sta(162)
#endif

    icnt_recv = 0
    do lrank = 0, nmrank - 1
       if ((y2z_recv(1,lrank) /= 0)) then
          call mpi_irecv(wk_recv1(1,y2z_recv(2,lrank)), y2z_recv(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_r(icnt_recv), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_CD_Direct_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 55, ierr)
           endif
          icnt_recv = icnt_recv + 1
       endif
    enddo

    icnt_send = 0
    do lrank = 0, nmrank - 1
       if ((y2z_send(1,lrank) /= 0)) then
          call mpi_isend(wk_send1(1,y2z_send(2,lrank)), y2z_send(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_s(icnt_send), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_CD_Direct_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 56, ierr)
           endif
          icnt_send = icnt_send + 1
       endif
    enddo

#ifdef __TIMER_DO__
  call timer_sta(163)
#endif
#ifdef __TIMER_DO__
  call timer_end(163)
#endif

    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_CD_Direct_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 57, ierr)
     endif

    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_CD_Direct_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 58, ierr)
     endif

#endif
!endif ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_end(162)
#endif

#ifdef __TIMER_DO__
  call timer_sta(164)
#endif
! (y,x,z) -> (z,x,y)
    nx = xyz_fftcd_z(2,1)-xyz_fftcd_z(1,1)+1
    ny = xyz_fftcd_z(2,2)-xyz_fftcd_z(1,2)+1
    nz = xyz_fftcd_z(2,3)-xyz_fftcd_z(1,3)+1
    if (kimg == 1) then
       do ib = 1, iesize
          do i = 1, nx/2
             do j = 1, ny
                do l = 1, fftcd_Y_z_dim
                   nis = nis_fftcd_Y_z(l)
                   nie = nie_fftcd_Y_z(l)
                   nnz = nie - nis + 1
                   iadd0 = nx/2*ny*nnz*(ib-1)*2
                   do k = nis, nie
                      iadd1 = iadd0 + (j+(i-1)*ny+(k-nis)*ny*nx/2)*2
                      iadd = (k+(i-1)*nz+(j-1)*nz*nx/2)*2
                      wk_afft_l(iadd-1,ib) = wk_recv1(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv1(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end do
    else
       do ib = 1, iesize
          do i = 1, nx
             do j = 1, ny
                do l = 1, fftcd_Y_z_dim
                   nis = nis_fftcd_Y_z(l)
                   nie = nie_fftcd_Y_z(l)
                   nnz = nie - nis + 1
                   iadd0 = nx*ny*nnz*(ib-1)*2
                   do k = nis, nie
                      iadd1 = iadd0+(j+(i-1)*ny+(k-nis)*ny*nx)*2
                      iadd  = (k+(i-1)*nz+(j-1)*nz*nx)*2
                      wk_afft_l(iadd-1,ib) = wk_recv1(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv1(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end do
    end if
#ifdef __TIMER_DO__
  call timer_end(164)
#endif

!
! Z-axis (x-y div)
!
#ifdef __TIMER_DO__
  call timer_sta(258)
#endif
    if(kimg==1) then
       do ib = 1, iesize
          call DM_V1DMCFT(wk_afft_l(1,ib),nz,nz,nx*ny/2,-1,ierr)
       end do
!!     call DM_V1DMCFT(wk_afft_l(1,ib),nz,nz,nx*ny/2*iesize,-1,ierr)
    else
       do ib = 1, iesize
          call DM_V1DMCFT(wk_afft_l(1,ib),nz,nz,nx*ny,1,ierr)
       end do
!!     call DM_V1DMCFT(wk_afft_l(1,1),nz,nz,nx*ny*iesize,1,ierr)
    end if
#ifdef __TIMER_DO__
  call timer_end(258)
#endif

#ifdef __TIMER_DO__
  call timer_sta(165)
#endif
    if (kimg == 1) then
       do ib = 1, iesize
          do j = 1, ny
             do i = 1, nx/2
                do l = 1, fftcd_X_z_dim
                   nis = nis_fftcd_X_z(l)
                   nie = nie_fftcd_X_z(l)
                   nnz = nie - nis + 1
                   iadd0 = nx/2*ny*nnz*(ib-1)*2
                   do k = nis, nie
                      iadd  = (k+(i-1)*nz+(j-1)*nz*nx/2)*2
                      iadd1 = iadd0 + (k-nis+1+(i-1)*nnz+(j-1)*nnz*nx/2)*2
                      wk_send2(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send2(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
    else
       do ib = 1, iesize
          do j = 1, ny
             do i = 1, nx
                do l = 1, fftcd_X_z_dim
                   nis = nis_fftcd_X_z(l)
                   nie = nie_fftcd_X_z(l)
                   nnz = nie - nis + 1
                   iadd0 = nx*ny*nnz*(ib-1)*2
                   do k = nis, nie
                      iadd  = (k+(i-1)*nz+(j-1)*nz*nx)*2
                      iadd1 = iadd0+(k-nis+1+(i-1)*nnz+(j-1)*nnz*nx)*2
                      wk_send2(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send2(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
    end if
#ifdef __TIMER_DO__
  call timer_end(165)
#endif

#ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_fftcd_xz_world)
  call timer_sta(166)
#endif

    call MPI_ALLTOALL(wk_send2, z2x_smax*kimg*iesize, mpi_double_precision,   &
   &                  wk_recv2, z2x_rmax*kimg*iesize, mpi_double_precision,   &
   &                                                  mpi_fftcd_xz_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_CD_Direct_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 60, ierr)
     endif

#ifdef __TIMER_DO__
  call timer_sta(167)
#endif
#ifdef __TIMER_DO__
  call timer_end(167)
#endif

#else
!else ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_comm)
  call timer_sta(166)
#endif

    icnt_recv = 0
    do lrank = 0, nmrank - 1
       if ((z2x_recv(1,lrank) /= 0)) then
          call mpi_irecv(wk_recv2(1,z2x_recv(2,lrank)), z2x_recv(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_r(icnt_recv), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_CD_Direct_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 60, ierr)
           endif
          icnt_recv = icnt_recv + 1
       endif
    enddo

    icnt_send = 0
    do lrank = 0, nmrank - 1
       if ((z2x_send(1,lrank) /= 0)) then
          call mpi_isend(wk_send2(1,z2x_send(2,lrank)), z2x_send(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_s(icnt_send), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_CD_Direct_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 61, ierr)
           endif
          icnt_send = icnt_send + 1
       endif
    enddo

#ifdef __TIMER_DO__
  call timer_sta(167)
#endif
#ifdef __TIMER_DO__
  call timer_end(167)
#endif

    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_CD_Direct_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 62, ierr)
     endif

    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_CD_Direct_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 63, ierr)
     endif

#endif
!endif ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_end(166)
#endif

#ifdef __TIMER_DO__
  call timer_sta(168)
#endif
    nx = xyz_fftcd_x(2,1)-xyz_fftcd_x(1,1)+1
    ny = xyz_fftcd_x(2,2)-xyz_fftcd_x(1,2)+1
    nz = xyz_fftcd_x(2,3)-xyz_fftcd_x(1,3)+1
!  (z,x,y) -> (x,y,z)
    if (kimg == 1) then
       do l = 1, fftcd_Z_x_dim
          nis = nis_fftcd_Z_x(l)
          nie = nie_fftcd_Z_x(l)
          nnx = nie - nis + 1
          if(mod(nis,2)>0) then
            nisx = nis/2 + 1
          else
            nisx = nis/2
          endif
          niex = nie/2
          do ib = 1, iesize
             iadd0 = nnx/2*ny*nz*(ib-1)*2
             do k = 1, nz
                do j = 1, ny
                   do i = nisx, niex
                      iadd1 = iadd0 + (k+(i-nisx)*nz+(j-1)*nz*nnx/2)*2
!XX                   iadd  = (k+(i-1)*nz+(j-1)*nz*nx/2)*2
                      iadd  = (i+(j-1)*nx/2+(k-1)*nx/2*ny)*2
                      wk_afft_l(iadd-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end do
    else
       do ib = 1, iesize
          do k = 1, nz
             do j = 1, ny
                do l = 1, fftcd_Z_x_dim
                   nis = nis_fftcd_Z_x(l)
                   nie = nie_fftcd_Z_x(l)
                   nnx = nie - nis + 1
                   iadd0 = nnx*ny*nz*(ib-1)*2
                   do i = nis, nie
                      iadd1 = iadd0+(k+(i-nis)*nz+(j-1)*nz*nnx)*2
                      iadd  = (i+(j-1)*nx+(k-1)*nx*ny)*2
                      wk_afft_l(iadd-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end do
    end if
#ifdef __TIMER_DO__
  call timer_end(168)
#endif

!
! X-axis (y-z div)
!
#ifdef __TIMER_DO__
  call timer_sta(259)
#endif
    if(kimg==1) then
! (x,y,z) <- (z,x,y)
!XX    wk_afft_y = wk_afft_l
!XX    do ib = 1, iesize
!XX       do k = 1, nz
!XX          do j = 1, ny
!XX             do i = 1, nx/2
!XX                iadd0 = (k+(i-1)*nz+(j-1)*nz*nx/2)*2
!XX                iadd1 = (i+(j-1)*nx/2+(k-1)*nx/2*ny)*2
!XX                wk_afft_l(iadd1-1,ib) = wk_afft_y(iadd0-1,ib)
!XX                wk_afft_l(iadd1  ,ib) = wk_afft_y(iadd0  ,ib)
!XX             enddo
!XX          enddo
!XX       enddo
!XX    enddo

!.     do ib = 1, iesize
!.        do k = 1, nz
!.           do j = 1, ny
!.              call DM_V1DRCF2(wk_afft_l(1+nx*(j-1)+nx*ny*(k-1),ib),nx-2,      &
!.             &                wk_afft_l(1+nx*(j-1)+nx*ny*(k-1),ib),1,-1,ierr)
!.           enddo
!.        enddo
!.     enddo
       nsize(1:3) = (/nx-2,ny,nz/)
       isin(1:3)  = (/-1,0,0/)
!OCL INDEPENDENT[ DVMRF2 ]
!OCL PARALLEL_STRONG
       do ib = 1, iesize
          call DVMRF2(wk_afft_l(1,ib),nsize,3,isin,-1,ierr)
       end do
    else
       do ib = 1, iesize
          call DM_V1DMCFT(wk_afft_l(1,ib),nx,nx,ny*nz,1,ierr)
       end do
!      call DM_V1DMCFT(wk_afft_l(1,1),nx,nx,ny*nz*iesize,1,ierr)
    end if
#ifdef __TIMER_DO__
  call timer_end(259)
#endif

1000 continue

!F  if (allocated(req_r))       deallocate(req_r)
!F  if (allocated(req_s))       deallocate(req_s)
!F  if (allocated(sta_r))       deallocate(sta_r)
!F  if (allocated(sta_s))       deallocate(sta_s)
!F  if (allocated(wk_recvcnt))    deallocate(wk_recvcnt)
!F  if (allocated(wk_sendcnt))    deallocate(wk_sendcnt)
!   if (allocated(wk_recvrank))   deallocate(wk_recvrank)
!   if (allocated(wk_sendrank))   deallocate(wk_sendrank)
!F  if (allocated(wk_afft_y))   deallocate(wk_afft_y)
    if (allocated(wk_mp_fft_y)) deallocate(wk_mp_fft_y)

#ifdef __TIMER_SUB__
    call timer_end(109)
#endif
  end subroutine m_FFT_CD_Direct_3D
!------------------------------------------------------------------------------
#endif
!endif ifndef FFT_USE_SSL2

#ifndef FFT_USE_SSL2
!------------------------------------------------------------------------------
  subroutine m_FFT_CD_Inverse_3D(nfout,wk_afft_l, wk_size, iesize)
    integer, intent(in) :: nfout, wk_size, iesize
    real(kind=DP), dimension(wk_size*kimg,iesize), intent(inout) :: wk_afft_l

    real(kind=DP), allocatable, dimension(:) :: wk_allfft
    integer,      allocatable, dimension(:)   :: wk_recvdsp
    integer, parameter :: FFTW_MEASURE=0
    integer, parameter :: FFTW_ESTIMATE=64
    integer, parameter :: FFTW_RANK=1
    integer            :: FFTW_FLAG
    integer,dimension(1) :: NFFTW3, NEMBED, NEREAL
    integer :: nx, ny, nz, i, j, k, l, ri, ib, iadd, lrank, lx
    integer :: irank, itag, icnt_send, icnt_recv
    integer,save, allocatable, dimension(:)   :: req_r, req_s
    integer,save, allocatable, dimension(:,:) :: sta_r, sta_s
    integer,save, allocatable, dimension(:)   :: wk_recvcnt, wk_sendcnt
!F  real(kind=DP),allocatable, dimension(:,:) :: wk_recv, wk_send
    real(kind=DP),allocatable, dimension(:,:),save :: wk_recv1, wk_send1
    real(kind=DP),allocatable, dimension(:,:),save :: wk_recv2, wk_send2
    integer,      allocatable, dimension(:)   :: wk_mp_fft_y
    real(kind=DP),allocatable, dimension(:,:),save  :: wk_afft_y
    real(kind=DP),allocatable, dimension(:,:,:) :: wk_gather
    integer :: mpi_comm, myrank, nmrank
    integer(kind=8) :: plan
    integer :: max_x, max_y, max_z
    integer,save :: max_elm = 0, savesize = 0

    integer,save, allocatable, dimension(:,:) :: x2z_recv, x2z_send, z2y_recv, z2y_send
    integer,save :: x2z_rrank, x2z_srank, z2y_rrank, z2y_srank
    integer,save :: x2z_rmax, x2z_smax, x2z_srmax, z2y_rmax, z2y_smax, z2y_srmax

    integer ::  nnx, nny, nnz, iadd0, iadd1, nis, nie, jadd, kadd
#ifdef FFT_ALLTOALL
    integer,dimension(2) :: max_send,max_recv
#endif
#ifdef __TIMER_SUB__
    call mpi_barrier(mpi_ke_world,ierr)
    call timer_sta(110)
#endif

    if (kimg == 1) then
       FFTW_FLAG = -1      ! FFTW_FORWARD
    else
       FFTW_FLAG = +1      ! FFTW_BACKWARD
    endif
    plan = 0

!fj!$$#ifdef CD_FFT_ALL
!fj!$$    mpi_comm = MPI_CommGroup
!fj!$$    myrank = mype
!fj!$$    nmrank = npes
!fj!$$#else
    mpi_comm = mpi_ke_world
    myrank = myrank_g
    nmrank = nrank_g
!fj!$$#endif
    itag = 10

    if (np_fftcd_x==0)then
       nx = 0
       nx = 0
       nz = 0
       go to 1000
    endif

    if (firstcall_cd_inverse_3d) then
       savesize = 0
#ifdef __TIMER_DO__
  call timer_sta(173)
#endif
       max_x = maxval(nel_fftcd_x(:))
       max_y = maxval(nel_fftcd_y(:))
       max_z = maxval(nel_fftcd_z(:))
       max_elm = max(max_x,max_y,max_z)

       if(allocated(req_r)) deallocate(req_r)
       if(allocated(req_s)) deallocate(req_s)
       if(allocated(sta_r)) deallocate(sta_r)
       if(allocated(sta_s)) deallocate(sta_s)
       if(allocated(wk_recvcnt)) deallocate(wk_recvcnt)
       if(allocated(wk_sendcnt)) deallocate(wk_sendcnt)
       allocate(req_r(0:nmrank-1), stat=ierr)
       allocate(req_s(0:nmrank-1), stat=ierr)
       allocate(sta_r(MPI_STATUS_SIZE,0:nmrank-1), stat=ierr)
       allocate(sta_s(MPI_STATUS_SIZE,0:nmrank-1), stat=ierr)
       allocate(wk_recvcnt(0:nmrank-1), stat=ierr)
       allocate(wk_sendcnt(0:nmrank-1), stat=ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_FFT_CD_Inverse_3D :  Not allocate '
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 68, ierr)
        endif

       if(allocated(x2z_recv)) deallocate(x2z_recv)
       if(allocated(x2z_send)) deallocate(x2z_send)
       allocate(x2z_recv(2,0:nmrank-1))
       allocate(x2z_send(2,0:nmrank-1))
       x2z_send = 0
       x2z_recv = 0
       wk_recvcnt = 0
       wk_sendcnt = 0
       do i = 1, nel_fftcd_z(myrank)
          irank = map_fftcd_x(mp_fftcd_z(i)) - 1
          wk_recvcnt(irank) = wk_recvcnt(irank) + 1
       enddo
       do i = 1, nel_fftcd_x(myrank)
          irank = map_fftcd_z(mp_fftcd_x(i)) - 1
          wk_sendcnt(irank) = wk_sendcnt(irank) + 1
       enddo
       k = 0
       do i = 0, nmrank - 1
          if(wk_recvcnt(i) /= 0) then
             k = k + 1
             x2z_recv(1,i) = wk_recvcnt(i)
             x2z_recv(2,i) = k
          endif
       enddo
       x2z_rrank = k
       k = 0
       do i = 0, nmrank - 1
          if(wk_sendcnt(i) /= 0) then
             k = k + 1
             x2z_send(1,i) = wk_sendcnt(i)
             x2z_send(2,i) = k
          endif
       enddo
       x2z_srank = k
       x2z_rmax = maxval(wk_recvcnt)
       x2z_smax = maxval(wk_sendcnt)
       x2z_srmax = max(x2z_rmax,x2z_smax)

       if(allocated(z2y_recv)) deallocate(z2y_recv)
       if(allocated(z2y_send)) deallocate(z2y_send)
       allocate(z2y_recv(2,0:nmrank-1))
       allocate(z2y_send(2,0:nmrank-1))
       z2y_send = 0
       z2y_recv = 0
       wk_recvcnt = 0
       wk_sendcnt = 0

       nx = xyz_fftcd_y(2,1)-xyz_fftcd_y(1,1)+1
       ny = xyz_fftcd_y(2,2)-xyz_fftcd_y(1,2)+1
       nz = xyz_fftcd_y(2,3)-xyz_fftcd_y(1,3)+1
       allocate(wk_mp_fft_y(nx*ny*nz) ,stat=ierr)
       do k = 0, nz-1
          do j = 0, ny-1
             do i = 0, nx-1
                wk_mp_fft_y(i+j*nx+k*nx*ny+1) = mp_fftcd_y(i+k*nx+j*nx*nz+1)
             enddo
          enddo
       enddo

       do i = 1, nel_fftcd_y(myrank)
          irank = map_fftcd_z(wk_mp_fft_y(i)) - 1
          wk_recvcnt(irank) = wk_recvcnt(irank) + 1
       enddo
       do i = 1, nel_fftcd_z(myrank)
          irank = map_fftcd_y(mp_fftcd_z(i)) - 1
          wk_sendcnt(irank) = wk_sendcnt(irank) + 1
       enddo
       k = 0
       do i = 0, nmrank - 1
          if(wk_recvcnt(i) /= 0) then
             k = k + 1
             z2y_recv(1,i) = wk_recvcnt(i)
             z2y_recv(2,i) = k
          endif
       enddo
       z2y_rrank = k
       k = 0
      do i = 0, nmrank - 1
          if(wk_sendcnt(i) /= 0) then
             k = k + 1
             z2y_send(1,i) = wk_sendcnt(i)
             z2y_send(2,i) = k
          endif
       enddo
       z2y_srank = k
       z2y_rmax = maxval(wk_recvcnt)
       z2y_smax = maxval(wk_sendcnt)
       z2y_srmax = max(z2y_rmax,z2y_smax)

       deallocate(wk_mp_fft_y)

#ifdef FFT_ALLTOALL
       max_send(1) = x2z_rmax
       max_send(2) = x2z_smax
       call mpi_allreduce(max_send,max_recv,2,mpi_integer,mpi_max,mpi_fftcd_xz_world,ierr)
       x2z_rmax = max_recv(1)
       x2z_smax = max_recv(2)
       max_send(1) = z2y_rmax
       max_send(2) = z2y_smax
       call mpi_allreduce(max_send,max_recv,2,mpi_integer,mpi_max,mpi_fftcd_zy_world,ierr)
       z2y_rmax = max_recv(1)
       z2y_smax = max_recv(2)
#endif

       firstcall_cd_inverse_3d = .false.
#ifdef __TIMER_DO__
  call timer_end(173)
#endif
    endif

#ifdef FFT_ALLTOALL
    if (iesize /= savesize) then
#else
    if (iesize > savesize) then
#endif
#ifdef __TIMER_DO__
  call timer_sta(174)
#endif
       if (allocated(wk_recv1)) deallocate(wk_recv1)
       if (allocated(wk_send1)) deallocate(wk_send1)
       if (allocated(wk_recv2)) deallocate(wk_recv2)
       if (allocated(wk_send2)) deallocate(wk_send2)
       if (allocated(wk_afft_y)) deallocate(wk_afft_y)
       allocate(wk_recv1(x2z_rmax*kimg*iesize,x2z_rrank), stat=ierr)
       allocate(wk_send1(x2z_smax*kimg*iesize,x2z_srank), stat=ierr)
       allocate(wk_recv2(z2y_rmax*kimg*iesize,z2y_rrank), stat=ierr)
       allocate(wk_send2(z2y_smax*kimg*iesize,z2y_srank), stat=ierr)
       allocate(wk_afft_y(max_elm*kimg,iesize) ,stat=ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_FFT_CD_Inverse_3D :  Not allocate '
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 69, ierr)
        endif
       savesize = iesize
#ifdef __TIMER_DO__
  call timer_end(174)
#endif
    end if

!
! X-axis (y-z div)
!
    nx = xyz_fftcd_x(2,1)-xyz_fftcd_x(1,1)+1
    ny = xyz_fftcd_x(2,2)-xyz_fftcd_x(1,2)+1
    nz = xyz_fftcd_x(2,3)-xyz_fftcd_x(1,3)+1
    if(kimg==1) then
       NFFTW3(1) = fft_box_size_CD_3D(1,1)
       NEREAL(1) = fft_box_size_CD_3D(1,1)
       NEMBED(1) = fft_box_size_CD_3D(1,1)/2
       call dfftw_plan_many_dft_r2c(plan, FFTW_RANK, NFFTW3, ny*nz,   &
      &                             wk_afft_l, NEREAL, 1, nx,          &
      &                             wk_afft_l, NEMBED, 1, nx/2,        &
      &                             FFTW_ESTIMATE)
    else
       NFFTW3(1) = fft_box_size_CD_3D(1,1)
       NEMBED(1) = nx
       NEMBED(1) = fft_box_size_CD_3D(1,0)
       call dfftw_plan_many_dft    (plan, FFTW_RANK, NFFTW3, ny*nz,   &
      &                             wk_afft_l, NEMBED, 1, nx,          &
      &                             wk_afft_l, NEMBED, 1, nx,          &
      &                             FFTW_FLAG, FFTW_ESTIMATE)
    endif
#ifdef __TIMER_DO__
  call timer_sta(260)
#endif
!OCL INDEPENDENT[dfftw_execute_dft_r2c,dfftw_execute_dft]
!OCL PARALLEL_STRONG
    do ib = 1, iesize
       if (kimg==1) then
          call dfftw_execute_dft_r2c(plan,wk_afft_l(1,ib),wk_afft_l(1,ib))
       else
          call dfftw_execute_dft    (plan,wk_afft_l(1,ib),wk_afft_l(1,ib))
       endif
    enddo
#ifdef __TIMER_DO__
  call timer_end(260)
#endif
    call dfftw_destroy_plan(plan)

#ifdef __TIMER_DO__
  call timer_sta(175)
#endif
    do l = 1, fftcd_Z_x_dim
       nis = nis_fftcd_Z_x(l)
       nie = nie_fftcd_Z_x(l)
       nnx = nie - nis + 1
       if (kimg == 1) then
          do ib = 1, iesize
             iadd0 = nnx*ny*nz*(ib-1)
             do k = 1, nz
                do j = 1, ny
                   do i = nis, nie
                      wk_send1(iadd0+(i-nis+1+(j-1)*nnx+(k-1)*nnx*ny),l) = wk_afft_l(i+(j-1)*nx+(k-1)*nx*ny,ib)
                   end do
                end do
             end do
          end do
       else
          do ib = 1, iesize
             iadd0 = nnx*ny*nz*(ib-1)*2
             do k = 1, nz
                jadd = (k-1)*nnx*ny
                do j = 1, ny
                   kadd = (j-1)*nnx+jadd
                   do i = nis, nie
                      iadd1 = iadd0+(i-nis+1+kadd)*2
                      iadd  = (i+(j-1)*nx+(k-1)*nx*ny)*2
                      wk_send1(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send1(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end if
    end do
#ifdef __TIMER_DO__
  call timer_end(175)
#endif

#ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_fftcd_xz_world)
  call timer_sta(176)
#endif

    call MPI_ALLTOALL(wk_send1, x2z_smax*kimg*iesize, mpi_double_precision,   &
   &                  wk_recv1, x2z_rmax*kimg*iesize, mpi_double_precision,   &
   &                                                  mpi_fftcd_xz_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_CD_Inverse_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 70, ierr)
     endif

#ifdef __TIMER_DO__
  call timer_sta(177)
#endif
#ifdef __TIMER_DO__
  call timer_end(177)
#endif

#else
!else ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_comm)
  call timer_sta(176)
#endif

    icnt_recv = 0
    do lrank = 0, nmrank - 1
       if ((lrank /= myrank) .and. (x2z_recv(1,lrank) /= 0)) then
          call mpi_irecv(wk_recv1(1,x2z_recv(2,lrank)), x2z_recv(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_r(icnt_recv), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_CD_Inverse_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 70, ierr)
           endif
          icnt_recv = icnt_recv + 1
       endif
    enddo

    icnt_send = 0
    do lrank = 0, nmrank - 1
       if ((lrank /= myrank) .and. (x2z_send(1,lrank) /= 0)) then
          call mpi_isend(wk_send1(1,x2z_send(2,lrank)), x2z_send(1,lrank)*kimg*iesize, &
                         mpi_double_precision, lrank, itag, mpi_comm, req_s(icnt_send), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_CD_Inverse_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 71, ierr)
           endif
          icnt_send = icnt_send + 1
       endif
    enddo

#ifdef __TIMER_DO__
  call timer_sta(177)
#endif
    do i = 1, x2z_recv(1,myrank)*kimg*iesize
       wk_recv1(i,x2z_recv(2,myrank)) = wk_send1(i,x2z_send(2,myrank))
    enddo
#ifdef __TIMER_DO__
  call timer_end(177)
#endif

    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_CD_Inverse_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 72, ierr)
     endif
    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_CD_Inverse_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 73, ierr)
     endif
#endif
!endif ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_end(176)
#endif

#ifdef __TIMER_DO__
  call timer_sta(178)
#endif
    nx = xyz_fftcd_z(2,1)-xyz_fftcd_z(1,1)+1
    ny = xyz_fftcd_z(2,2)-xyz_fftcd_z(1,2)+1
    nz = xyz_fftcd_z(2,3)-xyz_fftcd_z(1,3)+1
    do l = 1, fftcd_X_z_dim
       nis = nis_fftcd_X_z(l)
       nie = nie_fftcd_X_z(l)
       nnz = nie - nis + 1
       if (kimg == 1) then
          do ib = 1, iesize
             iadd0 = nx*ny*nnz*(ib-1)
             do k = nis, nie
                do j = 1, ny
                   do i = 1, nx
                      wk_afft_l(i+(j-1)*nx+(k-1)*nx*ny,ib) = wk_recv1(iadd0+i+(j-1)*nx+(k-nis)*nx*ny,l)
                   end do
                end do
             end do
          end do
       else
          do ib = 1, iesize
             iadd0 = nx*ny*nnz*(ib-1)*2
             do k = nis, nie
                do j = 1, ny
                   do i = 1, nx
                      iadd1 = iadd0+(i+(j-1)*nx+(k-nis)*nx*ny)*2
                      iadd  = (i+(j-1)*nx+(k-1)*nx*ny)*2
                      wk_afft_l(iadd-1,ib) = wk_recv1(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv1(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end if
    end do
#ifdef __TIMER_DO__
  call timer_end(178)
#endif

!
! Z-axis (x-y div)
!
    nx = xyz_fftcd_z(2,1)-xyz_fftcd_z(1,1)+1
    ny = xyz_fftcd_z(2,2)-xyz_fftcd_z(1,2)+1
    nz = xyz_fftcd_z(2,3)-xyz_fftcd_z(1,3)+1

    if(kimg==1) then
       NFFTW3(1) = fft_box_size_CD_3D(3,1)
       NEMBED(1) = fft_box_size_CD_3D(3,0)
       call dfftw_plan_many_dft    (plan, FFTW_RANK, NFFTW3, nx*ny/2, &
      &                             wk_afft_l, NEMBED, nx*ny/2, 1,     &
      &                             wk_afft_l, NEMBED, nx*ny/2, 1,     &
      &                             FFTW_FLAG, FFTW_ESTIMATE)
    else
       NFFTW3(1) = fft_box_size_CD_3D(3,1)
       NEMBED(1) = nz
       NEMBED(1) = fft_box_size_CD_3D(3,0)
       call dfftw_plan_many_dft    (plan, FFTW_RANK, NFFTW3, nx*ny,   &
      &                             wk_afft_l, NEMBED, nx*ny, 1,       &
      &                             wk_afft_l, NEMBED, nx*ny, 1,       &
      &                             FFTW_FLAG, FFTW_ESTIMATE)
    endif
#ifdef __TIMER_DO__
  call timer_sta(261)
#endif
!OCL INDEPENDENT[dfftw_execute_dft]
!OCL PARALLEL_STRONG
    do ib = 1, iesize
       call dfftw_execute_dft(plan,wk_afft_l(1,ib),wk_afft_l(1,ib))
    end do
#ifdef __TIMER_DO__
  call timer_end(261)
#endif
    call dfftw_destroy_plan(plan)

    nx = xyz_fftcd_y(2,1)-xyz_fftcd_y(1,1)+1
    ny = xyz_fftcd_y(2,2)-xyz_fftcd_y(1,2)+1
    nz = xyz_fftcd_y(2,3)-xyz_fftcd_y(1,3)+1
!XX allocate(wk_mp_fft_y(nx*ny*nz) ,stat=ierr)
!XX  if (ierr /= 0) then
!XX     write(nfout,*)' m_FFT_CD_Inverse_3D :  Not allocate '
!XX     call flush(nfout)
!XX     call mpi_abort(mpi_comm_world, 74, ierr)
!XX  endif
#ifdef __TIMER_DO__
  call timer_sta(179)
#endif
!XX do k = 0, nz-1
!XX    do j = 0, ny-1
!XX       do i = 0, nx-1
!XX          wk_mp_fft_y(i+j*nx+k*nx*ny+1) = mp_fftcd_y(i+k*nx+j*nx*nz+1)
!XX       enddo
!XX    enddo
!XX enddo
#ifdef __TIMER_DO__
  call timer_end(179)
#endif

!F  allocate(wk_recv2(z2y_rmax*kimg*iesize,z2y_rrank), stat=ierr)
!F  allocate(wk_send2(z2y_smax*kimg*iesize,z2y_srank), stat=ierr)
!F   if (ierr /= 0) then
!F      write(nfout,*)' m_FFT_CD_Inverse_3D :  Not allocate '
!F      call flush(nfout)
!F      call mpi_abort(mpi_comm_world, 75, ierr)
!F   endif

#ifdef __TIMER_DO__
  call timer_sta(180)
#endif
    nx = xyz_fftcd_z(2,1)-xyz_fftcd_z(1,1)+1
    ny = xyz_fftcd_z(2,2)-xyz_fftcd_z(1,2)+1
    nz = xyz_fftcd_z(2,3)-xyz_fftcd_z(1,3)+1
    do l = 1, fftcd_Y_z_dim
       nis = nis_fftcd_Y_z(l)
       nie = nie_fftcd_Y_z(l)
       nnz = nie - nis + 1
       if (kimg == 1) then
          do ib = 1, iesize
             iadd0 = nx*ny*nnz*(ib-1)
             do k = nis, nie
                do j = 1, ny
                   do i = 1, nx
                      wk_send2(iadd0+(i+(j-1)*nx+(k-nis)*nx*ny),l) = wk_afft_l(i+(j-1)*nx+(k-1)*nx*ny,ib)
                   end do
                end do
             end do
          end do
       else
          do ib = 1, iesize
             iadd0 = nx*ny*nnz*(ib-1)*2
             do k = nis, nie
                jadd = (k-nis)*nx*ny
                do j = 1, ny
                   kadd = (j-1)*nx+jadd
                   do i = 1, nx
                      iadd1 = iadd0+(i+kadd)*2
                      iadd  = (i+(j-1)*nx+(k-1)*nx*ny)*2
                      wk_send2(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send2(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end if
    end do
#ifdef __TIMER_DO__
  call timer_end(180)
#endif

#ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_fftcd_zy_world)
  call timer_sta(181)
#endif

    call MPI_ALLTOALL(wk_send2, z2y_smax*kimg*iesize, mpi_double_precision,   &
   &                  wk_recv2, z2y_rmax*kimg*iesize, mpi_double_precision,   &
   &                                                  mpi_fftcd_zy_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_CD_Inverse_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 76, ierr)
     endif

#ifdef __TIMER_DO__
  call timer_sta(182)
#endif
#ifdef __TIMER_DO__
  call timer_end(182)
#endif

#else
!else ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_comm)
  call timer_sta(181)
#endif

    icnt_recv = 0
    do lrank = 0, nmrank - 1
       if ((lrank /= myrank) .and. (z2y_recv(1,lrank) /= 0)) then
          call mpi_irecv(wk_recv2(1,z2y_recv(2,lrank)), z2y_recv(1,lrank)*kimg*iesize, &
                         mpi_double_precision, lrank, itag, mpi_comm, req_r(icnt_recv), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_CD_FFT_Inverse_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 76, ierr)
           endif
          icnt_recv = icnt_recv + 1
       endif
    enddo

    icnt_send = 0
    do lrank = 0, nmrank - 1
       if ((lrank /= myrank) .and. (z2y_send(1,lrank) /= 0)) then
          call mpi_isend(wk_send2(1,z2y_send(2,lrank)), z2y_send(1,lrank)*kimg*iesize, &
                         mpi_double_precision, lrank, itag, mpi_comm, req_s(icnt_send), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_CD_FFT_Inverse_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 77, ierr)
           endif
          icnt_send = icnt_send + 1
       endif
    enddo

#ifdef __TIMER_DO__
  call timer_sta(182)
#endif
    do i = 1, z2y_recv(1,myrank)*kimg*iesize
       wk_recv2(i,z2y_recv(2,myrank)) = wk_send2(i,z2y_send(2,myrank))
    enddo
#ifdef __TIMER_DO__
  call timer_end(182)
#endif

    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_CD_Inverse_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 78, ierr)
     endif
    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_CD_Inverse_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 79, ierr)
     endif

#endif
!endif ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_end(181)
#endif

#ifdef __TIMER_DO__
  call timer_sta(183)
#endif
    nx = xyz_fftcd_y(2,1)-xyz_fftcd_y(1,1)+1
    ny = xyz_fftcd_y(2,2)-xyz_fftcd_y(1,2)+1
    nz = xyz_fftcd_y(2,3)-xyz_fftcd_y(1,3)+1
    do l = 1, fftcd_Z_y_dim
       nis = nis_fftcd_Z_y(l)
       nie = nie_fftcd_Z_y(l)
       nny = nie - nis + 1
       if (kimg == 1) then
          do ib = 1, iesize
             iadd0 = nx*nny*nz*(ib-1)
             do k = 1, nz
                do j = nis, nie
                   do i = 1, nx
                      wk_afft_y(i+(j-1)*nx+(k-1)*nx*ny,ib) = wk_recv2(iadd0+i+(j-nis)*nx+(k-1)*nx*nny,l)
                   end do
                end do
             end do
          end do
       else
          do ib = 1, iesize
             iadd0 = nx*nny*nz*(ib-1)*2
             do k = 1, nz
                do j = nis, nie
                   do i = 1, nx
                      iadd1 = iadd0+(i+(j-nis)*nx+(k-1)*nx*nny)*2
                      iadd  = (i+(j-1)*nx+(k-1)*nx*ny)*2
                      wk_afft_y(iadd-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_y(iadd  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end if
    end do
#ifdef __TIMER_DO__
  call timer_end(183)
#endif

#ifdef __TIMER_DO__
  call timer_sta(184)
#endif
    if (kimg == 1) then
       do ib = 1, iesize
          do k = 0, nz-1
             do j = 0, ny-1
                do i = 0, nx-1
                   do ri = 1, kimg
                      wk_afft_l((i+k*nx+j*nx*nz+1)*kimg+(ri-kimg),ib) = wk_afft_y((i+j*nx+k*nx*ny+1)*kimg+(ri-kimg),ib)
                   enddo
                enddo
             enddo
          enddo
       enddo
    else
       do ib = 1, iesize
          do k = 0, nz-1
             do j = 0, ny-1
                do i = 0, nx-1
                   wk_afft_l((i+k*nx+j*nx*nz+1)*2-1,ib) = wk_afft_y((i+j*nx+k*nx*ny+1)*2-1,ib)
                   wk_afft_l((i+k*nx+j*nx*nz+1)*2  ,ib) = wk_afft_y((i+j*nx+k*nx*ny+1)*2  ,ib)
                enddo
             enddo
          enddo
       enddo
    endif
#ifdef __TIMER_DO__
  call timer_end(184)
#endif

!
! Y-axis (z-x div)
!
    if(kimg==1) then
       NFFTW3(1) = fft_box_size_CD_3D(2,1)
       NEMBED(1) = fft_box_size_CD_3D(2,1)
       call dfftw_plan_many_dft    (plan, FFTW_RANK, NFFTW3, nz*nx/2, &
      &                             wk_afft_l, NEMBED, nz*nx/2, 1,    &
      &                             wk_afft_l, NEMBED, nz*nx/2, 1,    &
      &                             FFTW_FLAG, FFTW_ESTIMATE)
    else
       NFFTW3(1) = fft_box_size_CD_3D(2,1)
       NEMBED(1) = ny
       NEMBED(1) = fft_box_size_CD_3D(2,1)
       lx = fft_box_size_CD_3D(1,1)
       call dfftw_plan_many_dft    (plan, FFTW_RANK, NFFTW3, nz*nx,   &
      &                             wk_afft_l, NEMBED, nz*nx, 1,      &
      &                             wk_afft_l, NEMBED, nz*nx, 1,      &
      &                             FFTW_FLAG, FFTW_ESTIMATE)
    end if
#ifdef __TIMER_DO__
  call timer_sta(262)
#endif
!OCL INDEPENDENT[dfftw_execute_dft]
!OCL PARALLEL_STRONG
    do ib = 1,iesize
       call dfftw_execute_dft(plan,wk_afft_l(1,ib),wk_afft_l(1,ib))
    end do
#ifdef __TIMER_DO__
  call timer_end(262)
#endif
    call dfftw_destroy_plan(plan)

!XX!   do ib = 1, iesize
!XX!      do k = 0, nz-1
!XX!         do j = 0, ny-1
!XX!            do i = 0, nx-1
!XX!   do ri = 1, kimg
!XX!               wk_afft_y((i+j*nx+k*nx*ny+1)*kimg+(ri-kimg),ib) = wk_afft_l((i+k*nx+j*nx*nz+1)*kimg+(ri-kimg),ib)
!XX!   enddo
!XX!            enddo
!XX!         enddo
!XX!      enddo
!XX!   enddo
!XX!   do ib = 1, iesize
!XX!      do k = 0, nz-1
!XX!         do j = 0, ny-1
!XX!            do i = 0, nx-1
!XX!   do ri = 1, kimg
!XX!               wk_afft_l((i+j*nx+k*nx*ny+1)*kimg+(ri-kimg),ib) = wk_afft_y((i+j*nx+k*nx*ny+1)*kimg+(ri-kimg),ib)
!XX!   enddo
!XX!            enddo
!XX!         enddo
!XX!      enddo
!XX!   enddo

1000 continue

!F  if (allocated(req_r))       deallocate(req_r)
!F  if (allocated(req_s))       deallocate(req_s)
!F  if (allocated(sta_r))       deallocate(sta_r)
!F  if (allocated(sta_s))       deallocate(sta_s)
!F  if (allocated(wk_recvcnt))    deallocate(wk_recvcnt)
!F  if (allocated(wk_sendcnt))    deallocate(wk_sendcnt)
!F  if (allocated(wk_afft_y))   deallocate(wk_afft_y)
    if (allocated(wk_mp_fft_y)) deallocate(wk_mp_fft_y)

#ifdef __TIMER_SUB__
    call timer_end(110)
#endif
  end subroutine m_FFT_CD_Inverse_3D
!------------------------------------------------------------------------------
#else
!else ifndef FFT_USE_SSL2
!------------------------------------------------------------------------------
! For Fujitsu SSL2
  subroutine m_FFT_CD_Inverse_3D(nfout,wk_afft_l, wk_size, iesize)
    integer, intent(in) :: nfout, wk_size, iesize
    real(kind=DP), dimension(wk_size*kimg,iesize), intent(inout) :: wk_afft_l

    real(kind=DP), allocatable, dimension(:) :: wk_allfft
    integer,      allocatable, dimension(:)   :: wk_recvdsp
    integer, parameter :: FFTW_MEASURE=0
    integer, parameter :: FFTW_ESTIMATE=64
    integer, parameter :: FFTW_RANK=1
    integer            :: FFTW_FLAG
    integer,dimension(1) :: NFFTW3, NEMBED, NEREAL
    integer :: nx, ny, nz, i, j, k, l, ri, ib, iadd, lrank, lx
    integer :: irank, itag, icnt_send, icnt_recv
    integer,save, allocatable, dimension(:)   :: req_r, req_s
    integer,save, allocatable, dimension(:,:) :: sta_r, sta_s
    integer,save, allocatable, dimension(:)   :: wk_recvcnt, wk_sendcnt
!F  real(kind=DP),allocatable, dimension(:,:) :: wk_recv, wk_send
    real(kind=DP),allocatable, dimension(:,:),save :: wk_recv1, wk_send1
    real(kind=DP),allocatable, dimension(:,:),save :: wk_recv2, wk_send2
    integer,      allocatable, dimension(:)   :: wk_mp_fft_y
    real(kind=DP),allocatable, dimension(:,:),save  :: wk_afft_y
    real(kind=DP),allocatable, dimension(:,:,:) :: wk_gather
    integer :: mpi_comm, myrank, nmrank
    integer(kind=8) :: plan
    integer :: max_x, max_y, max_z
    integer,save :: max_elm = 0, savesize = 0

    integer,save, allocatable, dimension(:,:) :: x2z_recv, x2z_send, z2y_recv, z2y_send
    integer,save :: x2z_rrank, x2z_srank, z2y_rrank, z2y_srank
    integer,save :: x2z_rmax, x2z_smax, x2z_srmax, z2y_rmax, z2y_smax, z2y_srmax

    integer ::  nnx, nny, nnz, iadd0, iadd1, nis, nie, jadd, kadd
    integer,dimension(3) :: nsize, isin
#ifdef FFT_ALLTOALL
    integer,dimension(2) :: max_send,max_recv
#endif
#ifdef __TIMER_SUB__
    call mpi_barrier(mpi_ke_world,ierr)
    call timer_sta(110)
#endif

    if (kimg == 1) then
       FFTW_FLAG = -1      ! FFTW_FORWARD
    else
       FFTW_FLAG = +1      ! FFTW_BACKWARD
    endif
    plan = 0

!fj!$$#ifdef CD_FFT_ALL
!fj!$$    mpi_comm = MPI_CommGroup
!fj!$$    myrank = mype
!fj!$$    nmrank = npes
!fj!$$#else
    mpi_comm = mpi_ke_world
    myrank = myrank_g
    nmrank = nrank_g
!fj!$$#endif
    itag = 10

    if (np_fftcd_x==0)then
       nx = 0
       nx = 0
       nz = 0
       go to 1000
    endif

    if (firstcall_cd_inverse_3d) then
       savesize = 0
#ifdef __TIMER_DO__
  call timer_sta(173)
#endif
       max_x = maxval(nel_fftcd_x(:))
       max_y = maxval(nel_fftcd_y(:))
       max_z = maxval(nel_fftcd_z(:))
       max_elm = max(max_x,max_y,max_z)

       if(allocated(req_r)) deallocate(req_r)
       if(allocated(req_s)) deallocate(req_s)
       if(allocated(sta_r)) deallocate(sta_r)
       if(allocated(sta_s)) deallocate(sta_s)
       if(allocated(wk_recvcnt)) deallocate(wk_recvcnt)
       if(allocated(wk_sendcnt)) deallocate(wk_sendcnt)
       allocate(req_r(0:nmrank-1), stat=ierr)
       allocate(req_s(0:nmrank-1), stat=ierr)
       allocate(sta_r(MPI_STATUS_SIZE,0:nmrank-1), stat=ierr)
       allocate(sta_s(MPI_STATUS_SIZE,0:nmrank-1), stat=ierr)
       allocate(wk_recvcnt(0:nmrank-1), stat=ierr)
       allocate(wk_sendcnt(0:nmrank-1), stat=ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_FFT_CD_Inverse_3D :  Not allocate '
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 68, ierr)
        endif

       if(allocated(x2z_recv)) deallocate(x2z_recv)
       if(allocated(x2z_send)) deallocate(x2z_send)
       allocate(x2z_recv(2,0:nmrank-1))
       allocate(x2z_send(2,0:nmrank-1))
       x2z_send = 0
       x2z_recv = 0
       wk_recvcnt = 0
       wk_sendcnt = 0
       do i = 1, nel_fftcd_z(myrank)
          irank = map_fftcd_x(mp_fftcd_z(i)) - 1
          wk_recvcnt(irank) = wk_recvcnt(irank) + 1
       enddo
       do i = 1, nel_fftcd_x(myrank)
          irank = map_fftcd_z(mp_fftcd_x(i)) - 1
          wk_sendcnt(irank) = wk_sendcnt(irank) + 1
       enddo
       k = 0
       do i = 0, nmrank - 1
          if(wk_recvcnt(i) /= 0) then
             k = k + 1
             x2z_recv(1,i) = wk_recvcnt(i)
             x2z_recv(2,i) = k
          endif
       enddo
       x2z_rrank = k
       k = 0
       do i = 0, nmrank - 1
          if(wk_sendcnt(i) /= 0) then
             k = k + 1
             x2z_send(1,i) = wk_sendcnt(i)
             x2z_send(2,i) = k
          endif
       enddo
       x2z_srank = k
       x2z_rmax = maxval(wk_recvcnt)
       x2z_smax = maxval(wk_sendcnt)
       x2z_srmax = max(x2z_rmax,x2z_smax)

       if(allocated(z2y_recv)) deallocate(z2y_recv)
       if(allocated(z2y_send)) deallocate(z2y_send)
       allocate(z2y_recv(2,0:nmrank-1))
       allocate(z2y_send(2,0:nmrank-1))
       z2y_send = 0
       z2y_recv = 0
       wk_recvcnt = 0
       wk_sendcnt = 0

       nx = xyz_fftcd_y(2,1)-xyz_fftcd_y(1,1)+1
       ny = xyz_fftcd_y(2,2)-xyz_fftcd_y(1,2)+1
       nz = xyz_fftcd_y(2,3)-xyz_fftcd_y(1,3)+1
       allocate(wk_mp_fft_y(nx*ny*nz) ,stat=ierr)
       do k = 0, nz-1
          do j = 0, ny-1
             do i = 0, nx-1
                wk_mp_fft_y(i+j*nx+k*nx*ny+1) = mp_fftcd_y(i+k*nx+j*nx*nz+1)
             enddo
          enddo
       enddo

       do i = 1, nel_fftcd_y(myrank)
          irank = map_fftcd_z(wk_mp_fft_y(i)) - 1
          wk_recvcnt(irank) = wk_recvcnt(irank) + 1
       enddo
       do i = 1, nel_fftcd_z(myrank)
          irank = map_fftcd_y(mp_fftcd_z(i)) - 1
          wk_sendcnt(irank) = wk_sendcnt(irank) + 1
       enddo
       k = 0
       do i = 0, nmrank - 1
          if(wk_recvcnt(i) /= 0) then
             k = k + 1
             z2y_recv(1,i) = wk_recvcnt(i)
             z2y_recv(2,i) = k
          endif
       enddo
       z2y_rrank = k
       k = 0
      do i = 0, nmrank - 1
          if(wk_sendcnt(i) /= 0) then
             k = k + 1
             z2y_send(1,i) = wk_sendcnt(i)
             z2y_send(2,i) = k
          endif
       enddo
       z2y_srank = k
       z2y_rmax = maxval(wk_recvcnt)
       z2y_smax = maxval(wk_sendcnt)
       z2y_srmax = max(z2y_rmax,z2y_smax)

       deallocate(wk_mp_fft_y)

#ifdef FFT_ALLTOALL
       max_send(1) = x2z_rmax
       max_send(2) = x2z_smax
       call mpi_allreduce(max_send,max_recv,2,mpi_integer,mpi_max,mpi_fftcd_xz_world,ierr)
       x2z_rmax = max_recv(1)
       x2z_smax = max_recv(2)
       max_send(1) = z2y_rmax
       max_send(2) = z2y_smax
       call mpi_allreduce(max_send,max_recv,2,mpi_integer,mpi_max,mpi_fftcd_zy_world,ierr)
       z2y_rmax = max_recv(1)
       z2y_smax = max_recv(2)
#endif

       firstcall_cd_inverse_3d = .false.
#ifdef __TIMER_DO__
  call timer_end(173)
#endif
    endif

#ifdef FFT_ALLTOALL
    if (iesize /= savesize) then
#else
    if (iesize > savesize) then
#endif
#ifdef __TIMER_DO__
  call timer_sta(174)
#endif
       if (allocated(wk_recv1)) deallocate(wk_recv1)
       if (allocated(wk_send1)) deallocate(wk_send1)
       if (allocated(wk_recv2)) deallocate(wk_recv2)
       if (allocated(wk_send2)) deallocate(wk_send2)
       if (allocated(wk_afft_y)) deallocate(wk_afft_y)
       allocate(wk_recv1(x2z_rmax*kimg*iesize,x2z_rrank), stat=ierr)
       allocate(wk_send1(x2z_smax*kimg*iesize,x2z_srank), stat=ierr)
       allocate(wk_recv2(z2y_rmax*kimg*iesize,z2y_rrank), stat=ierr)
       allocate(wk_send2(z2y_smax*kimg*iesize,z2y_srank), stat=ierr)
       allocate(wk_afft_y(max_elm*kimg,iesize) ,stat=ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_FFT_CD_Inverse_3D :  Not allocate '
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 69, ierr)
        endif
       savesize = iesize
#ifdef __TIMER_DO__
  call timer_end(174)
#endif
    end if
!
! X-axis (y-z div)
!
    nx = xyz_fftcd_x(2,1)-xyz_fftcd_x(1,1)+1
    ny = xyz_fftcd_x(2,2)-xyz_fftcd_x(1,2)+1
    nz = xyz_fftcd_x(2,3)-xyz_fftcd_x(1,3)+1
#ifdef __TIMER_DO__
  call timer_sta(260)
#endif
    if (kimg==1) then
!.     do ib = 1, iesize
!.        do k = 1, nz
!.           do j = 1, ny
!.              call DM_V1DRCF2(wk_afft_l(1+nx*(j-1)+nx*ny*(k-1),ib),nx-2,      &
!.             &                wk_afft_l(1+nx*(j-1)+nx*ny*(k-1),ib),1,1,ierr)
!.           enddo
!.        enddo
!.     enddo
       nsize(1:3) = (/nx-2,ny,nz/)
       isin(1:3)  = (/1,0,0/)
!OCL INDEPENDENT[ DVMRF2 ]
!OCL PARALLEL_STRONG
       do ib = 1, iesize
          call DVMRF2(wk_afft_l(1,ib),nsize,3,isin,1,ierr)
       end do
    else
       do ib = 1, iesize
          call DM_V1DMCFT(wk_afft_l(1,ib),nx,nx,ny*nz,-1,ierr)
       end do
!      call DM_V1DMCFT(wk_afft_l(1,1),nx,nx,ny*nz*iesize,-1,ierr)
    end if
#ifdef __TIMER_DO__
  call timer_end(260)
#endif

#ifdef __TIMER_DO__
  call timer_sta(175)
#endif
    do l = 1, fftcd_Z_x_dim
       nis = nis_fftcd_Z_x(l)
       nie = nie_fftcd_Z_x(l)
       nnx = nie - nis + 1
       if (kimg == 1) then
          do ib = 1, iesize
             iadd0 = nnx*ny*nz*(ib-1)
             do k = 1, nz
                do j = 1, ny
                   do i = nis, nie
                      wk_send1(iadd0+(i-nis+1+(j-1)*nnx+(k-1)*nnx*ny),l) = wk_afft_l(i+(j-1)*nx+(k-1)*nx*ny,ib)
                   end do
                end do
             end do
          end do
       else
          do ib = 1, iesize
             iadd0 = nnx*ny*nz*(ib-1)*2
             do k = 1, nz
                do j = 1, ny
                   do i = nis, nie
                      iadd  = (i+(j-1)*nx+(k-1)*nx*ny)*2
                      iadd1 = iadd0+(i-nis+1+(j-1)*nnx+(k-1)*nnx*ny)*2
!fj --------------------
!                     do ri = 1, kimg
!                       wk_send1(iadd1+(ri-kimg),l) = wk_afft_l(iadd+(ri-kimg),ib)
!                     end do
                      wk_send1(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send1(iadd1  ,l) = wk_afft_l(iadd  ,ib)
!fj --------------------
                   end do
                end do
             end do
          end do
       end if
    end do
#ifdef __TIMER_DO__
  call timer_end(175)
#endif

#ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_fftcd_xz_world)
  call timer_sta(176)
#endif

    call MPI_ALLTOALL(wk_send1, x2z_smax*kimg*iesize, mpi_double_precision,   &
   &                  wk_recv1, x2z_rmax*kimg*iesize, mpi_double_precision,   &
   &                                                  mpi_fftcd_xz_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_CD_Inverse_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 70, ierr)
     endif

#ifdef __TIMER_DO__
  call timer_sta(177)
#endif
#ifdef __TIMER_DO__
  call timer_end(177)
#endif

#else
!else ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_comm)
  call timer_sta(176)
#endif

    icnt_recv = 0
    do lrank = 0, nmrank - 1
       if ((x2z_recv(1,lrank) /= 0)) then
          call mpi_irecv(wk_recv1(1,x2z_recv(2,lrank)), x2z_recv(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_r(icnt_recv), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_CD_Inverse_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 70, ierr)
           endif
          icnt_recv = icnt_recv + 1
       endif
    enddo

    icnt_send = 0
    do lrank = 0, nmrank - 1
       if ((x2z_send(1,lrank) /= 0)) then
          call mpi_isend(wk_send1(1,x2z_send(2,lrank)), x2z_send(1,lrank)*kimg*iesize, &
                         mpi_double_precision, lrank, itag, mpi_comm, req_s(icnt_send), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_CD_Inverse_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 71, ierr)
           endif
          icnt_send = icnt_send + 1
       endif
    enddo

#ifdef __TIMER_DO__
  call timer_sta(177)
#endif
#ifdef __TIMER_DO__
  call timer_end(177)
#endif

    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_CD_Inverse_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 72, ierr)
     endif
    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_CD_Inverse_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 73, ierr)
     endif
#endif
!endif ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_end(176)
#endif

#ifdef __TIMER_DO__
  call timer_sta(178)
#endif
    nx = xyz_fftcd_z(2,1)-xyz_fftcd_z(1,1)+1
    ny = xyz_fftcd_z(2,2)-xyz_fftcd_z(1,2)+1
    nz = xyz_fftcd_z(2,3)-xyz_fftcd_z(1,3)+1
    if (kimg == 1) then
       do ib = 1, iesize
          do j = 1, ny
             do i = 1, nx/2
                do l = 1, fftcd_X_z_dim
                   nis = nis_fftcd_X_z(l)
                   nie = nie_fftcd_X_z(l)
                   nnz = nie - nis + 1
                   iadd0 = nx/2*ny*nnz*(ib-1)*2
                   do k = nis, nie
                      iadd1 = iadd0 + (i+(j-1)*nx/2+(k-nis)*nx/2*ny)*2
                      iadd  = (k+(i-1)*nz+(j-1)*nz*nx/2)*2
                      wk_afft_l(iadd-1,ib) = wk_recv1(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv1(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end do
    else
       do ib = 1, iesize
          do j = 1, ny
             do i = 1, nx
                do l = 1, fftcd_X_z_dim
                   nis = nis_fftcd_X_z(l)
                   nie = nie_fftcd_X_z(l)
                   nnz = nie - nis + 1
                   iadd0 = nx*ny*nnz*(ib-1)*2
                   do k = nis, nie
                      iadd1 = iadd0+(i+(j-1)*nx+(k-nis)*nx*ny)*2
                      iadd  = (k+(i-1)*nz+(j-1)*nz*nx)*2
                      wk_afft_l(iadd-1,ib) = wk_recv1(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv1(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end do
    end if
#ifdef __TIMER_DO__
  call timer_end(178)
#endif

!
! Z-axis (x-y div)
!
#ifdef __TIMER_DO__
  call timer_sta(261)
#endif
    if (kimg == 1) then
       do ib = 1, iesize
          call DM_V1DMCFT(wk_afft_l(1,ib),nz,nz,nx*ny/2,1,ierr)
       end do
!.     call DM_V1DMCFT(wk_afft_l(1,ib),nz,nz,nx*ny/2*iesize,-1,ierr)
    else
       do ib = 1, iesize
          call DM_V1DMCFT(wk_afft_l(1,ib),nz,nz,nx*ny,-1,ierr)
       end do
!.     call DM_V1DMCFT(wk_afft_l(1,1),nz,nz,nx*ny*iesize,-1,ierr)
    end if
#ifdef __TIMER_DO__
  call timer_end(261)
#endif

    nx = xyz_fftcd_y(2,1)-xyz_fftcd_y(1,1)+1
    ny = xyz_fftcd_y(2,2)-xyz_fftcd_y(1,2)+1
    nz = xyz_fftcd_y(2,3)-xyz_fftcd_y(1,3)+1
!XX allocate(wk_mp_fft_y(nx*ny*nz) ,stat=ierr)
!XX  if (ierr /= 0) then
!XX     write(nfout,*)' m_FFT_CD_Inverse_3D :  Not allocate '
!XX     call flush(nfout)
!XX     call mpi_abort(mpi_comm_world, 74, ierr)
!XX  endif
#ifdef __TIMER_DO__
  call timer_sta(179)
#endif
!XX do k = 0, nz-1
!XX    do j = 0, ny-1
!XX       do i = 0, nx-1
!XX          wk_mp_fft_y(i+j*nx+k*nx*ny+1) = mp_fftcd_y(i+k*nx+j*nx*nz+1)
!XX       enddo
!XX    enddo
!XX enddo
#ifdef __TIMER_DO__
  call timer_end(179)
#endif

!F  allocate(wk_recv2(z2y_rmax*kimg*iesize,z2y_rrank), stat=ierr)
!F  allocate(wk_send2(z2y_smax*kimg*iesize,z2y_srank), stat=ierr)
!F   if (ierr /= 0) then
!F      write(nfout,*)' m_FFT_CD_Inverse_3D :  Not allocate '
!F      call flush(nfout)
!F      call mpi_abort(mpi_comm_world, 75, ierr)
!F   endif

#ifdef __TIMER_DO__
  call timer_sta(180)
#endif
    nx = xyz_fftcd_z(2,1)-xyz_fftcd_z(1,1)+1
    ny = xyz_fftcd_z(2,2)-xyz_fftcd_z(1,2)+1
    nz = xyz_fftcd_z(2,3)-xyz_fftcd_z(1,3)+1
    if (kimg == 1) then
       do ib = 1, iesize
          do j = 1, ny
             do i = 1, nx/2
                do l = 1, fftcd_Y_z_dim
                   nis = nis_fftcd_Y_z(l)
                   nie = nie_fftcd_Y_z(l)
                   nnz = nie - nis + 1
                   iadd0 = nx/2*ny*nnz*(ib-1)*2
                   do k = nis, nie
                      iadd  = (k+(i-1)*nz+(j-1)*nz*nx/2)*2
                      iadd1 = iadd0 + (k-nis+1+(i-1)*nnz+(j-1)*nnz*nx/2)*2
                      wk_send2(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send2(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
    else
       do ib = 1, iesize
          do j = 1, ny
             do i = 1, nx
                do l = 1, fftcd_Y_z_dim
                   nis = nis_fftcd_Y_z(l)
                   nie = nie_fftcd_Y_z(l)
                   nnz = nie - nis + 1
                   iadd0 = nx*ny*nnz*(ib-1)*2
                   do k = nis, nie
                      iadd  = (k+(i-1)*nz+(j-1)*nz*nx)*2
                      iadd1 = iadd0+(k-nis+1+(i-1)*nnz+(j-1)*nnz*nx)*2
!fj --------------------
!                     do ri = 1, kimg
!                        wk_send2(iadd1+(ri-kimg),l) = wk_afft_l(iadd+(ri-kimg),ib)
!                     end do
                      wk_send2(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send2(iadd1  ,l) = wk_afft_l(iadd  ,ib)
!fj --------------------
                   end do
                end do
             end do
          end do
       end do
    end if
#ifdef __TIMER_DO__
  call timer_end(180)
#endif

#ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_fftcd_zy_world)
  call timer_sta(181)
#endif

    call MPI_ALLTOALL(wk_send2, z2y_smax*kimg*iesize, mpi_double_precision,   &
   &                  wk_recv2, z2y_rmax*kimg*iesize, mpi_double_precision,   &
   &                                                  mpi_fftcd_zy_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_CD_Inverse_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 76, ierr)
     endif

#ifdef __TIMER_DO__
  call timer_sta(182)
#endif
#ifdef __TIMER_DO__
  call timer_end(182)
#endif

#else
!else ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_comm)
  call timer_sta(181)
#endif

    icnt_recv = 0
    do lrank = 0, nmrank - 1
       if ((z2y_recv(1,lrank) /= 0)) then
          call mpi_irecv(wk_recv2(1,z2y_recv(2,lrank)), z2y_recv(1,lrank)*kimg*iesize, &
                         mpi_double_precision, lrank, itag, mpi_comm, req_r(icnt_recv), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_CD_FFT_Inverse_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 76, ierr)
           endif
          icnt_recv = icnt_recv + 1
       endif
    enddo

    icnt_send = 0
    do lrank = 0, nmrank - 1
       if ((z2y_send(1,lrank) /= 0)) then
          call mpi_isend(wk_send2(1,z2y_send(2,lrank)), z2y_send(1,lrank)*kimg*iesize, &
                         mpi_double_precision, lrank, itag, mpi_comm, req_s(icnt_send), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_CD_FFT_Inverse_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 77, ierr)
           endif
          icnt_send = icnt_send + 1
       endif
    enddo

#ifdef __TIMER_DO__
  call timer_sta(182)
#endif
#ifdef __TIMER_DO__
  call timer_end(182)
#endif

    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_CD_Inverse_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 78, ierr)
     endif
    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_CD_Inverse_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 79, ierr)
     endif

#endif
!endif ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_end(181)
#endif

#ifdef __TIMER_DO__
  call timer_sta(183)
#endif
    nx = xyz_fftcd_y(2,1)-xyz_fftcd_y(1,1)+1
    ny = xyz_fftcd_y(2,2)-xyz_fftcd_y(1,2)+1
    nz = xyz_fftcd_y(2,3)-xyz_fftcd_y(1,3)+1
    if (kimg == 1) then
       do ib = 1, iesize
          do k = 1, nz
             do i = 1, nx/2
                do l = 1, fftcd_Z_y_dim
                   nis = nis_fftcd_Z_y(l)
                   nie = nie_fftcd_Z_y(l)
                   nny = nie - nis + 1
                   iadd0 = nx/2*nny*nz*(ib-1)*2
                   do j = nis, nie
                      iadd1 = iadd0 + (k+(i-1)*nz+(j-nis)*nz*nx/2)*2
                      iadd  = (j+(i-1)*ny+(k-1)*ny*nx/2)*2
                      wk_afft_y(iadd-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_y(iadd  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end do
    else
       do ib = 1, iesize
          do k = 1, nz
             do i = 1, nx
                do l = 1, fftcd_Z_y_dim
                   nis = nis_fftcd_Z_y(l)
                   nie = nie_fftcd_Z_y(l)
                   nny = nie - nis + 1
                   iadd0 = nx*nny*nz*(ib-1)*2
                   do j = nis, nie
                      iadd1 = iadd0+(k+(i-1)*nz+(j-nis)*nz*nx)*2
                      iadd  = (j+(i-1)*ny+(k-1)*ny*nx)*2
                      wk_afft_y(iadd-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_y(iadd  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end do
    end if
#ifdef __TIMER_DO__
  call timer_end(183)
#endif


!
! Y-axis (z-x div)
!
#ifdef __TIMER_DO__
  call timer_sta(262)
#endif
    if (kimg == 1) then
       do ib = 1, iesize
          call DM_V1DMCFT(wk_afft_y(1,ib),ny,ny,nx/2*nz,1,ierr)
       end do
!.     call DM_V1DMCFT(wk_afft_y(1,1),ny,ny,nx/2*ny*iesize,1,ierr)
    else
       do ib = 1, iesize
          call DM_V1DMCFT(wk_afft_y(1,ib),ny,ny,nx*nz,-1,ierr)
       end do
!.     call DM_V1DMCFT(wk_afft_l(1,1),nz,nz,nx*ny*iesize,-1,ierr)
    end if
#ifdef __TIMER_DO__
  call timer_end(262)
#endif

#ifdef __TIMER_DO__
  call timer_sta(184)
#endif
! (x,z,y) <- (y,x,z)
    if (kimg == 1) then
       do ib = 1, iesize
          do k = 1, nz
             do j = 1, ny
                do i = 1, nx/2
                   iadd0 = (j+(i-1)*ny+(k-1)*ny*nx/2)*2
                   iadd1 = (i+(k-1)*nx/2+(j-1)*nx/2*nz)*2
                   wk_afft_l(iadd1-1,ib) = wk_afft_y(iadd0-1,ib)
                   wk_afft_l(iadd1  ,ib) = wk_afft_y(iadd0  ,ib)
                enddo
             enddo
          enddo
       enddo
    else
       do ib = 1, iesize
          do j = 0, ny-1
             do k = 0, nz-1
                do i = 0, nx-1
                   wk_afft_l((i+k*nx+j*nx*nz+1)*2-1,ib) = wk_afft_y((j+i*ny+k*ny*nx+1)*2-1,ib)
                   wk_afft_l((i+k*nx+j*nx*nz+1)*2  ,ib) = wk_afft_y((j+i*ny+k*ny*nx+1)*2  ,ib)
                enddo
             enddo
          enddo
       enddo
    end if
#ifdef __TIMER_DO__
  call timer_end(184)
#endif

1000 continue

!F  if (allocated(req_r))       deallocate(req_r)
!F  if (allocated(req_s))       deallocate(req_s)
!F  if (allocated(sta_r))       deallocate(sta_r)
!F  if (allocated(sta_s))       deallocate(sta_s)
!F  if (allocated(wk_recvcnt))    deallocate(wk_recvcnt)
!F  if (allocated(wk_sendcnt))    deallocate(wk_sendcnt)
!F  if (allocated(wk_afft_y))   deallocate(wk_afft_y)
    if (allocated(wk_mp_fft_y)) deallocate(wk_mp_fft_y)

#ifdef __TIMER_SUB__
    call timer_end(110)
#endif
  end subroutine m_FFT_CD_Inverse_3D
!------------------------------------------------------------------------------
#endif
!endif ifndef FFT_USE_SSL2

#ifdef FFT_USE_SSL2
!------------------------------------------------------------------------------
  subroutine m_FFT_Direct_XYZ_3D(nfout,wk_afft_l, wk_size, iesize)
    integer, intent(in) :: nfout, wk_size, iesize
    real(kind=DP), dimension(wk_size*kimg,iesize), intent(inout) :: wk_afft_l

    integer,      allocatable, dimension(:)   :: wk_recvdsp
    integer, parameter :: FFTW_MEASURE=0
    integer, parameter :: FFTW_ESTIMATE=64
    integer, parameter :: FFTW_RANK=1
    integer            :: FFTW_FLAG
    integer,dimension(FFTW_RANK) :: NFFTW3, NEMBED, NEREAL
    integer :: nx, ny, nz, i, j, k, l, ri, ib, iadd, lrank
    integer :: irank, itag, icnt_send, icnt_recv
    integer,save, allocatable, dimension(:)   :: req_r, req_s
    integer,save, allocatable, dimension(:,:) :: sta_r, sta_s
    integer,save, allocatable, dimension(:)   :: wk_recvcnt, wk_sendcnt
    real(kind=DP),allocatable, dimension(:,:),save :: wk_recv1, wk_send1
    real(kind=DP),allocatable, dimension(:,:),save :: wk_recv2, wk_send2
    real(kind=DP),allocatable, dimension(:,:,:)   :: wk_gather
    integer :: mpi_comm, myrank, nmrank
    integer(kind=8),save :: planx1 = 0, planx2 = 0
    integer(kind=8),save :: plany1 = 0, plany2 = 0
    integer(kind=8),save :: planz1 = 0, planz2 = 0
    integer(kind=8),save :: planx1_1d = 0, planx2_1d = 0
    integer(kind=8),save :: plany1_1d = 0, plany2_1d = 0
    integer(kind=8),save :: planz1_1d = 0, planz2_1d = 0
    integer(kind=8) :: plan
    integer :: max_x, max_y, max_z
    integer,save :: max_elm = 0, savesize = 0, expo_x = 0, expo_y = 0, expo_z = 0

    integer,save, allocatable, dimension(:,:) :: z2y_recv, z2y_send, y2x_recv, y2x_send
    integer,save :: z2y_rrank, z2y_srank, y2x_rrank, y2x_srank
    integer,save :: z2y_rmax, z2y_smax, z2y_srmax, y2x_rmax, y2x_smax, y2x_srmax

    integer ::  nnx, nny, nnz, iadd0, iadd1, nis, nie, jadd, kadd, nisx, niex, nxp, nyp, nzp

    integer,dimension(3) :: nsize, isin
    real(kind=DP),allocatable, dimension(:),save :: wwx,wwy,wwz
    integer :: isw
#ifdef FFT_ALLTOALL
    integer,dimension(2) :: max_send,max_recv
#endif

#ifdef __TIMER_SUB__
    call mpi_barrier(mpi_ke_world,ierr)
    call timer_sta(105)
#endif

    if (kimg == 1) then
       FFTW_FLAG = +1      ! FFTW_BACKWARD
    else
       FFTW_FLAG = -1      ! FFTW_FORWARD
    endif
    plan = 0

    mpi_comm = mpi_ke_world
    myrank = myrank_g
    nmrank = nrank_g
    itag = 10

    if (np_fft_x==0)then
       nx = 0
       nx = 0
       nz = 0
       go to 1000
    endif

    if (firstcall_direct_xyz_3d) then
       savesize = 0
#ifdef __TIMER_DO__
  call timer_sta(118)
#endif
       max_x = maxval(nel_fft_x(:))
       max_y = maxval(nel_fft_y(:))
       max_z = maxval(nel_fft_z(:))
       max_elm = max(max_x,max_y,max_z)

       if(allocated(req_r)) deallocate(req_r)
       if(allocated(req_s)) deallocate(req_s)
       if(allocated(sta_r)) deallocate(sta_r)
       if(allocated(sta_s)) deallocate(sta_s)
       if(allocated(wk_recvcnt)) deallocate(wk_recvcnt)
       if(allocated(wk_sendcnt)) deallocate(wk_sendcnt)
       allocate(req_r(0:nmrank-1), stat=ierr)
       allocate(req_s(0:nmrank-1), stat=ierr)
       allocate(sta_r(MPI_STATUS_SIZE,0:nmrank-1), stat=ierr)
       allocate(sta_s(MPI_STATUS_SIZE,0:nmrank-1), stat=ierr)
       allocate(wk_recvcnt(0:nmrank-1), stat=ierr)
       allocate(wk_sendcnt(0:nmrank-1), stat=ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_FFT_Direct_XYZ_3D :  Not allocate '
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 51, ierr)
        endif

       if(allocated(z2y_recv)) deallocate(z2y_recv)
       if(allocated(z2y_send)) deallocate(z2y_send)
       allocate(z2y_recv(2,0:nmrank-1))
       allocate(z2y_send(2,0:nmrank-1))
       z2y_send = 0
       z2y_recv = 0
       wk_recvcnt = 0
       wk_sendcnt = 0
       do i = 1, nel_fft_y(myrank)
          irank = map_fft_z(mp_fft_y(i)) - 1
          wk_recvcnt(irank) = wk_recvcnt(irank) + 1
       enddo
       do i = 1, nel_fft_z(myrank)
          irank = map_fft_y(mp_fft_z(i)) - 1
          wk_sendcnt(irank) = wk_sendcnt(irank) + 1
       enddo
       k = 0
       do i = 0, nmrank - 1
          if(wk_recvcnt(i) /= 0) then
             k = k + 1
             z2y_recv(1,i) = wk_recvcnt(i)
             z2y_recv(2,i) = k
          endif
       enddo
       z2y_rrank = k
       k = 0
       do i = 0, nmrank - 1
          if(wk_sendcnt(i) /= 0) then
             k = k + 1
             z2y_send(1,i) = wk_sendcnt(i)
             z2y_send(2,i) = k
          endif
       enddo
       z2y_srank = k
       z2y_rmax = maxval(wk_recvcnt)
       z2y_smax = maxval(wk_sendcnt)
       z2y_srmax = max(z2y_rmax,z2y_smax)

       if(allocated(y2x_recv)) deallocate(y2x_recv)
       if(allocated(y2x_send)) deallocate(y2x_send)
       allocate(y2x_recv(2,0:nmrank-1))
       allocate(y2x_send(2,0:nmrank-1))
       y2x_send = 0
       y2x_recv = 0
       wk_recvcnt = 0
       wk_sendcnt = 0
       do i = 1, nel_fft_x(myrank)
          irank = map_fft_y(mp_fft_x(i)) - 1
          wk_recvcnt(irank) = wk_recvcnt(irank) + 1
       enddo
       do i = 1, nel_fft_y(myrank)
          irank = map_fft_x(mp_fft_y(i)) - 1
          wk_sendcnt(irank) = wk_sendcnt(irank) + 1
       enddo
       k = 0
       do i = 0, nmrank - 1
          if(wk_recvcnt(i) /= 0) then
             k = k + 1
             y2x_recv(1,i) = wk_recvcnt(i)
             y2x_recv(2,i) = k
          endif
       enddo
       y2x_rrank = k
       k = 0
      do i = 0, nmrank - 1
          if(wk_sendcnt(i) /= 0) then
             k = k + 1
             y2x_send(1,i) = wk_sendcnt(i)
             y2x_send(2,i) = k
          endif
       enddo
       y2x_srank = k
       y2x_rmax = maxval(wk_recvcnt)
       y2x_smax = maxval(wk_sendcnt)
       y2x_srmax = max(y2x_rmax,y2x_smax)

       nx = xyz_fft_x(2,1)-xyz_fft_x(1,1)+1
       ny = xyz_fft_y(2,2)-xyz_fft_y(1,2)+1
       nz = xyz_fft_z(2,3)-xyz_fft_z(1,3)+1
       if ((nx==64).or.(nx==128).or.(nx==256).or.(nx==512).or.(nx==1024).or.(nx==2048)) then
          expo_x = 1
       end if
       if ((ny==64).or.(ny==128).or.(ny==256).or.(ny==512).or.(ny==1024).or.(ny==2048)) then
          expo_y = 1
       end if
       if ((nz==64).or.(nz==128).or.(nz==256).or.(nz==512).or.(nz==1024).or.(nz==2048)) then
          expo_z = 1
       end if
#ifdef FFT_ALLTOALL
       max_send(1) = z2y_rmax
       max_send(2) = z2y_smax
       call mpi_allreduce(max_send,max_recv,2,mpi_integer,mpi_max,mpi_fft_yz_world,ierr)
       z2y_rmax = max_recv(1)
       z2y_smax = max_recv(2)
       max_send(1) = y2x_rmax
       max_send(2) = y2x_smax
       call mpi_allreduce(max_send,max_recv,2,mpi_integer,mpi_max,mpi_fft_xy_world,ierr)
       y2x_rmax = max_recv(1)
       y2x_smax = max_recv(2)
#endif
       if (ipri > 1) then
          write(nfout,'("m_FFT_Direct_XYZ_3D   --   myrank_g=",i4)') myrank_g
          write(nfout,'("z2y_send")')
          write(nfout,'(10(i8,", "))') (z2y_send(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (z2y_send(2,i),i=0,nmrank-1)
          write(nfout,'("z2y_recv")')
          write(nfout,'(10(i8,", "))') (z2y_recv(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (z2y_recv(2,i),i=0,nmrank-1)
          write(nfout,'("y2x_send")')
          write(nfout,'(10(i8,", "))') (y2x_send(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (y2x_send(2,i),i=0,nmrank-1)
          write(nfout,'("y2x_recv")')
          write(nfout,'(10(i8,", "))') (y2x_recv(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (y2x_recv(2,i),i=0,nmrank-1)
       endif

       firstcall_direct_xyz_3d = .false.
#ifdef __TIMER_DO__
  call timer_end(118)
#endif
    endif

#ifdef FFT_ALLTOALL
    if (iesize /= savesize) then
#else
    if (iesize > savesize) then
#endif
#ifdef __TIMER_DO__
  call timer_sta(120)
#endif
       if (allocated(wk_recv1)) deallocate(wk_recv1)
       if (allocated(wk_send1)) deallocate(wk_send1)
       if (allocated(wk_recv2)) deallocate(wk_recv2)
       if (allocated(wk_send2)) deallocate(wk_send2)
       allocate(wk_recv1(z2y_rmax*kimg*iesize,z2y_rrank), stat=ierr)
       allocate(wk_send1(z2y_smax*kimg*iesize,z2y_srank), stat=ierr)
       allocate(wk_recv2(y2x_rmax*kimg*iesize,y2x_rrank), stat=ierr)
       allocate(wk_send2(y2x_smax*kimg*iesize,y2x_srank), stat=ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_FFT_Direct_XYZ_3D :  Not allocate '
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 52, ierr)
        endif
       savesize = iesize
#ifdef __TIMER_DO__
  call timer_end(120)
#endif
    end if

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP          SHARED(     nfout,wk_afft_l,iesize,req_r,req_s,sta_r,sta_s,    &
!$OMP                 wk_recv1,wk_send1,wk_recv2,wk_send2,                    &
!$OMP                 z2y_recv, z2y_send, y2x_recv, y2x_send,itag,            &
!$OMP                 xyz_fft_y,xyz_fft_z,xyz_fft_x,fft_Z_x_dim,nis_fft_Z_x,  &
!$OMP                 nie_fft_Z_x,fft_X_z_dim,nis_fft_X_z,nie_fft_X_z,        &
!$OMP                 fft_Y_z_dim,nis_fft_Y_z,nie_fft_Y_z,fft_Z_y_dim,        &
!$OMP                 nie_fft_Z_y,nis_fft_Z_y,                                &
!$OMP                 nie_fft_X_y,nis_fft_X_y,fft_X_y_dim,                    &
!$OMP                 nie_fft_Y_x,nis_fft_Y_x,fft_Y_x_dim,                    &
#ifdef FFT_ALLTOALL
!$OMP                 mpi_fft_yz_world, z2y_rmax, z2y_smax,                   &
!$OMP                 mpi_fft_xy_world, y2x_rmax, y2x_smax,                   &
#endif
!$OMP                 mpi_comm,myrank,nmrank,kimg,expo_x,expo_y,expo_z    )   &
!$OMP          PRIVATE(nx,ny,nz,i,j,k,l,ri,ib,iadd,lrank,icnt_send,icnt_recv, &
!$OMP                  nxp,nyp,nzp,nisx,niex,nsize,isin,ierr,wwx,wwy,wwz,isw, &
!$OMP                  irank,nnx,nny,nnz,iadd0,iadd1,nis,nie,jadd,kadd  )
    nx = xyz_fft_z(2,1)-xyz_fft_z(1,1)+1
    ny = xyz_fft_z(2,2)-xyz_fft_z(1,2)+1
    nz = xyz_fft_z(2,3)-xyz_fft_z(1,3)+1
    nzp = nz
!
! Z-axis (x-y div)
!
#ifdef __TIMER_DO__
  call timer_sta(251)
#endif
    if(kimg==1) then
       nsize(1:3) = (/nz,(nx/2),ny/)
       isin(1:3)  = (/-1,0,0/)
!$OMP DO
       do ib = 1, iesize
          call DVMCF2(wk_afft_l(1,ib),nsize,3,isin,ierr)
!x     end do
! (z,x,y) -> (z,x,y)
!xOMP DO
!x     do ib = 1, iesize
          do j = 1, ny
             do i = 1, nx/2
                do l = 1, fft_Y_z_dim
                   nis = nis_fft_Y_z(l)
                   nie = nie_fft_Y_z(l)
                   nnz = nie - nis + 1
                   iadd0 = nx/2*ny*nnz*(ib-1)*2
                   do k = nis, nie
                      iadd  = (k+(i-1)*nzp+(j-1)*nzp*nx/2)*2
                      iadd1 = iadd0+((k-nis+1)+(i-1)*nnz+(j-1)*nnz*nx/2)*2
                      wk_send1(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send1(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
    else
     if (expo_z==1) then
        allocate(wwz((2*nz+70)*kimg))
     endif
     if (iesize > 1) then
       nsize(1:3) = (/nz,nx,ny/)
       isin(1:3)  = (/1,0,0/)
       isw = 1
!$OMP DO
       do ib = 1, iesize
          if (expo_z==1) then
             do j = 1, ny
                do i = 1, nx
                   call DVCFM1(wk_afft_l((1+(i-1)*nz+(j-1)*nz*nx)*2-1,ib),nz,isw,1,wwz,ierr)
!x                 call DVMCF2(wk_afft_l((1+(i-1)*nz+(j-1)*nz*nx)*2-1,ib),nz,1,1,ierr)
                end do
             end do
          else
             call DVMCF2(wk_afft_l(1,ib),nsize,3,isin,ierr)
          end if
!x     end do
! (z,x,y) -> (z,x,y)
!xOMP DO
!x     do ib = 1, iesize
          do j = 1, ny
             do i = 1, nx
                do l = 1, fft_Y_z_dim
                   nis = nis_fft_Y_z(l)
                   nie = nie_fft_Y_z(l)
                   nnz = nie - nis + 1
                   iadd0 = nx*ny*nnz*(ib-1)*2
                   do k = nis, nie
                      iadd  = (k+(i-1)*nzp+(j-1)*nzp*nx)*2
                      iadd1 = iadd0+((k-nis+1)+(i-1)*nnz+(j-1)*nnz*nx)*2
                      wk_send1(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send1(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
     else
       isw = 1
       do ib = 1, iesize
          if (expo_z==1) then
!$OMP DO
             do j = 1, ny
                do i = 1, nx
                   call DVCFM1(wk_afft_l((1+(i-1)*nz+(j-1)*nz*nx)*2-1,ib),nz,isw,1,wwz,ierr)
                end do
             end do
          else
!$OMP DO
             do j = 1, ny
                do i = 1, nx
                   call DVMCF2(wk_afft_l((1+(i-1)*nz+(j-1)*nz*nx)*2-1,ib),nz,1,1,ierr)
                end do
             end do
          end if
       end do
! (z,x,y) -> (z,x,y)
       do ib = 1, iesize
!$OMP DO
          do j = 1, ny
             do i = 1, nx
                do l = 1, fft_Y_z_dim
                   nis = nis_fft_Y_z(l)
                   nie = nie_fft_Y_z(l)
                   nnz = nie - nis + 1
                   iadd0 = nx*ny*nnz*(ib-1)*2
                   do k = nis, nie
                      iadd  = (k+(i-1)*nzp+(j-1)*nzp*nx)*2
                      iadd1 = iadd0+((k-nis+1)+(i-1)*nnz+(j-1)*nnz*nx)*2
                      wk_send1(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send1(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
     end if
     if (expo_z==1) then
        deallocate(wwz)
     end if
    end if
#ifdef __TIMER_DO__
  call timer_end(251)
#endif

#ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_fft_yz_world)
  call timer_sta(123)
#endif
!$OMP MASTER
    call MPI_ALLTOALL(wk_send1, z2y_smax*kimg*iesize, mpi_double_precision,   &
   &                  wk_recv1, z2y_rmax*kimg*iesize, mpi_double_precision,   &
   &                                                  mpi_fft_yz_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Direct_XYZ_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 55, ierr)
     endif
!$OMP END MASTER
!$OMP BARRIER
#ifdef __TIMER_DO__
  call timer_sta(124)
#endif
#ifdef __TIMER_DO__
  call timer_end(124)
#endif

#else
!else ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_comm)
  call timer_sta(123)
#endif
!$OMP MASTER
    icnt_recv = 0
    lrank = myrank + 1
    if (lrank > (nmrank-1)) lrank = 0
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank-1)) lrank = 0
       if ((y2z_recv(1,lrank) /= 0)) then
          call mpi_irecv(wk_recv1(1,z2y_recv(2,lrank)), z2y_recv(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_r(icnt_recv), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Direct_XYZ_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 55, ierr)
           endif
          icnt_recv = icnt_recv + 1
       endif
    enddo

    icnt_send = 0
    lrank = myrank
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank-1)) lrank = 0
       if ((y2z_send(1,lrank) /= 0)) then
          call mpi_isend(wk_send1(1,z2y_send(2,lrank)), z2y_send(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_s(icnt_send), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Direct_XYZ_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 56, ierr)
           endif
          icnt_send = icnt_send + 1
       endif
    enddo

#ifdef __TIMER_DO__
  call timer_sta(124)
#endif
#ifdef __TIMER_DO__
  call timer_end(124)
#endif

    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Direct_XYZ_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 57, ierr)
     endif

    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Direct_XYZ_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 58, ierr)
     endif
!$OMP END MASTER
!$OMP BARRIER
#endif
!endif ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_end(123)
#endif

    nx = xyz_fft_y(2,1)-xyz_fft_y(1,1)+1
    ny = xyz_fft_y(2,2)-xyz_fft_y(1,2)+1
    nz = xyz_fft_y(2,3)-xyz_fft_y(1,3)+1
    nyp = ny
!
! Y-axis (z-x div)
!
#ifdef __TIMER_DO__
  call timer_sta(252)
#endif
    if(kimg==1) then
       nsize(1:3) = (/ny,nz,(nx/2)/)
       isin(1:3)  = (/-1,0,0/)
! (z,x,y) -> (y,z,x)
!$OMP DO
       do ib = 1, iesize
          do i = 1, nx/2
             do k = 1, nz
                do l = 1, fft_Z_y_dim
                   nis = nis_fft_Z_y(l)
                   nie = nie_fft_Z_y(l)
                   nny = nie - nis + 1
                   iadd0 = nx/2*nny*nz*(ib-1)*2
                   do j = nis, nie
                      iadd1 = iadd0 + (k+(i-1)*nz+(j-nis)*nz*nx/2)*2
                      iadd = (j+(k-1)*nyp+(i-1)*nyp*nz)*2
                      wk_afft_l(iadd-1,ib) = wk_recv1(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv1(iadd1  ,l)
                   end do
                end do
             end do
          end do
!x     end do
!xOMP DO
!x     do ib = 1, iesize
          call DVMCF2(wk_afft_l(1,ib),nsize,3,isin,ierr)
!x     end do
! (z,x,y) -> (y,z,x)
!xOMP DO
!x     do ib = 1, iesize
          do i = 1, nx/2
             do k = 1, nz
                do l = 1, fft_X_y_dim
                   nis = nis_fft_X_y(l)
                   nie = nie_fft_X_y(l)
                   nny = nie - nis + 1
                   iadd0 = nx/2*nny*nz*(ib-1)*2
                   do j = nis, nie
                      iadd  = (j+(k-1)*nyp+(i-1)*nyp*nz)*2
                      iadd1 = iadd0 + (j-nis+1+(k-1)*nny+(i-1)*nny*nz)*2
                      wk_send2(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send2(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
    else
     if (expo_y==1) then
        allocate(wwy((2*ny+70)*kimg))
     end if
     if (iesize > 1) then
       nsize(1:3) = (/ny,nz,nx/)
       isin(1:3)  = (/1,0,0/)
       isw = 1
! (y,z,x) -> (y,z,x)
!$OMP DO
       do ib = 1, iesize
          do i = 1, nx
             do k = 1, nz
                do l = 1, fft_Z_y_dim
                   nis = nis_fft_Z_y(l)
                   nie = nie_fft_Z_y(l)
                   nny = nie - nis + 1
                   iadd0 = nx*nny*nz*(ib-1)*2
                   do j = nis, nie
                      iadd1 = iadd0+(k+(i-1)*nz+(j-nis)*nz*nx)*2
                      iadd  = (j+(k-1)*nyp+(i-1)*nyp*nz)*2
                      wk_afft_l(iadd-1,ib) = wk_recv1(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv1(iadd1  ,l)
                   end do
                end do
             end do
          end do
!x     end do
!xOMP DO
!x     do ib = 1, iesize
          if (expo_y==1) then
             do i = 1, nx
                do k = 1, nz
                   call DVCFM1(wk_afft_l((1+(k-1)*ny+(i-1)*ny*nz)*2-1,ib),ny,isw,1,wwy,ierr)
!x                 call DVMCF2(wk_afft_l((1+(k-1)*ny+(i-1)*ny*nz)*2-1,ib),ny,1,1,ierr)
                end do
             end do
          else
             call DVMCF2(wk_afft_l(1,ib),nsize,3,isin,ierr)
          end if
!x     end do
! (y,z,x) -> (y,z,x)
!xOMP DO
!x     do ib = 1, iesize
          do i = 1, nx
             do k = 1, nz
                do l = 1, fft_X_y_dim
                   nis = nis_fft_X_y(l)
                   nie = nie_fft_X_y(l)
                   nny = nie - nis + 1
                   iadd0 = nx*nny*nz*(ib-1)*2
                   do j = nis, nie
                      iadd  = (j+(k-1)*nyp+(i-1)*nyp*nz)*2
                      iadd1 = iadd0+(j-nis+1+(k-1)*nny+(i-1)*nny*nz)*2
                      wk_send2(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send2(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
     else
       isw = 1
! (y,z,x) -> (y,z,x)
       do ib = 1, iesize
!$OMP DO
          do i = 1, nx
             do k = 1, nz
                do l = 1, fft_Z_y_dim
                   nis = nis_fft_Z_y(l)
                   nie = nie_fft_Z_y(l)
                   nny = nie - nis + 1
                   iadd0 = nx*nny*nz*(ib-1)*2
                   do j = nis, nie
                      iadd1 = iadd0+(k+(i-1)*nz+(j-nis)*nz*nx)*2
                      iadd  = (j+(k-1)*nyp+(i-1)*nyp*nz)*2
                      wk_afft_l(iadd-1,ib) = wk_recv1(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv1(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end do
       do ib = 1, iesize
          if (expo_y==1) then
!$OMP DO
             do i = 1, nx
                do k = 1, nz
                   call DVCFM1(wk_afft_l((1+(k-1)*ny+(i-1)*ny*nz)*2-1,ib),ny,isw,1,wwy,ierr)
                end do
             end do
          else
!$OMP DO
             do i = 1, nx
                do k = 1, nz
                   call DVMCF2(wk_afft_l((1+(k-1)*ny+(i-1)*ny*nz)*2-1,ib),ny,1,1,ierr)
                end do
             end do
          end if
       end do
! (y,z,x) -> (y,z,x)
       do ib = 1, iesize
!$OMP DO
          do i = 1, nx
             do k = 1, nz
                do l = 1, fft_X_y_dim
                   nis = nis_fft_X_y(l)
                   nie = nie_fft_X_y(l)
                   nny = nie - nis + 1
                   iadd0 = nx*nny*nz*(ib-1)*2
                   do j = nis, nie
                      iadd  = (j+(k-1)*nyp+(i-1)*nyp*nz)*2
                      iadd1 = iadd0+(j-nis+1+(k-1)*nny+(i-1)*nny*nz)*2
                      wk_send2(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send2(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
     end if
     if (expo_y==1) then
        deallocate(wwy)
     end if
    end if
#ifdef __TIMER_DO__
  call timer_end(252)
#endif

#ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_fft_xy_world)
  call timer_sta(127)
#endif
!$OMP MASTER
    call MPI_ALLTOALL(wk_send2, y2x_smax*kimg*iesize, mpi_double_precision,   &
   &                  wk_recv2, y2x_rmax*kimg*iesize, mpi_double_precision,   &
   &                                                  mpi_fft_xy_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Direct_XYZ_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 60, ierr)
     endif
!$OMP END MASTER
!$OMP BARRIER
#ifdef __TIMER_DO__
  call timer_sta(128)
#endif
#ifdef __TIMER_DO__
  call timer_end(128)
#endif

#else
!else ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_comm)
  call timer_sta(127)
#endif
!$OMP MASTER
    icnt_recv = 0
    lrank = myrank + 1
    if (lrank > (nmrank-1)) lrank = 0
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank-1)) lrank = 0
       if ((z2x_recv(1,lrank) /= 0)) then
          call mpi_irecv(wk_recv2(1,z2x_recv(2,lrank)), z2x_recv(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_r(icnt_recv), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Direct_XYZ_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 60, ierr)
           endif
          icnt_recv = icnt_recv + 1
       endif
    enddo

    icnt_send = 0
    lrank = myrank
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank-1)) lrank = 0
       if ((z2x_send(1,lrank) /= 0)) then
          call mpi_isend(wk_send2(1,z2x_send(2,lrank)), z2x_send(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_s(icnt_send), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Direct_XYZ_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 61, ierr)
           endif
          icnt_send = icnt_send + 1
       endif
    enddo

#ifdef __TIMER_DO__
  call timer_sta(128)
#endif
#ifdef __TIMER_DO__
  call timer_end(128)
#endif

    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Direct_XYZ_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 62, ierr)
     endif

    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Direct_XYZ_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 63, ierr)
     endif
!$OMP END MASTER
!$OMP BARRIER
#endif
!endif ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_end(127)
#endif

!   wk_afft_l = 0.0d0

    nx = xyz_fft_x(2,1)-xyz_fft_x(1,1)+1
    ny = xyz_fft_x(2,2)-xyz_fft_x(1,2)+1
    nz = xyz_fft_x(2,3)-xyz_fft_x(1,3)+1
!
! X-axis (y-z div)
!
#ifdef __TIMER_DO__
  call timer_sta(253)
#endif
    if(kimg==1) then
       nsize(1:3) = (/nx-2,ny,nz/)
       isin(1:3)  = (/-1,0,0/)
! (y,z,x) -> (x,y,z)
!$OMP DO
       do ib = 1, iesize
          do k = 1, nz
             do j = 1, ny
                do l = 1, fft_Y_x_dim
                   nis = nis_fft_Y_x(l)
                   nie = nie_fft_Y_x(l)
                   nnx = nie - nis + 1
                   if(mod(nis,2)>0) then
                     nisx = nis/2 + 1
                   else
                     nisx = nis/2
                   endif
                   niex = nie/2
                   iadd0 = nnx/2*ny*nz*(ib-1)*2
                   do i = nisx, niex
                      iadd1 = iadd0 + (j+(k-1)*ny+(i-nisx)*ny*nz)*2
                      iadd  = (i+(j-1)*nx/2+(k-1)*nx/2*ny)*2
                      wk_afft_l(iadd-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
          end do
!x     end do
!xOMP DO
!x     do ib = 1, iesize
          call DVMRF2(wk_afft_l(1,ib),nsize,3,isin,-1,ierr)
       end do
    else
     if (expo_x==1) then
        allocate(wwx((2*nx+70)*kimg))
     end if
     if (iesize > 1) then
       nsize(1:3) = (/nx,ny,nz/)
       isin(1:3)  = (/1,0,0/)
       isw = 1
! (y,z,x) -> (x,y,z)
!$OMP DO
       do ib = 1, iesize
          do k = 1, nz
             do j = 1, ny
                do l = 1, fft_Y_x_dim
                   nis = nis_fft_Y_x(l)
                   nie = nie_fft_Y_x(l)
                   nnx = nie - nis + 1
                   iadd0 = nnx*ny*nz*(ib-1)*2
                   do i = nis, nie
                      iadd1 = iadd0+(j+(k-1)*ny+(i-nis)*ny*nz)*2
                      iadd  = (i+(j-1)*nx+(k-1)*nx*ny)*2
                      wk_afft_l(iadd-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
          end do
!x     end do
!xOMP DO
!x     do ib = 1, iesize
          if (expo_x==1) then
             do k = 1, nz
                do j = 1, ny
                   call DVCFM1(wk_afft_l((1+(j-1)*nx+(k-1)*nx*ny)*2-1,ib),nx,isw,1,wwx,ierr)
!x                 call DVMCF2(wk_afft_l((1+(j-1)*nx+(k-1)*nx*ny)*2-1,ib),nx,1,1,ierr)
                end do
             end do
          else
             call DVMCF2(wk_afft_l(1,ib),nsize,3,isin,ierr)
          end if
       end do
     else
       isw = 1
! (y,z,x) -> (x,y,z)
       do ib = 1, iesize
!$OMP DO
          do k = 1, nz
             do j = 1, ny
                do l = 1, fft_Y_x_dim
                   nis = nis_fft_Y_x(l)
                   nie = nie_fft_Y_x(l)
                   nnx = nie - nis + 1
                   iadd0 = nnx*ny*nz*(ib-1)*2
                   do i = nis, nie
                      iadd1 = iadd0+(j+(k-1)*ny+(i-nis)*ny*nz)*2
                      iadd  = (i+(j-1)*nx+(k-1)*nx*ny)*2
                      wk_afft_l(iadd-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end do
       do ib = 1, iesize
          if (expo_x==1) then
!$OMP DO
             do k = 1, nz
                do j = 1, ny
                   call DVCFM1(wk_afft_l((1+(j-1)*nx+(k-1)*nx*ny)*2-1,ib),nx,isw,1,wwx,ierr)
                end do
             end do
          else
!$OMP DO
             do k = 1, nz
                do j = 1, ny
                   call DVMCF2(wk_afft_l((1+(j-1)*nx+(k-1)*nx*ny)*2-1,ib),nx,1,1,ierr)
                end do
             end do
          end if
       end do
     end if
     if (expo_x==1) then
        deallocate(wwx)
     end if
    end if
#ifdef __TIMER_DO__
  call timer_end(253)
#endif
!$OMP END PARALLEL

1000 continue

#ifdef __TIMER_SUB__
    call timer_end(105)
#endif
  end subroutine m_FFT_Direct_XYZ_3D
!------------------------------------------------------------------------------
#else
!else #ifdef FFT_USE_SSL2
#ifdef FFTW_NOSTRIDE
!------------------------------------------------------------------------------
  subroutine m_FFT_Direct_XYZ_3D(nfout,wk_afft_l, wk_size, iesize)
    integer, intent(in) :: nfout, wk_size, iesize
    real(kind=DP), dimension(wk_size*kimg,iesize), intent(inout) :: wk_afft_l

    integer, parameter :: FFTW_MEASURE=0
    integer, parameter :: FFTW_ESTIMATE=64
    integer, parameter :: FFTW_RANK=1
    integer            :: FFTW_FLAG
    integer,dimension(FFTW_RANK) :: NFFTW3, NEMBED, NEREAL
    integer :: nx, ny, nz, i, j, k, l, ri, ib, iadd, lrank
    integer :: irank, itag, icnt_send, icnt_recv
    integer,save, allocatable, dimension(:)   :: req_r, req_s
    integer,save, allocatable, dimension(:,:) :: sta_r, sta_s
    integer,save, allocatable, dimension(:)   :: wk_recvcnt, wk_sendcnt
    real(kind=DP),allocatable, dimension(:,:),save :: wk_recv1, wk_send1
    real(kind=DP),allocatable, dimension(:,:),save :: wk_recv2, wk_send2
    integer :: mpi_comm, myrank, nmrank
    integer(kind=8),save :: planx1 = 0, planx2 = 0
    integer(kind=8),save :: plany1 = 0, plany2 = 0
    integer(kind=8),save :: planz1 = 0, planz2 = 0
    integer(kind=8),save :: planx1_1d = 0, planx2_1d = 0
    integer(kind=8),save :: plany1_1d = 0, plany2_1d = 0
    integer(kind=8),save :: planz1_1d = 0, planz2_1d = 0
    integer(kind=8) :: plan
    integer :: max_x, max_y, max_z
    integer,save :: max_elm = 0, savesize = 0

    integer,save, allocatable, dimension(:,:) :: z2y_recv, z2y_send, y2x_recv, y2x_send
    integer,save :: z2y_rrank, z2y_srank, y2x_rrank, y2x_srank
    integer,save :: z2y_rmax, z2y_smax, z2y_srmax, y2x_rmax, y2x_smax, y2x_srmax

    integer ::  nnx, nny, nnz, iadd0, iadd1, nis, nie, jadd, kadd, nisx, niex, nxp, nyp, nzp

    integer,dimension(3) :: nsize, isin
#ifdef FFT_ALLTOALL
    integer,dimension(2) :: max_send,max_recv
#endif

#ifdef __TIMER_SUB__
    call mpi_barrier(mpi_ke_world,ierr)
    call timer_sta(105)
#endif

    if (kimg == 1) then
       FFTW_FLAG = +1      ! FFTW_BACKWARD
    else
       FFTW_FLAG = -1      ! FFTW_FORWARD
    endif
    plan = 0

!!  mpi_comm = mpi_kg_world
!!  myrank = myrank_e
!!  nmrank = nrank_e
    mpi_comm = mpi_ke_world
    myrank = myrank_g
    nmrank = nrank_g
    itag = 10

    if (np_fft_x==0)then
       nx = 0
       nx = 0
       nz = 0
       go to 1000
    endif

    if (firstcall_direct_xyz_3d) then
       savesize = 0
#ifdef __TIMER_DO__
  call timer_sta(118)
#endif
       max_x = maxval(nel_fft_x(:))
       max_y = maxval(nel_fft_y(:))
       max_z = maxval(nel_fft_z(:))
       max_elm = max(max_x,max_y,max_z)

       if(allocated(req_r)) deallocate(req_r)
       if(allocated(req_s)) deallocate(req_s)
       if(allocated(sta_r)) deallocate(sta_r)
       if(allocated(sta_s)) deallocate(sta_s)
       if(allocated(wk_recvcnt)) deallocate(wk_recvcnt)
       if(allocated(wk_sendcnt)) deallocate(wk_sendcnt)
       allocate(req_r(0:nmrank-1), stat=ierr)
       allocate(req_s(0:nmrank-1), stat=ierr)
       allocate(sta_r(MPI_STATUS_SIZE,0:nmrank-1), stat=ierr)
       allocate(sta_s(MPI_STATUS_SIZE,0:nmrank-1), stat=ierr)
       allocate(wk_recvcnt(0:nmrank-1), stat=ierr)
       allocate(wk_sendcnt(0:nmrank-1), stat=ierr)
!      allocate(wk_recvrank(0:nmrank-1), stat=ierr)
!      allocate(wk_sendrank(0:nmrank-1), stat=ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_FFT_Direct_XYZ_3D :  Not allocate '
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 51, ierr)
        endif

       if(allocated(z2y_recv)) deallocate(z2y_recv)
       if(allocated(z2y_send)) deallocate(z2y_send)
       allocate(z2y_recv(2,0:nmrank-1))
       allocate(z2y_send(2,0:nmrank-1))
       z2y_send = 0
       z2y_recv = 0
       wk_recvcnt = 0
       wk_sendcnt = 0
       do i = 1, nel_fft_y(myrank)
          irank = map_fft_z(mp_fft_y(i)) - 1
          wk_recvcnt(irank) = wk_recvcnt(irank) + 1
       enddo
       do i = 1, nel_fft_z(myrank)
          irank = map_fft_y(mp_fft_z(i)) - 1
          wk_sendcnt(irank) = wk_sendcnt(irank) + 1
       enddo
       k = 0
       do i = 0, nmrank - 1
          if(wk_recvcnt(i) /= 0) then
             k = k + 1
             z2y_recv(1,i) = wk_recvcnt(i)
             z2y_recv(2,i) = k
          endif
       enddo
       z2y_rrank = k
       k = 0
       do i = 0, nmrank - 1
          if(wk_sendcnt(i) /= 0) then
             k = k + 1
             z2y_send(1,i) = wk_sendcnt(i)
             z2y_send(2,i) = k
          endif
       enddo
       z2y_srank = k
       z2y_rmax = maxval(wk_recvcnt)
       z2y_smax = maxval(wk_sendcnt)
       z2y_srmax = max(z2y_rmax,z2y_smax)

       if(allocated(y2x_recv)) deallocate(y2x_recv)
       if(allocated(y2x_send)) deallocate(y2x_send)
       allocate(y2x_recv(2,0:nmrank-1))
       allocate(y2x_send(2,0:nmrank-1))
       y2x_send = 0
       y2x_recv = 0
       wk_recvcnt = 0
       wk_sendcnt = 0
       do i = 1, nel_fft_x(myrank)
          irank = map_fft_y(mp_fft_x(i)) - 1
          wk_recvcnt(irank) = wk_recvcnt(irank) + 1
       enddo
       do i = 1, nel_fft_y(myrank)
          irank = map_fft_x(mp_fft_y(i)) - 1
          wk_sendcnt(irank) = wk_sendcnt(irank) + 1
       enddo
       k = 0
       do i = 0, nmrank - 1
          if(wk_recvcnt(i) /= 0) then
             k = k + 1
             y2x_recv(1,i) = wk_recvcnt(i)
             y2x_recv(2,i) = k
          endif
       enddo
       y2x_rrank = k
       k = 0
      do i = 0, nmrank - 1
          if(wk_sendcnt(i) /= 0) then
             k = k + 1
             y2x_send(1,i) = wk_sendcnt(i)
             y2x_send(2,i) = k
          endif
       enddo
       y2x_srank = k
       y2x_rmax = maxval(wk_recvcnt)
       y2x_smax = maxval(wk_sendcnt)
       y2x_srmax = max(y2x_rmax,y2x_smax)

#ifdef FFT_ALLTOALL
       max_send(1) = z2y_rmax
       max_send(2) = z2y_smax
       call mpi_allreduce(max_send,max_recv,2,mpi_integer,mpi_max,mpi_fft_yz_world,ierr)
       z2y_rmax = max_recv(1)
       z2y_smax = max_recv(2)
       max_send(1) = y2x_rmax
       max_send(2) = y2x_smax
       call mpi_allreduce(max_send,max_recv,2,mpi_integer,mpi_max,mpi_fft_xy_world,ierr)
       y2x_rmax = max_recv(1)
       y2x_smax = max_recv(2)
#endif
       if (ipri > 1) then
          write(nfout,'("m_FFT_Direct_XYZ_3D   --   myrank_g=",i4)') myrank_g
          write(nfout,'("z2y_send")')
          write(nfout,'(10(i8,", "))') (z2y_send(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (z2y_send(2,i),i=0,nmrank-1)
          write(nfout,'("z2y_recv")')
          write(nfout,'(10(i8,", "))') (z2y_recv(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (z2y_recv(2,i),i=0,nmrank-1)
          write(nfout,'("y2x_send")')
          write(nfout,'(10(i8,", "))') (y2x_send(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (y2x_send(2,i),i=0,nmrank-1)
          write(nfout,'("y2x_recv")')
          write(nfout,'(10(i8,", "))') (y2x_recv(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (y2x_recv(2,i),i=0,nmrank-1)
       endif

#ifdef __TIMER_DO__
  call timer_sta(119)
#endif
       nx = xyz_fft_z(2,1)-xyz_fft_z(1,1)+1
       ny = xyz_fft_z(2,2)-xyz_fft_z(1,2)+1
       nz = xyz_fft_z(2,3)-xyz_fft_z(1,3)+1
!      if(kimg==1) then
          NFFTW3(1) = nz
          NEMBED(1) = nz
          call dfftw_plan_many_dft    (planz1,     1, nz, nx*ny/2,  &
         &                             wk_afft_l, nz,  1,      nz,  &
         &                             wk_afft_l, nz,  1,      nz,  &
         &                             FFTW_FLAG, FFTW_MEASURE  )
!      else
          NFFTW3(1) = nz
          NEMBED(1) = nz
          call dfftw_plan_many_dft    (planz2,     1, nz, nx*ny,  &
         &                             wk_afft_l, nz,  1,    nz,  &
         &                             wk_afft_l, nz,  1,    nz,  &
         &                             FFTW_FLAG, FFTW_MEASURE )
          call dfftw_plan_dft_1d(planz2_1d, nz, wk_afft_l, wk_afft_l, FFTW_FLAG, FFTW_MEASURE)
!      endif

       nx = xyz_fft_y(2,1)-xyz_fft_y(1,1)+1
       ny = xyz_fft_y(2,2)-xyz_fft_y(1,2)+1
       nz = xyz_fft_y(2,3)-xyz_fft_y(1,3)+1
!      if(kimg==1) then
          call dfftw_plan_many_dft    (plany1   ,  1, ny, nz*nx/2,  &
         &                             wk_afft_l, ny,  1,      ny,  &
         &                             wk_afft_l, ny,  1,      ny,  &
         &                             FFTW_FLAG, FFTW_MEASURE )
!      else
          call dfftw_plan_many_dft    (plany2   ,  1, ny, nz*nx,   &
         &                             wk_afft_l, ny,  1,    ny,   &
         &                             wk_afft_l, ny,  1,    ny,   &
         &                             FFTW_FLAG, FFTW_MEASURE )
          call dfftw_plan_dft_1d(plany2_1d, ny, wk_afft_l, wk_afft_l, FFTW_FLAG, FFTW_MEASURE)
!      end if

       nx = xyz_fft_x(2,1)-xyz_fft_x(1,1)+1
       ny = xyz_fft_x(2,2)-xyz_fft_x(1,2)+1
       nz = xyz_fft_x(2,3)-xyz_fft_x(1,3)+1
!      if(kimg==1) then
          call dfftw_plan_many_dft_c2r(planx1   ,   1, nx-2, ny*nz,   &
         &                             wk_afft_l, nx,     1,  nx/2,   &
         &                             wk_afft_l, nx,     1,    nx,   &
         &                             FFTW_MEASURE )
!      else
          call dfftw_plan_many_dft    (planx2   ,  1, nx, ny*nz,   &
         &                             wk_afft_l, nx,  1,    nx,   &
         &                             wk_afft_l, nx,  1,    nx,   &
         &                             FFTW_FLAG, FFTW_MEASURE )
          call dfftw_plan_dft_1d(planx2_1d, nx, wk_afft_l, wk_afft_l, FFTW_FLAG, FFTW_MEASURE)
!      endif
#ifdef __TIMER_DO__
  call timer_end(119)
#endif

       firstcall_direct_xyz_3d = .false.
       go to 1000
#ifdef __TIMER_DO__
  call timer_end(118)
#endif
    endif

#ifdef FFT_ALLTOALL
    if (iesize /= savesize) then
#else
    if (iesize > savesize) then
#endif
#ifdef __TIMER_DO__
  call timer_sta(120)
#endif
       if (allocated(wk_recv1)) deallocate(wk_recv1)
       if (allocated(wk_send1)) deallocate(wk_send1)
       if (allocated(wk_recv2)) deallocate(wk_recv2)
       if (allocated(wk_send2)) deallocate(wk_send2)
       allocate(wk_recv1(z2y_rmax*kimg*iesize,z2y_rrank), stat=ierr)
       allocate(wk_send1(z2y_smax*kimg*iesize,z2y_srank), stat=ierr)
       allocate(wk_recv2(y2x_rmax*kimg*iesize,y2x_rrank), stat=ierr)
       allocate(wk_send2(y2x_smax*kimg*iesize,y2x_srank), stat=ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_FFT_Direct_XYZ_3D :  Not allocate '
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 52, ierr)
        endif
       savesize = iesize
#ifdef __TIMER_DO__
  call timer_end(120)
#endif
    end if

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP          SHARED(     nfout,wk_afft_l,iesize,req_r,req_s,sta_r,sta_s,    &
!$OMP                 wk_recv1,wk_send1,wk_recv2,wk_send2,                    &
!$OMP                 planx1,planx2,planz1,planz2,plany1,plany2,              &
!$OMP                 planx1_1d,planx2_1d,planz1_1d,planz2_1d,plany1_1d,plany2_1d, &
!$OMP                 z2y_recv, z2y_send, y2x_recv, y2x_send,itag,            &
!$OMP                 xyz_fft_y,xyz_fft_z,xyz_fft_x,fft_Z_x_dim,nis_fft_Z_x,  &
!$OMP                 nie_fft_Z_x,fft_X_z_dim,nis_fft_X_z,nie_fft_X_z,        &
!$OMP                 fft_Y_z_dim,nis_fft_Y_z,nie_fft_Y_z,fft_Z_y_dim,        &
!$OMP                 nie_fft_Z_y,nis_fft_Z_y,                                &
!$OMP                 nie_fft_X_y,nis_fft_X_y,fft_X_y_dim,                    &
!$OMP                 nie_fft_Y_x,nis_fft_Y_x,fft_Y_x_dim,                    &
#ifdef FFT_ALLTOALL
!$OMP                 mpi_fft_yz_world, z2y_rmax, z2y_smax,                   &
!$OMP                 mpi_fft_xy_world, y2x_rmax, y2x_smax,                   &
#endif
!$OMP                 mpi_comm,myrank,nmrank,kimg,ierr                    )   &
!$OMP          PRIVATE(nx,ny,nz,i,j,k,l,ri,ib,iadd,lrank,icnt_send,icnt_recv, &
!$OMP                  nxp,nyp,nzp,nisx,niex,                                  &
!$OMP                  irank,nnx,nny,nnz,iadd0,iadd1,nis,nie,jadd,kadd,plan   )
    nx = xyz_fft_z(2,1)-xyz_fft_z(1,1)+1
    ny = xyz_fft_z(2,2)-xyz_fft_z(1,2)+1
    nz = xyz_fft_z(2,3)-xyz_fft_z(1,3)+1
    nzp = nz
#ifdef __TIMER_DO__
  call timer_sta(251)
#endif
    if (kimg == 1) then
!
! Z-axis (x-y div)
!
!$OMP DO
       do ib = 1, iesize
          call dfftw_execute_dft(planz1,wk_afft_l(1,ib),wk_afft_l(1,ib))
!x     end do
! (z,x,y) -> (z,x,y)
!xOMP DO
!x     do ib = 1, iesize
          do j = 1, ny
             do i = 1, nx/2
                do l = 1, fft_Y_z_dim
                   nis = nis_fft_Y_z(l)
                   nie = nie_fft_Y_z(l)
                   nnz = nie - nis + 1
                   iadd0 = nx/2*ny*nnz*(ib-1)*2
                   do k = nis, nie
                      iadd  = (k+(i-1)*nzp+(j-1)*nzp*nx/2)*2
                      iadd1 = iadd0+((k-nis+1)+(i-1)*nnz+(j-1)*nnz*nx/2)*2
                      wk_send1(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send1(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
    else
!
! Z-axis (x-y div)
!
     if (iesize > 1) then
!$OMP DO
       do ib = 1, iesize
          call dfftw_execute_dft(planz2,wk_afft_l(1,ib),wk_afft_l(1,ib))
!x     end do
! (z,x,y) -> (z,x,y)
!xOMP DO
!x     do ib = 1, iesize
          do j = 1, ny
             do i = 1, nx
                do l = 1, fft_Y_z_dim
                   nis = nis_fft_Y_z(l)
                   nie = nie_fft_Y_z(l)
                   nnz = nie - nis + 1
                   iadd0 = nx*ny*nnz*(ib-1)*2
                   do k = nis, nie
                      iadd  = (k+(i-1)*nzp+(j-1)*nzp*nx)*2
                      iadd1 = iadd0+((k-nis+1)+(i-1)*nnz+(j-1)*nnz*nx)*2
                      wk_send1(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send1(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
     else
       do ib = 1, iesize
!$OMP DO
          do j = 1, ny
             do i = 1, nx
                iadd = (1+(i-1)*nzp+(j-1)*nzp*nx)*2-1
                call dfftw_execute_dft(planz2_1d,wk_afft_l(iadd,ib),wk_afft_l(iadd,ib))
             end do
!x        end do
!x     end do
! (z,x,y) -> (z,x,y)
!x     do ib = 1, iesize
!x        do j = 1, ny
             do i = 1, nx
                do l = 1, fft_Y_z_dim
                   nis = nis_fft_Y_z(l)
                   nie = nie_fft_Y_z(l)
                   nnz = nie - nis + 1
                   iadd0 = nx*ny*nnz*(ib-1)*2
                   do k = nis, nie
                      iadd  = (k+(i-1)*nzp+(j-1)*nzp*nx)*2
                      iadd1 = iadd0+((k-nis+1)+(i-1)*nnz+(j-1)*nnz*nx)*2
                      wk_send1(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send1(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
     end if
    end if
#ifdef __TIMER_DO__
  call timer_end(251)
#endif

#ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_fft_yz_world)
  call timer_sta(123)
#endif

!$OMP MASTER
    call MPI_ALLTOALL(wk_send1, z2y_smax*kimg*iesize, mpi_double_precision,   &
   &                  wk_recv1, z2y_rmax*kimg*iesize, mpi_double_precision,   &
   &                                                  mpi_fft_yz_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Direct_XYZ_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 55, ierr)
     endif
!$OMP END MASTER
!$OMP BARRIER

#ifdef __TIMER_DO__
  call timer_sta(124)
#endif
#ifdef __TIMER_DO__
  call timer_end(124)
#endif

#else
!else ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_comm)
  call timer_sta(123)
#endif

!$OMP MASTER
    icnt_recv = 0
    lrank = myrank + 1
    if (lrank > (nmrank-1)) lrank = 0
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank-1)) lrank = 0
       if ((y2z_recv(1,lrank) /= 0)) then
          call mpi_irecv(wk_recv1(1,z2y_recv(2,lrank)), z2y_recv(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_r(icnt_recv), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Direct_XYZ_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 55, ierr)
           endif
          icnt_recv = icnt_recv + 1
       endif
    enddo

    icnt_send = 0
    lrank = myrank
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank-1)) lrank = 0
       if ((y2z_send(1,lrank) /= 0)) then
          call mpi_isend(wk_send1(1,z2y_send(2,lrank)), z2y_send(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_s(icnt_send), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Direct_XYZ_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 56, ierr)
           endif
          icnt_send = icnt_send + 1
       endif
    enddo

#ifdef __TIMER_DO__
  call timer_sta(124)
#endif
#ifdef __TIMER_DO__
  call timer_end(124)
#endif

    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Direct_XYZ_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 57, ierr)
     endif

    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Direct_XYZ_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 58, ierr)
     endif
!$OMP END MASTER
!$OMP BARRIER

#endif
!endif ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_end(123)
#endif

    nx = xyz_fft_y(2,1)-xyz_fft_y(1,1)+1
    ny = xyz_fft_y(2,2)-xyz_fft_y(1,2)+1
    nz = xyz_fft_y(2,3)-xyz_fft_y(1,3)+1
    nyp = ny
#ifdef __TIMER_DO__
  call timer_sta(252)
#endif
    if (kimg == 1) then
! (z,x,y) -> (y,z,x)
!$OMP DO
       do ib = 1, iesize
          do i = 1, nx/2
             do k = 1, nz
                do l = 1, fft_Z_y_dim
                   nis = nis_fft_Z_y(l)
                   nie = nie_fft_Z_y(l)
                   nny = nie - nis + 1
                   iadd0 = nx/2*nny*nz*(ib-1)*2
                   do j = nis, nie
                      iadd1 = iadd0 + (k+(i-1)*nz+(j-nis)*nz*nx/2)*2
                      iadd = (j+(k-1)*nyp+(i-1)*nyp*nz)*2
                      wk_afft_l(iadd-1,ib) = wk_recv1(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv1(iadd1  ,l)
                   end do
                end do
             end do
          end do
!x     end do
!
! Y-axis (z-x div)
!
!xOMP DO
!x     do ib = 1, iesize
          call dfftw_execute_dft(plany1,wk_afft_l(1,ib),wk_afft_l(1,ib))
!x     end do
! (y,z,x) -> (y,z,x)
!xOMP DO
!x     do ib = 1, iesize
          do i = 1, nx/2
             do k = 1, nz
                do l = 1, fft_X_y_dim
                   nis = nis_fft_X_y(l)
                   nie = nie_fft_X_y(l)
                   nny = nie - nis + 1
                   iadd0 = nx/2*nny*nz*(ib-1)*2
                   do j = nis, nie
                      iadd  = (j+(k-1)*nyp+(i-1)*nyp*nz)*2
                      iadd1 = iadd0 + (j-nis+1+(k-1)*nny+(i-1)*nny*nz)*2
                      wk_send2(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send2(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
    else
     if (iesize > 1) then
! (z,x,y) -> (y,z,x)
!$OMP DO
       do ib = 1, iesize
          do i = 1, nx
             do k = 1, nz
                do l = 1, fft_Z_y_dim
                   nis = nis_fft_Z_y(l)
                   nie = nie_fft_Z_y(l)
                   nny = nie - nis + 1
                   iadd0 = nx*nny*nz*(ib-1)*2
                   do j = nis, nie
                      iadd1 = iadd0+(k+(i-1)*nz+(j-nis)*nz*nx)*2
                      iadd  = (j+(k-1)*nyp+(i-1)*nyp*nz)*2
                      wk_afft_l(iadd-1,ib) = wk_recv1(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv1(iadd1  ,l)
                   end do
                end do
             end do
          end do
!x     end do
!
! Y-axis (z-x div)
!
!xOMP DO
!x     do ib = 1, iesize
          call dfftw_execute_dft(plany2,wk_afft_l(1,ib),wk_afft_l(1,ib))
!x     end do
! (y,z,x) -> (y,z,x)
!xOMP DO
!x       do ib = 1, iesize
          do i = 1, nx
             do k = 1, nz
                do l = 1, fft_X_y_dim
                   nis = nis_fft_X_y(l)
                   nie = nie_fft_X_y(l)
                   nny = nie - nis + 1
                   iadd0 = nx*nny*nz*(ib-1)*2
                   do j = nis, nie
                      iadd  = (j+(k-1)*nyp+(i-1)*nyp*nz)*2
                      iadd1 = iadd0+(j-nis+1+(k-1)*nny+(i-1)*nny*nz)*2
                      wk_send2(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send2(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
     else
! (z,x,y) -> (y,z,x)
       do ib = 1, iesize
!$OMP DO
          do i = 1, nx
             do k = 1, nz
                do l = 1, fft_Z_y_dim
                   nis = nis_fft_Z_y(l)
                   nie = nie_fft_Z_y(l)
                   nny = nie - nis + 1
                   iadd0 = nx*nny*nz*(ib-1)*2
                   do j = nis, nie
                      iadd1 = iadd0+(k+(i-1)*nz+(j-nis)*nz*nx)*2
                      iadd  = (j+(k-1)*nyp+(i-1)*nyp*nz)*2
                      wk_afft_l(iadd-1,ib) = wk_recv1(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv1(iadd1  ,l)
                   end do
                end do
             end do
!o        end do
!o     end do
!
! Y-axis (z-x div)
!
!o     do ib = 1, iesize
!o        do i = 1, nx
             do k = 1, nz
              iadd = (1+(k-1)*nyp+(i-1)*nyp*nz)*2-1
              call dfftw_execute_dft(plany2_1d,wk_afft_l(iadd,ib),wk_afft_l(iadd,ib))
             end do
!o        end do
!o     end do
! (y,z,x) -> (y,z,x)
!o     do ib = 1, iesize
!o        do i = 1, nx
             do k = 1, nz
                do l = 1, fft_X_y_dim
                   nis = nis_fft_X_y(l)
                   nie = nie_fft_X_y(l)
                   nny = nie - nis + 1
                   iadd0 = nx*nny*nz*(ib-1)*2
                   do j = nis, nie
                      iadd  = (j+(k-1)*nyp+(i-1)*nyp*nz)*2
                      iadd1 = iadd0+(j-nis+1+(k-1)*nny+(i-1)*nny*nz)*2
                      wk_send2(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send2(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
     end if
    end if
#ifdef __TIMER_DO__
  call timer_end(252)
#endif

#ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_fft_xy_world)
  call timer_sta(127)
#endif

!$OMP MASTER
    call MPI_ALLTOALL(wk_send2, y2x_smax*kimg*iesize, mpi_double_precision,   &
   &                  wk_recv2, y2x_rmax*kimg*iesize, mpi_double_precision,   &
   &                                                  mpi_fft_xy_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Direct_XYZ_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 60, ierr)
     endif
!$OMP END MASTER
!$OMP BARRIER

#ifdef __TIMER_DO__
  call timer_sta(128)
#endif
#ifdef __TIMER_DO__
  call timer_end(128)
#endif

#else
!else ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_comm)
  call timer_sta(127)
#endif

!$OMP MASTER
    icnt_recv = 0
    lrank = myrank + 1
    if (lrank > (nmrank-1)) lrank = 0
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank-1)) lrank = 0
       if ((z2x_recv(1,lrank) /= 0)) then
          call mpi_irecv(wk_recv2(1,z2x_recv(2,lrank)), z2x_recv(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_r(icnt_recv), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Direct_XYZ_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 60, ierr)
           endif
          icnt_recv = icnt_recv + 1
       endif
    enddo

    icnt_send = 0
    lrank = myrank
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank-1)) lrank = 0
       if ((z2x_send(1,lrank) /= 0)) then
          call mpi_isend(wk_send2(1,z2x_send(2,lrank)), z2x_send(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_s(icnt_send), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Direct_XYZ_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 61, ierr)
           endif
          icnt_send = icnt_send + 1
       endif
    enddo

#ifdef __TIMER_DO__
  call timer_sta(128)
#endif
#ifdef __TIMER_DO__
  call timer_end(128)
#endif

    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Direct_XYZ_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 62, ierr)
     endif

    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Direct_XYZ_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 63, ierr)
     endif
!$OMP END MASTER
!$OMP BARRIER
#endif
!endif ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_end(127)
#endif

    nx = xyz_fft_x(2,1)-xyz_fft_x(1,1)+1
    ny = xyz_fft_x(2,2)-xyz_fft_x(1,2)+1
    nz = xyz_fft_x(2,3)-xyz_fft_x(1,3)+1
    nxp = nx
#ifdef __TIMER_DO__
  call timer_sta(253)
#endif
    if (kimg == 1) then
! (y,z,x) -> (x,y,z)
!$OMP DO
       do ib = 1, iesize
          do k = 1, nz
             do j = 1, ny
                do l = 1, fft_Y_x_dim
                   nis = nis_fft_Y_x(l)
                   nie = nie_fft_Y_x(l)
                   nnx = nie - nis + 1
                   if(mod(nis,2)>0) then
                     nisx = nis/2 + 1
                   else
                     nisx = nis/2
                   endif
                   niex = nie/2
                   iadd0 = nnx/2*ny*nz*(ib-1)*2
                   do i = nisx, niex
                      iadd1 = iadd0 + (j+(k-1)*ny+(i-nisx)*ny*nz)*2
                      iadd  = (i+(j-1)*nx/2+(k-1)*nx/2*ny)*2
                      wk_afft_l(iadd-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
          end do
!x     end do
!
! X-axis (y-z div)
!
!xOMP DO
!x     do ib = 1, iesize
          call dfftw_execute_dft_c2r(planx1,wk_afft_l(1,ib),wk_afft_l(1,ib))
       end do
    else
     if (iesize > 1) then
! (y,z,x) -> (x,y,z)
!$OMP DO
       do ib = 1, iesize
          do k = 1, nz
             do j = 1, ny
                do l = 1, fft_Y_x_dim
                   nis = nis_fft_Y_x(l)
                   nie = nie_fft_Y_x(l)
                   nnx = nie - nis + 1
                   iadd0 = nnx*ny*nz*(ib-1)*2
                   do i = nis, nie
                      iadd1 = iadd0+(j+(k-1)*ny+(i-nis)*ny*nz)*2
                      iadd  = (i+(j-1)*nxp+(k-1)*nxp*ny)*2
                      wk_afft_l(iadd-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
          end do
!x     end do
!
! X-axis (y-z div)
!
!xOMP DO
!x     do ib = 1, iesize
          call dfftw_execute_dft    (planx2,wk_afft_l(1,ib),wk_afft_l(1,ib))
       end do
     else
! (y,z,x) -> (x,y,z)
       do ib = 1, iesize
!$OMP DO
          do k = 1, nz
             do j = 1, ny
                do l = 1, fft_Y_x_dim
                   nis = nis_fft_Y_x(l)
                   nie = nie_fft_Y_x(l)
                   nnx = nie - nis + 1
                   iadd0 = nnx*ny*nz*(ib-1)*2
                   do i = nis, nie
                      iadd1 = iadd0+(j+(k-1)*ny+(i-nis)*ny*nz)*2
                      iadd  = (i+(j-1)*nxp+(k-1)*nxp*ny)*2
                      wk_afft_l(iadd-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
!o        end do
!o     end do
!
! X-axis (y-z div)
!
!o     do ib = 1, iesize
!o        do k = 1, nz
             do j = 1, ny
                iadd = (1+(j-1)*nxp+(k-1)*nxp*ny)*2-1
                call dfftw_execute_dft(planx2_1d,wk_afft_l(iadd,ib),wk_afft_l(iadd,ib))
             end do
          end do
       end do
     end if
    end if
#ifdef __TIMER_DO__
  call timer_end(253)
#endif

!$OMP END PARALLEL

1000 continue

#ifdef __TIMER_SUB__
    call timer_end(105)
#endif
  end subroutine m_FFT_Direct_XYZ_3D
!------------------------------------------------------------------------------
#else
!else #ifdef FFTW_NOSTRIDE
!------------------------------------------------------------------------------
  subroutine m_FFT_Direct_XYZ_3D(nfout,wk_afft_l, wk_size, iesize)
#ifdef KMATH_FFT3D
    use m_Const_Parameters, only : PAI2
    use kmath_fft3d_mod, only : KMATH_FFT3D_Transform
#endif
    integer, intent(in) :: nfout, wk_size, iesize
    real(kind=DP), dimension(wk_size*kimg,iesize), intent(inout) :: wk_afft_l

    integer, parameter :: FFTW_MEASURE=0
    integer, parameter :: FFTW_ESTIMATE=64
    integer, parameter :: FFTW_RANK=1
    integer            :: FFTW_FLAG
    integer,dimension(FFTW_RANK) :: NFFTW3, NEMBED, NEREAL
    integer :: nx, ny, nz, i, j, k, l, ri, ib, iadd, lrank
    integer :: irank, itag, icnt_send, icnt_recv
    integer,save, allocatable, dimension(:)   :: req_r, req_s
    integer,save, allocatable, dimension(:,:) :: sta_r, sta_s
    integer,save, allocatable, dimension(:)   :: wk_recvcnt, wk_sendcnt
    real(kind=DP),allocatable, dimension(:,:),save :: wk_recv1, wk_send1
    real(kind=DP),allocatable, dimension(:,:),save :: wk_recv2, wk_send2
    integer :: mpi_comm, myrank, nmrank
    integer(kind=8),save :: planx1 = 0, planx2 = 0
    integer(kind=8),save :: plany1 = 0, plany2 = 0
    integer(kind=8),save :: planz1 = 0, planz2 = 0
    integer(kind=8),save :: planx1_1d = 0, planx2_1d = 0
    integer(kind=8),save :: plany1_1d = 0, plany2_1d = 0
    integer(kind=8),save :: planz1_1d = 0, planz2_1d = 0
    integer(kind=8) :: plan
    integer :: max_x, max_y, max_z
    integer,save :: max_elm = 0, savesize = 0

    integer,save, allocatable, dimension(:,:) :: z2y_recv, z2y_send, y2x_recv, y2x_send
    integer,save :: z2y_rrank, z2y_srank, y2x_rrank, y2x_srank
    integer,save :: z2y_rmax, z2y_smax, z2y_srmax, y2x_rmax, y2x_smax, y2x_srmax

    integer ::  nnx, nny, nnz, iadd0, iadd1, nis, nie, jadd, kadd, nisx, niex, nxp, nyp, nzp

    integer,dimension(3) :: nsize, isin
#ifdef FFT_ALLTOALL
    integer,dimension(2) :: max_send,max_recv
#endif
#ifdef KMATH_FFT3D
    complex(kind=DP), allocatable, dimension(:) :: wk_afft_l_in, wk_afft_l_out
    integer :: icount
    integer, dimension(3) :: nsize2
    integer :: nx1,ny1,nz1,nx2,ny2,nz2
#endif
    integer :: id_sname = -1

    call tstatc0_begin('m_FFT_Direct ',id_sname)
#ifdef __TIMER_SUB__
    call mpi_barrier(mpi_ke_world,ierr)
    call timer_sta(105)
#endif

#ifdef KMATH_FFT3D
    if(sw_kmath_fft3d==ON .and. .not.firstcall_direct_xyz_3d) then
      nsize2 = fft_box_size_WF(:,1)

      nx1 = nsize2(1)/nproc_fft3d_direct(1)
      ny1 = nsize2(2)/nproc_fft3d_direct(2)
      nz1 = nsize2(3)/nproc_fft3d_direct(3)
      if (MOD(nsize2(1),nproc_fft3d_direct(1)) /= 0) &
        nx1 = nx1 + 1
      if (MOD(nsize2(2),nproc_fft3d_direct(2)) /= 0) &
        ny1 = ny1 + 1
      if (MOD(nsize2(3),nproc_fft3d_direct(3)) /= 0) &
        nz1 = nz1 + 1

      nx2 = nsize2(1)/nproc_fft3d_direct(2)
      if (MOD(nsize2(1),nproc_fft3d_direct(2)) /= 0) &
        nx2 = nx2 + 1
      ny2 = nsize2(2)/nproc_fft3d_direct(1)
      if (MOD(nsize2(2),nproc_fft3d_direct(1)) /= 0) &
        ny2 = ny2 + 1
      if (MOD(ny2,nproc_fft3d_direct(3)) /= 0) then
        ny2 = ny2/nproc_fft3d_direct(3) + 1
      else
        ny2 = ny2/nproc_fft3d_direct(3)
      end if
      nz2 = nsize2(3)

      nx = xyz_fft_z(2,1)-xyz_fft_z(1,1)+1
      ny = xyz_fft_z(2,2)-xyz_fft_z(1,2)+1
      nz = xyz_fft_z(2,3)-xyz_fft_z(1,3)+1
      allocate(wk_afft_l_in (wk_size)) ; wk_afft_l_in  = cmplx(0.d0,0.d0)
      allocate(wk_afft_l_out(wk_size)) ; wk_afft_l_out = cmplx(0.d0,0.d0)
      do k=1, nz
        do j=1, ny
          do l=1, nx
             iadd0 = k + nz*(l-1) + nz*nx*(j-1)
             iadd1 = l + nx*(j-1) + nx*ny*(k-1)
             wk_afft_l_in(iadd0) = cmplx(wk_afft_l(2*iadd1-1,1),-wk_afft_l(2*iadd1,1))
          enddo
        enddo
      enddo
      call KMATH_FFT3D_Transform(kmath3d_handle_wf_direct, wk_afft_l_in, wk_afft_l_out,.true.)
      nx = xyz_fft_x(2,1)-xyz_fft_x(1,1)+1
      ny = xyz_fft_x(2,2)-xyz_fft_x(1,2)+1
      nz = xyz_fft_x(2,3)-xyz_fft_x(1,3)+1
      !do k=1, nz
      !  do j=1, ny
      !    do l=1, nx
      !       iadd1 = l + nx2*(j-1) + nx2*ny2*(k-1)
      !       wk_afft_l_in(iadd1) = wk_afft_l_out(iadd1) * wk_size
      !    enddo
      !  enddo
      !enddo
      do i=1,wk_size
        wk_afft_l_in(i) = wk_afft_l_out(i)*product(fft_box_size_WF(:,1))
      enddo

      do i=1,wk_size
        wk_afft_l(2*i-1,1) =  real (wk_afft_l_in(i))
        wk_afft_l(2*i,1)   = -aimag(wk_afft_l_in(i))
      enddo
      deallocate(wk_afft_l_in)
      deallocate(wk_afft_l_out)
      call tstatc0_end(id_sname)
      return
    endif
#endif
    if (kimg == 1) then
       FFTW_FLAG = +1      ! FFTW_BACKWARD
    else
       FFTW_FLAG = -1      ! FFTW_FORWARD
    endif
    plan = 0

!!  mpi_comm = mpi_kg_world
!!  myrank = myrank_e
!!  nmrank = nrank_e
    mpi_comm = mpi_ke_world
    myrank = myrank_g
    nmrank = nrank_g
    itag = 10

    if (np_fft_x==0)then
       nx = 0
       nx = 0
       nz = 0
       go to 1000
    endif

    if (firstcall_direct_xyz_3d) then
#ifdef __FAPP__
    call fapp_start('fftwf_direct_firstcall',1,1)
#endif
       if(planz1    /= 0) call dfftw_destroy_plan(planz1)
       if(planz1_1d /= 0) call dfftw_destroy_plan(planz1_1d)
       if(planz2    /= 0) call dfftw_destroy_plan(planz2)
       if(planz2_1d /= 0) call dfftw_destroy_plan(planz2_1d)
       if(plany1    /= 0) call dfftw_destroy_plan(plany1)
       if(plany2    /= 0) call dfftw_destroy_plan(plany2)
       if(planx1    /= 0) call dfftw_destroy_plan(planx1)
       if(planx1_1d /= 0) call dfftw_destroy_plan(planx1_1d)
       if(planx2    /= 0) call dfftw_destroy_plan(planx2)
       if(planx2_1d /= 0) call dfftw_destroy_plan(planx2_1d)
       savesize = 0
#ifdef __TIMER_DO__
  call timer_sta(118)
#endif
       max_x = maxval(nel_fft_x(:))
       max_y = maxval(nel_fft_y(:))
       max_z = maxval(nel_fft_z(:))
       max_elm = max(max_x,max_y,max_z)

       if(allocated(req_r)) deallocate(req_r)
       if(allocated(req_s)) deallocate(req_s)
       if(allocated(sta_r)) deallocate(sta_r)
       if(allocated(sta_s)) deallocate(sta_s)
       if(allocated(wk_recvcnt)) deallocate(wk_recvcnt)
       if(allocated(wk_sendcnt)) deallocate(wk_sendcnt)
       allocate(req_r(0:nmrank-1), stat=ierr)
       allocate(req_s(0:nmrank-1), stat=ierr)
       allocate(sta_r(MPI_STATUS_SIZE,0:nmrank-1), stat=ierr)
       allocate(sta_s(MPI_STATUS_SIZE,0:nmrank-1), stat=ierr)
       allocate(wk_recvcnt(0:nmrank-1), stat=ierr)
       allocate(wk_sendcnt(0:nmrank-1), stat=ierr)
!      allocate(wk_recvrank(0:nmrank-1), stat=ierr)
!      allocate(wk_sendrank(0:nmrank-1), stat=ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_FFT_Direct_XYZ_3D :  Not allocate '
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 51, ierr)
        endif

       if(allocated(z2y_recv)) deallocate(z2y_recv)
       if(allocated(z2y_send)) deallocate(z2y_send)
       allocate(z2y_recv(2,0:nmrank-1))
       allocate(z2y_send(2,0:nmrank-1))
       z2y_send = 0
       z2y_recv = 0
       wk_recvcnt = 0
       wk_sendcnt = 0
       do i = 1, nel_fft_y(myrank)
          irank = map_fft_z(mp_fft_y(i)) - 1
          wk_recvcnt(irank) = wk_recvcnt(irank) + 1
       enddo
       do i = 1, nel_fft_z(myrank)
          irank = map_fft_y(mp_fft_z(i)) - 1
          wk_sendcnt(irank) = wk_sendcnt(irank) + 1
       enddo
       k = 0
       do i = 0, nmrank - 1
          if(wk_recvcnt(i) /= 0) then
             k = k + 1
             z2y_recv(1,i) = wk_recvcnt(i)
             z2y_recv(2,i) = k
          endif
       enddo
       z2y_rrank = k
       k = 0
       do i = 0, nmrank - 1
          if(wk_sendcnt(i) /= 0) then
             k = k + 1
             z2y_send(1,i) = wk_sendcnt(i)
             z2y_send(2,i) = k
          endif
       enddo
       z2y_srank = k
       z2y_rmax = maxval(wk_recvcnt)
       z2y_smax = maxval(wk_sendcnt)
       z2y_srmax = max(z2y_rmax,z2y_smax)

       if(allocated(y2x_recv)) deallocate(y2x_recv)
       if(allocated(y2x_send)) deallocate(y2x_send)
       allocate(y2x_recv(2,0:nmrank-1))
       allocate(y2x_send(2,0:nmrank-1))
       y2x_send = 0
       y2x_recv = 0
       wk_recvcnt = 0
       wk_sendcnt = 0
       do i = 1, nel_fft_x(myrank)
          irank = map_fft_y(mp_fft_x(i)) - 1
          wk_recvcnt(irank) = wk_recvcnt(irank) + 1
       enddo
       do i = 1, nel_fft_y(myrank)
          irank = map_fft_x(mp_fft_y(i)) - 1
          wk_sendcnt(irank) = wk_sendcnt(irank) + 1
       enddo
       k = 0
       do i = 0, nmrank - 1
          if(wk_recvcnt(i) /= 0) then
             k = k + 1
             y2x_recv(1,i) = wk_recvcnt(i)
             y2x_recv(2,i) = k
          endif
       enddo
       y2x_rrank = k
       k = 0
      do i = 0, nmrank - 1
          if(wk_sendcnt(i) /= 0) then
             k = k + 1
             y2x_send(1,i) = wk_sendcnt(i)
             y2x_send(2,i) = k
          endif
       enddo
       y2x_srank = k
       y2x_rmax = maxval(wk_recvcnt)
       y2x_smax = maxval(wk_sendcnt)
       y2x_srmax = max(y2x_rmax,y2x_smax)

#ifdef FFT_ALLTOALL
       max_send(1) = z2y_rmax
       max_send(2) = z2y_smax
       call mpi_allreduce(max_send,max_recv,2,mpi_integer,mpi_max,mpi_fft_yz_world,ierr)
       z2y_rmax = max_recv(1)
       z2y_smax = max_recv(2)
       max_send(1) = y2x_rmax
       max_send(2) = y2x_smax
       call mpi_allreduce(max_send,max_recv,2,mpi_integer,mpi_max,mpi_fft_xy_world,ierr)
       y2x_rmax = max_recv(1)
       y2x_smax = max_recv(2)
#endif
       if (ipri > 1) then
          write(nfout,'("m_FFT_Direct_XYZ_3D   --   myrank_g=",i4)') myrank_g
          write(nfout,'("z2y_send")')
          write(nfout,'(10(i8,", "))') (z2y_send(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (z2y_send(2,i),i=0,nmrank-1)
          write(nfout,'("z2y_recv")')
          write(nfout,'(10(i8,", "))') (z2y_recv(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (z2y_recv(2,i),i=0,nmrank-1)
          write(nfout,'("y2x_send")')
          write(nfout,'(10(i8,", "))') (y2x_send(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (y2x_send(2,i),i=0,nmrank-1)
          write(nfout,'("y2x_recv")')
          write(nfout,'(10(i8,", "))') (y2x_recv(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (y2x_recv(2,i),i=0,nmrank-1)
       endif

#ifdef __TIMER_DO__
  call timer_sta(119)
#endif
       nx = xyz_fft_z(2,1)-xyz_fft_z(1,1)+1
       ny = xyz_fft_z(2,2)-xyz_fft_z(1,2)+1
       nz = xyz_fft_z(2,3)-xyz_fft_z(1,3)+1
!      if(kimg==1) then
          call dfftw_plan_many_dft    (planz1,     1,      nz, nx*ny/2,  &
         &                             wk_afft_l, nz, nx*ny/2,       1,  &
         &                             wk_afft_l, nz, nx*ny/2,       1,  &
         &                             FFTW_FLAG, FFTW_MEASURE )
          call dfftw_plan_many_dft    (planz1_1d,  1,      nz, nx/2   , &
         &                             wk_afft_l, nz, nx*ny/2,       1, &
         &                             wk_afft_l, nz, nx*ny/2,       1, &
         &                             FFTW_FLAG, FFTW_MEASURE )
!      else
          call dfftw_plan_many_dft    (planz2,     1,    nz, nx*ny,  &
         &                             wk_afft_l, nz, nx*ny,     1,  &
         &                             wk_afft_l, nz, nx*ny,     1,  &
         &                             FFTW_FLAG, FFTW_MEASURE )
          call dfftw_plan_many_dft    (planz2_1d,  1,    nz, nx   ,  &
         &                             wk_afft_l, nz, nx*ny,     1,  &
         &                             wk_afft_l, nz, nx*ny,     1,  &
         &                             FFTW_FLAG, FFTW_MEASURE )
!      endif

       nx = xyz_fft_y(2,1)-xyz_fft_y(1,1)+1
       ny = xyz_fft_y(2,2)-xyz_fft_y(1,2)+1
       nz = xyz_fft_y(2,3)-xyz_fft_y(1,3)+1
!      if(kimg==1) then
          call dfftw_plan_many_dft    (plany1   ,  1,   ny, nx/2,  &
         &                             wk_afft_l, ny, nx/2,    1,  &
         &                             wk_afft_l, ny, nx/2,    1,  &
         &                             FFTW_FLAG, FFTW_MEASURE )
!      else
          call dfftw_plan_many_dft    (plany2   ,  1, ny, nx,   &
         &                             wk_afft_l, ny, nx,  1,   &
         &                             wk_afft_l, ny, nx,  1,   &
         &                             FFTW_FLAG, FFTW_MEASURE )
!      end if

       nx = xyz_fft_x(2,1)-xyz_fft_x(1,1)+1
       ny = xyz_fft_x(2,2)-xyz_fft_x(1,2)+1
       nz = xyz_fft_x(2,3)-xyz_fft_x(1,3)+1
!      if(kimg==1) then
          call dfftw_plan_many_dft_c2r(planx1   ,   1, nx-2, ny*nz,   &
         &                             wk_afft_l, nx,     1,  nx/2,   &
         &                             wk_afft_l, nx,     1,    nx,   &
         &                             FFTW_MEASURE )
!?        call dfftw_plan_many_dft_c2r(planx1_1d,   1,   nx-2, ny   , &
!?       &                             wk_afft_l,  nx-2,    1,    nx, &
!?       &                             wk_afft_l, (nx-2)/2, 1,  nx/2, &
!?       &                             FFTW_MEASURE )
          call dfftw_plan_dft_c2r_1d(planx1_1d, nx-2, wk_afft_l, wk_afft_l, FFTW_FLAG, FFTW_MEASURE)
!      else
          call dfftw_plan_many_dft    (planx2   ,  1, nx, ny*nz,   &
         &                             wk_afft_l, nx,  1,    nx,   &
         &                             wk_afft_l, nx,  1,    nx,   &
         &                             FFTW_FLAG, FFTW_MEASURE )
          call dfftw_plan_many_dft    (planx2_1d,  1, nx, ny   ,   &
         &                             wk_afft_l, nx,  1,    nx,   &
         &                             wk_afft_l, nx,  1,    nx,   &
         &                             FFTW_FLAG, FFTW_MEASURE )
!      endif
#ifdef __TIMER_DO__
  call timer_end(119)
#endif

       firstcall_direct_xyz_3d = .false.
#ifdef __FAPP__
    call fapp_stop('fftwf_direct_firstcall',1,1)
#endif
       go to 1000
#ifdef __TIMER_DO__
  call timer_end(118)
#endif
    endif

#ifdef FFT_ALLTOALL
    if (iesize /= savesize) then
#else
    if (iesize > savesize) then
#endif
#ifdef __TIMER_DO__
  call timer_sta(120)
#endif
       if (allocated(wk_recv1)) deallocate(wk_recv1)
       if (allocated(wk_send1)) deallocate(wk_send1)
       if (allocated(wk_recv2)) deallocate(wk_recv2)
       if (allocated(wk_send2)) deallocate(wk_send2)
       allocate(wk_recv1(z2y_rmax*kimg*iesize,z2y_rrank), stat=ierr)
       allocate(wk_send1(z2y_smax*kimg*iesize,z2y_srank), stat=ierr)
       allocate(wk_recv2(y2x_rmax*kimg*iesize,y2x_rrank), stat=ierr)
       allocate(wk_send2(y2x_smax*kimg*iesize,y2x_srank), stat=ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_FFT_Direct_XYZ_3D :  Not allocate '
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 52, ierr)
        endif
       savesize = iesize
#ifdef __TIMER_DO__
  call timer_end(120)
#endif
    end if

    if (np_fft_x==0)then
       nx = 0
       nx = 0
       nz = 0
       go to 1000
    endif

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP          SHARED(     nfout,wk_afft_l,iesize,req_r,req_s,sta_r,sta_s,    &
!$OMP                 wk_recv1,wk_send1,wk_recv2,wk_send2,                    &
!$OMP                 planx1,planx2,planz1,planz2,plany1,plany2,              &
!$OMP                 planx1_1d,planx2_1d,planz1_1d,planz2_1d,plany1_1d,plany2_1d, &
!$OMP                 z2y_recv, z2y_send, y2x_recv, y2x_send,itag,            &
!$OMP                 xyz_fft_y,xyz_fft_z,xyz_fft_x,fft_Z_x_dim,nis_fft_Z_x,  &
!$OMP                 nie_fft_Z_x,fft_X_z_dim,nis_fft_X_z,nie_fft_X_z,        &
!$OMP                 fft_Y_z_dim,nis_fft_Y_z,nie_fft_Y_z,fft_Z_y_dim,        &
!$OMP                 nie_fft_Z_y,nis_fft_Z_y,                                &
!$OMP                 nie_fft_X_y,nis_fft_X_y,fft_X_y_dim,                    &
!$OMP                 nie_fft_Y_x,nis_fft_Y_x,fft_Y_x_dim,                    &
#ifdef FFT_ALLTOALL
!$OMP                 mpi_fft_yz_world, z2y_rmax, z2y_smax,                   &
!$OMP                 mpi_fft_xy_world, y2x_rmax, y2x_smax,                   &
#endif
!$OMP                 mpi_comm,myrank,nmrank,kimg,ierr                    )   &
!$OMP          PRIVATE(nx,ny,nz,i,j,k,l,ri,ib,iadd,lrank,icnt_send,icnt_recv, &
!$OMP                  nxp,nyp,nzp,nisx,niex,                                  &
!$OMP                  irank,nnx,nny,nnz,iadd0,iadd1,nis,nie,jadd,kadd,plan   )
!
! Z-axis (x-y div)
!
    nx = xyz_fft_z(2,1)-xyz_fft_z(1,1)+1
    ny = xyz_fft_z(2,2)-xyz_fft_z(1,2)+1
    nz = xyz_fft_z(2,3)-xyz_fft_z(1,3)+1
    nzp = nz

#ifdef __FAPP__
    call fapp_start('fftwf_direct_fftw1',1,1)
#endif
#ifdef __TIMER_DO__
  call timer_sta(251)
#endif
    if (kimg == 1) then
     if (iesize > 1) then
!$OMP DO
       do ib = 1, iesize
          call dfftw_execute_dft(planz1,wk_afft_l(1,ib),wk_afft_l(1,ib))
!x     end do
!xOMP DO
!x     do ib = 1, iesize
          do l = 1, fft_Y_z_dim
             nis = nis_fft_Y_z(l)
             nie = nie_fft_Y_z(l)
             nnz = nie - nis + 1
             iadd0 = nx/2*ny*nnz*(ib-1)*2
             do k = nis, nie
                do j = 1, ny
                   do i = 1, nx/2
                      iadd  = (i+(j-1)*nx/2+(k-1)*nx*ny/2)*2
                      iadd1 = iadd0+(i+(j-1)*nx/2+(k-nis)*nx*ny/2)*2
                      wk_send1(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send1(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
     else
       do ib = 1, iesize
!$OMP DO
          do j = 1, ny
             iadd = (1+(j-1)*nx/2)*2-1
             call dfftw_execute_dft(planz1_1d,wk_afft_l(iadd,ib),wk_afft_l(iadd,ib))
          end do
       end do
       do ib = 1, iesize
!$OMP DO
          do l = 1, fft_Y_z_dim
             nis = nis_fft_Y_z(l)
             nie = nie_fft_Y_z(l)
             nnz = nie - nis + 1
             iadd0 = nx/2*ny*nnz*(ib-1)*2
             do k = nis, nie
                do j = 1, ny
                   do i = 1, nx/2
                      iadd  = (i+(j-1)*nx/2+(k-1)*nx*ny/2)*2
                      iadd1 = iadd0+(i+(j-1)*nx/2+(k-nis)*nx*ny/2)*2
                      wk_send1(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send1(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
     end if
    else
     if (iesize > 1) then
!$OMP DO
       do ib = 1, iesize
          call dfftw_execute_dft(planz2,wk_afft_l(1,ib),wk_afft_l(1,ib))
!x     end do
!xOMP DO
!x     do ib = 1, iesize
          do l = 1, fft_Y_z_dim
             nis = nis_fft_Y_z(l)
             nie = nie_fft_Y_z(l)
             nnz = nie - nis + 1
             iadd0 = nx*ny*nnz*(ib-1)*2
             do k = nis, nie
                do j = 1, ny
                   do i = 1, nx
                      iadd  = (i+(j-1)*nx+(k-1)*nx*ny)*2
                      iadd1 = iadd0+(i+(j-1)*nx+(k-nis)*nx*ny)*2
                      wk_send1(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send1(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
     else
       do ib = 1, iesize
!$OMP DO
          do j = 1, ny
             iadd = (1+(j-1)*nx)*2-1
             call dfftw_execute_dft(planz2_1d,wk_afft_l(iadd,ib),wk_afft_l(iadd,ib))
          enddo
!o     end do
!o     do ib = 1, iesize
!$OMP DO
          do l = 1, fft_Y_z_dim
             nis = nis_fft_Y_z(l)
             nie = nie_fft_Y_z(l)
             nnz = nie - nis + 1
             iadd0 = nx*ny*nnz*(ib-1)*2
             do k = nis, nie
                do j = 1, ny
                   do i = 1, nx
                      iadd  = (i+(j-1)*nx+(k-1)*nx*ny)*2
                      iadd1 = iadd0+(i+(j-1)*nx+(k-nis)*nx*ny)*2
                      wk_send1(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send1(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
     end if
    end if
#ifdef __FAPP__
    call fapp_stop('fftwf_direct_fftw1',1,1)
#endif
#ifdef __TIMER_DO__
  call timer_end(251)
#endif

#ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_fft_yz_world)
  call timer_sta(123)
#endif

#ifdef __FAPP__
    call fapp_start('fftwf_direct_comm1',1,1)
#endif
!$OMP MASTER
    call MPI_ALLTOALL(wk_send1, z2y_smax*kimg*iesize, mpi_double_precision,   &
   &                  wk_recv1, z2y_rmax*kimg*iesize, mpi_double_precision,   &
   &                                                  mpi_fft_yz_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Direct_XYZ_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 55, ierr)
     endif
!$OMP END MASTER
!$OMP BARRIER

#ifdef __TIMER_DO__
  call timer_sta(124)
#endif
#ifdef __TIMER_DO__
  call timer_end(124)
#endif

#else
!else ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_comm)
  call timer_sta(123)
#endif

!$OMP MASTER
    icnt_recv = 0
    lrank = myrank + 1
    if (lrank > (nmrank-1)) lrank = 0
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank-1)) lrank = 0
       if ((y2z_recv(1,lrank) /= 0)) then
          call mpi_irecv(wk_recv1(1,z2y_recv(2,lrank)), z2y_recv(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_r(icnt_recv), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Direct_XYZ_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 55, ierr)
           endif
          icnt_recv = icnt_recv + 1
       endif
    enddo

    icnt_send = 0
    lrank = myrank
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank-1)) lrank = 0
       if ((y2z_send(1,lrank) /= 0)) then
          call mpi_isend(wk_send1(1,z2y_send(2,lrank)), z2y_send(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_s(icnt_send), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Direct_XYZ_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 56, ierr)
           endif
          icnt_send = icnt_send + 1
       endif
    enddo

#ifdef __TIMER_DO__
  call timer_sta(124)
#endif
#ifdef __TIMER_DO__
  call timer_end(124)
#endif

    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Direct_XYZ_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 57, ierr)
     endif

    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Direct_XYZ_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 58, ierr)
     endif
!$OMP END MASTER
!$OMP BARRIER

#endif
!endif ifdef FFT_ALLTOALL

#ifdef __FAPP__
    call fapp_stop('fftwf_direct_comm1',1,1)
#endif

#ifdef __TIMER_COMM__
  call timer_end(123)
#endif

#ifdef __FAPP__
    call fapp_start('fftwf_direct_fftw2',1,1)
#endif
    nx = xyz_fft_y(2,1)-xyz_fft_y(1,1)+1
    ny = xyz_fft_y(2,2)-xyz_fft_y(1,2)+1
    nz = xyz_fft_y(2,3)-xyz_fft_y(1,3)+1
    nyp = ny
#ifdef __TIMER_DO__
  call timer_sta(252)
#endif
    if (kimg == 1) then
     if (iesize > 1) then
!$OMP DO
       do ib = 1, iesize
          do k = 1, nz
             do l = 1, fft_Z_y_dim
                nis = nis_fft_Z_y(l)
                nie = nie_fft_Z_y(l)
                nny = nie - nis + 1
                iadd0 = nx/2*nny*nz*(ib-1)*2
                do j = nis, nie
                   do i = 1, nx/2
                      iadd1 = iadd0 + (i+(j-nis)*nx/2+(k-1)*nx*nny/2)*2
                      iadd = (i+(j-1)*nx/2+(k-1)*nx*nyp/2)*2
                      wk_afft_l(iadd-1,ib) = wk_recv1(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv1(iadd1  ,l)
                   end do
                end do
             end do
          end do
!x     end do
!
! Y-axis (z-x div)
!
!xOMP DO
!x     do ib = 1, iesize
          do k = 1, nz
             call dfftw_execute_dft(plany1,wk_afft_l((nx/2*nyp*(k-1)+1)*2-1,ib),wk_afft_l((nx/2*nyp*(k-1)+1)*2-1,ib))
          end do
!x     end do
!xOMP DO
!x     do ib = 1, iesize
          do k = 1, nz
             do l = 1, fft_X_y_dim
                nis = nis_fft_X_y(l)
                nie = nie_fft_X_y(l)
                nny = nie - nis + 1
                iadd0 = nx/2*nny*nz*(ib-1)*2
                do j = nis, nie
                   do i = 1, nx/2
                      iadd  = (i+(j-1)*nx/2+(k-1)*nx*nyp/2)*2
                      iadd1 = iadd0 + (i+(j-nis)*nx/2+(k-1)*nx*nny/2)*2
                      wk_send2(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send2(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
     else
       do ib = 1, iesize
!$OMP DO
          do k = 1, nz
             do l = 1, fft_Z_y_dim
                nis = nis_fft_Z_y(l)
                nie = nie_fft_Z_y(l)
                nny = nie - nis + 1
                iadd0 = nx/2*nny*nz*(ib-1)*2
                do j = nis, nie
                   do i = 1, nx/2
                      iadd1 = iadd0 + (i+(j-nis)*nx/2+(k-1)*nx*nny/2)*2
                      iadd = (i+(j-1)*nx/2+(k-1)*nx*nyp/2)*2
                      wk_afft_l(iadd-1,ib) = wk_recv1(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv1(iadd1  ,l)
                   end do
                end do
             end do
          end do
!o     end do
!
! Y-axis (z-x div)
!
!o     do ib = 1, iesize
!$OMP DO
          do k = 1, nz
             call dfftw_execute_dft(plany1,wk_afft_l((nx/2*nyp*(k-1)+1)*2-1,ib),wk_afft_l((nx/2*nyp*(k-1)+1)*2-1,ib))
          end do
!o     end do
!o     do ib = 1, iesize
!$OMP DO
          do k = 1, nz
             do l = 1, fft_X_y_dim
                nis = nis_fft_X_y(l)
                nie = nie_fft_X_y(l)
                nny = nie - nis + 1
                iadd0 = nx/2*nny*nz*(ib-1)*2
                do j = nis, nie
                   do i = 1, nx/2
                      iadd  = (i+(j-1)*nx/2+(k-1)*nx*nyp/2)*2
                      iadd1 = iadd0 + (i+(j-nis)*nx/2+(k-1)*nx*nny/2)*2
                      wk_send2(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send2(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
     end if
    else
     if (iesize > 1) then
!$OMP DO
       do ib = 1, iesize
          do k = 1, nz
             do l = 1, fft_Z_y_dim
                nis = nis_fft_Z_y(l)
                nie = nie_fft_Z_y(l)
                nny = nie - nis + 1
                iadd0 = nx*nny*nz*(ib-1)*2
                do j = nis, nie
                   do i = 1, nx
                      iadd1 = iadd0+(i+(j-nis)*nx+(k-1)*nx*nny)*2
                      iadd  = (i+(j-1)*nx+(k-1)*nx*nyp)*2
                      wk_afft_l(iadd-1,ib) = wk_recv1(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv1(iadd1  ,l)
                   end do
                end do
             end do
          end do
!x     end do
!
! Y-axis (z-x div)
!
!xOMP DO
!x     do ib = 1, iesize
          do k = 1, nz
             call dfftw_execute_dft(plany2,wk_afft_l((nx*nyp*(k-1)+1)*2-1,ib),wk_afft_l((nx*nyp*(k-1)+1)*2-1,ib))
          end do
!x     end do
!xOMP DO
!x     do ib = 1, iesize
          do k = 1, nz
             do l = 1, fft_X_y_dim
                nis = nis_fft_X_y(l)
                nie = nie_fft_X_y(l)
                nny = nie - nis + 1
                iadd0 = nx*nny*nz*(ib-1)*2
                do j = nis, nie
                   do i = 1, nx
                      iadd  = (i+(j-1)*nx+(k-1)*nx*nyp)*2
                      iadd1 = iadd0+(i+(j-nis)*nx+(k-1)*nx*nny)*2
                      wk_send2(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send2(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
     else
       do ib = 1, iesize
!$OMP DO
          do k = 1, nz
             do l = 1, fft_Z_y_dim
                nis = nis_fft_Z_y(l)
                nie = nie_fft_Z_y(l)
                nny = nie - nis + 1
                iadd0 = nx*nny*nz*(ib-1)*2
                do j = nis, nie
                   do i = 1, nx
                      iadd1 = iadd0+(i+(j-nis)*nx+(k-1)*nx*nny)*2
                      iadd  = (i+(j-1)*nx+(k-1)*nx*nyp)*2
                      wk_afft_l(iadd-1,ib) = wk_recv1(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv1(iadd1  ,l)
                   end do
                end do
             end do
!o        end do
!o     end do
!
! Y-axis (z-x div)
!
!o     do ib = 1, iesize
!o        do k = 1, nz
             call dfftw_execute_dft(plany2,wk_afft_l((nx*nyp*(k-1)+1)*2-1,ib),wk_afft_l((nx*nyp*(k-1)+1)*2-1,ib))
!o        end do
!o     end do
!o     do ib = 1, iesize
!o        do k = 1, nz
             do l = 1, fft_X_y_dim
                nis = nis_fft_X_y(l)
                nie = nie_fft_X_y(l)
                nny = nie - nis + 1
                iadd0 = nx*nny*nz*(ib-1)*2
                do j = nis, nie
                   do i = 1, nx
                      iadd  = (i+(j-1)*nx+(k-1)*nx*nyp)*2
                      iadd1 = iadd0+(i+(j-nis)*nx+(k-1)*nx*nny)*2
                      wk_send2(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send2(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
     end if
    end if
#ifdef __TIMER_DO__
  call timer_end(252)
#endif
#ifdef __FAPP__
    call fapp_stop('fftwf_direct_fftw2',1,1)
#endif

#ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_fft_xy_world)
  call timer_sta(127)
#endif

#ifdef __FAPP__
    call fapp_start('fftwf_direct_comm2',1,1)
#endif
!$OMP MASTER
    call MPI_ALLTOALL(wk_send2, y2x_smax*kimg*iesize, mpi_double_precision,   &
   &                  wk_recv2, y2x_rmax*kimg*iesize, mpi_double_precision,   &
   &                                                  mpi_fft_xy_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Direct_XYZ_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 60, ierr)
     endif
!$OMP END MASTER
!$OMP BARRIER

#ifdef __TIMER_DO__
  call timer_sta(128)
#endif
#ifdef __TIMER_DO__
  call timer_end(128)
#endif

#else
!else ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_comm)
  call timer_sta(127)
#endif

!$OMP MASTER
    icnt_recv = 0
    lrank = myrank + 1
    if (lrank > (nmrank-1)) lrank = 0
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank-1)) lrank = 0
       if ((z2x_recv(1,lrank) /= 0)) then
          call mpi_irecv(wk_recv2(1,z2x_recv(2,lrank)), z2x_recv(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_r(icnt_recv), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Direct_XYZ_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 60, ierr)
           endif
          icnt_recv = icnt_recv + 1
       endif
    enddo

    icnt_send = 0
    lrank = myrank
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank-1)) lrank = 0
       if ((z2x_send(1,lrank) /= 0)) then
          call mpi_isend(wk_send2(1,z2x_send(2,lrank)), z2x_send(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_s(icnt_send), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Direct_XYZ_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 61, ierr)
           endif
          icnt_send = icnt_send + 1
       endif
    enddo

#ifdef __TIMER_DO__
  call timer_sta(128)
#endif
#ifdef __TIMER_DO__
  call timer_end(128)
#endif

    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Direct_XYZ_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 62, ierr)
     endif

    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Direct_XYZ_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 63, ierr)
     endif
!$OMP END MASTER
!$OMP BARRIER
#endif

#ifdef __FAPP__
    call fapp_stop('fftwf_direct_comm2',1,1)
#endif
!endif ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_end(127)
#endif

#ifdef __FAPP__
    call fapp_start('fftwf_direct_fftw3',1,1)
#endif
    nx = xyz_fft_x(2,1)-xyz_fft_x(1,1)+1
    ny = xyz_fft_x(2,2)-xyz_fft_x(1,2)+1
    nz = xyz_fft_x(2,3)-xyz_fft_x(1,3)+1
#ifdef __TIMER_DO__
  call timer_sta(253)
#endif
    if (kimg == 1) then
     if (iesize > 1) then
!$OMP DO
       do ib = 1, iesize
          do k = 1, nz
             do j = 1, ny
                do l = 1, fft_Y_x_dim
                   nis = nis_fft_Y_x(l)
                   nie = nie_fft_Y_x(l)
                   nnx = nie - nis + 1
                   if(mod(nis,2)>0) then
                     nisx = nis/2 + 1
                   else
                     nisx = nis/2
                   endif
                   niex = nie/2
                   iadd0 = nnx/2*ny*nz*(ib-1)*2
                   do i = nisx, niex
                      iadd1 = iadd0 + (i-nisx+1+(j-1)*nnx/2+(k-1)*nnx/2*ny)*2
                      iadd  = (i+(j-1)*nx/2+(k-1)*nx*ny/2)*2
                      wk_afft_l(iadd-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
          end do
!x     end do
!
! X-axis (y-z div)
!
!x     do ib = 1, iesize
          call dfftw_execute_dft_c2r(planx1,wk_afft_l(1,ib),wk_afft_l(1,ib))
       end do
     else
       do ib = 1, iesize
!$OMP DO
          do k = 1, nz
             do j = 1, ny
                do l = 1, fft_Y_x_dim
                   nis = nis_fft_Y_x(l)
                   nie = nie_fft_Y_x(l)
                   nnx = nie - nis + 1
                   if(mod(nis,2)>0) then
                     nisx = nis/2 + 1
                   else
                     nisx = nis/2
                   endif
                   niex = nie/2
                   iadd0 = nnx/2*ny*nz*(ib-1)*2
                   do i = nisx, niex
                      iadd1 = iadd0 + (i-nisx+1+(j-1)*nnx/2+(k-1)*nnx/2*ny)*2
                      iadd  = (i+(j-1)*nx/2+(k-1)*nx*ny/2)*2
                      wk_afft_l(iadd-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
          end do
!o     end do
!
! X-axis (y-z div)
!
!o     do ib = 1, iesize
!$OMP DO
          do k = 1, nz
             do j = 1, ny
                iadd = (1+(j-1)*nx/2+(k-1)*nx/2*ny)*2-1
                call dfftw_execute_dft_c2r(planx1_1d,wk_afft_l(iadd,ib),wk_afft_l(iadd,ib))
             end do
          end do
       end do
     end if
    else
     if (iesize > 1) then
!$OMP DO
       do ib = 1, iesize
          do k = 1, nz
             do j = 1, ny
                do l = 1, fft_Y_x_dim
                   nis = nis_fft_Y_x(l)
                   nie = nie_fft_Y_x(l)
                   nnx = nie - nis + 1
                   iadd0 = nnx*ny*nz*(ib-1)*2
                   do i = nis, nie
                      iadd1 = iadd0+(i-nis+1+(j-1)*nnx+(k-1)*nnx*ny)*2
                      iadd  = (i+(j-1)*nx+(k-1)*nx*ny)*2
                      wk_afft_l(iadd-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
          end do
!x     end do
!
! X-axis (y-z div)
!
!xOMP DO
!x     do ib = 1, iesize
          call dfftw_execute_dft    (planx2,wk_afft_l(1,ib),wk_afft_l(1,ib))
       end do
     else
       do ib = 1, iesize
!$OMP DO
          do k = 1, nz
             do j = 1, ny
                do l = 1, fft_Y_x_dim
                   nis = nis_fft_Y_x(l)
                   nie = nie_fft_Y_x(l)
                   nnx = nie - nis + 1
                   iadd0 = nnx*ny*nz*(ib-1)*2
                   do i = nis, nie
                      iadd1 = iadd0+(i-nis+1+(j-1)*nnx+(k-1)*nnx*ny)*2
                      iadd  = (i+(j-1)*nx+(k-1)*nx*ny)*2
                      wk_afft_l(iadd-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
!o        end do
!o     end do
!
! X-axis (y-z div)
!
!o     do ib = 1, iesize
!o        do k = 1, nz
             iadd = (1+(k-1)*nx*ny)*2-1
             call dfftw_execute_dft(planx2_1d,wk_afft_l(iadd,ib),wk_afft_l(iadd,ib))
          end do
       end do
     end if
    end if
!      do i=1,wk_size
!        write(30+mype,'(i,2f20.10)') i,wk_afft_l(2*i-1,1),wk_afft_l(2*i,1)
!      enddo
!      stop
#ifdef __FAPP__
    call fapp_stop('fftwf_direct_fftw3',1,1)
#endif
#ifdef __TIMER_DO__
  call timer_end(253)
#endif

!$OMP END PARALLEL

1000 continue

#ifdef __TIMER_SUB__
    call timer_end(105)
#endif
    call tstatc0_end(id_sname)

  end subroutine m_FFT_Direct_XYZ_3D
!------------------------------------------------------------------------------
#endif
!endif #ifdef FFTW_NOSTRIDE
#endif
!endif #ifdef FFT_USE_SSL2

#ifdef FFT_USE_SSL2
!------------------------------------------------------------------------------
  subroutine m_FFT_Inverse_XYZ_3D(nfout,wk_afft_l, wk_size, iesize)
    integer, intent(in) :: nfout, wk_size, iesize
    real(kind=DP), dimension(wk_size*kimg,iesize), intent(inout) :: wk_afft_l

    integer, parameter :: FFTW_MEASURE=0
    integer, parameter :: FFTW_ESTIMATE=64
    integer, parameter :: FFTW_RANK=1
    integer            :: FFTW_FLAG
    integer,dimension(1) :: NFFTW3, NEMBED, NEREAL
    integer :: nx, ny, nz, i, j, k, l, ri, ib, iadd, lrank
    integer :: irank, itag, icnt_send, icnt_recv
    integer,save, allocatable, dimension(:)   :: req_r, req_s
    integer,save, allocatable, dimension(:,:) :: sta_r, sta_s
    integer,save, allocatable, dimension(:)   :: wk_recvcnt, wk_sendcnt
    real(kind=DP),allocatable, dimension(:,:),save :: wk_recv1, wk_send1
    real(kind=DP),allocatable, dimension(:,:),save :: wk_recv2, wk_send2
    integer,      allocatable, dimension(:)   :: wk_mp_fft_y
    integer :: mpi_comm, myrank, nmrank
    integer(kind=8),save :: planx1 = 0, planx2 = 0
    integer(kind=8),save :: plany1 = 0, plany2 = 0
    integer(kind=8),save :: planz1 = 0, planz2 = 0
    integer(kind=8),save :: planx1_1d = 0, planx2_1d = 0
    integer(kind=8),save :: plany1_1d = 0, plany2_1d = 0
    integer(kind=8),save :: planz1_1d = 0, planz2_1d = 0
    integer(kind=8) :: plan
    integer :: max_x, max_y, max_z
    integer,save :: max_elm = 0, savesize = 0, expo_x = 0, expo_y = 0, expo_z = 0

    integer,save, allocatable, dimension(:,:) :: x2y_recv, x2y_send, y2z_recv, y2z_send
    integer,save :: x2y_rrank, x2y_srank, y2z_rrank, y2z_srank
    integer,save :: x2y_rmax, x2y_smax, x2y_srmax, y2z_rmax, y2z_smax, y2z_srmax

    integer ::  nnx, nny, nnz, iadd0, iadd1, nis, nie, jadd, kadd, nxp, nyp, nzp
    integer                  , dimension(3) :: isin,nsize
    real(kind=DP),allocatable, dimension(:),save :: wwx,wwy,wwz
    integer :: isw
#ifdef FFT_ALLTOALL
    integer,dimension(2) :: max_send,max_recv
#endif

#ifdef __TIMER_SUB__
    call timer_barrier(mpi_ke_world)
    call timer_sta(106)
#endif

    if (kimg == 1) then
       FFTW_FLAG = -1      ! FFTW_FORWARD
    else
       FFTW_FLAG = +1      ! FFTW_BACKWARD
    endif

    mpi_comm = mpi_ke_world
    myrank = myrank_g
    nmrank = nrank_g
    itag = 10

    if (np_fft_x==0)then
       nx = 0
       nx = 0
       nz = 0
       go to 2000
    endif

    if (firstcall_inverse_xyz_3d) then
       savesize = 0
#ifdef __TIMER_DO__
  call timer_sta(134)
#endif
       max_x = maxval(nel_fft_x(:))
       max_y = maxval(nel_fft_y(:))
       max_z = maxval(nel_fft_z(:))
       max_elm = max(max_x,max_y,max_z)
       if(allocated(req_r)) deallocate(req_r)
       if(allocated(req_s)) deallocate(req_s)
       if(allocated(sta_r)) deallocate(sta_r)
       if(allocated(sta_s)) deallocate(sta_s)
       if(allocated(wk_recvcnt)) deallocate(wk_recvcnt)
       if(allocated(wk_sendcnt)) deallocate(wk_sendcnt)
       allocate(req_r(0:nmrank-1), stat=ierr)
       allocate(req_s(0:nmrank-1), stat=ierr)
       allocate(sta_r(MPI_STATUS_SIZE,0:nmrank-1), stat=ierr)
       allocate(sta_s(MPI_STATUS_SIZE,0:nmrank-1), stat=ierr)
       allocate(wk_recvcnt(0:nmrank-1), stat=ierr)
       allocate(wk_sendcnt(0:nmrank-1), stat=ierr)
       if (ierr /= 0) then
          write(nfout,*)' m_FFT_Inverse_3D :  Not allocate '
          call flush(nfout)
          call mpi_abort(mpi_comm_world, 68, ierr)
       endif

       if(allocated(x2y_recv)) deallocate(x2y_recv)
       if(allocated(x2y_send)) deallocate(x2y_send)
       allocate(x2y_recv(2,0:nmrank-1))
       allocate(x2y_send(2,0:nmrank-1))
       x2y_send = 0
       x2y_recv = 0
       wk_recvcnt = 0
       wk_sendcnt = 0
       do i = 1, nel_fft_y(myrank)
          irank = map_fft_x(mp_fft_y(i)) - 1
          wk_recvcnt(irank) = wk_recvcnt(irank) + 1
       enddo
       do i = 1, nel_fft_x(myrank)
          irank = map_fft_y(mp_fft_x(i)) - 1
          wk_sendcnt(irank) = wk_sendcnt(irank) + 1
       enddo
       k = 0
       do i = 0, nmrank - 1
          if(wk_recvcnt(i) /= 0) then
             k = k + 1
             x2y_recv(1,i) = wk_recvcnt(i)
             x2y_recv(2,i) = k
          endif
       enddo
       x2y_rrank = k
       k = 0
       do i = 0, nmrank - 1
          if(wk_sendcnt(i) /= 0) then
             k = k + 1
             x2y_send(1,i) = wk_sendcnt(i)
             x2y_send(2,i) = k
          endif
       enddo
       x2y_srank = k
       x2y_rmax = maxval(wk_recvcnt)
       x2y_smax = maxval(wk_sendcnt)
       x2y_srmax = max(x2y_rmax,x2y_smax)

       if(allocated(y2z_recv)) deallocate(y2z_recv)
       if(allocated(y2z_send)) deallocate(y2z_send)
       allocate(y2z_recv(2,0:nmrank-1))
       allocate(y2z_send(2,0:nmrank-1))
       y2z_send = 0
       y2z_recv = 0
       wk_recvcnt = 0
       wk_sendcnt = 0

       nx = xyz_fft_y(2,1)-xyz_fft_y(1,1)+1
       ny = xyz_fft_y(2,2)-xyz_fft_y(1,2)+1
       nz = xyz_fft_y(2,3)-xyz_fft_y(1,3)+1
       allocate(wk_mp_fft_y(nx*ny*nz) ,stat=ierr)
       do k = 0, nz-1
          do j = 0, ny-1
             do i = 0, nx-1
                wk_mp_fft_y(i+j*nx+k*nx*ny+1) = mp_fft_y(i+k*nx+j*nx*nz+1)
             enddo
          enddo
       enddo

       do i = 1, nel_fft_z(myrank)
          irank = map_fft_y(mp_fft_z(i)) - 1
          wk_recvcnt(irank) = wk_recvcnt(irank) + 1
       enddo
       do i = 1, nel_fft_y(myrank)
          irank = map_fft_z(wk_mp_fft_y(i)) - 1
          wk_sendcnt(irank) = wk_sendcnt(irank) + 1
       enddo
       k = 0
       do i = 0, nmrank - 1
          if(wk_recvcnt(i) /= 0) then
             k = k + 1
             y2z_recv(1,i) = wk_recvcnt(i)
             y2z_recv(2,i) = k
          endif
       enddo
       y2z_rrank = k
       k = 0
      do i = 0, nmrank - 1
          if(wk_sendcnt(i) /= 0) then
             k = k + 1
             y2z_send(1,i) = wk_sendcnt(i)
             y2z_send(2,i) = k
          endif
       enddo
       y2z_srank = k
       y2z_rmax = maxval(wk_recvcnt)
       y2z_smax = maxval(wk_sendcnt)
       y2z_srmax = max(y2z_rmax,y2z_smax)

       deallocate(wk_mp_fft_y)

       nx = xyz_fft_x(2,1)-xyz_fft_x(1,1)+1
       ny = xyz_fft_y(2,2)-xyz_fft_y(1,2)+1
       nz = xyz_fft_z(2,3)-xyz_fft_z(1,3)+1
       if ((nx==64).or.(nx==128).or.(nx==256).or.(nx==512).or.(nx==1024).or.(nx==2048)) then
          expo_x = 1
       end if
       if ((ny==64).or.(ny==128).or.(ny==256).or.(ny==512).or.(ny==1024).or.(ny==2048)) then
          expo_y = 1
       end if
       if ((nz==64).or.(nz==128).or.(nz==256).or.(nz==512).or.(nz==1024).or.(nz==2048)) then
          expo_z = 1
       end if
#ifdef FFT_ALLTOALL
       max_send(1) = x2y_rmax
       max_send(2) = x2y_smax
       call mpi_allreduce(max_send,max_recv,2,mpi_integer,mpi_max,mpi_fft_xy_world,ierr)
       x2y_rmax = max_recv(1)
       x2y_smax = max_recv(2)
       max_send(1) = y2z_rmax
       max_send(2) = y2z_smax
       call mpi_allreduce(max_send,max_recv,2,mpi_integer,mpi_max,mpi_fft_yz_world,ierr)
       y2z_rmax = max_recv(1)
       y2z_smax = max_recv(2)
#endif
       if (ipri > 1) then
          write(nfout,'("m_FFT_Inverse_XYZ_3D   --   myrank_g=",i4)') myrank_g
          write(nfout,'("x2y_send")')
          write(nfout,'(10(i8,", "))') (x2y_send(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (x2y_send(2,i),i=0,nmrank-1)
          write(nfout,'("x2y_recv")')
          write(nfout,'(10(i8,", "))') (x2y_recv(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (x2y_recv(2,i),i=0,nmrank-1)
          write(nfout,'("y2z_send")')
          write(nfout,'(10(i8,", "))') (y2z_send(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (y2z_send(2,i),i=0,nmrank-1)
          write(nfout,'("y2z_recv")')
          write(nfout,'(10(i8,", "))') (y2z_recv(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (y2z_recv(2,i),i=0,nmrank-1)
          call flush(nfout)
       endif

       firstcall_inverse_xyz_3d = .false.
#ifdef __TIMER_DO__
  call timer_end(134)
#endif
    endif

#ifdef FFT_ALLTOALL
    if (iesize /= savesize) then
#else
    if (iesize > savesize) then
#endif
#ifdef __TIMER_ETC__
  call timer_sta(264)
#endif
       if (allocated(wk_recv1)) deallocate(wk_recv1)
       if (allocated(wk_send1)) deallocate(wk_send1)
       if (allocated(wk_recv2)) deallocate(wk_recv2)
       if (allocated(wk_send2)) deallocate(wk_send2)
       allocate(wk_recv1(x2y_rmax*kimg*iesize,x2y_rrank), stat=ierr)
       allocate(wk_send1(x2y_smax*kimg*iesize,x2y_srank), stat=ierr)
       allocate(wk_recv2(y2z_rmax*kimg*iesize,y2z_rrank), stat=ierr)
       allocate(wk_send2(y2z_smax*kimg*iesize,y2z_srank), stat=ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_FFT_Inverse_3D :  Not allocate '
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 69, ierr)
        endif
       savesize = iesize
#ifdef __TIMER_ETC__
  call timer_end(264)
#endif
    end if

 2000 continue
#ifdef __TIMER_SUB__
    call mpi_barrier(mpi_ke_world, ierr)
    call timer_sta(1498)
#endif
    if (np_fft_x==0)then
       nx = 0
       nx = 0
       nz = 0
       go to 1000
    endif

#ifdef __TIMER_SUB__
    call timer_sta(1499)
#endif
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP          SHARED(     nfout,wk_afft_l,iesize,req_r,req_s,sta_r,sta_s,    &
!$OMP                 wk_recv1,wk_send1,wk_recv2,wk_send2,                    &
!$OMP                 x2y_recv, x2y_send, y2z_recv, y2z_send,itag,            &
!$OMP                 xyz_fft_y,xyz_fft_z,xyz_fft_x,fft_Z_x_dim,nis_fft_Z_x,  &
!$OMP                 nie_fft_Z_x,fft_X_z_dim,nis_fft_X_z,nie_fft_X_z,        &
!$OMP                 fft_Y_z_dim,nis_fft_Y_z,nie_fft_Y_z,fft_Z_y_dim,        &
!$OMP                 nie_fft_Z_y,nis_fft_Z_y,                                &
!$OMP                 nie_fft_X_y,nis_fft_X_y,fft_X_y_dim,                    &
!$OMP                 nie_fft_Y_x,nis_fft_Y_x,fft_Y_x_dim,                    &
#ifdef FFT_ALLTOALL
!$OMP                 mpi_fft_xy_world, x2y_rmax, x2y_smax,                   &
!$OMP                 mpi_fft_yz_world, y2z_rmax, y2z_smax,                   &
#endif
!$OMP                 mpi_comm,myrank,nmrank,kimg,expo_x,expo_y,expo_z    )   &
!$OMP          PRIVATE(nx,ny,nz,i,j,k,l,ri,ib,iadd,lrank,icnt_send,icnt_recv, &
!$OMP                  nxp,nyp,nzp,nsize,isin,ierr,                           &
!$OMP                                              wwx,wwy,wwz,isw,           &
!$OMP                  irank,nnx,nny,nnz,iadd0,iadd1,nis,nie,jadd,kadd        )
!
! X-axis (y-z div)
!
    nx = xyz_fft_x(2,1)-xyz_fft_x(1,1)+1
    ny = xyz_fft_x(2,2)-xyz_fft_x(1,2)+1
    nz = xyz_fft_x(2,3)-xyz_fft_x(1,3)+1
    nxp = nx

#ifdef __TIMER_DO__
  call timer_sta(254)
#endif
    if (kimg==1) then
       nsize(1:3) = (/nx-2,ny,nz/)
       isin(1:3)  = (/1,0,0/)
!$OMP DO
       do ib = 1, iesize
          call DVMRF2(wk_afft_l(1,ib),nsize,3,isin,1,ierr)
       end do
!$OMP DO
       do l = 1, fft_Y_x_dim
          nis = nis_fft_Y_x(l)
          nie = nie_fft_Y_x(l)
          nnx = nie - nis + 1
          do ib = 1, iesize
             iadd0 = nnx*ny*nz*(ib-1)
             do k = 1, nz
                do j = 1, ny
                   do i = nis, nie
                      wk_send1(iadd0+(i-nis+1+(j-1)*nnx+(k-1)*nnx*ny),l) = wk_afft_l(i+(j-1)*nx+(k-1)*nx*ny,ib)
                   end do
                end do
             end do
          end do
       end do
    else
     if (expo_x==1) then
        allocate(wwx((2*nx+70)*kimg))
     end if
     if (iesize > 1) then
       nsize(1:3) = (/nx,ny,nz/)
       isin(1:3)  = (/-1,0,0/)
       isw = 1
!$OMP DO
       do ib = 1, iesize
          if (expo_x==1) then
             do k = 1, nz
                do j = 1, ny
                   call DVCFM1(wk_afft_l((1+(j-1)*nx+(k-1)*nx*ny)*2-1,ib),nx,isw,-1,wwx,ierr)
!x                 call DVMCF2(wk_afft_l((1+(j-1)*nx+(k-1)*nx*ny)*2-1,ib),nx,1,-1,ierr)
                end do
             end do
          else
             call DVMCF2(wk_afft_l(1,ib),nsize,3,isin,ierr)
          end if
!x     end do
!xOMP DO
!x     do ib = 1, iesize
          do k = 1, nz
             do j = 1, ny
                do l = 1, fft_Y_x_dim
                   nis = nis_fft_Y_x(l)
                   nie = nie_fft_Y_x(l)
                   nnx = nie - nis + 1
                   iadd0 = nnx*ny*nz*(ib-1)*2
                   do i = nis, nie
                      iadd  = (i+(j-1)*nxp+(k-1)*nxp*ny)*2
                      iadd1 = iadd0+(i-nis+1+(j-1)*nnx+(k-1)*nnx*ny)*2
                      wk_send1(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send1(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
     else
       isw = 1
       do ib = 1, iesize
          if (expo_x==1) then
!$OMP DO
             do k = 1, nz
                do j = 1, ny
                   call DVCFM1(wk_afft_l((1+(j-1)*nx+(k-1)*nx*ny)*2-1,ib),nx,isw,-1,wwx,ierr)
                end do
             end do
          else
!$OMP DO
             do k = 1, nz
                do j = 1, ny
                   call DVMCF2(wk_afft_l((1+(j-1)*nx+(k-1)*nx*ny)*2-1,ib),nx,1,-1,ierr)
                end do
             end do
          end if
       end do
       do ib = 1, iesize
!$OMP DO
          do k = 1, nz
             do j = 1, ny
                do l = 1, fft_Y_x_dim
                   nis = nis_fft_Y_x(l)
                   nie = nie_fft_Y_x(l)
                   nnx = nie - nis + 1
                   iadd0 = nnx*ny*nz*(ib-1)*2
                   do i = nis, nie
                      iadd  = (i+(j-1)*nxp+(k-1)*nxp*ny)*2
                      iadd1 = iadd0+(i-nis+1+(j-1)*nnx+(k-1)*nnx*ny)*2
                      wk_send1(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send1(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
     end if
     if (expo_x==1) then
        deallocate(wwx)
     end if
    end if
#ifdef __TIMER_DO__
  call timer_end(254)
#endif

#ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_fft_xy_world)
  call timer_sta(137)
#endif
!$OMP MASTER
    call MPI_ALLTOALL(wk_send1, x2y_smax*kimg*iesize, mpi_double_precision,   &
   &                  wk_recv1, x2y_rmax*kimg*iesize, mpi_double_precision,   &
   &                                                  mpi_fft_xy_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Inverse_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 70, ierr)
     endif
!$OMP END MASTER
!$OMP BARRIER
#ifdef __TIMER_DO__
  call timer_sta(138)
#endif
#ifdef __TIMER_DO__
  call timer_end(138)
#endif

#else
!else ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_comm)
  call timer_sta(137)
#endif
!$OMP MASTER
    icnt_recv = 0
    lrank = myrank + 1
    if (lrank > (nmrank -1)) lrank = 0
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank -1)) lrank = 0
       if ((x2y_recv(1,lrank) /= 0)) then
          call mpi_irecv(wk_recv1(1,x2y_recv(2,lrank)), x2y_recv(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_r(icnt_recv), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Inverse_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 70, ierr)
           endif
          icnt_recv = icnt_recv + 1
       endif
    enddo

    icnt_send = 0
    lrank = myrank
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank -1)) lrank = 0
       if ((x2y_send(1,lrank) /= 0)) then
          call mpi_isend(wk_send1(1,x2y_send(2,lrank)), x2y_send(1,lrank)*kimg*iesize, &
                         mpi_double_precision, lrank, itag, mpi_comm, req_s(icnt_send), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Inverse_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 71, ierr)
           endif
          icnt_send = icnt_send + 1
       endif
    enddo

#ifdef __TIMER_DO__
  call timer_sta(138)
#endif

#ifdef __TIMER_DO__
  call timer_end(138)
#endif

    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Inverse_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 72, ierr)
     endif

    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Inverse_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 73, ierr)
     endif
!$OMP END MASTER
!$OMP BARRIER
#endif
!endif ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_end(137)
#endif

    nx = xyz_fft_y(2,1)-xyz_fft_y(1,1)+1
    ny = xyz_fft_y(2,2)-xyz_fft_y(1,2)+1
    nz = xyz_fft_y(2,3)-xyz_fft_y(1,3)+1
    nyp = ny
!
! Y-axis (z-x div)
!
#ifdef __TIMER_DO__
  call timer_sta(255)
#endif
    if (kimg == 1) then
       nsize(1:3) = (/ny,nz,(nx/2)/)
       isin(1:3)  = (/1,0,0/)
! (x,y,z) -> (y,z,x)
!$OMP DO
       do ib = 1, iesize
          do i = 1, nx/2
             do k = 1, nz
                do l = 1, fft_X_y_dim
                   nis = nis_fft_X_y(l)
                   nie = nie_fft_X_y(l)
                   nny = nie - nis + 1
                   iadd0 = nx/2*nny*nz*(ib-1)*2
                   do j = nis, nie
                      iadd1 = iadd0 + (i+(j-nis)*nx/2+(k-1)*nx/2*nny)*2
                      iadd  = (j+(k-1)*nyp+(i-1)*nyp*nz)*2
                      wk_afft_l(iadd-1,ib) = wk_recv1(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv1(iadd1  ,l)
                   end do
                end do
             end do
          end do
!x     end do
!xOMP DO
!x     do ib = 1, iesize
          call DVMCF2(wk_afft_l(1,ib),nsize,3,isin,ierr)
!x     end do
! (y,z,x) -> (y,z,x)
!xOMP DO
!x     do ib = 1, iesize
          do i = 1, nx/2
             do k = 1, nz
                do l = 1, fft_Z_y_dim
                   nis = nis_fft_Z_y(l)
                   nie = nie_fft_Z_y(l)
                   nny = nie - nis + 1
                   iadd0 = nx/2*nny*nz*(ib-1)*2
                   do j = nis, nie
                      iadd  = (j+(k-1)*nyp+(i-1)*nyp*nz)*2
                      iadd1 = iadd0 + (j-nis+1+(k-1)*nny+(i-1)*nny*nz)*2
                      wk_send2(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send2(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
    else
     if (expo_y==1) then
        allocate(wwy((2*ny+70)*kimg))
     end if
     if (iesize > 1) then
       nsize(1:3) = (/ny,nz,nx/)
       isin(1:3)  = (/-1,0,0/)
       isw = 1
! (x,y,z) -> (y,z,x)
!$OMP DO
       do ib = 1, iesize
          do i = 1, nx
             do k = 1, nz
                do l = 1, fft_X_y_dim
                   nis = nis_fft_X_y(l)
                   nie = nie_fft_X_y(l)
                   nny = nie - nis + 1
                   iadd0 = nx*nny*nz*(ib-1)*2
                   do j = nis, nie
                      iadd1 = iadd0+(i+(j-nis)*nx+(k-1)*nx*nny)*2
                      iadd  = (j+(k-1)*nyp+(i-1)*nyp*nz)*2
                      wk_afft_l(iadd-1,ib) = wk_recv1(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv1(iadd1  ,l)
                   end do
                end do
             end do
          end do
!x     end do
!xOMP DO
!x     do ib = 1, iesize
          if (expo_y==1) then
             do i = 1, nx
                do k = 1, nz
                   call DVCFM1(wk_afft_l((1+(k-1)*ny+(i-1)*ny*nz)*2-1,ib),ny,isw,-1,wwy,ierr)
!x                 call DVMCF2(wk_afft_l((1+(k-1)*ny+(i-1)*ny*nz)*2-1,ib),ny,1,-1,ierr)
                end do
             end do
          else
             call DVMCF2(wk_afft_l(1,ib),nsize,3,isin,ierr)
          end if
!x     end do
! (y,z,x) -> (y,z,x)
!xOMP DO
!x     do ib = 1, iesize
          do i = 1, nx
             do k = 1, nz
                do l = 1, fft_Z_y_dim
                   nis = nis_fft_Z_y(l)
                   nie = nie_fft_Z_y(l)
                   nny = nie - nis + 1
                   iadd0 = nx*nny*nz*(ib-1)*2
                   do j = nis, nie
                      iadd  = (j+(k-1)*nyp+(i-1)*nyp*nz)*2
                      iadd1 = iadd0+(j-nis+1+(k-1)*nny+(i-1)*nny*nz)*2
                      wk_send2(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send2(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
     else
       isw = 1
! (x,y,z) -> (y,z,x)
!$OMP DO
       do ib = 1, iesize
          do i = 1, nx
             do k = 1, nz
                do l = 1, fft_X_y_dim
                   nis = nis_fft_X_y(l)
                   nie = nie_fft_X_y(l)
                   nny = nie - nis + 1
                   iadd0 = nx*nny*nz*(ib-1)*2
                   do j = nis, nie
                      iadd1 = iadd0+(i+(j-nis)*nx+(k-1)*nx*nny)*2
                      iadd  = (j+(k-1)*nyp+(i-1)*nyp*nz)*2
                      wk_afft_l(iadd-1,ib) = wk_recv1(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv1(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end do
       do ib = 1, iesize
          if (expo_y==1) then
!$OMP DO
             do i = 1, nx
                do k = 1, nz
                   call DVCFM1(wk_afft_l((1+(k-1)*ny+(i-1)*ny*nz)*2-1,ib),ny,isw,-1,wwy,ierr)
                end do
             end do
          else
!$OMP DO
             do i = 1, nx
                do k = 1, nz
                   call DVMCF2(wk_afft_l((1+(k-1)*ny+(i-1)*ny*nz)*2-1,ib),ny,1,-1,ierr)
                end do
             end do
          end if
       end do
! (y,z,x) -> (y,z,x)
       do ib = 1, iesize
!$OMP DO
          do i = 1, nx
             do k = 1, nz
                do l = 1, fft_Z_y_dim
                   nis = nis_fft_Z_y(l)
                   nie = nie_fft_Z_y(l)
                   nny = nie - nis + 1
                   iadd0 = nx*nny*nz*(ib-1)*2
                   do j = nis, nie
                      iadd  = (j+(k-1)*nyp+(i-1)*nyp*nz)*2
                      iadd1 = iadd0+(j-nis+1+(k-1)*nny+(i-1)*nny*nz)*2
                      wk_send2(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send2(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
     end if
     if (expo_y==1) then
        deallocate(wwy)
     end if
    end if
#ifdef __TIMER_DO__
  call timer_end(255)
#endif

#ifdef __TIMER_DO__
  call timer_sta(140)
#endif
#ifdef __TIMER_DO__
  call timer_end(140)
#endif

#ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_fft_yz_world)
  call timer_sta(142)
#endif
!$OMP MASTER
    call MPI_ALLTOALL(wk_send2, y2z_smax*kimg*iesize, mpi_double_precision,   &
   &                  wk_recv2, y2z_rmax*kimg*iesize, mpi_double_precision,   &
   &                                                  mpi_fft_yz_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Inverse_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 76, ierr)
     endif
!$OMP END MASTER
!$OMP BARRIER
#ifdef __TIMER_DO__
  call timer_sta(143)
#endif
#ifdef __TIMER_DO__
  call timer_end(143)
#endif

#else
!else ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_comm)
  call timer_sta(142)
#endif
!$OMP MASTER
    icnt_recv = 0
    lrank = myrank + 1
    if (lrank > (nmrank-1)) lrank = 0
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank-1)) lrank = 0
       if ((z2y_recv(1,lrank) /= 0)) then
          call mpi_irecv(wk_recv2(1,z2y_recv(2,lrank)), z2y_recv(1,lrank)*kimg*iesize, &
                         mpi_double_precision, lrank, itag, mpi_comm, req_r(icnt_recv), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Inverse_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 76, ierr)
           endif
          icnt_recv = icnt_recv + 1
       endif
    enddo

    icnt_send = 0
    lrank = myrank
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank -1)) lrank = 0
       if ((z2y_send(1,lrank) /= 0)) then
          call mpi_isend(wk_send2(1,z2y_send(2,lrank)), z2y_send(1,lrank)*kimg*iesize, &
                         mpi_double_precision, lrank, itag, mpi_comm, req_s(icnt_send), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Inverse_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 77, ierr)
           endif
          icnt_send = icnt_send + 1
       endif
    enddo

#ifdef __TIMER_DO__
  call timer_sta(143)
#endif
#ifdef __TIMER_DO__
  call timer_end(143)
#endif

    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Inverse_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 78, ierr)
     endif

    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Inverse_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 79, ierr)
     endif
!$OMP END MASTER
!$OMP BARRIER
#endif
!endif ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_end(142)
#endif

    nx = xyz_fft_z(2,1)-xyz_fft_z(1,1)+1
    ny = xyz_fft_z(2,2)-xyz_fft_z(1,2)+1
    nz = xyz_fft_z(2,3)-xyz_fft_z(1,3)+1
    nzp = nz
!
! Z-axis (x-y div)
!
#ifdef __TIMER_DO__
  call timer_sta(256)
#endif
    if (kimg == 1) then
       nsize(1:3) = (/nz,nx/2,ny/)
       isin(1:3)  = (/1,0,0/)
! (y,z,x) -> (z,x,y)
!$OMP DO
       do ib = 1, iesize
          do j = 1, ny
             do i = 1, nx/2
                do l = 1, fft_Y_z_dim
                   nis = nis_fft_Y_z(l)
                   nie = nie_fft_Y_z(l)
                   nnz = nie - nis + 1
                   iadd0 = nx/2*ny*nnz*(ib-1)*2
                   do k = nis, nie
                      iadd1 = iadd0 + (j+(k-nis)*ny+(i-1)*ny*nnz)*2
                      iadd  = (k+(i-1)*nzp+(j-1)*nzp*nx/2)*2
                      wk_afft_l(iadd-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
          end do
!x     end do
!xOMP DO
!x     do ib = 1, iesize
          call DVMCF2(wk_afft_l(1,ib),nsize,3,isin,ierr)
       end do
    else
     if (expo_z==1) then
        allocate(wwz((2*nz+70)*kimg))
     end if
     if (iesize > 1) then
       nsize(1:3) = (/nz,nx,ny/)
       isin(1:3)  = (/-1,0,0/)
       isw = 1
! (y,z,x) -> (z,x,y)
!$OMP DO
       do ib = 1, iesize
          do j = 1, ny
             do i = 1, nx
                do l = 1, fft_Y_z_dim
                   nis = nis_fft_Y_z(l)
                   nie = nie_fft_Y_z(l)
                   nnz = nie - nis + 1
                   iadd0 = nx*ny*nnz*(ib-1)*2
                   do k = nis, nie
                      iadd1 = iadd0+(j+(k-nis)*ny+(i-1)*ny*nnz)*2
                      iadd  = (k+(i-1)*nzp+(j-1)*nzp*nx)*2
                      wk_afft_l(iadd-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
          end do
!x     end do
!xOMP DO
!x     do ib = 1, iesize
          if (expo_z==1) then
             do j = 1, ny
                do i = 1, nx
                   call DVCFM1(wk_afft_l((1+(i-1)*nz+(j-1)*nz*nx)*2-1,ib),nz,isw,-1,wwz,ierr)
!x                 call DVMCF2(wk_afft_l((1+(i-1)*nz+(j-1)*nz*nx)*2-1,ib),nz,1,-1,ierr)
                end do
             end do
          else
             call DVMCF2(wk_afft_l(1,ib),nsize,3,isin,ierr)
          end if
       end do
     else
       isw = 1
! (y,z,x) -> (z,x,y)
       do ib = 1, iesize
!$OMP DO
          do j = 1, ny
             do i = 1, nx
                do l = 1, fft_Y_z_dim
                   nis = nis_fft_Y_z(l)
                   nie = nie_fft_Y_z(l)
                   nnz = nie - nis + 1
                   iadd0 = nx*ny*nnz*(ib-1)*2
                   do k = nis, nie
                      iadd1 = iadd0+(j+(k-nis)*ny+(i-1)*ny*nnz)*2
                      iadd  = (k+(i-1)*nzp+(j-1)*nzp*nx)*2
                      wk_afft_l(iadd-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end do
       do ib = 1, iesize
          if (expo_z==1) then
!$OMP DO
             do j = 1, ny
                do i = 1, nx
                   call DVCFM1(wk_afft_l((1+(i-1)*nz+(j-1)*nz*nx)*2-1,ib),nz,isw,-1,wwz,ierr)
                end do
             end do
          else
!$OMP DO
             do j = 1, ny
                do i = 1, nx
                   call DVMCF2(wk_afft_l((1+(i-1)*nz+(j-1)*nz*nx)*2-1,ib),nz,1,-1,ierr)
                end do
             end do
          end if
       end do
     end if
     if (expo_z==1) then
        deallocate(wwz)
     end if
    end if
#ifdef __TIMER_DO__
    call timer_end(256)
#endif

#ifdef __TIMER_DO__
  call timer_sta(145)
#endif
#ifdef __TIMER_DO__
  call timer_end(145)
#endif

#ifdef __TIMER_SUB__
  call timer_end(1498)
#endif

#ifdef __TIMER_SUB__
  call timer_end(1499)
#endif
!$OMP END PARALLEL
1000 continue

#ifdef __TIMER_SUB__
    call timer_end(106)
#endif
  end subroutine m_FFT_Inverse_XYZ_3D
!------------------------------------------------------------------------------
#else
!else #ifdef FFT_USE_SSL2
#ifdef FFTW_NOSTRIDE
!------------------------------------------------------------------------------
  subroutine m_FFT_Inverse_XYZ_3D(nfout,wk_afft_l, wk_size, iesize)
    integer, intent(in) :: nfout, wk_size, iesize
    real(kind=DP), dimension(wk_size*kimg,iesize), intent(inout) :: wk_afft_l

    integer, parameter :: FFTW_MEASURE=0
    integer, parameter :: FFTW_ESTIMATE=64
    integer, parameter :: FFTW_RANK=1
    integer            :: FFTW_FLAG
    integer,dimension(1) :: NFFTW3, NEMBED, NEREAL
    integer :: nx, ny, nz, i, j, k, l, ri, ib, iadd, lrank
    integer :: irank, itag, icnt_send, icnt_recv
    integer,save, allocatable, dimension(:)   :: req_r, req_s
    integer,save, allocatable, dimension(:,:) :: sta_r, sta_s
    integer,save, allocatable, dimension(:)   :: wk_recvcnt, wk_sendcnt
    real(kind=DP),allocatable, dimension(:,:),save :: wk_recv1, wk_send1
    real(kind=DP),allocatable, dimension(:,:),save :: wk_recv2, wk_send2
    integer,      allocatable, dimension(:)   :: wk_mp_fft_y
    integer :: mpi_comm, myrank, nmrank
    integer(kind=8),save :: planx1 = 0, planx2 = 0
    integer(kind=8),save :: plany1 = 0, plany2 = 0
    integer(kind=8),save :: planz1 = 0, planz2 = 0
    integer(kind=8),save :: planx1_1d = 0, planx2_1d = 0
    integer(kind=8),save :: plany1_1d = 0, plany2_1d = 0
    integer(kind=8),save :: planz1_1d = 0, planz2_1d = 0
    integer(kind=8) :: plan
    integer :: max_x, max_y, max_z
    integer,save :: max_elm = 0, savesize = 0

    integer,save, allocatable, dimension(:,:) :: x2y_recv, x2y_send, y2z_recv, y2z_send
    integer,save :: x2y_rrank, x2y_srank, y2z_rrank, y2z_srank
    integer,save :: x2y_rmax, x2y_smax, x2y_srmax, y2z_rmax, y2z_smax, y2z_srmax

    integer ::  nnx, nny, nnz, iadd0, iadd1, nis, nie, jadd, kadd, nxp, nyp, nzp
    integer                  , dimension(3) :: isin,nsize
#ifdef FFT_ALLTOALL
    integer,dimension(2) :: max_send,max_recv
#endif

#ifdef __TIMER_SUB__
    call timer_barrier(mpi_ke_world)
    call timer_sta(106)
#endif

    if (kimg == 1) then
       FFTW_FLAG = -1      ! FFTW_FORWARD
    else
       FFTW_FLAG = +1      ! FFTW_BACKWARD
    endif

    mpi_comm = mpi_ke_world
    myrank = myrank_g
    nmrank = nrank_g
    itag = 10

    if (np_fft_x==0)then
       nx = 0
       nx = 0
       nz = 0
!      go to 1000
       go to 2000
    endif

    if (firstcall_inverse_xyz_3d) then
       savesize = 0
#ifdef __TIMER_DO__
  call timer_sta(134)
#endif
       max_x = maxval(nel_fft_x(:))
       max_y = maxval(nel_fft_y(:))
       max_z = maxval(nel_fft_z(:))
       max_elm = max(max_x,max_y,max_z)
       if(allocated(req_r)) deallocate(req_r)
       if(allocated(req_s)) deallocate(req_s)
       if(allocated(sta_r)) deallocate(sta_r)
       if(allocated(sta_s)) deallocate(sta_s)
       if(allocated(wk_recvcnt)) deallocate(wk_recvcnt)
       if(allocated(wk_sendcnt)) deallocate(wk_sendcnt)
       allocate(req_r(0:nmrank-1), stat=ierr)
       allocate(req_s(0:nmrank-1), stat=ierr)
       allocate(sta_r(MPI_STATUS_SIZE,0:nmrank-1), stat=ierr)
       allocate(sta_s(MPI_STATUS_SIZE,0:nmrank-1), stat=ierr)
       allocate(wk_recvcnt(0:nmrank-1), stat=ierr)
       allocate(wk_sendcnt(0:nmrank-1), stat=ierr)
       if (ierr /= 0) then
          write(nfout,*)' m_FFT_Inverse_3D :  Not allocate '
          call flush(nfout)
          call mpi_abort(mpi_comm_world, 68, ierr)
       endif

       if(allocated(x2y_recv)) deallocate(x2y_recv)
       if(allocated(x2y_send)) deallocate(x2y_send)
       allocate(x2y_recv(2,0:nmrank-1))
       allocate(x2y_send(2,0:nmrank-1))
       x2y_send = 0
       x2y_recv = 0
       wk_recvcnt = 0
       wk_sendcnt = 0
       do i = 1, nel_fft_y(myrank)
          irank = map_fft_x(mp_fft_y(i)) - 1
          wk_recvcnt(irank) = wk_recvcnt(irank) + 1
       enddo
       do i = 1, nel_fft_x(myrank)
          irank = map_fft_y(mp_fft_x(i)) - 1
          wk_sendcnt(irank) = wk_sendcnt(irank) + 1
       enddo
       k = 0
       do i = 0, nmrank - 1
          if(wk_recvcnt(i) /= 0) then
             k = k + 1
             x2y_recv(1,i) = wk_recvcnt(i)
             x2y_recv(2,i) = k
          endif
       enddo
       x2y_rrank = k
       k = 0
       do i = 0, nmrank - 1
          if(wk_sendcnt(i) /= 0) then
             k = k + 1
             x2y_send(1,i) = wk_sendcnt(i)
             x2y_send(2,i) = k
          endif
       enddo
       x2y_srank = k
       x2y_rmax = maxval(wk_recvcnt)
       x2y_smax = maxval(wk_sendcnt)
       x2y_srmax = max(x2y_rmax,x2y_smax)

       if(allocated(y2z_recv)) deallocate(y2z_recv)
       if(allocated(y2z_send)) deallocate(y2z_send)
       allocate(y2z_recv(2,0:nmrank-1))
       allocate(y2z_send(2,0:nmrank-1))
       y2z_send = 0
       y2z_recv = 0
       wk_recvcnt = 0
       wk_sendcnt = 0

       nx = xyz_fft_y(2,1)-xyz_fft_y(1,1)+1
       ny = xyz_fft_y(2,2)-xyz_fft_y(1,2)+1
       nz = xyz_fft_y(2,3)-xyz_fft_y(1,3)+1
       allocate(wk_mp_fft_y(nx*ny*nz) ,stat=ierr)
       do k = 0, nz-1
          do j = 0, ny-1
             do i = 0, nx-1
                wk_mp_fft_y(i+j*nx+k*nx*ny+1) = mp_fft_y(i+k*nx+j*nx*nz+1)
             enddo
          enddo
       enddo

       do i = 1, nel_fft_z(myrank)
          irank = map_fft_y(mp_fft_z(i)) - 1
          wk_recvcnt(irank) = wk_recvcnt(irank) + 1
       enddo
       do i = 1, nel_fft_y(myrank)
          irank = map_fft_z(wk_mp_fft_y(i)) - 1
          wk_sendcnt(irank) = wk_sendcnt(irank) + 1
       enddo
       k = 0
       do i = 0, nmrank - 1
          if(wk_recvcnt(i) /= 0) then
             k = k + 1
             y2z_recv(1,i) = wk_recvcnt(i)
             y2z_recv(2,i) = k
          endif
       enddo
       y2z_rrank = k
       k = 0
      do i = 0, nmrank - 1
          if(wk_sendcnt(i) /= 0) then
             k = k + 1
             y2z_send(1,i) = wk_sendcnt(i)
             y2z_send(2,i) = k
          endif
       enddo
       y2z_srank = k
       y2z_rmax = maxval(wk_recvcnt)
       y2z_smax = maxval(wk_sendcnt)
       y2z_srmax = max(y2z_rmax,y2z_smax)

       deallocate(wk_mp_fft_y)

#ifdef FFT_ALLTOALL
       max_send(1) = x2y_rmax
       max_send(2) = x2y_smax
       call mpi_allreduce(max_send,max_recv,2,mpi_integer,mpi_max,mpi_fft_xy_world,ierr)
       x2y_rmax = max_recv(1)
       x2y_smax = max_recv(2)
       max_send(1) = y2z_rmax
       max_send(2) = y2z_smax
       call mpi_allreduce(max_send,max_recv,2,mpi_integer,mpi_max,mpi_fft_yz_world,ierr)
       y2z_rmax = max_recv(1)
       y2z_smax = max_recv(2)
#endif

       if (ipri > 1) then
          write(nfout,'("m_FFT_Inverse_XYZ_3D   --   myrank_g=",i4)') myrank_g
          write(nfout,'("x2y_send")')
          write(nfout,'(10(i8,", "))') (x2y_send(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (x2y_send(2,i),i=0,nmrank-1)
          write(nfout,'("x2y_recv")')
          write(nfout,'(10(i8,", "))') (x2y_recv(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (x2y_recv(2,i),i=0,nmrank-1)
          write(nfout,'("y2z_send")')
          write(nfout,'(10(i8,", "))') (y2z_send(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (y2z_send(2,i),i=0,nmrank-1)
          write(nfout,'("y2z_recv")')
          write(nfout,'(10(i8,", "))') (y2z_recv(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (y2z_recv(2,i),i=0,nmrank-1)
          call flush(nfout)
       endif

#ifdef __TIMER_DO__
  call timer_sta(135)
#endif
       nx = xyz_fft_x(2,1)-xyz_fft_x(1,1)+1
       ny = xyz_fft_x(2,2)-xyz_fft_x(1,2)+1
       nz = xyz_fft_x(2,3)-xyz_fft_x(1,3)+1
       if(kimg==1) then
          call dfftw_plan_many_dft_r2c(planx1   ,        1, nx-2, ny*nz,  &
         &                             wk_afft_l,     nx-2,    1,    nx,  &
         &                             wk_afft_l, (nx-2)/2,    1,  nx/2,  &
         &                             FFTW_MEASURE )
       else
          call dfftw_plan_many_dft    (planx2   ,  1, nx, ny*nz,  &
         &                             wk_afft_l, nx,  1,    nx,  &
         &                             wk_afft_l, nx,  1,    nx,  &
         &                             FFTW_FLAG, FFTW_MEASURE )
          call dfftw_plan_dft_1d(planx2_1d, nx, wk_afft_l, wk_afft_l, FFTW_FLAG, FFTW_MEASURE)
       endif

       nx = xyz_fft_y(2,1)-xyz_fft_y(1,1)+1
       ny = xyz_fft_y(2,2)-xyz_fft_y(1,2)+1
       nz = xyz_fft_y(2,3)-xyz_fft_y(1,3)+1
       if(kimg==1) then
          call dfftw_plan_many_dft    (plany1   ,  1, ny, nz*nx/2,  &
         &                             wk_afft_l, ny,  1,      ny,  &
         &                             wk_afft_l, ny,  1,      ny,  &
         &                             FFTW_FLAG, FFTW_MEASURE )
       else
          call dfftw_plan_many_dft    (plany2   ,  1, ny, nz*nx,  &
         &                             wk_afft_l, ny,  1,    ny,  &
         &                             wk_afft_l, ny,  1,    ny,  &
         &                             FFTW_FLAG, FFTW_MEASURE )
          call dfftw_plan_dft_1d(plany2_1d, ny, wk_afft_l, wk_afft_l, FFTW_FLAG, FFTW_MEASURE)
       end if

       nx = xyz_fft_z(2,1)-xyz_fft_z(1,1)+1
       ny = xyz_fft_z(2,2)-xyz_fft_z(1,2)+1
       nz = xyz_fft_z(2,3)-xyz_fft_z(1,3)+1
       if(kimg==1) then
          call dfftw_plan_many_dft    (planz1   ,  1, nz, nx*ny/2,  &
         &                             wk_afft_l, nz,  1,      nz,  &
         &                             wk_afft_l, nz,  1,      nz,  &
         &                             FFTW_FLAG, FFTW_MEASURE )
       else
          call dfftw_plan_many_dft    (planz2   ,  1, nz, nx*ny,  &
         &                             wk_afft_l, nz,  1,    nz,  &
         &                             wk_afft_l, nz,  1,    nz,  &
         &                             FFTW_FLAG, FFTW_MEASURE )
          call dfftw_plan_dft_1d(planz2_1d, nz, wk_afft_l, wk_afft_l, FFTW_FLAG, FFTW_MEASURE)
       endif
#ifdef __TIMER_DO__
  call timer_end(135)
#endif

       firstcall_inverse_xyz_3d = .false.
       go to 1000
#ifdef __TIMER_DO__
  call timer_end(134)
#endif
    endif

#ifdef FFT_ALLTOALL
    if (iesize /= savesize) then
#else
    if (iesize > savesize) then
#endif
#ifdef __TIMER_ETC__
  call timer_sta(264)
#endif
       if (allocated(wk_recv1)) deallocate(wk_recv1)
       if (allocated(wk_send1)) deallocate(wk_send1)
       if (allocated(wk_recv2)) deallocate(wk_recv2)
       if (allocated(wk_send2)) deallocate(wk_send2)
       allocate(wk_recv1(x2y_rmax*kimg*iesize,x2y_rrank), stat=ierr)
       allocate(wk_send1(x2y_smax*kimg*iesize,x2y_srank), stat=ierr)
       allocate(wk_recv2(y2z_rmax*kimg*iesize,y2z_rrank), stat=ierr)
       allocate(wk_send2(y2z_smax*kimg*iesize,y2z_srank), stat=ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_FFT_Inverse_3D :  Not allocate '
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 69, ierr)
        endif
       savesize = iesize
#ifdef __TIMER_ETC__
  call timer_end(264)
#endif
    end if

 2000 continue
#ifdef __TIMER_SUB__
    call mpi_barrier(mpi_ke_world, ierr)
    call timer_sta(1498)
#endif
    if (np_fft_x==0)then
       nx = 0
       nx = 0
       nz = 0
       go to 1000
    endif
#ifdef __TIMER_SUB__
    call timer_sta(1499)
#endif
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP          SHARED(     nfout,wk_afft_l,iesize,req_r,req_s,sta_r,sta_s,    &
!$OMP                 wk_recv1,wk_send1,wk_recv2,wk_send2,                    &
!$OMP                 planx1,planx2,planz1,planz2,plany1,plany2,              &
!$OMP                 planx1_1d,planx2_1d,planz1_1d,planz2_1d,plany1_1d,plany2_1d, &
!$OMP                 x2y_recv, x2y_send, y2z_recv, y2z_send,itag,            &
!$OMP                 xyz_fft_y,xyz_fft_z,xyz_fft_x,fft_Z_x_dim,nis_fft_Z_x,  &
!$OMP                 nie_fft_Z_x,fft_X_z_dim,nis_fft_X_z,nie_fft_X_z,        &
!$OMP                 fft_Y_z_dim,nis_fft_Y_z,nie_fft_Y_z,fft_Z_y_dim,        &
!$OMP                 nie_fft_Z_y,nis_fft_Z_y,                                &
!$OMP                 nie_fft_X_y,nis_fft_X_y,fft_X_y_dim,                    &
!$OMP                 nie_fft_Y_x,nis_fft_Y_x,fft_Y_x_dim,                    &
#ifdef FFT_ALLTOALL
!$OMP                 mpi_fft_xy_world, x2y_rmax, x2y_smax,                   &
!$OMP                 mpi_fft_yz_world, y2z_rmax, y2z_smax,                   &
#endif
!$OMP                 mpi_comm,myrank,nmrank,kimg,ierr                    )   &
!$OMP          PRIVATE(nx,ny,nz,i,j,k,l,ri,ib,iadd,lrank,icnt_send,icnt_recv, &
!$OMP                  nxp,nyp,nzp,                                           &
!$OMP                  irank,nnx,nny,nnz,iadd0,iadd1,nis,nie,jadd,kadd,plan   )
!
! X-axis (y-z div)
!
    nx = xyz_fft_x(2,1)-xyz_fft_x(1,1)+1
    ny = xyz_fft_x(2,2)-xyz_fft_x(1,2)+1
    nz = xyz_fft_x(2,3)-xyz_fft_x(1,3)+1
    nxp = nx
#ifdef __TIMER_DO__
  call timer_sta(254)
#endif
    if (kimg == 1) then
!$OMP DO
       do ib = 1, iesize
          call dfftw_execute_dft_r2c(planx1,wk_afft_l(1,ib),wk_afft_l(1,ib))
       enddo
!$OMP DO
       do l = 1, fft_Y_x_dim
          nis = nis_fft_Y_x(l)
          nie = nie_fft_Y_x(l)
          nnx = nie - nis + 1
          do ib = 1, iesize
             iadd0 = nnx*ny*nz*(ib-1)
             do k = 1, nz
                do j = 1, ny
                   do i = nis, nie
                      wk_send1(iadd0+(i-nis+1+(j-1)*nnx+(k-1)*nnx*ny),l) = wk_afft_l(i+(j-1)*nx+(k-1)*nx*ny,ib)
                   end do
                end do
             end do
          end do
       end do
    else
     if (iesize > 1) then
!$OMP DO
       do ib = 1, iesize
          call dfftw_execute_dft    (planx2,wk_afft_l(1,ib),wk_afft_l(1,ib))
!x     enddo
!xOMP DO
!x     do ib = 1, iesize
          do k = 1, nz
             do j = 1, ny
                do l = 1, fft_Y_x_dim
                   nis = nis_fft_Y_x(l)
                   nie = nie_fft_Y_x(l)
                   nnx = nie - nis + 1
                   iadd0 = nnx*ny*nz*(ib-1)*2
                   do i = nis, nie
                      iadd  = (i+(j-1)*nxp+(k-1)*nxp*ny)*2
                      iadd1 = iadd0+(i-nis+1+(j-1)*nnx+(k-1)*nnx*ny)*2
                      wk_send1(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send1(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
     else
       do ib = 1, iesize
!$OMP DO
          do k = 1, nz
             do j = 1, ny
                iadd = (1+(j-1)*nxp+(k-1)*nxp*ny)*2-1
                call dfftw_execute_dft(planx2_1d,wk_afft_l(iadd,ib),wk_afft_l(iadd,ib))
             end do
!o        end do
!o     enddo
!o     do ib = 1, iesize
!o        do k = 1, nz
             do j = 1, ny
                do l = 1, fft_Y_x_dim
                   nis = nis_fft_Y_x(l)
                   nie = nie_fft_Y_x(l)
                   nnx = nie - nis + 1
                   iadd0 = nnx*ny*nz*(ib-1)*2
                   do i = nis, nie
                      iadd  = (i+(j-1)*nxp+(k-1)*nxp*ny)*2
                      iadd1 = iadd0+(i-nis+1+(j-1)*nnx+(k-1)*nnx*ny)*2
                      wk_send1(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send1(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
     end if
    end if
#ifdef __TIMER_DO__
  call timer_end(254)
#endif

#ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_fft_xy_world)
  call timer_sta(137)
#endif

!$OMP MASTER
    call MPI_ALLTOALL(wk_send1, x2y_smax*kimg*iesize, mpi_double_precision,   &
   &                  wk_recv1, x2y_rmax*kimg*iesize, mpi_double_precision,   &
   &                                                  mpi_fft_xy_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Inverse_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 70, ierr)
     endif
!$OMP END MASTER
!$OMP BARRIER

#ifdef __TIMER_DO__
  call timer_sta(138)
#endif
#ifdef __TIMER_DO__
  call timer_end(138)
#endif

#else
!else ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_comm)
  call timer_sta(137)
#endif

!$OMP MASTER
    icnt_recv = 0
    lrank = myrank + 1
    if (lrank > (nmrank -1)) lrank = 0
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank -1)) lrank = 0
       if ((x2y_recv(1,lrank) /= 0)) then
          call mpi_irecv(wk_recv1(1,x2y_recv(2,lrank)), x2y_recv(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_r(icnt_recv), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Inverse_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 70, ierr)
           endif
          icnt_recv = icnt_recv + 1
       endif
    enddo

    icnt_send = 0
    lrank = myrank
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank -1)) lrank = 0
       if ((x2y_send(1,lrank) /= 0)) then
          call mpi_isend(wk_send1(1,x2y_send(2,lrank)), x2y_send(1,lrank)*kimg*iesize, &
                         mpi_double_precision, lrank, itag, mpi_comm, req_s(icnt_send), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Inverse_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 71, ierr)
           endif
          icnt_send = icnt_send + 1
       endif
    enddo

#ifdef __TIMER_DO__
  call timer_sta(138)
#endif

#ifdef __TIMER_DO__
  call timer_end(138)
#endif

    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Inverse_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 72, ierr)
     endif

    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Inverse_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 73, ierr)
     endif
!$OMP END MASTER
!$OMP BARRIER

#endif
!endif ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_end(137)
#endif

    nx = xyz_fft_y(2,1)-xyz_fft_y(1,1)+1
    ny = xyz_fft_y(2,2)-xyz_fft_y(1,2)+1
    nz = xyz_fft_y(2,3)-xyz_fft_y(1,3)+1
    nyp = ny
#ifdef __TIMER_DO__
  call timer_sta(255)
#endif
    if (kimg == 1) then
! (x,y,z) -> (y,z,x)
!$OMP DO
       do ib = 1, iesize
          do i = 1, nx/2
             do k = 1, nz
                do l = 1, fft_X_y_dim
                   nis = nis_fft_X_y(l)
                   nie = nie_fft_X_y(l)
                   nny = nie - nis + 1
                   iadd0 = nx/2*nny*nz*(ib-1)*2
                   do j = nis, nie
                      iadd1 = iadd0 + (i+(j-nis)*nx/2+(k-1)*nx/2*nny)*2
                      iadd  = (j+(k-1)*nyp+(i-1)*nyp*nz)*2
                      wk_afft_l(iadd-1,ib) = wk_recv1(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv1(iadd1  ,l)
                   end do
                end do
             end do
          end do
!x     end do
!
! Y-axis (z-x div)
!
!xOMP DO
!x     do ib = 1, iesize
          call dfftw_execute_dft(plany1,wk_afft_l(1,ib),wk_afft_l(1,ib))
!x     end do
! (y,z,x) -> (y,z,x)
!xOMP DO
!x     do ib = 1, iesize
          do i = 1, nx/2
             do k = 1, nz
                do l = 1, fft_Z_y_dim
                   nis = nis_fft_Z_y(l)
                   nie = nie_fft_Z_y(l)
                   nny = nie - nis + 1
                   iadd0 = nx/2*nny*nz*(ib-1)*2
                   do j = nis, nie
                      iadd  = (j+(k-1)*nyp+(i-1)*nyp*nz)*2
                      iadd1 = iadd0 + (j-nis+1+(k-1)*nny+(i-1)*nny*nz)*2
                      wk_send2(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send2(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
    else
     if (iesize > 1) then
! (x,y,z) -> (y,z,x)
!$OMP DO
       do ib = 1, iesize
          do i = 1, nx
             do k = 1, nz
                do l = 1, fft_X_y_dim
                   nis = nis_fft_X_y(l)
                   nie = nie_fft_X_y(l)
                   nny = nie - nis + 1
                   iadd0 = nx*nny*nz*(ib-1)*2
                   do j = nis, nie
                      iadd1 = iadd0+(i+(j-nis)*nx+(k-1)*nx*nny)*2
                      iadd  = (j+(k-1)*nyp+(i-1)*nyp*nz)*2
                      wk_afft_l(iadd-1,ib) = wk_recv1(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv1(iadd1  ,l)
                   end do
                end do
             end do
          end do
!x     end do
!
! Y-axis (z-x div)
!
!xOMP DO
!x     do ib = 1, iesize
          call dfftw_execute_dft(plany2,wk_afft_l(1,ib),wk_afft_l(1,ib))
!x     end do
! (y,z,x) -> (y,z,x)
!xOMP DO
!x     do ib = 1, iesize
          do i = 1, nx
             do k = 1, nz
                do l = 1, fft_Z_y_dim
                   nis = nis_fft_Z_y(l)
                   nie = nie_fft_Z_y(l)
                   nny = nie - nis + 1
                   iadd0 = nx*nny*nz*(ib-1)*2
                   do j = nis, nie
                      iadd  = (j+(k-1)*nyp+(i-1)*nyp*nz)*2
                      iadd1 = iadd0+(j-nis+1+(k-1)*nny+(i-1)*nny*nz)*2
                      wk_send2(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send2(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
     else
! (x,y,z) -> (y,z,x)
       do ib = 1, iesize
!$OMP DO
          do i = 1, nx
             do k = 1, nz
                do l = 1, fft_X_y_dim
                   nis = nis_fft_X_y(l)
                   nie = nie_fft_X_y(l)
                   nny = nie - nis + 1
                   iadd0 = nx*nny*nz*(ib-1)*2
                   do j = nis, nie
                      iadd1 = iadd0+(i+(j-nis)*nx+(k-1)*nx*nny)*2
                      iadd  = (j+(k-1)*nyp+(i-1)*nyp*nz)*2
                      wk_afft_l(iadd-1,ib) = wk_recv1(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv1(iadd1  ,l)
                   end do
                end do
             end do
!o        end do
!o     end do
!
! Y-axis (z-x div)
!
!o     do ib = 1, iesize
!o        do i = 1, nx
             do k = 1, nz
                iadd = (1+(k-1)*nyp+(i-1)*nyp*nz)*2-1
                call dfftw_execute_dft(plany2_1d,wk_afft_l(iadd,ib),wk_afft_l(iadd,ib))
             end do
!o        end do
!o     end do
! (y,z,x) -> (y,z,x)
!o     do ib = 1, iesize
!o        do i = 1, nx
             do k = 1, nz
                do l = 1, fft_Z_y_dim
                   nis = nis_fft_Z_y(l)
                   nie = nie_fft_Z_y(l)
                   nny = nie - nis + 1
                   iadd0 = nx*nny*nz*(ib-1)*2
                   do j = nis, nie
                      iadd  = (j+(k-1)*nyp+(i-1)*nyp*nz)*2
                      iadd1 = iadd0+(j-nis+1+(k-1)*nny+(i-1)*nny*nz)*2
                      wk_send2(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send2(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
     end if
    end if

#ifdef __TIMER_DO__
  call timer_end(255)
#endif

#ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_fft_yz_world)
  call timer_sta(142)
#endif

!$OMP MASTER
    call MPI_ALLTOALL(wk_send2, y2z_smax*kimg*iesize, mpi_double_precision,   &
   &                  wk_recv2, y2z_rmax*kimg*iesize, mpi_double_precision,   &
   &                                                  mpi_fft_yz_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Inverse_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 76, ierr)
     endif
!$OMP END MASTER
!$OMP BARRIER

#ifdef __TIMER_DO__
  call timer_sta(143)
#endif
#ifdef __TIMER_DO__
  call timer_end(143)
#endif

#else
!else ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_comm)
  call timer_sta(142)
#endif

!$OMP MASTER
    icnt_recv = 0
    lrank = myrank + 1
    if (lrank > (nmrank-1)) lrank = 0
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank-1)) lrank = 0
       if ((z2y_recv(1,lrank) /= 0)) then
          call mpi_irecv(wk_recv2(1,z2y_recv(2,lrank)), z2y_recv(1,lrank)*kimg*iesize, &
                         mpi_double_precision, lrank, itag, mpi_comm, req_r(icnt_recv), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Inverse_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 76, ierr)
           endif
          icnt_recv = icnt_recv + 1
       endif
    enddo

    icnt_send = 0
    lrank = myrank
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank -1)) lrank = 0
       if ((z2y_send(1,lrank) /= 0)) then
          call mpi_isend(wk_send2(1,z2y_send(2,lrank)), z2y_send(1,lrank)*kimg*iesize, &
                         mpi_double_precision, lrank, itag, mpi_comm, req_s(icnt_send), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Inverse_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 77, ierr)
           endif
          icnt_send = icnt_send + 1
       endif
    enddo

#ifdef __TIMER_DO__
  call timer_sta(143)
#endif
#ifdef __TIMER_DO__
  call timer_end(143)
#endif

    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Inverse_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 78, ierr)
     endif

    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Inverse_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 79, ierr)
     endif
!$OMP END MASTER
!$OMP BARRIER

#endif
!endif ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_end(142)
#endif

    nx = xyz_fft_z(2,1)-xyz_fft_z(1,1)+1
    ny = xyz_fft_z(2,2)-xyz_fft_z(1,2)+1
    nz = xyz_fft_z(2,3)-xyz_fft_z(1,3)+1
    nzp = nz
#ifdef __TIMER_DO__
  call timer_sta(256)
#endif
    if (kimg == 1) then
! (y,z,x) -> (z,x,y)
!$OMP DO
       do ib = 1, iesize
          do j = 1, ny
             do i = 1, nx/2
                do l = 1, fft_Y_z_dim
                   nis = nis_fft_Y_z(l)
                   nie = nie_fft_Y_z(l)
                   nnz = nie - nis + 1
                   iadd0 = nx/2*ny*nnz*(ib-1)*2
                   do k = nis, nie
                      iadd1 = iadd0 + (j+(k-nis)*ny+(i-1)*ny*nnz)*2
                      iadd  = (k+(i-1)*nzp+(j-1)*nzp*nx/2)*2
                      wk_afft_l(iadd-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
          end do
!x     end do
!
! Z-axis (x-y div)
!
!xOMP DO
!x     do ib = 1, iesize
          call dfftw_execute_dft(planz1,wk_afft_l(1,ib),wk_afft_l(1,ib))
       end do
    else
     if (iesize > 1) then
! (y,z,x) -> (z,x,y)
!$OMP DO
       do ib = 1, iesize
          do j = 1, ny
             do i = 1, nx
                do l = 1, fft_Y_z_dim
                   nis = nis_fft_Y_z(l)
                   nie = nie_fft_Y_z(l)
                   nnz = nie - nis + 1
                   iadd0 = nx*ny*nnz*(ib-1)*2
                   do k = nis, nie
                      iadd1 = iadd0+(j+(k-nis)*ny+(i-1)*ny*nnz)*2
                      iadd  = (k+(i-1)*nzp+(j-1)*nzp*nx)*2
                      wk_afft_l(iadd-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
          end do
!x     end do
!
! Z-axis (x-y div)
!
!xOMP DO
!x     do ib = 1, iesize
          call dfftw_execute_dft(planz2,wk_afft_l(1,ib),wk_afft_l(1,ib))
       end do
     else
! (y,z,x) -> (z,x,y)
       do ib = 1, iesize
!$OMP DO
          do j = 1, ny
             do i = 1, nx
                do l = 1, fft_Y_z_dim
                   nis = nis_fft_Y_z(l)
                   nie = nie_fft_Y_z(l)
                   nnz = nie - nis + 1
                   iadd0 = nx*ny*nnz*(ib-1)*2
                   do k = nis, nie
                      iadd1 = iadd0+(j+(k-nis)*ny+(i-1)*ny*nnz)*2
                      iadd  = (k+(i-1)*nzp+(j-1)*nzp*nx)*2
                      wk_afft_l(iadd-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
!o        end do
!o     end do
!
! Z-axis (x-y div)
!
!o     do ib = 1, iesize
!o        do j = 1, ny
             do i = 1, nx
                iadd = (1+(i-1)*nzp+(j-1)*nzp*nx)*2-1
                call dfftw_execute_dft(planz2_1d,wk_afft_l(iadd,ib),wk_afft_l(iadd,ib))
             end do
          end do
       end do
     end if
    end if
#ifdef __TIMER_DO__
  call timer_end(256)
#endif

!$OMP END PARALLEL

#ifdef __TIMER_SUB__
  call timer_end(1498)
#endif

#ifdef __TIMER_SUB__
  call timer_end(1499)
#endif

1000 continue

#ifdef __TIMER_SUB__
    call timer_end(106)
#endif
  end subroutine m_FFT_Inverse_XYZ_3D
!------------------------------------------------------------------------------
#else
!else #ifdef FFTW_NOSTRIDE
!------------------------------------------------------------------------------
  subroutine m_FFT_Inverse_XYZ_3D(nfout,wk_afft_l, wk_size, iesize)
#ifdef KMATH_FFT3D
    use m_Const_Parameters, only : PAI2
    use kmath_fft3d_mod, only : KMATH_FFT3D_Transform
#endif
    integer, intent(in) :: nfout, wk_size, iesize
    real(kind=DP), dimension(wk_size*kimg,iesize), intent(inout) :: wk_afft_l

    integer, parameter :: FFTW_MEASURE=0
    integer, parameter :: FFTW_ESTIMATE=64
    integer, parameter :: FFTW_RANK=1
    integer            :: FFTW_FLAG
    integer,dimension(1) :: NFFTW3, NEMBED, NEREAL
    integer :: nx, ny, nz, i, j, k, l, ri, ib, iadd, lrank
    integer :: irank, itag, icnt_send, icnt_recv
    integer,save, allocatable, dimension(:)   :: req_r, req_s
    integer,save, allocatable, dimension(:,:) :: sta_r, sta_s
    integer,save, allocatable, dimension(:)   :: wk_recvcnt, wk_sendcnt
    real(kind=DP),allocatable, dimension(:,:),save :: wk_recv1, wk_send1
    real(kind=DP),allocatable, dimension(:,:),save :: wk_recv2, wk_send2
    integer,      allocatable, dimension(:)   :: wk_mp_fft_y
    integer :: mpi_comm, myrank, nmrank
    integer(kind=8),save :: planx1 = 0, planx2 = 0
    integer(kind=8),save :: plany1 = 0, plany2 = 0
    integer(kind=8),save :: planz1 = 0, planz2 = 0
    integer(kind=8),save :: planx1_1d = 0, planx2_1d = 0
    integer(kind=8),save :: plany1_1d = 0, plany2_1d = 0
    integer(kind=8),save :: planz1_1d = 0, planz2_1d = 0
    integer(kind=8) :: plan
    integer :: max_x, max_y, max_z
    integer,save :: max_elm = 0, savesize = 0

    integer,save, allocatable, dimension(:,:) :: x2y_recv, x2y_send, y2z_recv, y2z_send
    integer,save :: x2y_rrank, x2y_srank, y2z_rrank, y2z_srank
    integer,save :: x2y_rmax, x2y_smax, x2y_srmax, y2z_rmax, y2z_smax, y2z_srmax

    integer ::  nnx, nny, nnz, iadd0, iadd1, nis, nie, jadd, kadd, nxp, nyp, nzp
    integer                  , dimension(3) :: isin,nsize
#ifdef FFT_ALLTOALL
    integer,dimension(2) :: max_send,max_recv
#endif
#ifdef KMATH_FFT3D
    complex(kind=DP), allocatable, dimension(:) :: wk_afft_l_in, wk_afft_l_out
    integer :: icount
    integer, dimension(3) :: nsize2
    integer :: nx1,ny1,nz1,nx2,ny2,nz2
#endif
    integer :: id_sname = -1

    call tstatc0_begin('m_FFT_Inverse ',id_sname)
#ifdef __TIMER_SUB__
    call timer_barrier(mpi_ke_world)
    call timer_sta(106)
#endif
    if (kimg == 1) then
       FFTW_FLAG = -1      ! FFTW_FORWARD
    else
       FFTW_FLAG = +1      ! FFTW_BACKWARD
    endif

#ifdef KMATH_FFT3D
    if(sw_kmath_fft3d==ON .and. .not.firstcall_inverse_xyz_3d) then
      nsize2 = fft_box_size_WF(:,1)

      nx1 = nsize2(1)/nproc_fft3d_inverse(1)
      ny1 = nsize2(2)/nproc_fft3d_inverse(2)
      nz1 = nsize2(3)/nproc_fft3d_inverse(3)
      if (MOD(nsize2(1),nproc_fft3d_inverse(1)) /= 0) &
        nx1 = nx1 + 1
      if (MOD(nsize2(2),nproc_fft3d_inverse(2)) /= 0) &
        ny1 = ny1 + 1
      if (MOD(nsize2(3),nproc_fft3d_inverse(3)) /= 0) &
        nz1 = nz1 + 1

      nx2 = nsize2(1)/nproc_fft3d_inverse(2)
      if (MOD(nsize2(1),nproc_fft3d_inverse(2)) /= 0) &
        nx2 = nx2 + 1
      ny2 = nsize2(2)/nproc_fft3d_inverse(1)
      if (MOD(nsize2(2),nproc_fft3d_inverse(1)) /= 0) &
        ny2 = ny2 + 1
      if (MOD(ny2,nproc_fft3d_inverse(3)) /= 0) then
        ny2 = ny2/nproc_fft3d_inverse(3) + 1
      else
        ny2 = ny2/nproc_fft3d_inverse(3)
      end if
      nz2 = nsize2(3)

      nx = xyz_fft_x(2,1)-xyz_fft_x(1,1)+1
      ny = xyz_fft_x(2,2)-xyz_fft_x(1,2)+1
      nz = xyz_fft_x(2,3)-xyz_fft_x(1,3)+1
      allocate(wk_afft_l_in (nx1*ny1*nz1));wk_afft_l_in=cmplx(0.d0,0.d0)
      allocate(wk_afft_l_out(nx2*ny2*nz2));wk_afft_l_out=cmplx(0.d0,0.d0)
      do k=1, nz
        do j=1, ny
          do l=1, nx
             iadd0 = l + nx1*(j-1) + nx1*ny1*(k-1)
             iadd1 = l + nx*(j-1) + nx*ny*(k-1)
             wk_afft_l_in(iadd1) = cmplx(wk_afft_l(2*iadd1-1,1),-wk_afft_l(2*iadd1,1))
          enddo
        enddo
      enddo
      call KMATH_FFT3D_Transform(kmath3d_handle_wf_inverse, wk_afft_l_in, wk_afft_l_out,.false.)
      nx = xyz_fft_z(2,1)-xyz_fft_z(1,1)+1
      ny = xyz_fft_z(2,2)-xyz_fft_z(1,2)+1
      nz = xyz_fft_z(2,3)-xyz_fft_z(1,3)+1
      do k=1, nz
        do j=1, ny
          do l=1, nx
             iadd0 = k + nz*(l-1) + nz*nx*(j-1)
             iadd1 = l + nx2*(j-1) + nx2*ny2*(k-1)
             wk_afft_l_in(iadd1) = wk_afft_l_out(iadd0)
          enddo
        enddo
      enddo

      do i=1,wk_size
        wk_afft_l(2*i-1,1) =  real (wk_afft_l_in(i))
        wk_afft_l(2*i,1)   = -aimag(wk_afft_l_in(i))
      enddo
      deallocate(wk_afft_l_in)
      deallocate(wk_afft_l_out)
      call tstatc0_end(id_sname)
      return
    endif
#endif

    mpi_comm = mpi_ke_world
    myrank = myrank_g
    nmrank = nrank_g
    itag = 10

    if (np_fft_x==0)then
       nx = 0
       nx = 0
       nz = 0
!      go to 1000
       go to 2000
    endif

    if (firstcall_inverse_xyz_3d) then
#ifdef __FAPP__
    call fapp_start('fftwf_inverse_firstcall',1,1)
#endif
       if(planx1    /= 0) call dfftw_destroy_plan(planx1)
       if(planx1_1d /= 0) call dfftw_destroy_plan(planx1_1d)
       if(planx2    /= 0) call dfftw_destroy_plan(planx2)
       if(planx2_1d /= 0) call dfftw_destroy_plan(planx2_1d)
       if(plany1    /= 0) call dfftw_destroy_plan(plany1)
       if(plany2    /= 0) call dfftw_destroy_plan(plany2)
       if(planz1    /= 0) call dfftw_destroy_plan(planz1)
       if(planz1_1d /= 0) call dfftw_destroy_plan(planz1_1d)
       if(planz2    /= 0) call dfftw_destroy_plan(planz2)
       if(planz2_1d /= 0) call dfftw_destroy_plan(planz2_1d)
       savesize = 0
#ifdef __TIMER_DO__
  call timer_sta(134)
#endif
       max_x = maxval(nel_fft_x(:))
       max_y = maxval(nel_fft_y(:))
       max_z = maxval(nel_fft_z(:))
       max_elm = max(max_x,max_y,max_z)
       if(allocated(req_r)) deallocate(req_r)
       if(allocated(req_s)) deallocate(req_s)
       if(allocated(sta_r)) deallocate(sta_r)
       if(allocated(sta_s)) deallocate(sta_s)
       if(allocated(wk_recvcnt)) deallocate(wk_recvcnt)
       if(allocated(wk_sendcnt)) deallocate(wk_sendcnt)
       allocate(req_r(0:nmrank-1), stat=ierr)
       allocate(req_s(0:nmrank-1), stat=ierr)
       allocate(sta_r(MPI_STATUS_SIZE,0:nmrank-1), stat=ierr)
       allocate(sta_s(MPI_STATUS_SIZE,0:nmrank-1), stat=ierr)
       allocate(wk_recvcnt(0:nmrank-1), stat=ierr)
       allocate(wk_sendcnt(0:nmrank-1), stat=ierr)
       if (ierr /= 0) then
          write(nfout,*)' m_FFT_Inverse_3D :  Not allocate '
          call flush(nfout)
          call mpi_abort(mpi_comm_world, 68, ierr)
       endif

       if(allocated(x2y_recv)) deallocate(x2y_recv)
       if(allocated(x2y_send)) deallocate(x2y_send)
       allocate(x2y_recv(2,0:nmrank-1))
       allocate(x2y_send(2,0:nmrank-1))
       x2y_send = 0
       x2y_recv = 0
       wk_recvcnt = 0
       wk_sendcnt = 0
       do i = 1, nel_fft_y(myrank)
          irank = map_fft_x(mp_fft_y(i)) - 1
          wk_recvcnt(irank) = wk_recvcnt(irank) + 1
       enddo
       do i = 1, nel_fft_x(myrank)
          irank = map_fft_y(mp_fft_x(i)) - 1
          wk_sendcnt(irank) = wk_sendcnt(irank) + 1
       enddo
       k = 0
       do i = 0, nmrank - 1
          if(wk_recvcnt(i) /= 0) then
             k = k + 1
             x2y_recv(1,i) = wk_recvcnt(i)
             x2y_recv(2,i) = k
          endif
       enddo
       x2y_rrank = k
       k = 0
       do i = 0, nmrank - 1
          if(wk_sendcnt(i) /= 0) then
             k = k + 1
             x2y_send(1,i) = wk_sendcnt(i)
             x2y_send(2,i) = k
          endif
       enddo
       x2y_srank = k
       x2y_rmax = maxval(wk_recvcnt)
       x2y_smax = maxval(wk_sendcnt)
       x2y_srmax = max(x2y_rmax,x2y_smax)

       if(allocated(y2z_recv)) deallocate(y2z_recv)
       if(allocated(y2z_send)) deallocate(y2z_send)
       allocate(y2z_recv(2,0:nmrank-1))
       allocate(y2z_send(2,0:nmrank-1))
       y2z_send = 0
       y2z_recv = 0
       wk_recvcnt = 0
       wk_sendcnt = 0

       nx = xyz_fft_y(2,1)-xyz_fft_y(1,1)+1
       ny = xyz_fft_y(2,2)-xyz_fft_y(1,2)+1
       nz = xyz_fft_y(2,3)-xyz_fft_y(1,3)+1
       allocate(wk_mp_fft_y(nx*ny*nz) ,stat=ierr)
       do k = 0, nz-1
          do j = 0, ny-1
             do i = 0, nx-1
                wk_mp_fft_y(i+j*nx+k*nx*ny+1) = mp_fft_y(i+k*nx+j*nx*nz+1)
             enddo
          enddo
       enddo

       do i = 1, nel_fft_z(myrank)
          irank = map_fft_y(mp_fft_z(i)) - 1
          wk_recvcnt(irank) = wk_recvcnt(irank) + 1
       enddo
       do i = 1, nel_fft_y(myrank)
          irank = map_fft_z(wk_mp_fft_y(i)) - 1
          wk_sendcnt(irank) = wk_sendcnt(irank) + 1
       enddo
       k = 0
       do i = 0, nmrank - 1
          if(wk_recvcnt(i) /= 0) then
             k = k + 1
             y2z_recv(1,i) = wk_recvcnt(i)
             y2z_recv(2,i) = k
          endif
       enddo
       y2z_rrank = k
       k = 0
      do i = 0, nmrank - 1
          if(wk_sendcnt(i) /= 0) then
             k = k + 1
             y2z_send(1,i) = wk_sendcnt(i)
             y2z_send(2,i) = k
          endif
       enddo
       y2z_srank = k
       y2z_rmax = maxval(wk_recvcnt)
       y2z_smax = maxval(wk_sendcnt)
       y2z_srmax = max(y2z_rmax,y2z_smax)

       deallocate(wk_mp_fft_y)

#ifdef FFT_ALLTOALL
       max_send(1) = x2y_rmax
       max_send(2) = x2y_smax
       call mpi_allreduce(max_send,max_recv,2,mpi_integer,mpi_max,mpi_fft_xy_world,ierr)
       x2y_rmax = max_recv(1)
       x2y_smax = max_recv(2)
       max_send(1) = y2z_rmax
       max_send(2) = y2z_smax
       call mpi_allreduce(max_send,max_recv,2,mpi_integer,mpi_max,mpi_fft_yz_world,ierr)
       y2z_rmax = max_recv(1)
       y2z_smax = max_recv(2)
#endif

       if (ipri > 1) then
          write(nfout,'("m_FFT_Inverse_XYZ_3D   --   myrank_g=",i4)') myrank_g
          write(nfout,'("x2y_send")')
          write(nfout,'(10(i8,", "))') (x2y_send(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (x2y_send(2,i),i=0,nmrank-1)
          write(nfout,'("x2y_recv")')
          write(nfout,'(10(i8,", "))') (x2y_recv(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (x2y_recv(2,i),i=0,nmrank-1)
          write(nfout,'("y2z_send")')
          write(nfout,'(10(i8,", "))') (y2z_send(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (y2z_send(2,i),i=0,nmrank-1)
          write(nfout,'("y2z_recv")')
          write(nfout,'(10(i8,", "))') (y2z_recv(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (y2z_recv(2,i),i=0,nmrank-1)
          call flush(nfout)
       endif

#ifdef __TIMER_DO__
  call timer_sta(135)
#endif
       nx = xyz_fft_x(2,1)-xyz_fft_x(1,1)+1
       ny = xyz_fft_x(2,2)-xyz_fft_x(1,2)+1
       nz = xyz_fft_x(2,3)-xyz_fft_x(1,3)+1
       if(kimg==1) then
          call dfftw_plan_many_dft_r2c(planx1   ,        1, nx-2, ny*nz, &
         &                             wk_afft_l,     nx-2,    1,    nx, &
         &                             wk_afft_l, (nx-2)/2,    1,  nx/2, &
         &                             FFTW_MEASURE )
          call dfftw_plan_many_dft_r2c(planx1_1d,        1, nx-2, ny   , &
         &                             wk_afft_l,     nx-2,    1,    nx, &
         &                             wk_afft_l, (nx-2)/2,    1,  nx/2, &
         &                             FFTW_MEASURE )
!?        call dfftw_plan_dft_r2c_1d(planx1_1d, nx-2, wk_afft_l, wk_afft_l, FFTW_FLAG, FFTW_MEASURE)
       else
          call dfftw_plan_many_dft    (planx2   ,  1, nx, ny*nz, &
         &                             wk_afft_l, nx,  1,    nx, &
         &                             wk_afft_l, nx,  1,    nx, &
         &                             FFTW_FLAG, FFTW_MEASURE )
          call dfftw_plan_many_dft    (planx2_1d,  1, nx, ny   , &
         &                             wk_afft_l, nx,  1,    nx, &
         &                             wk_afft_l, nx,  1,    nx, &
         &                             FFTW_FLAG, FFTW_MEASURE )
       endif

       nx = xyz_fft_y(2,1)-xyz_fft_y(1,1)+1
       ny = xyz_fft_y(2,2)-xyz_fft_y(1,2)+1
       nz = xyz_fft_y(2,3)-xyz_fft_y(1,3)+1
       if(kimg==1) then
          call dfftw_plan_many_dft    (plany1   ,  1,   ny, nx/2, &
         &                             wk_afft_l, ny, nx/2,    1, &
         &                             wk_afft_l, ny, nx/2,    1, &
         &                             FFTW_FLAG, FFTW_MEASURE )
       else
          call dfftw_plan_many_dft    (plany2   ,  1, ny, nx, &
         &                             wk_afft_l, ny, nx,  1, &
         &                             wk_afft_l, ny, nx,  1, &
         &                             FFTW_FLAG, FFTW_MEASURE )
       end if

       nx = xyz_fft_z(2,1)-xyz_fft_z(1,1)+1
       ny = xyz_fft_z(2,2)-xyz_fft_z(1,2)+1
       nz = xyz_fft_z(2,3)-xyz_fft_z(1,3)+1
       if(kimg==1) then
          call dfftw_plan_many_dft    (planz1   ,  1,      nz, nx*ny/2, &
         &                             wk_afft_l, nz, nx*ny/2,       1, &
         &                             wk_afft_l, nz, nx*ny/2,       1, &
         &                             FFTW_FLAG, FFTW_MEASURE )
          call dfftw_plan_many_dft    (planz1_1d,  1,      nz, nx/2   , &
         &                             wk_afft_l, nz, nx*ny/2,       1, &
         &                             wk_afft_l, nz, nx*ny/2,       1, &
         &                             FFTW_FLAG, FFTW_MEASURE )
       else
          call dfftw_plan_many_dft    (planz2   ,  1,    nz, nx*ny, &
         &                             wk_afft_l, nz, nx*ny,     1, &
         &                             wk_afft_l, nz, nx*ny,     1, &
         &                             FFTW_FLAG, FFTW_MEASURE )
          call dfftw_plan_many_dft    (planz2_1d,  1, nz,    nx   , &
         &                             wk_afft_l, nz, nx*ny,     1, &
         &                             wk_afft_l, nz, nx*ny,     1, &
         &                             FFTW_FLAG, FFTW_MEASURE )
       endif
#ifdef __TIMER_DO__
  call timer_end(135)
#endif

       firstcall_inverse_xyz_3d = .false.
#ifdef __FAPP__
    call fapp_stop('fftwf_inverse_firstcall',1,1)
#endif
       go to 1000
#ifdef __TIMER_DO__
  call timer_end(134)
#endif
    endif

#ifdef FFT_ALLTOALL
    if (iesize /= savesize) then
#else
    if (iesize > savesize) then
#endif
#ifdef __TIMER_ETC__
  call timer_sta(264)
#endif
       if (allocated(wk_recv1)) deallocate(wk_recv1)
       if (allocated(wk_send1)) deallocate(wk_send1)
       if (allocated(wk_recv2)) deallocate(wk_recv2)
       if (allocated(wk_send2)) deallocate(wk_send2)
       allocate(wk_recv1(x2y_rmax*kimg*iesize,x2y_rrank), stat=ierr)
       allocate(wk_send1(x2y_smax*kimg*iesize,x2y_srank), stat=ierr)
       allocate(wk_recv2(y2z_rmax*kimg*iesize,y2z_rrank), stat=ierr)
       allocate(wk_send2(y2z_smax*kimg*iesize,y2z_srank), stat=ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_FFT_Inverse_3D :  Not allocate '
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 69, ierr)
        endif
       savesize = iesize
#ifdef __TIMER_ETC__
  call timer_end(264)
#endif
    end if

 2000 continue
#ifdef __TIMER_SUB__
    call mpi_barrier(mpi_ke_world, ierr)
    call timer_sta(1498)
#endif
    if (np_fft_x==0)then
       nx = 0
       nx = 0
       nz = 0
       go to 1000
    endif

#ifdef __TIMER_SUB__
    call timer_sta(1499)
#endif
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP          SHARED(     nfout,wk_afft_l,iesize,req_r,req_s,sta_r,sta_s,    &
!$OMP                 wk_recv1,wk_send1,wk_recv2,wk_send2,                    &
!$OMP                 planx1,planx2,planz1,planz2,plany1,plany2,              &
!$OMP                 planx1_1d,planz1_1d,planx2_1d,planz2_1d,                &
!$OMP                 x2y_recv, x2y_send, y2z_recv, y2z_send,itag,            &
!$OMP                 xyz_fft_y,xyz_fft_z,xyz_fft_x,fft_Z_x_dim,nis_fft_Z_x,  &
!$OMP                 nie_fft_Z_x,fft_X_z_dim,nis_fft_X_z,nie_fft_X_z,        &
!$OMP                 fft_Y_z_dim,nis_fft_Y_z,nie_fft_Y_z,fft_Z_y_dim,        &
!$OMP                 nie_fft_Z_y,nis_fft_Z_y,                                &
!$OMP                 nie_fft_X_y,nis_fft_X_y,fft_X_y_dim,                    &
!$OMP                 nie_fft_Y_x,nis_fft_Y_x,fft_Y_x_dim,                    &
#ifdef FFT_ALLTOALL
!$OMP                 mpi_fft_xy_world, x2y_rmax, x2y_smax,                   &
!$OMP                 mpi_fft_yz_world, y2z_rmax, y2z_smax,                   &
#endif
!$OMP                 mpi_comm,myrank,nmrank,kimg,ierr                    )   &
!$OMP          PRIVATE(nx,ny,nz,i,j,k,l,ri,ib,iadd,lrank,icnt_send,icnt_recv, &
!$OMP                  nxp,nyp,nzp,                                           &
!$OMP                  irank,nnx,nny,nnz,iadd0,iadd1,nis,nie,jadd,kadd,plan   )

    nx = xyz_fft_x(2,1)-xyz_fft_x(1,1)+1
    ny = xyz_fft_x(2,2)-xyz_fft_x(1,2)+1
    nz = xyz_fft_x(2,3)-xyz_fft_x(1,3)+1
    nxp = nx
#ifdef __TIMER_DO__
  call timer_sta(254)
#endif
#ifdef __FAPP__
    call fapp_start('fftwf_inverse_fftw1',1,1)
#endif
    if (kimg == 1) then
!
! X-axis (y-z div)
!
     if (iesize > 1) then
!$OMP DO
       do ib = 1, iesize
          call dfftw_execute_dft_r2c(planx1,wk_afft_l(1,ib),wk_afft_l(1,ib))
!o     enddo
!oOMP DO
!o     do ib = 1, iesize
          do k = 1, nz
             do j = 1, ny
                do l = 1, fft_Y_x_dim
                   nis = nis_fft_Y_x(l)
                   nie = nie_fft_Y_x(l)
                   nnx = nie - nis + 1
                   iadd0 = nnx*ny*nz*(ib-1)
                   do i = nis, nie
                      wk_send1(iadd0+(i-nis+1+(j-1)*nnx+(k-1)*nnx*ny),l) = wk_afft_l(i+(j-1)*nx+(k-1)*nx*ny,ib)
                   end do
                end do
             end do
          end do
       end do
     else
       do ib = 1, iesize
!$OMP DO
          do k = 1, nz
!?           do j = 1, ny
!?              iadd = (1+(j-1)*nx/2+(k-1)*nx*ny/2)*2-1
!?              call dfftw_execute_dft_r2c(planx1_1d,wk_afft_l(iadd,ib),wk_afft_l(iadd,ib))
!?           end do
             iadd = (1+(k-1)*nx*ny/2)*2-1
             call dfftw_execute_dft_r2c(planx1_1d,wk_afft_l(iadd,ib),wk_afft_l(iadd,ib))
!o        end do
!o     end do
!o     do ib = 1, iesize
!oOMP DO
!o        do k = 1, nz
             do j = 1, ny
                do l = 1, fft_Y_x_dim
                   nis = nis_fft_Y_x(l)
                   nie = nie_fft_Y_x(l)
                   nnx = nie - nis + 1
                   iadd0 = nnx*ny*nz*(ib-1)
                   do i = nis, nie
                      wk_send1(iadd0+(i-nis+1+(j-1)*nnx+(k-1)*nnx*ny),l) = wk_afft_l(i+(j-1)*nx+(k-1)*nx*ny,ib)
                   end do
                end do
             end do
          end do
       end do
     end if
    else
!
! X-axis (y-z div)
!
     if (iesize > 1) then
!$OMP DO
       do ib = 1, iesize
          call dfftw_execute_dft    (planx2,wk_afft_l(1,ib),wk_afft_l(1,ib))
!x     enddo
!x!$OMP DO
!x     do ib = 1, iesize
          do k = 1, nz
             do j = 1, ny
                do l = 1, fft_Y_x_dim
                   nis = nis_fft_Y_x(l)
                   nie = nie_fft_Y_x(l)
                   nnx = nie - nis + 1
                   iadd0 = nnx*ny*nz*(ib-1)*2
                   do i = nis, nie
                      iadd  = (i+(j-1)*nxp+(k-1)*nxp*ny)*2
                      iadd1 = iadd0+(i-nis+1+(j-1)*nnx+(k-1)*nnx*ny)*2
                      wk_send1(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send1(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
     else
       do ib = 1, iesize
!$OMP DO
          do k = 1, nz
             iadd = (1+(k-1)*nxp*ny)*2-1
             call dfftw_execute_dft(planx2_1d,wk_afft_l(iadd,ib),wk_afft_l(iadd,ib))
!o        end do
!o     enddo
!o     do ib = 1, iesize
!o        do k = 1, nz
             do j = 1, ny
                do l = 1, fft_Y_x_dim
                   nis = nis_fft_Y_x(l)
                   nie = nie_fft_Y_x(l)
                   nnx = nie - nis + 1
                   iadd0 = nnx*ny*nz*(ib-1)*2
                   do i = nis, nie
                      iadd  = (i+(j-1)*nxp+(k-1)*nxp*ny)*2
                      iadd1 = iadd0+(i-nis+1+(j-1)*nnx+(k-1)*nnx*ny)*2
                      wk_send1(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send1(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
     end if
    end if
#ifdef __FAPP__
    call fapp_stop('fftwf_inverse_fftw1',1,1)
#endif
#ifdef __TIMER_DO__
  call timer_end(254)
#endif

#ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_fft_xy_world)
  call timer_sta(137)
#endif

!$OMP MASTER
    call MPI_ALLTOALL(wk_send1, x2y_smax*kimg*iesize, mpi_double_precision,   &
   &                  wk_recv1, x2y_rmax*kimg*iesize, mpi_double_precision,   &
   &                                                  mpi_fft_xy_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Inverse_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 70, ierr)
     endif
!$OMP END MASTER
!$OMP BARRIER

#ifdef __TIMER_DO__
  call timer_sta(138)
#endif
#ifdef __TIMER_DO__
  call timer_end(138)
#endif

#else
!else ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_comm)
  call timer_sta(137)
#endif

!$OMP MASTER
    icnt_recv = 0
    lrank = myrank + 1
    if (lrank > (nmrank -1)) lrank = 0
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank -1)) lrank = 0
       if ((x2y_recv(1,lrank) /= 0)) then
          call mpi_irecv(wk_recv1(1,x2y_recv(2,lrank)), x2y_recv(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_r(icnt_recv), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Inverse_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 70, ierr)
           endif
          icnt_recv = icnt_recv + 1
       endif
    enddo

    icnt_send = 0
    lrank = myrank
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank -1)) lrank = 0
       if ((x2y_send(1,lrank) /= 0)) then
          call mpi_isend(wk_send1(1,x2y_send(2,lrank)), x2y_send(1,lrank)*kimg*iesize, &
                         mpi_double_precision, lrank, itag, mpi_comm, req_s(icnt_send), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Inverse_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 71, ierr)
           endif
          icnt_send = icnt_send + 1
       endif
    enddo

#ifdef __TIMER_DO__
  call timer_sta(138)
#endif

#ifdef __TIMER_DO__
  call timer_end(138)
#endif

    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Inverse_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 72, ierr)
     endif

    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Inverse_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 73, ierr)
     endif
!$OMP END MASTER
!$OMP BARRIER

#endif
!endif ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_end(137)
#endif

    nx = xyz_fft_y(2,1)-xyz_fft_y(1,1)+1
    ny = xyz_fft_y(2,2)-xyz_fft_y(1,2)+1
    nz = xyz_fft_y(2,3)-xyz_fft_y(1,3)+1
    nyp = ny
#ifdef __TIMER_DO__
  call timer_sta(255)
#endif
#ifdef __FAPP__
    call fapp_start('fftwf_inverse_fftw2',1,1)
#endif
    if (kimg == 1) then
     if (iesize > 1) then
!$OMP DO
       do ib = 1, iesize
          do k = 1, nz
             do l = 1, fft_X_y_dim
                nis = nis_fft_X_y(l)
                nie = nie_fft_X_y(l)
                nny = nie - nis + 1
                iadd0 = nx/2*nny*nz*(ib-1)*2
                do j = nis, nie
                   do i = 1, nx/2
                      iadd1 = iadd0 + (i+(j-nis)*nx/2+(k-1)*nx*nny/2)*2
                      iadd  = (i+(j-1)*nx/2+(k-1)*nx*nyp/2)*2
                      wk_afft_l(iadd-1,ib) = wk_recv1(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv1(iadd1  ,l)
                   end do
                end do
             end do
          end do
!x     end do
!
! Y-axis (z-x div)
!
!xOMP DO
!x     do ib = 1, iesize
          do k = 1, nz
             call dfftw_execute_dft(plany1,wk_afft_l((nx/2*nyp*(k-1)+1)*2-1,ib),wk_afft_l((nx/2*nyp*(k-1)+1)*2-1,ib))
          end do
!x     end do
!xOMP DO
!x     do ib = 1, iesize
          do k = 1, nz
             do l = 1, fft_Z_y_dim
                nis = nis_fft_Z_y(l)
                nie = nie_fft_Z_y(l)
                nny = nie - nis + 1
                iadd0 = nx/2*nny*nz*(ib-1)*2
                do j = nis, nie
                   do i = 1, nx/2
                      iadd  = (i+(j-1)*nx/2+(k-1)*nx*nyp/2)*2
                      iadd1 = iadd0 + (i+(j-nis)*nx/2+(k-1)*nx*nny/2)*2
                      wk_send2(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send2(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
     else
       do ib = 1, iesize
!$OMP DO
          do k = 1, nz
             do l = 1, fft_X_y_dim
                nis = nis_fft_X_y(l)
                nie = nie_fft_X_y(l)
                nny = nie - nis + 1
                iadd0 = nx/2*nny*nz*(ib-1)*2
                do j = nis, nie
                   do i = 1, nx/2
                      iadd1 = iadd0 + (i+(j-nis)*nx/2+(k-1)*nx*nny/2)*2
                      iadd  = (i+(j-1)*nx/2+(k-1)*nx*nyp/2)*2
                      wk_afft_l(iadd-1,ib) = wk_recv1(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv1(iadd1  ,l)
                   end do
                end do
             end do
          end do
!x     end do
!
! Y-axis (z-x div)
!
!x     do ib = 1, iesize
!$OMP DO
          do k = 1, nz
             call dfftw_execute_dft(plany1,wk_afft_l((nx/2*nyp*(k-1)+1)*2-1,ib),wk_afft_l((nx/2*nyp*(k-1)+1)*2-1,ib))
          end do
!x     end do
!x     do ib = 1, iesize
!$OMP DO
          do k = 1, nz
             do l = 1, fft_Z_y_dim
                nis = nis_fft_Z_y(l)
                nie = nie_fft_Z_y(l)
                nny = nie - nis + 1
                iadd0 = nx/2*nny*nz*(ib-1)*2
                do j = nis, nie
                   do i = 1, nx/2
                      iadd  = (i+(j-1)*nx/2+(k-1)*nx*nyp/2)*2
                      iadd1 = iadd0 + (i+(j-nis)*nx/2+(k-1)*nx*nny/2)*2
                      wk_send2(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send2(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
     end if
    else
     if (iesize > 1) then
!$OMP DO
       do ib = 1, iesize
          do k = 1, nz
             do l = 1, fft_X_y_dim
                nis = nis_fft_X_y(l)
                nie = nie_fft_X_y(l)
                nny = nie - nis + 1
                iadd0 = nx*nny*nz*(ib-1)*2
                do j = nis, nie
                   do i = 1, nx
                      iadd1 = iadd0+(i+(j-nis)*nx+(k-1)*nx*nny)*2
                      iadd  = (i+(j-1)*nx+(k-1)*nx*nyp)*2
                      wk_afft_l(iadd-1,ib) = wk_recv1(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv1(iadd1  ,l)
                   end do
                end do
             end do
          end do
!x     end do
!
! Y-axis (z-x div)
!
!xOMP DO
!x     do ib = 1, iesize
          do k = 1, nz
             call dfftw_execute_dft(plany2,wk_afft_l((nx*nyp*(k-1)+1)*2-1,ib),wk_afft_l((nx*nyp*(k-1)+1)*2-1,ib))
          end do
!x     end do
!xOMP DO
!x     do ib = 1, iesize
          do k = 1, nz
             do l = 1, fft_Z_y_dim
                nis = nis_fft_Z_y(l)
                nie = nie_fft_Z_y(l)
                nny = nie - nis + 1
                iadd0 = nx*nny*nz*(ib-1)*2
                do j = nis, nie
                   do i = 1, nx
                      iadd  = (i+(j-1)*nx+(k-1)*nx*nyp)*2
                      iadd1 = iadd0+(i+(j-nis)*nx+(k-1)*nx*nny)*2
                      wk_send2(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send2(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
     else
       do ib = 1, iesize
!$OMP DO
          do k = 1, nz
             do l = 1, fft_X_y_dim
                nis = nis_fft_X_y(l)
                nie = nie_fft_X_y(l)
                nny = nie - nis + 1
                iadd0 = nx*nny*nz*(ib-1)*2
                do j = nis, nie
                   do i = 1, nx
                      iadd1 = iadd0+(i+(j-nis)*nx+(k-1)*nx*nny)*2
                      iadd  = (i+(j-1)*nx+(k-1)*nx*nyp)*2
                      wk_afft_l(iadd-1,ib) = wk_recv1(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv1(iadd1  ,l)
                   end do
                end do
             end do
!o        end do
!o     end do
!
! Y-axis (z-x div)
!
!o     do ib = 1, iesize
!o        do k = 1, nz
             call dfftw_execute_dft(plany2,wk_afft_l((nx*nyp*(k-1)+1)*2-1,ib),wk_afft_l((nx*nyp*(k-1)+1)*2-1,ib))
!o        end do
!o     end do
!o     do ib = 1, iesize
!o        do k = 1, nz
             do l = 1, fft_Z_y_dim
                nis = nis_fft_Z_y(l)
                nie = nie_fft_Z_y(l)
                nny = nie - nis + 1
                iadd0 = nx*nny*nz*(ib-1)*2
                do j = nis, nie
                   do i = 1, nx
                      iadd  = (i+(j-1)*nx+(k-1)*nx*nyp)*2
                      iadd1 = iadd0+(i+(j-nis)*nx+(k-1)*nx*nny)*2
                      wk_send2(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send2(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
     end if
    end if
#ifdef __FAPP__
    call fapp_stop('fftwf_inverse_fftw2',1,1)
#endif
#ifdef __TIMER_DO__
  call timer_end(255)
#endif

#ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_fft_yz_world)
  call timer_sta(142)
#endif

!$OMP MASTER
    call MPI_ALLTOALL(wk_send2, y2z_smax*kimg*iesize, mpi_double_precision,   &
   &                  wk_recv2, y2z_rmax*kimg*iesize, mpi_double_precision,   &
   &                                                  mpi_fft_yz_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Inverse_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 76, ierr)
     endif
!$OMP END MASTER
!$OMP BARRIER

#ifdef __TIMER_DO__
  call timer_sta(143)
#endif
#ifdef __TIMER_DO__
  call timer_end(143)
#endif

#else
!else ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_comm)
  call timer_sta(142)
#endif

!$OMP MASTER
    icnt_recv = 0
    lrank = myrank + 1
    if (lrank > (nmrank-1)) lrank = 0
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank-1)) lrank = 0
       if ((z2y_recv(1,lrank) /= 0)) then
          call mpi_irecv(wk_recv2(1,z2y_recv(2,lrank)), z2y_recv(1,lrank)*kimg*iesize, &
                         mpi_double_precision, lrank, itag, mpi_comm, req_r(icnt_recv), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Inverse_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 76, ierr)
           endif
          icnt_recv = icnt_recv + 1
       endif
    enddo

    icnt_send = 0
    lrank = myrank
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank -1)) lrank = 0
       if ((z2y_send(1,lrank) /= 0)) then
          call mpi_isend(wk_send2(1,z2y_send(2,lrank)), z2y_send(1,lrank)*kimg*iesize, &
                         mpi_double_precision, lrank, itag, mpi_comm, req_s(icnt_send), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Inverse_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 77, ierr)
           endif
          icnt_send = icnt_send + 1
       endif
    enddo

#ifdef __TIMER_DO__
  call timer_sta(143)
#endif
#ifdef __TIMER_DO__
  call timer_end(143)
#endif

    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Inverse_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 78, ierr)
     endif

    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Inverse_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 79, ierr)
     endif
!$OMP END MASTER
!$OMP BARRIER

#endif
!endif ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_end(142)
#endif

    nx = xyz_fft_z(2,1)-xyz_fft_z(1,1)+1
    ny = xyz_fft_z(2,2)-xyz_fft_z(1,2)+1
    nz = xyz_fft_z(2,3)-xyz_fft_z(1,3)+1
    nzp = nz
#ifdef __TIMER_DO__
  call timer_sta(256)
#endif
#ifdef __FAPP__
    call fapp_start('fftwf_inverse_fftw3',1,1)
#endif
    if (kimg == 1) then
     if (iesize > 1) then
!$OMP DO
       do ib = 1, iesize
          do l = 1, fft_Y_z_dim
             nis = nis_fft_Y_z(l)
             nie = nie_fft_Y_z(l)
             nnz = nie - nis + 1
             iadd0 = nx/2*ny*nnz*(ib-1)*2
             do k = nis, nie
                do j = 1, ny
                   do i = 1, nx/2
                      iadd1 = iadd0 + (i+(j-1)*nx/2+(k-nis)*nx*ny/2)*2
                      iadd  = (i+(j-1)*nx/2+(k-1)*nx*ny/2)*2
                      wk_afft_l(iadd-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
          end do
!x     end do
!
! Z-axis (x-y div)
!
!xOMP DO
!x     do ib = 1, iesize
          call dfftw_execute_dft(planz1,wk_afft_l(1,ib),wk_afft_l(1,ib))
       end do
     else
       do ib = 1, iesize
!$OMP DO
          do l = 1, fft_Y_z_dim
             nis = nis_fft_Y_z(l)
             nie = nie_fft_Y_z(l)
             nnz = nie - nis + 1
             iadd0 = nx/2*ny*nnz*(ib-1)*2
             do k = nis, nie
                do j = 1, ny
                   do i = 1, nx/2
                      iadd1 = iadd0 + (i+(j-1)*nx/2+(k-nis)*nx*ny/2)*2
                      iadd  = (i+(j-1)*nx/2+(k-1)*nx*ny/2)*2
                      wk_afft_l(iadd-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end do
!
! Z-axis (x-y div)
!
       do ib = 1, iesize
!$OMP DO
          do j = 1, ny
             iadd = (1+(j-1)*nx/2)*2-1
             call dfftw_execute_dft(planz1_1d,wk_afft_l(iadd,ib),wk_afft_l(iadd,ib))
          end do
       end do
     end if
    else
     if (iesize > 1) then
!$OMP DO
       do ib = 1, iesize
          do l = 1, fft_Y_z_dim
             nis = nis_fft_Y_z(l)
             nie = nie_fft_Y_z(l)
             nnz = nie - nis + 1
             iadd0 = nx*ny*nnz*(ib-1)*2
             do k = nis, nie
                do j = 1, ny
                   do i = 1, nx
                      iadd1 = iadd0+(i+(j-1)*nx+(k-nis)*nx*ny)*2
                      iadd  = (i+(j-1)*nx+(k-1)*nx*ny)*2
                      wk_afft_l(iadd-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
          end do
!x     end do
!
! Z-axis (x-y div)
!
!xOMP DO
!x     do ib = 1, iesize
          call dfftw_execute_dft(planz2,wk_afft_l(1,ib),wk_afft_l(1,ib))
       end do
     else
       do ib = 1, iesize
!$OMP DO
          do l = 1, fft_Y_z_dim
             nis = nis_fft_Y_z(l)
             nie = nie_fft_Y_z(l)
             nnz = nie - nis + 1
             iadd0 = nx*ny*nnz*(ib-1)*2
             do k = nis, nie
                do j = 1, ny
                   do i = 1, nx
                      iadd1 = iadd0+(i+(j-1)*nx+(k-nis)*nx*ny)*2
                      iadd  = (i+(j-1)*nx+(k-1)*nx*ny)*2
                      wk_afft_l(iadd-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end do
!
! Z-axis (x-y div)
!
       do ib = 1, iesize
!$OMP DO
          do j = 1, ny
             iadd = (1+(j-1)*nx)*2-1
             call dfftw_execute_dft(planz2_1d,wk_afft_l(iadd,ib),wk_afft_l(iadd,ib))
          end do
       end do
     end if
    end if
#ifdef __FAPP__
    call fapp_stop('fftwf_inverse_fftw3',1,1)
#endif
#ifdef __TIMER_DO__
  call timer_end(256)
#endif

!$OMP END PARALLEL

#ifdef __TIMER_SUB__
  call timer_end(1498)
#endif

#ifdef __TIMER_SUB__
  call timer_end(1499)
#endif

1000 continue

#ifdef __TIMER_SUB__
    call timer_end(106)
#endif

    call tstatc0_end(id_sname)
  end subroutine m_FFT_Inverse_XYZ_3D

  subroutine m_FFT_Inverse_XYZ_3D_oo_place(nfout, wk_afft_l_input, wk_afft_l, wk_size, iesize)
    integer, intent(in) :: nfout, wk_size, iesize
    real(kind=DP), dimension(wk_size*kimg,iesize), intent(in) :: wk_afft_l_input
    real(kind=DP), dimension(wk_size*kimg,iesize), intent(out) :: wk_afft_l

    integer, parameter :: FFTW_MEASURE=0
    integer, parameter :: FFTW_ESTIMATE=64
    integer, parameter :: FFTW_RANK=1
    integer            :: FFTW_FLAG
    integer,dimension(1) :: NFFTW3, NEMBED, NEREAL
    integer :: nx, ny, nz, i, j, k, l, ri, ib, iadd, lrank
    integer :: irank, itag, icnt_send, icnt_recv
    integer,save, allocatable, dimension(:)   :: req_r, req_s
    integer,save, allocatable, dimension(:,:) :: sta_r, sta_s
    integer,save, allocatable, dimension(:)   :: wk_recvcnt, wk_sendcnt
    real(kind=DP),allocatable, dimension(:,:),save :: wk_recv1, wk_send1
    real(kind=DP),allocatable, dimension(:,:),save :: wk_recv2, wk_send2
    integer,      allocatable, dimension(:)   :: wk_mp_fft_y
    integer :: mpi_comm, myrank, nmrank
    integer(kind=8),save :: planx1 = 0, planx2 = 0
    integer(kind=8),save :: plany1 = 0, plany2 = 0
    integer(kind=8),save :: planz1 = 0, planz2 = 0
    integer(kind=8),save :: planx1_1d = 0, planx2_1d = 0
    integer(kind=8),save :: plany1_1d = 0, plany2_1d = 0
    integer(kind=8),save :: planz1_1d = 0, planz2_1d = 0
    integer(kind=8) :: plan
    integer :: max_x, max_y, max_z
    integer,save :: max_elm = 0, savesize = 0

    integer,save, allocatable, dimension(:,:) :: x2y_recv, x2y_send, y2z_recv, y2z_send
    integer,save :: x2y_rrank, x2y_srank, y2z_rrank, y2z_srank
    integer,save :: x2y_rmax, x2y_smax, x2y_srmax, y2z_rmax, y2z_smax, y2z_srmax

    integer ::  nnx, nny, nnz, iadd0, iadd1, nis, nie, jadd, kadd, nxp, nyp, nzp
    integer                  , dimension(3) :: isin,nsize
#ifdef FFT_ALLTOALL
    integer,dimension(2) :: max_send,max_recv
#endif
    integer :: id_sname = -1

    call tstatc0_begin('m_FFT_Inverse ',id_sname)

    if (kimg == 1) then
       FFTW_FLAG = -1      ! FFTW_FORWARD
    else
       FFTW_FLAG = +1      ! FFTW_BACKWARD
    endif

    mpi_comm = mpi_ke_world
    myrank = myrank_g
    nmrank = nrank_g
    itag = 10

    if (np_fft_x==0)then
       nx = 0
       nx = 0
       nz = 0
!      go to 1000
       go to 2000
    endif

    if (firstcall_inverse_xyz_3d_oop) then
       if(planx1    /= 0) call dfftw_destroy_plan(planx1)
       if(planx1_1d /= 0) call dfftw_destroy_plan(planx1_1d)
       if(planx2    /= 0) call dfftw_destroy_plan(planx2)
       if(planx2_1d /= 0) call dfftw_destroy_plan(planx2_1d)
       if(plany1    /= 0) call dfftw_destroy_plan(plany1)
       if(plany2    /= 0) call dfftw_destroy_plan(plany2)
       if(planz1    /= 0) call dfftw_destroy_plan(planz1)
       if(planz1_1d /= 0) call dfftw_destroy_plan(planz1_1d)
       if(planz2    /= 0) call dfftw_destroy_plan(planz2)
       if(planz2_1d /= 0) call dfftw_destroy_plan(planz2_1d)
       savesize = 0
       max_x = maxval(nel_fft_x(:))
       max_y = maxval(nel_fft_y(:))
       max_z = maxval(nel_fft_z(:))
       max_elm = max(max_x,max_y,max_z)
       if(allocated(req_r)) deallocate(req_r)
       if(allocated(req_s)) deallocate(req_s)
       if(allocated(sta_r)) deallocate(sta_r)
       if(allocated(sta_s)) deallocate(sta_s)
       if(allocated(wk_recvcnt)) deallocate(wk_recvcnt)
       if(allocated(wk_sendcnt)) deallocate(wk_sendcnt)
       allocate(req_r(0:nmrank-1), stat=ierr)
       allocate(req_s(0:nmrank-1), stat=ierr)
       allocate(sta_r(MPI_STATUS_SIZE,0:nmrank-1), stat=ierr)
       allocate(sta_s(MPI_STATUS_SIZE,0:nmrank-1), stat=ierr)
       allocate(wk_recvcnt(0:nmrank-1), stat=ierr)
       allocate(wk_sendcnt(0:nmrank-1), stat=ierr)
       if (ierr /= 0) then
          write(nfout,*)' m_FFT_Inverse_3D :  Not allocate '
          call flush(nfout)
          call mpi_abort(mpi_comm_world, 68, ierr)
       endif

       if(allocated(x2y_recv)) deallocate(x2y_recv)
       if(allocated(x2y_send)) deallocate(x2y_send)
       allocate(x2y_recv(2,0:nmrank-1))
       allocate(x2y_send(2,0:nmrank-1))
       x2y_send = 0
       x2y_recv = 0
       wk_recvcnt = 0
       wk_sendcnt = 0
       do i = 1, nel_fft_y(myrank)
          irank = map_fft_x(mp_fft_y(i)) - 1
          wk_recvcnt(irank) = wk_recvcnt(irank) + 1
       enddo
       do i = 1, nel_fft_x(myrank)
          irank = map_fft_y(mp_fft_x(i)) - 1
          wk_sendcnt(irank) = wk_sendcnt(irank) + 1
       enddo
       k = 0
       do i = 0, nmrank - 1
          if(wk_recvcnt(i) /= 0) then
             k = k + 1
             x2y_recv(1,i) = wk_recvcnt(i)
             x2y_recv(2,i) = k
          endif
       enddo
       x2y_rrank = k
       k = 0
       do i = 0, nmrank - 1
          if(wk_sendcnt(i) /= 0) then
             k = k + 1
             x2y_send(1,i) = wk_sendcnt(i)
             x2y_send(2,i) = k
          endif
       enddo
       x2y_srank = k
       x2y_rmax = maxval(wk_recvcnt)
       x2y_smax = maxval(wk_sendcnt)
       x2y_srmax = max(x2y_rmax,x2y_smax)

       if(allocated(y2z_recv)) deallocate(y2z_recv)
       if(allocated(y2z_send)) deallocate(y2z_send)
       allocate(y2z_recv(2,0:nmrank-1))
       allocate(y2z_send(2,0:nmrank-1))
       y2z_send = 0
       y2z_recv = 0
       wk_recvcnt = 0
       wk_sendcnt = 0

       nx = xyz_fft_y(2,1)-xyz_fft_y(1,1)+1
       ny = xyz_fft_y(2,2)-xyz_fft_y(1,2)+1
       nz = xyz_fft_y(2,3)-xyz_fft_y(1,3)+1
       allocate(wk_mp_fft_y(nx*ny*nz) ,stat=ierr)
       do k = 0, nz-1
          do j = 0, ny-1
             do i = 0, nx-1
                wk_mp_fft_y(i+j*nx+k*nx*ny+1) = mp_fft_y(i+k*nx+j*nx*nz+1)
             enddo
          enddo
       enddo

       do i = 1, nel_fft_z(myrank)
          irank = map_fft_y(mp_fft_z(i)) - 1
          wk_recvcnt(irank) = wk_recvcnt(irank) + 1
       enddo
       do i = 1, nel_fft_y(myrank)
          irank = map_fft_z(wk_mp_fft_y(i)) - 1
          wk_sendcnt(irank) = wk_sendcnt(irank) + 1
       enddo
       k = 0
       do i = 0, nmrank - 1
          if(wk_recvcnt(i) /= 0) then
             k = k + 1
             y2z_recv(1,i) = wk_recvcnt(i)
             y2z_recv(2,i) = k
          endif
       enddo
       y2z_rrank = k
       k = 0
      do i = 0, nmrank - 1
          if(wk_sendcnt(i) /= 0) then
             k = k + 1
             y2z_send(1,i) = wk_sendcnt(i)
             y2z_send(2,i) = k
          endif
       enddo
       y2z_srank = k
       y2z_rmax = maxval(wk_recvcnt)
       y2z_smax = maxval(wk_sendcnt)
       y2z_srmax = max(y2z_rmax,y2z_smax)

       deallocate(wk_mp_fft_y)

#ifdef FFT_ALLTOALL
       max_send(1) = x2y_rmax
       max_send(2) = x2y_smax
       call mpi_allreduce(max_send,max_recv,2,mpi_integer,mpi_max,mpi_fft_xy_world,ierr)
       x2y_rmax = max_recv(1)
       x2y_smax = max_recv(2)
       max_send(1) = y2z_rmax
       max_send(2) = y2z_smax
       call mpi_allreduce(max_send,max_recv,2,mpi_integer,mpi_max,mpi_fft_yz_world,ierr)
       y2z_rmax = max_recv(1)
       y2z_smax = max_recv(2)
#endif

       if (ipri > 1) then
          write(nfout,'("m_FFT_Inverse_XYZ_3D   --   myrank_g=",i4)') myrank_g
          write(nfout,'("x2y_send")')
          write(nfout,'(10(i8,", "))') (x2y_send(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (x2y_send(2,i),i=0,nmrank-1)
          write(nfout,'("x2y_recv")')
          write(nfout,'(10(i8,", "))') (x2y_recv(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (x2y_recv(2,i),i=0,nmrank-1)
          write(nfout,'("y2z_send")')
          write(nfout,'(10(i8,", "))') (y2z_send(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (y2z_send(2,i),i=0,nmrank-1)
          write(nfout,'("y2z_recv")')
          write(nfout,'(10(i8,", "))') (y2z_recv(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (y2z_recv(2,i),i=0,nmrank-1)
          call flush(nfout)
       endif

       nx = xyz_fft_x(2,1)-xyz_fft_x(1,1)+1
       ny = xyz_fft_x(2,2)-xyz_fft_x(1,2)+1
       nz = xyz_fft_x(2,3)-xyz_fft_x(1,3)+1
       if(kimg==1) then
          call dfftw_plan_many_dft_r2c(planx1   ,        1, nx-2, ny*nz, &
         &                             wk_afft_l_input,  nx-2,    1,    nx, &
         &                             wk_afft_l, (nx-2)/2,    1,  nx/2, &
         &                             FFTW_MEASURE )
          call dfftw_plan_many_dft_r2c(planx1_1d,        1, nx-2, ny   , &
         &                             wk_afft_l_input,     nx-2,    1,    nx, &
         &                             wk_afft_l, (nx-2)/2,    1,  nx/2, &
         &                             FFTW_MEASURE )
!?        call dfftw_plan_dft_r2c_1d(planx1_1d, nx-2, wk_afft_l, wk_afft_l, FFTW_FLAG, FFTW_MEASURE)
       else
          call dfftw_plan_many_dft    (planx2   ,  1, nx, ny*nz, &
         &                             wk_afft_l_input, nx,  1,    nx, &
         &                             wk_afft_l, nx,  1,    nx, &
         &                             FFTW_FLAG, FFTW_MEASURE )
          call dfftw_plan_many_dft    (planx2_1d,  1, nx, ny   , &
         &                             wk_afft_l_input, nx,  1,    nx, &
         &                             wk_afft_l, nx,  1,    nx, &
         &                             FFTW_FLAG, FFTW_MEASURE )
       endif

       nx = xyz_fft_y(2,1)-xyz_fft_y(1,1)+1
       ny = xyz_fft_y(2,2)-xyz_fft_y(1,2)+1
       nz = xyz_fft_y(2,3)-xyz_fft_y(1,3)+1
       if(kimg==1) then
          call dfftw_plan_many_dft    (plany1   ,  1,   ny, nx/2, &
         &                             wk_afft_l, ny, nx/2,    1, &
         &                             wk_afft_l, ny, nx/2,    1, &
         &                             FFTW_FLAG, FFTW_MEASURE )
       else
          call dfftw_plan_many_dft    (plany2   ,  1, ny, nx, &
         &                             wk_afft_l, ny, nx,  1, &
         &                             wk_afft_l, ny, nx,  1, &
         &                             FFTW_FLAG, FFTW_MEASURE )
       end if

       nx = xyz_fft_z(2,1)-xyz_fft_z(1,1)+1
       ny = xyz_fft_z(2,2)-xyz_fft_z(1,2)+1
       nz = xyz_fft_z(2,3)-xyz_fft_z(1,3)+1
       if(kimg==1) then
          call dfftw_plan_many_dft    (planz1   ,  1,      nz, nx*ny/2, &
         &                             wk_afft_l, nz, nx*ny/2,       1, &
         &                             wk_afft_l, nz, nx*ny/2,       1, &
         &                             FFTW_FLAG, FFTW_MEASURE )
          call dfftw_plan_many_dft    (planz1_1d,  1,      nz, nx/2   , &
         &                             wk_afft_l, nz, nx*ny/2,       1, &
         &                             wk_afft_l, nz, nx*ny/2,       1, &
         &                             FFTW_FLAG, FFTW_MEASURE )
       else
          call dfftw_plan_many_dft    (planz2   ,  1,    nz, nx*ny, &
         &                             wk_afft_l, nz, nx*ny,     1, &
         &                             wk_afft_l, nz, nx*ny,     1, &
         &                             FFTW_FLAG, FFTW_MEASURE )
          call dfftw_plan_many_dft    (planz2_1d,  1, nz,    nx   , &
         &                             wk_afft_l, nz, nx*ny,     1, &
         &                             wk_afft_l, nz, nx*ny,     1, &
         &                             FFTW_FLAG, FFTW_MEASURE )
       endif

       firstcall_inverse_xyz_3d_oop = .false.
!       go to 1000
    endif

#ifdef FFT_ALLTOALL
    if (iesize /= savesize) then
#else
    if (iesize > savesize) then
#endif
       if (allocated(wk_recv1)) deallocate(wk_recv1)
       if (allocated(wk_send1)) deallocate(wk_send1)
       if (allocated(wk_recv2)) deallocate(wk_recv2)
       if (allocated(wk_send2)) deallocate(wk_send2)
       allocate(wk_recv1(x2y_rmax*kimg*iesize,x2y_rrank), stat=ierr)
       allocate(wk_send1(x2y_smax*kimg*iesize,x2y_srank), stat=ierr)
       allocate(wk_recv2(y2z_rmax*kimg*iesize,y2z_rrank), stat=ierr)
       allocate(wk_send2(y2z_smax*kimg*iesize,y2z_srank), stat=ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_FFT_Inverse_3D :  Not allocate '
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 69, ierr)
        endif
       savesize = iesize
    end if

 2000 continue
    if (np_fft_x==0)then
       nx = 0
       nx = 0
       nz = 0
       go to 1000
    endif

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP          SHARED(     nfout,wk_afft_l_input,wk_afft_l,iesize,req_r,req_s,sta_r,sta_s,    &
!$OMP                 wk_recv1,wk_send1,wk_recv2,wk_send2,                    &
!$OMP                 planx1,planx2,planz1,planz2,plany1,plany2,              &
!$OMP                 planx1_1d,planz1_1d,planx2_1d,planz2_1d,                &
!$OMP                 x2y_recv, x2y_send, y2z_recv, y2z_send,itag,            &
!$OMP                 xyz_fft_y,xyz_fft_z,xyz_fft_x,fft_Z_x_dim,nis_fft_Z_x,  &
!$OMP                 nie_fft_Z_x,fft_X_z_dim,nis_fft_X_z,nie_fft_X_z,        &
!$OMP                 fft_Y_z_dim,nis_fft_Y_z,nie_fft_Y_z,fft_Z_y_dim,        &
!$OMP                 nie_fft_Z_y,nis_fft_Z_y,                                &
!$OMP                 nie_fft_X_y,nis_fft_X_y,fft_X_y_dim,                    &
!$OMP                 nie_fft_Y_x,nis_fft_Y_x,fft_Y_x_dim,                    &
#ifdef FFT_ALLTOALL
!$OMP                 mpi_fft_xy_world, x2y_rmax, x2y_smax,                   &
!$OMP                 mpi_fft_yz_world, y2z_rmax, y2z_smax,                   &
#endif
!$OMP                 mpi_comm,myrank,nmrank,kimg,ierr                    )   &
!$OMP          PRIVATE(nx,ny,nz,i,j,k,l,ri,ib,iadd,lrank,icnt_send,icnt_recv, &
!$OMP                  nxp,nyp,nzp,                                           &
!$OMP                  irank,nnx,nny,nnz,iadd0,iadd1,nis,nie,jadd,kadd,plan   )

    nx = xyz_fft_x(2,1)-xyz_fft_x(1,1)+1
    ny = xyz_fft_x(2,2)-xyz_fft_x(1,2)+1
    nz = xyz_fft_x(2,3)-xyz_fft_x(1,3)+1
    nxp = nx
#ifdef __TIMER_DO__
  call timer_sta(254)
#endif
#ifdef __FAPP__
    call fapp_start('fftwf_inverse_fftw1',1,1)
#endif
    if (kimg == 1) then
!
! X-axis (y-z div)
!
     if (iesize > 1) then
!$OMP DO
       do ib = 1, iesize
          call dfftw_execute_dft_r2c(planx1,wk_afft_l_input(1,ib),wk_afft_l(1,ib))
!o     enddo
!oOMP DO
!o     do ib = 1, iesize
          do k = 1, nz
             do j = 1, ny
                do l = 1, fft_Y_x_dim
                   nis = nis_fft_Y_x(l)
                   nie = nie_fft_Y_x(l)
                   nnx = nie - nis + 1
                   iadd0 = nnx*ny*nz*(ib-1)
                   do i = nis, nie
                      wk_send1(iadd0+(i-nis+1+(j-1)*nnx+(k-1)*nnx*ny),l) = wk_afft_l(i+(j-1)*nx+(k-1)*nx*ny,ib)
                   end do
                end do
             end do
          end do
       end do
     else
       do ib = 1, iesize
!$OMP DO
          do k = 1, nz
!?           do j = 1, ny
!?              iadd = (1+(j-1)*nx/2+(k-1)*nx*ny/2)*2-1
!?              call dfftw_execute_dft_r2c(planx1_1d,wk_afft_l(iadd,ib),wk_afft_l(iadd,ib))
!?           end do
             iadd = (1+(k-1)*nx*ny/2)*2-1
             call dfftw_execute_dft_r2c(planx1_1d,wk_afft_l_input(iadd,ib),wk_afft_l(iadd,ib))
!o        end do
!o     end do
!o     do ib = 1, iesize
!oOMP DO
!o        do k = 1, nz
             do j = 1, ny
                do l = 1, fft_Y_x_dim
                   nis = nis_fft_Y_x(l)
                   nie = nie_fft_Y_x(l)
                   nnx = nie - nis + 1
                   iadd0 = nnx*ny*nz*(ib-1)
                   do i = nis, nie
                      wk_send1(iadd0+(i-nis+1+(j-1)*nnx+(k-1)*nnx*ny),l) = wk_afft_l(i+(j-1)*nx+(k-1)*nx*ny,ib)
                   end do
                end do
             end do
          end do
       end do
     end if
    else
!
! X-axis (y-z div)
!
     if (iesize > 1) then
!$OMP DO
       do ib = 1, iesize
          call dfftw_execute_dft    (planx2,wk_afft_l_input(1,ib),wk_afft_l(1,ib))
!x     enddo
!x!$OMP DO
!x     do ib = 1, iesize
          do k = 1, nz
             do j = 1, ny
                do l = 1, fft_Y_x_dim
                   nis = nis_fft_Y_x(l)
                   nie = nie_fft_Y_x(l)
                   nnx = nie - nis + 1
                   iadd0 = nnx*ny*nz*(ib-1)*2
                   do i = nis, nie
                      iadd  = (i+(j-1)*nxp+(k-1)*nxp*ny)*2
                      iadd1 = iadd0+(i-nis+1+(j-1)*nnx+(k-1)*nnx*ny)*2
                      wk_send1(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send1(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
     else
       do ib = 1, iesize
!$OMP DO
          do k = 1, nz
             iadd = (1+(k-1)*nxp*ny)*2-1
             call dfftw_execute_dft(planx2_1d,wk_afft_l_input(iadd,ib),wk_afft_l(iadd,ib))
!o        end do
!o     enddo
!o     do ib = 1, iesize
!o        do k = 1, nz
             do j = 1, ny
                do l = 1, fft_Y_x_dim
                   nis = nis_fft_Y_x(l)
                   nie = nie_fft_Y_x(l)
                   nnx = nie - nis + 1
                   iadd0 = nnx*ny*nz*(ib-1)*2
                   do i = nis, nie
                      iadd  = (i+(j-1)*nxp+(k-1)*nxp*ny)*2
                      iadd1 = iadd0+(i-nis+1+(j-1)*nnx+(k-1)*nnx*ny)*2
                      wk_send1(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send1(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
     end if
    end if

#ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_fft_xy_world)
  call timer_sta(137)
#endif

!$OMP MASTER
    call MPI_ALLTOALL(wk_send1, x2y_smax*kimg*iesize, mpi_double_precision,   &
   &                  wk_recv1, x2y_rmax*kimg*iesize, mpi_double_precision,   &
   &                                                  mpi_fft_xy_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Inverse_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 70, ierr)
     endif
!$OMP END MASTER
!$OMP BARRIER

#else
!else ifdef FFT_ALLTOALL

!$OMP MASTER
    icnt_recv = 0
    lrank = myrank + 1
    if (lrank > (nmrank -1)) lrank = 0
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank -1)) lrank = 0
       if ((x2y_recv(1,lrank) /= 0)) then
          call mpi_irecv(wk_recv1(1,x2y_recv(2,lrank)), x2y_recv(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_r(icnt_recv), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Inverse_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 70, ierr)
           endif
          icnt_recv = icnt_recv + 1
       endif
    enddo

    icnt_send = 0
    lrank = myrank
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank -1)) lrank = 0
       if ((x2y_send(1,lrank) /= 0)) then
          call mpi_isend(wk_send1(1,x2y_send(2,lrank)), x2y_send(1,lrank)*kimg*iesize, &
                         mpi_double_precision, lrank, itag, mpi_comm, req_s(icnt_send), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Inverse_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 71, ierr)
           endif
          icnt_send = icnt_send + 1
       endif
    enddo

    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Inverse_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 72, ierr)
     endif

    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Inverse_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 73, ierr)
     endif
!$OMP END MASTER
!$OMP BARRIER

#endif
!endif ifdef FFT_ALLTOALL

    nx = xyz_fft_y(2,1)-xyz_fft_y(1,1)+1
    ny = xyz_fft_y(2,2)-xyz_fft_y(1,2)+1
    nz = xyz_fft_y(2,3)-xyz_fft_y(1,3)+1
    nyp = ny

    if (kimg == 1) then
     if (iesize > 1) then
!$OMP DO
       do ib = 1, iesize
          do k = 1, nz
             do l = 1, fft_X_y_dim
                nis = nis_fft_X_y(l)
                nie = nie_fft_X_y(l)
                nny = nie - nis + 1
                iadd0 = nx/2*nny*nz*(ib-1)*2
                do j = nis, nie
                   do i = 1, nx/2
                      iadd1 = iadd0 + (i+(j-nis)*nx/2+(k-1)*nx*nny/2)*2
                      iadd  = (i+(j-1)*nx/2+(k-1)*nx*nyp/2)*2
                      wk_afft_l(iadd-1,ib) = wk_recv1(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv1(iadd1  ,l)
                   end do
                end do
             end do
          end do
!x     end do
!
! Y-axis (z-x div)
!
!xOMP DO
!x     do ib = 1, iesize
          do k = 1, nz
             call dfftw_execute_dft(plany1,wk_afft_l((nx/2*nyp*(k-1)+1)*2-1,ib),wk_afft_l((nx/2*nyp*(k-1)+1)*2-1,ib))
          end do
!x     end do
!xOMP DO
!x     do ib = 1, iesize
          do k = 1, nz
             do l = 1, fft_Z_y_dim
                nis = nis_fft_Z_y(l)
                nie = nie_fft_Z_y(l)
                nny = nie - nis + 1
                iadd0 = nx/2*nny*nz*(ib-1)*2
                do j = nis, nie
                   do i = 1, nx/2
                      iadd  = (i+(j-1)*nx/2+(k-1)*nx*nyp/2)*2
                      iadd1 = iadd0 + (i+(j-nis)*nx/2+(k-1)*nx*nny/2)*2
                      wk_send2(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send2(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
     else
       do ib = 1, iesize
!$OMP DO
          do k = 1, nz
             do l = 1, fft_X_y_dim
                nis = nis_fft_X_y(l)
                nie = nie_fft_X_y(l)
                nny = nie - nis + 1
                iadd0 = nx/2*nny*nz*(ib-1)*2
                do j = nis, nie
                   do i = 1, nx/2
                      iadd1 = iadd0 + (i+(j-nis)*nx/2+(k-1)*nx*nny/2)*2
                      iadd  = (i+(j-1)*nx/2+(k-1)*nx*nyp/2)*2
                      wk_afft_l(iadd-1,ib) = wk_recv1(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv1(iadd1  ,l)
                   end do
                end do
             end do
          end do
!x     end do
!
! Y-axis (z-x div)
!
!x     do ib = 1, iesize
!$OMP DO
          do k = 1, nz
             call dfftw_execute_dft(plany1,wk_afft_l((nx/2*nyp*(k-1)+1)*2-1,ib),wk_afft_l((nx/2*nyp*(k-1)+1)*2-1,ib))
          end do
!x     end do
!x     do ib = 1, iesize
!$OMP DO
          do k = 1, nz
             do l = 1, fft_Z_y_dim
                nis = nis_fft_Z_y(l)
                nie = nie_fft_Z_y(l)
                nny = nie - nis + 1
                iadd0 = nx/2*nny*nz*(ib-1)*2
                do j = nis, nie
                   do i = 1, nx/2
                      iadd  = (i+(j-1)*nx/2+(k-1)*nx*nyp/2)*2
                      iadd1 = iadd0 + (i+(j-nis)*nx/2+(k-1)*nx*nny/2)*2
                      wk_send2(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send2(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
     end if
    else
     if (iesize > 1) then
!$OMP DO
       do ib = 1, iesize
          do k = 1, nz
             do l = 1, fft_X_y_dim
                nis = nis_fft_X_y(l)
                nie = nie_fft_X_y(l)
                nny = nie - nis + 1
                iadd0 = nx*nny*nz*(ib-1)*2
                do j = nis, nie
                   do i = 1, nx
                      iadd1 = iadd0+(i+(j-nis)*nx+(k-1)*nx*nny)*2
                      iadd  = (i+(j-1)*nx+(k-1)*nx*nyp)*2
                      wk_afft_l(iadd-1,ib) = wk_recv1(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv1(iadd1  ,l)
                   end do
                end do
             end do
          end do
!x     end do
!
! Y-axis (z-x div)
!
!xOMP DO
!x     do ib = 1, iesize
          do k = 1, nz
             call dfftw_execute_dft(plany2,wk_afft_l((nx*nyp*(k-1)+1)*2-1,ib),wk_afft_l((nx*nyp*(k-1)+1)*2-1,ib))
          end do
!x     end do
!xOMP DO
!x     do ib = 1, iesize
          do k = 1, nz
             do l = 1, fft_Z_y_dim
                nis = nis_fft_Z_y(l)
                nie = nie_fft_Z_y(l)
                nny = nie - nis + 1
                iadd0 = nx*nny*nz*(ib-1)*2
                do j = nis, nie
                   do i = 1, nx
                      iadd  = (i+(j-1)*nx+(k-1)*nx*nyp)*2
                      iadd1 = iadd0+(i+(j-nis)*nx+(k-1)*nx*nny)*2
                      wk_send2(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send2(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
     else
       do ib = 1, iesize
!$OMP DO
          do k = 1, nz
             do l = 1, fft_X_y_dim
                nis = nis_fft_X_y(l)
                nie = nie_fft_X_y(l)
                nny = nie - nis + 1
                iadd0 = nx*nny*nz*(ib-1)*2
                do j = nis, nie
                   do i = 1, nx
                      iadd1 = iadd0+(i+(j-nis)*nx+(k-1)*nx*nny)*2
                      iadd  = (i+(j-1)*nx+(k-1)*nx*nyp)*2
                      wk_afft_l(iadd-1,ib) = wk_recv1(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv1(iadd1  ,l)
                   end do
                end do
             end do
!o        end do
!o     end do
!
! Y-axis (z-x div)
!
!o     do ib = 1, iesize
!o        do k = 1, nz
             call dfftw_execute_dft(plany2,wk_afft_l((nx*nyp*(k-1)+1)*2-1,ib),wk_afft_l((nx*nyp*(k-1)+1)*2-1,ib))
!o        end do
!o     end do
!o     do ib = 1, iesize
!o        do k = 1, nz
             do l = 1, fft_Z_y_dim
                nis = nis_fft_Z_y(l)
                nie = nie_fft_Z_y(l)
                nny = nie - nis + 1
                iadd0 = nx*nny*nz*(ib-1)*2
                do j = nis, nie
                   do i = 1, nx
                      iadd  = (i+(j-1)*nx+(k-1)*nx*nyp)*2
                      iadd1 = iadd0+(i+(j-nis)*nx+(k-1)*nx*nny)*2
                      wk_send2(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send2(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
     end if
    end if

#ifdef FFT_ALLTOALL

!$OMP MASTER
    call MPI_ALLTOALL(wk_send2, y2z_smax*kimg*iesize, mpi_double_precision,   &
   &                  wk_recv2, y2z_rmax*kimg*iesize, mpi_double_precision,   &
   &                                                  mpi_fft_yz_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Inverse_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 76, ierr)
     endif
!$OMP END MASTER
!$OMP BARRIER

#else
!else ifdef FFT_ALLTOALL

!$OMP MASTER
    icnt_recv = 0
    lrank = myrank + 1
    if (lrank > (nmrank-1)) lrank = 0
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank-1)) lrank = 0
       if ((z2y_recv(1,lrank) /= 0)) then
          call mpi_irecv(wk_recv2(1,z2y_recv(2,lrank)), z2y_recv(1,lrank)*kimg*iesize, &
                         mpi_double_precision, lrank, itag, mpi_comm, req_r(icnt_recv), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Inverse_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 76, ierr)
           endif
          icnt_recv = icnt_recv + 1
       endif
    enddo

    icnt_send = 0
    lrank = myrank
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank -1)) lrank = 0
       if ((z2y_send(1,lrank) /= 0)) then
          call mpi_isend(wk_send2(1,z2y_send(2,lrank)), z2y_send(1,lrank)*kimg*iesize, &
                         mpi_double_precision, lrank, itag, mpi_comm, req_s(icnt_send), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_Inverse_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 77, ierr)
           endif
          icnt_send = icnt_send + 1
       endif
    enddo

    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Inverse_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 78, ierr)
     endif

    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_Inverse_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 79, ierr)
     endif
!$OMP END MASTER
!$OMP BARRIER

#endif
!endif ifdef FFT_ALLTOALL

    nx = xyz_fft_z(2,1)-xyz_fft_z(1,1)+1
    ny = xyz_fft_z(2,2)-xyz_fft_z(1,2)+1
    nz = xyz_fft_z(2,3)-xyz_fft_z(1,3)+1
    nzp = nz

    if (kimg == 1) then
     if (iesize > 1) then
!$OMP DO
       do ib = 1, iesize
          do l = 1, fft_Y_z_dim
             nis = nis_fft_Y_z(l)
             nie = nie_fft_Y_z(l)
             nnz = nie - nis + 1
             iadd0 = nx/2*ny*nnz*(ib-1)*2
             do k = nis, nie
                do j = 1, ny
                   do i = 1, nx/2
                      iadd1 = iadd0 + (i+(j-1)*nx/2+(k-nis)*nx*ny/2)*2
                      iadd  = (i+(j-1)*nx/2+(k-1)*nx*ny/2)*2
                      wk_afft_l(iadd-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
          end do
!x     end do
!
! Z-axis (x-y div)
!
!xOMP DO
!x     do ib = 1, iesize
          call dfftw_execute_dft(planz1,wk_afft_l(1,ib),wk_afft_l(1,ib))
       end do
     else
       do ib = 1, iesize
!$OMP DO
          do l = 1, fft_Y_z_dim
             nis = nis_fft_Y_z(l)
             nie = nie_fft_Y_z(l)
             nnz = nie - nis + 1
             iadd0 = nx/2*ny*nnz*(ib-1)*2
             do k = nis, nie
                do j = 1, ny
                   do i = 1, nx/2
                      iadd1 = iadd0 + (i+(j-1)*nx/2+(k-nis)*nx*ny/2)*2
                      iadd  = (i+(j-1)*nx/2+(k-1)*nx*ny/2)*2
                      wk_afft_l(iadd-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end do
!
! Z-axis (x-y div)
!
       do ib = 1, iesize
!$OMP DO
          do j = 1, ny
             iadd = (1+(j-1)*nx/2)*2-1
             call dfftw_execute_dft(planz1_1d,wk_afft_l(iadd,ib),wk_afft_l(iadd,ib))
          end do
       end do
     end if
    else
     if (iesize > 1) then
!$OMP DO
       do ib = 1, iesize
          do l = 1, fft_Y_z_dim
             nis = nis_fft_Y_z(l)
             nie = nie_fft_Y_z(l)
             nnz = nie - nis + 1
             iadd0 = nx*ny*nnz*(ib-1)*2
             do k = nis, nie
                do j = 1, ny
                   do i = 1, nx
                      iadd1 = iadd0+(i+(j-1)*nx+(k-nis)*nx*ny)*2
                      iadd  = (i+(j-1)*nx+(k-1)*nx*ny)*2
                      wk_afft_l(iadd-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
          end do
!x     end do
!
! Z-axis (x-y div)
!
!xOMP DO
!x     do ib = 1, iesize
          call dfftw_execute_dft(planz2,wk_afft_l(1,ib),wk_afft_l(1,ib))
       end do
     else
       do ib = 1, iesize
!$OMP DO
          do l = 1, fft_Y_z_dim
             nis = nis_fft_Y_z(l)
             nie = nie_fft_Y_z(l)
             nnz = nie - nis + 1
             iadd0 = nx*ny*nnz*(ib-1)*2
             do k = nis, nie
                do j = 1, ny
                   do i = 1, nx
                      iadd1 = iadd0+(i+(j-1)*nx+(k-nis)*nx*ny)*2
                      iadd  = (i+(j-1)*nx+(k-1)*nx*ny)*2
                      wk_afft_l(iadd-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end do
!
! Z-axis (x-y div)
!
       do ib = 1, iesize
!$OMP DO
          do j = 1, ny
             iadd = (1+(j-1)*nx)*2-1
             call dfftw_execute_dft(planz2_1d,wk_afft_l(iadd,ib),wk_afft_l(iadd,ib))
          end do
       end do
     end if
    end if

!$OMP END PARALLEL

1000 continue

    call tstatc0_end(id_sname)
  end subroutine m_FFT_Inverse_XYZ_3D_oo_place
!------------------------------------------------------------------------------
#endif
!endif #ifdef FFTW_NOSTRIDE
#endif
!endif #ifdef FFT_USE_SSL2

#ifdef MPI_FFTW

  subroutine m_FFT_init_mpifftw()
    integer(C_INTPTR_T) :: i, j, alloc_local, l, m, n, ln, ln_off, lm, lm_off, ll
    type(C_PTR) :: cdata
    l  = fft_box_size_WF(1,0)
    m  = fft_box_size_WF(2,0)
    n  = fft_box_size_WF(3,0)
    ll = fft_box_size_WF(1,1)
    if(kimg==2) then
      alloc_local = fftw_mpi_local_size_3d(n,m,l,mpi_ke_world,ln,ln_off)
      cdata = fftw_alloc_complex(alloc_local)
      call c_f_pointer(cdata, bfft_mpifftw, [l,m,ln])
      if (n /= m) then
        alloc_local = fftw_mpi_local_size_3d(m,n,l,mpi_ke_world,lm,lm_off)
        cdata = fftw_alloc_complex(alloc_local)
        call c_f_pointer(cdata, afft_mpifftw, [l,n,lm])
      else
        afft_mpifftw => bfft_mpifftw
      endif

      if(.not.distfftw_inverse_wf_flg) then
        plan_distfftw_inverse_wf = fftw_mpi_plan_dft_3d(n, m, l, bfft_mpifftw, afft_mpifftw, mpi_ke_world, &
                                                        FFTW_BACKWARD, FFTW_MEASURE+FFTW_MPI_TRANSPOSED_OUT)
      endif
      if(.not.distfftw_direct_wf_flg) then
        plan_distfftw_direct_wf  = fftw_mpi_plan_dft_3d(m, n, l, afft_mpifftw, bfft_mpifftw, mpi_ke_world, &
                                                        FFTW_FORWARD,  FFTW_MEASURE+FFTW_MPI_TRANSPOSED_OUT)
      endif
    else
      alloc_local = fftw_mpi_local_size_3d(m,n,l/2,mpi_ke_world,lm,lm_off)
      cdata = fftw_alloc_complex(alloc_local)
      call c_f_pointer(cdata, afft_mpifftw_kimg1, [l/2,n,lm])
      alloc_local = fftw_mpi_local_size_3d(n,m,l,mpi_ke_world,ln,ln_off)
      cdata = fftw_alloc_real(alloc_local)
      call c_f_pointer(cdata, bfft_mpifftw_kimg1, [l,m,ln])

      if(.not.distfftw_inverse_wf_flg) then
        plan_distfftw_inverse_wf = fftw_mpi_plan_dft_r2c_3d(n, m, ll, bfft_mpifftw_kimg1, afft_mpifftw_kimg1, mpi_ke_world, &
                                                            FFTW_MEASURE+FFTW_MPI_TRANSPOSED_OUT)
      endif
      if(.not.distfftw_direct_wf_flg) then
        plan_distfftw_direct_wf  = fftw_mpi_plan_dft_c2r_3d(m, n, ll, afft_mpifftw_kimg1, bfft_mpifftw_kimg1, mpi_ke_world, &
                                                            FFTW_MEASURE+FFTW_MPI_TRANSPOSED_OUT)
      endif
    endif

    distfftw_inverse_wf_flg = .true.
    distfftw_direct_wf_flg  = .true.
  end subroutine m_FFT_init_mpifftw

  subroutine m_FFT_Inverse_MPI_FFTW(nfout)
    integer, intent(in) :: nfout
    integer :: id_sname=-1
    call tstatc0_begin('m_FFT_Inverse_MPI_FFTW ',id_sname)
    if(kimg==2) then
      call fftw_mpi_execute_dft(plan_distfftw_inverse_wf,bfft_mpifftw,afft_mpifftw)
    else
      call fftw_mpi_execute_dft_r2c(plan_distfftw_inverse_wf,bfft_mpifftw_kimg1,afft_mpifftw_kimg1)
    endif
    call tstatc0_end(id_sname)
  end subroutine m_FFT_Inverse_MPI_FFTW

  subroutine m_FFT_Direct_MPI_FFTW(nfout)
    integer, intent(in) :: nfout
    integer :: id_sname=-1
    call tstatc0_begin('m_FFT_Direct_MPI_FFTW ',id_sname)
    if(kimg==2) then
      call fftw_mpi_execute_dft(plan_distfftw_direct_wf,afft_mpifftw, bfft_mpifftw)
    else
      call fftw_mpi_execute_dft_c2r(plan_distfftw_direct_wf,afft_mpifftw_kimg1, bfft_mpifftw_kimg1)
    endif
    call tstatc0_end(id_sname)
  end subroutine m_FFT_Direct_MPI_FFTW
#endif

!------------------------------------------------------------------------------
  subroutine m_FFT_CD_Direct_XYZ_3D(nfout,wk_afft_l, wk_size, iesize)
    integer, intent(in) :: nfout, wk_size, iesize
    real(kind=DP), dimension(wk_size*kimg,iesize), intent(inout) :: wk_afft_l

    integer, parameter :: FFTW_MEASURE=0
    integer, parameter :: FFTW_ESTIMATE=64
    integer, parameter :: FFTW_RANK=1
    integer            :: FFTW_FLAG
    integer,dimension(FFTW_RANK) :: NFFTW3, NEMBED, NEREAL
    integer :: nx, ny, nz, i, j, k, l, ri, ib, iadd, lrank
    integer :: irank, itag, icnt_send, icnt_recv
    integer,save, allocatable, dimension(:)   :: req_r, req_s
    integer,save, allocatable, dimension(:,:) :: sta_r, sta_s
    integer,save, allocatable, dimension(:)   :: wk_recvcnt, wk_sendcnt
    real(kind=DP),allocatable, dimension(:,:),save :: wk_recv1, wk_send1
    real(kind=DP),allocatable, dimension(:,:),save :: wk_recv2, wk_send2
    integer :: mpi_comm, myrank, nmrank
    integer(kind=8),save :: planx1 = 0, planx2 = 0
    integer(kind=8),save :: plany1 = 0, plany2 = 0
    integer(kind=8),save :: planz1 = 0, planz2 = 0
    integer(kind=8) :: plan
    integer :: max_x, max_y, max_z
    integer,save :: max_elm = 0, savesize = 0

    integer,save, allocatable, dimension(:,:) :: z2y_recv, z2y_send, y2x_recv, y2x_send
    integer,save :: z2y_rrank, z2y_srank, y2x_rrank, y2x_srank
    integer,save :: z2y_rmax, z2y_smax, z2y_srmax, y2x_rmax, y2x_smax, y2x_srmax

    integer ::  nnx, nny, nnz, iadd0, iadd1, nis, nie, jadd, kadd, nisx, niex, nxp, nyp, nzp

    integer,dimension(3) :: nsize, isin
#ifdef FFT_ALLTOALL
    integer,dimension(2) :: max_send,max_recv
#endif

#ifdef __TIMER_SUB__
    call mpi_barrier(mpi_ke_world,ierr)
    call timer_sta(109)
#endif

    if (kimg == 1) then
       FFTW_FLAG = +1      ! FFTW_BACKWARD
    else
       FFTW_FLAG = -1      ! FFTW_FORWARD
    endif
    plan = 0

    mpi_comm = mpi_ke_world
    myrank = myrank_g
    nmrank = nrank_g
    itag = 10

    if (np_fftcd_x==0)then
       nx = 0
       nx = 0
       nz = 0
       go to 1000
    endif

    if (firstcall_cd_direct_xyz_3d) then
       savesize = 0
#ifdef __TIMER_DO__
  call timer_sta(157)
#endif
       max_x = maxval(nel_fftcd_x(:))
       max_y = maxval(nel_fftcd_y(:))
       max_z = maxval(nel_fftcd_z(:))
       max_elm = max(max_x,max_y,max_z)

       if(allocated(req_r)) deallocate(req_r)
       if(allocated(req_s)) deallocate(req_s)
       if(allocated(sta_r)) deallocate(sta_r)
       if(allocated(sta_s)) deallocate(sta_s)
       if(allocated(wk_recvcnt)) deallocate(wk_recvcnt)
       if(allocated(wk_sendcnt)) deallocate(wk_sendcnt)
       allocate(req_r(0:nmrank-1), stat=ierr)
       allocate(req_s(0:nmrank-1), stat=ierr)
       allocate(sta_r(MPI_STATUS_SIZE,0:nmrank-1), stat=ierr)
       allocate(sta_s(MPI_STATUS_SIZE,0:nmrank-1), stat=ierr)
       allocate(wk_recvcnt(0:nmrank-1), stat=ierr)
       allocate(wk_sendcnt(0:nmrank-1), stat=ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_FFT_CD_Direct_3D :  Not allocate '
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 51, ierr)
        endif

       if(allocated(z2y_recv)) deallocate(z2y_recv)
       if(allocated(z2y_send)) deallocate(z2y_send)
       allocate(z2y_recv(2,0:nmrank-1))
       allocate(z2y_send(2,0:nmrank-1))
       z2y_send = 0
       z2y_recv = 0
       wk_recvcnt = 0
       wk_sendcnt = 0
       do i = 1, nel_fftcd_y(myrank)
          irank = map_fftcd_z(mp_fftcd_y(i)) - 1
          wk_recvcnt(irank) = wk_recvcnt(irank) + 1
       enddo
       do i = 1, nel_fftcd_z(myrank)
          irank = map_fftcd_y(mp_fftcd_z(i)) - 1
          wk_sendcnt(irank) = wk_sendcnt(irank) + 1
       enddo
       k = 0
       do i = 0, nmrank - 1
          if(wk_recvcnt(i) /= 0) then
             k = k + 1
             z2y_recv(1,i) = wk_recvcnt(i)
             z2y_recv(2,i) = k
          endif
       enddo
       z2y_rrank = k
       k = 0
       do i = 0, nmrank - 1
          if(wk_sendcnt(i) /= 0) then
             k = k + 1
             z2y_send(1,i) = wk_sendcnt(i)
             z2y_send(2,i) = k
          endif
       enddo
       z2y_srank = k
       z2y_rmax = maxval(wk_recvcnt)
       z2y_smax = maxval(wk_sendcnt)
       z2y_srmax = max(z2y_rmax,z2y_smax)

       if(allocated(y2x_recv)) deallocate(y2x_recv)
       if(allocated(y2x_send)) deallocate(y2x_send)
       allocate(y2x_recv(2,0:nmrank-1))
       allocate(y2x_send(2,0:nmrank-1))
       y2x_send = 0
       y2x_recv = 0
       wk_recvcnt = 0
       wk_sendcnt = 0
       do i = 1, nel_fftcd_x(myrank)
          irank = map_fftcd_y(mp_fftcd_x(i)) - 1
          wk_recvcnt(irank) = wk_recvcnt(irank) + 1
       enddo
       do i = 1, nel_fftcd_y(myrank)
          irank = map_fftcd_x(mp_fftcd_y(i)) - 1
          wk_sendcnt(irank) = wk_sendcnt(irank) + 1
       enddo
       k = 0
       do i = 0, nmrank - 1
          if(wk_recvcnt(i) /= 0) then
             k = k + 1
             y2x_recv(1,i) = wk_recvcnt(i)
             y2x_recv(2,i) = k
          endif
       enddo
       y2x_rrank = k
       k = 0
      do i = 0, nmrank - 1
          if(wk_sendcnt(i) /= 0) then
             k = k + 1
             y2x_send(1,i) = wk_sendcnt(i)
             y2x_send(2,i) = k
          endif
       enddo
       y2x_srank = k
       y2x_rmax = maxval(wk_recvcnt)
       y2x_smax = maxval(wk_sendcnt)
       y2x_srmax = max(y2x_rmax,y2x_smax)

#ifdef FFT_ALLTOALL
       max_send(1) = z2y_rmax
       max_send(2) = z2y_smax
       call mpi_allreduce(max_send,max_recv,2,mpi_integer,mpi_max,mpi_fftcd_yz_world,ierr)
       z2y_rmax = max_recv(1)
       z2y_smax = max_recv(2)
       max_send(1) = y2x_rmax
       max_send(2) = y2x_smax
       call mpi_allreduce(max_send,max_recv,2,mpi_integer,mpi_max,mpi_fftcd_xy_world,ierr)
       y2x_rmax = max_recv(1)
       y2x_smax = max_recv(2)
#endif
       if (ipri > 1) then
          write(nfout,'("m_FFT_CD_Direct_XYZ_3D   --   myrank_g=",i4)') myrank_g
          write(nfout,'("z2y_send")')
          write(nfout,'(10(i8,", "))') (z2y_send(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (z2y_send(2,i),i=0,nmrank-1)
          write(nfout,'("z2y_recv")')
          write(nfout,'(10(i8,", "))') (z2y_recv(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (z2y_recv(2,i),i=0,nmrank-1)
          write(nfout,'("y2x_send")')
          write(nfout,'(10(i8,", "))') (y2x_send(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (y2x_send(2,i),i=0,nmrank-1)
          write(nfout,'("y2x_recv")')
          write(nfout,'(10(i8,", "))') (y2x_recv(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (y2x_recv(2,i),i=0,nmrank-1)
       endif

       nx = xyz_fftcd_z(2,1)-xyz_fftcd_z(1,1)+1
       ny = xyz_fftcd_z(2,2)-xyz_fftcd_z(1,2)+1
       nz = xyz_fftcd_z(2,3)-xyz_fftcd_z(1,3)+1

       firstcall_cd_direct_xyz_3d = .false.
#ifdef __TIMER_DO__
  call timer_end(157)
#endif
    endif

#ifdef FFT_ALLTOALL
    if (iesize /= savesize) then
#else
    if (iesize > savesize) then
#endif
#ifdef __TIMER_DO__
  call timer_sta(158)
#endif
       if (allocated(wk_recv1)) deallocate(wk_recv1)
       if (allocated(wk_send1)) deallocate(wk_send1)
       if (allocated(wk_recv2)) deallocate(wk_recv2)
       if (allocated(wk_send2)) deallocate(wk_send2)
       allocate(wk_recv1(z2y_rmax*kimg*iesize,z2y_rrank), stat=ierr)
       allocate(wk_send1(z2y_smax*kimg*iesize,z2y_srank), stat=ierr)
       allocate(wk_recv2(y2x_rmax*kimg*iesize,y2x_rrank), stat=ierr)
       allocate(wk_send2(y2x_smax*kimg*iesize,y2x_srank), stat=ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_FFT_Direct_3D :  Not allocate '
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 52, ierr)
        endif
       savesize = iesize
#ifdef __TIMER_DO__
  call timer_end(158)
#endif
    end if

    nx = xyz_fftcd_z(2,1)-xyz_fftcd_z(1,1)+1
    ny = xyz_fftcd_z(2,2)-xyz_fftcd_z(1,2)+1
    nz = xyz_fftcd_z(2,3)-xyz_fftcd_z(1,3)+1
    nzp = nz
!
! Z-axis (x-y div)
!
#ifdef __TIMER_DO__
  call timer_sta(257)
#endif

#ifdef FFT_USE_SSL2
    if(kimg==1) then
       nsize(1:3) = (/nz,(nx/2),ny/)
       isin(1:3)  = (/-1,0,0/)
    else
       nsize(1:3) = (/nz,nx,ny/)
       isin(1:3)  = (/1,0,0/)
    end if
!OCL INDEPENDENT[ DVMCF2 ]
!OCL PARALLEL_STRONG
    do ib = 1, iesize
       call DVMCF2(wk_afft_l(1,ib),nsize,3,isin,ierr)
    end do
#else
    if(kimg==1) then
       NFFTW3(1) = fft_box_size_CD_3D(3,1)
       NEMBED(1) = fft_box_size_CD_3D(3,1)
       call dfftw_plan_many_dft    (plan,           1, NFFTW3, nx*ny/2, &
      &                             wk_afft_l, NEMBED,      1,      nz, &
      &                             wk_afft_l, NEMBED,      1,      nz, &
      &                             FFTW_FLAG, FFTW_ESTIMATE)
    else
       NFFTW3(1) = fft_box_size_CD_3D(3,1)
       NEMBED(1) = fft_box_size_CD_3D(3,0)
       call dfftw_plan_many_dft    (plan,           1, NFFTW3, nx*ny,   &
      &                             wk_afft_l, NEMBED,      1,    nz,   &
      &                             wk_afft_l, NEMBED,      1,    nz,   &
      &                             FFTW_FLAG, FFTW_ESTIMATE)
    end if
!OCL INDEPENDENT[dfftw_execute_dft]
!OCL PARALLEL_STRONG
    do ib = 1, iesize
       call dfftw_execute_dft(plan,wk_afft_l(1,ib),wk_afft_l(1,ib))
    end do
    call dfftw_destroy_plan(plan)
#endif

#ifdef __TIMER_DO__
  call timer_end(257)
#endif

#ifdef __TIMER_DO__
  call timer_sta(160)
#endif
#ifdef __TIMER_DO__
  call timer_end(160)
#endif

#ifdef __TIMER_DO__
  call timer_sta(161)
#endif
! (z,x,y) -> (z,x,y)
    if (kimg == 1) then
!OCL NORECURRENCE
!OCL PARALLEL_STRONG
       do ib = 1, iesize
!OCL SERIAL
          do j = 1, ny
!OCL SERIAL
             do i = 1, nx/2
!OCL SERIAL
                do l = 1, fftcd_Y_z_dim
                   nis = nis_fftcd_Y_z(l)
                   nie = nie_fftcd_Y_z(l)
                   nnz = nie - nis + 1
                   iadd0 = nx/2*ny*nnz*(ib-1)*2
!OCL SERIAL
                   do k = nis, nie
                      iadd  = (k+(i-1)*nzp+(j-1)*nzp*nx/2)*2
                      iadd1 = iadd0+((k-nis+1)+(i-1)*nnz+(j-1)*nnz*nx/2)*2
                      wk_send1(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send1(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
    else
!OCL NORECURRENCE
!OCL PARALLEL_STRONG
       do ib = 1, iesize
!OCL SERIAL
          do j = 1, ny
!OCL SERIAL
             do i = 1, nx
!OCL SERIAL
                do l = 1, fftcd_Y_z_dim
                   nis = nis_fftcd_Y_z(l)
                   nie = nie_fftcd_Y_z(l)
                   nnz = nie - nis + 1
                   iadd0 = nx*ny*nnz*(ib-1)*2
!OCL SERIAL
                   do k = nis, nie
                      iadd  = (k+(i-1)*nzp+(j-1)*nzp*nx)*2
                      iadd1 = iadd0+((k-nis+1)+(i-1)*nnz+(j-1)*nnz*nx)*2
                      wk_send1(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send1(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
    end if
#ifdef __TIMER_DO__
  call timer_end(161)
#endif

#ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_fftcd_yz_world)
  call timer_sta(162)
#endif

    call MPI_ALLTOALL(wk_send1, z2y_smax*kimg*iesize, mpi_double_precision,   &
   &                  wk_recv1, z2y_rmax*kimg*iesize, mpi_double_precision,   &
   &                                                  mpi_fftcd_yz_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_CD_Direct_XYZ_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 55, ierr)
     endif

#ifdef __TIMER_DO__
  call timer_sta(163)
#endif
#ifdef __TIMER_DO__
  call timer_end(163)
#endif

#else
!else ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_comm)
  call timer_sta(162)
#endif

    icnt_recv = 0
    lrank = myrank + 1
    if (lrank > (nmrank-1)) lrank = 0
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank-1)) lrank = 0
       if ((y2z_recv(1,lrank) /= 0)) then
          call mpi_irecv(wk_recv1(1,z2y_recv(2,lrank)), z2y_recv(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_r(icnt_recv), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_CD_Direct_XYZ_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 55, ierr)
           endif
          icnt_recv = icnt_recv + 1
       endif
    enddo

    icnt_send = 0
    lrank = myrank
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank-1)) lrank = 0
       if ((y2z_send(1,lrank) /= 0)) then
          call mpi_isend(wk_send1(1,z2y_send(2,lrank)), z2y_send(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_s(icnt_send), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_CD_Direct_XYZ_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 56, ierr)
           endif
          icnt_send = icnt_send + 1
       endif
    enddo

#ifdef __TIMER_DO__
  call timer_sta(163)
#endif
#ifdef __TIMER_DO__
  call timer_end(163)
#endif

    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_CD_Direct_XYZ_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 57, ierr)
     endif

    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_CD_Direct_XYZ_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 58, ierr)
     endif

#endif
!endif ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_end(162)
#endif

#ifdef __TIMER_DO__
  call timer_sta(164)
#endif
! (z,x,y) -> (y,z,x)
    nx = xyz_fftcd_y(2,1)-xyz_fftcd_y(1,1)+1
    ny = xyz_fftcd_y(2,2)-xyz_fftcd_y(1,2)+1
    nz = xyz_fftcd_y(2,3)-xyz_fftcd_y(1,3)+1
    nyp = ny
    if (kimg == 1) then
       do ib = 1, iesize
          do i = 1, nx/2
             do k = 1, nz
                do l = 1, fftcd_Z_y_dim
                   nis = nis_fftcd_Z_y(l)
                   nie = nie_fftcd_Z_y(l)
                   nny = nie - nis + 1
                   iadd0 = nx/2*nny*nz*(ib-1)*2
                   do j = nis, nie
                      iadd1 = iadd0 + (k+(i-1)*nz+(j-nis)*nz*nx/2)*2
                      iadd = (j+(k-1)*nyp+(i-1)*nyp*nz)*2
                      wk_afft_l(iadd-1,ib) = wk_recv1(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv1(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end do
    else
       do ib = 1, iesize
          do i = 1, nx
             do k = 1, nz
                do l = 1, fftcd_Z_y_dim
                   nis = nis_fftcd_Z_y(l)
                   nie = nie_fftcd_Z_y(l)
                   nny = nie - nis + 1
                   iadd0 = nx*nny*nz*(ib-1)*2
                   do j = nis, nie
                      iadd1 = iadd0+(k+(i-1)*nz+(j-nis)*nz*nx)*2
                      iadd  = (j+(k-1)*nyp+(i-1)*nyp*nz)*2
                      wk_afft_l(iadd-1,ib) = wk_recv1(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv1(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end do
    end if
#ifdef __TIMER_DO__
  call timer_end(164)
#endif
!
! Y-axis (z-x div)
!
#ifdef __TIMER_DO__
  call timer_sta(258)
#endif

#ifdef FFT_USE_SSL2
    if(kimg==1) then
       nsize(1:3) = (/ny,nz,(nx/2)/)
       isin(1:3)  = (/-1,0,0/)
    else
       nsize(1:3) = (/ny,nz,nx/)
       isin(1:3)  = (/1,0,0/)
    end if
!OCL INDEPENDENT[ DVMCF2 ]
!OCL PARALLEL_STRONG
    do ib = 1, iesize
       call DVMCF2(wk_afft_l(1,ib),nsize,3,isin,ierr)
    end do
#else
    if(kimg==1) then
       NFFTW3(1) = fft_box_size_CD_3D(2,1)
       NEMBED(1) = fft_box_size_CD_3D(2,1)
       call dfftw_plan_many_dft    (plan,           1, NFFTW3, nz*nx/2, &
      &                             wk_afft_l, NEMBED,      1,      ny, &
      &                             wk_afft_l, NEMBED,      1,      ny, &
      &                             FFTW_FLAG, FFTW_ESTIMATE)
    else
       NFFTW3(1) = fft_box_size_CD_3D(2,1)
       NEMBED(1) = fft_box_size_CD_3D(2,0)
       call dfftw_plan_many_dft    (plan,           1, NFFTW3, nz*nx,  &
      &                             wk_afft_l, NEMBED,      1,    ny,  &
      &                             wk_afft_l, NEMBED,      1,    ny,  &
      &                             FFTW_FLAG, FFTW_ESTIMATE)
    endif
!OCL INDEPENDENT[dfftw_execute_dft]
!OCL PARALLEL_STRONG
    do ib = 1, iesize
       call dfftw_execute_dft(plan,wk_afft_l(1,ib),wk_afft_l(1,ib))
    end do
    call dfftw_destroy_plan(plan)
#endif

#ifdef __TIMER_DO__
  call timer_end(258)
#endif

#ifdef __TIMER_DO__
  call timer_sta(165)
#endif
! (y,z,x) -> (y,z,x)
    if (kimg == 1) then
!OCL NORECURRENCE
!OCL PARALLEL_STRONG
       do ib = 1, iesize
!OCL SERIAL
          do i = 1, nx/2
!OCL SERIAL
             do k = 1, nz
!OCL SERIAL
                do l = 1, fftcd_X_y_dim
                   nis = nis_fftcd_X_y(l)
                   nie = nie_fftcd_X_y(l)
                   nny = nie - nis + 1
                   iadd0 = nx/2*nny*nz*(ib-1)*2
!OCL SERIAL
                   do j = nis, nie
                      iadd  = (j+(k-1)*nyp+(i-1)*nyp*nz)*2
                      iadd1 = iadd0 + (j-nis+1+(k-1)*nny+(i-1)*nny*nz)*2
                      wk_send2(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send2(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
    else
!OCL NORECURRENCE
!OCL PARALLEL_STRONG
       do ib = 1, iesize
!OCL SERIAL
          do i = 1, nx
!OCL SERIAL
             do k = 1, nz
!OCL SERIAL
                do l = 1, fftcd_X_y_dim
                   nis = nis_fftcd_X_y(l)
                   nie = nie_fftcd_X_y(l)
                   nny = nie - nis + 1
                   iadd0 = nx*nny*nz*(ib-1)*2
!OCL SERIAL
                   do j = nis, nie
                      iadd  = (j+(k-1)*nyp+(i-1)*nyp*nz)*2
                      iadd1 = iadd0+(j-nis+1+(k-1)*nny+(i-1)*nny*nz)*2
                      wk_send2(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send2(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
    end if
#ifdef __TIMER_DO__
  call timer_end(165)
#endif

#ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_fftcd_xy_world)
  call timer_sta(166)
#endif

    call MPI_ALLTOALL(wk_send2, y2x_smax*kimg*iesize, mpi_double_precision,   &
   &                  wk_recv2, y2x_rmax*kimg*iesize, mpi_double_precision,   &
   &                                                  mpi_fftcd_xy_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_CD_Direct_XYZ_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 60, ierr)
     endif

#ifdef __TIMER_DO__
  call timer_sta(167)
#endif
#ifdef __TIMER_DO__
  call timer_end(167)
#endif

#else
!else ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_comm)
  call timer_sta(166)
#endif

    icnt_recv = 0
    lrank = myrank + 1
    if (lrank > (nmrank-1)) lrank = 0
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank-1)) lrank = 0
       if ((z2x_recv(1,lrank) /= 0)) then
          call mpi_irecv(wk_recv2(1,z2x_recv(2,lrank)), z2x_recv(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_r(icnt_recv), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_CD_Direct_XYZ_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 60, ierr)
           endif
          icnt_recv = icnt_recv + 1
       endif
    enddo

    icnt_send = 0
    lrank = myrank
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank-1)) lrank = 0
       if ((z2x_send(1,lrank) /= 0)) then
          call mpi_isend(wk_send2(1,z2x_send(2,lrank)), z2x_send(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_s(icnt_send), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_CD_Direct_XYZ_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 61, ierr)
           endif
          icnt_send = icnt_send + 1
       endif
    enddo

#ifdef __TIMER_DO__
  call timer_sta(167)
#endif
#ifdef __TIMER_DO__
  call timer_end(167)
#endif

    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_CD_Direct_XYZ_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 62, ierr)
     endif

    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_CD_Direct_XYZ_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 63, ierr)
     endif
#endif
!endif ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_end(166)
#endif

#ifdef __TIMER_DO__
  call timer_sta(168)
#endif
    wk_afft_l = 0.0d0

    nx = xyz_fftcd_x(2,1)-xyz_fftcd_x(1,1)+1
    ny = xyz_fftcd_x(2,2)-xyz_fftcd_x(1,2)+1
    nz = xyz_fftcd_x(2,3)-xyz_fftcd_x(1,3)+1
! (y,z,x) -> (x,y,z)
    if (kimg == 1) then
       do ib = 1, iesize
          do k = 1, nz
             do j = 1, ny
                do l = 1, fftcd_Y_x_dim
                   nis = nis_fftcd_Y_x(l)
                   nie = nie_fftcd_Y_x(l)
                   nnx = nie - nis + 1
                   if(mod(nis,2)>0) then
                     nisx = nis/2 + 1
                   else
                     nisx = nis/2
                   endif
                   niex = nie/2
                   iadd0 = nnx/2*ny*nz*(ib-1)*2
                   do i = nisx, niex
                      iadd1 = iadd0 + (j+(k-1)*ny+(i-nisx)*ny*nz)*2
                      iadd  = (i+(j-1)*nx/2+(k-1)*nx/2*ny)*2
                      wk_afft_l(iadd-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end do
    else
       do ib = 1, iesize
          do k = 1, nz
             do j = 1, ny
                do l = 1, fftcd_Y_x_dim
                   nis = nis_fftcd_Y_x(l)
                   nie = nie_fftcd_Y_x(l)
                   nnx = nie - nis + 1
                   iadd0 = nnx*ny*nz*(ib-1)*2
                   do i = nis, nie
                      iadd1 = iadd0+(j+(k-1)*ny+(i-nis)*ny*nz)*2
                      iadd  = (i+(j-1)*nx+(k-1)*nx*ny)*2
                      wk_afft_l(iadd-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end do
    end if
#ifdef __TIMER_DO__
  call timer_end(168)
#endif
!
! X-axis (y-z div)
!
#ifdef __TIMER_DO__
  call timer_sta(259)
#endif

#ifdef FFT_USE_SSL2
    if(kimg==1) then
       nsize(1:3) = (/nx-2,ny,nz/)
       isin(1:3)  = (/-1,0,0/)
!OCL INDEPENDENT[ DVMRF2 ]
!OCL PARALLEL_STRONG
       do ib = 1, iesize
          call DVMRF2(wk_afft_l(1,ib),nsize,3,isin,-1,ierr)
       end do
    else
       nsize(1:3) = (/nx,ny,nz/)
       isin(1:3)  = (/1,0,0/)
!OCL INDEPENDENT[ DVMCF2 ]
!OCL PARALLEL_STRONG
       do ib = 1, iesize
          call DVMCF2(wk_afft_l(1,ib),nsize,3,isin,ierr)
       end do
    end if
#else
    if(kimg==1) then
       NFFTW3(1) = nx - 2
       NEMBED(1) = nx
       NEREAL(1) = nx
       call dfftw_plan_many_dft_c2r(plan,           1, NFFTW3, ny*nz,  &
      &                             wk_afft_l, NEMBED,      1,  nx/2,  &
      &                             wk_afft_l, NEREAL,      1,    nx,  &
      &                             FFTW_ESTIMATE)
    else
       NFFTW3(1) = fft_box_size_CD_3D(1,1)
       NEMBED(1) = fft_box_size_CD_3D(1,1)
       call dfftw_plan_many_dft    (plan,           1, NFFTW3, ny*nz,  &
      &                             wk_afft_l, NEMBED,      1,    nx,  &
      &                             wk_afft_l, NEMBED,      1,    nx,  &
      &                             FFTW_FLAG, FFTW_ESTIMATE)
    endif
!OCL INDEPENDENT[dfftw_execute_dft_r2c,dfftw_execute_dft]
!OCL PARALLEL_STRONG
    do ib = 1, iesize
       if (kimg==1) then
          call dfftw_execute_dft_c2r(plan,wk_afft_l(1,ib),wk_afft_l(1,ib))
       else
          call dfftw_execute_dft    (plan,wk_afft_l(1,ib),wk_afft_l(1,ib))
       endif
    end do
    call dfftw_destroy_plan(plan)
#endif
#ifdef __TIMER_DO__
  call timer_end(259)
#endif

1000 continue

#ifdef __TIMER_SUB__
    call timer_end(109)
#endif
  end subroutine m_FFT_CD_Direct_XYZ_3D
!------------------------------------------------------------------------------

  subroutine m_FFT_CD_Inverse_XYZ_3D(nfout,wk_afft_l, wk_size, iesize)
    integer, intent(in) :: nfout, wk_size, iesize
    real(kind=DP), dimension(wk_size*kimg,iesize), intent(inout) :: wk_afft_l

    integer, parameter :: FFTW_MEASURE=0
    integer, parameter :: FFTW_ESTIMATE=64
    integer, parameter :: FFTW_RANK=1
    integer            :: FFTW_FLAG
    integer,dimension(1) :: NFFTW3, NEMBED, NEREAL
    integer :: nx, ny, nz, i, j, k, l, ri, ib, iadd, lrank
    integer :: irank, itag, icnt_send, icnt_recv
    integer,save, allocatable, dimension(:)   :: req_r, req_s
    integer,save, allocatable, dimension(:,:) :: sta_r, sta_s
    integer,save, allocatable, dimension(:)   :: wk_recvcnt, wk_sendcnt
    real(kind=DP),allocatable, dimension(:,:),save :: wk_recv1, wk_send1
    real(kind=DP),allocatable, dimension(:,:),save :: wk_recv2, wk_send2
    integer,      allocatable, dimension(:)   :: wk_mp_fft_y
    integer :: mpi_comm, myrank, nmrank
    integer(kind=8),save :: planx1 = 0, planx2 = 0
    integer(kind=8),save :: plany1 = 0, plany2 = 0
    integer(kind=8),save :: planz1 = 0, planz2 = 0
    integer(kind=8) :: plan
    integer :: max_x, max_y, max_z
    integer,save :: max_elm = 0, savesize = 0

    integer,save, allocatable, dimension(:,:) :: x2y_recv, x2y_send, y2z_recv, y2z_send
    integer,save :: x2y_rrank, x2y_srank, y2z_rrank, y2z_srank
    integer,save :: x2y_rmax, x2y_smax, x2y_srmax, y2z_rmax, y2z_smax, y2z_srmax

    integer ::  nnx, nny, nnz, iadd0, iadd1, nis, nie, jadd, kadd, nxp, nyp, nzp
    integer                  , dimension(3) :: isin,nsize
#ifdef FFT_ALLTOALL
    integer,dimension(2) :: max_send,max_recv
#endif

#ifdef __TIMER_SUB__
    call timer_barrier(mpi_ke_world)
    call timer_sta(110)
#endif

    if (kimg == 1) then
       FFTW_FLAG = -1      ! FFTW_FORWARD
    else
       FFTW_FLAG = +1      ! FFTW_BACKWARD
    endif

    mpi_comm = mpi_ke_world
    myrank = myrank_g
    nmrank = nrank_g
    itag = 10

    if (np_fftcd_x==0)then
       nx = 0
       nx = 0
       nz = 0
       go to 1000
    endif

    if (firstcall_cd_inverse_xyz_3d) then
       savesize = 0
#ifdef __TIMER_DO__
  call timer_sta(173)
#endif
       max_x = maxval(nel_fftcd_x(:))
       max_y = maxval(nel_fftcd_y(:))
       max_z = maxval(nel_fftcd_z(:))
       max_elm = max(max_x,max_y,max_z)
       if(allocated(req_r)) deallocate(req_r)
       if(allocated(req_s)) deallocate(req_s)
       if(allocated(sta_r)) deallocate(sta_r)
       if(allocated(sta_s)) deallocate(sta_s)
       if(allocated(wk_recvcnt)) deallocate(wk_recvcnt)
       if(allocated(wk_sendcnt)) deallocate(wk_sendcnt)
       allocate(req_r(0:nmrank-1), stat=ierr)
       allocate(req_s(0:nmrank-1), stat=ierr)
       allocate(sta_r(MPI_STATUS_SIZE,0:nmrank-1), stat=ierr)
       allocate(sta_s(MPI_STATUS_SIZE,0:nmrank-1), stat=ierr)
       allocate(wk_recvcnt(0:nmrank-1), stat=ierr)
       allocate(wk_sendcnt(0:nmrank-1), stat=ierr)
       if (ierr /= 0) then
          write(nfout,*)' m_FFT_Inverse_3D :  Not allocate '
          call flush(nfout)
          call mpi_abort(mpi_comm_world, 68, ierr)
       endif

       if(allocated(x2y_recv)) deallocate(x2y_recv)
       if(allocated(x2y_send)) deallocate(x2y_send)
       allocate(x2y_recv(2,0:nmrank-1))
       allocate(x2y_send(2,0:nmrank-1))
       x2y_send = 0
       x2y_recv = 0
       wk_recvcnt = 0
       wk_sendcnt = 0
       do i = 1, nel_fftcd_y(myrank)
          irank = map_fftcd_x(mp_fftcd_y(i)) - 1
          wk_recvcnt(irank) = wk_recvcnt(irank) + 1
       enddo
       do i = 1, nel_fftcd_x(myrank)
          irank = map_fftcd_y(mp_fftcd_x(i)) - 1
          wk_sendcnt(irank) = wk_sendcnt(irank) + 1
       enddo
       k = 0
       do i = 0, nmrank - 1
          if(wk_recvcnt(i) /= 0) then
             k = k + 1
             x2y_recv(1,i) = wk_recvcnt(i)
             x2y_recv(2,i) = k
          endif
       enddo
       x2y_rrank = k
       k = 0
       do i = 0, nmrank - 1
          if(wk_sendcnt(i) /= 0) then
             k = k + 1
             x2y_send(1,i) = wk_sendcnt(i)
             x2y_send(2,i) = k
          endif
       enddo
       x2y_srank = k
       x2y_rmax = maxval(wk_recvcnt)
       x2y_smax = maxval(wk_sendcnt)
       x2y_srmax = max(x2y_rmax,x2y_smax)

       if(allocated(y2z_recv)) deallocate(y2z_recv)
       if(allocated(y2z_send)) deallocate(y2z_send)
       allocate(y2z_recv(2,0:nmrank-1))
       allocate(y2z_send(2,0:nmrank-1))
       y2z_send = 0
       y2z_recv = 0
       wk_recvcnt = 0
       wk_sendcnt = 0

       nx = xyz_fftcd_y(2,1)-xyz_fftcd_y(1,1)+1
       ny = xyz_fftcd_y(2,2)-xyz_fftcd_y(1,2)+1
       nz = xyz_fftcd_y(2,3)-xyz_fftcd_y(1,3)+1
       allocate(wk_mp_fft_y(nx*ny*nz) ,stat=ierr)
       do k = 0, nz-1
          do j = 0, ny-1
             do i = 0, nx-1
                wk_mp_fft_y(i+j*nx+k*nx*ny+1) = mp_fftcd_y(i+k*nx+j*nx*nz+1)
             enddo
          enddo
       enddo

       do i = 1, nel_fftcd_z(myrank)
          irank = map_fftcd_y(mp_fftcd_z(i)) - 1
          wk_recvcnt(irank) = wk_recvcnt(irank) + 1
       enddo
       do i = 1, nel_fftcd_y(myrank)
          irank = map_fftcd_z(wk_mp_fft_y(i)) - 1
          wk_sendcnt(irank) = wk_sendcnt(irank) + 1
       enddo
       k = 0
       do i = 0, nmrank - 1
          if(wk_recvcnt(i) /= 0) then
             k = k + 1
             y2z_recv(1,i) = wk_recvcnt(i)
             y2z_recv(2,i) = k
          endif
       enddo
       y2z_rrank = k
       k = 0
      do i = 0, nmrank - 1
          if(wk_sendcnt(i) /= 0) then
             k = k + 1
             y2z_send(1,i) = wk_sendcnt(i)
             y2z_send(2,i) = k
          endif
       enddo
       y2z_srank = k
       y2z_rmax = maxval(wk_recvcnt)
       y2z_smax = maxval(wk_sendcnt)
       y2z_srmax = max(y2z_rmax,y2z_smax)

       deallocate(wk_mp_fft_y)

#ifdef FFT_ALLTOALL
       max_send(1) = x2y_rmax
       max_send(2) = x2y_smax
       call mpi_allreduce(max_send,max_recv,2,mpi_integer,mpi_max,mpi_fftcd_xy_world,ierr)
       x2y_rmax = max_recv(1)
       x2y_smax = max_recv(2)
       max_send(1) = y2z_rmax
       max_send(2) = y2z_smax
       call mpi_allreduce(max_send,max_recv,2,mpi_integer,mpi_max,mpi_fftcd_yz_world,ierr)
       y2z_rmax = max_recv(1)
       y2z_smax = max_recv(2)
#endif
       if (ipri > 1) then
          write(nfout,'("m_FFT_CD_Inverse_XYZ_3D   --   myrank_g=",i4)') myrank_g
          write(nfout,'("x2y_send")')
          write(nfout,'(10(i8,", "))') (x2y_send(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (x2y_send(2,i),i=0,nmrank-1)
          write(nfout,'("x2y_recv")')
          write(nfout,'(10(i8,", "))') (x2y_recv(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (x2y_recv(2,i),i=0,nmrank-1)
          write(nfout,'("y2z_send")')
          write(nfout,'(10(i8,", "))') (y2z_send(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (y2z_send(2,i),i=0,nmrank-1)
          write(nfout,'("y2z_recv")')
          write(nfout,'(10(i8,", "))') (y2z_recv(1,i),i=0,nmrank-1)
          write(nfout,'(10(i8,", "))') (y2z_recv(2,i),i=0,nmrank-1)
          call flush(nfout)
       endif

       firstcall_cd_inverse_xyz_3d = .false.
#ifdef __TIMER_DO__
  call timer_end(173)
#endif
    endif

#ifdef FFT_ALLTOALL
    if (iesize /= savesize) then
#else
    if (iesize > savesize) then
#endif
#ifdef __TIMER_DO__
  call timer_sta(174)
#endif
       if (allocated(wk_recv1)) deallocate(wk_recv1)
       if (allocated(wk_send1)) deallocate(wk_send1)
       if (allocated(wk_recv2)) deallocate(wk_recv2)
       if (allocated(wk_send2)) deallocate(wk_send2)
       allocate(wk_recv1(x2y_rmax*kimg*iesize,x2y_rrank), stat=ierr)
       allocate(wk_send1(x2y_smax*kimg*iesize,x2y_srank), stat=ierr)
       allocate(wk_recv2(y2z_rmax*kimg*iesize,y2z_rrank), stat=ierr)
       allocate(wk_send2(y2z_smax*kimg*iesize,y2z_srank), stat=ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_FFT_CD_Inverse_XYZ_3D :  Not allocate '
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 69, ierr)
        endif
       savesize = iesize
#ifdef __TIMER_DO__
  call timer_end(174)
#endif
    end if
!
! X-axis (y-z div)
!
    nx = xyz_fftcd_x(2,1)-xyz_fftcd_x(1,1)+1
    ny = xyz_fftcd_x(2,2)-xyz_fftcd_x(1,2)+1
    nz = xyz_fftcd_x(2,3)-xyz_fftcd_x(1,3)+1
    nxp = nx
#ifdef __TIMER_DO__
  call timer_sta(260)
#endif

#ifdef FFT_USE_SSL2
    if (kimg==1) then
       nsize(1:3) = (/nx-2,ny,nz/)
       isin(1:3)  = (/1,0,0/)
!OCL INDEPENDENT[ DVMRF2 ]
!OCL PARALLEL_STRONG
       do ib = 1, iesize
          call DVMRF2(wk_afft_l(1,ib),nsize,3,isin,1,ierr)
       end do
    else
       nsize(1:3) = (/nx,ny,nz/)
       isin(1:3)  = (/-1,0,0/)
!OCL INDEPENDENT[ DVMCF2 ]
!OCL PARALLEL_STRONG
       do ib = 1, iesize
          call DVMCF2(wk_afft_l(1,ib),nsize,3,isin,ierr)
       end do
    end if
#else
    if(kimg==1) then
       NFFTW3(1) = fft_box_size_CD_3D(1,1)
       NEREAL(1) = fft_box_size_CD_3D(1,1)
       NEMBED(1) = fft_box_size_CD_3D(1,1)/2
       call dfftw_plan_many_dft_r2c(plan,           1, NFFTW3, ny*nz,  &
      &                             wk_afft_l, NEREAL, 1,         nx,  &
      &                             wk_afft_l, NEMBED, 1,       nx/2,  &
      &                             FFTW_ESTIMATE)
    else
       NFFTW3(1) = fft_box_size_CD_3D(1,1)
       NEMBED(1) = fft_box_size_CD_3D(1,0)
       call dfftw_plan_many_dft    (plan,           1, NFFTW3, ny*nz,  &
      &                             wk_afft_l, NEMBED,      1,    nx,  &
      &                             wk_afft_l, NEMBED,      1,    nx,  &
      &                             FFTW_FLAG, FFTW_ESTIMATE)
    endif
!OCL INDEPENDENT[dfftw_execute_dft_r2c,dfftw_execute_dft]
!OCL PARALLEL_STRONG
    do ib = 1, iesize
       if (kimg==1) then
          call dfftw_execute_dft_r2c(plan,wk_afft_l(1,ib),wk_afft_l(1,ib))
       else
          call dfftw_execute_dft    (plan,wk_afft_l(1,ib),wk_afft_l(1,ib))
       endif
    enddo
    call dfftw_destroy_plan(plan)
#endif

#ifdef __TIMER_DO__
  call timer_end(260)
#endif

#ifdef __TIMER_DO__
  call timer_sta(175)
#endif
    do l = 1, fftcd_Y_x_dim
       nis = nis_fftcd_Y_x(l)
       nie = nie_fftcd_Y_x(l)
       nnx = nie - nis + 1
       if (kimg == 1) then
          do ib = 1, iesize
             iadd0 = nnx*ny*nz*(ib-1)
             do k = 1, nz
                do j = 1, ny
                   do i = nis, nie
                      wk_send1(iadd0+(i-nis+1+(j-1)*nnx+(k-1)*nnx*ny),l) = wk_afft_l(i+(j-1)*nx+(k-1)*nx*ny,ib)
                   end do
                end do
             end do
          end do
       else
          do ib = 1, iesize
             iadd0 = nnx*ny*nz*(ib-1)*2
             do k = 1, nz
                do j = 1, ny
                   do i = nis, nie
                      iadd  = (i+(j-1)*nxp+(k-1)*nxp*ny)*2
                      iadd1 = iadd0+(i-nis+1+(j-1)*nnx+(k-1)*nnx*ny)*2
                      wk_send1(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send1(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end if
    end do
#ifdef __TIMER_DO__
  call timer_end(175)
#endif

#ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_fftcd_xy_world)
  call timer_sta(176)
#endif

    call MPI_ALLTOALL(wk_send1, x2y_smax*kimg*iesize, mpi_double_precision,   &
   &                  wk_recv1, x2y_rmax*kimg*iesize, mpi_double_precision,   &
   &                                                  mpi_fftcd_xy_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_CD_Inverse_XYZ_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 70, ierr)
     endif

#ifdef __TIMER_DO__
  call timer_sta(177)
#endif
#ifdef __TIMER_DO__
  call timer_end(177)
#endif

#else
!else ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_comm)
  call timer_sta(178)
#endif

    icnt_recv = 0
    lrank = myrank + 1
    if (lrank > (nmrank -1)) lrank = 0
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank -1)) lrank = 0
       if ((x2y_recv(1,lrank) /= 0)) then
          call mpi_irecv(wk_recv1(1,x2y_recv(2,lrank)), x2y_recv(1,lrank)*kimg*iesize, &
         &               mpi_double_precision, lrank, itag, mpi_comm, req_r(icnt_recv), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_CD_Inverse_XYZ_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 70, ierr)
           endif
          icnt_recv = icnt_recv + 1
       endif
    enddo

    icnt_send = 0
    lrank = myrank
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank -1)) lrank = 0
       if ((x2y_send(1,lrank) /= 0)) then
          call mpi_isend(wk_send1(1,x2y_send(2,lrank)), x2y_send(1,lrank)*kimg*iesize, &
                         mpi_double_precision, lrank, itag, mpi_comm, req_s(icnt_send), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_CD_Inverse_XYZ_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 71, ierr)
           endif
          icnt_send = icnt_send + 1
       endif
    enddo

#ifdef __TIMER_DO__
  call timer_sta(177)
#endif

#ifdef __TIMER_DO__
  call timer_end(177)
#endif

    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_CD_Inverse_XYZ_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 72, ierr)
     endif

    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_CD_Inverse_XYZ_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 73, ierr)
     endif

#endif
!endif ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_end(176)
#endif

#ifdef __TIMER_DO__
  call timer_sta(178)
#endif
    nx = xyz_fftcd_y(2,1)-xyz_fftcd_y(1,1)+1
    ny = xyz_fftcd_y(2,2)-xyz_fftcd_y(1,2)+1
    nz = xyz_fftcd_y(2,3)-xyz_fftcd_y(1,3)+1
    nyp = ny
! (x,y,z) -> (y,z,x)
    if (kimg == 1) then
       do ib = 1, iesize
          do i = 1, nx/2
             do k = 1, nz
                do l = 1, fftcd_X_y_dim
                   nis = nis_fftcd_X_y(l)
                   nie = nie_fftcd_X_y(l)
                   nny = nie - nis + 1
                   iadd0 = nx/2*nny*nz*(ib-1)*2
                   do j = nis, nie
                      iadd1 = iadd0 + (i+(j-nis)*nx/2+(k-1)*nx/2*nny)*2
                      iadd  = (j+(k-1)*nyp+(i-1)*nyp*nz)*2
                      wk_afft_l(iadd-1,ib) = wk_recv1(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv1(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end do
    else
       do ib = 1, iesize
          do i = 1, nx
             do k = 1, nz
                do l = 1, fftcd_X_y_dim
                   nis = nis_fftcd_X_y(l)
                   nie = nie_fftcd_X_y(l)
                   nny = nie - nis + 1
                   iadd0 = nx*nny*nz*(ib-1)*2
                   do j = nis, nie
                      iadd1 = iadd0+(i+(j-nis)*nx+(k-1)*nx*nny)*2
                      iadd  = (j+(k-1)*nyp+(i-1)*nyp*nz)*2
                      wk_afft_l(iadd-1,ib) = wk_recv1(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv1(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end do
    end if
#ifdef __TIMER_DO__
  call timer_end(178)
#endif
!
! Y-axis (z-x div)
!
#ifdef __TIMER_DO__
  call timer_sta(261)
#endif

#ifdef FFT_USE_SSL2
    if (kimg == 1) then
       nsize(1:3) = (/ny,nz,(nx/2)/)
       isin(1:3)  = (/1,0,0/)
    else
       nsize(1:3) = (/ny,nz,nx/)
       isin(1:3)  = (/-1,0,0/)
    end if
!OCL INDEPENDENT[ DVMCF2 ]
!OCL PARALLEL_STRONG
    do ib = 1, iesize
       call DVMCF2(wk_afft_l(1,ib),nsize,3,isin,ierr)
    end do
#else
    if(kimg==1) then
       NFFTW3(1) = fft_box_size_CD_3D(2,1)
       NEMBED(1) = fft_box_size_CD_3D(2,1)
       call dfftw_plan_many_dft    (plan,           1, NFFTW3, nz*nx/2, &
      &                             wk_afft_l, NEMBED,      1,      ny, &
      &                             wk_afft_l, NEMBED,      1,      ny, &
      &                             FFTW_FLAG, FFTW_ESTIMATE)
    else
       NFFTW3(1) = fft_box_size_CD_3D(2,1)
       NEMBED(1) = fft_box_size_CD_3D(2,1)
       lx = fft_box_size_CD_3D(1,1)
       call dfftw_plan_many_dft    (plan,           1, NFFTW3, nz*nx,  &
      &                             wk_afft_l, NEMBED,      1,    ny,  &
      &                             wk_afft_l, NEMBED,      1,    ny,  &
      &                             FFTW_FLAG, FFTW_ESTIMATE)
    end if
!OCL INDEPENDENT[dfftw_execute_dft]
!OCL PARALLEL_STRONG
    do ib = 1,iesize
       call dfftw_execute_dft(plan,wk_afft_l(1,ib),wk_afft_l(1,ib))
    end do
    call dfftw_destroy_plan(plan)
#endif

#ifdef __TIMER_DO__
  call timer_end(261)
#endif

#ifdef __TIMER_DO__
  call timer_sta(179)
#endif
#ifdef __TIMER_DO__
  call timer_end(179)
#endif

#ifdef __TIMER_DO__
  call timer_sta(180)
#endif
! (y,z,x) -> (y,z,x)
    if (kimg == 1) then
       do ib = 1, iesize
          do i = 1, nx/2
             do k = 1, nz
                do l = 1, fftcd_Z_y_dim
                   nis = nis_fftcd_Z_y(l)
                   nie = nie_fftcd_Z_y(l)
                   nny = nie - nis + 1
                   iadd0 = nx/2*nny*nz*(ib-1)*2
                   do j = nis, nie
                      iadd  = (j+(k-1)*nyp+(i-1)*nyp*nz)*2
                      iadd1 = iadd0 + (j-nis+1+(k-1)*nny+(i-1)*nny*nz)*2
                      wk_send2(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send2(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
    else
!OCL PARALLEL_STRONG
!OCL NORECURRENCE
       do ib = 1, iesize
!OCL SERIAL
          do i = 1, nx
!OCL SERIAL
             do k = 1, nz
!OCL SERIAL
                do l = 1, fftcd_Z_y_dim
                   nis = nis_fftcd_Z_y(l)
                   nie = nie_fftcd_Z_y(l)
                   nny = nie - nis + 1
                   iadd0 = nx*nny*nz*(ib-1)*2
!OCL SERIAL
                   do j = nis, nie
                      iadd  = (j+(k-1)*nyp+(i-1)*nyp*nz)*2
                      iadd1 = iadd0+(j-nis+1+(k-1)*nny+(i-1)*nny*nz)*2
                      wk_send2(iadd1-1,l) = wk_afft_l(iadd-1,ib)
                      wk_send2(iadd1  ,l) = wk_afft_l(iadd  ,ib)
                   end do
                end do
             end do
          end do
       end do
    end if

#ifdef __TIMER_DO__
  call timer_end(180)
#endif

#ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_fftcd_yz_world)
  call timer_sta(181)
#endif

    call MPI_ALLTOALL(wk_send2, y2z_smax*kimg*iesize, mpi_double_precision,   &
   &                  wk_recv2, y2z_rmax*kimg*iesize, mpi_double_precision,   &
   &                                                  mpi_fftcd_yz_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_CD_Inverse_XYZ_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 76, ierr)
     endif

#ifdef __TIMER_DO__
  call timer_sta(182)
#endif
#ifdef __TIMER_DO__
  call timer_end(182)
#endif

#else
!else ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_barrier(mpi_comm)
  call timer_sta(181)
#endif

    icnt_recv = 0
    lrank = myrank + 1
    if (lrank > (nmrank-1)) lrank = 0
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank-1)) lrank = 0
       if ((z2y_recv(1,lrank) /= 0)) then
          call mpi_irecv(wk_recv2(1,z2y_recv(2,lrank)), z2y_recv(1,lrank)*kimg*iesize, &
                         mpi_double_precision, lrank, itag, mpi_comm, req_r(icnt_recv), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_CD_Inverse_XYZ_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 76, ierr)
           endif
          icnt_recv = icnt_recv + 1
       endif
    enddo

    icnt_send = 0
    lrank = myrank
    do i = 0, nmrank - 1
       lrank = lrank + 1
       if (lrank > (nmrank -1)) lrank = 0
       if ((z2y_send(1,lrank) /= 0)) then
          call mpi_isend(wk_send2(1,z2y_send(2,lrank)), z2y_send(1,lrank)*kimg*iesize, &
                         mpi_double_precision, lrank, itag, mpi_comm, req_s(icnt_send), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_FFT_CD_Inverse_XYZ_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 77, ierr)
           endif
          icnt_send = icnt_send + 1
       endif
    enddo

#ifdef __TIMER_DO__
  call timer_sta(182)
#endif
#ifdef __TIMER_DO__
  call timer_end(182)
#endif

    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_CD_Inverse_XYZ_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 78, ierr)
     endif

    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_FFT_CD_Inverse_XYZ_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 79, ierr)
     endif

#endif
!endif ifdef FFT_ALLTOALL

#ifdef __TIMER_COMM__
  call timer_end(181)
#endif

#ifdef __TIMER_DO__
  call timer_sta(183)
#endif
    nx = xyz_fftcd_z(2,1)-xyz_fftcd_z(1,1)+1
    ny = xyz_fftcd_z(2,2)-xyz_fftcd_z(1,2)+1
    nz = xyz_fftcd_z(2,3)-xyz_fftcd_z(1,3)+1
    nzp = nz

! (y,z,x) -> (z,x,y)
    if (kimg == 1) then
       do ib = 1, iesize
          do j = 1, ny
             do i = 1, nx/2
                do l = 1, fftcd_Y_z_dim
                   nis = nis_fftcd_Y_z(l)
                   nie = nie_fftcd_Y_z(l)
                   nnz = nie - nis + 1
                   iadd0 = nx/2*ny*nnz*(ib-1)*2
                   do k = nis, nie
                      iadd1 = iadd0 + (j+(k-nis)*ny+(i-1)*ny*nnz)*2
                      iadd  = (k+(i-1)*nzp+(j-1)*nzp*nx/2)*2
                      wk_afft_l(iadd-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end do
    else
       do ib = 1, iesize
          do j = 1, ny
             do i = 1, nx
                do l = 1, fftcd_Y_z_dim
                   nis = nis_fftcd_Y_z(l)
                   nie = nie_fftcd_Y_z(l)
                   nnz = nie - nis + 1
                   iadd0 = nx*ny*nnz*(ib-1)*2
                   do k = nis, nie
                      iadd1 = iadd0+(j+(k-nis)*ny+(i-1)*ny*nnz)*2
                      iadd  = (k+(i-1)*nzp+(j-1)*nzp*nx)*2
                      wk_afft_l(iadd-1,ib) = wk_recv2(iadd1-1,l)
                      wk_afft_l(iadd  ,ib) = wk_recv2(iadd1  ,l)
                   end do
                end do
             end do
          end do
       end do
    end if
#ifdef __TIMER_DO__
  call timer_end(183)
#endif
!
! Z-axis (x-y div)
!
#ifdef __TIMER_DO__
  call timer_sta(262)
#endif

#ifdef FFT_USE_SSL2
    if (kimg == 1) then
       nsize(1:3) = (/nz,nx/2,ny/)
       isin(1:3)  = (/1,0,0/)
    else
       nsize(1:3) = (/nz,nx,ny/)
       isin(1:3)  = (/-1,0,0/)
    end if
!OCL INDEPENDENT[ DVMCF2 ]
!OCL PARALLEL_STRONG
    do ib = 1, iesize
       call DVMCF2(wk_afft_l(1,ib),nsize,3,isin,ierr)
    end do
#else
    if(kimg==1) then
       NFFTW3(1) = fft_box_size_CD_3D(3,1)
       NEMBED(1) = fft_box_size_CD_3D(3,0)
       call dfftw_plan_many_dft    (plan,           1, NFFTW3, nx*ny/2, &
      &                             wk_afft_l, NEMBED,      1,      nz, &
      &                             wk_afft_l, NEMBED,      1,      nz, &
      &                             FFTW_FLAG, FFTW_ESTIMATE)
    else
       NFFTW3(1) = fft_box_size_CD_3D(3,1)
       NEMBED(1) = fft_box_size_CD_3D(3,0)
       call dfftw_plan_many_dft    (plan,           1, NFFTW3, nx*ny,  &
      &                             wk_afft_l, NEMBED,      1,    nz,  &
      &                             wk_afft_l, NEMBED,      1,    nz,  &
      &                             FFTW_FLAG, FFTW_ESTIMATE)
    endif
!OCL INDEPENDENT[dfftw_execute_dft]
!OCL PARALLEL_STRONG
    do ib = 1, iesize
       call dfftw_execute_dft(plan,wk_afft_l(1,ib),wk_afft_l(1,ib))
    end do
    call dfftw_destroy_plan(plan)
#endif

#ifdef __TIMER_DO__
  call timer_end(262)
#endif

#ifdef __TIMER_DO__
  call timer_sta(184)
#endif
#ifdef __TIMER_DO__
  call timer_end(184)
#endif

1000 continue

#ifdef __TIMER_SUB__
    call timer_end(110)
#endif
  end subroutine m_FFT_CD_Inverse_XYZ_3D
!------------------------------------------------------------------------------

  subroutine m_FFT_WF_XYZ_3D(electron_or_positron,nfout,afft_l,lsize,iesize,inverse_or_direct)
    integer, intent(in) :: electron_or_positron, nfout, lsize, iesize, inverse_or_direct
    real(kind=DP), dimension(lsize*kimg,iesize), intent(inout) :: afft_l
    integer               :: nsize, ierr

#ifdef __TIMER_SUB__
    call timer_sta(104)
#endif
    if(inverse_or_direct == DIRECT) then
       call m_FFT_Direct_XYZ_3D (nfout, afft_l, lsize, iesize)
    else ! INVERSE
       call m_FFT_Inverse_XYZ_3D(nfout, afft_l, lsize, iesize)
    end if
#ifdef __TIMER_SUB__
    call timer_end(104)
#endif

  end subroutine m_FFT_WF_XYZ_3D
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  subroutine m_FFT_W_Vlocal_W_XYZ_3D(electron_or_positron,afft,bfft,lsize,ibesize,s)
    integer, intent(in) ::                             electron_or_positron
    integer, intent(in) ::                             lsize, ibesize
    real(kind=DP), intent(in),    dimension(lsize*kimg) :: afft
    real(kind=DP), intent(inout), dimension(lsize*kimg,ibesize) :: bfft
    real(kind=DP), intent(out),   dimension(ibesize)    :: s

    integer i, id, j, k, ip, nl, nm, nn, nd2, nd3, nlh, idh, ib, ierr, nx, ny, nz
    real(kind=DP) :: prd
#ifdef __TIMER_SUB__
    call timer_sta(108)
#endif

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

    if (np_fft_z == 0) then
       nx = 0
       ny = 0
       nz = 0
    else
       nx = xyz_fft_z(2,1)-xyz_fft_z(1,1)+1
       ny = xyz_fft_z(2,2)-xyz_fft_z(1,2)+1
       nz = xyz_fft_z(2,3)-xyz_fft_z(1,3)+1
    endif

    s = 0.d0
    if(kimg == 1) then

       do ib = 1, ibesize
#ifdef __TIMER_DO__
  call timer_sta(152)
#endif
          do i = 1, nx*nz*ny-1, 2
             bfft(i,ib) = afft(i)*(bfft(i,ib)**2 + bfft(i+1,ib)**2)
          end do
#ifdef __TIMER_DO__
  call timer_end(152)
#endif
#ifdef __TIMER_DO__
  call timer_sta(153)
#endif
          do i = 1, nx*nz*ny-1, 2
             s(ib) = s(ib) + bfft(i,ib)
          end do
#ifdef __TIMER_DO__
  call timer_end(153)
#endif
          s(ib) = s(ib) + s(ib)
#ifdef __TIMER_DO__
  call timer_sta(154)
#endif
          if(xyz_fft_z(1,1) == 1) then
#ifdef FFTW_STRIDE
             do k = 1, nz
                do j = 1, ny
                   do i = 1, 1
                      s(ib) = s(ib) - bfft((i+(j-1)*nx/2+(k-1)*nx*ny/2)*2-1,ib)
                   end do
                end do
             end do
#else
             do j = 1, ny
                do i = 1, 1
                   do k = 1, nz
                      s(ib) = s(ib) - bfft((k+(i-1)*nz+(j-1)*nz*nx/2)*2-1,ib)
                   end do
                end do
             end do
#endif
          endif
#ifdef __TIMER_DO__
  call timer_end(154)
#endif
#ifdef __TIMER_DO__
  call timer_sta(155)
#endif
          if(xyz_fft_z(2,1) == id) then
#ifdef FFTW_STRIDE
             do k = 1, nz
                do j = 1, ny
                   do i = nx/2, nx/2
                      s(ib) = s(ib) - bfft((i+(j-1)*nx/2+(k-1)*nx*ny/2)*2-1,ib)
                   end do
                end do
             end do
#else
             do j = 1, ny
                do i = nx/2, nx/2
                   do k = 1, nz
                      s(ib) = s(ib) - bfft((k+(i-1)*nz+(j-1)*nz*nx/2)*2-1,ib)
                   end do
                end do
             end do
#endif
          endif
#ifdef __TIMER_DO__
  call timer_end(155)
#endif
       end do

    else if(kimg == 2) then
#ifdef __TIMER_DO__
  call timer_sta(156)
#endif
       do ib = 1, ibesize
          do i = 1, nx*nz*ny*kimg-1, 2
             s(ib) = s(ib) + afft(i)*(bfft(i,ib)**2 + bfft(i+1,ib)**2)
          end do
       end do
#ifdef __TIMER_DO__
  call timer_end(156)
#endif
    end if
#ifdef __TIMER_SUB__
    call timer_end(108)
#endif
  end subroutine m_FFT_W_Vlocal_W_XYZ_3D
!------------------------------------------------------------------------------
#ifdef MPI_FFTW
!------------------------------------------------------------------------------
  subroutine m_FFT_W_Vlocal_W_mpifftw3d(lx,local_n,lz,afft,ibesize,s)
    integer(C_INTPTR_T), intent(in)                        :: lx,local_n,lz
    integer, intent(in)                                    :: ibesize
    real(kind=DP), intent(in),    dimension(lx,lz,local_n) :: afft
    real(kind=DP), intent(out),   dimension(ibesize)       :: s

    real(kind=DP) :: prd
    integer :: i1,i2,mm,ib
    integer :: id_sname=-1
    integer :: ix,iy,iz,lxh
    real(kind=DP), allocatable, dimension(:) :: bfft
    call tstatc0_begin('m_FFT_W_Vlocal_W_mpifftw ',id_sname)
    if(kimg==2) then
      s(1) = 0.d0
      do iy=1,local_n
        do iz=1,lz
          do ix=1,lx
!            i1 = (iy-1)*lx*lz+(iz-1)*lx+ix
            s(1) = s(1) + afft(ix,iz,iy)*(real(afft_mpifftw(ix,iz,iy))**2+aimag(afft_mpifftw(ix,iz,iy))**2)
          enddo
        enddo
      enddo
    else
      lxh = lx/2
      allocate(bfft(local_n*lz*lxh));bfft=0.d0
      do iy=1,local_n
        do iz=1,lz
          do ix=1,lxh
            i1 = (iy-1)*lxh*lz+(iz-1)*lxh+ix
            bfft(i1) = afft(ix,iz,iy)*(real(afft_mpifftw_kimg1(ix,iz,iy))**2+aimag(afft_mpifftw_kimg1(ix,iz,iy))**2)
          enddo
        enddo
      enddo
      s(1) = 0.d0
      do iy=1,local_n
        do iz=1,lz
          do ix=1,lxh
            i1 = (iy-1)*lxh*lz+(iz-1)*lxh+ix
            s(1) = s(1) + bfft(i1)
          enddo
        enddo
      enddo
      s(1) = s(1)+s(1)
      do iy=1,local_n
        do iz=1,lz
          ix = 1
          i1 = (iy-1)*lxh*lz+(iz-1)*lxh+ix
          s(1) = s(1) - bfft(i1)
          ix = lxh
          i1 = (iy-1)*lxh*lz+(iz-1)*lxh+ix
          s(1) = s(1) - bfft(i1)
        enddo
      enddo
      deallocate(bfft)
    endif
    call tstatc0_end(id_sname)
  end subroutine m_FFT_W_Vlocal_W_mpifftw3d

  subroutine m_FFT_W_Vlocal_W_mpifftw(afft,lsize,ibesize,s)
    integer, intent(in) ::                             lsize, ibesize
    real(kind=DP), intent(in),    dimension(lsize*kimg) :: afft
    real(kind=DP), intent(out),   dimension(ibesize)    :: s

    real(kind=DP) :: prd
    integer :: i1,i2,mm,ib
    integer(C_INTPTR_T)  :: local_n, local_n_offset, alloc_local, lx, ly, lz, mx,my,mz, lxh,mmy
    integer :: id_sname=-1
    integer :: ix,iy,iz
    real(kind=DP), allocatable, dimension(:) :: bfft
    call tstatc0_begin('m_FFT_W_Vlocal_W_mpifftw ',id_sname)
    lx = fft_box_size_WF(1,0)
    ly = fft_box_size_WF(2,0)
    lz = fft_box_size_WF(3,0)
    if(kimg==2) then
      alloc_local = fftw_mpi_local_size_3d(ly,lz,lx,mpi_ke_world,local_n,local_n_offset)
      s(1) = 0.d0
      do iy=1,local_n
        do iz=1,lz
          do ix=1,lx
            i1 = (iy-1)*lx*lz+(iz-1)*lx+ix
            s(1) = s(1) + afft(2*i1-1)*(real(afft_mpifftw(ix,iz,iy))**2+aimag(afft_mpifftw(ix,iz,iy))**2)
          enddo
        enddo
      enddo
    else
      lxh = lx/2
      alloc_local = fftw_mpi_local_size_3d(ly,lz,lxh,mpi_ke_world,local_n,local_n_offset)
      allocate(bfft(local_n*lz*lxh));bfft=0.d0
      do iy=1,local_n
        do iz=1,lz
          do ix=1,lxh
            i1 = (iy-1)*lxh*lz+(iz-1)*lxh+ix
            bfft(i1) = afft(2*i1-1)*(real(afft_mpifftw_kimg1(ix,iz,iy))**2+aimag(afft_mpifftw_kimg1(ix,iz,iy))**2)
          enddo
        enddo
      enddo
      s(1) = 0.d0
      do iy=1,local_n
        do iz=1,lz
          do ix=1,lxh
            i1 = (iy-1)*lxh*lz+(iz-1)*lxh+ix
            s(1) = s(1) + bfft(i1)
          enddo
        enddo
      enddo
      s(1) = s(1)+s(1)
      do iy=1,local_n
        do iz=1,lz
          ix = 1
          i1 = (iy-1)*lxh*lz+(iz-1)*lxh+ix
          s(1) = s(1) - bfft(i1)
          ix = lxh
          i1 = (iy-1)*lxh*lz+(iz-1)*lxh+ix
          s(1) = s(1) - bfft(i1)
        enddo
      enddo
      deallocate(bfft)
    endif
    call tstatc0_end(id_sname)
  end subroutine m_FFT_W_Vlocal_W_mpifftw
!------------------------------------------------------------------------------
#endif

#ifdef FFT_3D_DIVISION
!------------------------------------------------------------------------------
  subroutine m_FFT_Direct_3DIV_3D(nfout,wk_afft_l, wk_size, iesize)
    integer, intent(in) :: nfout, wk_size, iesize
    real(kind=DP), dimension(wk_size*2,iesize), intent(inout) :: wk_afft_l

    integer :: nfft1, nfft2, nfft3, nx, ny, nz, ib, i, j, k, iadd, is, iixx
    integer :: kx1p, kx2p, kx3p, nw
    complex*16, allocatable, dimension(:),save     :: work

#ifdef __TIMER_SUB__
    call mpi_barrier(mpi_ke_world,ierr)
    call timer_sta(105)
#endif
       is = 1

       nfft1 = fft_box_size_WF(1,1)
       nfft2 = fft_box_size_WF(2,1)
       nfft3 = fft_box_size_WF(3,1)

       kx1p = fft_X_x_nel
       kx2p = fft_X_y_nel
       kx3p = fft_X_z_nel
       nw = kx1p*kx2p*kx3p
       if (firstcall_direct_3ddiv) then
          allocate(work(nw), stat=ierr)
           if (ierr /= 0) then
              write(nfout,'("*** ERROR alloc m_FFT_Direct_3DIV_3D : error code -> ",i8)') ierr
           end if
          firstcall_direct_3ddiv = .false.
       end if
!
       if (np_fft_x > 0) then
         nx = xyz_fft_x(2,1)-xyz_fft_x(1,1)+1
         ny = xyz_fft_x(2,2)-xyz_fft_x(1,2)+1
         nz = xyz_fft_x(2,3)-xyz_fft_x(1,3)+1
       else
         nx = 0
         ny = 0
         nz = 0
       end if

       do ib = 1, iesize
#ifdef __TIMER_DO__
  call timer_sta(252)
#endif
          call DS_V3DCFT3(wk_afft_l(1,ib), kx1p, kx2p, kx1p, kx2p, kx3p, nfft1, nfft2, nfft3,  &
         &                fft_X_x_dim, fft_X_y_dim, fft_X_z_dim, work, nw, is, mpi_fft_zy_world, ierr)
          if (ierr /= 0) then
             write(nfout,'("*** ERROR DS_V3DCFT3 SUBOUTINE : error code -> ",i8)') ierr
          end if
#ifdef __TIMER_DO__
  call timer_end(252)
#endif
       end do

#ifdef __TIMER_SUB__
    call timer_end(105)
#endif
  end subroutine m_FFT_Direct_3DIV_3D
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  subroutine m_FFT_Inverse_3DIV_3D(nfout,wk_afft_l, wk_size, iesize)
    integer, intent(in) :: nfout, wk_size, iesize
    real(kind=DP), dimension(wk_size*2,iesize), intent(inout) :: wk_afft_l

    integer :: nx, ny, nz, i, j, k, ib, iadd, is, iixx
    integer :: nfft1, nfft2, nfft3
    integer :: kx1p, kx2p, kx3p, kx1, kx2, nw
    complex*16, allocatable, dimension(:),save     :: work

#ifdef __TIMER_SUB__
    call timer_barrier(mpi_ke_world)
    call timer_sta(106)
#endif
       is = -1

       nfft1 = fft_box_size_WF(1,1)
       nfft2 = fft_box_size_WF(2,1)
       nfft3 = fft_box_size_WF(3,1)

       kx1p = fft_X_x_nel
       kx2p = fft_X_y_nel
       kx3p = fft_X_z_nel
       nw = kx1p*kx2p*kx3p
       if (firstcall_inverse_3ddiv) then
          allocate(work(nw), stat=ierr)
           if (ierr /= 0) then
              write(nfout,'("*** ERROR alloc m_FFT_Inverse_3DIV_3D : error code -> ",i8)') ierr
           end if
          firstcall_inverse_3ddiv = .false.
       end if

       if(np_fft_x > 0) then
         nx = xyz_fft_x(2,1)-xyz_fft_x(1,1)+1
         ny = xyz_fft_x(2,2)-xyz_fft_x(1,2)+1
         nz = xyz_fft_x(2,3)-xyz_fft_x(1,3)+1
       else
         nx = 0
         ny = 0
         nz = 0
       end if

       do ib = 1, iesize
#ifdef __TIMER_DO__
  call timer_sta(255)
#endif
          call DS_V3DCFT3(wk_afft_l(1,ib), kx1p, kx2p, kx1p, kx2p, kx3p, nfft1, nfft2, nfft3,  &
         &                fft_X_x_dim, fft_X_y_dim, fft_X_z_dim, work, nw, is, mpi_fft_zy_world, ierr)
          if (ierr /= 0) then
             write(nfout,'("*** ERROR DS_V3DCFT3 SUBOUTINE : error code -> ",i8)') ierr
          end if
#ifdef __TIMER_DO__
  call timer_end(255)
#endif
       end do

#ifdef __TIMER_SUB__
    call timer_end(106)
#endif
  end subroutine m_FFT_Inverse_3DIV_3D
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  subroutine m_FFT_WF_3DIV_3D(electron_or_positron,nfout,afft_l,lsize,iesize,inverse_or_direct)
    integer, intent(in) :: electron_or_positron, nfout, lsize, iesize, inverse_or_direct
    real(kind=DP), dimension(lsize*2,iesize), intent(inout) :: afft_l
    integer               :: nsize, ierr

#ifdef __TIMER_SUB__
    call timer_sta(104)
#endif
    if(inverse_or_direct == DIRECT) then
       call m_FFT_Direct_3DIV_3D (nfout, afft_l, lsize, iesize)
    else ! INVERSE
       call m_FFT_Inverse_3DIV_3D(nfout, afft_l, lsize, iesize)
    end if
#ifdef __TIMER_SUB__
    call timer_end(104)
#endif

  end subroutine m_FFT_WF_3DIV_3D
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  subroutine m_FFT_Vlocal_W_3DIV_3D(afft_l,bfft_l,lsize,ibsize,nfft_l)
    integer, intent(in) :: lsize, ibsize, nfft_l
    real(kind=DP), intent(in)   , dimension(lsize*2)        :: afft_l
    real(kind=DP), intent(inout), dimension(lsize*2,ibsize) :: bfft_l
    integer i, j, k, nx, ny, nz, kx1p, kx2p, kx3p, iadd, ib
#ifdef __TIMER_SUB__
    call timer_sta(107)
#endif

    if (np_fft_x == 0) then
       nx = 0
       ny = 0
       nz = 0
    else
       nx = xyz_fft_x(2,1)-xyz_fft_x(1,1)+1
       ny = xyz_fft_x(2,2)-xyz_fft_x(1,2)+1
       nz = xyz_fft_x(2,3)-xyz_fft_x(1,3)+1
    endif

    kx1p = fft_X_x_nel
    kx2p = fft_X_y_nel
    kx3p = fft_X_z_nel

#ifdef __TIMER_DO__
  call timer_sta(151)
#endif
!   do ib = 1, ibsize
!      do k = 1, nz
!        do j = 1, ny
!          do i = 1, nx*2-1, 2
!            iadd = i+kx1p*2*(j-1)+kx1p*2*kx2p*(k-1)
!            bfft_l(iadd  ,ib) = afft_l(iadd  ) * bfft_l(iadd  ,ib)
!            bfft_l(iadd+1,ib) = afft_l(iadd  ) * bfft_l(iadd+1,ib)
!          end do
!        end do
!      end do
!   end do
    do ib = 1, ibsize
       do k = 1, nz
         do j = 1, ny
           do i = 1, nx
             iadd = i+kx1p*(j-1)+kx1p*kx2p*(k-1)
             bfft_l(iadd*2-1,ib) = afft_l(iadd*2-1) * bfft_l(iadd*2-1,ib)
             bfft_l(iadd*2  ,ib) = afft_l(iadd*2-1) * bfft_l(iadd*2  ,ib)
           end do
         end do
       end do
    end do
!      do j = 1, ibsize
!         do i = 1, nfft_l
!            bfft_l(i*2-1,j) = afft_l(i*2-1) * bfft_l(i*2-1,j)
!            bfft_l(i*2  ,j) = afft_l(i*2-1) * bfft_l(i*2  ,j)
!         end do
!      end do
#ifdef __TIMER_DO__
  call timer_end(151)
#endif
#ifdef __TIMER_SUB__
    call timer_end(107)
#endif
  end subroutine m_FFT_Vlocal_W_3DIV_3D
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  subroutine m_FFT_W_Vlocal_W_3DIV_3D(electron_or_positron,afft,bfft,lsize,ibesize,s)
    integer, intent(in) ::                             electron_or_positron
    integer, intent(in) ::                             lsize, ibesize
    real(kind=DP), intent(in),    dimension(lsize*2)         :: afft
    real(kind=DP), intent(inout), dimension(lsize*2,ibesize) :: bfft
    real(kind=DP), intent(out),   dimension(ibesize)    :: s

    integer i, id, j, k, ip, nl, nm, nn, nd2, nd3, nlh, idh, ib, ierr, nx, ny, nz
    integer kx1p, kx2p, kx3p, ladd
!   real(kind=DP),dimension(ibesize) :: s, s_mpi
    real(kind=DP) :: prd
#ifdef __TIMER_SUB__
    call timer_sta(108)
#endif

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

    if (np_fft_x == 0) then
       nx = 0
       ny = 0
       nz = 0
    else
       nx = xyz_fft_x(2,1)-xyz_fft_x(1,1)+1
       ny = xyz_fft_x(2,2)-xyz_fft_x(1,2)+1
       nz = xyz_fft_x(2,3)-xyz_fft_x(1,3)+1
    endif

    kx1p = fft_X_x_nel
    kx2p = fft_X_y_nel
    kx3p = fft_X_z_nel

    s = 0.d0
#ifdef __TIMER_DO__
  call timer_sta(156)
#endif
    do ib = 1, ibesize
       do k = 1, nz
         do j = 1, ny
           do i = 1, nx*2-1, 2
             ladd = i+kx1p*2*(j-1)+kx1p*2*kx2p*(k-1)
             s(ib) = s(ib) + afft(ladd)*(bfft(ladd,ib)**2 + bfft(ladd+1,ib)**2)
           end do
         end do
       end do
    end do
#ifdef __TIMER_DO__
  call timer_end(156)
#endif

#ifdef __TIMER_SUB__
    call timer_end(108)
#endif
  end subroutine m_FFT_W_Vlocal_W_3DIV_3D
!------------------------------------------------------------------------------
#endif
!endif #ifdef FFT_3D_DIVISION

#ifdef FFT_3D_DIVISION_CD
!------------------------------------------------------------------------------
  subroutine m_FFT_CD_Direct_3DIV_3D(nfout,wk_afft_l, wk_size, iesize)
    integer, intent(in) :: nfout, wk_size, iesize
    real(kind=DP), dimension(wk_size*2,iesize), intent(inout) :: wk_afft_l

    integer :: nfft1, nfft2, nfft3, nx, ny, nz, ib, i, j, k, iadd, is, iixx
    integer :: kx1p, kx2p, kx3p, nw
    complex*16, allocatable, dimension(:,:,:) :: wk_fft
    complex*16, allocatable, dimension(:),save     :: work

#ifdef __TIMER_SUB__
    call mpi_barrier(mpi_ke_world,ierr)
    call timer_sta(109)
#endif
       is = 1

       if (kimg == 1) then
          nfft1 = fft_box_size_CD_3D(1,1)
          nfft2 = fft_box_size_CD_3D(2,1)
          nfft3 = fft_box_size_CD_3D(3,1)
       else
          nfft1 = fft_box_size_CD_3D(1,0)
          nfft2 = fft_box_size_CD_3D(2,0)
          nfft3 = fft_box_size_CD_3D(3,0)
       end if

       kx1p = fftcd_X_x_nel
       kx2p = fftcd_X_y_nel
       kx3p = fftcd_X_z_nel
       nw = kx1p*kx2p*kx3p
       if (firstcall_cd_direct_3ddiv) then
          allocate(work(nw), stat=ierr)
          firstcall_cd_direct_3ddiv = .false.
       end if
!
       if (np_fftcd_x > 0) then
         nx = xyz_fftcd_x(2,1)-xyz_fftcd_x(1,1)+1
         ny = xyz_fftcd_x(2,2)-xyz_fftcd_x(1,2)+1
         nz = xyz_fftcd_x(2,3)-xyz_fftcd_x(1,3)+1
       else
         nx = 0
         ny = 0
         nz = 0
       end if

       do ib = 1, iesize
#ifdef __TIMER_DO__
  call timer_sta(256)
#endif
#ifdef __TIMER_DO__
  call timer_end(256)
#endif
#ifdef __TIMER_DO__
  call timer_sta(257)
#endif
          call DS_V3DCFT3(wk_afft_l(1,ib), kx1p, kx2p, kx1p, kx2p, kx3p, nfft1, nfft2, nfft3,  &
         &                fftcd_X_x_dim, fftcd_X_y_dim, fftcd_X_z_dim, work, nw, is,           &
         &                mpi_fftcd_zy_world, ierr)
          if (ierr /= 0) then
             write(nfout,'("*** ERROR DS_V3DCFT3 SUBOUTINE : error code -> ",i8)') ierr
          end if
#ifdef __TIMER_DO__
  call timer_end(257)
#endif
#ifdef __TIMER_DO__
  call timer_sta(258)
#endif
#ifdef __TIMER_DO__
  call timer_sta(258)
#endif
       end do

#ifdef __TIMER_SUB__
    call timer_end(109)
#endif
  end subroutine m_FFT_CD_Direct_3DIV_3D
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  subroutine m_FFT_CD_Inverse_3DIV_3D(nfout,wk_afft_l, wk_size, iesize)
    integer, intent(in) :: nfout, wk_size, iesize
    real(kind=DP), dimension(wk_size*2,iesize), intent(inout) :: wk_afft_l

    integer :: nx, ny, nz, i, j, k, ib, iadd, is, iixx
    integer :: nfft1, nfft2, nfft3
    integer :: kx1p, kx2p, kx3p, kx1, kx2, nw
    complex*16, allocatable, dimension(:,:,:) :: wk_fft
    complex*16, allocatable, dimension(:),save     :: work

#ifdef __TIMER_SUB__
    call timer_barrier(mpi_ke_world)
    call timer_sta(110)
#endif
       is = -1

       if (kimg == 1) then
          nfft1 = fft_box_size_CD_3D(1,1)
          nfft2 = fft_box_size_CD_3D(2,1)
          nfft3 = fft_box_size_CD_3D(3,1)
       else
          nfft1 = fft_box_size_CD_3D(1,0)
          nfft2 = fft_box_size_CD_3D(2,0)
          nfft3 = fft_box_size_CD_3D(3,0)
       end if

       kx1p = fftcd_X_x_nel
       kx2p = fftcd_X_y_nel
       kx3p = fftcd_X_z_nel
       nw = kx1p*kx2p*kx3p
       if (firstcall_cd_inverse_3ddiv) then
          allocate(work(nw), stat=ierr)
          firstcall_cd_inverse_3ddiv = .false.
       end if

       if(np_fftcd_x > 0) then
         nx = xyz_fftcd_x(2,1)-xyz_fftcd_x(1,1)+1
         ny = xyz_fftcd_x(2,2)-xyz_fftcd_x(1,2)+1
         nz = xyz_fftcd_x(2,3)-xyz_fftcd_x(1,3)+1
       else
         nx = 0
         ny = 0
         nz = 0
       end if

       do ib = 1, iesize
#ifdef __TIMER_DO__
  call timer_sta(260)
#endif
#ifdef __TIMER_DO__
  call timer_end(260)
#endif
#ifdef __TIMER_DO__
  call timer_sta(261)
#endif
          call DS_V3DCFT3(wk_afft_l(1,ib), kx1p, kx2p, kx1p, kx2p, kx3p, nfft1, nfft2, nfft3,  &
         &                fftcd_X_x_dim, fftcd_X_y_dim, fftcd_X_z_dim, work, nw, is,           &
         &                mpi_fftcd_zy_world, ierr)
          if (ierr /= 0) then
             write(nfout,'("*** ERROR DS_V3DCFT3 SUBOUTINE : error code -> ",i8)') ierr
          end if
#ifdef __TIMER_DO__
  call timer_end(261)
#endif
#ifdef __TIMER_DO__
  call timer_sta(262)
#endif
#ifdef __TIMER_DO__
  call timer_end(262)
#endif
       end do

#ifdef __TIMER_SUB__
    call timer_end(110)
#endif
  end subroutine m_FFT_CD_Inverse_3DIV_3D
!------------------------------------------------------------------------------
#endif
!endif #ifdef FFT_3D_DIVISION_CD

!#endif PARA3D
#endif

end module m_FFT
