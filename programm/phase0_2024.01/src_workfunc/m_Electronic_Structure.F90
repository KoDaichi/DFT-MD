!================================================
!  Software name : STM
!  Module : m_Electronic_Structure
!  Subroutine(s) : m_ES_rd_EigenValues_etc, wd_EigenValues_etc,
!                  m_ES_rd_VLCs, m_ES_VLCS_in_Rspace_fine_mesh,
!                  m_ES_rd_WFs, m_ES_WF_in_Rspace, m_ES_WF_in_Rspace_fine_mesh
!  Author(s)     : Takahiro Yamasaki and Koichi Kato (June 7, 2004)
!
!  FURTHER MODIFICATION: Junichiro  Koga (June 24, 2004)
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

module m_Electronic_Structure
! $Id: m_Electronic_Structure.f90,v 1.1.1.1 2004/06/26 11:21:37 yamasaki Exp $
  use m_ArraySize_Parameters, only : knv3, kngp, kimg, kspin, kng1 &
       &, keg
  use m_Const_Parameters,     only : SP, DP, INVERSE, ON, OFF
  use m_Files,                only : nfout,nfzaj, nfvlc
  use m_PlaneWaveBasisSet,    only : igf, igfp_l, igfpf_l, nbase, iba
  use m_FFT,                  only : FFT_WF, FFT_WF_fine &
       &                            ,nfft,nfftp,nfftpf,nlpf
  implicit none

!  real(kind=DP), dimension(kng1,keg,kimg)  :: zaj_l ! wave functions
!  real(kind=DP), dimension(kngp,kimg,kspin):: vlc_l ! local potentials
  real(kind=DP), allocatable, dimension(:,:,:)  :: zaj_l ! wave functions
  real(kind=DP), allocatable, dimension(:,:,:)  :: vlc_l ! local potentials

  real(kind=SP), pointer, dimension(:,:)   :: wf_l
  real(kind=DP), pointer, dimension(:,:)   :: vl_l
  ! A work array for reading of wave functions

!  integer, dimension(keg,knv3)        :: neordr, nrvf_ordr
!  real(DP),dimension(keg,knv3)        :: eko, occup_l
  integer, allocatable, dimension(:,:)        :: neordr, nrvf_ordr
  real(DP),allocatable, dimension(:,:)        :: eko, occup_l

  real(DP)                            :: efermi
  integer                             :: is, izaj

contains
  subroutine m_ES_rd_EigenValues_etc(nfcntn_bin)
    integer, intent(in) :: nfcntn_bin

    integer :: i, k

    allocate(zaj_l(kng1,keg,kimg))
    allocate(vlc_l(kngp,kimg,kspin))
    allocate(neordr(keg,knv3))
    allocate(nrvf_ordr(keg,knv3))
    allocate(eko(keg,knv3))
    allocate(occup_l(keg,knv3))

    read(nfcntn_bin) neordr
    read(nfcntn_bin) nrvf_ordr

    read(nfcntn_bin) eko

    read(nfcntn_bin) occup_l
    read(nfcntn_bin) efermi
#ifdef DEBUG
    write(0,*) 'neordr',neordr
    write(0,*) 'nrvf_ordr',nrvf_ordr
    write(0,*) 'eko',eko
    write(0,*) 'occup_l',occup_l
    write(0,*) 'efermi',efermi
#endif
  end subroutine m_ES_rd_EigenValues_etc

  subroutine wd_EigenValues_etc(nfcntn_bin)
    integer, intent(in) :: nfcntn_bin
    write(nfcntn_bin) neordr,nrvf_ordr
    write(nfcntn_bin) eko
    write(nfcntn_bin) occup_l
    write(nfcntn_bin) efermi
  end subroutine wd_EigenValues_etc

  subroutine m_ES_rd_VLCs(nfvlc)
    integer, intent(in) :: nfvlc
!!$    integer       :: ri,ia

!!$    rewind nfvlc
!!$    allocate(vl_l(kngp,kimg))
!!$    read(nfvlc) vl_l
!!$    print *,' read finished vl_l '
!!$    do ri = 1, kimg
!!$       do ia = 1, kngp
!!$          vlc_l(ia, ri) = vl_l(ia, ri)
!!$       end do
!!$    end do
!!$    deallocate(vl_l)
    read(nfvlc) vlc_l
  end subroutine m_ES_rd_VLCs

  subroutine m_ES_VLC_in_Rspace_fine_mesh(is,bfft)
    integer, intent(in)   :: is
    real(kind=DP), intent(inout), dimension(nfftpf) :: bfft

    integer :: i,i1,ri
    bfft = 0.d0
    do ri = 1, kimg
       do i = 1, kngp
          i1 = kimg*igfpf_l(i) + (ri - kimg)
          bfft(i1) = vlc_l(i,ri,is)
       end do
    end do
    !!write(nfout,*) 'VLC before FFT'
    !!write(nfout,'(2d20.8)') ( bfft(2*i-1), bfft(2*i), i=1, nlpf)
    call FFT_WF_fine(bfft,INVERSE)    !-(m_FFT)
    ! V_{loc}(G)) -> V_{loc}(r)
  end subroutine m_ES_VLC_in_Rspace_fine_mesh

  subroutine m_ES_rd_WFs(nfzaj)
    integer, intent(in) :: nfzaj
    integer       :: ib,ri

!!    print *,'(" !D kng1, keg = ",2i5," <<m_ES_rd_WFs>>")',kng1,keg
!!    print *,' !D Reading zaj'
!!$#ifdef VDB
!!$    allocate(wf_l(kng1,keg))
!!$    do ri = 1, kimg
!!$       read(nfzaj) wf_l
!!$       do ib = 1, keg
!!$          zaj_l(1:kng1,ib,ri) = wf_l(1:kng1,ib)
!!$       end do
!!$    end do
!!$#else             
    allocate(wf_l(kng1,kimg))
    do ib = 1, keg
       read(nfzaj) wf_l
       do ri = 1, kimg
          zaj_l(1:kng1,ib,1:kimg) = wf_l(1:kng1,1:kimg)
       end do
    end do
!!$#endif
    deallocate(wf_l)
!   write(nfout,*) 'WF after READ, kng1 =',kng1
!   write(nfout,*) ( zaj_l(ia,32,1), ia=1, 2000)
  end subroutine m_ES_rd_WFs

  subroutine m_ES_WF_in_Rspace(ikp,ib,bfft)
    integer, intent(in)                           :: ikp, ib
    real(kind=DP), intent(inout), dimension(nfft) :: bfft

    integer :: i,i1,ri
    bfft = 0.d0
    do ri = 1, kimg
       do i = 1, iba(ikp)
          i1 = kimg*igf(nbase(i,ikp)) + (ri - kimg)
          bfft(i1) = zaj_l(i,ib,ri)
       end do
    end do
    call FFT_WF(bfft,INVERSE,ON)    !-(m_FFT)
    ! \Psi(G) -> \Phi(r)
  end subroutine m_ES_WF_in_Rspace

  subroutine m_ES_WF_in_Rspace_fine_mesh(ikp,ib,afft)
    integer, intent(in)                             :: ikp, ib
    real(kind=DP), intent(inout), dimension(nfftpf) :: afft

    integer :: i,i1,ri
    afft = 0.d0
    do ri = 1, kimg
       do i = 1, iba(ikp)
          i1 = kimg*igfpf_l(nbase(i,ikp)) + (ri - kimg)
          afft(i1) = zaj_l(i,ib,ri)
       end do
    end do
!   write(nfout,*) 'WF before FFT, iba(ikp)=', iba(ikp)
!   write(nfout,*) ( afft(i), i=1, 100)
    call FFT_WF_fine(afft,INVERSE)    !-(m_FFT)
    ! \Psi(G) -> \Phi(r)
  end subroutine m_ES_WF_in_Rspace_fine_mesh

end module m_Electronic_Structure
