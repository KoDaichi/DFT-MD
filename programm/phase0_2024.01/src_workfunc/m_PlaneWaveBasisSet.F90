!================================================
!  Software name : STM
!  Module : m_PlaneWaveBasisSet
!  Subroutine(s) : wd_planewavebasisset, m_pwBS_rd_data,
!                  m_pwBS_kintetic_energies,
!                  m_pwBS_set_FFT_mapfunctions
!  Author(s)     : Takahiro Yamasaki (June 7, 2004)
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

module m_PlaneWaveBasisSet
! $Id: m_PlaneWaveBasisSet.f90,v 1.1.1.1 2004/06/26 11:21:37 yamasaki Exp $
  use m_ArraySize_Parameters,only : kng,kng1,knv3,kngp
  use m_Const_Parameters,    only : DP, BEGIN, END, CARTS, BUCS, CRDTYP
  use m_Kpoints,             only : vkxyz
  use m_Control_Parameters,  only : n_fc, nc_z

  implicit none

!  integer, dimension(kngp,3)      :: ngabc
  integer, allocatable, dimension(:,:)      :: ngabc
!  real(kind=DP), dimension(kngp)  :: gr_l

!  integer, dimension(kng)         :: igf
!  integer, dimension(kngp)        :: igfp_l
!  integer, dimension(kngp)        :: igfpf_l
  integer, allocatable, dimension(:)         :: igf
  integer, allocatable, dimension(:)        :: igfp_l
  integer, allocatable, dimension(:)        :: igfpf_l

  integer                         :: kg, kgp
!!$  integer, dimension(3)           :: n_rGv
!!$  integer, dimension(3)           :: n_rGpv
  ! n_rG(p)v is the size of a rhombohedron which contains
  ! a sphere whose radius is 2*Gmax (Gmaxp) .

!  integer       :: nbase(kng1,knv3) ! PW Basis Set for each k-points
!  integer       :: iba(knv3)
  integer, allocatable, dimension(:,:)      :: nbase ! PW Basis Set for each k-points
  integer, allocatable, dimension(:)        :: iba

  real(kind=DP), dimension(6) :: ttr

contains

  subroutine wd_planewavebasisset_data(nfcntn_bin)
    integer, intent(in) :: nfcntn_bin
    write(nfcntn_bin) ngabc
    write(nfcntn_bin) igf
    write(nfcntn_bin) igfp_l
    write(nfcntn_bin) kg, kgp
    write(nfcntn_bin) nbase
    write(nfcntn_bin) iba
    write(nfcntn_bin) ttr
  end subroutine wd_planewavebasisset_data

  subroutine m_pwBS_rd_data(nfcntn_bin)
    integer, intent(in) :: nfcntn_bin

    allocate(ngabc(kngp,3))
    allocate(igf(kng))
    allocate(igfp_l(kngp))
    allocate(igfpf_l(kngp))
    allocate(nbase(kng1,knv3))
    allocate(iba(knv3))

    read(nfcntn_bin) ngabc
    read(nfcntn_bin) igf
    read(nfcntn_bin) igfp_l
    read(nfcntn_bin) kg, kgp
    read(nfcntn_bin) nbase
    read(nfcntn_bin) iba
    read(nfcntn_bin) ttr
 !   call calc_length_of_G
  end subroutine m_pwBS_rd_data
       
!!$  subroutine calc_length_of_G
!!$    integer :: i
!!$    do i = 1, kgp
!!$       gr_l(i) =      dsqrt(ttr(1)*ngabc(i,1)*ngabc(i,1) &
!!$            &             + ttr(2)*ngabc(i,2)*ngabc(i,2) &
!!$            &             + ttr(3)*ngabc(i,3)*ngabc(i,3) &
!!$            &             + ttr(4)*ngabc(i,1)*ngabc(i,2) &
!!$            &             + ttr(5)*ngabc(i,2)*ngabc(i,3) &
!!$            &             + ttr(6)*ngabc(i,3)*ngabc(i,1))
!!$    enddo
!!$  end subroutine calc_length_of_G

  subroutine m_pwBS_kinetic_energies(ik,vkxyz,ekin,nmp,nnp)
    integer, intent(in)       :: ik
    real(kind=DP), intent(in) :: vkxyz(knv3,3,CRDTYP)
    real(kind=DP), intent(out):: ekin(*) ! d(nmp*nnp)
    integer, intent(in)       :: nmp,nnp

    integer i, j, jj, k, kk
    real(kind=DP) :: gb,gc
    integer       :: id_sname = -1

    real*8 :: ttr2, ttr3, ttr5

    if(nc_z .eq. 1) then
       ttr2 = ttr(2); ttr3 = ttr(3); ttr5 = ttr(5)
    else if(nc_z .eq. 2) then
       ttr2 = ttr(1); ttr3 = ttr(3); ttr5 = ttr(6)
    else if(nc_z .eq. 3) then
       ttr2 = ttr(1); ttr3 = ttr(2); ttr5 = ttr(4)
    end if

    i = 0
    do j = 1, nmp
       jj = j
       if (jj > nmp/2 ) jj = jj -nmp
       do k = 1, nnp
          kk = k
          if ( kk > nnp/2 ) kk = kk -nnp
          gb = vkxyz(ik,n_fc(2),BUCS) + jj
          gc = vkxyz(ik,n_fc(3),BUCS) + kk
          i = i + 1
          ekin(i) = ( ttr2*gb*gb + ttr3*gc*gc  + ttr5*gb*gc )*0.5d0
!
       end do
    end do
  end subroutine m_pwBS_kinetic_energies

  subroutine m_pwBS_set_FFT_mapfunctions(fft_box_size_WF &
       &, fft_box_size_CD, fft_box_size_fine)
    integer, intent(in), dimension(0:3)    :: fft_box_size_WF
    integer, intent(in), dimension(0:3)    :: fft_box_size_CD
    integer, intent(in), dimension(0:3)    :: fft_box_size_fine
      
    integer id, i
    integer igf1, igf2, igf3
    integer               :: id_sname = -1

    id = fft_box_size_WF(0)
    do i = 1, kg
       igf1 = ngabc(i,1) + 1
       igf2 = ngabc(i,2) + 1
       igf3 = ngabc(i,3) + 1
       if(ngabc(i,1) <= -1) igf1 = igf1 + fft_box_size_WF(1)
       if(ngabc(i,2) <= -1) igf2 = igf2 + fft_box_size_WF(2)
       if(ngabc(i,3) <= -1) igf3 = igf3 + fft_box_size_WF(3)
       igf(i) = igf1 + (igf2-1)*id + (igf3-1)*id*fft_box_size_WF(2)
    enddo

    id = fft_box_size_CD(0)
    do i = 1, kgp
       igf1 = ngabc(i,1) + 1
       igf2 = ngabc(i,2) + 1
       igf3 = ngabc(i,3) + 1
       if(ngabc(i,1) <= -1) igf1 = igf1 + fft_box_size_CD(1)
       if(ngabc(i,2) <= -1) igf2 = igf2 + fft_box_size_CD(2)
       if(ngabc(i,3) <= -1) igf3 = igf3 + fft_box_size_CD(3)
       igfp_l(i) = igf1 + (igf2-1)*id + (igf3-1)*id*fft_box_size_CD(2)
    enddo

    id = fft_box_size_fine(0)
    do i = 1, kgp
!!$       igf1 = ngabc(i,1) + 1
!!$       igf2 = ngabc(i,2) + 1
!!$       igf3 = ngabc(i,3) + 1
       igf1 = ngabc(i,n_fc(1)) + 1
       igf2 = ngabc(i,n_fc(2)) + 1
       igf3 = ngabc(i,n_fc(3)) + 1
!!$       if(ngabc(i,1) <= -1) igf1 = igf1 + fft_box_size_fine(1)
!!$       if(ngabc(i,2) <= -1) igf2 = igf2 + fft_box_size_fine(2)
!!$       if(ngabc(i,3) <= -1) igf3 = igf3 + fft_box_size_fine(3)
       if(ngabc(i,n_fc(1)) <= -1) igf1 = igf1 + fft_box_size_fine(1)
       if(ngabc(i,n_fc(2)) <= -1) igf2 = igf2 + fft_box_size_fine(2)
       if(ngabc(i,n_fc(3)) <= -1) igf3 = igf3 + fft_box_size_fine(3)
       igfpf_l(i) = igf1 + (igf2-1)*id + (igf3-1)*id*fft_box_size_fine(2)
    enddo

  end subroutine m_pwBS_set_FFT_mapfunctions

end module m_PlaneWaveBasisSet







