!================================================
!  Software name : STM
!  Module : m_Charge_File
!  Subroutine(s) : m_CF_rd_Header_For_Cube
!  Author(s)     : Junichiro Koga (June 24, 2004)
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

module m_Charge_File
! $Id: m_Charge_File.F90,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $
  use m_Const_Parameters,     only : DP
  implicit none

  integer nchg_type
  integer :: natm2, natm
  real(kind=DP) ::  x_origin,y_origin,z_origin
  real(kind=DP), dimension(3) :: cell1,cell2,cell3
  real(kind=DP), allocatable, dimension(:) :: ival
  real(kind=DP), allocatable, dimension(:,:) :: atom_pos

  integer, dimension(3) :: fft_param
  integer, allocatable, dimension(:) :: iatomtype

contains

  subroutine m_CF_rd_Header_For_Cube(nfcntn_bin)
    integer, intent(in) :: nfcntn_bin
    integer :: i

    read(nfcntn_bin) natm2,x_origin,y_origin,z_origin, natm

    read(nfcntn_bin) fft_param(1), cell1(1:3)
    read(nfcntn_bin) fft_param(2), cell2(1:3)
    read(nfcntn_bin) fft_param(3), cell3(1:3)

    allocate(iatomtype(natm2))
    allocate(ival(natm2))
    allocate(atom_pos(natm2,3))

    do i=1,natm2
      read(nfcntn_bin) iatomtype(i), ival(i), atom_pos(i,1:3)
    enddo

  end subroutine m_CF_rd_Header_For_Cube

end module m_Charge_File
