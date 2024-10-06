!================================================
!  Software name : STM
!  Module : m_Kpoints
!  Subroutine(s) : m_Kp_rd_kv3, m_Kp_alloc, m_Kp_rd
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

module m_Kpoints
! $Id: m_Kpoints.f90,v 1.1.1.1 2004/06/26 11:21:37 yamasaki Exp $
  use m_ArraySize_Parameters, only : knv3
  use m_Const_Parameters, only     : DP, CARTS, BUCS, CRDTYP

  implicit none

  integer                                      :: kv3   ! #sampling k-points
  integer                                      :: nspin=1 ! nspin( 1 or 2 )
  real(kind=DP), allocatable, dimension(:,:,:) :: vkxyz ! sampling k-points

contains
!  subroutine m_Kp_rd_kv3(nfinp)
!    integer, intent(in) :: nfinp
!
!    read(nfinp,*) kv3, nspin
!    print *, ' kv3, nspin =', kv3, nspin
!  end subroutine m_Kp_rd_kv3

  subroutine m_Kp_set_nspin(nsp)
    integer, intent(in) :: nsp
    nspin = nsp
  end subroutine m_Kp_set_nspin

  subroutine m_Kp_alloc
    allocate(vkxyz(kv3,3,CRDTYP))
  end subroutine m_Kp_alloc

  subroutine m_Kp_rd_kpoints(nfinp)
    integer, intent(in) :: nfinp
    integer             :: i, j, k

    do i = 1, kv3 * nspin
      read(nfinp,*) j, (vkxyz(i,k,CARTS),k=1,3),(vkxyz(i, k, BUCS), k=1,3)
    end do
  end subroutine m_Kp_rd_kpoints

end module m_Kpoints
