!================================================
!  Software name : STM
!  Module : m_Crystal_Structure
!  Subroutine(s) : m_CS_rd_neg_and_zl
!  Author(s)     : Takahiro Yamasaki (June 7, 2004)
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

module m_Crystal_Structure
! $Id: m_Crystal_Structure.f90,v 1.1.1.1 2004/06/26 11:21:37 yamasaki Exp $
  implicit none
  integer :: neg
  real    :: zl
contains
  subroutine m_CS_rd_neg_and_zl(nfinp)
    integer, intent(in) :: nfinp
    read(nfinp,*) neg,zl
  end subroutine m_CS_rd_neg_and_zl
end module m_Crystal_Structure
