!================================================
!  Software name : STM
!  Module : m_Const_Parameters
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

module m_Const_Parameters
! $Id: m_Const_Parameters.f90,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $
  implicit none
  integer, parameter :: DP = kind(1.d0)
  integer, parameter :: SP = kind(1.0)
  integer, parameter :: CMPLDP = kind((1.d0, 0.d0))
  integer, parameter :: CMPLSP = kind((1.0 , 0.0 ))
  integer, parameter :: ON = 1, OFF = 0
  integer, parameter :: INVERSE = 1, DIRECT = 2
  integer, parameter :: BEGIN = 1, END = 2
  integer, parameter :: NODATA = 0, PUCV = 0, CARTS = 1, BUCS = 2
  integer, parameter :: CRDTYP = BUCS
  integer, parameter :: OLD = 0, NEW = 1, UNKNOWN = 2
  integer, parameter :: FORMATTED = 0, UNFORMATTED = 1
  integer, parameter :: check_file_name_on  = 1
  integer, parameter :: check_file_name_off = 0
  integer, parameter :: GENERAL = 100
  integer, parameter :: LOWER = 1, IN_BETWEEN = 0, HIGHER = 2
  integer, parameter :: NEGATIVE = -1
  real(kind=DP), parameter :: PAI  = 3.141592653589793238462D0
  real(kind=DP), parameter :: PAI2 = PAI*2.d0, PAI4 = PAI*4.d0
  real(kind=DP), parameter :: Hartree = 27.21139615d0
  real(kind=DP), parameter :: BOHR = 0.5291772480d0

end module m_Const_Parameters
