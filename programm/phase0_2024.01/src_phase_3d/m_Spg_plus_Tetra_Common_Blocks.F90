!-----------------------------------------------------------------------!
!     Module Spg_plus_Tetra/Common_Blocks ver.1.0
!
!     Author : T. Hamada
!   
!     December 20, 2021
!
!     This module is a dummy module containing global variales
!     used by subroutinnes of spg+tetra.F90. The variables are
!     the common block variables used by spg+tetra.F subroutines.
!     This module also contains logical variables and constants
!     used by the spg+tetra.F90 subroutines.
! --------------------------------------------------------------------
      module m_Spg_plus_Tetra_Common_Blocks
! nspg00
         character(len=5), dimension(48) :: schoen
         character(len=2), dimension(3,48) :: jones
! nspg0
         integer :: ng0_0
         integer, dimension(3,24) :: ieuler
         integer, dimension(3,3,48) :: irot
         integer, dimension(48,48) :: im0
         integer, dimension(48) :: iv0
         integer, dimension(3,48) :: ir1234
         real(kind=8), dimension(3,24) :: euler
         real(kind=8), dimension(3,3,48) :: rot

! nspg2
         integer :: ig02, iv02
         integer, dimension(3) :: ieule2, ieulv2
         integer, dimension(3,3) :: irotr2, irotk2
         real(kind=8), dimension(3) :: eeule2, eeulv2

! nspg03
         real(kind=8), dimension(3,3) :: tca, tac, tab, tba, tcb, tbc

! nspg04

         real(kind=8), dimension(3,3) :: grc, gkc, gra, gka, grb, gkb
! nspg06
         integer :: id0, il, ng0, ng1
         integer, dimension(48) :: ig01, iv1
         integer, dimension(48,48) :: im1
         integer, dimension(3,3,48) :: lra1, lsa1, lrb1, lsb1
         real(kind=8), dimension(3,48) :: ta1, tb1
! nspg07
         integer, dimension(48) :: inver1
         real(kind=8), dimension(3,48) :: euler1

! spg1
         integer, dimension(3,48) :: it
         integer, dimension(48,48) :: im
         integer, dimension(48) :: iv

! spg2
         integer :: il9, ng9
         integer, dimension(48) :: ig0
         integer, dimension(2,3,48) :: jv

! ntcntr
         integer :: ncounter, ncounter1, ncounter2, ncounter3

! ipr
         integer, dimension(10,2) :: ipa

! logical
         logical :: nstt3i_called = .false.
         logical :: nstt1i_flag = .false.
         logical :: nskpb0_s_called = .false.
         logical :: nskp00_called = .false.
         logical :: ka00_called = .false.
         logical  :: nskma0_called = .false.
         logical :: tet_sub_called =  .false.
         logical :: ip0_index_generated = .false.
         logical :: tbspg_called = .false.
         logical :: tbspg_kt_called = .false.
         logical :: setkp0_default_n_called = .false.
         logical :: setkp0_default_n_kt_called = .false.

! nsdos0 and nstt1i
         integer :: npx, npy, npz, np, ntet, ncub
         integer, allocatable, dimension(:) :: ip0_index1
         integer, dimension(2,2,2) :: ip0_index2
         real(kind=8) :: rntet_1
         real(kind=8) :: tetra_eps = 1.0d-4

! constant
         real(kind=8) :: factor_1_3 = 1.0d0/3.0d0
         real(kind=8) :: factor_1_12 = 1.0d0/12.0d0
         contains
         subroutine dummy
         end subroutine dummy
      end  module m_Spg_plus_Tetra_Common_Blocks
