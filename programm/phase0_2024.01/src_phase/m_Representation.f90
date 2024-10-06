!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 606 $)
!
!  MODULE: m_Representation
!
!  AUTHOR(S): T. Yamamoto   May/01/2004
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
module m_Representation
  use m_Const_Parameters,  only : DP, CMPLDP, zi
  implicit none

  integer, parameter :: nrep0 = 12
  integer, parameter :: nsym0 = 48

  integer, parameter :: nrep0_dblegrp = 16
  
contains

  subroutine get_irr_rep(point_group,nsym,nrep,rep,symbol_irrep,active)
    character(len=3), intent(in) :: point_group
    integer, intent(out) :: nsym
    integer, intent(out) :: nrep
    integer, intent(out) :: rep(nsym0,nrep0)
    character(len=3), intent(out) :: symbol_irrep(nrep0) 
    character(len=4), intent(out) :: active(nrep0) 

    rep(1:nsym0,1:nrep0) = 0
    symbol_irrep(1:nrep0) = 'NON'
    active(1:nrep0) = 'NON '

    select case(point_group)  
    case('Oh ') !  1 
       call set_chartab_oh(rep)
       call set_symbol_oh(nrep,symbol_irrep,active)
       nsym = 48
    case('O  ') !  2
       call set_chartab_o(rep)
       call set_symbol_o(nrep,symbol_irrep,active)
       nsym = 24
    case('Td ') !  3
       call set_chartab_td(rep)
       call set_symbol_td(nrep,symbol_irrep,active)
       nsym = 24
    case('Th ') !  4
       call set_chartab_th(rep)
       call set_symbol_th(nrep,symbol_irrep,active)
       nsym = 24
    case('T  ') !  5
       call set_chartab_t(rep)
       call set_symbol_t(nrep,symbol_irrep,active)
       nsym = 12
    case('D4h') !  6
       call set_chartab_d4h(rep)
       call set_symbol_d4h(nrep,symbol_irrep,active)
       nsym = 16
    case('D4 ') !  7
       call set_chartab_d4(rep)
       call set_symbol_d4(nrep,symbol_irrep,active)
       nsym = 8
    case('D2d') !  8 
       call set_chartab_d2d(rep)
       call set_symbol_d2d(nrep,symbol_irrep,active)
       nsym = 8
    case('C4v') !  9
       call set_chartab_c4v(rep)
       call set_symbol_c4v(nrep,symbol_irrep,active)
       nsym = 8
    case('C4h') ! 10
       call set_chartab_c4h(rep)
       call set_symbol_c4h(nrep,symbol_irrep,active)
       nsym = 8
    case('S4 ') ! 11
       call set_chartab_s4(rep)
       call set_symbol_s4(nrep,symbol_irrep,active)
       nsym = 4
    case('C4 ') ! 12
       call set_chartab_c4(rep)
       call set_symbol_c4(nrep,symbol_irrep,active)
       nsym = 4
    case('D2h') ! 13
       call set_chartab_d2h(rep)
       call set_symbol_d2h(nrep,symbol_irrep,active)
       nsym = 8
    case('D2 ') ! 14
       call set_chartab_d2(rep)
       call set_symbol_d2(nrep,symbol_irrep,active)
       nsym = 4
    case('C2v') ! 15
       call set_chartab_c2v(rep)
       call set_symbol_c2v(nrep,symbol_irrep,active)
       nsym = 4
    case('D6h') ! 16
       call set_chartab_d6h(rep)
       call set_symbol_d6h(nrep,symbol_irrep,active)
       nsym = 24
    case('D6 ') ! 17
       call set_chartab_d6(rep)
       call set_symbol_d6(nrep,symbol_irrep,active)
       nsym = 12
    case('D3h') ! 18
       call set_chartab_d3h(rep)
       call set_symbol_d3h(nrep,symbol_irrep,active)
       nsym = 12
    case('C6v') ! 19
       call set_chartab_c6v(rep)
       call set_symbol_c6v(nrep,symbol_irrep,active)
       nsym = 12
    case('C6h') ! 20
       call set_chartab_c6h(rep)
       call set_symbol_c6h(nrep,symbol_irrep,active)
       nsym = 12
    case('C3h') ! 21
       call set_chartab_c3h(rep)
       call set_symbol_c3h(nrep,symbol_irrep,active)
       nsym = 6
    case('C6 ') ! 22
       call set_chartab_c6(rep)
       call set_symbol_c6(nrep,symbol_irrep,active)
       nsym = 6
    case('D3d') ! 23
       call set_chartab_d3d(rep)
       call set_symbol_d3d(nrep,symbol_irrep,active)
       nsym = 12
    case('D3 ') ! 24
       call set_chartab_d3(rep)
       call set_symbol_d3(nrep,symbol_irrep,active)
       nsym = 6
    case('C3v') ! 25
       call set_chartab_c3v(rep)
       call set_symbol_c3v(nrep,symbol_irrep,active)
       nsym = 6
    case('S6 ') ! 26
       call set_chartab_s6(rep)
       call set_symbol_s6(nrep,symbol_irrep,active)
       nsym = 6
    case('C3 ') ! 27
       call set_chartab_c3(rep)
       call set_symbol_c3(nrep,symbol_irrep,active)
       nsym = 3
    case('C2h') ! 28
       call set_chartab_c2h(rep)
       call set_symbol_c2h(nrep,symbol_irrep,active)
       nsym = 4
    case('Cs ') ! 29
       call set_chartab_cs(rep)
       call set_symbol_cs(nrep,symbol_irrep,active)
       nsym = 2
    case('C2 ') ! 30
       call set_chartab_c2(rep)
       call set_symbol_c2(nrep,symbol_irrep,active)
       nsym = 2
    case('Ci ') ! 30
       call set_chartab_ci(rep)
       call set_symbol_ci(nrep,symbol_irrep,active)
       nsym = 2
    case('C1 ') ! 32
       call set_chartab_c1(rep)
       call set_symbol_c1(nrep,symbol_irrep,active)
       nsym = 1
    end select

  contains

    subroutine set_chartab_oh(rep)
      integer, intent(out) :: rep(nsym0,nrep0)

      ! local varibales
      integer :: i
      integer, parameter :: nsym = 48
      integer, parameter :: nsym2 = nsym/2
      integer, parameter :: nrep = 10
      integer, parameter :: nrep2 = nrep/2

      call set_chartab_o(rep)
      do i=1,nsym2
         rep(i,nrep2+1:nrep)=rep(i,1:nrep2)
      end do
      do i=1,nsym2
         rep(i+nsym2,1:nrep2)=rep(i,1:nrep2)
      end do
      do i=nsym2+1,nsym
         rep(i,nrep2+1:nrep)=-rep(i,1:nrep2)
      end do

    end subroutine set_chartab_oh

    subroutine set_symbol_oh(nrep,symbol_irrep,active)
      integer, intent(out) :: nrep
      character(len=3), intent(out) :: symbol_irrep(10)
      character(len=4), intent(out) :: active(10)

      nrep = 10
      symbol_irrep(1) = 'A1g'; active(1) = 'R   '
      symbol_irrep(2) = 'A2g' 
      symbol_irrep(3) = 'Eg '; active(3) = 'R   '
      symbol_irrep(4) = 'T1g'
      symbol_irrep(5) = 'T2g'; active(5) = 'R   '
      symbol_irrep(6) = 'A1u'
      symbol_irrep(7) = 'A2u'
      symbol_irrep(8) = 'Eu '
      symbol_irrep(9) = 'T1u'; active(9) = 'IR  '
      symbol_irrep(10) = 'T2u'
    end subroutine set_symbol_oh

    subroutine set_chartab_o(rep)
      integer, intent(out) :: rep(nsym0,nrep0)

      ! local varibales
      integer, parameter :: nsym = 24

      rep(1:nsym,1)=1 ! A1

      rep(1:12,2)=1; rep(13:nsym,2)=-1 ! A2

      rep(1:4,3)=2; rep(5:12,3)=-1; rep(13:nsym,3)=0 ! E

      rep(1,4)=3; rep(2:4,4)=-1; rep(5:12,4)=0 ! T1
      rep(13:18,4)=-1; rep(19:nsym,4)=1          ! T1 (cont.)

      rep(1,5)=3; rep(2:4,5)=-1; rep(5:12,5)=0; rep(13:18,5)=1; rep(19:nsym,5)=-1 ! T2

    end subroutine set_chartab_o

    subroutine set_symbol_o(nrep,symbol_irrep,active)
      integer, intent(out) :: nrep
      character(len=3), intent(out) :: symbol_irrep(5)
      character(len=4), intent(out) :: active(5)
      nrep = 5
      symbol_irrep(1) = 'A1 '; active(1)='R   '
      symbol_irrep(2) = 'A2 '
      symbol_irrep(3) = 'E  '; active(3)='R   '
      symbol_irrep(4) = 'T1 '; active(4)='IR  '
      symbol_irrep(5) = 'T2 '; active(5)='R   '
    end subroutine set_symbol_o

    subroutine set_chartab_td(rep)
      integer, intent(out) :: rep(nsym0,nrep0)
      call set_chartab_o(rep)
    end subroutine set_chartab_td

    subroutine set_symbol_td(nrep,symbol_irrep,active)
      integer, intent(out) :: nrep
      character(len=3), intent(out) :: symbol_irrep(5)
      character(len=4), intent(out) :: active(5)
      call set_symbol_o(nrep,symbol_irrep,active)
      symbol_irrep(5) = 'T2 '; active(5)='IR&R'
    end subroutine set_symbol_td

    subroutine set_chartab_th(rep)
      integer, intent(out) :: rep(nsym0,nrep0)

      ! local varibales
      integer :: i
      integer, parameter :: nsym = 24
      integer, parameter :: nsym2 = nsym/2
      integer, parameter :: nrep = 6
      integer, parameter :: nrep2 = nrep/2

      call set_chartab_t(rep)
      do i=1,nsym2
         rep(i,nrep2+1:nrep)=rep(i,1:nrep2)
      end do
      do i=1,nsym2
         rep(i+nsym2,1:nrep2)=rep(i,1:nrep2)
      end do
      do i=nsym2+1,nsym
         rep(i,nrep2+1:nrep)=-rep(i,1:nrep2)
      end do

    end subroutine set_chartab_th

    subroutine set_symbol_th(nrep,symbol_irrep,active)
      integer, intent(out) :: nrep
      character(len=3), intent(out) :: symbol_irrep(6)
      character(len=4), intent(out) :: active(6)

      nrep = 6
      symbol_irrep(1) = 'Ag '; active(1)='R   '
      symbol_irrep(2) = 'Eg '; active(2)='R   '
      symbol_irrep(3) = 'Tg '; active(3)='R   '
      symbol_irrep(4) = 'Au '
      symbol_irrep(5) = 'Eu '
      symbol_irrep(6) = 'Tu '; active(6)='IR  '

    end subroutine set_symbol_th

    subroutine set_chartab_t(rep)
      integer, intent(out) :: rep(nsym0,nrep0)

      ! local varibales
      integer, parameter :: nsym = 12

      rep(1:nsym,1)=1 ! A

      rep(1:4,2)=2; rep(5:nsym,2)=-1 ! E

      rep(1,3)=3; rep(2:4,3)=-1; rep(5:nsym,3)=0 ! T

    end subroutine set_chartab_t

    subroutine set_symbol_t(nrep,symbol_irrep,active)
      integer, intent(out) :: nrep
      character(len=3), intent(out) :: symbol_irrep(3)
      character(len=4), intent(out) :: active(3)

      nrep = 3
      symbol_irrep(1) = 'A  '; active(1) = 'R   '
      symbol_irrep(2) = 'E  '; active(2) = 'R   '
      symbol_irrep(3) = 'T  '; active(3) = 'IR&R'

    end subroutine set_symbol_t

    subroutine set_chartab_d4h(rep)
      integer, intent(out) :: rep(nsym0,nrep0)

      ! local varibales
      integer :: i
      integer, parameter :: nsym = 16
      integer, parameter :: nsym2 = nsym/2
      integer, parameter :: nrep = 10
      integer, parameter :: nrep2 = nrep/2

      call set_chartab_d4(rep)
      do i=1,nsym2
         rep(i,nrep2+1:nrep)=rep(i,1:nrep2)
      end do
      do i=1,nsym2
         rep(i+nsym2,1:nrep2)=rep(i,1:nrep2)
      end do
      do i=nsym2+1,nsym
         rep(i,nrep2+1:nrep)=-rep(i,1:nrep2)
      end do

    end subroutine set_chartab_d4h

    subroutine set_symbol_d4h(nrep,symbol_irrep,active)
      integer, intent(out) :: nrep
      character(len=3), intent(out) :: symbol_irrep(10)
      character(len=4), intent(out) :: active(10)

      nrep = 10
      symbol_irrep(1) = 'A1g'; active(1) = 'R   '
      symbol_irrep(2) = 'A2g'
      symbol_irrep(3) = 'B1g'; active(3) = 'R   '
      symbol_irrep(4) = 'B2g'; active(4) = 'R   '
      symbol_irrep(5) = 'Eg '; active(5) = 'R   '
      symbol_irrep(6) = 'A1u'
      symbol_irrep(7) = 'A2u'; active(7) = 'IR  '
      symbol_irrep(8) = 'B1u'
      symbol_irrep(9) = 'B2u'
      symbol_irrep(10) = 'Eu '; active(10) = 'IR'
    end subroutine set_symbol_d4h

    subroutine set_chartab_d4(rep)
      integer, intent(out) :: rep(nsym0,nrep0)

      ! local varibales
      integer, parameter :: nsym = 8

      rep(1:nsym,1)=1 ! A1

      rep(1,2)=1; rep(2:3,2)=-1; rep(4,2)=1; rep(5:6,2)=-1; rep(7:nsym,2)=1 ! A2

      rep(1:4,3)=1; rep(5:nsym,3)=-1 ! B1

      rep(1,4)=1; rep(2:3,4)=-1; rep(4:6,4)=1; rep(7:nsym,4)=-1 ! B2

      rep(1,5)=2; rep(2:3,5)=0; rep(4,5)=-2; rep(6:nsym,5)=0 ! E

    end subroutine set_chartab_d4

    subroutine set_symbol_d4(nrep,symbol_irrep,active)
      integer, intent(out) :: nrep
      character(len=3), intent(out) :: symbol_irrep(5)
      character(len=4), intent(out) :: active(5)

      nrep = 5
      symbol_irrep(1) = 'A1 '; active(1) = 'R   '
      symbol_irrep(2) = 'A2 '; active(2) = 'IR  '
      symbol_irrep(3) = 'B1 '; active(3) = 'R   '
      symbol_irrep(4) = 'B2 '; active(4) = 'R   '
      symbol_irrep(5) = 'E  '; active(5) = 'IR&R'

    end subroutine set_symbol_d4

    subroutine set_chartab_d2d(rep)
      integer, intent(out) :: rep(nsym0,nrep0)

      ! local varibales
      integer, parameter :: nsym = 8

      rep(1:nsym,1)=1 ! A1

      rep(1:2,2)=1; rep(3:6,2)=-1; rep(7:nsym,2)=1 ! A2

      rep(1:4,3)=1; rep(5:nsym,3)=-1 ! B1

      rep(1:2,4)=1; rep(3:4,4)=-1; rep(5:6,4)=1; rep(7:nsym,4)=-1 ! B2

      rep(1,5)=2; rep(2,5)=-2; rep(3:nsym,5)=0 ! E

    end subroutine set_chartab_d2d

    subroutine set_symbol_d2d(nrep,symbol_irrep,active)
      integer, intent(out) :: nrep
      character(len=3), intent(out) :: symbol_irrep(5)
      character(len=4), intent(out) :: active(5)

      nrep = 5
      symbol_irrep(1) = 'A1 '; active(1) = 'R   '
      symbol_irrep(2) = 'A2 '
      symbol_irrep(3) = 'B1 '; active(3) = 'R   '
      symbol_irrep(4) = 'B2 '; active(4) = 'IR&R'
      symbol_irrep(5) = 'E  '; active(5) = 'IR&R'

    end subroutine set_symbol_d2d

    subroutine set_chartab_c4v(rep)
      integer, intent(out) :: rep(nsym0,nrep0)

      ! local varibales
      integer, parameter :: nsym = 8

      rep(1:nsym,1)=1 ! A1

      rep(1:4,2)=1; rep(5:nsym,2)=-1 ! A2

      rep(1:2,3)=1; rep(3:4,3)=-1; rep(5:6,3)=1; rep(7:nsym,3)=-1 ! B1

      rep(1:2,4)=1; rep(3:4,4)=-1; rep(5:6,4)=-1; rep(7:nsym,4)=1 ! B2

      rep(1,5)=2; rep(2,5)=-2; rep(3:nsym,5)=0 ! E

    end subroutine set_chartab_c4v

    subroutine set_symbol_c4v(nrep,symbol_irrep,active)
      integer, intent(out) :: nrep
      character(len=3), intent(out) :: symbol_irrep(5)
      character(len=4), intent(out) :: active(5)

      nrep = 5
      symbol_irrep(1) = 'A1 '; active(1) = 'IR&R'
      symbol_irrep(2) = 'A2 '
      symbol_irrep(3) = 'B1 '; active(3) = 'R   '
      symbol_irrep(4) = 'B2 '; active(4) = 'R   '
      symbol_irrep(5) = 'E  '; active(5) = 'IR&R'

    end subroutine set_symbol_c4v

    subroutine set_chartab_c4h(rep)
      integer, intent(out) :: rep(nsym0,nrep0)

      integer :: i
      integer, parameter :: nsym = 8
      integer, parameter :: nsym2 = nsym/2
      integer, parameter :: nrep = 6
      integer, parameter :: nrep2 = nrep/2

      call set_chartab_c4(rep)
      do i=1,nsym2
         rep(i,nrep2+1:nrep)=rep(i,1:nrep2)
      end do
      do i=1,nsym2
         rep(i+nsym2,1:nrep2)=rep(i,1:nrep2)
      end do
      do i=nsym2+1,nsym
         rep(i,nrep2+1:nrep)=-rep(i,1:nrep2)
      end do

    end subroutine set_chartab_c4h

    subroutine set_symbol_c4h(nrep,symbol_irrep,active)
      integer, intent(out) :: nrep
      character(len=3), intent(out) :: symbol_irrep(6)
      character(len=4), intent(out) :: active(6)

      nrep = 6
      symbol_irrep(1) = 'Ag '; active(1) = 'R   '
      symbol_irrep(2) = 'Bg '; active(2) = 'R   '
      symbol_irrep(3) = 'Eg '; active(3) = 'R   '
      symbol_irrep(4) = 'Au '; active(4) = 'IR  '
      symbol_irrep(5) = 'Bu '
      symbol_irrep(6) = 'Eu '; active(6) = 'IR  '

    end subroutine set_symbol_c4h

    subroutine set_chartab_s4(rep)
      integer, intent(out) :: rep(nsym0,nrep0)

      ! local varibales
      integer, parameter :: nsym = 4

      rep(1:nsym,1)=1 ! A

      rep(1:2,2)=1; rep(3:nsym,2)=-1 ! B

      rep(1,3)=2; rep(2,3)=-2; rep(3:nsym,3)=0 ! E

    end subroutine set_chartab_s4

    subroutine set_symbol_s4(nrep,symbol_irrep,active)
      integer, intent(out) :: nrep
      character(len=3), intent(out) :: symbol_irrep(3)
      character(len=4), intent(out) :: active(3)

      nrep = 3
      symbol_irrep(1) = 'A  '; active(1) = 'R   '
      symbol_irrep(2) = 'B  '; active(2) = 'IR&R'
      symbol_irrep(3) = 'E  '; active(3) = 'IR&R'

    end subroutine set_symbol_s4

    subroutine set_chartab_c4(rep)
      integer, intent(out) :: rep(nsym0,nrep0)

      ! local varibales
      integer, parameter :: nsym = 4

      rep(1:nsym,1)=1 ! A

      rep(1:2,2)=1; rep(3:nsym,2)=-1 ! B

      rep(1,3)=2; rep(2,3)=-2; rep(3:nsym,3)=0 ! E

    end subroutine set_chartab_c4

    subroutine set_symbol_c4(nrep,symbol_irrep,active)
      integer, intent(out) :: nrep
      character(len=3), intent(out) :: symbol_irrep(3)
      character(len=4), intent(out) :: active(3)

      nrep = 3
      call set_symbol_s4(nrep,symbol_irrep,active)
      active(1) = 'IR&R'
      active(2) = 'R   '

    end subroutine set_symbol_c4

    subroutine set_chartab_d2h(rep)
      integer, intent(out) :: rep(nsym0,nrep0)

      integer :: i
      integer, parameter :: nsym = 8
      integer, parameter :: nsym2 = nsym/2
      integer, parameter :: nrep = 8
      integer, parameter :: nrep2 = nrep/2

      call set_chartab_d2(rep)
      do i=1,nsym2
         rep(i,nrep2+1:nrep)=rep(i,1:nrep2)
      end do
      do i=1,nsym2
         rep(i+nsym2,1:nrep2)=rep(i,1:nrep2)
      end do
      do i=nsym2+1,nsym
         rep(i,nrep2+1:nrep)=-rep(i,1:nrep2)
      end do

    end subroutine set_chartab_d2h

    subroutine set_symbol_d2h(nrep,symbol_irrep,active)
      integer, intent(out) :: nrep
      character(len=3), intent(out) :: symbol_irrep(8)
      character(len=4), intent(out) :: active(8)

      nrep = 8
      symbol_irrep(1) = 'Ag '; active(1) = 'R  '
      symbol_irrep(2) = 'B1g'; active(2) = 'R  '
      symbol_irrep(3) = 'B2g'; active(3) = 'R  '
      symbol_irrep(4) = 'B3g'; active(4) = 'R  '
      symbol_irrep(5) = 'Au '
      symbol_irrep(6) = 'B1u'; active(6) = 'IR '
      symbol_irrep(7) = 'B2u'; active(7) = 'IR '
      symbol_irrep(8) = 'B3u'; active(8) = 'IR '

    end subroutine set_symbol_d2h

    subroutine set_chartab_d2(rep)
      integer, intent(out) :: rep(nsym0,nrep0)

      ! local varibales
      integer, parameter :: nsym = 4

      rep(1:nsym,1)=1 ! A

      rep(1,2)=1; rep(2:3,2)=-1; rep(nsym,2)=1 ! B1

      rep(1,3)=1; rep(2,3)=-1; rep(3,3)=1; rep(4,3)=-1 ! B2

      rep(1:2,4)=1; rep(3:nsym,4)=-1 ! B3

    end subroutine set_chartab_d2

    subroutine set_symbol_d2(nrep,symbol_irrep,active)
      integer, intent(out) :: nrep
      character(len=3), intent(out) :: symbol_irrep(4)
      character(len=4), intent(out) :: active(4)

      nrep = 4
      symbol_irrep(1) = 'A  '; active(1) = 'R   '
      symbol_irrep(2) = 'B1 '; active(2) = 'IR&R'
      symbol_irrep(3) = 'B2 '; active(3) = 'IR&R'
      symbol_irrep(4) = 'B3 '; active(4) = 'IR&R'

    end subroutine set_symbol_d2

    subroutine set_chartab_c2v(rep)
      integer, intent(out) :: rep(nsym0,nrep0)

      ! local varibales
      integer, parameter :: nsym = 4

      rep(1:nsym,1)=1 ! A1

      rep(1:2,2)=1; rep(3:nsym,2)=-1 ! A2

      rep(1,3)=1; rep(2:3,3)=-1; rep(nsym,3)=1 ! B1

      rep(1,4)=1; rep(2,4)=-1; rep(3,4)=1; rep(nsym,4)=-1 ! B2

    end subroutine set_chartab_c2v

    subroutine set_symbol_c2v(nrep,symbol_irrep,active)
      integer, intent(out) :: nrep
      character(len=3), intent(out) :: symbol_irrep(4)
      character(len=4), intent(out) :: active(4)

      nrep = 4
      symbol_irrep(1) = 'A1 '; active(1) = 'IR&R'
      symbol_irrep(2) = 'A2 '; active(2) = 'R   '
      symbol_irrep(3) = 'B1 '; active(3) = 'IR&R'
      symbol_irrep(4) = 'B2 '; active(4) = 'IR&R'

    end subroutine set_symbol_c2v

    subroutine set_chartab_d6h(rep)
      integer, intent(out) :: rep(nsym0,nrep0)

      integer :: i
      integer, parameter :: nsym = 24
      integer, parameter :: nsym2 = nsym/2
      integer, parameter :: nrep = 12
      integer, parameter :: nrep2 = nrep/2

      call set_chartab_d6(rep)
      do i=1,nsym2
         rep(i,nrep2+1:nrep)=rep(i,1:nrep2)
      end do
      do i=1,nsym2
         rep(i+nsym2,1:nrep2)=rep(i,1:nrep2)
      end do
      do i=nsym2+1,nsym
         rep(i,nrep2+1:nrep)=-rep(i,1:nrep2)
      end do

    end subroutine set_chartab_d6h

    subroutine set_symbol_d6h(nrep,symbol_irrep,active)
      integer, intent(out) :: nrep
      character(len=3), intent(out) :: symbol_irrep(12)
      character(len=4), intent(out) :: active(12)

      nrep = 12
      symbol_irrep(1) = 'A1g'; active(1) = 'R   '
      symbol_irrep(2) = 'A2g'
      symbol_irrep(3) = 'B1g'
      symbol_irrep(4) = 'B2g'
      symbol_irrep(5) = 'E1g'; active(5) = 'R   '
      symbol_irrep(6) = 'E2g'; active(6) = 'R   '
      symbol_irrep(7) = 'A1u'
      symbol_irrep(8) = 'A2u'; active(8) = 'IR  '
      symbol_irrep(9) = 'B1u'
      symbol_irrep(10) = 'B2u'
      symbol_irrep(11) = 'E1u'; active(11) = 'IR  '
      symbol_irrep(12) = 'E2u'

    end subroutine set_symbol_d6h

    subroutine set_chartab_d6(rep)
      integer, intent(out) :: rep(nsym0,nrep0)

      ! local varibales
      integer, parameter :: nsym = 12

      rep(1:nsym,1)=1 ! A1

      rep(1:6,2)=1; rep(7:nsym,2)=-1 ! A2

      rep(1,3)=1; rep(2,3)=-1; rep(3,3)=1; rep(4,3)=-1 ! B1
      rep(5,3)=1; rep(6,3)=-1; rep(7:9,3)=-1; rep(10:nsym,3)=1 ! B1 (cont.)

      rep(1,4)=1; rep(2,4)=-1; rep(3,4)=1; rep(4,4)=-1 ! B2
      rep(5,4)=1; rep(6,4)=-1; rep(7:9,4)=1; rep(10:nsym,4)=-1 ! B2 (cont.)

      rep(1,5)=2; rep(2,5)=1; rep(3,5)=-1; rep(4,5)=-2 ! E1
      rep(5,5)=-1; rep(6,5)=1; rep(7:nsym,5)=0 ! E1 (cont.)

      rep(1,6)=2; rep(2:3,6)=-1; rep(4,6)=2 ! E2
      rep(5:6,6)=-1; rep(7:nsym,6)=0 ! E2 (cont.)

    end subroutine set_chartab_d6

    subroutine set_symbol_d6(nrep,symbol_irrep,active)
      integer, intent(out) :: nrep
      character(len=3), intent(out) :: symbol_irrep(6)
      character(len=4), intent(out) :: active(6)

      nrep = 6
      symbol_irrep(1) = 'A1 '; active(1) = 'R   '
      symbol_irrep(2) = 'A2 '; active(2) = 'IR  '
      symbol_irrep(3) = 'B1 '
      symbol_irrep(4) = 'B2 '
      symbol_irrep(5) = 'E1 '; active(3) = 'IR&R'
      symbol_irrep(6) = 'E2 '; active(4) = 'R   '

    end subroutine set_symbol_d6

    subroutine set_chartab_d3h(rep)
      integer, intent(out) :: rep(nsym0,nrep0)

      ! local varibales
      integer, parameter :: nsym = 12

      rep(1:nsym,1)=1 ! A1

      rep(1:3,2)=1; rep(4:6,2)=-1; rep(7:9,2)=1; rep(10:nsym,2)=-1 ! A2

      rep(1:6,3)=1; rep(7:nsym,3)=-1 ! B1

      rep(1:3,4)=1; rep(4:9,4)=-1; rep(10:nsym,4)=1 ! B2

      rep(1,5)=2; rep(2:3,5)=-1; rep(4:6,5)=0; rep(7,5)=1 ! E1
      rep(8,5)=-2; rep(9,5)=1; rep(10:nsym,5)=0 ! E1 (cont.)

      rep(1,6)=2; rep(2:3,6)=-1; rep(4:6,6)=0; rep(7,6)=-1 ! E2
      rep(8,6)=2; rep(9,6)=-1; rep(10:nsym,6)=0 ! E2 (cont.)

    end subroutine set_chartab_d3h

    subroutine set_symbol_d3h(nrep,symbol_irrep,active)
      integer, intent(out) :: nrep
      character(len=3), intent(out) :: symbol_irrep(6)
      character(len=4), intent(out) :: active(6)

      nrep = 6
      symbol_irrep(1) = "A1'"; active(1) = 'R   '
      symbol_irrep(2) = "A2'"
      symbol_irrep(3) = 'A1"'
      symbol_irrep(4) = 'A2"'; active(4) = 'IR  '
      symbol_irrep(5) = 'E" '; active(5) = 'R   '         !! ASMS 2015/03/17
      symbol_irrep(6) = "E' "; active(6) = 'IR&R'         !! ASMS 2015/03/17

    end subroutine set_symbol_d3h

    subroutine set_chartab_c6v(rep)
      integer, intent(out) :: rep(nsym0,nrep0)

      ! local varibales
      integer, parameter :: nsym = 12

      call set_chartab_d6(rep)

    end subroutine set_chartab_c6v

    subroutine set_symbol_c6v(nrep,symbol_irrep,active)
      integer, intent(out) :: nrep
      character(len=3), intent(out) :: symbol_irrep(6)
      character(len=4), intent(out) :: active(6)

      nrep = 6
      symbol_irrep(1) = 'A1 '; active(1) = 'IR&R'
      symbol_irrep(2) = 'A2 '
      symbol_irrep(3) = 'B1 '
      symbol_irrep(4) = 'B2 '
      symbol_irrep(5) = 'E1 '; active(5) = 'IR&R'
      symbol_irrep(6) = 'E2 '; active(6) = 'R   '

    end subroutine set_symbol_c6v

    subroutine set_chartab_c6h(rep)
      integer, intent(out) :: rep(nsym0,nrep0)

      integer :: i
      integer, parameter :: nsym = 12
      integer, parameter :: nsym2 = nsym/2
      integer, parameter :: nrep = 8
      integer, parameter :: nrep2 = nrep/2

      call set_chartab_c6(rep)
      do i=1,nsym2
         rep(i,nrep2+1:nrep)=rep(i,1:nrep2)
      end do
      do i=1,nsym2
         rep(i+nsym2,1:nrep2)=rep(i,1:nrep2)
      end do
      do i=nsym2+1,nsym
         rep(i,nrep2+1:nrep)=-rep(i,1:nrep2)
      end do

    end subroutine set_chartab_c6h

    subroutine set_symbol_c6h(nrep,symbol_irrep,active)
      integer, intent(out) :: nrep
      character(len=3), intent(out) :: symbol_irrep(8)
      character(len=4), intent(out) :: active(8)

      nrep = 8
      symbol_irrep(1) = 'Ag '; active(1) = 'R   '
      symbol_irrep(2) = 'Bg '
      symbol_irrep(3) = 'E2g'; active(3) = 'R   '
      symbol_irrep(4) = 'E1g'; active(4) = 'R   '
      symbol_irrep(5) = 'Au '; active(5) = 'IR  '
      symbol_irrep(6) = 'Bu '
      symbol_irrep(7) = 'E2u'
      symbol_irrep(8) = 'E1u'; active(8) = 'IR  '

    end subroutine set_symbol_c6h

    subroutine set_chartab_c3h(rep)
      integer, intent(out) :: rep(nsym0,nrep0)

      ! local varibales
      integer, parameter :: nsym = 6

      rep(1:nsym,1)=1 ! A

      rep(1:3,2)=1; rep(4:nsym,2)=-1 ! B

      rep(1,3)=2; rep(2:3,3)=-1; rep(4,3)=-1; ! E2
      rep(5,3)=2; rep(nsym,3)=-1  ! E2

      rep(1,4)=2; rep(2:3,4)=-1; rep(4,4)=1 ! E1
      rep(5,4)=-2; rep(nsym,4)=1 ! E1

    end subroutine set_chartab_c3h

    subroutine set_symbol_c3h(nrep,symbol_irrep,active)
      integer, intent(out) :: nrep
      character(len=3), intent(out) :: symbol_irrep(4)
      character(len=4), intent(out) :: active(4)

      nrep = 4
      symbol_irrep(1) = "A' "; active(1) = 'R   '
      symbol_irrep(2) = 'A" '; active(2) = 'IR  '
      symbol_irrep(3) = "E' "; active(3) = 'IR&R'
      symbol_irrep(4) = 'E" '; active(4) = 'IR&R'

    end subroutine set_symbol_c3h

    subroutine set_chartab_c6(rep)
      integer, intent(out) :: rep(nsym0,nrep0)

      ! local varibales
      integer, parameter :: nsym = 6

      rep(1:nsym,1)=1 ! A

      rep(1,2)=1; rep(2,2)=-1; rep(3,2)=1; rep(4,2)=-1 ! B
      rep(5,2)=1; rep(nsym,2)=-1 ! B

      rep(1,3)=2; rep(2:3,3)=-1; rep(4,3)=2; rep(5:nsym,3)=-1 ! E2

      rep(1,4)=2; rep(2,4)=1; rep(3,4)=-1; rep(4,4)=-2 ! E1
      rep(5,4)=-1; rep(nsym,4)=1 ! E1

    end subroutine set_chartab_c6

    subroutine set_symbol_c6(nrep,symbol_irrep,active)
      integer, intent(out) :: nrep
      character(len=3), intent(out) :: symbol_irrep(4)
      character(len=4), intent(out) :: active(4)

      nrep = 4
      symbol_irrep(1) = 'A  '; active(1) = 'IR&R'
      symbol_irrep(2) = 'B  '
      symbol_irrep(3) = 'E2 '; active(3) = 'R   '
      symbol_irrep(4) = 'E1 '; active(4) = 'IR&R'

    end subroutine set_symbol_c6

    subroutine set_chartab_d3d(rep)
      integer, intent(out) :: rep(nsym0,nrep0)

      integer :: i
      integer, parameter :: nsym = 12 
      integer, parameter :: nsym2 = nsym/2
      integer, parameter :: nrep = 6
      integer, parameter :: nrep2 = nrep/2

      call set_chartab_d3(rep)
      do i=1,nsym2
         rep(i,nrep2+1:nrep)=rep(i,1:nrep2)
      end do
      do i=1,nsym2
         rep(i+nsym2,1:nrep2)=rep(i,1:nrep2)
      end do
      do i=nsym2+1,nsym
         rep(i,nrep2+1:nrep)=-rep(i,1:nrep2)
      end do

    end subroutine set_chartab_d3d

    subroutine set_symbol_d3d(nrep,symbol_irrep,active)
      integer, intent(out) :: nrep
      character(len=3), intent(out) :: symbol_irrep(6)
      character(len=4), intent(out) :: active(6)

      nrep = 6
      symbol_irrep(1) = 'A1g'; active(1) = 'R   '
      symbol_irrep(2) = 'A2g'
      symbol_irrep(3) = 'Eg '; active(3) = 'R   '
      symbol_irrep(4) = 'A1u'
      symbol_irrep(5) = 'A2u'; active(5) = 'IR  '
      symbol_irrep(6) = 'Eu '; active(6) = 'IR  '

    end subroutine set_symbol_d3d

    subroutine set_chartab_d3(rep)
      integer, intent(out) :: rep(nsym0,nrep0)

      ! local varibales
      integer, parameter :: nsym = 6

      rep(1:nsym,1)=1 ! A1

      rep(1:3,2)=1; rep(4:nsym,2)=-1 ! A2

      rep(1,3)=2; rep(2:3,3)=-1; rep(4:nsym,3)=0 ! E

    end subroutine set_chartab_d3

    subroutine set_symbol_d3(nrep,symbol_irrep,active)
      integer, intent(out) :: nrep
      character(len=3), intent(out) :: symbol_irrep(3)
      character(len=4), intent(out) :: active(3)

      nrep = 3
      symbol_irrep(1) = 'A1 '; active(1) = 'R   '
      symbol_irrep(2) = 'A2 '; active(2) = 'IR  '
      symbol_irrep(3) = 'E  '; active(3) = 'IR&R'

    end subroutine set_symbol_d3

    subroutine set_chartab_c3v(rep)
      integer, intent(out) :: rep(nsym0,nrep0)

      ! local varibales
      integer, parameter :: nsym = 6

      rep(1:nsym,1)=1 ! A1

      rep(1:3,2)=1; rep(4:nsym,2)=-1 ! A2

      rep(1,3)=2; rep(2:3,3)=-1; rep(4:nsym,3)=0 ! E

    end subroutine set_chartab_c3v

    subroutine set_symbol_c3v(nrep,symbol_irrep,active)
      integer, intent(out) :: nrep
      character(len=3), intent(out) :: symbol_irrep(3)
      character(len=4), intent(out) :: active(3)

      nrep = 3
      symbol_irrep(1) = 'A1 '; active(1) = 'IR&R'
      symbol_irrep(2) = 'A2 '
      symbol_irrep(3) = 'E  '; active(3) = 'IR&R'

    end subroutine set_symbol_c3v

    subroutine set_chartab_s6(rep)
      integer, intent(out) :: rep(nsym0,nrep0)

      integer :: i
      integer, parameter :: nsym = 6
      integer, parameter :: nsym2 = nsym/2
      integer, parameter :: nrep = 4
      integer, parameter :: nrep2 = nrep/2

      call set_chartab_c3(rep)
      do i=1,nsym2
         rep(i,nrep2+1:nrep)=rep(i,1:nrep2)
      end do
      do i=1,nsym2
         rep(i+nsym2,1:nrep2)=rep(i,1:nrep2)
      end do
      do i=nsym2+1,nsym
         rep(i,nrep2+1:nrep)=-rep(i,1:nrep2)
      end do

    end subroutine set_chartab_s6

    subroutine set_symbol_s6(nrep,symbol_irrep,active)
      integer, intent(out) :: nrep
      character(len=3), intent(out) :: symbol_irrep(4)
      character(len=4), intent(out) :: active(4)

      nrep = 4
      symbol_irrep(1) = 'Ag '; active(1) = 'R   '
      symbol_irrep(2) = 'Eg '; active(2) = 'R   '
      symbol_irrep(3) = 'Au '; active(3) = 'IR  '
      symbol_irrep(4) = 'Eu '; active(4) = 'IR  '

    end subroutine set_symbol_s6

    subroutine set_chartab_c3(rep)
      integer, intent(out) :: rep(nsym0,nrep0)

      ! local varibales
      integer, parameter :: nsym = 3

      rep(1:nsym,1)=1 ! A1

      rep(1,2)=2; rep(2:nsym,2)=-1 ! E

    end subroutine set_chartab_c3

    subroutine set_symbol_c3(nrep,symbol_irrep,active)
      integer, intent(out) :: nrep
      character(len=3), intent(out) :: symbol_irrep(2)
      character(len=4), intent(out) :: active(2)

      nrep = 2
      symbol_irrep(1) = 'A  '; active(1) = 'IR&R'
      symbol_irrep(2) = 'E  '; active(2) = 'IR&R'

    end subroutine set_symbol_c3

    subroutine set_chartab_c2h(rep)
      integer, intent(out) :: rep(nsym0,nrep0)

      integer :: i
      integer, parameter :: nsym = 4
      integer, parameter :: nsym2 = nsym/2
      integer, parameter :: nrep = 4
      integer, parameter :: nrep2 = nrep/2

      call set_chartab_c2(rep)
      do i=1,nsym2
         rep(i,nrep2+1:nrep)=rep(i,1:nrep2)
      end do
      do i=1,nsym2
         rep(i+nsym2,1:nrep2)=rep(i,1:nrep2)
      end do
      do i=nsym2+1,nsym
         rep(i,nrep2+1:nrep)=-rep(i,1:nrep2)
      end do

    end subroutine set_chartab_c2h

    subroutine set_symbol_c2h(nrep,symbol_irrep,active)
      integer, intent(out) :: nrep
      character(len=3), intent(out) :: symbol_irrep(4)
      character(len=4), intent(out) :: active(4)

      nrep = 4
      symbol_irrep(1) = 'Ag '; active(1) = 'R   '
      symbol_irrep(2) = 'Bg '; active(2) = 'R   '
      symbol_irrep(3) = 'Au '; active(3) = 'IR  '
      symbol_irrep(4) = 'Bu '; active(4) = 'IR  '

    end subroutine set_symbol_c2h

    subroutine set_chartab_cs(rep)
      integer, intent(out) :: rep(nsym0,nrep0)

      call set_chartab_c2(rep)

    end subroutine set_chartab_cs

    subroutine set_symbol_cs(nrep,symbol_irrep,active)
      integer, intent(out) :: nrep
      character(len=3), intent(out) :: symbol_irrep(2)
      character(len=4), intent(out) :: active(2)

      nrep = 2
      symbol_irrep(1) = "A' "; active(1) = 'IR&R'
      symbol_irrep(2) = 'A" '; active(2) = 'IR&R'

    end subroutine set_symbol_cs

    subroutine set_chartab_c2(rep)
      integer, intent(out) :: rep(nsym0,nrep0)

      ! local varibales
      integer, parameter :: nsym = 2

      rep(1:nsym,1)=1 ! A

      rep(1,2)=1; rep(nsym,2)=-1 ! B

    end subroutine set_chartab_c2

    subroutine set_symbol_c2(nrep,symbol_irrep,active)
      integer, intent(out) :: nrep
      character(len=3), intent(out) :: symbol_irrep(2)
      character(len=4), intent(out) :: active(2)

      nrep = 2
      symbol_irrep(1) = 'A  '; active(1) = 'IR&R'
      symbol_irrep(2) = 'B  '; active(2) = 'IR&R'

    end subroutine set_symbol_c2

    subroutine set_chartab_ci(rep)
      integer, intent(out) :: rep(nsym0,nrep0)

      ! local varibales
      integer, parameter :: nsym = 2

      rep(1:nsym,1)=1 ! Ag 
      rep(1,2)=1; rep(nsym,2)=-1 ! Au

    end subroutine set_chartab_ci

    subroutine set_symbol_ci(nrep,symbol_irrep,active)
      integer, intent(out) :: nrep
      character(len=3), intent(out) :: symbol_irrep(2)
      character(len=4), intent(out) :: active(2)

      nrep = 2
      symbol_irrep(1) = 'A  '; active(1) = 'R   '
      symbol_irrep(2) = 'B  '; active(2) = 'IR  '

    end subroutine set_symbol_ci

    subroutine set_chartab_c1(rep)
      integer, intent(out) :: rep(nsym0,nrep0)

      ! local varibales
      integer, parameter :: nsym = 1

      rep(1:nsym,1)=1 ! A

    end subroutine set_chartab_c1

    subroutine set_symbol_c1(nrep,symbol_irrep,active)
      integer, intent(out) :: nrep
      character(len=3), intent(out) :: symbol_irrep(1)
      character(len=4), intent(out) :: active(1)

      nrep = 1
      symbol_irrep(1) = 'A  '; active(1) = 'IR&R'

    end subroutine set_symbol_c1

  end subroutine get_irr_rep

! =======================================================================
!  Reference: 
!  https://www.cryst.ehu.es/cgi-bin/cryst/programs/representations_point.pl?tipogrupo=dbg
!
! =======================================================================
  subroutine get_irr_rep_dblegrp( point_group, nsym, nrep, rep, symbol_irrep )
    character(len=3), intent(in) :: point_group
    integer, intent(out) :: nsym, nrep
    complex(kind=CMPLDP), intent(out) :: rep( nsym0, 2, nrep0_dblegrp )
    character(len=4), intent(out) :: symbol_irrep( nrep0_dblegrp ) 

    integer :: i, j
    
    rep = 0.0d0
    symbol_irrep( 1:nrep0_dblegrp ) = 'NON'

    select case(point_group)  
    case('Oh ') !  1 
       call set_chartab_oh( rep )
       call set_symbol_oh( nrep, symbol_irrep )
       nsym = 48
    case('O  ') !  2
       call set_chartab_o( rep )
       call set_symbol_o( nrep, symbol_irrep )
       nsym = 24
    case('Td ') !  3
       call set_chartab_td( rep )
       call set_symbol_td( nrep, symbol_irrep )
       nsym = 24
    case('S4 ') !  11
       call set_chartab_s4( rep )
       call set_symbol_s4( nrep, symbol_irrep )
       nsym = 4
    case('D2h') ! 13
       call set_chartab_d2h( rep )
       call set_symbol_d2h( nrep, symbol_irrep )
       nsym = 8
    case('D2 ') ! 14
       call set_chartab_d2( rep )
       call set_symbol_d2( nrep, symbol_irrep )
       nsym = 4
    case('C2v') ! 15
       call set_chartab_c2v( rep )
       call set_symbol_c2v( nrep, symbol_irrep )
       nsym = 4
    end select

    return
    
    Do i=1, nrep
       write(*,*) "irep = ", i,  symbol_irrep(i)
       Do j=1, nsym
          write(*,'(I5,2F8.4)') j, rep(j,1,i), rep(j,2,i)
       End do
       write(*,*)
    End Do
    stop
    
  contains

    subroutine set_chartab_oh( rep )
      complex(kind=CMPLDP), intent(out) :: rep( nsym0, 2, nrep0_dblegrp )

      integer, parameter :: nsym = 48
      integer, parameter :: nsym2 = nsym/2
      integer, parameter :: nrep = 16
      integer, parameter :: nrep2 = nrep/2
      integer :: i
      
      call set_chartab_o(rep)
      Do i=nsym2+1, nsym
         rep( i, :, 1:nrep2 ) = rep(i-nsym2,:,1:nrep2) 
      End Do
      Do i=nrep2, 1, -1
         rep( 1:nsym2, :, 2*i )        =  rep( 1:nsym2, :, i )
         rep( nsym2+1:nsym, :, 2*i )   = -rep( 1:nsym2, :, i )
         rep( 1:nsym2, :, 2*i-1 )      =  rep( 1:nsym2, :, i )
         rep( nsym2+1:nsym, :, 2*i-1 ) =  rep( 1:nsym2, :, i )
      End Do
    end subroutine set_chartab_oh

    subroutine set_symbol_oh( nrep, symbol_irrep )
      integer, intent(out) :: nrep
      character(len=4), intent(out) :: symbol_irrep(16)

      nrep = 16
      symbol_irrep(1)  = 'GM1+';      symbol_irrep(2)  = 'GM1-' 
      symbol_irrep(3)  = 'GM2+';      symbol_irrep(4)  = 'GM2-'
      symbol_irrep(5)  = 'GM3+';      symbol_irrep(6)  = 'GM3-'
      symbol_irrep(7)  = 'GM4+';      symbol_irrep(8)  = 'GM4-'
      symbol_irrep(9)  = 'GM5+';      symbol_irrep(10) = 'GM5-'
      symbol_irrep(11) = 'GM6+';      symbol_irrep(12) = 'GM6-'
      symbol_irrep(13) = 'GM7+';      symbol_irrep(14) = 'GM7-'
      symbol_irrep(15) = 'GM8+';      symbol_irrep(16) = 'GM8-'
    end subroutine set_symbol_oh

    subroutine set_chartab_o( rep )
      complex(kind=CMPLDP), intent(out) :: rep( nsym0, 2, nrep0_dblegrp )

      integer, parameter :: nsym = 24
! GM1 = A1
      rep( 1:nsym, :, 1 ) = 1    
! GM2 = A2
      rep(     1, :, 2 ) =  1           ! E
      rep(  2: 4, :, 2 ) =  1           ! C2xyz
      rep(  5:12, :, 2 ) =  1           ! C3
      rep( 13:18, :, 2 ) = -1           ! C2a,..
      rep( 19:24, :, 2 ) = -1           ! C4...
! GM3 = E
      rep(     1, :, 3 ) =  2           ! E
      rep(  2: 4, :, 3 ) =  2           ! C2xyz
      rep(  5:12, :, 3 ) = -1           ! C3
      rep( 13:18, :, 3 ) =  0           ! C2a,..
      rep( 19:24, :, 3 ) =  0           ! C4...
! GM4 = T1
      rep(     1, :, 4 ) =  3           ! E
      rep(  2: 4, :, 4 ) = -1           ! C2xyz
      rep(  5:12, :, 4 ) =  0           ! C3
      rep( 13:18, :, 4 ) = -1           ! C2a,..
      rep( 19:24, :, 4 ) =  1           ! C4...
! GM5 = T2
      rep(     1, :, 5 ) =  3           ! E
      rep(  2: 4, :, 5 ) = -1           ! C2xyz
      rep(  5:12, :, 5 ) =  0           ! C3
      rep( 13:18, :, 5 ) =  1           ! C2a,..
      rep( 19:24, :, 5 ) = -1           ! C4...
! GM6 = E1bar
      rep(     1, 1, 6 ) =  2           ! E
      rep(     1, 2, 6 ) = -2           ! E *R
      rep(  2: 4, :, 6 ) =  0           ! C2xyz
      rep(  5:12, 1, 6 ) =  1           ! C3
      rep(  5:12, 2, 6 ) = -1           ! C3 *R
      rep( 13:18, :, 6 ) =  0           ! C2a,..
      rep( 19:24, 1, 6 ) =  sqrt(2.)    ! C4...
      rep( 19:24, 2, 6 ) = -sqrt(2.)    ! C4... *R
! GM7 = E2bar
      rep(     1, 1, 7 ) =  2           ! E
      rep(     1, 2, 7 ) = -2           ! E *R
      rep(  2: 4, :, 7 ) =  0           ! C2xyz
      rep(  5:12, 1, 7 ) =  1           ! C3
      rep(  5:12, 2, 7 ) = -1           ! C3 *R
      rep( 13:18, :, 7 ) =  0           ! C2a,..
      rep( 19:24, 1, 7 ) = -sqrt(2.)    ! C4...
      rep( 19:24, 2, 7 ) =  sqrt(2.)    ! C4... *R
! GM8 = Fbar
      rep(     1, 1, 8 ) =  4           ! E
      rep(     1, 2, 8 ) = -4           ! E *R
      rep(  2: 4, :, 8 ) =  0           ! C2xyz
      rep(  5:12, 1, 8 ) = -1           ! C3
      rep(  5:12, 2, 8 ) =  1           ! C3 *R
      rep( 13:18, :, 8 ) =  0           ! C2a,..
      rep( 19:24, :, 8 ) =  0           ! C4...

    end subroutine set_chartab_o

    subroutine set_symbol_o( nrep, symbol_irrep )
      integer, intent(out) :: nrep
      character(len=4), intent(out) :: symbol_irrep(8)

      nrep = 8
      symbol_irrep(1) = 'GM1 ';      symbol_irrep(2) = 'GM2 '
      symbol_irrep(3) = 'GM3 ';      symbol_irrep(4) = 'GM4 '
      symbol_irrep(5) = 'GM5 ';      symbol_irrep(6) = 'GM6 '
      symbol_irrep(7) = 'GM7 ';      symbol_irrep(8) = 'GM8 '
    end subroutine set_symbol_o

    subroutine set_chartab_td( rep )
      complex(kind=CMPLDP), intent(out) :: rep( nsym0, 2, nrep0_dblegrp )

      call set_chartab_o( rep )
    end subroutine set_chartab_td

    subroutine set_symbol_td( nrep, symbol_irrep )
      integer, intent(out) :: nrep
      character(len=4), intent(out) :: symbol_irrep(8)

      nrep = 8
      symbol_irrep(1) = 'GM1 ';      symbol_irrep(2) = 'GM2 '
      symbol_irrep(3) = 'GM3 ';      symbol_irrep(4) = 'GM4 '
      symbol_irrep(5) = 'GM5 ';      symbol_irrep(6) = 'GM6 '
      symbol_irrep(7) = 'GM7 ';      symbol_irrep(8) = 'GM8 '

    end subroutine set_symbol_td

    subroutine set_chartab_s4( rep )
      complex(kind=CMPLDP), intent(out) :: rep( nsym0, 2, nrep0_dblegrp )

      integer, parameter :: nsym = 4
      complex(kind=CMPLDP) :: ztmp1, ztmp2
      
      ztmp1 = ( 1.d0 +zi )/sqrt(2.d0)
      ztmp2 = ( 1.d0 -zi )/sqrt(2.d0)
      
! GM1 = A
      rep(     1, :, 1 ) =  1           ! E
      rep(     4, :, 1 ) =  1           ! C2   0  0  1
      rep(    21, :, 1 ) =  1           ! C4+  0  0  1
      rep(    24, :, 1 ) =  1           ! C4-  0  0  1
! GM2 = B
      rep(     1, :, 2 ) =  1           ! E
      rep(     4, :, 2 ) =  1           ! C2   0  0  1
      rep(    21, :, 2 ) = -1           ! C4+  0  0  1
      rep(    24, :, 2 ) = -1           ! C4-  0  0  1
! GM3 = 2E
      rep(     1, :, 3 ) =  1           ! E
      rep(     4, :, 3 ) = -1           ! C2   0  0  1
      rep(    21, :, 3 ) =  zi          ! C4+  0  0  1
      rep(    24, :, 3 ) = -zi          ! C4-  0  0  1
! GM4 = 1E
      rep(     1, :, 4 ) =  1           ! E
      rep(     4, :, 4 ) = -1           ! C2   0  0  1
      rep(    21, :, 4 ) = -zi          ! C4+  0  0  1
      rep(    24, :, 4 ) =  zi          ! C4-  0  0  1
! GM5 = 2Ebar1
      rep(     1, 1, 5 ) =  1           ! E
      rep(     4, 1, 5 ) = -zi          ! C2   0  0  1
      rep(    21, 1, 5 ) =  ztmp2       ! C4+  0  0  1
      rep(    24, 1, 5 ) =  ztmp1       ! C4-  0  0  1
      rep(     1, 2, 5 ) = -1           ! E *R
      rep(     4, 2, 5 ) =  zi          ! C2   0  0  1 *R
      rep(    21, 2, 5 ) = -ztmp2       ! C4+  0  0  1 *R
      rep(    24, 2, 5 ) = -ztmp1       ! C4-  0  0  1 *R
! GM6 = 1Ebar1
      rep(     1, 1, 6 ) =  1           ! E
      rep(     4, 1, 6 ) =  zi          ! C2   0  0  1
      rep(    21, 1, 6 ) =  ztmp1       ! C4+  0  0  1
      rep(    24, 1, 6 ) =  ztmp2       ! C4-  0  0  1
      rep(     1, 2, 6 ) = -1           ! E *R
      rep(     4, 2, 6 ) = -zi          ! C2   0  0  1 *R
      rep(    21, 2, 6 ) = -ztmp1       ! C4+  0  0  1 *R
      rep(    24, 2, 6 ) = -ztmp2       ! C4-  0  0  1 *R
! GM7 = 2Ebar2
      rep(     1, 1, 7 ) =  1           ! E
      rep(     4, 1, 7 ) = -zi          ! C2   0  0  1
      rep(    21, 1, 7 ) = -ztmp2       ! C4+  0  0  1
      rep(    24, 1, 7 ) = -ztmp1       ! C4-  0  0  1
      rep(     1, 2, 7 ) = -1           ! E *R
      rep(     4, 2, 7 ) =  zi          ! C2   0  0  1 *R
      rep(    21, 2, 7 ) =  ztmp2       ! C4+  0  0  1 *R
      rep(    24, 2, 7 ) =  ztmp1       ! C4-  0  0  1 *R
! GM8 = 1Ebar2
      rep(     1, 1, 8 ) =  1           ! E
      rep(     4, 1, 8 ) =  zi          ! C2   0  0  1
      rep(    21, 1, 8 ) = -ztmp1       ! C4+  0  0  1
      rep(    24, 1, 8 ) = -ztmp2       ! C4-  0  0  1
      rep(     1, 2, 8 ) = -1           ! E *R
      rep(     4, 2, 8 ) = -zi          ! C2   0  0  1 *R
      rep(    21, 2, 8 ) =  ztmp1       ! C4+  0  0  1 *R
      rep(    24, 2, 8 ) =  ztmp2       ! C4-  0  0  1 *R
    end subroutine set_chartab_s4

    subroutine set_symbol_s4( nrep, symbol_irrep )
      integer, intent(out) :: nrep
      character(len=4), intent(out) :: symbol_irrep(8)

      nrep = 8
      symbol_irrep(1) = 'GM1 ';      symbol_irrep(2) = 'GM2 '
      symbol_irrep(3) = 'GM3 ';      symbol_irrep(4) = 'GM4 '
      symbol_irrep(5) = 'GM5 ';      symbol_irrep(6) = 'GM6 '
      symbol_irrep(7) = 'GM7 ';      symbol_irrep(8) = 'GM8 '

    end subroutine set_symbol_s4

    subroutine set_chartab_d2h( rep )
      complex(kind=CMPLDP), intent(out) :: rep( nsym0, 2, nrep0_dblegrp )

      integer, parameter :: nsym = 8
      integer, parameter :: nsym2 = nsym/2
      integer, parameter :: nrep = 10
      integer, parameter :: nrep2 = nrep/2
      integer :: i
      
      call set_chartab_d2(rep)
      Do i=nsym2+1, nsym
         rep( i, :, 1:nrep2 ) = rep(i-nsym2,:,1:nrep2) 
      End Do
      Do i=nrep2, 1, -1
         rep( 1:nsym2, :, 2*i )        =  rep( 1:nsym2, :, i )
         rep( nsym2+1:nsym, :, 2*i )   = -rep( 1:nsym2, :, i )
         rep( 1:nsym2, :, 2*i-1 )      =  rep( 1:nsym2, :, i )
         rep( nsym2+1:nsym, :, 2*i-1 ) =  rep( 1:nsym2, :, i )
      End Do
    end subroutine set_chartab_d2h

    subroutine set_symbol_d2h( nrep, symbol_irrep )
      integer, intent(out) :: nrep
      character(len=4), intent(out) :: symbol_irrep(10)

      nrep = 10
      symbol_irrep(1)  = 'GM1+';      symbol_irrep(2)  = 'GM1-' 
      symbol_irrep(3)  = 'GM2+';      symbol_irrep(4)  = 'GM2-'
      symbol_irrep(5)  = 'GM3+';      symbol_irrep(6)  = 'GM3-'
      symbol_irrep(7)  = 'GM4+';      symbol_irrep(8)  = 'GM4-'
      symbol_irrep(9)  = 'GM5+';      symbol_irrep(10) = 'GM5-'
    end subroutine set_symbol_d2h

    subroutine set_chartab_d2( rep )
      complex(kind=CMPLDP), intent(out) :: rep( nsym0, 2, nrep0_dblegrp )

      integer, parameter :: nsym = 4
      
! GM1 = A1
      rep(     1, :, 1 ) =  1           ! E
      rep(     2, 1, 1 ) =  1           ! C2   1  0  0
      rep(     3, 1, 1 ) =  1           ! C2   0  1  0
      rep(     4, 1, 1 ) =  1           ! C2   0  0  1
! GM2 = B2
      rep(     1, :, 2 ) =  1           ! E
      rep(     2, 1, 2 ) = -1           ! C2   1  0  0
      rep(     3, 1, 2 ) =  1           ! C2   0  1  0
      rep(     4, 1, 2 ) = -1           ! C2   0  0  1
! GM3 = B1
      rep(     1, :, 3 ) =  1           ! E
      rep(     2, 1, 3 ) = -1           ! C2   1  0  0
      rep(     3, 1, 3 ) = -1           ! C2   0  1  0
      rep(     4, 1, 3 ) =  1           ! C2   0  0  1
! GM4 = B3
      rep(     1, :, 4 ) =  1           ! E
      rep(     2, 1, 4 ) =  1           ! C2   1  0  0
      rep(     3, 1, 4 ) = -1           ! C2   0  1  0
      rep(     4, 1, 4 ) = -1           ! C2   0  0  1
! GM5 = Ebar
      rep(     1, 1, 5 ) =  2           ! E
      rep(     1, 2, 5 ) = -2           ! E *R
    end subroutine set_chartab_d2

    subroutine set_symbol_d2( nrep, symbol_irrep )
      integer, intent(out) :: nrep
      character(len=4), intent(out) :: symbol_irrep(5)

      nrep = 5
      symbol_irrep(1) = 'GM1 ';      symbol_irrep(2) = 'GM2 '
      symbol_irrep(3) = 'GM3 ';      symbol_irrep(4) = 'GM4 '
      symbol_irrep(5) = 'GM5 '; 
    end subroutine set_symbol_d2

    subroutine set_chartab_c2v( rep )
      complex(kind=CMPLDP), intent(out) :: rep( nsym0, 2, nrep0_dblegrp )

      integer, parameter :: nsym = 4
      
! GM1 = A1
      rep(     1, :, 1 ) =  1           ! E
      rep(     4, :, 1 ) =  1           ! C2   0  0  1
      rep( 24 +2, :, 1 ) =  1           ! C2   1  0  0 *Inv
      rep( 25 +2, :, 1 ) =  1           ! C2   0  1  0 *Inv
! GM2 = B1
      rep(     1, :, 2 ) =  1           ! E
      rep(     4, :, 2 ) = -1           ! C2   0  0  1
      rep( 24 +2, :, 2 ) = -1           ! C2   1  0  0 *Inv
      rep( 25 +2, :, 2 ) =  1           ! C2   0  1  0 *Inv
! GM3 = A2
      rep(     1, :, 3 ) =  1           ! E
      rep(     4, :, 3 ) = -1           ! C2   0  0  1
      rep( 24 +2, :, 3 ) =  1           ! C2   1  0  0 *Inv
      rep( 25 +2, :, 3 ) = -1           ! C2   0  1  0 *Inv
! GM4 = B2
      rep(     1, :, 4 ) =  1           ! E
      rep(     4, :, 4 ) =  1           ! C2   0  0  1
      rep( 24 +2, :, 4 ) = -1           ! C2   1  0  0 *Inv
      rep( 25 +2, :, 4 ) = -1           ! C2   0  1  0 *Inv
! GM5 = barE
      rep(     1, 1, 5 ) =  2           ! E
      rep(     1, 2, 5 ) = -2           ! R
    end subroutine set_chartab_c2v

    subroutine set_symbol_c2v( nrep, symbol_irrep )
      integer, intent(out) :: nrep
      character(len=4), intent(out) :: symbol_irrep(5)

      nrep = 5
      symbol_irrep(1) = 'GM1 ';      symbol_irrep(2) = 'GM2 '
      symbol_irrep(3) = 'GM3 ';      symbol_irrep(4) = 'GM4 '
      symbol_irrep(5) = 'GM5 '; 
    end subroutine set_symbol_c2v

  end subroutine get_irr_rep_dblegrp

end module m_Representation

