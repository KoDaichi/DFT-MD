!================================================
!  Software name : STM
!  Subroutine(s) : rsreal, tsort, cmpstr, skip_lines, getttr, hpsort,
!                  error
!  Author(s)     : Koichi Kato and Takahiro Yamasaki (June 7, 2004)
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
!
! $Id: bottom_Subroutines.F90,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $
!                           @(#)bottom_Subroutines.F90 1.1 01/03/13 00:34:19
subroutine rsreal(n, obj)
  integer, intent(in) :: n
  integer, parameter :: DP = kind(1.d0)
  real(kind=DP), intent(out), dimension(n) :: obj

  integer :: i
  do i = 1, n
     obj(i) = 0.d0
  enddo
end subroutine rsreal

subroutine tsort(ecpu,msbrnm,iorder)
  implicit none
  integer, parameter :: DP = kind(1.d0)
  integer      , intent(in)    :: msbrnm
  real(kind=DP), intent(in), dimension(msbrnm) :: ecpu
  integer      , intent(out),dimension(msbrnm) :: iorder

  integer :: i,ib,ia,it
  do i = 1, msbrnm
     iorder(i) = i
  enddo
  do ib = 2, msbrnm
     do ia = 1, ib-1
        if(ecpu(iorder(ia)).lt.ecpu(iorder(ib))) then
           it = iorder(ia)
           iorder(ia) = iorder(ib)
           iorder(ib) = it
        endif
     enddo
  enddo
end subroutine tsort

subroutine cmpstr(name1,name2,nlen,ident)
  implicit none
  integer, intent(in) :: nlen
  character*(*), intent(in) :: name1, name2
  logical, intent(out) ::  ident

  character :: c
  integer   :: ndif,i

  ndif = ichar('A') - ichar('a')

  ident = .true.
  do i = 1, nlen
     c = name2(i:i)
     if(ichar(c).ge.ichar('a').and.ichar(c).le.ichar('z')) then
        c = char(ichar(c) + ndif)
     endif
     if(name1(i:i) .ne. c) goto 1001
  enddo
  goto 1002
1001 ident = .false.
1002 continue
end subroutine cmpstr

subroutine skip_lines(nfinp, nline)
  integer, intent(in) :: nfinp
  integer, intent(in) :: nline
  integer :: i
  do i = 1, nline
     read(nfinp,*)
  enddo
end subroutine skip_lines

subroutine getttr (rltv,ttr)
  implicit none
  integer, parameter :: DP = kind(1.d0)
  real(kind=DP), intent(in), dimension(3,3) :: rltv
  real(kind=DP), intent(out),dimension(6)   :: ttr

  integer i, j, jj

  call rsreal(6,ttr)
  do j = 1, 3
     jj = mod(j+1,3)
     if(jj.eq.0) jj = 3
     do i = 1, 3
        ttr(j) = ttr(j) + rltv(i,j)*rltv(i,j)
        ttr(j+3) = ttr(j+3) + 2.d0*rltv(i,j)*rltv(i,jj)
     enddo
  enddo
end subroutine getttr


subroutine hpsort(kn,n,bxyz,br)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!     THIS SUBROUTINE SORTS THE ELEMENTS OF ARRAY BR IN ASCENDING
!     ORDER. THE CODE IS BASED UPON HEAPSORT ALGORITHM WHICH HAS
!     N*LOG(N) SCALING BEHAVIOUR, PARTICULARLY FAVOURABLE FOR LARGE
!     VECTORS BR
!                              STEFAN BL"UGEL, ISSP, NOV 1989
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  implicit none
  integer, intent(in) :: kn, n
  real(kind=kind(1.d0)), intent(inout), dimension(kn,3) :: bxyz
  real(kind=kind(1.d0)), intent(inout), dimension(kn)   :: br
  integer   l,i,ii,iheap
  real*8    brr,bxx,byy,bzz

!=====> BUILD-UP OF THE HEAP
!
!-----> LOOP OVER ALL HIERACHICAL LEVELS OF THE HEAP TREE

  Loop_L :  do L = N/2 , 1 , -1
     brr = br(l)
     bxx = bxyz(l,1)
     byy = bxyz(l,2)
     bzz = bxyz(l,3)
     I   = L
     II  = L + L

!-----> GO DOWN ALL THE HIERACHICAL LEVELS OF THE HEAP TREE

20   IF ( II <= N ) THEN

!-----> COMPARE NEIGHBOURING ELEMENTS
        IF ( II < N ) THEN
           IF ( BR(II) < BR(II+1) ) II = II + 1
        ENDIF

!-----> COMPARE THE ELEMENTS OF TWO HIRACHICAL LEVELS
!       PROMOTE LARGER ONE, DEMOTE SMALLER ONE
        IF ( BRR < BR(II) ) THEN
           BR(I) = BR(II)
           bxyz(i,1) = bxyz(ii,1)
           bxyz(i,2) = bxyz(ii,2)
           bxyz(i,3) = bxyz(ii,3)
           I     = II
           II    = II + II
        ELSE

!-----> THIS PART OF THE TREE IS ORDERED , STOP BY PUTTING II=N+1
           II    = N + 1
        END IF
        GO TO 20
     END IF
!-----> PUT ELEMENTS IN THE PROPER SLOT
     br(i) = brr
     bxyz(i,1) = bxx
     bxyz(i,2) = byy
     bxyz(i,3) = bzz
  enddo Loop_L

!=====> NOW COLLECT ALL ELEMENTS FROM THE HEAP

  Loop_IHEAP : do IHEAP = N , 1 , -1

     brr = br(iheap)
     bxx = bxyz(iheap,1)
     byy = bxyz(iheap,2)
     bzz = bxyz(iheap,3)

!-----> THE FIRST ELEMENT IS ALWAYS THE LARGEST
     br(iheap) = br(1)
     bxyz(iheap,1) = bxyz(1,1)
     bxyz(iheap,2) = bxyz(1,2)
     bxyz(iheap,3) = bxyz(1,3)
!-----> NOW LARGEST ELEMENT OF ALL BR(I) WITH 1<=I<=IHEAP IS STORED

     I   = 1
     II  = 2

!-----> NOW GENERATE LARGEST ELEMENT OF BR(I) WITH 1<=I<=IHEAP-1

40   IF ( II <= IHEAP - 1 ) THEN

!-----> COMPARE NEIGHBOURING ELEMENTS
        IF ( II < IHEAP - 1 ) THEN
           IF ( BR(II) < BR(II+1) ) II = II + 1
        ENDIF

!-----> PROMOTE EACH ELEMENT OF THE TWIG OF BR UNTIL BRR > BR(I)
        IF ( BRR < BR(II) ) THEN
           br(i) = br(ii)
           bxyz(i,1) = bxyz(ii,1)
           bxyz(i,2) = bxyz(ii,2)
           bxyz(i,3) = bxyz(ii,3)
           I     = II
           II    = II + II
        ELSE

!-----> THIS PART OF THE TREE IS PROMOTED , STOP BY PUTTING II=IHEAP+1
           II    = IHEAP + 1
        END IF
        GO TO 40
     END IF
!-----> PUT ELEMENTS IN THE PROPER SLOT
     br(i) = brr
     bxyz(i,1) = bxx
     bxyz(i,2) = byy
     bxyz(i,3) = bzz
  enddo Loop_IHEAP
end subroutine hpsort

!!$subroutine eqivvl(value)
!!$c                by T.Yamasaki (JRCAT,ATP)
!!$  use m_Const_Parameters, only : DP
!!$  implicit none
!!$  real(kind=DP), intent(inout) :: value
!!$
!!$  logical flag
!!$#ifdef VPP
!!$#ifdef PARA
!!$  include 'nproc.inc'
!!$!xocl processor pe(npe)
!!$!xocl subprocessor pes(npe)=pe(1:npe)
!!$  include 'para_p.inc'
!!$  include 'inddef_p.inc'
!!$#endif
!!$#endif
!!$
!!$  flag = .false.
!!$!xocl spread do/ind_keg
!!$!xocl index 1
!!$  flag = .true.
!!$!xocl end spread
!!$!xocl broadcast(value) (flag)
!!$end subroutine eqivvl

subroutine error
  write(6,*) ' ! error '
end subroutine error


