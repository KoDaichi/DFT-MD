#ifdef HIUX
*OPTION MP(P(0))
#endif
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: heap_sort, swap
!
!  AUTHOR: Takenori Yamamotoi   June/01/2003
!  
!  Contact address :  Phase System Consortium
!                     E-mail: phase_system@nims.go.jp URL https://azuma.nims.go.jp
!  
!
!
!=======================================================================
!
! coded by Takenori YAMAMOTO (Univ. Tokyo) for gen_rgrid in m_SpecialKpoints
! at Jun. 1, 2003
!
subroutine  heap_sort(m,n,a,key)
! $Id: heap_sort.F90 570 2017-04-21 20:34:50Z yamasaki $
!
! Heap sorting
! a(M,N) : array of M-componet vectors
! key    : key componet for sorting 
!
implicit none

integer(4), intent(in) :: m,n,key
real(8), intent(inout) :: a(m,n)

integer(4) :: i,j,k,itmp
real(8) :: tmp(m)

!*** initiation: heapify a ***
do i=2,n
  tmp(1:m)=a(1:m,i)
  j=i
10 itmp=j/2
  if(a(key,itmp).ge.tmp(key)) go to 20
  a(1:m,j)=a(1:m,itmp)
  j=itmp
  if(j.gt.1) go to 10
20 a(1:m,j)=tmp(1:m)
end do
!*** to be fully sort ***
do k=n-1,1,-1
  call swap(a(1,1),a(1,k+1))
  tmp(1:m)=a(1:m,1)
  j=1
  itmp=2
30 if(itmp.gt.k) go to 40
  if(itmp.lt.k) then
    if(a(key,itmp+1).gt.a(key,itmp)) itmp=itmp+1
  end if
  if(tmp(key).ge.a(key,itmp)) go to 40
  a(1:m,j)=a(1:m,itmp)
  j=itmp
  itmp=j*2
  go to 30
40 a(1:m,j)=tmp(1:m)
end do

return

contains

!*** swap ***
subroutine swap(a,b)

implicit none

real(8), intent(inout) :: a(m),b(m)
real(8) :: tmp(m)

tmp(1:m)=a(1:m)
a(1:m)=b(1:m)
b(1:m)=tmp(1:m)

return
end subroutine swap

end subroutine heap_sort
