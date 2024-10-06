!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: gsort, spsort, shellsort, sort_gvec_within_shell_by_heap,
!             sift_down, sort_gvec_heap, sift_down, swap_i, 
!             sort_gvec_quick, sort_gvec_simple
!
!  AUTHOR(S): T. Yamasaki, M. Okamoto,    August/20/2003
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

! $Id: b_PlaneWaveBasisSet.F90 570 2017-04-21 20:34:50Z yamasaki $

!#####################################################################

subroutine sort_gvec_simple(ttr,k,n,ngvec)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: sort_gvec_simple
!
!  AUTHOR(S): M. Okamoto   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

  use m_Const_Parameters, only : DP
  implicit none
  real(kind=DP),intent(in),dimension(6) :: ttr
  integer,intent(in)                    :: k, n
  integer,intent(inout),dimension(k,3)  :: ngvec
  real(kind=DP),allocatable :: norm_gvec(:)
  integer,allocatable :: idx_gvec(:), ngvec_tmp(:,:)
  integer :: i, j, ipos, jpos, ia, ib, ic
  real(kind=DP) :: sum, norm_gi, norm_gj
 !++++++++++++++++++++++++++++++++++++++++++++++++++
  allocate(norm_gvec(n),idx_gvec(n),ngvec_tmp(k,3))
 !++++++++++++++++++++++++++++++++++++++++++++++++++
  do i = 1,n
     ia = ngvec(i,1) ; ib = ngvec(i,2) ; ic = ngvec(i,3)
     sum = ttr(1)*ia*ia + ttr(2)*ib*ib + ttr(3)*ic*ic &
         + ttr(4)*ia*ib + ttr(5)*ib*ic + ttr(6)*ic*ia
     norm_gvec(i) = sqrt(sum)
     idx_gvec (i) = i
  end do
  do i = 1,n-1
     ipos    = idx_gvec(i)
     norm_gi = norm_gvec(ipos)
     do j = i+1,n
        jpos    = idx_gvec(j)
        norm_gj = norm_gvec(jpos)
        if (norm_gi > norm_gj) then
           idx_gvec(i) = jpos
           idx_gvec(j) = ipos
           ipos    = jpos
           norm_gi = norm_gj
        end if
     end do
  end do
  ngvec_tmp(:,:) = ngvec(:,:)
  do i = 1,n
     ngvec(i,:) = ngvec_tmp(idx_gvec(i),:)
  end do
 !+++++++++++++++++++++++++++++++++++++++++++
  deallocate(norm_gvec,idx_gvec,ngvec_tmp)
 !+++++++++++++++++++++++++++++++++++++++++++
end subroutine sort_gvec_simple


subroutine sort_gvec_quick(ttr,k,n,ngvec)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: sort_gvec_quick
!
!  AUTHOR(S): M. Okamoto   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

  use m_Const_Parameters, only : DP
  implicit none
  real(kind=DP),intent(in),dimension(6) :: ttr
  integer,intent(in)                    :: k, n
  integer,intent(inout),dimension(k,3)  :: ngvec
  real(kind=DP),allocatable :: norm_gvec(:)
  integer,allocatable :: idx_gvec(:), ngvec_tmp(:,:)
  integer :: i, j, ia, ib, ic
  real(kind=DP) :: sum, v0
  integer,parameter :: NN=15, NSTACK=50
  integer :: jstack, l, r, ipos, jpos, ll, kk, istack(NSTACK)
 !++++++++++++++++++++++++++++++++++++++++++++++++++
  allocate(norm_gvec(n),idx_gvec(n),ngvec_tmp(k,3))
 !++++++++++++++++++++++++++++++++++++++++++++++++++
  do i = 1,n
     ia = ngvec(i,1) ; ib = ngvec(i,2) ; ic = ngvec(i,3)
     sum = ttr(1)*ia*ia + ttr(2)*ib*ib + ttr(3)*ic*ic &
         + ttr(4)*ia*ib + ttr(5)*ib*ic + ttr(6)*ic*ia
     norm_gvec(i) = sqrt(sum)
     idx_gvec (i) = i
  end do
  jstack = 0
  l = 1
  r = n
  do
     if (r-1 < NN) then
        do j = l+1,r
           jpos = idx_gvec(j)
           v0 = norm_gvec(jpos)
           do i = j-1,l,-1
              ipos = idx_gvec(i)
              if (norm_gvec(ipos) <= v0) then
                 exit
              end if
              idx_gvec(i+1) = idx_gvec(i)
           end do
           idx_gvec(i+1) = jpos
        end do
        if (jstack == 0) then
           return
        end if
        r = istack(jstack  )
        l = istack(jstack-1)
        jstack = jstack-2
     else
        kk = (l+r)/2
        call swap_i(idx_gvec(kk),idx_gvec(l+1))
        if (norm_gvec(idx_gvec(l)) > norm_gvec(idx_gvec(r))) then
           call swap_i(idx_gvec(l),idx_gvec(r))
        end if
        if (norm_gvec(idx_gvec(l+1)) > norm_gvec(idx_gvec(r))) then
           call swap_i(idx_gvec(l+1),idx_gvec(r))
        end if
        if (norm_gvec(idx_gvec(l)) > norm_gvec(idx_gvec(l+1))) then
           call swap_i(idx_gvec(l),idx_gvec(l+1))
        end if
        i = l+1
        j = r
        ll = idx_gvec(l+1)
        v0 = norm_gvec(ll)
        do
           do
              i = i+1
              if (norm_gvec(idx_gvec(i)) >= v0) then
                 exit
              end if
           end do
           do
              j = j-1
              if (norm_gvec(idx_gvec(j)) <= v0) then
                 exit
              end if
           end do
           if (j < i) then
              exit
           end if
           call swap_i(idx_gvec(i),idx_gvec(j))
        end do
        idx_gvec(l+1) = idx_gvec(j)
        idx_gvec(j) = ll
        jstack = jstack+2
        if (jstack > NSTACK) then
           write(*,*) 'NSTACK too small in sort_gvec_quick'
           stop
        end if
        if (r-i+1 >= j-1) then
           istack(jstack  ) = r
           istack(jstack-1) = i
           r = j-1
        else
           istack(jstack  ) = j-1
           istack(jstack-1) = l
           l = i
        end if
     end if
  end do
  ngvec_tmp(:,:) = ngvec(:,:)
  do i = 1,n
     ngvec(i,:) = ngvec_tmp(idx_gvec(i),:)
  end do
 !+++++++++++++++++++++++++++++++++++++++++
  deallocate(norm_gvec,idx_gvec,ngvec_tmp)
 !+++++++++++++++++++++++++++++++++++++++++
end subroutine sort_gvec_quick


subroutine swap_i(a,b)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: swap_i
!
!  AUTHOR(S): M. Okamoto    August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

  integer,intent(inout) :: a,b
  integer :: tmp
  tmp = a ; a = b ; b = tmp
end subroutine swap_i


subroutine sort_gvec_heap(ttr,k,n,ngvec)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: 
!
!  AUTHOR(S): M. Okamoto   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

  use m_Const_Parameters, only : DP
  implicit none
  real(kind=DP),intent(in),dimension(6) :: ttr
  integer,intent(in)                    :: k, n
  integer,intent(inout),dimension(k,3)  :: ngvec
  real(kind=DP) :: sum
  real(kind=DP),allocatable :: norm_gvec(:)
  integer,allocatable :: idx_gvec(:), ngvec_tmp(:,:)
  integer :: i, j, ia, ib, ic
 !++++++++++++++++++++++++++++++++++++++++++++++++++
  allocate(norm_gvec(n),idx_gvec(n),ngvec_tmp(k,3))
 !++++++++++++++++++++++++++++++++++++++++++++++++++
  do i = 1,n
     ia = ngvec(i,1) ; ib = ngvec(i,2) ; ic = ngvec(i,3)
     sum = ttr(1)*ia*ia + ttr(2)*ib*ib + ttr(3)*ic*ic &
         + ttr(4)*ia*ib + ttr(5)*ib*ic + ttr(6)*ic*ia
     norm_gvec(i) = sqrt(sum)
     idx_gvec (i) = i
  end do
  do i = n/2,1,-1
     call sift_down(i,n)
  end do
  do i = n,2,-1
     call swap_i(idx_gvec(1),idx_gvec(i))
     call sift_down(1,i-1)
  end do
  ngvec_tmp(:,:) = ngvec(:,:)
  do i = 1,n
     ngvec(i,:) = ngvec_tmp(idx_gvec(i),:)
  end do
 !+++++++++++++++++++++++++++++++++++++++++
  deallocate(norm_gvec,idx_gvec,ngvec_tmp)
 !+++++++++++++++++++++++++++++++++++++++++
!....................................................................
  contains
  subroutine sift_down(l,r)
     integer,intent(in) :: l, r
     integer :: j, j0, l0
     real(kind=DP) :: v0
     l0 = idx_gvec(l)
     v0 = norm_gvec(l0)
     j0 = l
     j  = l+l
     do while (j <= r)
        if (j < r) then
           if (norm_gvec(idx_gvec(j)) < norm_gvec(idx_gvec(j+1))) then
              j = j+1
           end if
        end if
        if (v0 >= norm_gvec(idx_gvec(j))) then
           exit
        end if
        idx_gvec(j0) = idx_gvec(j)
        j0 = j
        j  = j+j
     end do
     idx_gvec(j0) = l0
  end subroutine sift_down
!....................................................................
end subroutine sort_gvec_heap

subroutine sort_gvec_heap2(ttr,k,n,ngvec,gr)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: 
!
!  AUTHOR(S): M. Okamoto   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

  use m_Const_Parameters, only : DP
  implicit none
  real(kind=DP),intent(in),dimension(6) :: ttr
  integer,intent(in)                    :: k, n
  integer,intent(inout),dimension(k,3)  :: ngvec
  real(kind=DP),intent(inout),dimension(k) ::gr

!!$  real(kind=DP) :: sum
  real(kind=DP),allocatable :: norm_gvec(:), gr_tmp(:)
  integer,allocatable :: idx_gvec(:), ngvec_tmp(:,:)
  integer :: i, j, ia, ib, ic
 !++++++++++++++++++++++++++++++++++++++++++++++++++
  allocate(norm_gvec(n),idx_gvec(n),ngvec_tmp(k,3),gr_tmp(k))
 !++++++++++++++++++++++++++++++++++++++++++++++++++
  do i = 1,n
!!$     ia = ngvec(i,1) ; ib = ngvec(i,2) ; ic = ngvec(i,3)
!!$     sum = ttr(1)*ia*ia + ttr(2)*ib*ib + ttr(3)*ic*ic &
!!$         + ttr(4)*ia*ib + ttr(5)*ib*ic + ttr(6)*ic*ia
!!$     norm_gvec(i) = sqrt(sum)
     norm_gvec(i) = gr(i)
     idx_gvec (i) = i
  end do
  do i = n/2,1,-1
     call sift_down(i,n)
  end do
  do i = n,2,-1
     call swap_i(idx_gvec(1),idx_gvec(i))
     call sift_down(1,i-1)
  end do
  ngvec_tmp(:,:) = ngvec(:,:)
  gr_tmp(:)      = gr(:)
  do i = 1,n
     ngvec(i,:) = ngvec_tmp(idx_gvec(i),:)
     gr(i)      = gr_tmp(idx_gvec(i))
  end do
 !+++++++++++++++++++++++++++++++++++++++++
  deallocate(norm_gvec,idx_gvec,ngvec_tmp,gr_tmp)
 !+++++++++++++++++++++++++++++++++++++++++
!....................................................................
  contains
  subroutine sift_down(l,r)
     integer,intent(in) :: l, r
     integer :: j, j0, l0
     real(kind=DP) :: v0
     l0 = idx_gvec(l)
     v0 = norm_gvec(l0)
     j0 = l
     j  = l+l
     do while (j <= r)
        if (j < r) then
           if (norm_gvec(idx_gvec(j)) < norm_gvec(idx_gvec(j+1))) then
              j = j+1
           end if
        end if
        if (v0 >= norm_gvec(idx_gvec(j))) then
           exit
        end if
        idx_gvec(j0) = idx_gvec(j)
        j0 = j
        j  = j+j
     end do
     idx_gvec(j0) = l0
  end subroutine sift_down
!....................................................................
end subroutine sort_gvec_heap2

subroutine sort_gvec_heap3(ttr,k,n,ngvec,gr)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE:
!
!  AUTHOR(S): M. Okamoto   August/20/2003
!
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)
!
!=======================================================================

  use m_Const_Parameters, only : DP
  implicit none
  real(kind=DP),intent(in),dimension(6) :: ttr
  integer,intent(in)                    :: k, n
  integer,intent(inout),dimension(k,3)  :: ngvec
  real(kind=DP),intent(inout),dimension(k) ::gr

  real(kind=DP),allocatable :: gr_tmp(:)
  integer,allocatable :: idx_gvec(:), ngvec_tmp(:,:)
  integer :: i, j, ia, ib, ic, itmp
 !++++++++++++++++++++++++++++++++++++++++++++++++++
  allocate(idx_gvec(n),ngvec_tmp(k,3),gr_tmp(k))
 !++++++++++++++++++++++++++++++++++++++++++++++++++
  do i = 1,n
     idx_gvec (i) = i
  end do
  do i = n/2,1,-1
     call sift_down(i,n)
  end do
  do i = n,2,-1
     itmp = idx_gvec(1)
     idx_gvec(1) = idx_gvec(i)
     idx_gvec(i) = itmp
     call sift_down(1,i-1)
  end do
  ngvec_tmp(:,:) = ngvec(:,:)
  gr_tmp(:)      = gr(:)
  do i = 1,n
     ngvec(i,:) = ngvec_tmp(idx_gvec(i),:)
     gr(i)      = gr_tmp(idx_gvec(i))
  end do
 !+++++++++++++++++++++++++++++++++++++++++
  deallocate(idx_gvec,ngvec_tmp,gr_tmp)
 !+++++++++++++++++++++++++++++++++++++++++
!....................................................................
  contains
  subroutine sift_down(l,r)
     integer,intent(in) :: l, r
     integer :: j, j0, l0
     real(kind=DP) :: v0
     l0 = idx_gvec(l)
     v0 = gr(l0)
     j0 = l
     j  = l+l
     do while (j <= r)
        if (j < r) then
           if (gr(idx_gvec(j)) < gr(idx_gvec(j+1))) then
              j = j+1
           end if
        end if
        if (v0 >= gr(idx_gvec(j))) then
           exit
        end if
        idx_gvec(j0) = idx_gvec(j)
        j0 = j
        j  = j+j
     end do
     idx_gvec(j0) = l0
  end subroutine sift_down
!....................................................................
end subroutine sort_gvec_heap3

! This is a replacement of spsort
subroutine sort_gvec_within_shell_by_heap( &
              n,imin,jmin,jmax,kmin,kmax,na,nb,nc)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: 
!
!  AUTHOR(S): M. Okamoto   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

  use m_Const_Parameters, only : DP
  implicit none
  integer,intent(in)        :: n, imin, jmin, jmax, kmin, kmax
  integer,intent(inout)     :: na(n), nb(n), nc(n)
  real(kind=DP),allocatable :: val(:)
  integer,allocatable       :: idx(:), ntmp(:,:)
  integer                   :: i, b, c
 !++++++++++++++++++++++++++++++++++
  allocate(val(n),idx(n),ntmp(n,3))
 !++++++++++++++++++++++++++++++++++
  b = jmax-jmin+1 ; c = kmax-kmin+1
  do i = 1,n
     val(i) = (na(i)-imin)*b*c + (nb(i)-jmin)*c + (nc(i)-kmin) + 1
     idx(i) = i
  end do
  do i = n/2,1,-1
     call sift_down(i,n)
  end do
  do i = n,2,-1
     call swap_i(idx(1),idx(i))
     call sift_down(1,i-1)
  end do
  ntmp(:,1) = na(:) ; ntmp(:,2) = nb(:) ; ntmp(:,3) = nc(:)
  do i = 1,n
     na(i) = ntmp(idx(i),1)
     nb(i) = ntmp(idx(i),2)
     nc(i) = ntmp(idx(i),3)
  end do
 !+++++++++++++++++++++++++
  deallocate(val,idx,ntmp)
 !+++++++++++++++++++++++++
!....................................................................
  contains
  subroutine sift_down(l,r)
     integer,intent(in) :: l, r
     integer :: j, j0, l0
     real(kind=DP) :: v0
     l0 = idx(l)
     v0 = val(l0)
     j0 = l
     j  = l+l
     do while (j <= r)
        if (j < r) then
           if (val(idx(j)) < val(idx(j+1))) then
              j = j+1
           end if
        end if
        if (v0 >= val(idx(j))) then
           exit
        end if
        idx(j0) = idx(j)
        j0 = j
        j  = j+j
     end do
     idx(j0) = l0
  end subroutine sift_down
!....................................................................
end subroutine sort_gvec_within_shell_by_heap


!#####################################################################


subroutine shellsort(nfout,ipri,ttr,m,n,nG)
! $Id: b_PlaneWaveBasisSet.F90 570 2017-04-21 20:34:50Z yamasaki $
  use m_Const_Parameters, only : DP
  implicit none
  integer, intent(in)                       :: nfout,ipri
  real(kind=DP), intent(in), dimension(6)   :: ttr
  integer, intent(in)                       :: m, n
  integer, intent(inout), dimension(m,3) :: nG

      
  integer :: i, k, counter, m_shellsize = 0, nshell = 1, j
  real(kind=DP), parameter :: DELTA = 1.d-7
  real(kind=DP) :: length_start, length
  integer       :: imin, imax, jmin, jmax, kmin, kmax

  if(ipri >= 2 ) then
     write(nfout,*) ' ----- shell sorting of G-vectors ----'
     write(nfout,9000)
9000 format(' #shell  shellsize    range     length   range-x   range-y   range-z')
  endif
  imin = nG(1,1); imax = imin
  jmin = nG(1,2); jmax = jmin
  kmin = nG(1,3); kmax = kmin
  length_start = ttr(1)*nG(1,1)*nG(1,1) + ttr(2)*nG(1,2)*nG(1,2) &
       &       + ttr(3)*nG(1,3)*nG(1,3) + ttr(4)*nG(1,1)*nG(1,2) &
       &       + ttr(5)*nG(1,2)*nG(1,3) + ttr(6)*nG(1,3)*nG(1,1)
  k = 2
1 continue
  shellloop: do i = k, n
     length =    ttr(1)*nG(i,1)*nG(i,1) + ttr(2)*nG(i,2)*nG(i,2) &
          &    + ttr(3)*nG(i,3)*nG(i,3) + ttr(4)*nG(i,1)*nG(i,2) &
          &    + ttr(5)*nG(i,2)*nG(i,3) + ttr(6)*nG(i,3)*nG(i,1)
     if(length > length_start + DELTA .or. i == n) then
        counter = i - k + 1
        if(i == n) counter = i - k + 2
        if(counter > m_shellsize) m_shellsize = counter
        exit shellloop
     else
        if(imin > nG(i,1)) imin = nG(i,1)
        if(imax < nG(i,1)) imax = nG(i,1)
        if(jmin > nG(i,2)) jmin = nG(i,2)
        if(jmax < nG(i,2)) jmax = nG(i,2)
        if(kmin > nG(i,3)) kmin = nG(i,3)
        if(kmax < nG(i,3)) kmax = nG(i,3)
     endif
  enddo shellloop
#ifdef _RPL_SPSORT_
 ! (M.Okamoto)
  call sort_gvec_within_shell_by_heap( &
              counter,imin,jmin,jmax,kmin,kmax, &
              nG(k-1,1),nG(k-1,2),nG(k-1,3))
#else
  call spsort(counter,imin,jmin,jmax,kmin,kmax &
       & ,nG(k-1,1),nG(k-1,2),nG(k-1,3))
#endif
  if(ipri >= 2) write(nfout,9001) nshell, counter,k-1,k-1+counter-1&
       &,dsqrt(length_start),imin,imax,jmin,jmax,kmin,kmax
  if(ipri >= 3) then
     do j = 1, counter
        write(nfout,*) k-1+j-1,nG(k-1+j-1,1),nG(k-1+j-1,2),nG(k-1+j-1,3)
     end do
  end if
  length_start = length
  k = i + 1
  imin = nG(i,1); imax = imin; jmin = nG(i,2); jmax = jmin
  kmin = nG(i,3); kmax = kmin
9001 format(' ',i6,'  ',i4,' (',i6,':',i6,')', f10.5, ' (',i4,':',i4,',',i4,':',i4,',',i4,':',i4,')')
  
  nshell = nshell + 1
  if(k <= n) goto 1
  if(ipri >= 2) write(nfout,*) ' ! m_shellsize = ',m_shellsize
end subroutine shellsort

subroutine shellsort2(nfout,ipri,ttr,m,n,nG,gr)
! $Id: b_PlaneWaveBasisSet.F90 570 2017-04-21 20:34:50Z yamasaki $
  use m_Const_Parameters, only : DP
  implicit none
  integer, intent(in)                       :: nfout,ipri
  real(kind=DP), intent(in), dimension(6)   :: ttr
  integer, intent(in)                       :: m, n
  integer, intent(inout), dimension(m,3) :: nG
  real(kind=DP), intent(inout), dimension(m) :: gr

      
  integer :: i, k, counter, m_shellsize = 0, nshell = 1, j
  real(kind=DP), parameter :: DELTA = 1.d-7
  real(kind=DP) :: length_start, length
  integer       :: imin, imax, jmin, jmax, kmin, kmax

  if(ipri >= 2 ) then
     write(nfout,*) ' ----- shell sorting of G-vectors ----'
     write(nfout,9000)
9000 format(' #shell  shellsize    range     length   range-x   range-y   range-z')
  endif
  imin = nG(1,1); imax = imin
  jmin = nG(1,2); jmax = jmin
  kmin = nG(1,3); kmax = kmin
!!$  length_start = ttr(1)*nG(1,1)*nG(1,1) + ttr(2)*nG(1,2)*nG(1,2) &
!!$       &       + ttr(3)*nG(1,3)*nG(1,3) + ttr(4)*nG(1,1)*nG(1,2) &
!!$       &       + ttr(5)*nG(1,2)*nG(1,3) + ttr(6)*nG(1,3)*nG(1,1)
  length_start = gr(1)
  k = 2
1 continue
  shellloop: do i = k, n
!!$     length =    ttr(1)*nG(i,1)*nG(i,1) + ttr(2)*nG(i,2)*nG(i,2) &
!!$          &    + ttr(3)*nG(i,3)*nG(i,3) + ttr(4)*nG(i,1)*nG(i,2) &
!!$          &    + ttr(5)*nG(i,2)*nG(i,3) + ttr(6)*nG(i,3)*nG(i,1)
     length =    gr(i)
     if(length > length_start + DELTA .or. i == n) then
        counter = i - k + 1
        if(i == n) counter = i - k + 2
        if(counter > m_shellsize) m_shellsize = counter
        exit shellloop
     else
        if(imin > nG(i,1)) imin = nG(i,1)
        if(imax < nG(i,1)) imax = nG(i,1)
        if(jmin > nG(i,2)) jmin = nG(i,2)
        if(jmax < nG(i,2)) jmax = nG(i,2)
        if(kmin > nG(i,3)) kmin = nG(i,3)
        if(kmax < nG(i,3)) kmax = nG(i,3)
     endif
  enddo shellloop
#ifdef _RPL_SPSORT_
 ! (M.Okamoto)
  call sort_gvec_within_shell_by_heap( &
              counter,imin,jmin,jmax,kmin,kmax, &
              nG(k-1,1),nG(k-1,2),nG(k-1,3))
#else
  call spsort(counter,imin,jmin,jmax,kmin,kmax &
       & ,nG(k-1,1),nG(k-1,2),nG(k-1,3))
#endif
  if(ipri >= 2) write(nfout,9001) nshell, counter,k-1,k-1+counter-1&
       &,dsqrt(length_start),imin,imax,jmin,jmax,kmin,kmax
  if(ipri >= 3) then
     do j = 1, counter
        write(nfout,*) k-1+j-1,nG(k-1+j-1,1),nG(k-1+j-1,2),nG(k-1+j-1,3)
     end do
  end if
  length_start = length
  k = i + 1
  imin = nG(i,1); imax = imin; jmin = nG(i,2); jmax = jmin
  kmin = nG(i,3); kmax = kmin
9001 format(' ',i6,'  ',i4,' (',i6,':',i6,')', f10.5, ' (',i4,':',i4,',',i4,':',i4,',',i4,':',i4,')')
  
  nshell = nshell + 1
  if(k <= n) goto 1
  if(ipri >= 2) write(nfout,*) ' ! m_shellsize = ',m_shellsize
end subroutine shellsort2

subroutine spsort(n,imin,jmin,jmax,kmin,kmax,na,nb,nc)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: spsort
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

  implicit none
  integer, intent(in) :: n, imin,jmin,jmax,kmin,kmax
  integer, intent(inout), dimension(n) :: na,nb,nc
!C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!C     THIS SUBROUTINE SORTS THE ELEMENTS OF ARRAY BR IN ASCENDING
!C     ORDER. THE CODE IS BASED UPON HEAPSORT ALGORITHM WHICH HAS
!C     N*LOG(N) SCALING BEHAVIOUR, PARTICULARLY FAVOURABLE FOR LARGE
!C     VECTORS BR
!C                              STEFAN BL"UGEL, ISSP, NOV 1989
!c     interface is changed by T.Yamasaki at 28th Jun. 1994
!c     changed by T. Yamsaki at 17th Apr.1998
!c
!C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  integer b, c
  integer  l,i,ii,iheap
  integer ia,ib,ic
  integer p1,p2,p3

  b = jmax-jmin+1
  c = kmax-kmin+1

!C=====> BUILD-UP OF THE HEAP

!-----> LOOP OVER ALL HIERACHICAL LEVELS OF THE HEAP TREE

  Loop_L : do L = N/2 , 1 , -1
     ia = na(L)
     ib = nb(L)
     ic = nc(L)
!!$     p1 = (ia-imin+1) + (ib-jmin)*a + (ic-kmin)*a*b
     p1 = (ia-imin)*b*c + (ib-jmin)*c + (ic-kmin+1)
     I   = L
     II  = L + L

!-----> GO DOWN ALL THE HIERACHICAL LEVELS OF THE HEAP TREE

20   IF ( II <= N ) THEN
!!$        p2 = (na(ii)-imin+1) + (nb(ii)-jmin)*a + (nc(ii)-kmin)*a*b
        p2 = (na(ii)-imin)*b*c + (nb(ii)-jmin)*c + (nc(ii)-kmin+1)
!-----> COMPARE NEIGHBOURING ELEMENTS
        IF ( II < N ) THEN
!!$           p3 = (na(ii+1)-imin+1) + (nb(ii+1)-jmin)*a + (nc(ii+1)-kmin)*a*b
           p3 = (na(ii+1)-imin)*b*c+(nb(ii+1)-jmin)*c+(nc(ii+1)-kmin+1)
           if(p2 < p3) then
              ii = ii + 1
              p2 = p3
           endif
        endif

!-----> COMPARE THE ELEMENTS OF TWO HIRACHICAL LEVELS
!       PROMOTE LARGER ONE, DEMOTE SMALLER ONE

        if( p1 < p2 ) then
           na(i) = na(ii)
           nb(i) = nb(ii)
           nc(i) = nc(ii)
           I     = II
           II    = II + II
        ELSE

!-----> THIS PART OF THE TREE IS ORDERED , STOP BY PUTTING II=N+1
           II    = N + 1
        END IF
        GO TO 20
     END IF
!-----> PUT ELEMENTS IN THE PROPER SLOT
     na(i) = ia
     nb(i) = ib
     nc(i) = ic
  enddo Loop_L

!=====> NOW COLLECT ALL ELEMENTS FROM THE HEAP

  Loop_IHEAP :  do IHEAP = N , 1 , -1

     ia = na(iheap)
     ib = nb(iheap)
     ic = nc(iheap)
!!$     p1  = (ia-imin+1) + (ib-jmin)*a + (ic-kmin)*a*b
     p1  = (ia-imin)*b*c + (ib-jmin)*c + (ic-kmin+1)

!-----> THE FIRST ELEMENT IS ALWAYS THE LARGEST
     na(iheap) = na(1)
     nb(iheap) = nb(1)
     nc(iheap) = nc(1)
!-----> NOW LARGEST ELEMENT OF ALL BR(I) WITH 1<=I<=IHEAP IS STORED

     I   = 1
     II  = 2

!-----> NOW GENERATE LARGEST ELEMENT OF BR(I) WITH 1<=I<=IHEAP-1

40   IF ( II <= IHEAP - 1 ) THEN

!!$        p2 = (na(ii)-imin+1) + (nb(ii)-jmin)*a + (nc(ii)-kmin)*a*b
        p2 = (na(ii)-imin)*b*c + (nb(ii)-jmin)*c + (nc(ii)-kmin+1)
!-----> COMPARE NEIGHBOURING ELEMENTS
        if ( ii < iheap - 1 ) then
!!$           p3 = (na(ii+1)-imin+1) + (nb(ii+1)-jmin)*a + (nc(ii+1)-kmin)*a*b
           p3 = (na(ii+1)-imin)*b*c+(nb(ii+1)-jmin)*c+(nc(ii+1)-kmin+1)
           if ( p2 < p3 ) then
              ii = ii + 1
              p2 = p3
           endif
        endif

!-----> PROMOTE EACH ELEMENT OF THE TWIG OF BR UNTIL BRR > BR(I)
        if ( p1 < p2 ) then
           na(i) = na(ii)
           nb(i) = nb(ii)
           nc(i) = nc(ii)
           I     = II
           II    = II + II
        ELSE

!-----> THIS PART OF THE TREE IS PROMOTED , STOP BY PUTTING II=IHEAP+1
           II    = IHEAP + 1
        END IF
        GO TO 40
     END IF
!-----> PUT ELEMENTS IN THE PROPER SLOT
     na(i) = ia
     nb(i) = ib
     nc(i) = ic
  enddo Loop_IHEAP
end subroutine spsort

!C-----CALL HPSOR2(KG,IGP,GX,GY,GZ,GR)
!
subroutine gsort(ttr,k,n,nG)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: gsort
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

!C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!C     THIS SUBROUTINE SORTS THE ELEMENTS OF ARRAY BR IN ASCENDING
!C     ORDER. THE CODE IS BASED UPON HEAPSORT ALGORITHM WHICH HAS
!C     N*LOG(N) SCALING BEHAVIOUR, PARTICULARLY FAVOURABLE FOR LARGE
!C     VECTORS BR
!C                              STEFAN BL"UGEL, ISSP, NOV 1989
!c     interface is changed by T.Yamasaki at 28th Jun. 1994
!c     changed as it accords with fortran90 by T. Yamsaki at 17th Apr.1998
!c
!C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  use m_Const_Parameters, only : DP
  implicit none
  real(kind=DP), intent(in), dimension(6)   :: ttr
  integer, intent(in)                       :: k, n
  integer, intent(inout), dimension(k,3)    :: nG

!      integer DELTA
!      parameter (DELTA = 1.d-10)
  integer  l,i,ii,iheap
  integer ia,ib,ic
  real*8  brr,brr2, brr3

!C=====> BUILD-UP OF THE HEAP

!-----> LOOP OVER ALL HIERACHICAL LEVELS OF THE HEAP TREE

  loop_L :  do L = N/2 , 1 , -1
     ia = nG(L,1)
     ib = nG(L,2)
     ic = nG(L,3)
     brr    = ttr(1)*ia*ia + ttr(2)*ib*ib + ttr(3)*ic*ic &
          & + ttr(4)*ia*ib + ttr(5)*ib*ic + ttr(6)*ic*ia
     I   = L
     II  = L + L

!-----> GO DOWN ALL THE HIERACHICAL LEVELS OF THE HEAP TREE

20   IF ( II <= N ) THEN
        brr2 =  ttr(1)*nG(ii,1)*nG(ii,1)+ttr(2)*nG(ii,2)*nG(ii,2) &
             &+ ttr(3)*nG(ii,3)*nG(ii,3)+ttr(4)*nG(ii,1)*nG(ii,2) &
             &+ ttr(5)*nG(ii,2)*nG(ii,3)+ttr(6)*nG(ii,3)*nG(ii,1)

!-----> COMPARE NEIGHBOURING ELEMENTS
        IF ( II < N ) THEN
           brr3 = &
                &  ttr(1)*nG(ii+1,1)*nG(ii+1,1)+ttr(2)*nG(ii+1,2)*nG(ii+1,2) &
                &+ ttr(3)*nG(ii+1,3)*nG(ii+1,3)+ttr(4)*nG(ii+1,1)*nG(ii+1,2) &
                &+ ttr(5)*nG(ii+1,2)*nG(ii+1,3)+ttr(6)*nG(ii+1,3)*nG(ii+1,1)
           if(brr2 < brr3) then
!             if(brr2  < brr3 - DELTA) then
              ii = ii + 1
              brr2 = brr3
           endif
        endif

!-----> COMPARE THE ELEMENTS OF TWO HIRACHICAL LEVELS
!       PROMOTE LARGER ONE, DEMOTE SMALLER ONE

!           if( brr < brr2 - DELTA ) then
        if( brr < brr2 ) then
           nG(i,1) = nG(ii,1)
           nG(i,2) = nG(ii,2)
           nG(i,3) = nG(ii,3)
           I     = II
           II    = II + II
        ELSE

!-----> THIS PART OF THE TREE IS ORDERED , STOP BY PUTTING II=N+1
           II    = N + 1
        END IF
        GO TO 20
     END IF
!-----> PUT ELEMENTS IN THE PROPER SLOT
     nG(i,1) = ia
     nG(i,2) = ib
     nG(i,3) = ic
  enddo Loop_L

!=====> NOW COLLECT ALL ELEMENTS FROM THE HEAP

  loop_IHEAP :  do IHEAP = N , 1 , -1

     ia = nG(iheap,1)
     ib = nG(iheap,2)
     ic = nG(iheap,3)
     brr  =   ttr(1)*ia*ia + ttr(2)*ib*ib + ttr(3)*ic*ic &
          & + ttr(4)*ia*ib + ttr(5)*ib*ic + ttr(6)*ic*ia

!-----> THE FIRST ELEMENT IS ALWAYS THE LARGEST
     nG(iheap,1) = nG(1,1)
     nG(iheap,2) = nG(1,2)
     nG(iheap,3) = nG(1,3)
!-----> NOW LARGEST ELEMENT OF ALL BR(I) WITH 1<=I<=IHEAP IS STORED

     I   = 1
     II  = 2

!-----> NOW GENERATE LARGEST ELEMENT OF BR(I) WITH 1<=I<=IHEAP-1

40   IF ( II <= IHEAP - 1 ) THEN

        brr2 =    ttr(1)*ng(ii,1)*nG(ii,1) + ttr(2)*nG(ii,2)*nG(ii,2) &
             &  + ttr(3)*nG(ii,3)*nG(ii,3) + ttr(4)*nG(ii,1)*nG(ii,2) &
             &  + ttr(5)*nG(ii,2)*nG(ii,3) + ttr(6)*nG(ii,3)*nG(ii,1)
!-----> COMPARE NEIGHBOURING ELEMENTS
        if ( ii < iheap - 1 ) then
           brr3 = &
                &  ttr(1)*nG(ii+1,1)*nG(ii+1,1)+ttr(2)*nG(ii+1,2)*nG(ii+1,2)&
                &+ ttr(3)*nG(ii+1,3)*nG(ii+1,3)+ttr(4)*nG(ii+1,1)*nG(ii+1,2)&
                &+ ttr(5)*nG(ii+1,2)*nG(ii+1,3)+ttr(6)*nG(ii+1,3)*nG(ii+1,1)
!             if ( brr2 < brr3 - DELTA ) then
           if ( brr2 < brr3 ) then
              ii = ii + 1
              brr2 = brr3
           endif
        endif

!-----> PROMOTE EACH ELEMENT OF THE TWIG OF BR UNTIL BRR > BR(I)
!           if ( brr < brr2 - DELTA ) then
        if ( brr < brr2 ) then
           nG(i,1) = nG(ii,1)
           nG(i,2) = nG(ii,2)
           nG(i,3) = nG(ii,3)
           I     = II
           II    = II + II
        ELSE

!-----> THIS PART OF THE TREE IS PROMOTED , STOP BY PUTTING II=IHEAP+1
           II    = IHEAP + 1
        END IF
        GO TO 40
     END IF
!-----> PUT ELEMENTS IN THE PROPER SLOT
     nG(i,1) = ia
     nG(i,2) = ib
     nG(i,3) = ic
  enddo Loop_IHEAP
end subroutine gsort

subroutine gsort2(ttr,k,n,nG,gr)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: gsort2
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

!C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!C     THIS SUBROUTINE SORTS THE ELEMENTS OF ARRAY BR IN ASCENDING
!C     ORDER. THE CODE IS BASED UPON HEAPSORT ALGORITHM WHICH HAS
!C     N*LOG(N) SCALING BEHAVIOUR, PARTICULARLY FAVOURABLE FOR LARGE
!C     VECTORS BR
!C                              STEFAN BL"UGEL, ISSP, NOV 1989
!c     interface is changed by T.Yamasaki at 28th Jun. 1994
!c     changed as it accords with fortran90 by T. Yamsaki at 17th Apr.1998
!c
!C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  use m_Const_Parameters, only : DP
  implicit none
  real(kind=DP), intent(in), dimension(6)   :: ttr
  integer, intent(in)                       :: k, n
  integer, intent(inout), dimension(k,3)    :: nG
  real(kind=DP), intent(inout), dimension(k) :: gr

!      integer DELTA
!      parameter (DELTA = 1.d-10)
  integer  l,i,ii,iheap
  integer ia,ib,ic
  real*8  brr,brr2, brr3

!C=====> BUILD-UP OF THE HEAP

!-----> LOOP OVER ALL HIERACHICAL LEVELS OF THE HEAP TREE

  loop_L :  do L = N/2 , 1 , -1
     ia = nG(L,1)
     ib = nG(L,2)
     ic = nG(L,3)
!!$     brr    = ttr(1)*ia*ia + ttr(2)*ib*ib + ttr(3)*ic*ic &
!!$          & + ttr(4)*ia*ib + ttr(5)*ib*ic + ttr(6)*ic*ia
     brr    = gr(L)
     I   = L
     II  = L + L

!-----> GO DOWN ALL THE HIERACHICAL LEVELS OF THE HEAP TREE

20   IF ( II <= N ) THEN
!!$        brr2 =  ttr(1)*nG(ii,1)*nG(ii,1)+ttr(2)*nG(ii,2)*nG(ii,2) &
!!$             &+ ttr(3)*nG(ii,3)*nG(ii,3)+ttr(4)*nG(ii,1)*nG(ii,2) &
!!$             &+ ttr(5)*nG(ii,2)*nG(ii,3)+ttr(6)*nG(ii,3)*nG(ii,1)
        brr2 =  gr(ii)
!-----> COMPARE NEIGHBOURING ELEMENTS
        IF ( II < N ) THEN
!!$           brr3 = &
!!$                &  ttr(1)*nG(ii+1,1)*nG(ii+1,1)+ttr(2)*nG(ii+1,2)*nG(ii+1,2) &
!!$                &+ ttr(3)*nG(ii+1,3)*nG(ii+1,3)+ttr(4)*nG(ii+1,1)*nG(ii+1,2) &
!!$                &+ ttr(5)*nG(ii+1,2)*nG(ii+1,3)+ttr(6)*nG(ii+1,3)*nG(ii+1,1)
           brr3 =  gr(ii+1)
           if(brr2 < brr3) then
!             if(brr2  < brr3 - DELTA) then
              ii = ii + 1
              brr2 = brr3
           endif
        endif

!-----> COMPARE THE ELEMENTS OF TWO HIRACHICAL LEVELS
!       PROMOTE LARGER ONE, DEMOTE SMALLER ONE

!           if( brr < brr2 - DELTA ) then
        if( brr < brr2 ) then
           nG(i,1) = nG(ii,1)
           nG(i,2) = nG(ii,2)
           nG(i,3) = nG(ii,3)
           gr(i)   = gr(ii)
           I     = II
           II    = II + II
        ELSE

!-----> THIS PART OF THE TREE IS ORDERED , STOP BY PUTTING II=N+1
           II    = N + 1
        END IF
        GO TO 20
     END IF
!-----> PUT ELEMENTS IN THE PROPER SLOT
     nG(i,1) = ia
     nG(i,2) = ib
     nG(i,3) = ic
     gr(i)   = brr
  enddo Loop_L

!=====> NOW COLLECT ALL ELEMENTS FROM THE HEAP

  loop_IHEAP :  do IHEAP = N , 1 , -1

     ia = nG(iheap,1)
     ib = nG(iheap,2)
     ic = nG(iheap,3)
!!$     brr  =   ttr(1)*ia*ia + ttr(2)*ib*ib + ttr(3)*ic*ic &
!!$          & + ttr(4)*ia*ib + ttr(5)*ib*ic + ttr(6)*ic*ia
     brr  =   gr(iheap)
!-----> THE FIRST ELEMENT IS ALWAYS THE LARGEST
     nG(iheap,1) = nG(1,1)
     nG(iheap,2) = nG(1,2)
     nG(iheap,3) = nG(1,3)
     gr(iheap)   = gr(1)
!-----> NOW LARGEST ELEMENT OF ALL BR(I) WITH 1<=I<=IHEAP IS STORED

     I   = 1
     II  = 2

!-----> NOW GENERATE LARGEST ELEMENT OF BR(I) WITH 1<=I<=IHEAP-1

40   IF ( II <= IHEAP - 1 ) THEN

!!$        brr2 =    ttr(1)*ng(ii,1)*nG(ii,1) + ttr(2)*nG(ii,2)*nG(ii,2) &
!!$             &  + ttr(3)*nG(ii,3)*nG(ii,3) + ttr(4)*nG(ii,1)*nG(ii,2) &
!!$             &  + ttr(5)*nG(ii,2)*nG(ii,3) + ttr(6)*nG(ii,3)*nG(ii,1)
        brr2 =    gr(ii)
!-----> COMPARE NEIGHBOURING ELEMENTS
        if ( ii < iheap - 1 ) then
!!$           brr3 = &
!!$                &  ttr(1)*nG(ii+1,1)*nG(ii+1,1)+ttr(2)*nG(ii+1,2)*nG(ii+1,2)&
!!$                &+ ttr(3)*nG(ii+1,3)*nG(ii+1,3)+ttr(4)*nG(ii+1,1)*nG(ii+1,2)&
!!$                &+ ttr(5)*nG(ii+1,2)*nG(ii+1,3)+ttr(6)*nG(ii+1,3)*nG(ii+1,1)
           brr3 = gr(ii+1)
!             if ( brr2 < brr3 - DELTA ) then
           if ( brr2 < brr3 ) then
              ii = ii + 1
              brr2 = brr3
           endif
        endif

!-----> PROMOTE EACH ELEMENT OF THE TWIG OF BR UNTIL BRR > BR(I)
!           if ( brr < brr2 - DELTA ) then
        if ( brr < brr2 ) then
           nG(i,1) = nG(ii,1)
           nG(i,2) = nG(ii,2)
           nG(i,3) = nG(ii,3)
           gr(i)   = gr(ii)
           I     = II
           II    = II + II
        ELSE

!-----> THIS PART OF THE TREE IS PROMOTED , STOP BY PUTTING II=IHEAP+1
           II    = IHEAP + 1
        END IF
        GO TO 40
     END IF
!-----> PUT ELEMENTS IN THE PROPER SLOT
     nG(i,1) = ia
     nG(i,2) = ib
     nG(i,3) = ic
     gr(i)   = brr
  enddo Loop_IHEAP
end subroutine gsort2


