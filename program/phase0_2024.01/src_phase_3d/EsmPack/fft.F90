! Copyright (c) 2012, Minoru Otani <minoru.otani@aist.go.jp> 
! 
! Permission is hereby granted, free of charge, to any person 
! obtaining a copy of this software and associated documentation 
! files (the "Software"), to deal in the Software without restriction, 
! including without limitation the rights to use, copy, modify, merge, 
! publish, distribute, sublicense, and/or sell copies of the Software, 
! and to permit persons to whom the Software is furnished to do so, 
! subject to the following conditions:
 
! The above copyright notice and this permission notice shall be 
! included in all copies or substantial portions of the Software.
 
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
! OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
! NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
! HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
! WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
! DEALINGS IN THE SOFTWARE.

Subroutine cft_1z(c, nsl, nz, ldz, isign, cout)
  Implicit None
#ifndef SX
#ifdef __NEC__
  Include "aslfftw3.f"
#else
#include "fftw3.f"
#endif
#endif

  Integer,    Intent(In)  :: nsl, nz, ldz, isign
  Complex(8), Intent(In)  :: c(0:nsl*ldz-1)
  Complex(8), Intent(Out) :: cout(0:nsl*ldz-1)
  
#ifndef SX
  Integer(8) :: plan
  Integer :: FFTW_WARD, isl

  If ( isign < 0 ) Then
     FFTW_WARD = FFTW_FORWARD
  Else
     FFTW_WARD = FFTW_BACKWARD
  End If

  Do isl = 0, nsl-1
     Call dfftw_plan_dft_1d ( plan, nz, c(isl*ldz), cout(isl*ldz), FFTW_WARD, &
          FFTW_ESTIMATE )
     Call dfftw_execute ( plan )
     Call dfftw_destroy_plan ( plan )
  End Do

  If( isign < 0 ) Then
     cout(:) = cout(:) / nz
  End If
#else
  Integer, save :: isw=0
  Integer :: n,ierr,isl
  Real(8),allocatable,dimension(:) :: rpart,ipart
  Integer,save :: ifax(20)
  Real(8),allocatable,dimension(:),save :: trigs
  Real(8),allocatable,dimension(:) :: wk
  allocate(rpart(0:nsl*ldz-1));rpart(:) = dble(c(:))
  allocate(ipart(0:nsl*ldz-1));ipart(:) = dimag(c(:))
  n = nsl*ldz
  if(isw==0) allocate(trigs(2*n))
  allocate(wk(2*n))

  if(isw.eq.0)then
     do isl=0,nsl-1
        call dfc1fb(nz,rpart(isl*ldz),ipart(isl*ldz),nz,isw,ifax,trigs,wk,ierr)
     enddo
     isw=1
  else
     if(isign<0)then
        isw = 1
     else
        isw = -1
     endif
     do isl=0,nsl-1
        call dfc1bf(nz,rpart(isl*ldz),ipart(isl*ldz),nz,isw,ifax,trigs,wk,ierr)
     enddo
  endif
  cout(:) = dcmplx(rpart(:),ipart(:))
  if(isign<0) cout(:) = cout(:)/nz
  deallocate(rpart)
  deallocate(ipart)
  deallocate(wk)
#endif


End Subroutine cft_1z

