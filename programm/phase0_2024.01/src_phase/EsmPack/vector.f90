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

!-----------------------------------------------------------------
!
! Making reciprocal lattice vectors from lattice vectors
!
Subroutine MakeRecLatVec( A, B )
  Implicit None
  Real(8), intent(in ) :: A(3,3)
  Real(8), intent(out) :: B(3,3)
  Real(8), Parameter :: PI2 = 2.0d0 * 3.141592653589793238d0
  Real(8), external :: InnerProd

  Call OuterProd( A(1,2), A(1,3), B(1,1) )
  Call OuterProd( A(1,3), A(1,1), B(1,2) )
  Call OuterProd( A(1,1), A(1,2), B(1,3) )
  
  B(:,1) = B(:,1) / InnerProd( A(1,1), B(1,1) )
  B(:,2) = B(:,2) / InnerProd( A(1,2), B(1,2) )
  B(:,3) = B(:,3) / InnerProd( A(1,3), B(1,3) )

  Return
End Subroutine MakeRecLatVec

!-----------------------------------------------------------------
!
! Outer product of 3D vector
!
Subroutine OuterProd( A, B, C )
  Implicit None
  Real(8), Intent(in ) :: A(3), B(3)
  Real(8), Intent(out) :: C(3)

  C(1) = A(2) * B(3) - A(3) * B(2)
  C(2) = A(3) * B(1) - A(1) * B(3)
  C(3) = A(1) * B(2) - A(2) * B(1)

  Return
End Subroutine OuterProd

!-----------------------------------------------------------------
!
! Inner product of 3D vector
!
Function InnerProd( A, B )
  Implicit None
  Real(8), Intent(in ) :: A(3), B(3)
  Real(8) :: InnerProd

  InnerProd = Sum( A(:) * B(:) )

  Return
End Function InnerProd

