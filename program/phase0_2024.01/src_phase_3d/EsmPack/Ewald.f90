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
!-------------------------------------------------------------------

!
! Calculates Ewald energy R-space terms.
!
Subroutine Ewald( nAtom, tau, charge, at, alat, alpha, bg, ewaldr )
  Implicit None
  Integer, Intent(in) :: nAtom            ! Number of Atoms
  Real(8), Intent(in) :: tau(3,nAtom)     ! Coordinates of Atoms
  Real(8), Intent(in) :: charge(nAtom)    ! Charge of Atoms
  Real(8), Intent(in) :: at(3,3)         ! Lattice Vectors
  Real(8), Intent(in) :: alat             ! Lattice Parameter
  Real(8), Intent(in) :: alpha            ! Convergence Factor
  Real(8), Intent(in) :: bg(3,3)          ! Reciprocal Lattice Vectors

  Real(8), Intent(out) :: ewaldr          ! Real Space Ewald Energy

  Real(8), external :: Dnrm2

  ! Local Variables
  Real(8) :: rmax, rmax2, dtau(3), t(3), t2, r
  Integer :: iA, iB, i, j, k, ni, nj, nk
  rmax = 4.0d0 / sqrt (alpha) / alat
  ni = int( Dnrm2(3, bg(1,1), 1) * rmax ) + 2
  nj = int( Dnrm2(3, bg(1,2), 1) * rmax ) + 2
  nk = int( Dnrm2(3, bg(1,3), 1) * rmax ) + 2

  rmax2 = rmax * rmax
  ewaldr = 0.0d0
  Do iA = 1, nAtom
  Do iB = 1, nAtom
     dtau(:) = tau(:,iA) - tau(:,iB)
     Do i = -ni, ni
     Do j = -nj, nj
     Do k = -nk, nk
        t(:) = dtau(:) + i * at(:,1) + j * at(:,2) + k * at(:,3)
        t2 = Sum( t(:) * t(:) )
        If( t2 > rmax2 ) Cycle
        If( t2 < 1.0d-10 ) Cycle

        r = Sqrt(t2) * alat
        ewaldr = ewaldr + charge(iA) * charge(iB) &
             * Erfc( Sqrt(alpha) * r ) / r
        
     End Do
     End Do
     End Do
  End Do
  End Do

  Return
End Subroutine Ewald

