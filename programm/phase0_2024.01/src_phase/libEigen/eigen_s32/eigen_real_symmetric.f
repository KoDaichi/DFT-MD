!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  n   : dimension of matrix
!  a   : matrix
!  lda : size of array "a"
!  w   : array for eigen value
!  z   : array for eigen vector
!  ldz : size of array "z"
!  d   : array to store a main diagonal element
!  e   : array to store a sub diagonal element
!  m0  : coefficient of blocking
!  ifl : switch flag ( 0:eigen-value and eigen-vector )
!                    ( 1:eigen-value only )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine eigen_real_symmetric(n, a, lda, w, z, ldz, m0, ifl)
      implicit none

      integer, intent(in)    :: n, lda, ldz, m0, ifl
      real(8), intent(inout) :: a(lda,*), w(*), z(ldz,*)

      include 'switch_ndim.h'

      if(n<ndim_switch) then
        call eigen_s(n, a, lda, w, z, ldz, m0, ifl)
      else
        call eigen_sx(n, a, lda, w, z, ldz, m0, ifl)
      endif

      end subroutine

