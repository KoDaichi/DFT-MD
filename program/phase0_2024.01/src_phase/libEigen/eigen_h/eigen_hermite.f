!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  n   : dimension of matrix
!  ar  : matrix(real part)
!  ai  : matrix(imaginary part)
!  w   : array for eigen value
!  lda : size of array "a"
!  m0  : coefficient of blocking
!  ifl : switch flag ( 0:eigen-value and eigen-vector )
!                    ( 1:eigen-value only )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine eigen_h(n,ar,ai,w,lda,m0,ifl)
!
      use communication_h, only : eigen_init
     &               , eigen_free
!
      implicit none
      integer, intent(in) :: n, lda, m0, ifl
      real(8), intent(inout) :: ar(lda,*), ai(lda,*), w(*)
      real(8), pointer :: zr(:), zi(:), d(:), e(:), e2(:), tau(:)
      integer :: nb, nm, nma, nmz, nmz1, nmw, nmw1, info, ierr
      integer :: nprow, npcol, i_1, istat
      real(8) :: hs0, hs1
!
      include 'mpif.h'
      include 'trd.h'
!
      call eigen_init(2)
      NPROW = size_of_col
      NPCOL = size_of_row
      call eigen_free(0)
!
      NB  = 48
      nmz = ((n-1)/NPROW+1)
      nmz = ((nmz-1)/NB+1)*NB+1
      nmw = ((n-1)/NPCOL+1)
      nmw = ((nmw-1)/NB+1)*NB+1

      nm=(n/2)*2+1

      nmz1=nm
      nmw1=((n-1)/nprocs+1)

      allocate(
     &               zr(MAX(nmz*nmw+nm+1,nmz1*nmw1)),
     &               zi(MAX(nmz*nmw+nm+1,nmz1*nmw1)),
     &               d(nm),e(nm),e2(nm),tau(2*nm),
     &               stat=istat)
      if(istat.ne.0) then
           if(myrank==1) print*,"Memory exhausted"
           call mpi_abort(mpi_comm_eigen, 1, ierr)
      endif
*
         hs0=MPI_Wtime()
         nma=lda
!----
         call eigen_hrd(n,ar,ai,nma,d,e,e2,tau,nm,m0,zr,zi)
!----
         e2(1:n)=e(1:n)
         w(1:n)=d(1:n)
!----
         call eigen_dch(n,w,e2(2),zr,nmz1,info)
!----
         nma=lda
         call eigen_hbk(ar,ai,nma, zr,zi,nmz1, tau,nm,  n)
!----
#ifdef TIMER
         hs1=MPI_Wtime()
         if(myrank==1)then
            print*," "
            print*,"Total Execution Time of eigen_h =",hs1-hs0,"(sec)"
         endif
#endif
*
         do i_1=1,(n-1)/nprocs+1
           ar(1:n,i_1)=zr(1+nmz1*(i_1-1):n+nmz1*(i_1-1))
           ai(1:n,i_1)=zi(1+nmz1*(i_1-1):n+nmz1*(i_1-1))
         enddo

         deallocate(zr,zi,d,e,e2,tau)
*
      return
      end

