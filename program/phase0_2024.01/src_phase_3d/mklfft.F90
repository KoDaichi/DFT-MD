!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: mklfft
!
!  AUTHOR(S): T. Yamamoto   March/8/2006
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
#ifdef MKLFFT
module mklfft
!  $Id: mklfft.F90 570 2017-04-21 20:34:50Z yamasaki $
  use MKL_DFTI
  implicit none

  type mklfft_rdh
    TYPE(DFTI_DESCRIPTOR), POINTER :: dh_r1
    TYPE(DFTI_DESCRIPTOR), POINTER :: dh_c2
    integer :: worksize
    integer :: ldx,ldy,ldz
  end type

contains

  subroutine init_mkl_rfft(nx,ny,nz,ldx,ldy,ldz,rdh)
    integer, intent(in) :: nx,ny,nz,ldx,ldy,ldz
    type(mklfft_rdh), intent(inout) :: rdh

    integer :: strides(3),Status
    integer :: nc2(2)

    rdh%worksize = ldx*ldy*ldz/2
    rdh%ldx = ldx
    rdh%ldy = ldy
    rdh%ldz = ldz
    nc2(1)=nz; nc2(2)=ny

    ! REAL X FFT
    Status = DftiCreateDescriptor(rdh%dh_r1, &
     &       DFTI_DOUBLE,DFTI_REAL,1,nx)
    Status = DftiSetValue(rdh%dh_r1, &
     &       DFTI_NUMBER_OF_TRANSFORMS,ldy*nz)
    Status = DftiSetValue(rdh%dh_r1, &
     &       DFTI_INPUT_DISTANCE,ldx)
    Status = DftiCommitDescriptor(rdh%dh_r1)
    ! COMPLEX YZ FFT
    Status = DftiCreateDescriptor(rdh%dh_c2, &
     &       DFTI_DOUBLE,DFTI_COMPLEX,2,nc2)
    Status = DftiSetValue(rdh%dh_c2, &
    !! &       DFTI_NUMBER_OF_TRANSFORMS,nx)
     &       DFTI_NUMBER_OF_TRANSFORMS,nx/2+1)
    strides(1) = 0
    strides(2) = 1
    strides(3) = ldz
    Status = DftiSetValue(rdh%dh_c2,DFTI_INPUT_STRIDES,strides)
    Status = DftiSetValue(rdh%dh_c2,DFTI_INPUT_DISTANCE,ldz*ldy)
    Status = DftiCommitDescriptor(rdh%dh_c2)
  end subroutine init_mkl_rfft

  subroutine mkl_rfft(isw,rdh,a)
    integer, intent(in) :: isw
    type(mklfft_rdh), intent(in) :: rdh
    real(8), intent(inout) :: a(*)

    complex(8) :: work(rdh%worksize)
    integer :: ldx2,ldy,ldz

    ldx2 = rdh%ldx/2
    ldy = rdh%ldy
    ldz = rdh%ldz

    if(isw == -1) then
       call compute_forward(a,work)
    else if(isw == 1) then
       call compute_backward(a,work)
    else
       call phase_error_with_msg(nfout,'mklfft: isw is not 1 nor -1.',__LINE__,__FILE__)
    end if
  contains
    subroutine compute_forward(ar,ac)
      real(8), intent(inout) :: ar(*)
      complex(8), intent(inout) :: ac(*)
      integer :: Status
      Status = DftiComputeForward(rdh%dh_r1,ar)
      call trans_rc(ar,ac)
      Status = DftiComputeForward(rdh%dh_c2,ac)
      call trans_cr(ac,ar)
    end subroutine compute_forward
    subroutine compute_backward(ar,ac)
      real(8), intent(inout) :: ar(*)
      complex(8), intent(inout) :: ac(*)
      integer :: Status
      call trans_rc(ar,ac)
      Status = DftiComputeBackward(rdh%dh_c2,ac)
      call trans_cr(ac,ar)
      Status = DftiComputeBackward(rdh%dh_r1,ar)
    end subroutine compute_backward
    subroutine trans_rc(r,c)
      real(8), intent(in)  :: r(*)
      complex(8), intent(out) :: c(*)

      integer :: i,j,k
      integer :: ir,ic

      do i=1,ldx2
      do j=1,ldy
      do k=1,ldz
         ir = i + ldx2*((j-1) + ldy*(k-1))
         ic = k + ldz*((j-1) + ldy*(i-1))
         c(ic) = dcmplx(r(2*ir-1),r(2*ir))
      end do
      end do
      end do
    end subroutine trans_rc
    subroutine trans_cr(c,r)
      complex(8), intent(in) :: c(*)
      real(8), intent(out)  :: r(*)

      integer :: i,j,k
      integer :: ir,ic

      do i=1,ldx2
      do j=1,ldy
      do k=1,ldz
         ir = i + ldx2*((j-1) + ldy*(k-1))
         ic = k + ldz*((j-1) + ldy*(i-1))
         r(2*ir-1) = dble(c(ic))
         r(2*ir) = dimag(c(ic))
      end do
      end do
      end do
    end subroutine trans_cr
  end subroutine mkl_rfft

end module mklfft
#else
subroutine mklfft_dummy
end 
#endif
