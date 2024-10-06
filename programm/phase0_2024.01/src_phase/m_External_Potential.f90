
module m_External_Potential
  use m_Control_Parameters, only : sw_external_potential
  use m_Parallelization, only : MPI_CommGroup

  implicit none

  ! external point charge
  integer :: natom_ex = 0
  real(8), allocatable, dimension(:,:) :: pos_ex
  real(8), allocatable, dimension(:) :: charge_ex
  real(8), allocatable, dimension(:) :: coulomb_rc
  real(8), allocatable, dimension(:) :: coulomb_rn
  ! external potential (electrostatic potential)
  real(8), allocatable, dimension(:,:) :: espot
  real(8), allocatable, dimension(:,:,:) :: espot_g

  ! temporary
  real(8), dimension(:,:), pointer :: mesh_coord
  real(8) :: mesh_delta(3,3)
  real(8), dimension(:), pointer :: charge_r

contains

  subroutine m_EP_alloc_point_charge
    implicit none

    if( .not. allocated(pos_ex) ) allocate(pos_ex(natom_ex,3))
    if( .not. allocated(charge_ex) ) allocate(charge_ex(natom_ex))
    if( .not. allocated(coulomb_rc) ) allocate(coulomb_rc(natom_ex))
    if( .not. allocated(coulomb_rn) ) allocate(coulomb_rn(natom_ex))

  end subroutine m_EP_alloc_point_charge

  subroutine m_EP_dealloc_point_charge
    implicit none

    if( allocated(pos_ex) ) deallocate(pos_ex)
    if( allocated(charge_ex) ) deallocate(charge_ex)
    if( allocated(coulomb_rc) ) deallocate(coulomb_rc)
    if( allocated(coulomb_rn) ) deallocate(coulomb_rn)

  end subroutine m_EP_dealloc_point_charge

  subroutine m_EP_set_point_charge(n,pos,q,rc,rn)
    implicit none
    integer n
    real(8) pos(n,3)
    real(8) q(n)
    real(8) rc(n)
    real(8) rn(n)
    integer i

    natom_ex = n
    call m_EP_alloc_point_charge

    pos_ex = pos
    charge_ex = q
    coulomb_rc = rc
    coulomb_rn = rn

  end subroutine m_EP_set_point_charge

! === "phase" repository merge: To make 3D_Parallel by tkato 2011/11/10 ========
! ==============================================================================
  subroutine m_EP_espot_by_external_charge

    use m_PlaneWaveBasisSet,    only : igfp_l
    use m_Crystal_Structure,    only : altv, univol
    use m_Control_Parameters,   only : nspin, kimg
    use m_Charge_Density,       only : chgq_l
    use m_FFT, only : nfftp, fft_box_size_CD, &
                      m_FFT_CD_direct, m_FFT_alloc_CD_box, m_FFT_dealloc_CD_box
    use m_Parallelization,      only : ista_kngp,iend_kngp, npes, mype, ierr &
                                     , ista_fftp,iend_fftp,ista_fftph,iend_fftph &
                                     , nis_fftp, nie_fftp, nel_fftp, idisp_fftp, mp_fftp
    use m_Files, only : nfout
    use mpi

    implicit none
!    include 'mpif.h'

    integer i, j, k, l, n, ik, iloop
    integer idp, nlp, nmp, nnp, nd2p
    integer ispin
    integer nn, ip
    integer, allocatable, dimension(:)  :: inx, jnx, knx
    real(8) rinplw
    real(8) r(3), r0(3), dist, rc, rn
    real(8) vpot, vpot1, vpot2
    real(8), allocatable, dimension(:) :: afft
    real(8), allocatable, dimension(:) :: afft_mpi1

    if(mype==0) write(6,*) 'electrostatic potential (point charge)...'

    ! allocate
    call m_FFT_alloc_CD_box()

    if( .not. allocated( espot ) ) then
       allocate( espot(ista_fftph:iend_fftph,nspin ) )
       espot = 0.0d0
    end if
    if( .not. allocated( espot_g ) ) then
       allocate( espot_g(ista_kngp:iend_kngp,kimg,nspin))
       espot_g = 0.d0
    end if

    allocate(afft(ista_fftp:iend_fftp))
    allocate(afft_mpi1(nfftp))
    !allocate(inx(ista_fftp:iend_fftp))
    !allocate(jnx(ista_fftp:iend_fftp))
    !allocate(knx(ista_fftp:iend_fftp))
    allocate(inx(1:nfftp))
    allocate(jnx(1:nfftp))
    allocate(knx(1:nfftp))

    ! map FFT mesh

    idp  = fft_box_size_CD(1,0)
    nlp  = fft_box_size_CD(1,1)
    nmp  = fft_box_size_CD(2,1)
    nnp  = fft_box_size_CD(3,1)
    nd2p = fft_box_size_CD(2,0)

    !do n = ista_fftp, iend_fftp
    do n = 1, nfftp
       nn = (n+kimg-1)/kimg
       i  = mod(nn,idp)
       j  = mod((nn-1)/idp,nd2p) + 1
       k  = (nn - (j-1)*idp - i)/(idp*nd2p) + 1
       inx(n) = i - 1;
       jnx(n) = j - 1;
       knx(n) = k - 1;
    end do

    espot = 0.d0
    do i = 1, natom_ex

       write(6,'(i5,3f15.8,3f10.5)') i, pos_ex(i,:), charge_ex(i), coulomb_rc(i), coulomb_rn(i)

       do k = ista_fftph, iend_fftph
          ip = 2*k-1
          r(1:3) = altv(1:3,1)*dble(inx(ip))/dble(nlp)+altv(1:3,2)*dble(jnx(ip))/dble(nmp)&
               +altv(1:3,3)*knx(ip)/dble(nnp)

          r0(1:3) = pos_ex(i,1:3)
          dist = sqrt(dot_product(r-r0, r-r0))
          rc = coulomb_rc(i)
          rn = coulomb_rn(i)
          vpot1 = 0.d0
          vpot2 = 0.d0
          if( dist > 0.d0 ) then
             do n = 0, int(rn)
                if( n /= 0 ) vpot1 = vpot1 + (rc**dble(rn-n))*(dist**dble(n-1))
                vpot2 = vpot2 + (rc**dble(rn-n))*(dist**dble(n))
             end do
          else
             vpot1 = rc**(rn-1.d0)
             vpot2 = rc**rn
          end if
          vpot = vpot1 / vpot2
          espot(k,1:nspin)=espot(k,1:nspin)-charge_ex(i)*vpot
       end do
    end do

    afft_mpi1 = 0.d0
    do iloop = 1, nspin
       ! copy R-space potential to FFT mesh
       afft = 0.d0
       do i = ista_fftph, iend_fftph
          afft(2*i-1) = espot(i,iloop)
       end do

       ! FFTl to G-space
       call m_FFT_CD_direct(nfout,afft)

       if(npes >= 2) then
          call mpi_allgatherv(afft,nel_fftp(mype),mpi_double_precision &  ! MPI
               &       ,afft_mpi1,nel_fftp,idisp_fftp,mpi_double_precision,MPI_CommGroup,ierr)
       else
          afft_mpi1 = afft
       end if

       rinplw = 1.d0/product(fft_box_size_CD(1:3,1))
       do ik = 1, kimg
          do i = ista_kngp, iend_kngp  !for mpi
             ip = (igfp_l(i)-1)*kimg + ik
             espot_g(i,ik,iloop)   = afft_mpi1(ip)*rinplw
          end do
       end do
    end do

    ! deallocate
    deallocate( inx )
    deallocate( jnx )
    deallocate( knx )
    deallocate(afft)
    deallocate(afft_mpi1)
    call m_FFT_dealloc_CD_box()

  end subroutine m_EP_espot_by_external_charge
! === "phase" repository merge: To make 3D_Parallel by tkato 2011/11/10 ========
! ==============================================================================

end module m_External_Potential

! === "phase" repository merge: To make 3D_Parallel by tkato 2011/11/10 ========
! ==============================================================================
subroutine charge_density_realspace

  use m_PlaneWaveBasisSet,    only : igfp_l
  use m_Crystal_Structure,    only : altv, univol
  use m_Control_Parameters,   only : nspin, kimg
  use m_Charge_Density,       only : chgq_l
  use m_FFT, only : nfftp, fft_box_size_CD, &
                    m_FFT_CD_direct, m_FFT_CD_inverse_c, m_FFT_check_of_negative_CD, &
                    m_FFT_alloc_CD_box, m_FFT_dealloc_CD_box
  use m_Parallelization,      only : ista_kngp,iend_kngp, npes, mype, ierr &
                                   , npes_cdfft &
                                   , ista_fftp,iend_fftp,ista_fftph,iend_fftph &
                                   , nis_fftp, nie_fftp, nel_fftp, idisp_fftp, mp_fftp &
                                   , MPI_CommGroup
  use m_Files, only : nfout
  use m_External_Potential, only : mesh_coord, mesh_delta, charge_r
  use mpi

  implicit none

!  include 'mpif.h'

  integer i, j, k, l, n, ik, iloop
  integer idp, nlp, nmp, nnp, nd2p
  integer ispin
  integer nn, ip
  integer, allocatable, dimension(:)  :: inx, jnx, knx
  real(8) r(3)
  real(8), allocatable, dimension(:) :: afft
  real(8), allocatable, dimension(:) :: afft_mpi1
  real(8), allocatable, dimension(:) :: afft_mpi2,afft_mpi3
  real(8), allocatable, dimension(:) :: work

  ! allocate
  call m_FFT_alloc_CD_box()

  !ista_fftph = 1
  !iend_fftph = nfftp/2
  if( associated( mesh_coord ) ) deallocate(mesh_coord)
  allocate( mesh_coord(1:nfftp,3) )
  if(associated(charge_r)) deallocate(charge_r)
  allocate(charge_r(1:nfftp/2))
  charge_r = 0.0d0

  allocate(afft(ista_fftp:iend_fftp))
  allocate(afft_mpi1(nfftp))
  if(npes >= 2) then
     allocate(afft_mpi2(mp_fftp))
     allocate(afft_mpi3(mp_fftp))
  end if
  !  allocate(inx(ista_fftp:iend_fftp))
  !  allocate(jnx(ista_fftp:iend_fftp))
  !  allocate(knx(ista_fftp:iend_fftp))
  allocate(inx(1:nfftp))
  allocate(jnx(1:nfftp))
  allocate(knx(1:nfftp))

  ! map FFT mesh

  idp  = fft_box_size_CD(1,0)
  nlp  = fft_box_size_CD(1,1)
  nmp  = fft_box_size_CD(2,1)
  nnp  = fft_box_size_CD(3,1)
  nd2p = fft_box_size_CD(2,0)

  !  do n = ista_fftp, iend_fftp
  do n = 1, nfftp
     nn = (n+kimg-1)/kimg
     i  = mod(nn,idp)
     j  = mod((nn-1)/idp,nd2p) + 1
     k  = (nn - (j-1)*idp - i)/(idp*nd2p) + 1
     inx(n) = i - 1;
     jnx(n) = j - 1;
     knx(n) = k - 1;
  end do

  do k = 1, nfftp/2
     ip = 2*k-1
     r(1:3) = altv(1:3,1)*dble(inx(ip))/dble(nlp)+altv(1:3,2)*dble(jnx(ip))/dble(nmp)&
          +altv(1:3,3)*knx(ip)/dble(nnp)
     mesh_coord(ip,1:3) = r(1:3)
     mesh_delta(1:3,1) = altv(1:3,1)/dble(nlp)
     mesh_delta(1:3,2) = altv(1:3,2)/dble(nmp)
     mesh_delta(1:3,3) = altv(1:3,3)/dble(nnp)
  end do

  afft_mpi1 = 0.d0
  do iloop = 1, nspin
     ! map charge on FFT box
     do j = 1, kimg
        do i = ista_kngp, iend_kngp  !for mpi
           ip = (igfp_l(i)-1)*kimg + j
           afft_mpi1(ip) = afft_mpi1(ip) + chgq_l(i,j,iloop) !mpi
        end do
     end do

     if(npes >= 2) then
        call mpi_barrier(MPI_CommGroup,ierr)
        do j = 0, npes-1

           do i = nis_fftp(j),nie_fftp(j)
              afft_mpi2(i-nis_fftp(j)+1) = afft_mpi1(i)
           end do
           call mpi_allreduce(afft_mpi2,afft_mpi3,mp_fftp &
                &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)

           if(j == mype) then
              do i = ista_fftp, iend_fftp
                 afft(i) = afft_mpi3(i - ista_fftp + 1)
              end do
           end if
        end do
     else
        afft = afft_mpi1
     end if

     ! FFT G-space to R-space
     call m_FFT_CD_inverse_c(nfout,afft)
     !call m_FFT_check_of_negative_CD(afft,nfout,nspin,iloop)
     call m_FFT_check_of_negative_CD(npes_cdfft,ista_fftp,iend_fftp&
          & ,ista_fftph,iend_fftph,afft,nfout,nspin,iloop)

     ! copy afft to charge density(in R-space)
     do i = ista_fftph, iend_fftph
        charge_r(i) = afft(i*2-1)*univol/dble(idp*nd2p*nd2p)
     end do
  end do

  allocate(work(1:nfftp/2))
  call mpi_allreduce(charge_r,work,nfftp/2 &
       &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
  charge_r = work
  deallocate(work)

  ! deallocate
  deallocate( inx )
  deallocate( jnx )
  deallocate( knx )
  deallocate(afft)
  deallocate(afft_mpi1)
  if(npes >= 2) then
     deallocate(afft_mpi2)
     deallocate(afft_mpi3)
  end if
  call m_FFT_dealloc_CD_box()

end subroutine charge_density_realspace
! === "phase" repository merge: To make 3D_Parallel by tkato 2011/11/10 ========
! ==============================================================================


