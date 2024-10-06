program wan_iterp
  implicit none
  integer, parameter :: DP = 8
  real(kind=DP), parameter :: PAI  = 3.141592653589793238462D0
  real(kind=DP), parameter :: PAI2 = PAI*2.d0, PAI4 = PAI*4.d0
  real(kind=DP), parameter :: Hartree = 27.21139615d0

  integer :: num_wan, nrpts, num_band_kpt
  integer, allocatable :: ndegen(:), ipos(:,:)
  real(kind=DP) :: efermi
  real(kind=DP), allocatable :: kpt_band(:,:)
  complex(kind=DP), allocatable :: hmat(:,:,:)
  character(len=128) :: seedname

  if(iargc()<1) then
     stop 'wan_interp.x seedname'
  end if
  call getarg(1,seedname)

!  call read_file_nnkp
  call read_file_hr
  call read_file_kpoint
  call read_file_efermi
  call calc_energy_band

contains

  subroutine read_file_nnkp
    integer, parameter :: nfwannier = 10
    integer :: i, nkpt
    real(kind=DP) :: altv(3,3), rltv(3,3)
    real(kind=DP), allocatable :: kvec(:,:)

    open(nfwannier,file=trim(seedname)//".nnkp")
    do i=1,5
       read(nfwannier,*)
    end do
    read(nfwannier,*) altv(1:3,1)
    read(nfwannier,*) altv(1:3,2)
    read(nfwannier,*) altv(1:3,3)
    read(nfwannier,*);     read(nfwannier,*);     read(nfwannier,*)
    read(nfwannier,*) rltv(1:3,1)
    read(nfwannier,*) rltv(1:3,2)
    read(nfwannier,*) rltv(1:3,3)
    read(nfwannier,*);     read(nfwannier,*);     read(nfwannier,*)
    read(nfwannier,*) nkpt

    allocate(kvec(3,nkpt))
    do i=1,nkpt
       read(nfwannier,*) kvec(1:3,i)
    end do
    close(nfwannier)
    
    write(*,'("Lattice vectors:")')
    do i=1,3
       write(*,'(3(1x,f10.5))') altv(1:3,i)
    end do
    write(*,'("Reciprocal lattice vectors:")')
    do i=1,3
       write(*,'(3(1x,f10.5))') rltv(1:3,i)
    end do
    write(*,'("k-points: nkpt=",i5)') nkpt
    do i=1,nkpt
       write(*,'(3(1x,f10.5))') kvec(1:3,i)
    end do
!
    deallocate( kvec )

  end subroutine read_file_nnkp

  subroutine read_file_hr
    integer, parameter :: nfwannier = 10

    integer :: ir, j, i, id, jd
    real(kind=DP) :: hr, hi

    open(nfwannier,file=trim(seedname)//"_hr.dat")
    read(nfwannier,*)
    read(nfwannier,*) num_wan
    read(nfwannier,*) nrpts
    allocate(ndegen(nrpts))
    allocate(hmat(num_wan,num_wan,nrpts))
    allocate(ipos(3,nrpts))
    read(nfwannier,*) ndegen(1:nrpts)
    do ir=1,nrpts 
       do j=1,num_wan
          do i=1,num_wan
             read(nfwannier,*) ipos(1:3,ir),id,jd, hr, hi 
             hmat(i,j,ir) = cmplx(hr,hi)
          end do
       end do
    end do
    close(nfwannier)

    write(*,'("Hamiltonian matrix in the WF basis:")')
    do ir=1,nrpts 
       do j=1,num_wan
          do i=1,num_wan
             write(*,'(5i5,2f10.6)') ipos(1:3,ir),i,j,hmat(i,j,ir)
          end do
       end do
    end do
  end subroutine read_file_hr

  subroutine read_file_kpoint
    integer                    :: fn_number_of_words
    integer,     parameter     :: len_str = 80
    integer, parameter :: nfkpoint = 10

    integer :: ikx,iky,ikz,idk, nwei, i
    character(len=len_str)     :: str

    open(nfkpoint,file="kpoint.data")
    read(nfkpoint,'(a80)') str

    select case (fn_number_of_words(str))
    case (0)
       stop ' error in reading of the nfkpoint file'
    case (1)
       read(str,*) num_band_kpt
       nwei = num_band_kpt

    case (2)
       read(str,*) num_band_kpt, nwei
    end select

    write(*,*) 'num_band_kpt =', num_band_kpt

    allocate(kpt_band(3,num_band_kpt));  kpt_band = 0.0d0

    select case (fn_number_of_words(str))
    case (1)
       do i=1,num_band_kpt
          read(nfkpoint,*) kpt_band(1:3,i)
       end do
    case (2)
       do i=1,num_band_kpt
          read(nfkpoint,*) ikx,iky,ikz,idk
          kpt_band(1,i) = real(ikx)/real(idk)
          kpt_band(2,i) = real(iky)/real(idk)
          kpt_band(3,i) = real(ikz)/real(idk)
       end do
    end select

    close( nfkpoint )
  end subroutine read_file_kpoint
!
  subroutine read_file_efermi
    integer, parameter :: nfefermi = 10

    efermi = 0.d0
    open(nfefermi,file="nfefermi.data")
    read(nfefermi,*) efermi
    close(nfefermi)
  end subroutine read_file_efermi

  subroutine calc_energy_band
    integer, parameter :: nfeng = 10

    integer :: i, nspin
    real(kind=DP), allocatable :: eig(:)
    complex(kind=DP), allocatable :: umat(:,:)

    nspin = 1
    open(nfeng,file="nfenergy.data")
    write(nfeng,'(" num_kpoints = ",i6)') num_band_kpt
    write(nfeng,'(" num_bands   = ",i6)') num_wan
    write(nfeng,'(" nspin       = ",i6)') nspin
    write(nfeng,'(" Fermi energy level = ",f10.5)') efermi
    write(nfeng,*)
    
    allocate(eig(num_wan))
    allocate(umat(num_wan,num_wan))
    do i=1,num_band_kpt
       call solve_hmat_k(kpt_band(1,i),eig(1),umat(1,1))
       write(nfeng,'("=== energy_eigen_values ===")')
       write(nfeng,'(" ik = ",i6," (",3f10.6," )")') i,kpt_band(1:3,i)
       write(nfeng,'(4f18.10)') eig(1:num_wan)/Hartree
    end do
    close(nfeng)

    deallocate( hmat, umat, eig )
    
  end subroutine calc_energy_band

  subroutine solve_hmat_k(kpt,eig,umat)
    implicit none
    real(kind=DP), intent(in) :: kpt(3)
    real(kind=DP), intent(out) :: eig(num_wan)
    complex(kind=DP), intent(out) :: umat(num_wan,num_wan)
 
    integer :: i,j,ir
    real(kind=DP) :: ph
    complex(kind=DP) :: zf
    complex(kind=DP), allocatable :: hmat_k(:,:)

    character(len=1) :: JOBZ,UPLO
    integer :: lwork, lrwork, info
    complex(kind=DP),allocatable,dimension(:) :: work
    real(kind=DP),allocatable,dimension(:) :: rwork

    allocate(hmat_k(num_wan,num_wan))

    hmat_k = cmplx(0.d0,0.d0)
    do j=1,num_wan
       do i=1,num_wan
          do ir=1,nrpts
             ph = PAI2 * dot_product(kpt,real(ipos(1:3,ir),DP))
             zf = cmplx(cos(ph),sin(ph))/real(ndegen(ir),DP)
             hmat_k(i,j) = hmat_k(i,j) + zf * hmat(i,j,ir)
          end do
       end do
    end do

    JOBZ = 'V'
    UPLO = 'U'
    lwork = max(1,2*num_wan-1)
    lrwork = max(3*num_wan-2,1)

    allocate(work(lwork))
    allocate(rwork(lrwork))
    umat = hmat_k
    call zheev(JOBZ,UPLO,num_wan,umat,num_wan,eig,work,lwork,rwork,info)
    deallocate(work)
    deallocate(rwork)

    deallocate(hmat_k)
  end subroutine solve_hmat_k
  
end program wan_iterp

function fn_number_of_words(string)
  implicit none

  integer :: fn_number_of_words
  character(*),intent(in) :: string
  integer        :: ipos, lpos
  character(256) :: tmp_string
  tmp_string = string
  lpos       = len(tmp_string)
  fn_number_of_words = 0

  if (lpos == 0) goto 99

  tmp_string = adjustl(tmp_string)
  if (index(tmp_string,' ') == 1) goto 99

  do while (index(tmp_string,' ') > 1)
     fn_number_of_words = fn_number_of_words + 1
     ipos               = index(tmp_string,' ')
     if (ipos == 0) goto 99
     tmp_string = adjustl(tmp_string(ipos:lpos))
     tmp_string = adjustl(tmp_string)
  end do
99 continue

end function fn_number_of_words
