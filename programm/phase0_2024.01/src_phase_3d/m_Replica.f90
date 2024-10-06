
!   Multi Replica
!   Nudged Elastic Band Method
!   Nobutaka Nishikawa

module m_Replica

  use m_Ionic_System,  only : number_of_replicas
  use m_Const_Parameters,   only : DP

  implicit none

  type unit_cell_type

     real(kind=DP)  a, b, c, alpha, beta, gamma

  end type unit_cell_type

  type symm_type
    integer :: nopr
    real(kind=DP),pointer,dimension(:,:,:) :: op
    real(kind=DP),pointer,dimension(:,:,:) :: tau
  end type symm_type

  type image_type

     integer num_atom

     type(unit_cell_type)  unit_cell
     type(symm_type) symmetry

     real(kind=DP),pointer,dimension(:,:) :: pos, pos0, cps, cps0
     real(kind=DP),pointer,dimension(:,:) :: velocity, velocity0
     real(kind=DP),pointer,dimension(:,:) :: force, force0
     real(kind=DP),pointer,dimension(:,:) :: force_org, true_force, spring_force, tau

     real(kind=DP)  energy, energy0
     !real(8)  total_energy
     !real(8)  potential_energy

     real(8)  k_spring, k_spring0

     integer, pointer :: fix_flag(:)

     logical scf_convergence

  end type image_type

  type neb_cond_type

     integer  condition
     integer  ci_neb_start, ci_neb_end
     integer :: ci_index
     real(kind=8) :: ci_thres
     real(kind=8) :: sd_factor

     logical  climbing_image
     real(8)  k_spring_init, k_spring_min, k_spring_max

     logical k_constant
     logical k_variable
     logical k_damping
     real(8) k_damping_factor

     logical penalty_function

     integer convergence_condition
     real(8) convergence_threshold

     character(20) time_integral

     logical :: end0_energy_given
     logical :: end1_energy_given

     real(kind=8) :: dimer_thres

  end type neb_cond_type

  type neb_file_type

     character(100)  F_OUT
     character(100)  F_DYNM
     character(100)  F_ENF
     character(100)  F_CHGT
     character(100)  F_CNTN
     character(100)  F_CNTN_BIN
     character(100)  F_ZAJ
     character(100)  F_CHR

     character(100)  F_NEB_OUT

! ======================== KT_add ================= 13.0D
     character(100)  F_CNTN_BIN_PAW
! ================================================= 13.0D

     character(100)  F_OCCMAT

  end type neb_file_type

  type neb_type

     integer  number_of_images
     integer  max_iteration

     type(image_type), pointer :: image(:)
     type(neb_cond_type) :: cond
     type(neb_file_type) :: file

     integer  step
     real(kind=DP)  time
     real(kind=DP)  dt

     ! temporary
     real(kind=DP),pointer,dimension(:,:) :: energy

  end type neb_type

  type(neb_type)  neb
  integer nrank_r, mype_r
  integer, allocatable :: npes_image(:)

  logical :: first_replica_done = .false.
  real(kind=DP) :: curr_thres = 1e+30
  logical :: changed_to_cineb=.false.

  integer :: nopr_lowest
  real(kind=DP), allocatable, dimension(:,:,:) :: op_lowest,tau_lowest

  contains

  subroutine m_Replica_resolve_lowest_symmetry(nfneb)
    use m_Parallelization, only : mype
    integer, intent(in) :: nfneb
    integer :: i,j,k,l,m,n,currindex,currsize
    integer :: icount
    real(kind=DP), allocatable, dimension(:,:,:) :: opbuf,taubuf

    currsize = neb%image(1)%symmetry%nopr
    allocate(opbuf(3,3,currsize))
    allocate(taubuf(3,currsize,2))
    allocate(op_lowest(3,3,currsize));op_lowest=neb%image(1)%symmetry%op
    allocate(tau_lowest(3,currsize,2));tau_lowest=neb%image(1)%symmetry%tau
    do i=2,neb%number_of_images
       currindex=0
       do j=1,currsize
          kloop:do k=1,neb%image(i)%symmetry%nopr
             icount = 0
             do m=1,3
                if (epsilon_equals(tau_lowest(m,j,1),neb%image(i)%symmetry%tau(m,k,1))) then
                   icount = icount+1 
                endif
                if (epsilon_equals(tau_lowest(m,j,2),neb%image(i)%symmetry%tau(m,k,2))) then
                   icount = icount+1 
                endif
                do l=1,3
                   if (epsilon_equals(op_lowest(l,m,j),neb%image(i)%symmetry%op(l,m,k))) then
                      icount = icount+1 
                   endif
                enddo
             enddo
             if (icount==15) then
               currindex = currindex+1
               opbuf(:,:,currindex) = op_lowest(:,:,j)
               taubuf(:,currindex,:) = tau_lowest(:,j,:)
               exit kloop
             endif
          enddo kloop
       enddo
       currsize = currindex
       do j=1,currsize
          op_lowest(:,:,j) = opbuf(:,:,j)
          tau_lowest(:,j,:) = taubuf(:,j,:)
       enddo
    enddo
    nopr_lowest = currsize

    if(mype==0) then
      write(nfneb,*) 'nopr ',nopr_lowest
      do i=1,nopr_lowest
         write(nfneb,*) 'op no. ',i
         do j=1,3
            write(nfneb,'(3f10.5)') op_lowest(j,1:3,i)
         enddo
         write(nfneb,*) 'tau no. ',i
         write(nfneb,'(6f10.5)') tau_lowest(1:3,i,1),tau_lowest(1:3,i,2)
      enddo
    endif
    deallocate(opbuf)
    deallocate(taubuf)

    contains 

    logical function epsilon_equals(a,b)
        real(kind=DP), intent(in) :: a,b
        real(kind=DP) :: eps = 1.d-8
        epsilon_equals = abs(a-b)<eps
    end function epsilon_equals

  end subroutine m_Replica_resolve_lowest_symmetry

  subroutine m_Replica_finalize()
    integer :: i,j
    do i=1,neb%number_of_images
      deallocate(neb%image(i)%pos)
      deallocate(neb%image(i)%pos0)
      deallocate(neb%image(i)%cps)
      deallocate(neb%image(i)%cps0)
      deallocate(neb%image(i)%tau)
      deallocate(neb%image(i)%force)
      deallocate(neb%image(i)%force0)
      deallocate(neb%image(i)%force_org)
      deallocate(neb%image(i)%spring_force)
      deallocate(neb%image(i)%true_force)
      deallocate(neb%image(i)%fix_flag)
      deallocate(neb%image(i)%velocity)
      deallocate(neb%image(i)%velocity0)
      deallocate(neb%image(i)%symmetry%op)
      deallocate(neb%image(i)%symmetry%tau)
    enddo
!!    deallocate(neb%energy)
    deallocate(neb%image)
  end subroutine m_Replica_finalize

  logical function FirstReplicaDone()
    FirstReplicaDone = first_replica_done
  end function FirstReplicaDone

  function get_max_force(neb,i,mode) result(fmax)
    type(neb_type), intent(in) :: neb
    integer, intent(in) :: i
    character(len=*), intent(in) :: mode
    real(kind=8) :: ftmp,fmax
    integer :: j
    fmax = 0.d0
    do j=1, neb%image(i)%num_atom
      if(neb%image(i)%fix_flag(j)==0)cycle
      if (mode=='true') then
      ftmp = dsqrt(dot_product(neb%image(i)%true_force(j,:),neb%image(i)%true_force(j,:)))
      else if (mode == 'org') then
      ftmp = dsqrt(dot_product(neb%image(i)%force_org(j,:),neb%image(i)%force_org(j,:)))
      else
      ftmp = dsqrt(dot_product(neb%image(i)%force(j,:),neb%image(i)%force(j,:)))
      endif
      if(fmax < ftmp) fmax=ftmp
    enddo
  end function get_max_force

end module m_Replica
