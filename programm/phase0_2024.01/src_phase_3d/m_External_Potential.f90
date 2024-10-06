
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

end module m_External_Potential

! === "phase" repository merge: To make 3D_Parallel by tkato 2011/11/10 ========
! ==============================================================================


