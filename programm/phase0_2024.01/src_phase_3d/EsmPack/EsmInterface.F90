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
subroutine Esm_interface_map_parameters(nat_,upf_zp_,tau_,alat_,at_,nr1x_,nr2x_,nr3x_,  &
                                   &    esm_bc_,gamma_only_,ngm_,nspin_,mill_,nl_,nlm_, &
                                   &    esm_w_,esm_efield_,izwall_,z_wall_,bar_height_,bar_width_)
  use ESM_VARS
  implicit none
  integer, intent(in) :: nat_
  real(8), dimension(nat_), intent(in)   :: upf_zp_
  real(8), dimension(3,nat_), intent(in) :: tau_
  real(8), intent(in) :: alat_
  real(8), intent(in) :: at_(3,3)
  integer, intent(in) :: nr1x_,nr2x_,nr3x_
  character(len=3), intent(in) :: esm_bc_
  logical, intent(in) :: gamma_only_
  integer, intent(in) :: ngm_,nspin_
  integer, dimension(3,ngm_), intent(in) :: mill_
  integer, dimension(ngm_), intent(in) :: nl_,nlm_
  real(8), intent(in) :: esm_w_
  real(8), intent(in) :: esm_efield_
  integer, intent(in) :: izwall_
  real(8), intent(in) :: z_wall_
  real(8), intent(in) :: bar_height_
  real(8), intent(in) :: bar_width_

  Real(8) :: at12(3)
  Real(8), External :: InnerProd
  
  integer :: i,j
  Real(8) :: a1,a2,cc

  nat = nat_
  if(.not.allocated(tau)) allocate(tau(3,nat))
  if(.not.allocated(upf_zp)) allocate(upf_zp(nat))
  do i=1,nat
     do j=1,3
        tau(j,i) = tau_(j,i)
     enddo
     upf_zp(i) = upf_zp_(i)
  enddo
  alat = alat_
  at = at_
  nr1x = nr1x_
  nr2x = nr2x_
  nr3x = nr3x_
  esm_bc = esm_bc_
  gamma_only = gamma_only_
  ngm = ngm_
  nspin = nspin_
  if(.not.allocated(mill)) allocate(mill(3,ngm))
  if(.not.allocated(nl))   allocate(nl(ngm))
  if(.not.allocated(nlm))  allocate(nlm(ngm))

  mill(:,:) = mill_(:,:)
  nl(:)     = nl_(:)
  nlm(:)    = nlm_(:)

  esm_w      = esm_w_
  esm_efield = esm_efield_

  izwall     = izwall_
  z_wall     = z_wall_
  bar_height = bar_height_
  bar_width  = bar_width_

! calculate parameters
  Call OuterProd( at(1,1), at(1,2), at12 )
  omega = InnerProd( at12, at(1,3) ) * ( alat ** 3 )

!  gew = sqrt(pi/omega/(at(3,3)*alat))

  tpiba2 = ( tpi / alat ) ** 2
  !!$gcutm = dual * ecutwfc / tpiba2

  Call MakeRecLatVec( at, bg )

  nrxx = nr1x * nr2x * nr3x

  allocate(rhog_(ngm,nspin))

  call esm_ggen_2d()

  call esm_gen_vbar()

  return
end subroutine Esm_interface_map_parameters

subroutine Esm_interface_set_tau(tau_)
  use ESM_VARS
  implicit none
  real(kind=8), intent(in), dimension(3,nat) :: tau_
  integer :: i,j
  tau(:,:) = tau_(:,:)
end subroutine Esm_interface_set_tau

subroutine Esm_interface_set_communicator(comm)
  use ESM_VARS
  implicit none
  integer, intent(in) :: comm
  communicator = comm 
end subroutine Esm_interface_set_communicator

subroutine Esm_interface_set_npes_mype(np,mp)
  use ESM_VARS
  implicit none
  integer, intent(in) :: np,mp
  npes = np
  mype = mp
end subroutine Esm_interface_set_npes_mype
