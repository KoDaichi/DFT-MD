subroutine UpdateUeff()
  use m_Control_Parameters, only : printable, num_projectors,proj_attribute, nUeff, driver
  use m_Const_Parameters, only : DRIVER_URAMP
  use m_Files, only : nfout
  use m_IterationNumbers, only : iteration_uramp
  implicit none
  integer :: i
  logical :: Uramping
  if (.not.Uramping()) return
  if (printable) write(nfout,'(a,i5,a,i5)') ' !** Uramping no ',iteration_uramp,'/',nUeff
  do i=1,num_projectors
     proj_attribute(i)%Ueff = (iteration_uramp-1) * proj_attribute(i)%deltaUeff + proj_attribute(i)%initialUeff
     if(printable) write(nfout,'(a,i5,a,f10.5)') ' !** new Ueff for projector',i,':',proj_attribute(i)%Ueff
  enddo
end subroutine UpdateUeff

function InUeffRamping() result(ret)
  use m_Const_Parameters, only : DRIVER_URAMP
  use m_Control_Parameters, only : driver, nUeff
  use m_IterationNumbers, only : iteration_uramp
  implicit none
  logical :: ret
  if (driver /= DRIVER_URAMP) then
    ret = .false.
    return
  endif
  if (iteration_uramp >= nUeff) then
    ret = .false.
    return
  endif 
  ret = .true.
end function InUeffRamping

function Uramping() result(ret)
  use m_Const_Parameters, only : DRIVER_URAMP
  use m_Control_Parameters, only : driver, nUeff
  implicit none
  logical :: ret
  ret = driver == DRIVER_URAMP
end function Uramping

