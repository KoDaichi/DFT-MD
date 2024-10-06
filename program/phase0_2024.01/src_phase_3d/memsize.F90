subroutine memsize(nfout)
implicit none
integer, intent(in) :: nfout
integer :: ipid, imemsize
character*16 :: cpid, c1, c2
character*256 :: status_file
real(kind=8) :: mb
integer :: getpid

ipid = getpid()
write(cpid,*) ipid

write(status_file,*) '/proc/', trim(adjustl(cpid)), '/status'

open(1, file=trim(adjustl(status_file)), action='read')
rewind 1
do while(.true.)
   read(1,*,end=999) c1, c2
   if(c1(1:6) .eq. 'VmHWM:') then
      read(c2,*) imemsize
      exit
   end if
end do

999 continue

close(1)

mb = dble(imemsize) / 1024.d0 ! with MB

!print '(1x,a,f10.3)', 'Memory Usage[MB]: ', mb
write(nfout, '(1x,a,f10.3)') 'Memory Usage[MB]: ', mb

return
end subroutine memsize
