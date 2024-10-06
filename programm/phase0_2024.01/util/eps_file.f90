       program main
       implicit none
!      Program of epsilon file converesion
!      This program read 'eps.data' file containing dietric tensor
!      and write 'r_eps.data' and 'i_eps.data' files containing real and
!      imaginary parts of the dielectric tensor, respectively.
!      INPUT : eps.data   (file no.=1)
!      OUTPUT: r_eps.data (file no.=2)
!            : i_eps.data (file no.=3)
!      T. Hamada, Univ. Tokyo (2003.06.12)
!      iline : output control parameter
!            : write down data every iline
!            : default==1
!      T. Hamada, Univ. Tokyo (2004.1.12)
       integer :: ind, i, len
       integer :: iline,ilcount1,ilcount2
       character*80  :: file_name
       character*105 :: a
       real(kind=8) :: energy, eps_xx, eps_yy, eps_zz, eps_xy, eps_xz, eps_yz
       open(1,file='eps.data',status="old")
       open(2,file='r_eps.data',status="replace")
       open(3,file='i_eps.data',status="replace")

! read skipped line number -> not used at present
!       write(6,'("skip line = ?")')
!       read(5,'(i3)') iline
!       if(iline==0) iline=1
       iline = 1
       call read_write_header
       ilcount1=0
       ilcount2=0
     1 read(1,10,end=99) a
       ind=index(a,'(')
       if(ind==0) then
        backspace 1
        read(1,20) energy, eps_xx, eps_yy, eps_zz, eps_xy, eps_xz, eps_yz
        if(mod(ilcount1,iline)==0)  write(2,20) energy, eps_xx, eps_yy, eps_zz, eps_xy, eps_xz, eps_yz
        ilcount1=ilcount1+1
       else
        backspace 1
        read(1,30)  eps_xx, eps_yy, eps_zz, eps_xy, eps_xz, eps_yz
        if(mod(ilcount2,iline)==0)  write(3,20) energy, eps_xx, eps_yy, eps_zz, eps_xy, eps_xz, eps_yz
        ilcount2=ilcount2+1
       end if
       goto 1
    99 write(6,'("line number of r_eps.data = ",i5)') ilcount1
       write(6,'("line number of i_eps.data = ",i5)') ilcount2
       stop
    10 format(a105)
    20 format(3x,f10.5,6x,6(f10.5,5x))
    30 format(18x,6(1x,f10.5,4x))
       end program main
 
       subroutine read_write_header
       implicit none
       character*105 ::a
       integer :: i
       rewind 1
       rewind 2
       rewind 3
       do i=1,3
        read(1,10) a
        if(i==1.and.index(a,'Tensor')==0) then
         write(*,99)
         stop
        end if
        write(2,10) a
        write(3,10) a
       end do
    10 format(a105)
    99 format("eps.data file has no dielectric tensor component. STOP")
       end subroutine read_write_header
