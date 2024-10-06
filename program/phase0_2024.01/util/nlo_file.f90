       program main
       implicit none
!      Program of nlo susceptibility decomposition
!      This program read 'nlo.data' file containing the susceptibility data 
!      and write the tensor component specified to 'nlo_ijk.data'.
!      ijk is the tensor index
!      INPUT : nlo.data     (file no.=1)
!      OUTPUT: nlo_ijk.data (file no.=2)
!      T. Hamada, Univ. Tokyo (2005.9.1)
       integer       :: ind, i
       integer       :: icount
       character*12  :: file_name
       character*59  :: a
       character*3   :: ind_chi2
       real(kind=8)  :: energy, r_chi2, i_chi2, abs

       open(1,file = 'nlo.data',status="old")
       read(5,'(a3)') ind_chi2
       file_name = 'nlo_'//ind_chi2//'.data'
       open(2,file = file_name, status="replace")

       call read_write_header

    1  read(1,10,end=99) a
       if(index(a,ind_chi2)/=0) then
          read(1,10,end=99) a
          write(2,10) a
          icount=0
    2     read(1,10) a
          if(index(a,'*')==0) then
             backspace 1
             read(1,20) energy, r_chi2, i_chi2, abs
             write(2,20) energy, r_chi2, i_chi2, abs
             icount=icount+1
          else
             write(6,'("line number of ",a12 " = ",i5)') file_name, icount-1
             stop
          end if
          goto 2
        end if
       goto 1
    99 stop 
    10 format(a59)
    20 format(3x,f10.5,6x,3(f10.5,5x))
       end program main
 
       subroutine read_write_header
       implicit none
       character*59 :: a
       integer :: i
       rewind 1
       rewind 2
       read(1,10) a
       if(i==1.and.index(a,'Tensor')==0) then
          write(*,99)
          stop
       end if
       write(2,10) a
    10 format(a59)
    99 format("nlo.data file has no dielectric tensor component. STOP")
       end subroutine read_write_header
