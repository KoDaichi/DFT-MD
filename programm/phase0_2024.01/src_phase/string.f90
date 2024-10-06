module stringModule
  implicit none
  private
  public:: readLine, strSplit
  character(len=120),dimension(:),allocatable,public:: splitStr
contains

  function readLine(unit,line)
    logical:: readLine
    integer,intent(in) :: unit
    character(len=*),intent(out):: line
    integer:: ios,pos
    character(len=1024)::buf
    do 
      read(unit,'(A1024)',iostat=ios) buf
      if(ios/=0) then
        readLine = .true.
        return
      endif
      do                              ! replace TAB with SPACE
        pos = scan(buf,'	')
        if (pos == 0) exit
        buf(pos:pos) = ' '
      enddo
      pos = scan(buf,'#')       
      if(pos/=0) buf = buf(:pos-1)               ! erase comment
      if (len_trim(buf)/=0) then
        line = trim(buf)
        readLine = .false.
        return
      endif
    enddo
  end function readLine

  subroutine strSplit(str)
    character(len=*),intent(in):: str
    integer:: i, iarg, ist,ied,i1, splitCount
    logical:: space

    if(allocated(splitStr)) deallocate(splitStr)
    splitCount = 0

    if(len_trim(str)==0) return
    i1 = len_trim(str)

    space = .true.
    iarg = 0
    do i=1,len_trim(str)
      if(space) then
        if(str(i:i)/=' ') then
          space = .false.
          iarg = iarg + 1
        endif
      else
        if(str(i:i)==' ') space = .true.
      endif
    enddo
    splitCount = iarg
    allocate(splitStr(splitCount))
    
    space = .true.
    iarg = 1
    do i=1,len_trim(str)
      if(space) then
        if(str(i:i)/=' ') then
          ist = i
          space = .false.
        endif
      else
        if(str(i:i)==' ') then
          ied = i-1
          splitStr(iarg) = str(ist:ied)
          iarg = iarg + 1
          space = .true.
        endif
      endif
    enddo
    ied = len_trim(str)
    splitStr(iarg) = str(ist:ied)
  end subroutine strSplit

end module stringModule

!program test
!  use misc
!  use stringModule
!  character(len=120):: aLine
!  integer:: i
!  read(*,'(A)') aLine
!  call strSplit(aLine)
!  do i=1,splitCount
!    write(6,*) ml, i,"<",trim(splitStr(i)),">"
!  end do
!end program test