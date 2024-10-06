module m_db
  use m_Files,                     only : F_INP,m_Files_reopen_nfinp
  use stringModule,                only : strSplit,splitStr
  use m_Const_Parameters,          only : LOWER,NOCONV

  implicit none
  private

  public:: getIntDB_TB,getLogDB_TB,getRealDB_TB
  public:: getStringDB_TB
  public:: getIntDB,getRealVecDB,getRealDB
  public:: getStringDB,initDB

contains

  subroutine initDB
    call m_Files_reopen_nfinp(1)
  end subroutine initDB

  subroutine getIntDB(tag,retInt,defint)
    character(len=*),intent(in):: tag
    integer,intent(out):: retInt
    integer,intent(in),optional:: defInt

    integer:: f_getIntValue,f_selectTop,f_selectBlock
    integer:: iret,nsize
    logical:: br

    call strSplit(trim(tag))
    nsize=size(splitStr)
    iret = f_selectTop()
    br=.false.
    do iret=1,nsize-1
       if(f_selectBlock(trim(splitStr(iret)))/=0) then
          br=.true.
          exit
       endif
    end do

    if(f_getIntValue(trim(splitStr(nsize)),iret)==0.and..not.br) then
       retInt=iret
    else
       if(present(defint)) then
          retInt=defInt
       else
          write(6,'(3a)') 'Error in getIntDB :',tag,' is not found.'
          stop
       end if
    end if

  end subroutine getIntDB

  subroutine getIntDB_TB(tag,int_vec,def_vec)
    character(len=*),intent(in):: tag
    integer,dimension(:),intent(out):: int_vec
    integer,dimension(:),intent(in),optional :: def_vec
    integer:: iret,i,matSize,nsize
    integer:: f_selectTop,f_selectBlock
    integer:: f_selectFirstTableLine
    integer:: f_selectNextTableLine
    integer:: f_getIntValue

    matSize=size(int_vec)
    call strSplit(trim(tag))
    nsize=size(splitStr)
    iret = f_selectTop()
    do iret=1,nsize-1
       if(f_selectBlock(trim(splitStr(iret))) /= 0) go to 10
    end do

    i = 1
    do while(.true.)
       if (i == 1) then
          if(f_selectFirstTableLine() /= 0) then
             exit
          end if
       else
          if(f_selectNextTableLine() /= 0) then
             exit
          end if
       end if
       if(i > matSize) exit
       if( f_getIntValue(trim(splitStr(nsize)),iret) == 0) &
            int_vec(i) = iret
       i=i+1
    end do
!    if(i-1.ne.matSize) then
!       write(6,*) 'Error in getIntDB_TB:',i-1,matSize
!       stop
!    end if
    return

10  continue
    if(present(def_vec)) then
       int_vec(:)=def_vec(:)
    else
       write(6,'(3a)') 'Error in getIntDB_TB: ',tag,' is not found!'
       stop
    end if
 
  end subroutine getIntDB_TB

  subroutine getLogDB_TB(tag,log_vec,def_vec)
    character(len=*),intent(in):: tag
    logical,dimension(:),intent(out):: log_vec
    logical,dimension(:),intent(in),optional :: def_vec
    integer:: iret,i,matSize,nsize
    integer:: f_selectTop,f_selectBlock
    integer:: f_selectFirstTableLine
    integer:: f_selectNextTableLine
    integer:: f_getIntValue

    matSize=size(log_vec)
    call strSplit(trim(tag))
    nsize=size(splitStr)
    iret = f_selectTop()
    do iret=1,nsize-1
       if(f_selectBlock(trim(splitStr(iret))) /= 0) go to 10
    end do

    i = 1
    do while(.true.)
       if (i == 1) then
          if(f_selectFirstTableLine() /= 0) then
             exit
          end if
       else
          if(f_selectNextTableLine() /= 0) then
             exit
          end if
       end if
       if(i > matSize) exit
       if( f_getIntValue(trim(splitStr(nsize)),iret) == 0) then
          if(iret.eq.1) then
             log_vec(i)=.true.
          else if(iret.eq.0) then
             log_vec(i)=.false.
          else
             write(6,*) 'Error in getLogDB_TB: ',iret
             stop
          end if
       end if
       i=i+1
    end do
!    if(i-1.ne.matSize) then
!       write(6,*) 'Error in getIntDB_TB:',i-1,matSize
!       stop
!    end if
    return

10  continue
    if(present(def_vec)) then
       log_vec(:)=def_vec(:)
    else
       write(6,'(3a)') 'Error in getIntDB_TB: ',tag,' is not found!'
       stop
    end if
 
  end subroutine getLogDB_TB

  subroutine getRealDB_TB(tag,re_vec,def_vec)
    character(len=*),intent(in):: tag
    real(8),dimension(:),intent(out):: re_vec
    real(8),dimension(:),intent(in),optional :: def_vec
    integer:: iret,i,matSize,nsize
    integer:: f_selectTop,f_selectBlock
    integer:: f_selectFirstTableLine
    integer:: f_selectNextTableLine
    integer:: f_getRealValue
    real(8):: dret

    matSize=size(re_vec)
    call strSplit(trim(tag))
    nsize=size(splitStr)
    iret = f_selectTop()
    do iret=1,nsize-1
       if(f_selectBlock(trim(splitStr(iret))) /= 0) go to 10
    end do

    i = 1
    do while(.true.)
       if (i == 1) then
          if(f_selectFirstTableLine() /= 0) then
             exit
          end if
       else
          if(f_selectNextTableLine() /= 0) then
             exit
          end if
       end if
       if(i > matSize) exit
       if( f_getRealValue(trim(splitStr(nsize)),dret,'') == 0) &
            re_vec(i) = dret
       i=i+1
    end do
!    if(i-1.ne.matSize) then
!       write(6,*) 'Error in getIntDB_TB:',i-1,matSize
!       stop
!    end if
    return

10  continue
    if(present(def_vec)) then
       re_vec(:)=def_vec(:)
    else
       write(6,'(3a)') 'Error in getRealDB_TB: ',tag,' is not found!'
       stop
    end if
 
  end subroutine getRealDB_TB

  subroutine getStringDB_TB(tag,str_vec,def_vec)
    character(len=*),intent(in):: tag
    character(len=64),dimension(:),intent(out):: str_vec
    character(len=64),dimension(:),intent(in),optional :: def_vec
    integer:: iret,i,matSize,nsize
    integer:: f_selectTop,f_selectBlock
    integer:: f_selectFirstTableLine
    integer:: f_selectNextTableLine
    integer:: f_getStringValue
    character(len=64):: str_ret

    matSize=size(str_vec)
    call strSplit(trim(tag))
    nsize=size(splitStr)
    iret = f_selectTop()
    do iret=1,nsize-1
       if(f_selectBlock(trim(splitStr(iret))) /= 0) go to 10
    end do

    i = 1
    do while(.true.)
       if (i == 1) then
          if(f_selectFirstTableLine() /= 0) then
             exit
          end if
       else
          if(f_selectNextTableLine() /= 0) then
             exit
          end if
       end if
       if(i > matSize) exit
       if( f_getStringValue(trim(splitStr(nsize)),str_ret,NOCONV) == 0) &
            str_vec(i) = str_ret
       i=i+1
    end do
!    if(i-1.ne.matSize) then
!       write(6,*) 'Error in getIntDB_TB:',i-1,matSize
!       stop
!    end if
    return

10  continue
    if(present(def_vec)) then
       str_vec(:)=def_vec(:)
    else
       write(6,'(3a)') 'Error in getStringDB_TB: ',tag,' is not found!'
       stop
    end if
 
  end subroutine getStringDB_TB

  subroutine getRealVecDB(tag,retVec,defVec)
    character(len=*),intent(in):: tag
    real(8),intent(out):: retVec(3)
    real(8),intent(in),optional:: defVec(3)

    integer:: f_selectTop,f_selectBlock,getRealVectorValue 
	
    integer:: iret,nsize
    real(8):: tmp_vec(3)
    character(16) :: readunit
    logical :: br
!!$interface 
!!$	integer function getRealVectorValue(ttag,rret,unit)
!!$	!DEC$ ATTRIBUTES C, ALIAS: '_getrealvectorvalue' :: getRealVectorValue
!!$	!DEC$ ATTRIBUTES REFERENCE :: ttag
!!$	!DEC$ ATTRIBUTES REFERENCE :: unit
!!$	!DEC$ ATTRIBUTES REFERENCE :: rret
!!$	character*(*):: ttag
!!$	real(8):: rret(3)
!!$	character*(*):: unit	
!!$	end function getRealVectorValue
!!$end interface

    call strSplit(trim(tag))
    nsize=size(splitStr)
    iret = f_selectTop()
    br=.false.
    do iret=1,nsize-1
       if(f_selectBlock(trim(splitStr(iret)))/=0) then
          br=.true.
          exit
       endif
    end do
    if(getRealVectorValue(trim(splitStr(nsize))//char(0), &
                                  tmp_vec,readunit)==0.and..not.br) then
       retVec(1:3)=tmp_vec(1:3)
    else
       if(present(defVec)) then
          retVec(1:3)=defVec(1:3)
       else
          write(6,'(3a)') 'Error in getRealVecDB :',tag,' is not found.'
          stop
       end if
    end if

  end subroutine getRealVecDB

  subroutine getRealDB(tag,retRe,defRe)
    character(len=*),intent(in):: tag
    real(8),intent(out):: retRe
    real(8),intent(in),optional:: defRe

    integer:: f_getRealValue,f_selectTop,f_selectBlock
    integer:: iret,nsize
    real(8):: dret
    logical:: br

    call strSplit(trim(tag))
    nsize=size(splitStr)
    iret = f_selectTop()
    br=.false.
    do iret=1,nsize-1
       if(f_selectBlock(trim(splitStr(iret)))/=0) then
         br=.true.
         exit
       endif
    end do

    if(f_getRealValue(trim(splitStr(nsize)),dret,'')==0.and..not.br) then
       retRe=dret
    else
       if(present(defRe)) then
          retRe=defRe
       else
          write(6,'(3a)') 'Error in getRealDB :',tag,' is not found.'
          stop
       end if
    end if

  end subroutine getRealDB

  subroutine getStringDB(tag,retStr,defStr)
    character(len=*),intent(in):: tag
    character(len=*),intent(out):: retStr
    character(len=*),intent(in),optional:: defStr

    integer:: f_getStringValue,f_selectTop,f_selectBlock
    integer:: iret,nsize
    character(len=64):: str_ret
    logical:: br

    call strSplit(trim(tag))
    nsize=size(splitStr)
    iret = f_selectTop()
    br=.false.
    do iret=1,nsize-1
       if(f_selectBlock(trim(splitStr(iret)))/=0) then
          br=.true.
          exit
       endif
    end do

    if(f_getStringValue(trim(splitStr(nsize)),str_ret,LOWER)==0.and..not.br) then
       retStr=str_ret
    else
       if(present(defStr)) then
          retStr=defStr
       else
          write(6,'(3a)') 'Error in getStringDB :',tag,' is not found.'
          stop
       end if
    end if

  end subroutine getStringDB

end module m_db
