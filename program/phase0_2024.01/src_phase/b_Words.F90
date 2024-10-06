!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: skip2number, f_strcpy, skip2RHS, read_LHS_number, read_RHS_number,
!             rewind_to_tag0, f_strncmp, ch4prc, chnnm2, skip2alph, ch4mix, 
!             rd_alpha, rd_cutoff, rd_hownew, rd_nbxmix, rd_istrbr, 
!             check_mix_tag_name, check_decomp, nameck, skip2next, 
!             strncmp2, strncmp0, chchn0, isnumm, isnums, isalph, chnnm,
!             skip_to_subtagname, read_energy_unit, read_itagvalue, read_dtagvalue
!
!  FUNCTION:  skip_to_tagbegin, skip_to_tagbegin2, strncmp3
!
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  Contact address :  Phase System Consortium
!                     E-mail: phase_system@nims.go.jp URL https://azuma.nims.go.jp
!  
!
!
!=======================================================================
!
!     The original version of this set of the computer programs "PHASE"
!  was developed by the members of the Theory Group of Joint Research
!  Center for Atom Technology (JRCAT), based in Tsukuba, in the period
!  1993-2001.
!
!     Since 2002, this set has been tuned and new functions have been
!  added to it as a part of the national project "Frontier Simulation 
!  Software for Industrial Science (FSIS)",  which is supported by
!  the IT program of the Ministry of Education, Culture, Sports,
!  Science and Technology (MEXT) of Japan. 
!     Since 2006, this program set has been developed as a part of the
!  national project "Revolutionary Simulation Software (RSS21)", which
!  is supported by the next-generation IT program of MEXT of Japan.
!   Since 2013, this program set has been further developed centering on PHASE System
!  Consortium.
!   The activity of development of this program set has been supervised by Takahisa Ohno.
!
!                           @(#)b_Words.f90 1.3 03/02/19 00:46:15
!***************************************************************
!
!  Subroutines in this file "b_Words" are all coded by T. Yamasaki
!  in 1993 - 1999
!
! $Id: b_Words.f90 570 2017-04-21 20:34:50Z yamasaki $
logical function skip_to_tagbegin(nfinp,tag)
  integer, intent(in)          :: nfinp
  character(len=*), intent(in) :: tag
      
  integer, parameter :: len_str = 132
  character(len=len_str) :: str
  logical :: tf

  rewind nfinp
1 read(nfinp,'(a132)',end = 1001) str
  call strncmp2(str,len(tag),tag,len(tag),tf)
  if(.not.tf) goto 1
1001 continue

  skip_to_tagbegin = tf

end function skip_to_tagbegin

logical function skip_to_tagbegin2(nfinp,tag)
  integer, intent(in)         :: nfinp
  character(len=*),intent(in) :: tag

  integer :: ip
  integer, parameter :: len_str = 132
  character(len=len_str) :: str

  skip_to_tagbegin2 = .false.
  rewind nfinp
  ip = 0
1 read(nfinp,'(a132)',end=1001) str
  ip = index(str,tag)
!!$      write(6,'(i4,a76)') ip,str(1:76)
  if(ip <= 0) goto 1
1001 continue
  if(ip >= 1) skip_to_tagbegin2 = .true.
  write(6,'(" ip (skip_to_tagbegin2) = ",i4)') ip
end function skip_to_tagbegin2

subroutine read_dtagvalue(nfinp,tag_begin,tag,value)
  use m_Const_Parameters, only : DP
  integer, intent(in)          :: nfinp
  character(len=*), intent(in) :: tag_begin,tag
  real(kind=DP), intent(inout) :: value
      
  logical :: tf, tf2, skip_to_tagbegin
  integer :: nsize
  real(kind=DP), dimension(1) :: work
  
  integer, parameter :: len_str = 132
  character(len=len_str) :: str

  if(.not.skip_to_tagbegin(nfinp,tag_begin)) return

  tf = .false.
1 read(nfinp,'(a132)',end = 1001) str
  call strncmp2(str,3,"end",3,tf2)
  if(tf2) goto 1001
  call strncmp2(str,len_str,tag,len(tag),tf)
  if(.not.tf) goto 1
1001 continue

  if(tf) then
!!$         write(6,'(" -- tag is found -- ",a31)') tag
     call chnnm2(str,len_str,1,work,nsize)
     value = work(1)
  else
     write(6,'(" -- no tag value -- ",a31)') tag
  end if
!!$      write(6,'(20x,a31," = ",d20.12)') tag, value
  write(6,'(3x,a31," = ",d20.12)') tag, value

end subroutine read_dtagvalue

subroutine read_itagvalue(nfinp,tag_begin,tag,value)
  use m_Const_Parameters, only : DP
  integer, intent(in)          :: nfinp
  character(len=*), intent(in) :: tag_begin,tag
  integer, intent(inout)       :: value
      
  logical :: tf, tf2, skip_to_tagbegin
  integer :: nsize

  real(kind=DP), dimension(1)  :: work
  integer, parameter :: len_str = 132
  character(len=len_str) :: str

  if(.not.skip_to_tagbegin(nfinp,tag_begin)) return

  tf = .false.
1 read(nfinp,'(a132)',end = 1001) str
  call strncmp2(str,3,"end",3,tf2)
  if(tf2) goto 1001
  call strncmp2(str,len_str,tag,len(tag),tf)
  if(.not.tf) goto 1
1001 continue

  if(tf) then
     call chnnm2(str,len_str,1,work,nsize)
     value = nint(work(1))
  else
     write(6,'(" -- no tag value -- ",a31,"  default value is set")') tag
  end if
  write(6,'(3x,a31," = ",i8)') tag, value

end subroutine read_itagvalue

subroutine read_energy_unit(nfinp,tag_begin,tag,energy_unit)
  use m_Const_Parameters, only : HR_ENERGY_UNIT, EV_ENERGY_UNIT
  integer, intent(in) ::               nfinp
  character(len=*), intent(in) ::      tag_begin,tag
  integer, intent(out) ::              energy_unit

  logical :: tf
  integer :: ip, strncmp3
  integer, parameter :: len_str = 132
  character(len=len_str) :: str

  call skip_to_subtagname(nfinp, tag_begin,tag,len_str,str,tf)
  if(tf) then
     ip = strncmp3(str,"=")
     if(strncmp3(str(ip+1:len_str),"hartree") >=1 ) then
        energy_unit = HR_ENERGY_UNIT
        write(6,'(" !  energy_unit (dos)= HARTREE")')
     else if(strncmp3(str(ip+1:len_str),"eV") >= 1 ) then
        energy_unit = eV_ENERGY_UNIT
        write(6,'(" !  energy_unit (dos)= eV")')
     end if
  end if
end subroutine read_energy_unit

subroutine skip_to_subtagname(nfinp,tag_begin,tag,len_str,str,tf)
  integer, intent(in) ::               nfinp,len_str
  character(len=*), intent(in) ::      tag_begin,tag
  character(len=len_str),intent(out):: str
  logical, intent(out) ::              tf
      
  logical :: tf1, tf2, skip_to_tagbegin

  if(.not.skip_to_tagbegin(nfinp,tag_begin)) return

  tf1 = .false.
1 read(nfinp,'(a132)',end = 1001) str
  call strncmp2(str,3,"end",3,tf2)
  if(tf2) goto 1001
  call strncmp2(str,len_str,tag,len(tag),tf1)
  if(.not.tf1) goto 1
1001 continue
  tf = tf1
end subroutine skip_to_subtagname

subroutine chnnm(str,len_str,n,wka,incunt)
  use m_Const_Parameters, only : DP
  implicit none
  character(len=*), intent(in) :: str
  integer, intent(in) :: len_str, n
  real(kind=DP), intent(out), dimension(n) :: wka
  integer, intent(out) :: incunt
  
  logical tf
  character(len=1) c
!!$  integer, parameter :: maxclm = 10
  real(kind=DP) :: xnum1
  integer icount
  integer i,ic2
  
  xnum1 = 0.d0
  call rsreal(n,wka)
  icount = 0
  incunt = 0
  
1001 continue
  do i = icount+1, len_str
     icount = i
     c = str(i:i)
     if(c.eq.':') goto 1006
     call isalph(c,tf)
     if(tf) goto 1006
     call isnumm(c,tf)
     if(tf) then
        read(str(i:len_str),*) xnum1
!          call chtonm(str(i:len_str),maxclm,xnum1)
        goto 1002
     endif
  enddo
  goto 1006
1002 continue
  incunt = incunt + 1
  wka(incunt) = xnum1

  if(incunt.eq.n) goto 1006
  
  ic2 = icount + 1
  do i = icount+1, len_str
     if(str(i:i).eq.' ') then
        ic2 = i
        goto 1004
     else if(str(i:i).eq.':') then
        goto 1006
     endif
  enddo
1004 continue
  icount = ic2
  goto 1001
1006 continue
end subroutine chnnm

subroutine isalph(chr,trorfa)
  implicit none
  character, intent(in):: chr
  logical, intent(out) :: trorfa
  if((ichar(chr) >= ichar('a') .and. ichar(chr) <= ichar('z')) &
       &     .or.(ichar(chr) >= ichar('A') .and. ichar(chr) <= ichar('Z')) &
       &     ) then
     trorfa = .true.
  else
     trorfa = .false.
  endif
end subroutine isalph

subroutine isnums(chr,trorfa )
  implicit none
  character, intent(in) :: chr
  logical ,intent(out) :: trorfa

  if(( ichar(chr).ge.ichar('0').and.ichar(chr).le.ichar('9') ) &
       &    .or.chr.eq.'.'.or.chr.eq.'-' &
       &    .or.(chr.eq.'D'.or.chr.eq.'E'.or.chr.eq.'Q') &
       &    ) then
     trorfa = .true.
  else
     trorfa = .false.
  endif
end subroutine isnums

subroutine isnumm(chr,trorfa )
  implicit none
  character, intent(in) :: chr
  logical ,intent(out) :: trorfa
  
  if(( ichar(chr).ge.ichar('0').and.ichar(chr).le.ichar('9') ) &
       &    .or.chr.eq.'.'.or.chr.eq.'-' ) then
     trorfa = .true.
  else
     trorfa = .false.
  endif
end subroutine isnumm

subroutine chchn0(str,len_str,n,name,ns)
  implicit none
  integer, intent(in)         :: len_str,n
  character(len=len_str), intent(in) :: str
  character(len=n), intent(out)  :: name
  integer, intent(out)        :: ns
  
  integer, parameter :: DP = kind(1.d0)
  real(kind=DP) :: xnum1
  integer :: ip, i, icount, ndif
  character c
  logical tf

  ip = 0
  do i = 1, n
     icount = i
     call isalph(str(i:i),tf)
     if(tf) goto 2000
  enddo
  !stop ' >>> Invalid input data <<<'
  call phase_error_with_msg(6, ' >>> Invalid input data <<<',__LINE__,__FILE__)
2000 continue
  do i = icount, icount + n - 1
     icount = i
     if(str(i:i) /= ' '.and.str(i:i) /= ':') then
        ip     = ip + 1
        name(ip:ip) = str(i:i)
     else
        goto 2001
     endif
  enddo
2001 continue

  name(ip+1:n) = ' '
  !      do i = ip+1, n
  !        name(i:i) = ' '
  !      enddo

  ndif = ichar('a') - ichar('A')
  do i = 1, n
     c = name(i:i)
     if(ichar(c) >= ichar('A') .and. ichar(c) <= ichar('Z')) then
        name(i:i) = char(ichar(c) + ndif)
     endif
  enddo

  do i = icount+1, len_str
     icount = i
     c = str(i:i)
     if(c.eq.':') goto 2006
     call isnumm(c,tf)
     if(tf) then
        read(str(i:len_str),*) xnum1
!          call chtonm(str(i:len_str),maxclm,xnum1)
        goto 2002
     endif
  enddo
2002 continue

  do i = icount+1, len_str
     if(str(i:i) == ' ') then
        goto 2004
     else if(str(i:i) == ':') then
        goto 2006
     endif
  enddo
2004 continue
  
2006 continue
  ns   = xnum1
end subroutine chchn0

subroutine strncmp0(str1,str2,tf)
  implicit none
  character(len=*), intent(in) :: str1, str2
  logical, intent(out) ::         tf

  integer :: len1, lencmp
  integer :: ndif, i, icharAcap,icharZcap,ichara,icharz,icharx &
       &   , icharupper, icharlower
!!$  character c
!  ----------
  len1 = len(str1)

  icharAcap = ichar('A')
  icharZcap = ichar('Z')
  ichara    = ichar('a')
  icharz    = ichar('z')

  ndif = ichar('A') - ichar('a')

  lencmp = len(str2)
  if(lencmp > len1) lencmp = len1
  tf = .false.
  do i = 1, lencmp
     icharx = ichar(str2(i:i))
     if(icharAcap <= icharx .and. icharx <= icharZcap) then
        icharupper = icharx
        icharlower = icharx - ndif
     else if(ichara <= icharx .and. icharx <= icharz) then
        icharupper = icharx + ndif
        icharlower = icharx
     else
        icharupper = icharx; icharlower = icharx
     end if
     if(str1(i:i).ne.char(icharupper).and.str1(i:i).ne.char(icharlower)) goto 1000
  enddo
  tf = .true.
1000 continue
end subroutine strncmp0

subroutine strncmp2(str1,len_str,str2,n,tf)
  implicit none
  character(len=*), intent(in) :: str1, str2
  integer, intent(in)       :: len_str, n
  logical, intent(out)      :: tf

  integer :: ndif, istart, i2,i
  character c
!  ----------
  ndif = ichar('A') - ichar('a')

  tf = .false.
  do istart = 1, len_str - n + 1
     i2 = 0
     do i = istart, istart-1+n
        i2 = i2 + 1
        c = char(ichar(str2(i2:i2)) + ndif)
        if(str1(i:i).ne.str2(i2:i2).and.str1(i:i).ne.c) goto 1000
     enddo
     tf = .true.
1000 continue
     if(tf) return
  enddo

end subroutine strncmp2

subroutine strncmp2b(str1,len_str,str2,n,tf)
  implicit none
  character(len=*), intent(in) :: str1, str2
  integer, intent(in)       :: len_str, n
  logical, intent(out)      :: tf

  integer :: ndif, istart, i2,i
  character c
!  ----------
  ndif = ichar('A') - ichar('a')

  tf = .false.
  do istart = 1, len_str - n + 1
     i2 = 0
     do i = istart, istart-1+n
        i2 = i2 + 1
        c = char(ichar(str2(i2:i2)) + ndif)
        if(str1(i:i).ne.str2(i2:i2).and.str1(i:i).ne.c) goto 1000
     enddo
     tf = .true.
1000 continue
     if(tf) then
        if(istart >= 2) then
           c = str1(istart-1:istart-1)
           if(c.ne.' ') tf = .false.
        end if
        if(i <= len_str) then
           c = str1(i:i)
           if(c.ne.' ' .and. c.ne.'=') tf = .false.
        end if
     end if
     if(tf) return
  enddo

end subroutine strncmp2b

integer function strncmp3(str1,str2)
  implicit none
  character(len=*), intent(in) :: str1, str2
!!$  integer, intent(in)       :: len_str

  integer :: ip,n,len1
  logical :: tf
  integer :: ndif, istart, i2,i, icharAcap,icharZcap,ichara,icharz,icharx &
       &   , icharupper, icharlower
!!$  character c
!  ----------
  len1 = len(str1)
  n = len(str2)
  icharAcap = ichar('A')
  icharZcap = ichar('Z')
  ichara    = ichar('a')
  icharz    = ichar('z')

  ndif = ichar('A') - ichar('a')

  tf = .false.
  strncmp3 = 0
  ip = 0
  do istart = 1, len1 - n + 1
     i2 = 0
     do i = istart, istart-1+n
        i2 = i2 + 1
        icharx = ichar(str2(i2:i2))
        if(icharAcap <= icharx .and. icharx <= icharZcap) then
           icharupper = icharx
           icharlower = icharx - ndif
        else if(ichara <= icharx .and. icharx <= icharz) then
           icharupper = icharx + ndif
           icharlower = icharx
        else
           icharupper = icharx; icharlower = icharx
        end if
!!$        c = char(ichar(str2(i2:i2)) + ndif)
!!$        if(str1(i:i).ne.str2(i2:i2).and.str1(i:i).ne.c) goto 1000
        if(str1(i:i).ne.char(icharupper).and.str1(i:i).ne.char(icharlower)) goto 1000
     enddo
     tf = .true.
     ip = istart
     strncmp3 = ip
1000 continue
     if(tf) return
  enddo
end function strncmp3

subroutine skip2next(str,len_str,charc,ic)
  implicit none
  character(len=*), intent(in) :: str
  character,     intent(in)  ::   charc
  integer,       intent(in)  ::   len_str
  integer,       intent(out) ::   ic

  logical ::          tf, flag
  character(len=1) :: c, c0
  integer ::          ndif, i, ipnt, ipnt2

  ndif = ichar('A') - ichar('a')
  c0   = charc
  if(ichar(c0) >= ichar('a') .and. ichar(c0) <= ichar('z')) c0 = char(ichar(c0) + ndif)

  flag = .false.

  do i = 1, len_str
     c = str(i:i)
     call isalph(c,tf)
     if(.not.flag.and.tf) then
        ipnt = i
        if(ichar(c) >= ichar('a') .and. ichar(c) <= ichar('z')) then
           c = char(ichar(c) + ndif)
        endif
        if(c == c0) goto 1001
     endif
     if(tf) then
        flag = .true.
     else
        flag = .false.
     endif
  enddo
  stop ' no matching character'

1001 continue
  do i = ipnt+1, len_str
     ipnt2 = i
     if(str(i:i) == '=') goto 1002
  enddo
  stop ' no = character'
  
1002 continue
  do i = ipnt2+1, len_str
     call isnumm(str(i:i),tf)
     if(tf) then
        ic = i
        return
     else
        call isalph(str(i:i),tf)
        if(tf) then
           ic = i
           return
        endif
     endif
  enddo
  stop ' no alphabet or number'
end subroutine skip2next

subroutine nameck(name,n)
  implicit none
  character(len=*), intent(inout) :: name
  integer, intent(in)          :: n

  integer i, icount

  icount = 0
  do  i = 1, n
     if(name(i:i).eq.' ') then
        icount = i
        goto 2
     endif
  enddo

2 continue
  if(icount.ne.0) then
     do  i = icount+1, n
        name(i:i) = ' '
     enddo
  endif
end subroutine nameck

subroutine ch4mix(str,len_str,waymix,istrbr,nbxmix,hownew &
     & ,cutoff_mix,decomp_into_charge_spin,work,nsize)
  use m_Const_Parameters, only : DP, SIMPLE, BROYD1, BROYD2, DFP, PULAY &
       &, ANEW, RENEW, SMALL, MEDIUM, LARGE
  implicit none
  character(len=*), intent(in)  :: str
  integer, intent(in)        :: len_str
  integer, intent(out)       :: waymix, istrbr, nbxmix, hownew
  integer, intent(out)       :: cutoff_mix
  logical, intent(out)       :: decomp_into_charge_spin
  real(DP), intent(out)      :: work(*)
  integer, intent(out)       :: nsize

  logical       :: tf
  integer       :: ic
  real(kind=DP) :: alpha

  call check_decomp

  call skip2alph(str,len_str,ic)
  if(ic <= 0) then
     call chnnm(str,len_str,5,work,nsize)
     if(nsize < 4) stop ' shortage of data for mixing ratios(2)'
     waymix = SIMPLE
     return
  end if

  call check_mix_tag_name(ic) ! -> waymix

  if(waymix /= SIMPLE) then
     call rd_istrbr   ! -> istrbr
     call rd_nbxmix   ! -> nbxmix
     call rd_hownew   ! -> hownew
     call rd_alpha    ! -> alpha
     work(1:5) = alpha
     nsize = 5
     call rd_cutoff   ! -> cutoff_mix
  end if

contains
  subroutine check_decomp
    logical :: tf

    call strncmp2(str,len_str,'nodecomp',8,tf)
    if(tf) then
       decomp_into_charge_spin = .false.
    else
       call strncmp2(str,len_str,'decomp',6,tf)
       if(tf) then
          decomp_into_charge_spin = .true.
       else
          decomp_into_charge_spin = .false.
       end if
    end if
  end subroutine check_decomp

  subroutine check_mix_tag_name(ic)
    integer, intent(in) :: ic
    logical :: tf
    call strncmp2(str(ic:len_str),6,'simple',6,tf)
    if(tf) then
       waymix = SIMPLE
       call chnnm2(str,len_str,5,work,nsize)
    end if

    if(.not.tf) then
       call strncmp2(str(ic:len_str),6,'broyd1',6,tf)
       if(tf) waymix = BROYD1
    end if

    if(.not.tf) then
       call strncmp2(str(ic:len_str),8,'broyden1',8,tf)
       if(tf) waymix = BROYD1
    end if

    if(.not.tf) then
       call strncmp2(str(ic:len_str),6,'broyd2',6,tf)
       if(tf) waymix = BROYD2
    end if

    if(.not.tf) then
       call strncmp2(str(ic:len_str),8,'broyden2',8,tf)
       if(tf) waymix = BROYD2
    end if

    if(.not.tf) then
       call strncmp2(str(ic:len_str),8,'dfp',3,tf)
       if(tf) waymix = DFP
    end if

    if(.not.tf) then
       call strncmp2(str(ic:len_str),7,'bluegel',7,tf)
       if(tf) waymix = DFP
    end if

    if(.not.tf) then
       call strncmp2(str(ic:len_str),2,'bl',2,tf)
       if(tf) waymix = DFP
    end if

    if(.not.tf) then
       call strncmp2(str(ic:len_str),5 ,'pulay',5,tf)
       if(tf) waymix = PULAY
    end if

    if(.not.tf)  stop ' ! no appropriate tag name'
  end subroutine check_mix_tag_name

  subroutine rd_istrbr
    integer :: ic
    call skip2next(str,len_str,'I',ic)
    read(str(ic:len_str),*) istrbr
    if(istrbr <= 0) then
       print *, ' istrbr is reset as to be 1, because istrbr(=',istrbr,') is invalid.'
       istrbr = 1
    endif
  end subroutine rd_istrbr

  subroutine rd_nbxmix
    integer :: ic
     call skip2next(str,len_str,'N',ic)
     read(str(ic:len_str),*) nbxmix
     if(nbxmix <= 0) then
        print *, ' nbxmix is reset as to be 5, because nbxmix(=',nbxmix,') is invalid.'
        nbxmix = 5
     endif
   end subroutine rd_nbxmix

   subroutine rd_hownew
     integer :: ic
     call skip2next(str,len_str,'H',ic)
     call strncmp2(str(ic:len_str),5 ,'renew',5,tf)
     if(tf) hownew = RENEW
     if(.not.tf) then
        call strncmp2(str(ic:len_str),4,'anew',4,tf)
        if(tf) then
           hownew = ANEW
        else
           print *,' hownew is set as to be RENEW, because hownew is not given properly '
           hownew = RENEW
        end if
     end if
   end subroutine rd_hownew

   subroutine rd_cutoff
     logical :: tf
     call strncmp2(str,len_str,'cutoff',6,tf)
     if(tf) call skip2next(str,len_str,'C',ic)
     call strncmp2(str(ic:len_str),5,'small',5,tf)
     if(tf) then
        cutoff_mix = SMALL
        return
     endif

     call strncmp2(str(ic:len_str),6,'medium',6,tf)
     if(tf) then
        cutoff_mix = MEDIUM
        return
     end if

     call strncmp2(str(ic:len_str),5,'large',5,tf)
     if(tf) then
        cutoff_mix = LARGE
        return
     end if
   end subroutine rd_cutoff

   subroutine rd_alpha
     integer :: ic
     call skip2next(str,len_str,'A',ic)
     read(str(ic:len_str),*) alpha
     if(alpha < 0.d0 .or. alpha > 1.0) then
        print *, ' alpha is set as to be 0.3,'&
        &       , ' because alpha(=',alpha,') is not given properly'
        alpha = 0.3d0
     end if
   end subroutine rd_alpha
end subroutine ch4mix

subroutine skip2alph(str,len_str,ic)
  implicit none
  character(len=*), intent(in) :: str
  integer, intent(in)       :: len_str
  integer, intent(out)      :: ic

  integer :: i
  logical :: tf
  ic = 0
  do i = 1, len_str
     ic = i
     call isalph(str(i:i),tf)
     if(str(i:i).eq.':') goto 999
     if(tf) return
  enddo
999 ic = 0
  return
end subroutine skip2alph

subroutine chnnm2(str,len_str,n,wka,incunt)
  use m_Const_Parameters, only : DP
  implicit none
  character(len=*), intent(in) :: str
  integer, intent(in)       :: len_str,n
  real(DP), intent(out),dimension(n) :: wka
  integer, intent(out)      :: incunt

  logical  ::         tf
  character(len=1) :: c
  logical  ::         is_chr
  real(kind=DP) ::    xnum1
  integer  ::         icount, i, ic2

  xnum1 = 0.d0

  wka = 0.d0

  incunt = 0
  icount = 0
  is_chr = .false.
1001 continue
  do i = icount+1, len_str
     icount = i
     c = str(i:i)
     if(c == ':') goto 1006
     call isnumm(c,tf)
     if(.not.is_chr.and.tf) then
        read(str(i:len_str),*) xnum1
        goto 1002
     endif
     if(c == '='.or.c == ' ') then
        is_chr = .false.
     else
        is_chr = .true.
     endif
  end do
  goto 1006
1002 continue
  incunt = incunt + 1
  wka(incunt) = xnum1

  if(incunt == n) goto 1006

  ic2 = icount + 1
  do i = icount+1, len_str
     if(str(i:i) == ' '.or.str(i:i) == ',') then
        ic2 = i
        goto 1004
     else if(str(i:i) == ':') then
        goto 1006
     endif
  end do
1004 continue
  icount = ic2
  is_chr = .false.
  goto 1001

1006 continue
end subroutine chnnm2

subroutine ch4prc(str,len_str,c_precon,amix,bmix,precon_charge_only &
     &,isterm)
  use m_Const_Parameters, only : DP
  implicit none
  character(len=*),intent(in) :: str
  integer,intent(in) ::          len_str
  logical,intent(out) ::         c_precon, precon_charge_only
  real(DP), intent(out) ::       amix, bmix
  logical,intent(out) ::         isterm

  logical  :: tf
  integer  :: i, ic, icc

  c_precon = .false.
  amix = -1.d0
  bmix = -1.d0
  precon_charge_only = .false.
  isterm = .false.

  call skip2alph(str,len_str,ic)
  if(ic <= 0) return

  isterm = .true.
  call strncmp2(str(ic:len_str),2,'cp',2,tf)
  if(tf) then
     c_precon = .true.
  else
     c_precon = .false.
  endif

  if(c_precon) then
     call skip2next(str,len_str,'A',ic)
     read(str(ic:len_str),*) amix
     if(amix > 1.d0) then
        print *,' amix is reset as to be 1.0, because amix (=',amix,') is too large.'
        amix = 1.d0
     endif

     call skip2next(str,len_str,'B',ic)
     read(str(ic:len_str),*) bmix

     do i = ic, len_str
        icc = i
        call isalph(str(i:i),tf)
        if(tf) goto 2000
        if(str(i:i).eq.':') return
     enddo

2000 continue
     ic = icc

     write(6,*) ' !! precon_charge_only = ', precon_charge_only
     if(str(ic:ic).eq.'p'.or.str(ic:ic).eq.'P') then
        call f_strncmp(str(ic:len_str),'precon_charge_only',17,tf)
        if(tf) precon_charge_only = .true.
     endif
  endif
end subroutine ch4prc

subroutine f_strncmp(str1,str2,n,tf)
  implicit none
  character(len=*), intent(in) :: str1, str2
  integer, intent(in) ::          n
  logical, intent(out) ::         tf

  character :: c
  integer   :: ndif, i
  ndif = ichar('A') - ichar('a')

  tf = .false.
  do i = 1, n
     c = char(ichar(str2(i:i)) + ndif)
     if(str1(i:i)/=str2(i:i).and.str1(i:i)/=c) goto 1000
  enddo
  tf = .true.
1000 continue
end subroutine f_strncmp

subroutine rewind_to_tag0(nfinp, len_tag, tag, eof_reach, tag_is_found&
     &, str, len_str)
!                           @(#)b_Words.f90 1.3 98/12/24 02:06:35
  implicit none
  integer, intent(in) ::       nfinp
  integer, intent(in) ::       len_tag
  character(len=*), intent(in) :: tag
  logical, intent(out) ::      eof_reach, tag_is_found
  integer, intent(in) ::       len_str
  character(len=*) ::          str

  tag_is_found = .false.
  rewind(nfinp)
  do while(.not.tag_is_found)
     read(nfinp,'(a132)',end=1001) str
!!$     call strncmp2(str,len_str,tag,len_tag,tag_is_found)
     call strncmp2b(str,len_str,tag,len_tag,tag_is_found)
  enddo
  return
1001 tag_is_found = .false.
  eof_reach = .true.
  return
end subroutine rewind_to_tag0

subroutine read_RHS_number(str,len_str,values,n)
  use m_Const_Parameters, only : DP
  implicit none
  integer, intent(in)              :: len_str,n
  character(len=len_str), intent(out) :: str
  real(kind=DP),intent(out)        :: values(n)

  integer ic, incunt

  call skip2RHS(str,len_str,ic)
  call chnnm2(str(ic:len_str),len_str-ic+1,n,values,incunt)
  if(incunt < n) print *,' !D incunt(=',incunt,' < n(= ',n,')'
end subroutine read_RHS_number

subroutine read_LHS_number(str,len_str,values,n)
  use m_Const_Parameters, only : DP
  implicit none
  integer, intent(in)              :: len_str,n
  character(len=len_str), intent(out) :: str
  real(kind=DP),intent(out)        :: values(n)

  integer   :: ic, incunt

  call skip2number(str,len_str,ic)
  call chnnm(str(ic:len_str),len_str-ic+1,n,values,incunt)
  if(incunt < n) print *, ' !D incunt(=',incunt,' < n(= ',n,')'
end subroutine read_LHS_number

subroutine skip2RHS(str,len_str,ic)
  use m_Const_Parameters, only : DP
  implicit none
  integer, intent(in)       :: len_str
  integer, intent(out)      :: ic
  character(len=*), intent(in) :: str

  logical     ::  tf
  character(len=1) :: charc, c
  data charc/'='/
  integer     :: i, ipnt, ipnt2

  do i = 1, len_str
     c = str(i:i)
     if(c == charc) then
        ipnt = i
        tf = .true.
        goto 1001
     endif
  end do
  stop ' no matching character'

1001 continue

  ipnt2 = ipnt
  do i = ipnt2+1, len_str
     call isnumm(str(i:i),tf)
     if(tf) then
        ic = i
        return
     else
        call isalph(str(i:i),tf)
        if(tf) then
           ic = i
           return
        endif
     endif
  end do
  stop ' no alphabet or number'
end subroutine skip2RHS

subroutine f_strcpy(strin,len_strin,strout,len_strout)
  implicit none
  integer, intent(in)      :: len_strin, len_strout
  character(len=len_strin),  intent(in)  :: strin
  character(len=len_strout), intent(out) :: strout

  integer  :: i, len_cp

  if(len_strout < len_strin) then
     len_cp = len_strout
  else
     len_cp = len_strin
  endif

  do i = len_cp + 1, len_strout
     strout(i:i) = ' '
  enddo

  do i = 1, len_cp
     strout(i:i) = strin(i:i)
  enddo

end subroutine f_strcpy

subroutine skip2number(str,len_str,ic)
  implicit none
  integer,         intent(in)  ::    len_str
  integer,         intent(out) ::    ic
  character(len=len_str), intent(in)  :: str

  integer :: i
  logical :: tf
  ic = 0
  do i = 1, len_str
     ic = i
     call isnumm(str(i:i),tf)
     if(str(i:i) == ':') goto 999
     if(tf) return
  enddo
999 ic = 0
end subroutine skip2number

