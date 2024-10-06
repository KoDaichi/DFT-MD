!PHASEの関数、一部を抜粋
!Copyright (C) 2002-2014 NIMS
!Copyright (C) 2014 Hiroyuki Oshima

module m_phase_commons
  implicit none
#define SCHOENFLIES(MAIN,SUB) MAIN//SUB
contains
!m_CS_SpaceGroup::get_point_group_name
  subroutine get_point_group_name(nfout,system,nsym,nopr,pg_name,pg_name_i)
    integer, intent(in) :: nfout
    character(len=9), intent(in) :: system
    integer, intent(in) :: nsym
    integer, intent(in) :: nopr(nsym)
    character(*), intent(out) :: pg_name
    character(*), intent(out) :: pg_name_i

    logical :: fexist(48)
    integer :: g,i

    fexist(1:48) = .false.
    do i=1,nsym
       fexist(nopr(i)) = .true.
    end do

    pg_name = SCHOENFLIES('C','1')
    pg_name_i = '1'
    g = nsym
    if(system(1:1) == 'c' .or. system(1:1) == 'C') then
       if(fexist(25)) then
! IE
          g = g/2
       end if
       if(g==1) then
          if(fexist(25)) then
! IE
100          pg_name = SCHOENFLIES('C','i')
             pg_name_i = '-1'
          else
             pg_name = SCHOENFLIES('C','1')
             pg_name_i = '1'
          end if
       else if(g==2) then
          if(fexist(2).or.fexist(3).or.fexist(4).or. &
           & fexist(13).or.fexist(14).or.fexist(15).or. &
           & fexist(16).or.fexist(17).or.fexist(18)) then
! C2x or C2y or C2z or
! C2a or C2b or C2c or
! C2d or C2e or C2f
             if(fexist(25)) then
! IE
                pg_name = SCHOENFLIES('C','2h')
                pg_name_i = '2/m'
             else
                pg_name = SCHOENFLIES('C','2')
                pg_name_i = '2'
             end if
          else if(fexist(26).or.fexist(27).or.fexist(28).or. &
                & fexist(37).or.fexist(38).or.fexist(39).or. &
                & fexist(40).or.fexist(41).or.fexist(42)) then
! IC2x or IC2y or IC2z or
! IC2a or IC2b or IC2c or
! IC2d or IC2e or IC2f
             pg_name = SCHOENFLIES('C','s')
             pg_name_i = 'm'
          end if
       else if(g==3) then
          if(fexist(25)) then
! IE
             pg_name = SCHOENFLIES('S','6')
             pg_name_i = '-3'
          else
             pg_name = SCHOENFLIES('C','3')
             pg_name_i = '3'
          end if
       else if(g==4) then
          if((fexist(2).and.fexist(3).and.fexist(4)).or.&
           & (fexist(2).and.fexist(16).and.fexist(18)).or.&
           & (fexist(3).and.fexist(15).and.fexist(17)).or.&
           & (fexist(4).and.fexist(13).and.fexist(14))) then
! C2x and C2y and C2z or
! C2x and C2d and C2f or
! C2y and C2c and C2e or
! C2z and C2a and C2b
             if(fexist(25)) then
! IE
                pg_name = SCHOENFLIES('D','2h')
                pg_name_i ='mmm'
             else
                pg_name =SCHOENFLIES('D','2')
                pg_name_i ='222'
             end if
          else if(fexist(26).or.fexist(27).or.fexist(28).or. &
                & fexist(37).or.fexist(38).or.fexist(39).or. &
                & fexist(40).or.fexist(41).or.fexist(42).and.&
                & .not.(fexist(19).or.fexist(20).or.fexist(21))) then
! IC2x or IC2y or IC2z or
! IC2a or IC2b or IC2c or
! IC2d or IC2e or IC2f
             pg_name =SCHOENFLIES('C','2v')
             pg_name_i ='mm2'
          else if((fexist(19).and.fexist(22)).or. &
                & (fexist(20).and.fexist(23)).or. &
                & (fexist(21).and.fexist(24))) then
! (C4x+ and C4x-) or (C4y+ and C4y-) or (C4z+ and C4z-)
             if(fexist(25)) then
! IE
                pg_name =SCHOENFLIES('C','4h')
                pg_name_i ='4/m'
             else
                pg_name =SCHOENFLIES('C','4')
                pg_name_i ='4'
             end if
          else if((fexist(43).and.fexist(46)).or. &
                & (fexist(44).and.fexist(47)).or. &
                & (fexist(45).and.fexist(48))) then
! (IC4x+ and IC4x-) or (IC4y+ and IC4y-) or (IC4z+ and IC4z-)
             pg_name =SCHOENFLIES('S','4')
             pg_name_i ='-4'
          end if
       else if(g==6) then
          if((fexist(5).and.fexist(14).and.fexist(17).and.fexist(18)) .or. &
           & (fexist(6).and.fexist(14).and.fexist(15).and.fexist(16)) .or. &
           & (fexist(7).and.fexist(13).and.fexist(15).and.fexist(18)) .or. &
           & (fexist(8).and.fexist(13).and.fexist(16).and.fexist(17))) then
! C31+ and C2b and C2e and C2f
! C32+ and C2b and C2c and C2d
! C33+ and C2a and C2c and C2f
! C34+ and C2a and C2d and C2e
             if(fexist(25)) then
! IE
                pg_name =SCHOENFLIES('D','3d')
                pg_name_i ='-3m'
             else
                pg_name =SCHOENFLIES('D','3')
                pg_name_i ='32'
             end if
          else
             pg_name =SCHOENFLIES('C','3v')
             pg_name_i ='3m'
          end if
       else if(g==8) then
          if((fexist(19).and.fexist(22)).or. &
           & (fexist(20).and.fexist(23)).or. &
           & (fexist(21).and.fexist(24))) then
! (C4x+ and C4x-) or (C4y+ and C4y-) or (C4z+ and C4z-)
             if(fexist(37).or.fexist(38).or.fexist(39).or. &
              & fexist(40).or.fexist(41).or.fexist(42)) then
! IC2a or IC2b or IC2c or IC2d or IC2e or IC2f
                if(fexist(25)) then
! IE
                   pg_name =SCHOENFLIES('D','4h')
                   pg_name_i ='4/mmm'
                else
                   pg_name =SCHOENFLIES('C','4v')
                   pg_name_i ='4mm'
                end if
             else
                pg_name =SCHOENFLIES('D','4')
                pg_name_i ='422'
             end if
          else
             pg_name =SCHOENFLIES('D','2d')
             pg_name_i ='-42m'
          end if
       else if(g==12) then
          if(fexist(25)) then
! IE
             pg_name ='Th'
             pg_name_i ='m-3'
          else
             pg_name ='T'
             pg_name_i ='23'
          end if
       else if(g==24) then
          if(fexist(19).and.fexist(22).and. &
           & fexist(20).and.fexist(23).and. &
           & fexist(21).and.fexist(24)) then
! C4x+ and C4x- and C4y+ and C4y- and C4z+ and C4z-
             if(fexist(25)) then
! IE
                pg_name =SCHOENFLIES('O','h')
                pg_name_i ='m-3m'
             else
                pg_name ='O'
                pg_name_i ='432'
             end if
          else
             pg_name ='Td'
             pg_name_i ='-43m'
          end if
       end if
    else if(system(1:1) == 'h' .or. system(1:1) == 'H') then! hexagonal
       if(fexist(13)) then
! IE
          g = g/2
          if(g==1) go to 100 !add By Tomita
       else if(fexist(16)) then
! IC2
          g = g/2
          if(g==1) go to 200 !add by Tomita
       end if
       if(g==2) then
           if(fexist(7).or.fexist(8).or.fexist(9).or. &
           & fexist(10).or.fexist(11).or.fexist(12).or.fexist(4)) then
! C211 or C221 or C231 or C212 or C222 or C232 or C2
             if(fexist(13)) then
! IE
                pg_name =SCHOENFLIES('C','2h')
                pg_name_i ='2/m'
             else
! add by oshima 2013 10 20              
                IF(FEXIST(19).OR.FEXIST(20).OR.FEXIST(21) &
                     .OR.FEXIST(22).OR.FEXIST(23).OR.FEXIST(24) &
                     )THEN
                   pg_name = SCHOENFLIES('C','2v')  
                   pg_name_i = 'mm2'
                else
!ccccccccccccccccccccccccccc
                   pg_name =SCHOENFLIES('C','2')
                   pg_name_i ='2'
                end if
             end if
          else if(fexist(19).or.fexist(20).or.fexist(21).or. &
                & fexist(22).or.fexist(23).or.fexist(24).or.fexist(16)) then
! IC211 or IC221 or IC231 or IC212 or IC222 or IC232 or IC2
200          pg_name =SCHOENFLIES('C','s')
             pg_name_i ='m'
          end if
       else if(g==3) then
          if(fexist(13)) then
! IE
             pg_name =SCHOENFLIES('S','6')
             pg_name_i ='-3'
          else if(fexist(16)) then
! IC2
             pg_name =SCHOENFLIES('C','3h')
             pg_name_i ='-6'
          else
             pg_name =SCHOENFLIES('C','3')
             pg_name_i ='3'
          end if
       else if(g==4) then
          if(fexist(4).and.( &
            & (fexist(7).and.fexist(10)).or. &
            & (fexist(8).and.fexist(11)).or. &
            & (fexist(9).and.fexist(12)))) then
! C2 and ( C211 and C212 ) or ( C221 and C222 ) or ( C231 and C232 ) 
             if(fexist(13)) then
! IE
                pg_name =SCHOENFLIES('D','2h')
                pg_name_i ='mmm'
             else
                pg_name =SCHOENFLIES('D','2')
                pg_name_i ='222'
             end if
          else if( (fexist(4).and.( &
            & (fexist(19).and.fexist(22)).or. &
            & (fexist(20).and.fexist(23)).or. &
            & (fexist(21).and.fexist(24))) ) .or. &
            & (fexist(16).and.( &
            & (fexist(7).and.fexist(22)).or. &
            & (fexist(8).and.fexist(23)).or. &
            & (fexist(9).and.fexist(24)).or. &
            & (fexist(19).and.fexist(10)).or. &
            & (fexist(20).and.fexist(11)).or. &
            & (fexist(21).and.fexist(12))) ) ) then
! C2 and (IC211 and IC212) or (IC221 and IC222) or (IC231 and IC232) 
! or IC2 and (C211 and IC212) or (C221 and IC222) or (C231 and IC232)
!         or (IC211 and C212) or (IC221 and C222) or (IC231 and C232)
             pg_name =SCHOENFLIES('C','2v')
             pg_name_i ='mm2'
          end if
       else if(g==6) then
          if(fexist(2).and.fexist(6)) then
! C6+ and C6-
             if(fexist(13)) then
! IE
                pg_name =SCHOENFLIES('C','6h')
                pg_name_i ='6/m'
             else
                pg_name =SCHOENFLIES('C','6')
                pg_name_i ='6'
             end if
          else if((fexist(7).and.fexist(8).and.fexist(9)) .or. &
                & (fexist(10).and.fexist(11).and.fexist(12))) then
! (C211 and C221 and C231) or (C212 and C222 and C232)
             if(fexist(13)) then
! IE
                pg_name =SCHOENFLIES('D','3d')
                pg_name_i ='-3m'
             else if(fexist(16)) then
! IC2
                pg_name =SCHOENFLIES('D','3h')
                pg_name_i ='-6m2'
             else
                pg_name =SCHOENFLIES('D','3')
                pg_name_i ='32'
             end if
          else if((fexist(19).and.fexist(20).and.fexist(21)) .or. &
                & (fexist(22).and.fexist(23).and.fexist(24))) then
! (IC211 and IC221 and IC231) or (IC212 and IC222 and IC232)
             pg_name =SCHOENFLIES('C','3v')
             pg_name_i ='3m'
          end if
       else if(g==12) then
          if(fexist(7).and.fexist(8).and.fexist(9).and. &
           & fexist(10).and.fexist(11).and.fexist(12)) then
! C211 and C221 and C231 and C212 and C222 and C232
             if(fexist(13)) then
! IE
                pg_name =SCHOENFLIES('D','6h')
                pg_name_i ='6/mmm'
             else
                pg_name =SCHOENFLIES('D','6')
                pg_name_i ='622'
             end if
          else if(fexist(19).and.fexist(20).and.fexist(21).and. &
                & fexist(22).and.fexist(23).and.fexist(24)) then
! IC211 and IC221 and IC231 and IC212 and IC222 and IC232
             pg_name =SCHOENFLIES('C','6v')
             pg_name_i ='6mm'
          end if
       end if
    else
       go to 20
    end if
    if(g>1 .and. trim(pg_name) == SCHOENFLIES('C','1')) then
20     write(nfout,*)'get_point_group_name: error'
       write(nfout,*)'system=',system
       write(nfout,*)'nsym=',nsym
       write(nfout,*)'g=',g
       write(nfout,*)'fexist=',fexist
!      stop 'get_point_group_name: error'
    end if
  end subroutine get_point_group_name
end module m_phase_commons
