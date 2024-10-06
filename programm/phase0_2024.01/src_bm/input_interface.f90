!=======================================================================
!
!  PROGRAM  PHASE/0 2014.01 (rev.375)
!
!  "First-principles Electronic Structure Calculation Program"
!
!  FUNCTION:  getUnitId, setUnit, updateUnits, clearUnitFlag, setDefaultUnits,
!            realConvByUnit, f_openInputFile, f_closeInputFile, f_selectTop,
!        f_selectBlock, f_selectParentBlock, f_getIntValue, f_getIntVectorValue
!          f_getRealValue, f_getRealVectorValue, f_getStringValue,
!          f_selectFirstTableLine, f_selectNextTableLine, f_readUnitCell,
!          f_readSymmetryGenerator, f_readKPoints, fstreq
!
!  AUTHOR(S): K. Mae, T. Yamasaki   August/20/2003
!  
!  Contact address :  Phase System Consortium
!                     E-mail: phase_system@nims.go.jp URL https://azuma.nims.go.jp
!  
!
!
!=======================================================================
! $Id: input_interface.F90 238 2012-11-12 04:11:13Z yamasaki $
!
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$! interface functions
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$module input_interface
!!$
!!$        use input_tag
!!$        use unit_conv
!!$
!!$	integer, parameter :: FMAXTAGLEN = 256
!!$	integer, parameter :: FMAXVALLEN = 256
!!$	integer, parameter :: FMAXUNITLEN = 16
!!$	real(8), parameter :: PI = 3.14159265358979323846d0
!!$
!!$        type unitsystem
!!$            integer :: energy
!!$            integer :: length
!!$            integer :: time
!!$            integer :: mass
!!$            integer :: angle
!!$            integer :: temperature
!!$            integer :: velocity
!!$            integer :: force
!!$            integer :: pressure
!!$        end type unitsystem
!!$        type(unitsystem) :: usys_file, usys_ret
!!$
!!$        private getUnitNo
!!$        !interface
!!$        !    integer function getUnitNo( usys, unittypeid )
!!$        !        type(unitsystem), intent(in) :: usys
!!$        !        integer, intent(in) :: unittypeid
!!$        !    end function getUnitNo
!!$        !end interface
!!$contains


integer function getUnitId( usys, unittype )
    use m_Const_Parameters
    implicit none
    type(unitsystem), intent(in) :: usys
    integer, intent(in) :: unittype

    select case ( unittype )
    case ( TYPE_ENERGY )
        getUnitId = usys%energy
    case ( TYPE_LENGTH )
        getUnitId = usys%length
    case ( TYPE_TIME )
        getUnitId = usys%time
    case ( TYPE_MASS )
        getUnitId = usys%mass
    case ( TYPE_ANGLE )
        getUnitId = usys%angle
    case ( TYPE_TEMPERATURE )
        getUnitId = usys%temperature
    case ( TYPE_VELOCITY )
        getUnitId = usys%velocity
    case ( TYPE_FORCE )
        getUnitId = usys%force
    case ( TYPE_PRESSURE )
        getUnitId = usys%pressure
    end select

    return
end function getUnitId

integer function setUnit( usys, unitid, unittype )
    use m_Const_Parameters
    implicit none
    type(unitsystem), intent(inout) :: usys
    integer, intent(in) :: unitid, unittype

    select case ( unittype )
    case ( TYPE_ENERGY )
        if( usys%inpfg_energy ) then
            setUnit = -1
            return           
        end if
        usys%energy = unitid
        usys%inpfg_energy = .true.
    case ( TYPE_LENGTH )
        if( usys%inpfg_length ) then
            setUnit = -1
            return           
        end if
        usys%length = unitid
        usys%inpfg_length = .true.
    case ( TYPE_TIME )
        if( usys%inpfg_time ) then
            setUnit = -1
            return           
        end if
        usys%time = unitid
        usys%inpfg_time = .true.
    case ( TYPE_MASS )
        if( usys%inpfg_mass ) then
            setUnit = -1
            return           
        end if
        usys%mass = unitid
        usys%inpfg_mass = .true.
    case ( TYPE_ANGLE )
        if( usys%inpfg_angle ) then
            setUnit = -1
            return           
        end if
        usys%angle = unitid
        usys%inpfg_angle = .true.
    case ( TYPE_TEMPERATURE )
        if( usys%inpfg_temperature ) then
            setUnit = -1
            return           
        end if
        usys%temperature = unitid
        usys%inpfg_temperature = .true.
    case ( TYPE_VELOCITY )
        if( usys%inpfg_velocity ) then
            setUnit = -1
            return           
        end if
        usys%velocity = unitid
        usys%inpfg_velocity = .true.
    case ( TYPE_FORCE )
        if( usys%inpfg_force ) then
            setUnit = -1
            return           
        end if
        usys%force = unitid
        usys%inpfg_force = .true.
    case ( TYPE_PRESSURE )
        if( usys%inpfg_pressure ) then
            setUnit = -1
            return           
        end if
        usys%pressure = unitid
        usys%inpfg_pressure = .true.
    end select

    setUnit = 0
    return
end function setUnit

integer function updateUnits( usys )
    use m_Const_Parameters
    implicit none
    type(unitsystem), intent(inout) :: usys
    character(FMAXUNITLEN) :: unit1, unit2, unit3
    integer unitid, unittype, iret

    !!! velocity unit
    if( .not. usys%inpfg_velocity ) then
        unit1 = unit_list(usys%length)%name
        unit2 = unit_list(usys%time)%name
        unit3 = trim(unit1)//'/'//trim(unit2)
        iret = get_unit_id( unit3, unitid, unittype )
        if( iret == 0 .and. unittype == TYPE_VELOCITY ) then
            usys%velocity = unitid
        end if
    end if

    !!! force unit
    if( .not. usys%inpfg_force ) then
        unit1 = unit_list(usys%energy)%name
        unit2 = unit_list(usys%length)%name
        if( unit1 == 'j' .and. unit2 == 'm' ) then
            unit3 = 'n'
        else
            unit3 = trim(unit1)//'/'//trim(unit2)
        end if
        iret = get_unit_id( unit3, unitid, unittype )
        if( iret == 0 .and. unittype == TYPE_FORCE ) then
            usys%force = unitid
        end if
    end if

    !!! force pressure
    if( .not. usys%inpfg_pressure ) then
        unit1 = unit_list(usys%energy)%name
        unit2 = unit_list(usys%length)%name
        if( unit1 == 'j' .and. unit2 == 'm' ) then
            unit3 = 'pa'
        else
            unit3 = trim(unit1)//'/'//trim(unit2)//'3'
        end if
        iret = get_unit_id( unit3, unitid, unittype )
        if( iret == 0 .and. unittype == TYPE_PRESSURE ) then
            usys%pressure = unitid
        end if
    end if

    updateUnits = 0
    return
end function updateUnits


integer function clearUnitFlag( usys )
    use m_Const_Parameters
    implicit none
    type(unitsystem), intent(out) :: usys

    usys%inpfg_energy = .false.
    usys%inpfg_length = .false.
    usys%inpfg_time = .false.
    usys%inpfg_mass = .false.
    usys%inpfg_angle = .false.
    usys%inpfg_temperature = .false.
    usys%inpfg_velocity = .false.
    usys%inpfg_force = .false.
    usys%inpfg_pressure = .false.

    clearUnitFlag = 0
    return
end function clearUnitFlag

integer function setDefaultUnits()
    use m_Const_Parameters
    implicit none
    integer unitid, unittype, i, iret
    integer setUnit, clearUnitFlag

    do i = 1, numdefaultunits
        iret = get_unit_id( defaultunits(i), unitid, unittype )
        iret = setUnit( usys_file(1), unitid, unittype )
    end do
    iret = clearUnitFlag( usys_file(1) )
    setDefaultUnits = iret
    return
end function setDefaultUnits


integer function realConvByUnit( val_f, val_r, unit_f, unit_r )
        use m_Const_Parameters
	implicit none
	real(8), intent(in)     :: val_f
	real(8), intent(out)     :: val_r
	character(*), intent(in) :: unit_f
	character(*), intent(in) :: unit_r
        real(8) wkret
	character(FMAXUNITLEN) wk_unit_f, wk_unit_r
	integer unitid_r, unittype_r, unitid_f, unittype_f, iret
        integer getUnitId

	wk_unit_f = unit_f
	wk_unit_r = unit_r
	call f_strlower( wk_unit_f )
	call f_strlower( wk_unit_r )
	!!! print '( a, "=>", a )', unit_f, wk_unit_f
	!!!print '( a, "=>", a )', unit_r, wk_unit_r
	
        if( unit_r == '' ) then
            val_r = val_f
            realConvByUnit = 0
        else
	    iret = get_unit_id( wk_unit_r, unitid_r, unittype_r ) 
 !!!debug print '( "get_unit_id,1: wk_unit_r=[", a, "] unitid_r=",I4,"unittype_r=",I4 )', wk_unit_r, unitid_r, unittype_r
            if( iret < 0 ) then
	       print '( "The unit [", a, "] in the program is not in the unit list." )', trim(unit_r)
	       realConvByUnit = -1
	       return
	    end if    
	    if( wk_unit_f == '' ) then
	        unittype_f = unittype_r
	        unitid_f = getUnitId( usys_file(fblkdepth), unittype_f )
!!!debug print '( "getUnitId,1: unitid_f=", I4, " unittype_f=",I4," fblkdepth=",I4 )', unitid_f, unittype_f, fblkdepth
 	    else
	        iret = get_unit_id( wk_unit_f, unitid_f, unittype_f )
	        if( iret < 0 ) then
		    print '( "The unit [", a, "] in the file is not in the unit list." )', trim(unit_f)
		    realConvByUnit = -1
		    return    
		else if( unittype_r /= unittype_f ) then
		    print *, "unit category mismatch"
		    realConvByUnit = -2
		    return
		end if 
	    end if
!!!debug print '( "unit_conv_byid,1: unitid_f=", I4, " unitid_r=",I4 )', unitid_f, unitid_r
            realConvByUnit = unit_conv_byid( val_f, val_r, unitid_f, unitid_r )
        end if 
	return
end function realConvByUnit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! open and read input file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer function f_openInputFile( fname )

	implicit none
	character(*), intent(in) :: fname
	integer openInputFile, setDefaultUnits, iret
	f_openInputFile = openInputFile( trim(fname)//char(0) )
        iret = setDefaultUnits()
	return
end function f_openInputFile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! dealloc read buffer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer function f_closeInputFile()

	implicit none
	integer closeInputFile

	f_closeInputFile = closeInputFile()
	return;
end function f_closeInputFile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer function f_selectTop()
  use m_Const_Parameters, only : fblkdepth
	implicit none
	integer selectTop

	f_selectTop = selectTop()
	if( f_selectTop >= 0 ) then
	    fblkdepth = 1
	end if
	return
end function f_selectTop


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer function f_selectBlock( blocktag )
  use m_Const_Parameters
  
	implicit none
	character(*), intent(in) :: blocktag
	integer numunits, ret, i, unittype, unitid
	character(FMAXUNITLEN) :: unit
	integer selectBlock, getNumBlockUnits, getBlockUnit, setUnit, iret 
	integer clearUnitFlag, updateUnits, f_selectParentBlock 

	f_selectBlock = selectBlock( trim(blocktag)//char(0) )
	if( f_selectBlock /= 0 ) then
	    return
	end if

        fblkdepth = fblkdepth + 1
        usys_file(fblkdepth) = usys_file(fblkdepth-1)
        f_selectBlock = clearUnitFlag( usys_file(fblkdepth) ) 
	f_selectBlock = getNumBlockUnits( numunits )
	do i = 0, numunits-1
	    iret = getBlockUnit( unit, i )
            iret = get_unit_id( unit, unitid, unittype )
	    if( iret < 0 ) then
	        f_selectBlock  = iret
                iret = f_selectParentBlock()
                print '( "The unit [", a, "] in the block [", a, "]is not found in the unit list.")', trim(unit), trim(blocktag)
 		return
	    end if
	    iret = setUnit( usys_file(fblkdepth), unitid, unittype )
	    if( iret < 0 ) then
	        f_selectBlock  = iret
                iret = f_selectParentBlock()
                print '( "A unit with the same dimension as [", a, "] has been already given in the block [", a, "]." )', &
                trim(unit), trim(blocktag)
 		return
	    end if
	end do
        iret = updateUnits( usys_file(fblkdepth) )
	return
end function f_selectBlock

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer function f_selectParentBlock()
  use m_Const_Parameters, only : fblkdepth
	implicit none
	integer selectParentBlock

	f_selectParentBlock = selectParentBlock()
	if( f_selectParentBlock /= 0 ) then
	    return
	end if
        fblkdepth = fblkdepth - 1
	return;
end function f_selectParentBlock

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  get integer value at tag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer function f_getIntValue( tag, ret )

	implicit none
	character(*), intent(in) :: tag
	integer, intent(out)     :: ret
	integer getIntValue

	f_getIntValue = getIntValue( trim(tag)//char(0), ret )
	return
end function f_getIntValue

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  get integer vector value at tag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer function f_getIntVectorValue( tag, ret )

	implicit none
	character(*), intent(in) :: tag
	integer, dimension(:), intent(out)     :: ret
	integer getIntVectorValue

	f_getIntVectorValue = getIntVectorValue( trim(tag)//char(0), ret )
	return
end function f_getIntVectorValue

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  get real value at tag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer function f_getRealValue( tag, ret, unit_r )
  use m_Const_Parameters
	implicit none
	character(*), intent(in) :: tag
	real(8), intent(out)     :: ret
	character(*), intent(in) :: unit_r
        real(8) wkret
	character(FMAXUNITLEN) unit_f
	integer getRealValue, realConvByUnit, iret_unit, iret_val

        unit_f = ''
	iret_val = getRealValue( trim(tag)//char(0), wkret, unit_f )
        if( iret_val /= 0 ) then
           f_getRealValue = iret_val
           return
        end if

  	iret_unit = realConvByUnit( wkret, ret, unit_f, unit_r )
        f_getRealValue = iret_unit !!!iret_val

	return
end function f_getRealValue

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get real vector value at tag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer function f_getRealVectorValue( tag, ret, unit_r )
  use m_Const_Parameters
	implicit none
	character(*), intent(in)           :: tag
	real(8), dimension(:), intent(out) :: ret
	character(*), intent(in) :: unit_r
        real(8), dimension(3) :: wkret
	character(FMAXUNITLEN) unit_f
	integer i
	integer getRealVectorValue, realConvByUnit, iret_unit, iret_val

	iret_val = getRealVectorValue( trim(tag)//char(0), wkret, unit_f )
        if( iret_val /= 0 ) then
           f_getRealVectorValue = iret_val
            return
        end if

	do i = 1, 3
            iret_unit = realConvByUnit( wkret(i), ret(i), unit_f, unit_r )
	    if( iret_unit < 0 ) then
                f_getRealVectorValue = iret_unit
       	        return;
	    end if 
	end do
        f_getRealVectorValue = iret_val
	return
end function f_getRealVectorValue

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get string value at tag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer function f_getStringValue( tag, ret, convfg )

	implicit none
	character(*), intent(in)  :: tag
	character(*), intent(out) :: ret
	integer, intent(in)  :: convfg
	integer getStringValue

	!!ret = "\0"
	f_getStringValue = getStringValue( trim(tag)//char(0), ret, convfg )

	!! print *,"in f_getStringValue ret=", ret
	return
end function f_getStringValue

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! seek the first line of a table
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer function f_selectFirstTableLine()

	implicit none
	integer selectFirstTableLine

	f_selectFirstTableLine = selectFirstTableLine()
	return
end function f_selectFirstTableLine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! seek the next line in a table
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer function f_selectNextTableLine()

	implicit none
	integer selectNextTableLine

	f_selectNextTableLine = selectNextTableLine()
	return
end function f_selectNextTableLine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! intelligent reading of unit cell
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$integer function f_readUnitCell(avec,bvec,cvec,a,b,c,alpha,beta,gamma,unit &
!!$     & , tag_avector,tag_bvector,tag_cvector,tag_a,tag_b,tag_c,tag_alpha,tag_beta,tag_gamma)
integer function f_readUnitCell(avec,bvec,cvec,a,b,c,ca,cb,cc,ucinptype,unit &
     & , tag_avector,tag_bvector,tag_cvector,tag_a,tag_b,tag_c,tag_alpha,tag_beta,tag_gamma,printable)
  use m_Const_Parameters, only : DP, FMAXUNITLEN,PAI

	implicit none
	real(DP), dimension(3), intent(inout) :: avec
	real(DP), dimension(3), intent(out) :: bvec
	real(DP), dimension(3), intent(out) :: cvec
	real(DP), intent(out) :: a, b, c
	real(DP), intent(out) :: ca, cb, cc
        integer, intent(out) :: ucinptype
        character(*), intent(in) :: unit
        character(*), intent(in) :: tag_avector,tag_bvector,tag_cvector
        character(*), intent(in) :: tag_a,tag_b,tag_c,tag_alpha,tag_beta,tag_gamma
        logical, intent(in) :: printable

        character(FMAXUNITLEN) :: readunit
	real(DP) :: s, al, be, ga, sin_be, sin_ga, wk, alpha, beta,gamma
	real(DP), dimension(3) :: wvec
	integer getRealVectorValue, realConvByUnit
	integer getRealValue
	integer i

	f_readUnitCell = -1
	if( getRealVectorValue( trim(TAG_AVECTOR)//char(0), wvec, readunit ) == 0 ) then
!!$           if(printable) then
!!$              write(6,'("  wvec(i), avec(i) = ",2f8.4)') wvec(i),avec(i)
!!$              write(6,'(" !! readunit = ",a16, " unit = ",a16)') readunit,unit
!!$           end if
		do i = 1, 3
		    if( realConvByUnit( wvec(i), avec(i), readunit, unit ) < 0 ) then
       		        return;
		    end if
!!$                    if(printable) write(6,'("  wvec(i), avec(i) = ",2f8.4)') wvec(i),avec(i)
		end do
		if( getRealVectorValue( trim(TAG_BVECTOR)//char(0), wvec, readunit ) /= 0 ) then
   		    return;
		end if
		do i = 1, 3
		    if( realConvByUnit( wvec(i), bvec(i), readunit, unit ) < 0 ) then
       		        return;
		    end if 
		end do
		if( getRealVectorValue( trim(TAG_CVECTOR)//char(0), wvec, readunit ) /= 0 ) then
   		    return;
		end if
		do i = 1, 3
		    if( realConvByUnit( wvec(i), cvec(i), readunit, unit ) < 0 ) then
       		        return;
		    end if 
		end do

		a = sqrt( sum(avec(:)*avec(:)) )
		b = sqrt( sum(bvec(:)*bvec(:)) )
		c = sqrt( sum(cvec(:)*cvec(:)) )
!!$		alpha = acos( dot_product(bvec,cvec)/(b*c) )*180.d0/PAI
!!$		beta = acos( dot_product(cvec,avec)/(c*a) )*180.d0/PAI
!!$		gamma = acos( dot_product(avec,bvec)/(a*b) )*180.d0/PAI
		ca =  dot_product(bvec,cvec)/(b*c)
		cb =  dot_product(cvec,avec)/(c*a)
		cc =  dot_product(avec,bvec)/(a*b)
                ucinptype = 1
	else if( getRealValue( trim(TAG_A)//char(0), wk, readunit ) == 0 ) then
!!$                if(printable) write(6,'(" !!! readunit = ",a15)') trim(readunit)
	 	if( realConvByUnit( wk, a, readunit, unit ) < 0 ) then
   		    return;
		end if

		if( getRealValue( trim(TAG_B)//char(0), wk, readunit ) /= 0 ) then
		    return;
		end if
	 	if( realConvByUnit( wk, b, readunit, unit ) < 0 ) then
   		    return;
		end if

		if( getRealValue( trim(TAG_C)//char(0), wk, readunit ) /= 0 ) then
		    return;
		end if
	 	if( realConvByUnit( wk, c, readunit, unit ) < 0 ) then
   		    return;
		end if

		if( getRealValue( trim(TAG_ALPHA)//char(0), wk, readunit ) /= 0 ) then
		    return;
		end if
		alpha = wk
	 	! if( realConvByUnit( wk, alpha, readunit, unit ) < 0 ) then
   		!     return;
		! end if

		if( getRealValue( trim(TAG_BETA)//char(0), wk, readunit ) /= 0 ) then
		    return;
		end if
		beta = wk
	 	! if( realConvByUnit( wk, beta, readunit, unit ) < 0 ) then
   		!     return;
		! end if

		if( getRealValue( trim(TAG_GAMMA)//char(0), wk, readunit ) /= 0 ) then
		    return;
		end if
		gamma = wk
	 	! if( realConvByUnit( wk, gamma, readunit, unit ) < 0 ) then
   		!    return;
		! end if

		al = alpha*PAI/180.d0
		be = beta*PAI/180.d0
		ga = gamma*PAI/180.d0
!!$		sin_al = sin(al)
		ca = cos(al)
		sin_be = sin(be)
		cb = cos(be)
		sin_ga = sin(ga)
		cc = cos(ga)
		s = (cc*cb-ca)/sin_ga/cb
		avec(1) = a
		avec(2) = 0.d0
		avec(3) = 0.d0
		bvec(1) = b*cc
		bvec(2) = b*sin_ga
		bvec(3) = 0.d0
		cvec(1) = c*cb
		cvec(2) = c*abs(s)*cb
		cvec(3) = c*sqrt(sin_be*sin_be-s*s*cb*cb)
                ucinptype = 0
	else
                ucinptype = -1
		return
	end if
	
	f_readUnitCell = 0
	return
end function f_readUnitCell


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! intelligent reading of Symmetry Generator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! integer function f_readSymmetryGenerator(sys,op,tvec)  ! old version 
integer function f_readSymmetryGenerator(sys,op,txu,txd,tyu,tyd,tzu,tzd,tag_rotation,tag_tx,tag_ty,tag_tz,printable)
  use m_Const_Parameters, only : DP, FMAXVALLEN, num_d6h, num_oh, d6h_symbol, oh_symbol, LOWER, NOCONV

	implicit none

	integer, intent(in) :: sys
	integer, intent(out) :: op
	!!! real(8), dimension(3), intent(out) :: tvec  ! old version
        integer, intent(out) :: txu, txd, tyu, tyd, tzu, tzd
        character(*),intent(in) :: tag_rotation,tag_tx,tag_ty,tag_tz
        logical, intent(in) :: printable
	integer getIntValue, getStringValue, frac2intint            !!!frac2real
	character(FMAXVALLEN) opstr, tstr
	integer i, iret

	if( getIntValue( trim(TAG_ROTATION)//char(0), iret ) == 0 ) then
	    if( sys <= 0 .and. (iret <=0 .or. iret > num_d6h) ) then
                print '("In f_readSymmetryGenerator: op =",I4," (1 < op <  ",i2,")" )', iret,num_d6h
	        f_readSymmetryGenerator = -1
		return
	    else if( sys > 0 .and. (iret <=0 .or. iret > num_oh) ) then
                print '("In f_readSymmetryGenerator: op =",I4," (1 < op <  ",i2,")" )', iret,num_oh
	        f_readSymmetryGenerator = -1
		return
	    end if
	    op = iret
	else
	f_readSymmetryGenerator = getStringValue( trim(TAG_ROTATION)//char(0), opstr, LOWER )
	if( f_readSymmetryGenerator /= 0 ) then
		return
	end if
	if( sys <= 0 ) then	! for rhombo, hexagonal  !!! should be checked
		do i = 1, num_d6h
			if( opstr == d6h_symbol(i) ) then
				exit
			end if
		end do
		op = i
                if(op > num_d6h) goto 1001
	else
		do i = 1, num_oh
			if( opstr == oh_symbol(i) ) then
				exit
			end if
		end do
		op = i
                if(op > num_oh) goto 1001
	end if
	endif
	f_readSymmetryGenerator = getStringValue( trim(TAG_TX)//char(0), tstr, NOCONV  )
	if( f_readSymmetryGenerator /= 0 ) then
		return
	end if
	!!! f_readSymmetryGenerator = frac2real( trim(tstr)//char(0), tvec(1) )    ! old version
	f_readSymmetryGenerator = frac2intint( trim(tstr)//char(0), txu, txd )
	if( f_readSymmetryGenerator /= 0 ) then
		return
	end if
	f_readSymmetryGenerator = getStringValue( trim(TAG_TY)//char(0), tstr, NOCONV )
	if( f_readSymmetryGenerator /= 0 ) then
		return
	end if
	!!! f_readSymmetryGenerator = frac2real( trim(tstr)//char(0), tvec(2) )     ! old version
	f_readSymmetryGenerator = frac2intint( trim(tstr)//char(0), tyu, tyd )
	if( f_readSymmetryGenerator /= 0 ) then
		return
	end if
	f_readSymmetryGenerator = getStringValue( trim(TAG_TZ)//char(0), tstr, NOCONV )
	if( f_readSymmetryGenerator /= 0 ) then
		return
	end if
	!!! f_readSymmetryGenerator = frac2real( trim(tstr)//char(0), tvec(3) )     ! old version
	f_readSymmetryGenerator = frac2intint( trim(tstr)//char(0), tzu, tzd )
	if( f_readSymmetryGenerator /= 0 ) then
		return
	end if
	return
1001    continue
        if(printable) then
           write(6,'(" !** illegal operator = [ ",a4,"]")') opstr
           if(sys <= 0) then
              write(6,'(" !** op = ",i6," > num_d6h = ",i6)') op,num_d6h
           else
              write(6,'(" !** op = ",i6," > num_oh = ",i6)') op,num_oh
           end if
        end if
        stop ' !* illegal operation generator <f_readSymmetryGenerator>'

end function f_readSymmetryGenerator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! intelligent reading of K-Point Coordinates
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer function f_readKPoints(tag_kx,tag_ky,tag_kz,tag_k_denom,tag_k_weight &
     &                        ,kvec,weight)

  use m_Const_Parameters, only : DP

	implicit none
        character(*), intent(in) :: tag_kx, tag_ky, tag_kz, tag_k_denom, tag_k_weight
	real(DP), dimension(3), intent(out) :: kvec
	integer, intent(out) :: weight
	integer getIntValue, i, wk
	real(DP) a

	f_readKPoints = getIntValue(trim(TAG_KX)//char(0),wk)
        kvec(1) = dble(wk)
	if( f_readKPoints /= 0 ) then
		return
	end if 
	f_readKPoints = getIntValue(trim(TAG_KY)//char(0),wk)
        kvec(2) = dble(wk)
	if( f_readKPoints /= 0 ) then
		return
	end if 
	f_readKPoints = getIntValue(trim(TAG_KZ)//char(0),wk)
        kvec(3) = dble(wk)
	if( f_readKPoints /= 0 ) then
		return
	end if 
	f_readKPoints = getIntValue(trim(TAG_K_DENOM)//char(0),wk)
        a = dble(wk)
	if( f_readKPoints /= 0 ) then
		return
	end if

	kvec(:) = kvec(:)/a 

	f_readKPoints = getIntValue(trim(TAG_K_WEIGHT)//char(0),weight)
	if( f_readKPoints /= 0 ) then
		return
	end if 
	return

end function f_readKPoints

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  convert to lower case string
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine f_strlower( str )
    implicit none
    character(*), intent(inout) :: str
    integer i, slen

    slen = len_trim(str)
    do i = 1, slen
        if( str(i:i) >= 'A' .and. str(i:i) <= 'Z' ) then
	    str(i:i) = achar( iachar(str(i:i)) - iachar('A') + iachar('a') )
        end if
    end do

end subroutine f_strlower


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  check s1 == s2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
logical function fstreq( s1, s2 )
    implicit none
    character(*), intent(in) :: s1
    character(*), intent(in) :: s2
    integer len1, len2, i

    len1 = len_trim(s1)
    len2 = len_trim(s2)
   !! print '("len1=",i8,"s1=[",a,"]")', len1, s1
   !! print '("len2=",i8,"s2=[",a,"]")', len2, s2
    if( len1 /= len2 ) then
        fstreq = .false.
        return
    end if

    do i = 1, len1
        if( s1(i:i) /= s2(i:i) ) then
            fstreq = .false.
            return
        end if       
    end do	
    fstreq = .true.
    return
end function fstreq

!!$end module input_interface
