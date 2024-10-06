module m_PAW_Tecplot
    use m_Const_Parameters,     only : DP
    implicit none
    private
    
    public:: output_box_tecplot
    public:: output_box_tecplot2
    public:: unusedUnit
    
contains
    
    subroutine output_box_tecplot(v_data,vec,title,axis,unit,mode)
        real(kind=DP),intent(in),dimension(0:,0:,0:) :: v_data
        real(kind=DP),dimension(3,3),intent(in):: vec
        character(len=*),intent(in):: title,axis
        integer,intent(in):: unit,mode

        integer:: i,j,k,ie,je,ke
        real(DP):: x,y,z,xx,yy,zz,vec2(3,3)

        ie=size(v_data,1)-1
        je=size(v_data,2)-1
        ke=size(v_data,3)-1
        
        vec2=vec
        
        if(mode.eq.0) &
            write(unit,'(3a)') 'VARIABLES = "X", "Y", "Z", "',trim(axis),'"'
        write(unit,'(3a,i3,a,i3,a,i3,a)') &
           'ZONE T=',trim(title),' I=',ie+1,', J=',je+1,', K=',ke+1,', F=POINT'
        do k=0,ke
	        do j=0,je
	            do i=0,ie
                    if(ie.ne.0) xx=dble(i)/dble(ie)
                    if(je.ne.0) yy=dble(j)/dble(je)
                    if(ke.ne.0) zz=dble(k)/dble(ke)
	                x=vec2(1,1)*xx+vec2(1,2)*yy+vec2(1,3)*zz
	                y=vec2(2,1)*xx+vec2(2,2)*yy+vec2(2,3)*zz
                    z=vec2(3,1)*xx+vec2(3,2)*yy+vec2(3,3)*zz
	                write(unit,'(4e19.6)') x,y,z,v_data(i,j,k)
                end do
            end do
        end do
        return
    end subroutine output_box_tecplot
    
    subroutine output_box_tecplot2(r_data,v_data,title,axis,unit,mode)
        real(kind=DP),intent(in),dimension(1:,1:,1:,1:)   :: r_data
        real(kind=DP),intent(in),dimension(1:,1:,1:)      :: v_data
        character(len=*),intent(in):: title,axis
        integer,intent(in):: unit,mode

        integer:: i,j,k,ie,je,ke

        ie=size(v_data,1)
        je=size(v_data,2)
        ke=size(v_data,3)
        
        if(mode.eq.0) &
            write(unit,'(3a)') 'VARIABLES = "X", "Y", "Z", "',trim(axis),'"'
        write(unit,'(3a,i3,a,i3,a,i3,a)') &
           'ZONE T=',trim(title),' I=',ie,', J=',je,', K=',ke,', F=POINT'
        do k=1,ke
	        do j=1,je
	            do i=1,ie
	                write(unit,'(4e19.6)') r_data(i,j,k,1), &
	                                        r_data(i,j,k,2), &
	                                        r_data(i,j,k,3), &
	                                        v_data(i,j,k)
                end do
            end do
        end do
        return
    end subroutine output_box_tecplot2
    
    function unusedUnit()    
        integer:: unusedUnit
        integer:: unit
        logical:: opened
        do unit=99,10,-1
        inquire(unit,opened=opened)
        if(.not.opened) then
            unusedUnit = unit
            return
        endif
        enddo
        write(6,*) 'All devices are used!'
        stop
    end function unusedUnit


end module m_PAW_Tecplot