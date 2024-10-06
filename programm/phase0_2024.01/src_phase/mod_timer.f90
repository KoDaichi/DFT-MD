module mod_timer
implicit none
integer, private, parameter :: num_max_routines = 100 ! Number of max routines.
integer, private, parameter :: num_max_namelen = 100 ! Max length of routine name.
integer, private :: num_of_routines = 0 ! Counter for number of routines.

character(len=num_max_namelen), private, dimension(num_max_routines) :: t_name = "" ! Array for routine name.
real(kind=8), private, dimension(num_max_routines) :: t_value = 0.0d0, ts ! Timer values.
integer, private, dimension(num_max_routines) :: call_count = 0 ! Counter of routine call.
logical, private, dimension(num_max_routines) :: timer_start = .false. ! Flag for timer is started or not.

private :: wallclock

contains

   subroutine start_timer(subroutine_name)
      character(len=*), intent(in) :: subroutine_name
      integer :: nlen
      integer :: i

      nlen=len(subroutine_name)

      do i = 1, num_of_routines
         if(subroutine_name(1:nlen) == trim(t_name(i))) then
            if(timer_start(i)) then
               write(0,*) 'Timer for ', trim(subroutine_name), ' is already started!!!'
               return
            end if
            timer_start(i) = .true.
            call_count(i) = call_count(i) + 1
            call wallclock(ts(i))
            return
         end if
      end do

      num_of_routines = num_of_routines + 1
      t_name(num_of_routines) = subroutine_name(1:nlen)
      timer_start(num_of_routines) = .true.
      call_count(num_of_routines) = call_count(num_of_routines) + 1
      call wallclock(ts(num_of_routines))

      return
   end subroutine start_timer

   subroutine stop_timer(subroutine_name)
      character(len=*), intent(in) :: subroutine_name
      integer :: nlen
      integer :: i
      real(kind=8) :: te_tmp

      call wallclock(te_tmp)

      nlen=len(subroutine_name)

      do i = 1, num_of_routines
         if(subroutine_name(1:nlen) == trim(t_name(i))) then
            if(.not. timer_start(i)) then
               write(0,*) 'Timer for ', trim(subroutine_name), ' is NOT started!!!'
               return
            end if
            timer_start(i) = .false.
            t_value(i) = t_value(i) + (te_tmp - ts(i))
            return
         end if
      end do

      write(0,*) 'Timer for ', trim(subroutine_name), ' is NOT started!!!'

      return
   end subroutine stop_timer

   subroutine print_timer(mype)
      integer, intent(in) :: mype
      integer :: i
      character(len=31) :: p_name

      write(40000+mype,'(a)') ''
      write(40000+mype,'(a)') '[Timer Output]'
      write(40000+mype,'(a)') '+-------------------------------+----------+------------+'
      write(40000+mype,'(a)') '|Timer region                   |Called    |Elapsed     |'
      write(40000+mype,'(a)') '|                               |          |Time[s]     |'
      write(40000+mype,'(a)') '+-------------------------------+----------+------------+'
      do i = 1, num_of_routines
         if(timer_start(i)) then
            write(0,*) 'Timer for ', trim(t_name(i)), ' is NOT stopped!!!'
            return
         end if
         p_name = ''
         p_name = trim(t_name(i))
         write(40000+mype,'(a,a,i10,a,f12.3,a)') '|', p_name // '|', call_count(i), '|', t_value(i), '|'
      end do
      write(40000+mype,'(a)') '+-------------------------------+----------+------------+'

      return
   end subroutine print_timer

   subroutine wallclock(t)
      real(kind=8), intent(out) :: t
      integer(kind=8) :: c, c_rate

      call system_clock(c, c_rate)

      t = dble(c)/dble(c_rate)

      return
   end subroutine wallclock

end module mod_timer
