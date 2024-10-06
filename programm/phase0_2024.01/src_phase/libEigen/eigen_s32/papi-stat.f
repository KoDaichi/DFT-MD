      subroutine papi_initialize()
      implicit NONE
      include 'f90papi.h'

      integer   :: i, flags, res

      integer, parameter :: counter_max = 10
      integer*4 :: counter_list(counter_max)
      integer   :: event_set
      integer*8 :: i_results(counter_max)
      real*8    :: fp_results(counter_max)
      integer   :: PAPI_test
      integer   :: PAPI_dummy

      character*(PAPI_MAX_STR_LEN) :: evt_name(10)
      logical      :: PAPI_ready

      common /papi_common/ PAPI_ready, event_set, PAPI_test, PAPI_dummy,
     $                     counter_list, evt_name, i_results

      res = PAPI_VER_CURRENT
      call PAPIf_library_init(res)
      if(res/=PAPI_VER_CURRENT)then
         print*,"PAPI initialize failed",res
      endif
      PAPI_ready=(res==PAPI_VER_CURRENT)
      if(PAPI_ready)then
         counter_list(1)= PAPI_TLB_DM
         evt_name(1)    ="PAPI_TLB_DM"
         counter_list(2)= PAPI_L1_DCM
         evt_name(2)    ="PAPI_L1_DCM"
         counter_list(3)= PAPI_L2_DCM
         evt_name(3)    ="PAPI_L2_DCM"
         counter_list(4)= PAPI_TLB_TL
         evt_name(4)    ="PAPI_TLB_TL"
         counter_list(5)= PAPI_TLB_IM
         evt_name(5)    ="PAPI_TLB_IM"
         counter_list(6)= PAPI_RES_STL
         evt_name(6)    ="PAPI_RES_STL"

         PAPI_test = 2

         event_set = PAPI_NULL
         call PAPIf_create_eventset(event_set, res)
         if(res/=PAPI_OK)then
            print*,"PAPI event cannot be created",res
            stop
         endif
         do i=1,PAPI_test
            call PAPIf_query_event(counter_list(i), res)
            if(res/=PAPI_OK)then
               print*,"PAPI event is not supported"
               stop
            endif
            call PAPIf_add_event(event_set, counter_list(i), res)
            if(res/=PAPI_OK)then
               print*,"PAPI event cannot be added"
               stop
            endif
         enddo
         print*,"PAPI all events are initialized"
         PAPI_ready=(res==PAPI_OK)
      endif

      end subroutine

      subroutine papi_timer_start()

      include 'f90papi.h'

      integer   :: flags, res

      integer, parameter :: counter_max = 10
      integer*4 :: counter_list(counter_max)
      integer   :: event_set
      integer*8 :: i_results(counter_max)
      real*8    :: fp_results(counter_max)
      integer   :: PAPI_test
      integer   :: PAPI_dummy

      character*(PAPI_MAX_STR_LEN) :: evt_name(10)
      logical      :: PAPI_ready

      common /papi_common/ PAPI_ready, event_set, PAPI_test, PAPI_dummy,
     $                     counter_list, evt_name, i_results

         if(PAPI_ready)then
            call PAPIf_start(event_set, res)
            if(res/=PAPI_OK)then
               print*,"ERR PAPIf_start"
            endif
         endif

      end subroutine

      subroutine papi_timer_stop()

      include 'f90papi.h'

      integer   :: flags, res

      integer, parameter :: counter_max = 10
      integer*4 :: counter_list(counter_max)
      integer   :: event_set
      integer*8 :: i_results(counter_max)
      real*8    :: fp_results(counter_max)
      integer   :: PAPI_test
      integer   :: PAPI_dummy

      character*(PAPI_MAX_STR_LEN) :: evt_name(10)
      logical      :: PAPI_ready

      common /papi_common/ PAPI_ready, event_set, PAPI_test, PAPI_dummy,
     $                     counter_list, evt_name, i_results

         if(PAPI_ready)then
            call PAPIf_stop(event_set, i_results, res)
            if(res/=PAPI_OK)then
               print*,"ERR PAPIf_stop"
            endif
         endif

      end subroutine

      subroutine papi_timer_report()

      include 'f90papi.h'

      integer   :: flags, res

      integer, parameter :: counter_max = 10
      integer*4 :: counter_list(counter_max)
      integer   :: event_set
      integer*8 :: i_results(counter_max)
      real*8    :: fp_results(counter_max)
      integer   :: PAPI_test
      integer   :: PAPI_dummy

      character*(PAPI_MAX_STR_LEN) :: evt_name(10)
      logical      :: PAPI_ready

      common /papi_common/ PAPI_ready, event_set, PAPI_test, PAPI_dummy,
     $                     counter_list, evt_name, i_results

         if(PAPI_ready)then
            print*,">>>>> PAPI result start <<<<<"
            do icc=1,PAPI_test
               print*,evt_name(icc)(1:16),i_results(icc)
            enddo
            print*,">>>>> PAPI result end   <<<<<"
         endif

      end subroutine

