      module  eigen_devel
!$    use OMP_LIB
      use MPI, only : MPI_COMM_WORLD, MPI_THREAD_SINGLE
      implicit NONE


         public ! default is public

!$       integer                :: MPI_THREAD_MODE = MPI_THREAD_SINGLE
         character*(8), private :: Process_Grid_Major
         logical, private       :: Eigen_initialized_flag = .FALSE.

         integer                :: TRD_inod       = 0
         integer                :: TRD_nnod       = 0
         integer                :: TRD_COMM_WORLD = MPI_COMM_WORLD

         integer                :: x_inod, x_nnod, x_COMM_WORLD,
     &                             y_inod, y_nnod, y_COMM_WORLD,
     &                             z_inod, z_nnod, z_COMM_WORLD,
     &                             w_inod, w_nnod, w_COMM_WORLD,
     &                             n_common, diag_0, diag_1

         integer, pointer       :: p0_(:), q0_(:)

         integer, parameter     :: nsx = 480
         integer, parameter     :: nsm = 256 ! 64 +32
         integer, parameter     :: ns0 = nsm * nsm + 6
         integer, parameter     :: MBAND = 2

         real(8), pointer       :: rt_timer(:,:)
         integer                :: rt_ptr(0:16)

         real(8)                :: time_bcast, time_reduce,
     &                             time_redist, time_gather


      contains


         real(8) function eigen_get_wtime( ) result(ret)
!$       use OMP_LIB
         use MPI
         implicit NONE

#ifdef _OPENMP
!$          ret = omp_get_wtime( )
#else
            ret = MPI_Wtime( )
#endif

         return
         end function  eigen_get_wtime

         subroutine eigen_set_initialized( )

            Eigen_initialized_flag = .TRUE.

         return
         end subroutine eigen_set_initialized

         subroutine eigen_unset_initialized( )

            Eigen_initialized_flag = .FALSE.

         return
         end subroutine eigen_unset_initialized

         subroutine eigen_get_initialized( flag )
         logical :: flag

            flag = Eigen_initialized_flag

         return
         end subroutine eigen_get_initialized

         subroutine eigen_set_grid_major( Major )
         character*(*) Major

            Process_Grid_Major(1:1) = Major(1:1)

         return
         end subroutine eigen_set_grid_major

         subroutine eigen_get_grid_major( Major )
         character*(*) Major

            Major(1:1) = Process_Grid_Major(1:1)

         return
         end subroutine eigen_get_grid_major

         subroutine  eigen_timer_reset( )

            time_bcast  = 0D0
            time_reduce = 0D0
            time_redist = 0D0
            time_gather = 0D0

         return
         end subroutine eigen_timer_reset

         subroutine eigen_timer_print( )

#if TIMER_PRINT
            if ( TRD_inod == 1 ) then
               print*, "COMM_STAT"
               print*, "   BCAST  :: ", time_bcast
               print*, "   REDUCE :: ", time_reduce
               print*, "   REDIST :: ", time_redist
               print*, "   GATHER :: ", time_gather
            end if
#endif

         return
         end subroutine  eigen_timer_print

!$       subroutine sync_other_than_master( TRBK_lock )
!$       use OMP_LIB
!$       integer(kind=omp_lock_kind) :: TRBK_lock(1:2)
!$       integer(kind=omp_lock_kind) :: TMP_lock(1:2)
!$
!$          if ( omp_get_num_threads() == 1 ) return
!$          if ( omp_get_thread_num() == 0 ) return
!$
!$          TMP_lock(1) = 2**omp_get_num_threads() - 2
!$
!$          call omp_set_lock( TRBK_lock(1) )
!$          TRBK_lock(2) = ibset( TRBK_lock(2), omp_get_thread_num() )
!$          call omp_unset_lock( TRBK_lock(1) )
!$
!$          do
!$             call omp_set_lock( TRBK_lock(1) )
!$             TMP_lock(2) = TRBK_lock(2)
!$             call omp_unset_lock( TRBK_lock(1) )
!$             if ( TMP_lock(1) == TMP_lock(2) ) exit
!$          end do
!$
!$       end subroutine  sync_other_than_master


      end module  eigen_devel

