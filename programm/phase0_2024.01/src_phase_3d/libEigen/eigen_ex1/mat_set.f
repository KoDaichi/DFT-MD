       subroutine mat_set( n, a, nm, mtype )
       use MPI
       use eigen_libs
       implicit NONE

       integer, intent(in)   :: n, nm, mtype
       real(8)               :: a(1:nm,*)

! Parameters BLACS array descritor(the position of entry tags), etc
       INTEGER, PARAMETER           :: BLOCK_CYCLIC_2D = 1
       INTEGER, PARAMETER           :: DLEN_  = 9
       INTEGER, PARAMETER           :: DTYPE_ = 1
       INTEGER, PARAMETER           :: CTXT_  = 2
       INTEGER, PARAMETER           :: M_     = 3
       INTEGER, PARAMETER           :: N_     = 4
       INTEGER, PARAMETER           :: MB_    = 5
       INTEGER, PARAMETER           :: NB_    = 6
       INTEGER, PARAMETER           :: RSRC_  = 7
       INTEGER, PARAMETER           :: CSRC_  = 8
       INTEGER, PARAMETER           :: LLD_   = 9

       integer                      :: DESCA( DLEN_ )
       integer                      :: NPROW, NPCOL

       real(8), pointer             :: as(:,:)
       real(8)                      :: t

       integer                      :: nnod, x_nnod, y_nnod
       integer                      :: inod, x_inod, y_inod
       integer                      :: i, i_1, i_2, i_3
       integer                      :: j, j_1, j_2, j_3
       integer                      :: nx, info, ICTXT
       integer, pointer             :: iseed(:)


          call eigen_get_procs( nnod, x_nnod, y_nnod )
          call eigen_get_id   ( inod, x_inod, y_inod )

          j_2 = eigen_loop_start( 1, x_nnod, x_inod )
          j_3 = eigen_loop_end  ( n, x_nnod, x_inod )
          i_2 = eigen_loop_start( 1, y_nnod, y_inod )
          i_3 = eigen_loop_end  ( n, y_nnod, y_inod )

          if ( mtype == 0 ) then
             do i_1 = i_2, i_3
                i = eigen_translate_l2g( i_1, y_nnod, y_inod )
                do j_1 = j_2, j_3
                   j = eigen_translate_l2g( j_1, x_nnod, x_inod )
                   a(j_1, i_1) = (n+1-Max(n+1-i,n+1-j))*1.0D+00
                end do
             end do
          end if
          if ( mtype == 1 ) then
             do i_1 = i_2, i_3
                i = eigen_translate_l2g( i_1, y_nnod, y_inod )
                do j_1 = j_2, j_3
                   j = eigen_translate_l2g( j_1, x_nnod, x_inod )
                   if ( i == j ) then
                      a(j_1, i_1) = -7.2D0
                   else
                      a(j_1, i_1) = -3.0D0/(i-j)**2
                   end if
                end do
             end do
          end if
          if ( mtype == 2 ) then

             NPROW = x_nnod
             NPCOL = y_nnod

             ICTXT = eigen_get_blacs_context( )

             CALL DESCINIT( DESCA, n, n, 1, 1, 0, 0, ICTXT, nm, INFO )

             nx = (n-1)/NPCOL+1
             allocate( as(1:nm,nx) )

             CALL RANDOM_SEED( size = i )
             allocate( iseed(i) )
             iseed(1) = inod
             CALL RANDOM_SEED( put = iseed )
             deallocate( iseed )

             do i_1 = i_2, i_3
                do j_1 = j_2, j_3
                   CALL RANDOM_NUMBER( t )
                   as(j_1, i_1) = t
                   a (j_1, i_1) = t
                end do
             end do

             CALL PDTRAN( n, n,
     &                    1D0, as, 1, 1, DESCA, 1D0, a, 1, 1, DESCA )

             deallocate( as )

          end if


       return
       end

