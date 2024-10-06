       integer                :: TRD_inod,TRD_nnod,TRD_COMM_WORLD
       common /USEMPI/           TRD_inod,TRD_nnod,TRD_COMM_WORLD

       integer                :: x_inod, x_nnod, x_COMM_WORLD,
     &                           y_inod, y_nnod, y_COMM_WORLD,
     &                           p0_(1024),q0_(1024), n_common,
     &                           diag_0, diag_1
       common /CYCL2D/           x_inod, x_nnod, x_COMM_WORLD,
     &                           y_inod, y_nnod, y_COMM_WORLD,
     &                           p0_      ,q0_      , n_common,
     &                           diag_0, diag_1

       external               :: loop_c
       integer, external      :: loop_c_start, loop_c_end
       integer, external      :: loop_c_node
       integer, external      :: loop_c_g2l_depth
       integer, external      :: loop_c_l2g_depth

       integer, parameter :: nsx = 480
       integer, parameter :: nsm = 256 ! 64 +32

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  for DC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer          :: mpi_comm_eigen, ICTXTtoDC, MYROWtoDC,
     &                    MYCOLtoDC
      common /COMMtoDC/   mpi_comm_eigen, ICTXTtoDC, MYROWtoDC,
     &                    MYCOLtoDC
