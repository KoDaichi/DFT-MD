
      integer                :: nprocs,myrank
      common /usempi/           nprocs,myrank

      integer                :: my_col, size_of_col, mpi_comm_col,
     &                          my_row, size_of_row, mpi_comm_row,
     &                          p0_(1024),q0_(1024), n_common, diag_0
      common /cycl2d/           my_col, size_of_col, mpi_comm_col,
     &                          my_row, size_of_row, mpi_comm_row,
     &                          p0_      ,q0_      , n_common, diag_0
      integer                :: es_cluster,
     &                          es_icnod, es_ncnod, es_ccomm_world,
     &                          es_ignod, es_ngnod, es_gcomm_world
      common /escy2d/           es_cluster,
     &                          es_icnod, es_ncnod, es_ccomm_world,
     &                          es_ignod, es_ngnod, es_gcomm_world

!     integer, external      :: get_loop_start, get_loop_end
!     integer, external      :: get_owner_node
!     integer, external      :: translate_g2l
!     integer, external      :: translate_l2g

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  for DC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer          :: mpi_comm_eigen, ICTXTtoDC, MYROWtoDC,
     &                    MYCOLtoDC
      common /COMMtoDC/   mpi_comm_eigen, ICTXTtoDC, MYROWtoDC,
     &                    MYCOLtoDC
