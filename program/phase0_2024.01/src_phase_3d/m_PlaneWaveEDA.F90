module m_PlaneWaveEDA

#ifdef __EDA__
!  use m_Charge_Density,                   only  : chgq_l_EDA, chgq_l, hsr
  use m_Charge_Density,                   only  : chgq_l, hsr
  use m_Const_Parameters,                 only  : DIRECT, DP, EXC_ONLY, INVERSE, OFF, ON, Valence_plus_PC_Charge &
       &                                        , ELECTRON, GAMMA
  use m_Control_Parameters,               only  : af, kimg, nspin, neg
  use m_Crystal_Structure,                only  : altv, rltv, univol
  use m_Electronic_Structure,             only  : m_ES_WF_in_Rspace_3D, m_ES_WF_kin_in_Rspace &
       &                                        , occup_l, totch, zaj_l
  use m_FFT,                              only  : m_FFT_alloc_WF_work, m_FFT_dealloc_WF_work &
       &                                        , m_FFT_alloc_CD_box, m_FFT_dealloc_CD_box &
       &                                        , m_FFT_CD_direct0, m_FFT_CD_inverse0 &
       &                                        , m_FFT_WF &
       &                                        , fft_box_size_CD_nonpara, fft_box_size_WF, nfft, nfftp, nfftp_nonpara
  use m_Files,                            only  : nfout
  use m_Ionic_System,                     only  : cps, iatomn_EDA, ityp, iwei, natm, natm2, ntyp, pos, speciesname, zfm3_l &
       &                                        , eewald_per_atom, icount_EDA, zfm3_l_EDA
  use m_Kpoints,                          only  : kv3, vkxyz
  use m_Parallelization,                  only  : ista_e, iend_e, istep_e &
       &                                        , ista_kngp, iend_kngp &
       &                                        , map_k, map_e, map_z,  myrank_k, myrank_e, np_e, nrank_e, nrank_k, myrank_g
  use m_PlaneWaveBasisSet,                only  : gr_l, iba, igf, igfp_nonpara, kg, nbase, nbase_gamma, nbmx, ngabc, kg1
  use m_PseudoPotential,                  only  : dion, etot1, ilmt, ltp, mtp, psc_l &
       &                                        , epc_per_atom, PP_per_atom
  use m_XC_Potential_2D,                  only  : exc_on_a_grid, m_XC_cal_potential
  use m_Parallelization,                  only  : MPI_CommGroup, mype, npes &
       &                                        , nel_fft_x, nel_fft_y, nel_fft_z, mpi_chg_world, mpi_ke_world &
       &                                        , mpi_ge_world, mpi_kg_world, ista_g1k, iend_g1k, map_ek, myrank_ke
  use m_Timing,                           only  : tstatc0_begin, tstatc0_end, tstatc_wd, tstatc_init, tstatc_iter
  use m_Kpoints,                          only  : k_symmetry
  use mpi

  implicit none
!  include 'mpif.h'

! -----  for partitioning function  -----

  real(kind=DP),private,allocatable,dimension(:)       :: bragg_slater_radii
  real(kind=DP),private,dimension(3,3)                 :: rltv_inverse       ! <- inverse matrix of rltv
  real(kind=DP), allocatable, dimension(:,:)           :: rltv_t
  real(kind=DP),private,allocatable,dimension(:,:)     :: aij
  real(kind=DP),private,allocatable,dimension(:,:)     :: r_catom            ! <- internal vectors of atomic centers
  real(kind=DP),private,allocatable,dimension(:,:)     :: cps_new            ! <- real-space vectors of atomic centers

  real(kind=DP),private,allocatable,dimension(:)       :: grid_x_CD, grid_y_CD, grid_z_CD
  integer,private                                      :: ifxr_CD, ifyr_CD, ifzr_CD
  real(kind=DP),private                                :: grid_weight_CD
  real(kind=DP),private,allocatable,dimension(:,:,:,:) :: partition_function
  real(kind=DP),private,allocatable,dimension(:,:)     :: partition_function_on_a_grid
  real(kind=DP),allocatable,dimension(:)               :: afft
  real(kind=DP),allocatable,dimension(:)               :: bfft1
  real(kind=DP),allocatable,dimension(:)               :: bfft2
  real(kind=DP),allocatable,dimension(:)               :: density_on_a_grid
  real(kind=DP),allocatable,dimension(:)               :: ekin_on_a_grid
  real(kind=DP),allocatable,dimension(:)               :: ehartr_on_a_grid
  real(kind=DP),allocatable,dimension(:)               :: ehartr_on_a_grid2
  real(kind=DP),allocatable,dimension(:)               :: elocal_on_a_grid
  real(kind=DP),allocatable,dimension(:)               :: elocal_on_a_grid1
  real(kind=DP),allocatable,dimension(:)               :: elocal_on_a_grid2
  real(kind=DP),allocatable,dimension(:,:)             :: PP_on_a_grid
  real(kind=DP),allocatable,dimension(:)               :: PP_on_a_grid_wk

! -----  EDA  -----

  real(kind=DP),private,allocatable,dimension(:)       :: density_per_atom
  real(kind=DP),private,allocatable,dimension(:)       :: ekin_per_atom
  real(kind=DP),private,allocatable,dimension(:)       :: ehartr_per_atom
  real(kind=DP),private,allocatable,dimension(:)       :: exc_per_atom
  real(kind=DP),private,allocatable,dimension(:)       :: elocal_per_atom
  real(kind=DP),private,allocatable,dimension(:)       :: elocal_per_atom_from_electron
  real(kind=DP),private,allocatable,dimension(:)       :: elocal_from_nucleus1
  real(kind=DP),private,allocatable,dimension(:)       :: elocal_from_nucleus2
  real(kind=DP),private,allocatable,dimension(:)       :: enonlocal_per_atom
  real(kind=DP),private,allocatable,dimension(:)       :: epcx_per_atom
  real(kind=DP),private,allocatable,dimension(:)       :: etotal_per_atom

! -----  for FFT  -----

  real(kind=DP),allocatable,dimension(:)               :: wf_EDA        ! <- wavefunction in reciprocal space
  real(kind=DP),allocatable,dimension(:)               :: wf_kin_EDA    ! <- |k+G|^2 * |k+G>
  real(kind=DP),allocatable,dimension(:)               :: wf_hartr_EDA  ! <- 4pi*rho(G)/|G|^2 * |k+G>

! -----  Parameters of G-vectors  -----

  real(kind=DP),private,allocatable,dimension(:,:,:)   :: ekin_l
  real(kind=DP),private,allocatable,dimension(:,:,:)   :: ehartr_l
  real(kind=DP),private,allocatable,dimension(:,:)     :: elocal_l

! -----  Others  -----
  real(kind=DP),dimension(6)                           :: ttl

contains

  subroutine m_PW_EDA_allocation
    allocate(afft(nfftp_nonpara))
    allocate(bfft1(nfftp_nonpara))
    allocate(bfft2(nfftp_nonpara))
    call m_FFT_alloc_CD_box()
  end subroutine m_PW_EDA_allocation

  subroutine m_PW_EDA_deallocation
    deallocate(afft)
    deallocate(bfft1)
    deallocate(bfft2)
!   call m_FFT_dealloc_CD_box()
  end subroutine m_PW_EDA_deallocation

  subroutine execute_PlaneWaveEDA
    use m_IterationNumbers, only : iteration, first_iteration_of_this_job

    call tstatc_init
    call write_headermark_of_PW_EDA

       WRITE(6,*) 'before_Partitioning_Function'
    call before_Partitioning_Function                                         ! <- set Becke's space-partitioning function

       WRITE(6,*) 'set_partition_function'
    call set_partition_function

       WRITE(6,*) 'Each_Energy_Term'
    call Each_Energy_Term                                               ! <- calculating each energy term per atom

       WRITE(6,*) 'write_energies_per_atom'
    call write_energies_per_atom

    call write_endmark_of_PW_EDA
    call tstatc_iter(iteration, first_iteration_of_this_job)
    call tstatc_wd(iteration)

  end subroutine execute_PlaneWaveEDA

!-------------------------------------------------------------------------------
!   >>>>> preparing partitioning function >>>>>
!-------------------------------------------------------------------------------

  subroutine before_Partitioning_Function

!-------------------------------------------------------------------------------
!   1. set Bragg-Slater radii for all atoms
!       - subroutine set_bragg_slater_radii
!
!   2. set inverse matrix of rltv for calculating atom-center coordinates in G-space
!       - subroutine inverse_matrix
!
!   3. set aij: determined only by Bragg-Slater radii
!       - subroutine set_aij
!
!   4. set grid points in G-space, which are the same those used in FFT and
!      calculate partitioning function on each grid.
!       - subroutine set_grid_points
!-------------------------------------------------------------------------------

    allocate(bragg_slater_radii(137))
    call set_bragg_slater_radii

    call inverse_matrix(rltv,rltv_inverse)

    allocate(aij(natm2,natm2))
    call set_aij

    ifxr_CD = fft_box_size_CD_nonpara(1,0)
    ifyr_CD = fft_box_size_CD_nonpara(2,0)
    ifzr_CD = fft_box_size_CD_nonpara(3,0)

    grid_weight_CD = 1.d0/(dble(product(fft_box_size_CD_nonpara(1:3,0))))

    allocate(r_catom(natm2,3), cps_new(natm2,3), &
       &     grid_x_CD(ifxr_CD+1), grid_y_CD(ifyr_CD+1), grid_z_CD(ifzr_CD+1), &
       &     partition_function(ifxr_CD+1,ifyr_CD+1,ifzr_CD+1,natm2), partition_function_on_a_grid(natm2,nfftp_nonpara))

    call set_grid_points

 contains
      subroutine set_bragg_slater_radii

        real(kind=DP),dimension(137)               :: BSR_work

!-----------------------------------------------------------------------
!     Bragg-Slater radii for determining the relative size of the
!     polyhedra in the polyatomic integration scheme
!-----------------------------------------------------------------------

        data BSR_work  /0.35d0, 0.31d0, 1.45d0, 1.05d0, 0.85d0,  &
                     &  0.70d0, 0.65d0, 0.60d0, 0.50d0, 0.38d0,  &
                     &  1.80d0, 1.50d0, 1.25d0, 1.10d0, 1.00d0,  &
                     &  1.00d0, 1.00d0, 0.71d0, 2.20d0, 1.80d0,  &
                     &  1.60d0, 1.40d0, 1.35d0, 1.40d0, 1.40d0,  &
                     &  1.40d0, 1.35d0, 1.35d0, 1.35d0, 1.35d0,  &
                     &  1.30d0, 1.25d0, 1.15d0, 1.15d0, 1.15d0,  &
                     &  0.88d0, 2.35d0, 2.00d0, 1.80d0, 1.55d0,  &
                     &  1.45d0, 1.45d0, 1.35d0, 1.30d0, 1.35d0,  &
                     &  1.40d0, 1.60d0, 1.55d0, 1.55d0, 1.45d0,  &
                     &  1.45d0, 1.40d0, 1.40d0, 1.08d0, 2.60d0,  &
                     &  2.15d0, 1.95d0, 1.85d0, 1.85d0, 1.85d0,  &
                     &  1.85d0, 1.85d0, 1.85d0, 1.80d0, 1.75d0,  &
                     &  1.75d0, 1.75d0, 1.75d0, 1.75d0, 1.75d0,  &
                     &  1.75d0, 1.55d0, 1.45d0, 1.35d0, 1.35d0,  &
                     &  1.30d0, 1.35d0, 1.35d0, 1.35d0, 1.50d0,  &
                     &  1.90d0, 1.80d0, 1.60d0, 1.90d0, 1.27d0,  &
                     &  1.20d0, 2.60d0, 2.15d0, 1.95d0, 1.80d0,  &
                     &  1.75d0, 1.75d0, 1.75d0, 1.75d0, 1.75d0,  &
                     &  1.75d0, 1.75d0, 1.75d0, 1.75d0, 1.75d0,  &
                     &  1.75d0, 1.75d0, 1.75d0, 1.75d0, 1.75d0,  &
                     &  1.75d0, 1.75d0, 1.75d0, 1.75d0, 1.75d0,  &
                     &  1.75d0, 1.75d0, 1.75d0, 1.75d0, 1.75d0,  &
                     &  1.75d0, 1.75d0, 1.75d0, 1.75d0, 1.75d0,  &
                     &  1.75d0, 1.75d0, 1.75d0, 1.75d0, 1.75d0,  &
                     &  1.75d0, 1.75d0, 1.75d0, 1.75d0, 1.75d0,  &
                     &  1.75d0, 1.75d0, 1.75d0, 1.75d0, 1.75d0,  &
                     &  1.75d0, 1.75d0/

        bragg_slater_radii = BSR_work

      end subroutine set_bragg_slater_radii

      subroutine set_aij()

        real(kind=DP)                              :: aij_work
        integer                                    :: iatm, jatm
        real(kind=DP)                              :: chi, chi2

        IatmLoop: do iatm = 1, natm2
          JatmLoop: do jatm = 1, natm2

            chi  = bragg_slater_radii(iatomn_EDA(iatm))/bragg_slater_radii(iatomn_EDA(jatm))
            chi2 = (chi-1.d0)/(chi+1.d0)
            aij_work = chi2/(chi2*chi2-1.d0)

            if(aij_work > 0.5d0) then
              aij_work =  0.5d0
            elseif(aij_work < -0.5d0) then
              aij_work = -0.5d0
            endif

            aij(iatm,jatm) = aij_work

          enddo JatmLoop
        enddo IatmLoop

      end subroutine set_aij

      subroutine set_grid_points

        call calculate_ttl

        call internal_vector                                            ! <- internal vector of atomic centers

        call packuc_EDA                                                 ! <- packing atoms in the unit cell

        call determine_grid_points

      end subroutine set_grid_points

      subroutine calculate_ttl                                          ! <- subroutine calttl in b_Ldos_f77.F

        integer                                      :: i, j, jj

        ttl = 0.d0

        do j = 1, 3
          jj = mod(j+1,3)
          if(jj == 0) jj = 3

          do i = 1, 3
            ttl(j) = ttl(j) + altv(i,j)*altv(i,j)
            ttl(j+3) = ttl(j+3) + 2.d0*altv(i,j)*altv(i,jj)
          enddo
        enddo

      end subroutine calculate_ttl

      subroutine internal_vector

        real(kind=DP),dimension(3)                   :: r_catom_work
        real(kind=DP)                                :: pi2
        integer                                      :: i, j, k

!-----------------------------------------------------------------------
!     r_catom is internal vector of atomic centers
!-----------------------------------------------------------------------

        pi2 = 8.d0*datan(1.d0)

        do i = 1, natm2
          do j = 1, 3
            r_catom_work = 0.d0

            do k = 1, 3
              r_catom_work(j) = r_catom_work(j) + rltv(k,j)*cps(i,k)/pi2
            enddo
            r_catom(i,j) = r_catom_work(j)
          enddo
        enddo

      end subroutine internal_vector

      subroutine packuc_EDA                                             ! <- subroutine packuc in b_Ldos_f77.F

        integer                                      :: i, j, n

        do j = 1, 3
          do i = 1, natm2
            n = r_catom(i,j)
            r_catom(i,j) = r_catom(i,j) - n

            if(r_catom(i,j) < 0.d0) r_catom(i,j) = r_catom(i,j) + 1.d0
          enddo
        enddo

      end subroutine packuc_EDA

      subroutine determine_grid_points()

        integer                                      :: ix, iy, iz
        real(kind=DP)                                :: dx, dy, dz

        dx = 1.d0/dble(ifxr_CD)
        dy = 1.d0/dble(ifyr_CD)
        dz = 1.d0/dble(ifzr_CD)

        do ix = 1, ifxr_CD+1
          grid_x_CD(ix) = dx * (ix-1)
        enddo

        do iy = 1, ifyr_CD+1
          grid_y_CD(iy) = dy * (iy-1)
        enddo

        do iz = 1, ifzr_CD+1
          grid_z_CD(iz) = dz * (iz-1)
        enddo

      end subroutine determine_grid_points
  end subroutine before_Partitioning_Function

  subroutine set_partition_function()

    real(kind=DP)                                :: total
    integer                                      :: ix, iy, iz, j, mm
    integer                                      :: iix, iixmax, iiy, iiymax, iiz, iizmax &
                                                & , ijx, ijxmin, ijxmax &
                                                & , ijy, ijymin, ijymax &
                                                & , ijz, ijzmin, ijzmax
    integer, parameter                           :: cell_range = 3  ! cell_range must be odd number
    integer                                      :: cell_range2, center
    integer                                      :: next_to_center_x, next_to_center_y, next_to_center_z
    integer                                      :: icell1,  cell1_x,  cell1_y,  cell1_z &
                                                & , icellj,  cellj_x,  cellj_y,  cellj_z
    integer, dimension(natm2,3)                  :: edge_or_not
    integer                                      :: iatm, jatm, jmax
    integer                                      :: key
    integer                                      :: iteration
    real(kind=DP)                                :: ri, rj, rij
    real(kind=DP)                                :: xmuij, zmuij, f2, f4, cutij, cutji
    real(kind=DP)                                :: x_flag, y_flag, z_flag
    real(kind=DP), dimension(natm2,2,2,2)        :: pA
    real(kind=DP), dimension(3)                  :: mesh
    real(kind=DP), dimension(3)                  :: r_catom_wk_i, r_catom_wk_j
    integer                                      :: nsize, ig, nx, ny, nz, ierr
    integer                                      :: id_sname = -1
    call tstatc0_begin('set_partition_function ',id_sname,1)

    call determine_edge_or_not(edge_or_not)
    cell_range2 = cell_range/2
    center = cell_range**3/2 + 1

    next_to_center_x = cell_range2 + 2 + (cell_range2    )*cell_range + (cell_range2    )*cell_range**2
    next_to_center_y = cell_range2 + 1 + (cell_range2 + 1)*cell_range + (cell_range2    )*cell_range**2
    next_to_center_z = cell_range2 + 1 + (cell_range2    )*cell_range + (cell_range2 + 1)*cell_range**2

    mesh = 0.d0
    nx = ifxr_CD+1
    ny = ifyr_CD+1
    nz = ifzr_CD+1
    nsize = nx*ny*nz
    partition_function = 0.d0
    do ig=1, nsize
      if(mod(ig,npes)/=mype) cycle
      iz = (ig-1)/(nx*ny)+1
      mm = mod(ig,nx*ny)
      if(mm==0) mm=nx*ny
      iy = (mm-1)/nx+1
      ix = mod(mm,nx)
      if(ix==0) ix=nx

!    do iz = 1, ifzr_CD+1
      z_flag = 1.d0
      if(iz == 1 .or. iz == ifzr_CD+1) z_flag = 0.5d0

!    do iy = 1, ifyr_CD+1
      y_flag = 1.d0
      if(iy == 1 .or. iy == ifyr_CD+1) y_flag = 0.5d0

!    do ix = 1, ifxr_CD+1
      x_flag = 1.d0
      if(ix == 1 .or. ix == ifxr_CD+1) x_flag = 0.5d0

!-----------------------------------------------------------------------
!     grids in the center cell
!-----------------------------------------------------------------------
      do icell1 = 1, cell_range**3
        call flush(0)
        call cell_xyz(cell_range,icell1,cell1_x,cell1_y,cell1_z)

        mesh(1) = grid_x_CD(ix) + dble(cell1_x)
        mesh(2) = grid_y_CD(iy) + dble(cell1_y)
        mesh(3) = grid_z_CD(iz) + dble(cell1_z)

! ----- Setting atom i ( i must be in the center cell.)

        call set_pA(pA,edge_or_not)

        IatomLoop: do iatm = 1, natm2

          if(abs(aij(iatm,1)+1.d0) < 1.d-05) then
            pA(iatm,:,:,:) = 0.d0
            cycle IatomLoop
          endif

! ----- Distance between the grid point and the center of atom i

          do iiz = 1, edge_or_not(iatm,3)
            r_catom_wk_i(3) = r_catom(iatm,3) + dble(iiz - 1)
          do iiy = 1, edge_or_not(iatm,2)
            r_catom_wk_i(2) = r_catom(iatm,2) + dble(iiy - 1)
          do iix = 1, edge_or_not(iatm,1)
            r_catom_wk_i(1) = r_catom(iatm,1) + dble(iix - 1)

            call set_r(ri,r_catom_wk_i,mesh)

! ----- Setting atom j (i and j do not exist in the same cell.)

          do icellj = 1, cell_range**3
            if(icellj < center) then
              key = 1
              call cell_xyz(cell_range,icellj,cellj_x,cellj_y,cellj_z)
              jmax = natm2
            elseif(icellj == center) then
              key = 2
              cellj_x = 0; cellj_y = 0; cellj_z = 0

              if(product(edge_or_not(iatm,1:3)) == 1) then
                jmax = iatm - 1
              else
                jmax = iatm
              endif
            else
              key = 3
              call cell_xyz(cell_range,icellj,cellj_x,cellj_y,cellj_z)
              jmax = natm2
            endif

            JatomLoop: do jatm = 1, jmax

! ----- Avoid the case that atom j matches the edge atom of center cell

              ijxmin = 1; ijymin = 1; ijzmin = 1

              if(key == 3) then
                if(edge_or_not(jatm,1) == 2 .and. icellj == next_to_center_x) ijxmin = 2
                if(edge_or_not(jatm,2) == 2 .and. icellj == next_to_center_y) ijymin = 2
                if(edge_or_not(jatm,3) == 2 .and. icellj == next_to_center_z) ijzmin = 2
              endif

! ----- Consider edge atom only if the calculated cell is rightest.

              ijxmax = 1; ijymax = 1; ijzmax = 1

              if(key /= 2) then
                if(cellj_x == cell_range2) ijxmax = edge_or_not(jatm,1)
                if(cellj_y == cell_range2) ijymax = edge_or_not(jatm,2)
                if(cellj_z == cell_range2) ijzmax = edge_or_not(jatm,3)
              endif

! ----- Distance between the grid point and the center of atom i

              do ijz = ijzmin, ijzmax
                r_catom_wk_j(3) = r_catom(jatm,3) + dble(cellj_z) + dble(ijz - 1)
              do ijy = ijymin, ijymax
                r_catom_wk_j(2) = r_catom(jatm,2) + dble(cellj_y) + dble(ijy - 1)
              do ijx = ijxmin, ijxmax
                r_catom_wk_j(1) = r_catom(jatm,1) + dble(cellj_x) + dble(ijx - 1)

                if(abs(aij(iatm,jatm)-1.d0) < 1.d-05) then
                  cutij = 1.d0
                  cutji = 1.d0
                  go to 3
                endif

                call set_r(rj,r_catom_wk_j,mesh)
                call set_r(rij,r_catom_wk_i,r_catom_wk_j)  ! <- distance between the centers of atom i and j
                if(rij == 0.d0) cycle

! ----- Becke's fuzzy cell method for molecular grid quadrature

                zmuij = (ri-rj)/rij
                xmuij = zmuij + aij(iatm,jatm)*(1.d0-zmuij*zmuij)

                f4 = xmuij

                do iteration = 1, 4
                  f4 = f4*(1.5d0 - 0.5d0*f4*f4)
                enddo

                f2    = 0.5d0*f4
                cutij = 0.5d0-f2
                cutji = 0.5d0+f2

 3              pA(iatm,iix,iiy,iiz) = pA(iatm,iix,iiy,iiz)*cutij
                if(key == 2) then
                   pA(jatm,ijx,ijy,ijz) = pA(jatm,ijx,ijy,ijz)*cutji
                endif

              enddo
              enddo
              enddo
            enddo JatomLoop

          enddo

          enddo
          enddo
          enddo
        enddo IatomLoop
!        total = 0.d0
        do iatm = 1, natm2
          do iiz = 1, edge_or_not(iatm,3)
            do iiy = 1, edge_or_not(iatm,2)
              do iix = 1, edge_or_not(iatm,1)
                partition_function(ix,iy,iz,iatm) = partition_function(ix,iy,iz,iatm) &
                                                & + pA(iatm,iix,iiy,iiz)/dble(product(edge_or_not(iatm,1:3)))
!                total = total + pA(iatm,iix,iiy,iiz)/dble(product(edge_or_not(iatm,1:3)))
              enddo
            enddo
          enddo
        enddo

      enddo

!      total = 0.d0
      total = sum(partition_function(ix,iy,iz,1:natm2))
      partition_function(ix,iy,iz,:) = x_flag*y_flag*z_flag*partition_function(ix,iy,iz,:)/total

!    enddo
!    enddo
    enddo
    call mpi_allreduce(mpi_in_place, partition_function, nsize*natm2, mpi_double_precision, mpi_sum, MPI_CommGroup, ierr)
    call tstatc0_end(id_sname)
  end subroutine set_partition_function

  subroutine determine_edge_or_not(edge_or_not)

    integer, intent(out), dimension(natm2,3)   :: edge_or_not
    integer                                    :: iatm1, j

    do iatm1 = 1, natm2
      do j = 1, 3
        if(abs(r_catom(iatm1,j)) < 1.d-8) then
          edge_or_not(iatm1,j) = 2
        else
          edge_or_not(iatm1,j) = 1
        endif
      enddo
    enddo

  end subroutine determine_edge_or_not

  subroutine cell_xyz(cell_range,icell,cellx,celly,cellz)

    integer, intent(in)                          :: cell_range, icell
    integer, intent(out)                         :: cellx, celly, cellz
    integer                                      :: half_range
    integer                                      :: i0, i1, i2, i3, j, k

    half_range = cell_range/2 + 1

    i0 = mod(icell,cell_range)
    if(i0 .lt. 1) i0 = cell_range
    i1 = i0 - half_range

    j = (icell - i1)/cell_range
    i2 = mod(j,cell_range) + 1 - half_range

    k = (icell - i1 - (i2-1)*cell_range)/cell_range**2
    i3 = k + 1 - half_range

    cellx = i1; celly = i2; cellz = i3

  end subroutine cell_xyz

  subroutine set_pA(pA,edge_or_not)

    real(kind=DP), intent(out)                   :: pA(natm2,2,2,2)
    integer, intent(in)                          :: edge_or_not(natm2,3)
    integer                                      :: iatm1

    pA = 1.d0
    do iatm1 = 1, natm2
      if(edge_or_not(iatm1,1) == 2) then
        pA(iatm1,2,:,:) = pA(iatm1,2,:,:) * 1.d0
      else
        pA(iatm1,2,:,:) = pA(iatm1,2,:,:) * 0.d0
      endif
      if(edge_or_not(iatm1,2) == 2) then
        pA(iatm1,:,2,:) = pA(iatm1,:,2,:) * 1.d0
      else
        pA(iatm1,:,2,:) = pA(iatm1,:,2,:) * 0.d0
      endif
      if(edge_or_not(iatm1,3) == 2) then
        pA(iatm1,:,:,2) = pA(iatm1,:,:,2) * 1.d0
      else
        pA(iatm1,:,:,2) = pA(iatm1,:,:,2) * 0.d0
      endif
    enddo

  end subroutine set_pA

  subroutine set_r(r,point_1,point_2)

    real(kind=DP), intent(out)                   :: r
    real(kind=DP), intent(in)                    :: point_1(3), point_2(3)
    integer                                      :: j, k
    real(kind=DP), dimension(3)                  :: cps_1(3), cps_2(3)
    real(kind=DP)                                :: r_sqrd
    real(kind=DP)                                :: pi2

    pi2 = 8.d0*datan(1.d0)
    cps_1 = 0.d0; cps_2 = 0.d0

    do j = 1, 3
      do k = 1, 3
        cps_1(j) = cps_1(j) + pi2*point_1(k)*rltv_inverse(k,j)
        cps_2(j) = cps_2(j) + pi2*point_2(k)*rltv_inverse(k,j)
      enddo
    enddo

    r_sqrd = (cps_1(1) - cps_2(1))**2 &
         & + (cps_1(2) - cps_2(2))**2 &
         & + (cps_1(3) - cps_2(3))**2

    r = dsqrt(r_sqrd)

  end subroutine set_r

!------------------------------------------------------------------------------
!   >>>>> calculating each energy terms
!------------------------------------------------------------------------------

  subroutine Each_Energy_Term

!------------------------------------------------------------------------------
!   Etotal = Ekin + Ehartree + Exc + Elocal + Enonlocal
!
!   energies obtained in G-space: kinetic energy
!                                 hartree(coulomb) energy
!                                 local energy
!     - these terms are calculated in G-space and then
!       transformed into R-space expression by FFT.
!
!   energies obtained in R-space: exchange-correlation energy
!
!------------------------------------------------------------------------------

    allocate(density_per_atom(natm2)); density_per_atom = 0.d0
    call calculate_density_per_atom

    allocate(ekin_per_atom(natm2)); ekin_per_atom = 0.d0
    call calculate_kinetic_energy_per_atom

    allocate(ehartr_per_atom(natm2)); ehartr_per_atom = 0.d0
    call calculate_hartree_energy_per_atom

    allocate(exc_per_atom(natm2)); exc_per_atom = 0.d0
    call calculate_exchange_correlation_energy_per_atom

    allocate(elocal_per_atom(natm2)); elocal_per_atom = 0.d0
    call calculate_local_energy_per_atom

    allocate(enonlocal_per_atom(natm2)); enonlocal_per_atom = 0.d0
    call calculate_nonlocal_energy_per_atom

    allocate(epcx_per_atom(natm2)); epcx_per_atom = 0.d0
    call calculate_partial_core_xc_energy_per_atom

    allocate(etotal_per_atom(natm2)); etotal_per_atom = 0.d0
    call calculate_total_energy_per_atom

 contains
    subroutine calculate_density_per_atom

      integer                                    :: ispin, ib, ik, ifft, ri, iend, igp, ip, i1, iloop
      integer                                    :: iatm, ix, iy, iz
      real(kind=DP)                              :: fac, total
      integer                                    :: ierr
      integer                                    :: id_sname = -1
      call tstatc0_begin('calculate_density_per_atom ',id_sname,1)

!     allocate(density_on_a_grid(nfftp)); density_on_a_grid = 0.d0
        afft = 0.d0
!     call m_FFT_alloc_CD_box()

      do ispin = 1, nspin
        do ri = 1, kimg
          do igp = ista_kngp, iend_kngp
            ip = (igfp_nonpara(igp)-1)*kimg + ri
            afft(ip) = afft(ip) &
                                & + univol*chgq_l(igp,ri,ispin)
          enddo
        enddo
      enddo

      call mpi_allreduce(mpi_in_place, afft, nfftp_nonpara, mpi_double_precision &
      &                 ,mpi_sum, mpi_chg_world, ierr)
      call m_FFT_CD_inverse0(nfout,afft)
      total = 0.d0
      do igp = 1, nfftp_nonpara
        total = total + afft(igp)
      enddo

      call partition_into_atoms(afft,density_per_atom)

!     deallocate(density_on_a_grid)
      call m_FFT_dealloc_CD_box()

      call tstatc0_end(id_sname)

    end subroutine calculate_density_per_atom

    subroutine calculate_kinetic_energy_per_atom
      integer                                    :: ispin, ib, ik, ifft, ri, iend, igp, i, i1, ierr
      integer                                    :: ig
      real(kind=DP)                              :: fac
      integer                                    :: nfft
      real(kind=DP), allocatable, dimension(:,:) :: zajtmp
      integer(kind=DP)                           :: plan
      integer, parameter                         :: FFTW_MEASURE=0
      integer                                    :: id_sname = -1
      call tstatc0_begin('calculate_kinetic_energy_per_atom ',id_sname,1)

      allocate(zajtmp(kg1,kimg))
      call m_FFT_alloc_WF_work()

      nfft = product(fft_box_size_WF(1:3,0)) * kimg
      allocate(wf_EDA(nfft))
      allocate(wf_kin_EDA(nfft))

      allocate(ekin_on_a_grid(nfft)); ekin_on_a_grid = 0.d0
      allocate(ekin_l(ista_kngp:iend_kngp,kimg,nspin)); ekin_l = 0.d0

!-----------------------------------------------------------------------
!               1                i
!     Ekin(G) = -*sigma[|k+G|^2*c  ^2]
!               2  i,k           k+G
!-----------------------------------------------------------------------
      do ispin = 1, nspin, af+1
         ekin_on_a_grid = 0.0d0             ! added by ASMS 2024/03/19

         do ik = ispin, kv3+ispin-nspin, nspin
            if(map_k(ik) /= myrank_k) cycle                               ! MPI
            do ib = 1, neg
               if(map_e(ib) /= myrank_e) cycle
               !if(map_e(ib) /= myrank_e) cycle
               !if(map_ek(ib,ik) /= myrank_ke) cycle
               zajtmp = 0.d0
               do ig=ista_g1k(ik),iend_g1k(ik)
                  zajtmp(ig,:) = zaj_l(ig-ista_g1k(ik)+1,map_z(ib),ik,:)
               enddo
               call mpi_allreduce(mpi_in_place, zajtmp, kg1*kimg, mpi_double_precision, mpi_sum, mpi_ke_world, ierr)
               call map_zaj_to_fft_box(ik, kg1, nfft, zajtmp, wf_EDA)
               call m_FFT_WF(ELECTRON, nfout, wf_EDA, INVERSE, ON)
               call m_ES_WF_kin_in_Rspace(ik,ib,wf_kin_EDA,zajtmp)                  ! <- |k+G|^2*exp(-i(k+G)*r)
               do ifft = 1, nfft-1, 2
                  ekin_on_a_grid(ifft) = ekin_on_a_grid(ifft) &
                       &    + occup_l(map_z(ib),ik) &
                       &    *(wf_EDA(ifft)*wf_kin_EDA(ifft)+wf_EDA(ifft+1)*wf_kin_EDA(ifft+1))
               enddo
            enddo
         enddo
         
         call mpi_allreduce(mpi_in_place, ekin_on_a_grid, nfft, mpi_double_precision, &
              &                  mpi_sum, mpi_kg_world, ierr)
         call mpi_allreduce(mpi_in_place, ekin_on_a_grid, nfft, mpi_double_precision, &
              &                  mpi_sum, mpi_ge_world, ierr)
         call m_FFT_WF(ELECTRON,nfout,ekin_on_a_grid,DIRECT,OFF)
         
         fac = 2.d0/(kv3*product(fft_box_size_WF(1:3,1)))
         
         do ri = 1, kimg
            iend = iend_kngp
            if( iend_kngp > kg ) iend = kg
            if( ista_kngp <= iend ) then
               do igp = ista_kngp, iend                                    ! MPI
                  i1 = kimg*igf(igp) + (ri - kimg)
                  ekin_l(igp,ri,ispin) = ekin_on_a_grid(i1)*fac             ! <- Ekin(G)
               enddo
            endif
         enddo
      enddo
      
      call m_FFT_dealloc_WF_work()
      deallocate(wf_EDA); deallocate(wf_kin_EDA)
      deallocate(zajtmp)

      deallocate(ekin_on_a_grid)
!      allocate(ekin_on_a_grid(nfftp)); ekin_on_a_grid = 0.d0
      allocate(ekin_on_a_grid(nfftp_nonpara)); ekin_on_a_grid = 0.d0
      call m_FFT_alloc_CD_box()

      ekin_on_a_grid = 0.d0
      do ispin = 1, nspin
         do ri = 1, kimg
            do igp = ista_kngp, iend_kngp
               i1 = kimg*igfp_nonpara(igp) + (ri - kimg)
               ekin_on_a_grid(i1) = ekin_on_a_grid(i1) + ekin_l(igp,ri,ispin)
            enddo
         enddo
      enddo
      call mpi_allreduce(mpi_in_place, ekin_on_a_grid, nfftp_nonpara, mpi_double_precision, mpi_sum, &
      &                  mpi_ke_world, ierr)

      call m_FFT_CD_inverse0(nfout,ekin_on_a_grid)                    ! <- Ekin(G) --> Ekin(r)
      call partition_into_atoms(ekin_on_a_grid,ekin_per_atom)

      deallocate(ekin_l)
      deallocate(ekin_on_a_grid)
      call m_FFT_dealloc_CD_box()
      call tstatc0_end(id_sname)


    end subroutine calculate_kinetic_energy_per_atom

    subroutine map_zaj_to_fft_box(ik,kg1,nfft,psi_l,bfft)
      integer, intent(in):: ik,kg1,nfft
      real(kind=DP), dimension(kg1,kimg), intent(in) :: psi_l
      real(kind=DP), dimension(nfft), intent(out) :: bfft
      integer :: i,i1,ri, j, i2, ii
      integer :: id_sname2=-1
      call tstatc0_begin('map_zaj_to_fft_box ',id_sname2)

      bfft = 0.d0
      if(k_symmetry(ik) == GAMMA) then
         if(kimg == 1) then
            i1 = igf(1)
            bfft(i1) = psi_l(1,1)
#ifdef NEC_TUNE_SMP
!CDIR NODEP
#endif
            do ii = 2, iba(ik)
!!$             i = nbase(ii,1)
               i = nbase(ii,ik)
               i1 = igf(i)
               bfft(i1) = psi_l(ii,1)
               j = nbase_gamma(ii,2)
               i2 = igf(j)
               bfft(i2) =   psi_l(ii,1)
            end do
         else if(kimg == 2) then
            i1 = 2*igf(1) - 1
            bfft(i1)   = psi_l(1,1)
            bfft(i1+1) = psi_l(1,2)
#ifdef NEC_TUNE_SMP
!CDIR NODEP
#endif
            do ii = 2, iba(ik)
!!$             i = nbase(ii,1)
               i = nbase(ii,ik)
               i1 = 2*igf(i)-1
               bfft(i1  ) = psi_l(ii,1)
               bfft(i1+1) = psi_l(ii,2)
               j = nbase_gamma(ii,2)
               i2 = 2*igf(j)-1
               bfft(i2  ) = psi_l(ii,1)
               bfft(i2+1) = -psi_l(ii,2)
            end do
         end if
      else
#ifdef NEC_TUNE_SMP
!CDIR NOLOOPCHG
#endif
         do ri = 1, kimg
#ifdef NEC_TUNE_SMP
!CDIR NODEP
#endif
            do i = 1, iba(ik)
               i1 = kimg*igf(nbase(i,ik)) + (ri - kimg)
               bfft(i1) = psi_l(i,ri)   ! MPI
            end do
         end do
      end if
    end subroutine map_zaj_to_fft_box

    subroutine calculate_hartree_energy_per_atom

!-----------------------------------------------------------------------
!                   4pi            rho(G)*rho(-G)
!     Ehartree(G) = --- * V     *  --------------   (G /= 0)
!                    2     cell         |G|^2
!-----------------------------------------------------------------------

      integer                                    :: ispin, ri, igp, ip, ierr
      integer, parameter                         :: UP = 1, DOWN = 2
      real(kind=DP)                              :: pi
      integer                                    :: id_sname = -1
      call tstatc0_begin('calculate_hartree_energy_per_atom ',id_sname,1)
      allocate(density_on_a_grid(nfftp_nonpara)); density_on_a_grid = 0.d0
      allocate(ehartr_on_a_grid2(nfftp_nonpara)); ehartr_on_a_grid2 = 0.d0
      call m_FFT_alloc_CD_box()

      do ri = 1, kimg
        do igp = ista_kngp, iend_kngp
          ip = (igfp_nonpara(igp)-1)*kimg + ri
          if(nspin == 1) then
            density_on_a_grid(ip) = density_on_a_grid(ip) + chgq_l(igp,ri,UP)
            if(igp == 1) cycle
            ehartr_on_a_grid2(ip) = ehartr_on_a_grid2(ip) + chgq_l(igp,ri,UP)/gr_l(igp)**2
          elseif(nspin == 2) then
            density_on_a_grid(ip) = density_on_a_grid(ip) + (chgq_l(igp,ri,UP)+chgq_l(igp,ri,DOWN))
            if(igp == 1) cycle
            ehartr_on_a_grid2(ip) = ehartr_on_a_grid2(ip) + (chgq_l(igp,ri,UP)+chgq_l(igp,ri,DOWN))/gr_l(igp)**2
          endif
        enddo
      enddo

      call mpi_allreduce(mpi_in_place, ehartr_on_a_grid2, nfftp_nonpara, mpi_double_precision &
      &                 ,mpi_sum, mpi_chg_world, ierr)
      call mpi_allreduce(mpi_in_place, density_on_a_grid, nfftp_nonpara, mpi_double_precision &
      &                 ,mpi_sum, mpi_chg_world, ierr)

      call m_FFT_CD_inverse0(nfout,density_on_a_grid)
      call m_FFT_CD_inverse0(nfout,ehartr_on_a_grid2)
!      call m_FFT_CD_inverse_ascat(nfout,ehartr_on_a_grid2)
!     call m_FFT_CD_direct0(nfout,ehartr_on_a_grid2)

      allocate(ehartr_on_a_grid(nfftp_nonpara)); ehartr_on_a_grid = 0.d0

      pi = 4.d0*datan(1.d0)

      do ip = 1, nfftp_nonpara
        ehartr_on_a_grid(ip) = 4.d0*pi*univol*0.5d0*density_on_a_grid(ip)*ehartr_on_a_grid2(ip)
      enddo

      call partition_into_atoms(ehartr_on_a_grid,ehartr_per_atom)

      deallocate(ehartr_on_a_grid)
      deallocate(ehartr_on_a_grid2)
      deallocate(density_on_a_grid)
      call m_FFT_dealloc_CD_box()
      call tstatc0_end(id_sname)
    end subroutine calculate_hartree_energy_per_atom

    subroutine calculate_exchange_correlation_energy_per_atom
      integer :: id_sname = -1
      call tstatc0_begin('calculate_exchange_correlation_energy_per_atom ',id_sname,1)
      call m_XC_cal_potential(nfout,Valence_plus_PC_Charge,chgq_l,EXC_ONLY)
      call partition_into_atoms(exc_on_a_grid,exc_per_atom)                 ! exc_on_a_grid -> m_XC_Potential
      call tstatc0_end(id_sname)
    end subroutine calculate_exchange_correlation_energy_per_atom

    subroutine calculate_local_energy_per_atom

!-----------------------------------------------------------------------
!                          atom                     local
!   Elocal = V     * sigma sigma rho(G)*exp(-iGR )*V(G)   + totch*etot1
!             cell   G/=0    A                  A   A
!-----------------------------------------------------------------------

      integer                                    :: ri, ip, igp, ispin, it, iatm, i, ierr
      integer, parameter                         :: UP = 1, DOWN = 2
      integer                                    :: id_sname = -1

      call tstatc0_begin('calculate_local_energy_per_atom ',id_sname,1)

!-----------------------------------------------------------------------
!     partition totch*etot1
!-----------------------------------------------------------------------

      allocate(elocal_on_a_grid1(nfftp_nonpara)); elocal_on_a_grid1 = 0.d0
      call m_FFT_alloc_CD_box()

      do ispin = 1, nspin
        do ri = 1, kimg
          do igp = ista_kngp, iend_kngp
            ip = (igfp_nonpara(igp)-1)*kimg + ri
            elocal_on_a_grid1(ip) = elocal_on_a_grid1(ip) & ! from electron
                                & + 0.5d0*etot1*univol*chgq_l(igp,ri,ispin)
                                                      !~~~~~~~~~~~~~~~~~~~~~~~~
          enddo
        enddo
      enddo

      call mpi_allreduce(mpi_in_place, elocal_on_a_grid1, nfftp_nonpara, mpi_double_precision &
      &                 ,mpi_sum, mpi_chg_world, ierr)
      call m_FFT_CD_inverse0(nfout,elocal_on_a_grid1)

      allocate(elocal_from_nucleus1(natm2)); elocal_from_nucleus1 = 0.d0

      do iatm = 1, natm2
        elocal_from_nucleus1(iatm) = elocal_from_nucleus1(iatm) & ! from nuclei
                                 & + 0.5d0*totch*PP_per_atom(ityp(iatm))/dble(icount_EDA(ityp(iatm)))
!                                               ~~~~~~~~~~~~~~~~~~~~~~~~
      enddo

!-----------------------------------------------------------------------
!     partition main part
!-----------------------------------------------------------------------

      allocate(density_on_a_grid(nfftp_nonpara)); density_on_a_grid = 0.d0

      do ispin = 1, nspin
        do ri = 1, kimg
          do igp = ista_kngp, iend_kngp
            ip = (igfp_nonpara(igp)-1)*kimg + ri
            density_on_a_grid(ip) = density_on_a_grid(ip) + chgq_l(igp,ri,ispin)
          enddo
        enddo
      enddo
      call mpi_allreduce(mpi_in_place, density_on_a_grid, nfftp_nonpara, mpi_double_precision &
      &                 ,mpi_sum, mpi_chg_world, ierr)

      call m_FFT_CD_inverse0(nfout,density_on_a_grid)

      allocate(PP_on_a_grid_wk(nfftp_nonpara)); PP_on_a_grid_wk = 0.d0

! ---> ASMS mod 2024/03/19
!!!      allocate(PP_on_a_grid(nfftp_nonpara,natm2)); PP_on_a_grid = 0.d0
      allocate(elocal_on_a_grid2(nfftp)); elocal_on_a_grid2 = 0.d0
      allocate(elocal_from_nucleus2(natm2)); elocal_from_nucleus2 = 0.d0
! <--- ASMS mod 2024/03/19

      do iatm = 1, natm2
        it = ityp(iatm)
        PP_on_a_grid_wk = 0.d0
        do ri = 1, kimg
          do igp = ista_kngp, iend_kngp
            ip = (igfp_nonpara(igp)-1)*kimg + ri
           !PP_on_a_grid_wk(ip) = PP_on_a_grid_wk(ip) + (-1.d0)**(ri-1)*psc_l(igp,it)*zfm3_l_EDA(igp,iatm,ri)
            PP_on_a_grid_wk(ip) = PP_on_a_grid_wk(ip) + psc_l(igp,it)*zfm3_l_EDA(igp,iatm,ri)
          enddo
        enddo
!        call m_FFT_CD_inverse_ascat(nfout,PP_on_a_grid_wk)
        call mpi_allreduce(mpi_in_place, PP_on_a_grid_wk, nfftp_nonpara, mpi_double_precision &
        &                 ,mpi_sum, mpi_chg_world, ierr)
        call m_FFT_CD_inverse0(nfout,PP_on_a_grid_wk)

! ---> ASMS mod 2024/03/19 
!!        PP_on_a_grid(:,iatm) = PP_on_a_grid_wk(:)
        do i = 1, nfftp_nonpara
           elocal_on_a_grid2(i) = elocal_on_a_grid2(i) &
                & + 0.5d0*univol*density_on_a_grid(i)*PP_on_a_grid_wk(i)
           elocal_from_nucleus2(iatm) = elocal_from_nucleus2(iatm) &
                &           + 0.5d0*univol*density_on_a_grid(i)*PP_on_a_grid_wk(i)
        enddo
! <--- ASMS mod 2024/03/19 

      enddo

      deallocate(PP_on_a_grid_wk)

! ---> ASMS mod 2024/03/19 
!      allocate(elocal_on_a_grid2(nfftp_nonpara)); elocal_on_a_grid2 = 0.d0
!      allocate(elocal_from_nucleus2(natm2)); elocal_from_nucleus2 = 0.d0
!
!      do iatm = 1, natm
!        do i = 1, nfftp_nonpara
!          elocal_on_a_grid2(i) = elocal_on_a_grid2(i) &
!                             & + 0.5d0*univol*density_on_a_grid(i)*PP_on_a_grid(i,iatm)
!          elocal_from_nucleus2(iatm) = elocal_from_nucleus2(iatm) &
!                                   & + 0.5d0*univol*density_on_a_grid(i)*PP_on_a_grid(i,iatm)
!        enddo
!      enddo
!      deallocate(PP_on_a_grid)
! <--- ASMS mod 2024/03/19 

      deallocate(density_on_a_grid)
      allocate(elocal_on_a_grid(nfftp_nonpara)); elocal_on_a_grid = 0.d0
      elocal_on_a_grid = elocal_on_a_grid1 + elocal_on_a_grid2

      deallocate(elocal_on_a_grid1)
      deallocate(elocal_on_a_grid2)
      allocate(elocal_per_atom_from_electron(natm2)); elocal_per_atom_from_electron = 0.d0

      call partition_into_atoms(elocal_on_a_grid,elocal_per_atom_from_electron)

      elocal_per_atom = elocal_per_atom_from_electron + elocal_from_nucleus1 + &
      &                 elocal_from_nucleus2/dble(product(fft_box_size_CD_nonpara(1:3,0)))

      deallocate(elocal_per_atom_from_electron)
      deallocate(elocal_on_a_grid)
      deallocate(elocal_from_nucleus1)
      deallocate(elocal_from_nucleus2)
      call m_FFT_dealloc_CD_box()
      call tstatc0_end(id_sname)

    end subroutine calculate_local_energy_per_atom

    subroutine calculate_nonlocal_energy_per_atom

      integer                                    :: ispin, it,lmt1, lmt2, il1, im1, il2, im2, iatm
      real(kind=DP)                              :: fac
      integer                                    :: id_sname = -1

!-----------------------------------------------------------------------
!                 natom              A,ion  occ   A,i*    A,i
!     Enonlocal = sigma sigma sigma D      sigma f     * f
!                   A     lm   tt'   tt'l   i,k   tlmk    t'lmk
!
!          ion
!       - D    , f are defined by Vanderbilt's ultrasoft potential
!
!-----------------------------------------------------------------------

      call tstatc0_begin('calculate_nonlocal_energy_per_atom ',id_sname,1)

      do ispin = 1, nspin, af+1
        do it = 1, ntyp
          do lmt1 = 1, ilmt(it)
            il1 = ltp(lmt1,it); im1 = mtp(lmt1,it)
            do lmt2 = lmt1, ilmt(it)
              il2 = ltp(lmt2,it); im2 = mtp(lmt2,it)
              if( il1 /= il2 .or. im1 /= im2) cycle
              do iatm = 1, natm
                if(ityp(iatm) /= it) cycle
                fac = 2.d0*iwei(iatm); if(lmt1 == lmt2) fac = fac*0.5d0
                enonlocal_per_atom(iatm) = enonlocal_per_atom(iatm) &
                                       & + fac*dion(lmt1,lmt2,it)*hsr(iatm,lmt1,lmt2,ispin)
              end do
            end do
          end do
        end do
      end do

      enonlocal_per_atom = enonlocal_per_atom*(af+1)
      call tstatc0_end(id_sname)
    end subroutine calculate_nonlocal_energy_per_atom

    subroutine calculate_partial_core_xc_energy_per_atom

      integer                                    :: iatm
      integer :: id_sname = -1
      call tstatc0_begin('calculate_partial_core_xc_energy_per_atom ',id_sname,1)

      do iatm = 1, natm2
        epcx_per_atom(iatm) = -epc_per_atom(ityp(iatm))/dble(icount_EDA(ityp(iatm)))
      enddo
      call tstatc0_end(id_sname)
    end subroutine calculate_partial_core_xc_energy_per_atom

    subroutine calculate_total_energy_per_atom

      integer                                    :: iatm
      integer                                    :: id_sname = -1
      call tstatc0_begin('calculate_total_energy_per_atom ',id_sname,1)

      do iatm = 1, natm2
        etotal_per_atom(iatm) = etotal_per_atom(iatm)    &
                            & + ekin_per_atom(iatm)      &
                            & + ehartr_per_atom(iatm)    &
                            & + exc_per_atom(iatm)       &
                            & + elocal_per_atom(iatm)    &
                            & + enonlocal_per_atom(iatm) &
                            & + eewald_per_atom(iatm)    &
                            & + epcx_per_atom(iatm)
      enddo
      call tstatc0_end(id_sname)

    end subroutine calculate_total_energy_per_atom

  end subroutine Each_Energy_Term

!------------------------------------------------------------------------------
!   >>>>> partitioning subroutines >>>>>
!------------------------------------------------------------------------------

    subroutine partition_into_atoms(afft_on_a_grid,afft_per_atom)

      integer                                    :: id, mm, nl, nm, nn, nlhf, inew, jnew, knew
      integer                                    :: iatm, ix, iy, iz, ix2, iy2, iz2, ixnew, iynew, iznew, i1
      real(kind=DP),intent(in),dimension(nfftp_nonpara)  :: afft_on_a_grid
      real(kind=DP),intent(out),dimension(natm2) :: afft_per_atom
      integer :: id_sname = -1
      call tstatc0_begin('partition_into_atoms ',id_sname,1)

      id = fft_box_size_CD_nonpara(1,0)
      mm = fft_box_size_CD_nonpara(2,0)
      nl = fft_box_size_CD_nonpara(1,0)
      nm = fft_box_size_CD_nonpara(2,0)
      nn = fft_box_size_CD_nonpara(3,0)

      nlhf = id

      do iz2 = 1, ifzr_CD+1                                              ! ifzr_CD = nn
        iz = iz2
        if(iz2 == ifzr_CD+1) iz = 1

        do iy2 = 1, ifyr_CD+1                                            ! ifyr_CD = nm
          iy = iy2
          if(iy2 == ifyr_CD+1) iy = 1

          do ix2 = 1, ifxr_CD+1                                          ! ifxr_CD = nl
            ix = ix2
            if(ix2 == ifxr_CD+1) ix = 1

            if(kimg == 1 .and. ix > nlhf) then
              ixnew = id - ix
              iynew = nm + 2 - iy
              iznew = nn + 2 - iz

              if(iynew > nm) then
                iynew = iynew - nm
              endif
              if(iznew > nn) then
                iznew = iznew - nn
              endif

            else
              ixnew = ix; iynew = iy; iznew = iz
            end if

            i1 = ixnew + nlhf*mm*(iznew-1) + nlhf*(iynew-1)

            do iatm = 1, natm2
              afft_per_atom(iatm) = afft_per_atom(iatm) &
                                & + grid_weight_CD*partition_function(ix2,iy2,iz2,iatm) &
                                & * (afft_on_a_grid(i1*2-1)+afft_on_a_grid(i1*2))
            enddo
          enddo
        enddo
      enddo
      call tstatc0_end(id_sname)

    end subroutine partition_into_atoms

    subroutine write_energies_per_atom

      integer                                    :: iatm

      write(nfout,'("")')
      write(nfout,'(6x,"atom",8x,"density",8x,"kinetic",8x,"hartree",12x, &
      "exc",10x,"local",7x,"nonlocal",10x,"ewald",12x,"epc",10x,"total")')

      do iatm = 1, natm2
        write(nfout,'(i4,"  ",a4,9f15.7)') iatm, speciesname(ityp(iatm)), &
                                         & density_per_atom(iatm),        &
                                         & ekin_per_atom(iatm),           &
                                         & ehartr_per_atom(iatm),         &
                                         & exc_per_atom(iatm),            &
                                         & elocal_per_atom(iatm),         &
                                         & enonlocal_per_atom(iatm),      &
                                         & eewald_per_atom(iatm),         &
                                         & epcx_per_atom(iatm),           &
                                         & etotal_per_atom(iatm)
      enddo
      write(nfout,'("  total   ",9f15.7)') sum(density_per_atom(1:natm2)),   &
                                         & sum(ekin_per_atom(1:natm2)),      &
                                         & sum(ehartr_per_atom(1:natm2)),    &
                                         & sum(exc_per_atom(1:natm2)),       &
                                         & sum(elocal_per_atom(1:natm2)),    &
                                         & sum(enonlocal_per_atom(1:natm2)), &
                                         & sum(eewald_per_atom(1:natm2)),    &
                                         & sum(epcx_per_atom(1:natm2)),      &
                                         & sum(etotal_per_atom(1:natm2))

      write(nfout,'("")')

    end subroutine write_energies_per_atom

    subroutine inverse_matrix(matrix_in,matrix_out)

      real(kind=DP),intent(in),dimension(3,3)    :: matrix_in
      real(kind=DP),intent(out),dimension(3,3)   :: matrix_out
      real(kind=DP),dimension(3,6,4)             :: matrix_sweep
      real(kind=DP),dimension(3,6,4)             :: matrix_temp
      integer                                    :: i, j

!-----------------------------------------------------------------------
!     Due to symmetry, internal vectors of atomic centers are moved.
!     So inverse matrix of rltv is required  to obtain new real-space
!     coordinates of atomic centers.
!-----------------------------------------------------------------------

! -----  At first, matrix_sweep(:,i,1) (i=4,5,6) is identity matrix.  -----

      do i = 1, 3
        do j = 1, 3
          matrix_sweep(i,j,1) = matrix_in(i,j)
          if(i == j) then
            matrix_sweep(i,j+3,1) = 1.d0
          else
            matrix_sweep(i,j+3,1) = 0.d0
          endif
        enddo
      enddo

! -----  Sweep out matrix_in. In this procedure, matrix_sweeo becomes inverse matrix of matrix_in.  -----

      if(matrix_sweep(1,1,1) /= 0.d0) then
        matrix_sweep(1,:,2) = matrix_sweep(1,:,1)/matrix_sweep(1,1,1)
        matrix_sweep(2,:,2) = matrix_sweep(2,:,1) - matrix_sweep(2,1,1)/matrix_sweep(1,1,1)*matrix_sweep(1,:,1)
        matrix_sweep(3,:,2) = matrix_sweep(3,:,1) - matrix_sweep(3,1,1)/matrix_sweep(1,1,1)*matrix_sweep(1,:,1)

        if(matrix_sweep(2,2,2) /= 0.d0) then
          matrix_sweep(1,:,3) = matrix_sweep(1,:,2) - matrix_sweep(1,2,2)/matrix_sweep(2,2,2)*matrix_sweep(2,:,2)
          matrix_sweep(2,:,3) = matrix_sweep(2,:,2)/matrix_sweep(2,2,2)
          matrix_sweep(3,:,3) = matrix_sweep(3,:,2) - matrix_sweep(3,2,2)/matrix_sweep(2,2,2)*matrix_sweep(2,:,2)

          if(matrix_sweep(3,3,3) /= 0.d0) then
            matrix_sweep(1,:,4) = matrix_sweep(1,:,3) - matrix_sweep(1,3,3)/matrix_sweep(3,3,3)*matrix_sweep(3,:,3)
            matrix_sweep(2,:,4) = matrix_sweep(2,:,3) - matrix_sweep(2,3,3)/matrix_sweep(3,3,3)*matrix_sweep(3,:,3)
            matrix_sweep(3,:,4) = matrix_sweep(3,:,3)/matrix_sweep(3,3,3)

            go to 2
          else
            go to 1
          endif

        elseif(matrix_sweep(3,2,2) /= 0.d0) then

          matrix_temp(2,:,2)  = matrix_sweep(2,:,2)
          matrix_temp(3,:,2)  = matrix_sweep(3,:,2)
          matrix_sweep(2,:,2) = matrix_temp(3,:,2)
          matrix_sweep(3,:,2) = matrix_temp(2,:,2)

          matrix_sweep(1,:,3) = matrix_sweep(1,:,2) - matrix_sweep(1,2,2)/matrix_sweep(2,2,2)*matrix_sweep(2,:,2)
          matrix_sweep(2,:,3) = matrix_sweep(2,:,2)/matrix_sweep(2,2,2)
          matrix_sweep(3,:,3) = matrix_sweep(3,:,2) - matrix_sweep(3,2,2)/matrix_sweep(2,2,2)*matrix_sweep(2,:,2)

          if(matrix_sweep(3,3,3) /= 0.d0) then
            matrix_sweep(1,:,4) = matrix_sweep(1,:,3) - matrix_sweep(1,3,3)/matrix_sweep(3,3,3)*matrix_sweep(3,:,3)
            matrix_sweep(2,:,4) = matrix_sweep(2,:,3) - matrix_sweep(2,3,3)/matrix_sweep(3,3,3)*matrix_sweep(3,:,3)
            matrix_sweep(3,:,4) = matrix_sweep(3,:,3)/matrix_sweep(3,3,3)

            go to 2
          else
            go to 1
          endif
        endif

      elseif(matrix_sweep(2,1,1) /= 0.d0 .and. matrix_sweep(1,2,1) /= 0.d0) then

        matrix_temp(1,:,1)  = matrix_sweep(1,:,1)
        matrix_temp(2,:,1)  = matrix_sweep(2,:,1)
        matrix_sweep(1,:,1) = matrix_temp(2,:,1)
        matrix_sweep(2,:,1) = matrix_temp(1,:,1)

        matrix_sweep(1,:,2) = matrix_sweep(1,:,1)/matrix_sweep(1,1,1)
        matrix_sweep(2,:,2) = matrix_sweep(2,:,1) - matrix_sweep(2,1,1)/matrix_sweep(1,1,1)*matrix_sweep(1,:,1)
        matrix_sweep(3,:,2) = matrix_sweep(3,:,1) - matrix_sweep(3,1,1)/matrix_sweep(1,1,1)*matrix_sweep(1,:,1)

!       if(matrix_sweep(2,2,2) /= 0.d0) then
          matrix_sweep(1,:,3) = matrix_sweep(1,:,2) - matrix_sweep(1,2,2)/matrix_sweep(2,2,2)*matrix_sweep(2,:,2)
          matrix_sweep(2,:,3) = matrix_sweep(2,:,2)/matrix_sweep(2,2,2)
          matrix_sweep(3,:,3) = matrix_sweep(3,:,2) - matrix_sweep(3,2,2)/matrix_sweep(2,2,2)*matrix_sweep(2,:,2)

          if(matrix_sweep(3,3,3) /= 0.d0) then
            matrix_sweep(1,:,4) = matrix_sweep(1,:,3) - matrix_sweep(1,3,3)/matrix_sweep(3,3,3)*matrix_sweep(3,:,3)
            matrix_sweep(2,:,4) = matrix_sweep(2,:,3) - matrix_sweep(2,3,3)/matrix_sweep(3,3,3)*matrix_sweep(3,:,3)
            matrix_sweep(3,:,4) = matrix_sweep(3,:,3)/matrix_sweep(3,3,3)

            go to 2
          else
            go to 1
          endif

      elseif(matrix_sweep(3,1,1) /= 0.d0 .and. matrix_sweep(1,2,1) /= 0.d0) then

        matrix_temp(1,:,1)  = matrix_sweep(1,:,1)
        matrix_temp(3,:,1)  = matrix_sweep(3,:,1)
        matrix_sweep(1,:,1) = matrix_temp(3,:,1)
        matrix_sweep(3,:,1) = matrix_temp(1,:,1)

        matrix_sweep(1,:,2) = matrix_sweep(1,:,1)/matrix_sweep(1,1,1)
        matrix_sweep(2,:,2) = matrix_sweep(2,:,1) - matrix_sweep(2,1,1)/matrix_sweep(1,1,1)*matrix_sweep(1,:,1)
        matrix_sweep(3,:,2) = matrix_sweep(3,:,1) - matrix_sweep(3,1,1)/matrix_sweep(1,1,1)*matrix_sweep(1,:,1)

!       if(matrix_sweep(2,2,2) /= 0.d0) then
          matrix_sweep(1,:,3) = matrix_sweep(1,:,2) - matrix_sweep(1,2,2)/matrix_sweep(2,2,2)*matrix_sweep(2,:,2)
          matrix_sweep(2,:,3) = matrix_sweep(2,:,2)/matrix_sweep(2,2,2)
          matrix_sweep(3,:,3) = matrix_sweep(3,:,2) - matrix_sweep(3,2,2)/matrix_sweep(2,2,2)*matrix_sweep(2,:,2)

          if(matrix_sweep(3,3,3) /= 0.d0) then
            matrix_sweep(1,:,4) = matrix_sweep(1,:,3) - matrix_sweep(1,3,3)/matrix_sweep(3,3,3)*matrix_sweep(3,:,3)
            matrix_sweep(2,:,4) = matrix_sweep(2,:,3) - matrix_sweep(2,3,3)/matrix_sweep(3,3,3)*matrix_sweep(3,:,3)
            matrix_sweep(3,:,4) = matrix_sweep(3,:,3)/matrix_sweep(3,3,3)

            go to 2
          else
            go to 1
          endif

      endif

 1    write(nfout,'("Inverse matrix of this matrix does not exist")')

 2    do i = 1, 3
        matrix_out(:,i) = matrix_sweep(:,i+3,4)
      enddo

    end subroutine inverse_matrix

    subroutine write_headermark_of_PW_EDA

      write(nfout,'("")')
      write(nfout,'("=====   Start of Plane-wave EDA  (by ascat)   =====")')
      write(nfout,'("")')

    end subroutine write_headermark_of_PW_EDA

    subroutine write_endmark_of_PW_EDA

      write(nfout,'("")')
      write(nfout,'("=====   End of Plane-wave EDA  (by ascat)   =====")')
      write(nfout,'("")')

    end subroutine write_endmark_of_PW_EDA
#endif

end module m_PlaneWaveEDA
