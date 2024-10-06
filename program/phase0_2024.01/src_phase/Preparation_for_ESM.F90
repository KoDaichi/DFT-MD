#ifdef ENABLE_ESM_PACK
subroutine Preparation_for_ESM()
    use m_Const_Parameters, only : ON,BARE,PE1,PE2,DP
    use m_Control_Parameters, only : gmaxp,kimg,nspin,sw_esm,esm_z1,esm_izwall &
  &                                , esm_iexpot,esm_z_wall,esm_bar_height,esm_bar_width, esm_e_field &
  &                                , esm_fix_ef,esm_add_elec,esm_gps,esm_gpe,ipriesm,nspin,esm_bc,esm_w &
  &                                , esm_izwall, esm_z_wall, esm_bar_height, esm_bar_width, noncol
    use m_Ionic_System, only : natm,natm2,ityp,iatomn,ntyp,ival,cps
    use m_FFT, only : fft_box_size_CD
    use m_Crystal_Structure, only : altv
    use m_Parallelization, only : npes,mype, ista_kngp, iend_kngp, MPI_CommGroup
    use m_PlaneWaveBasisSet, only : ngabc,igfp_l,kgp
    use mpi

    implicit none
!    include 'mpif.h'
    character(len=3) :: esm_bc_c
    real(kind=DP), allocatable,dimension(:) :: ival_at
    real(kind=DP), allocatable, dimension(:,:) :: cps_tmp
    integer, allocatable, dimension(:,:) :: ngabc_esm
    integer, allocatable, dimension(:) :: igfp_l_esm
    integer :: i,j,ierr
    integer :: nspin_m

    if(sw_esm/=ON) return
    if(esm_bc==BARE)then
       esm_bc_c = 'bc1'
    else if (esm_bc == PE1) then
       esm_bc_c = 'bc2'
    else if (esm_bc == PE2) then
       esm_bc_c = 'bc3'
    endif
    allocate(ival_at(natm));ival_at=0.d0
    do i=1,natm
       ival_at(i) = ival(ityp(i))
    enddo
!    allocate(ngabc_esm(3,1:iend_kngp-ista_kngp+1));ngabc_esm=0
    allocate(ngabc_esm(3,1:kgp));ngabc_esm=0
    do i=1,3
       do j=ista_kngp,iend_kngp
          ngabc_esm(i,j) = ngabc(j,i)
       enddo
    enddo
    call mpi_allreduce(mpi_in_place,ngabc_esm,3*kgp,mpi_integer,mpi_sum,MPI_CommGroup,ierr)
    allocate(cps_tmp(3,natm))
    do i=1,3
       do j=1,natm
          cps_tmp(i,j) = cps(j,i)
       enddo
    enddo
    allocate(igfp_l_esm(1:kgp));igfp_l_esm=0
    do i=ista_kngp,iend_kngp
       igfp_l_esm(i) = igfp_l(i)
    enddo
    call mpi_allreduce(mpi_in_place,igfp_l_esm,kgp,mpi_integer,mpi_sum,MPI_CommGroup,ierr)

    if ( noncol ) then
       nspin_m = 1
    else
       nspin_m = nspin
    endif

    call Esm_interface_map_parameters(natm,ival_at,cps_tmp,1.0d0,altv, &
    & fft_box_size_CD(1,0),fft_box_size_CD(2,0),fft_box_size_CD(3,0),  &
    & esm_bc_c,.false.,kgp,nspin_m,ngabc_esm,                            &
    & igfp_l_esm,igfp_l_esm,esm_w,2.0d0*esm_e_field, esm_izwall,       &
    & esm_z_wall, esm_bar_height,esm_bar_width)
    call Esm_interface_set_communicator(MPI_CommGroup)
    call Esm_interface_set_npes_mype(npes,mype)
    deallocate(igfp_l_esm)
    deallocate(ival_at)
    deallocate(ngabc_esm)
    deallocate(cps_tmp)
end subroutine Preparation_for_ESM
#endif
