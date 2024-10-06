#define INDEX_L_ORDER_CONFLICT
!!$#define LDOShardpartsub3
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 633 $)
!
!  MODULE: m_Ldos
!
!  AUTHOR(S): T. Yamasaki   January/18/2004
!
!  Contact address :  Phase System Consortium
!                     E-mail: phase_system@nims.go.jp URL https://azuma.nims.go.jp
!
!
!
!=======================================================================
!
!     The original version of this set of the computer programs "PHASE"
!  was developed by the members of the Theory Group of Joint Research
!  Center for Atom Technology (JRCAT), based in Tsukuba, in the period
!  1993-2001.
!
!     Since 2002, this set has been tuned and new functions have been
!  added to it as a part of the national project "Frontier Simulation
!  Software for Industrial Science (FSIS)",  which is supported by
!  the IT program of the Ministry of Education, Culture, Sports,
!  Science and Technology (MEXT) of Japan.
!     Since 2006, this program set has been developed as a part of the
!  national project "Revolutionary Simulation Software (RSS21)", which
!  is supported by the next-generation IT program of MEXT of Japan.
!   Since 2013, this program set has been further developed centering on PHASE System
!  Consortium.
!   The activity of development of this program set has been supervised by Takahisa Ohno.
!
module m_Ldos
! $Id: m_Ldos.F90 633 2020-12-01 05:11:03Z jkoga $
  use m_Const_Parameters, only   : DP, CMPLDP, REGULAR_INTERVALS, BY_ATOMIC_POSITIONS &
       &                         , DELTA10, EXECUT, ON, SOFTPART, HARDPART, DIRECT, PAI2 &
       &                         , ALDOS, LAYERDOS, NO, OFF, ELECTRON, BOHR
  use m_Control_Parameters, only : kimg, sw_aldos, sw_layerdos, crtdst_aldos, naldos_from, naldos_to &
       &                         , slicing_way_winlay, deltaz_winlay, normal_axis_winlay &
       &                         , crtdst_winlay,crtdst_is_given, nspin, neg, ipridos,ekmode &
       &                         , af, printable, hardpart_subroutine, sw_rspace_ldos, sw_save_ldos_weight &
       &                         , integration_dimension_winlay, sw_checksum &
       &                         , m_CtrlP_check_naldos_range &
       &                         , sw_ac_mesh, acmesh_factor, ipriparallel, sw_modified_kpoint_increment
  use m_Files, only              : nfout, nfldos, nfwfk_local_decomp &
       &                         , m_Files_open_nfldos, m_Files_close_nfldos
  use m_Crystal_Structure, only :  altv, rltv, univol
  use m_Ionic_System,   only :     natm2,natm,cps,numlay,ityp,iwei,if_aldos,pos
  use m_Parallelization, only :    ista_e,iend_e,istep_e, np_e, map_z, npes, ierr, mype &
       &                         , map_k, myrank_k, map_ek, ista_kngp, iend_kngp, ista_k, iend_k &
       &                         , nrank_e, nrank_k, myrank_e, map_e, mpi_k_world &
       &                         , MPI_CommGroup, ista_atm2, iend_atm2, m_Parallel_init_mpi_atm2 &
       &                         , nis_kv3_ek
  use m_Kpoints, only :            kv3, kv3_ek
  use m_FFT, only :                nfft, nfftp,nfftp_nonpara, fft_box_size_WF, fft_box_size_CD &
       &                         , fft_box_size_CD_nonpara &
       &                         , m_FFT_alloc_CD_box, m_FFT_dealloc_CD_box &
       &                         , m_FFT_CD_inverse0,   m_FFT_WF &
       &                         , m_FFT_CD_inverse &
       &                         , m_FFT_alloc_WF_work, m_FFT_dealloc_WF_work, m_FFT_CD_direct0
  use m_PlaneWaveBasisSet, only :  kg, igf, ngabc, kgp, kg1
  use m_PseudoPotential, only :    modnrm
!$$  use m_Electronic_Structure,only :neordr, totch, efermi &
!$$       &                         , m_ES_WF_in_Rspace &
!$$       &                         , m_ES_wd_zaj_small_portion
  use m_Electronic_Structure,only :neordr, totch, efermi
!$$#ifndef PARA3D
  use m_Electronic_Structure,only :m_ES_WF_in_Rspace, m_ES_WF_in_Rspace_kt &
       &                         , m_ES_wd_zaj_small_portion, zaj_l
!$$#endif
  use m_Charge_Density, only :     chgq_l, chgq_enl &
       &                         , m_CD_hardpart_sub &
       &                         , m_CD_hardpart_sub2 &
       &                         , m_CD_map_chgqenl_to_fft_box &
       &                         , m_CD_map_valence_charge_to_fft_box &
       &                         , m_CD_map_fft_box_to_chgqenl &
       &                         , m_CD_restore_chgq &
       &                         , m_CD_cp_chgq_to_chgqo & !ASMS
       &                         , m_CD_map_chgq_to_fft_box &
       &                         , m_CD_set_ylm_enl_etc &
       &                         , m_CD_dealloc_ylm_enl_etc &
       &                         , m_CD_keep_retrieve_hsr

! ============================== added by K. Tagami ============== 11.0
  use m_Control_Parameters,    only : noncol, ndim_magmom, ndim_spinor, sw_ldos
  use m_Charge_Density, only :    m_CD_hardpart_sub_noncl, &
       &                          m_CD_hardpart_sub2_noncl, &
       &                          m_CD_map_chgqenl_to_fft_box_kt, &
       &                          m_CD_hardpart_sub2_rs
! ================================================================ 11.0

  use m_Control_Parameters,   only : num_extra_bands, sw_rspace_lband, &
       &                             sw_calc_wf_atom_decomposition, &
       &                             sw_calc_wf_layer_decomposition, sw_lband, &
       &                             sw_wd_only_specified_atoms, &
       &                             sw_band_unfolding, band_unfolding_active

  use m_Timing,              only : tstatc0_begin, tstatc0_end
  use mpi

  implicit none

  integer, private, target, allocatable, dimension(:,:,:) ::  mesh, meshp, meshl, meshpl

  integer ::  nlayer = 0
  integer ::  mlayer = 1
  integer ::  naldos = 0
  integer, private :: naldos_write = 0
  integer ::  maldos = 1
  integer, private ::  n_total_ldoscal = 0
  real(kind=DP),private,allocatable,dimension(:,:)  :: winlay ! d(mlayer,2)
  integer, private, allocatable, dimension(:)       :: nmeshlay, nmeshplay ! d(mlayer)
!  real(kind=DP),public,allocatable,dimension(:,:,:,:) :: weiwsc ! d(maldos,neg|np_e,1|ista_k:iend_k,nspin)
!  real(kind=DP),public,allocatable,dimension(:,:,:,:) :: weilay ! d(mlayer,neg|np_e,1|ista_k:iend_k,nspin)
  real(kind=DP),public,allocatable,dimension(:,:,:) :: weiwsc ! d(maldos,neg|np_e,1|ista_k:iend_k)
  real(kind=DP),public,allocatable,dimension(:,:,:) :: weilay ! d(mlayer,neg|np_e,1|ista_k:iend_k)
!  real(kind=DP),private,allocatable,dimension(:) ::    bfft   ! d(nfft)
  real(kind=DP),private :: maxhv, maxv
  real(kind=DP),private :: height
  integer,private,allocatable,dimension(:) :: if_aldos_full ! d(natm2)

  character(len("HARDPART")), private, dimension(2) :: tag_Hard_or_Soft = (/'SOFTPART','HARDPART'/)

! ============================= added by K. Tagami ===================== 11.0
  real(kind=DP),public,allocatable,dimension(:,:,:,:) :: weiwsc_noncl
  real(kind=DP),public,allocatable,dimension(:,:,:,:) :: weilay_noncl
! ====================================================================== 11.0

  integer, parameter :: ilen = 5

  integer, allocatable, target, dimension(:,:,:) :: ac_mesh,ac_mesh_cd
  integer, allocatable, target, dimension(:) :: nac_mesh,nac_mesh_cd

  real(kind=DP), allocatable, dimension(:,:,:,:) :: dos_weight_ek

  integer :: natom_decomp
  integer, allocatable :: atom_decomp_map(:)

!  include 'mpif.h'
  integer istatus(mpi_status_size)

contains
  integer function m_Ldos_what_is_n_total_ldos()
!!$    m_Ldos_what_is_n_total_ldos = n_total_ldos
    m_Ldos_what_is_n_total_ldos = naldos_write + mlayer
  end function m_Ldos_what_is_n_total_ldos

  subroutine m_Ldos_preparation()

    if(sw_aldos == ON .or. sw_calc_wf_atom_decomposition==ON ) then
       call set_naldos()      ! -> naldos, maldos
       call m_CtrlP_check_naldos_range(nfout,maldos)
!!$       call check_naldos_range()
       call set_if_aldos_full()
       call set_naldos_write()
       call alloc_mesh()
       call fillup_mesh(fillmode=ALDOS)
    end if
    if ( sw_calc_wf_atom_decomposition==ON ) call set_natom_decomp()

    call set_nlayer()
    if(sw_layerdos == ON .or. sw_calc_wf_layer_decomposition==ON ) then
       call alloc_winlay()
       call set_winlay()
       if(integration_dimension_winlay==3) then
          call alloc_meshl()
!!$          call fillup_meshl(fillmode=LAYERDOS)
          call fillup_mesh(fillmode=LAYERDOS)
       end if
!!$       call dealloc_winlay()
    end if
  end subroutine m_Ldos_preparation

  subroutine set_naldos()
    naldos = natm2
    maldos = naldos + 1
  end subroutine set_naldos

!!$  subroutine check_naldos_range()
!!$    if(naldos_from == 0 .and. naldos_to == 0) then
!!$       naldos_from = 1
!!$       naldos_to   = maldos
!!$    else
!!$       if(naldos_from < 1)     naldos_from = 1
!!$       if(naldos_from > maldos) naldos_from = maldos
!!$       if(naldos_to < 1)       naldos_to   = 1
!!$       if(naldos_to > maldos)   naldos_to   = maldos
!!$       if(naldos_to < naldos_from) naldos_to = naldos_from
!!$    end if
!!$    write(nfout,'(" !!ldos naldos_from         = ",i6," <<check_naldos_range>>")') naldos_from
!!$    write(nfout,'(" !!ldos naldos_to           = ",i6," <<check_naldos_range>>")') naldos_to
!!$  end subroutine check_naldos_range

  subroutine set_if_aldos_full()
    integer, allocatable, dimension(:) :: ip_atom
    integer :: i, nb

    if(.not.allocated(if_aldos_full)) allocate(if_aldos_full(naldos_from:naldos_to))
    if(.not.allocated(ip_atom)) allocate(ip_atom(natm2+1)); ip_atom = 0

    nb = natm
    do i = 1, natm
       ip_atom(i) = i
       if(iwei(i) == 2) then
          nb = nb + 1
          ip_atom(nb) = i
       end if
    end do

    do i = 1, natm2
       if(ip_atom(i) == 0) then
          if(ipridos>=1) write(nfout,'(" !!ldos i, ip_atom = ",i5,i5)') i,ip_atom(i)
          call phase_error_with_msg(nfout,' ip_atom is illegal <<set_if_aldos_full>>',__LINE__,__FILE__)
       end if
    end do

    if_aldos_full = OFF
    do i = naldos_from, naldos_to
       if(i <= natm2) then
          if_aldos_full(i) = if_aldos(ip_atom(i))
       else
          if_aldos_full(i) = ON
       end if
    end do
    
    nb = 0
    ip_atom = 0
    do i = naldos_from, naldos_to
       if(if_aldos_full(i) == ON) then
          nb = nb + 1
          ip_atom(nb) = i
       end if
    end do
    if(ipridos>=1) write(nfout,'(" !!ldos aldos_atoms = ",10i5)') (ip_atom(i),i=1,nb)

    deallocate(ip_atom)
  end subroutine set_if_aldos_full

  subroutine set_naldos_write()
    integer :: i
    naldos_write = 0
    do i = naldos_from, naldos_to
       if(if_aldos_full(i) == ON) naldos_write = naldos_write + 1
    end do
    if(ipridos>=1) write(nfout,'(" !!ldos naldos_write = ",i5)') naldos_write
!!$    naldos_write = naldos_to - naldos_from + 1
  end subroutine set_naldos_write

  subroutine set_natom_decomp()
    integer :: i

    natom_decomp = 0
    if ( sw_wd_only_specified_atoms == ON ) then
       do i = 1, natm2
          if(if_aldos_full(i) == ON) natom_decomp = natom_decomp +1
       end do
    else
       natom_decomp = natm2
    endif

    allocate( atom_decomp_map(natom_decomp+1) )

    natom_decomp = 0
    if ( sw_wd_only_specified_atoms == ON ) then
       do i = 1, natm2
          if(if_aldos_full(i) == ON) natom_decomp = natom_decomp +1
          atom_decomp_map( natom_decomp ) = i
       end do
    else
       do i = 1, natm2
          natom_decomp = natom_decomp +1
          atom_decomp_map( natom_decomp ) = i
       end do
    endif
    atom_decomp_map( natom_decomp +1 ) = natm2 +1

  end subroutine set_natom_decomp

  subroutine set_nlayer()
    integer :: i,j, io,jo
    real(kind=DP) :: h
    integer, allocatable, dimension(:) :: layer_index, layer_order

!    if(sw_layerdos == ON) then
    if (sw_layerdos == ON .or. sw_calc_wf_layer_decomposition==ON ) then
       call get_height(h)
       height = h
       if(slicing_way_winlay == REGULAR_INTERVALS) then
!!$       call get_height(h)
          if(kimg == 1) h = h*0.5d0
          nlayer = h/deltaz_winlay + (1-DELTA10)
          mlayer = nlayer + 1
       else if(slicing_way_winlay == BY_ATOMIC_POSITIONS) then
          allocate(layer_index(natm))
          nlayer = 1
          layer_index(nlayer) = numlay(1)
          do i = 2, natm
             do j = 1, nlayer
                if(layer_index(j) == numlay(i)) goto 1001
             end do
             nlayer = nlayer + 1
             layer_index(nlayer) = numlay(i)
1001         continue
          end do
          if(minval(layer_index(1:nlayer)) /= 1 .or. maxval(layer_index(1:nlayer)) /= nlayer) then
             allocate(layer_order(nlayer)); layer_order(:) = (/(i,i=1,nlayer)/)
             do i = 1, nlayer-1
                do j = i+1, nlayer
                   io = layer_order(i)
                   jo = layer_order(j)
                   if(layer_index(io) < layer_index(io)) then
                      layer_order(j) = io
                      layer_order(i) = jo
                   end if
                end do
             end do
             if(ipridos>=1) then
                write(nfout,'(" !!ldos i, layer_index, layer_order <<m_Ldos.set_nlayer>>")')
                do i = 1, nlayer
                   write(nfout,'(" !!ldos ",3i8)') i, layer_index(i), layer_order(i)
                end do
             end if
             do i = 1, natm
                do j = 1, nlayer
                   if(layer_index(j) == numlay(i)) then
                      numlay(i) = layer_order(j)
                      goto 1002
                   end if
                end do
1002            continue
             end do
             if(ipridos>=1) then
                write(nfout,'(" !!ldos #atom, numlay <<m_Ldos.set_nlayer>>")')
                do i = 1, natm
                   write(nfout,'(" !!ldos ",3i8)') i, numlay(i)
                end do
             end if
             deallocate(layer_order)
          end if
          deallocate(layer_index)
          mlayer = nlayer + 1
       end if
       if(ipridos>=1) write(nfout,'(" !!ldos nlayer, mlayer = ",2i8)') nlayer, mlayer
    else
       nlayer = 0
       mlayer = 0
    end if
  end subroutine set_nlayer

  subroutine alloc_mesh()
    integer :: id,nd2,nn
    id  = fft_box_size_WF(1,0)
    nd2 = fft_box_size_WF(2,0)
    nn  = fft_box_size_WF(3,1)
    if(.not.allocated(mesh)) allocate(mesh(id,nd2,nn))
    id  = fft_box_size_CD_nonpara(1,0)
    nd2 = fft_box_size_CD_nonpara(2,0)
    nn  = fft_box_size_CD(3,1)
    if(.not.allocated(meshp)) allocate(meshp(id,nd2,nn))
  end subroutine alloc_mesh

  subroutine alloc_meshl()
    integer :: id,nd2,nn
    id  = fft_box_size_WF(1,0)
    nd2 = fft_box_size_WF(2,0)
    nn  = fft_box_size_WF(3,1)
    if(.not.allocated(meshl)) allocate(meshl(id,nd2,nn))
    id  = fft_box_size_CD_nonpara(1,0)
    nd2 = fft_box_size_CD_nonpara(2,0)
    nn  = fft_box_size_CD(3,1)
    if(.not.allocated(meshpl)) allocate(meshpl(id,nd2,nn))
  end subroutine alloc_meshl

  subroutine m_Ldos_dealloc_mesh()
    deallocate(mesh)
    deallocate(meshp)
  end subroutine m_Ldos_dealloc_mesh

  subroutine alloc_winlay()
    if(allocated(winlay)) deallocate(winlay)
    if(allocated(nmeshlay)) deallocate(nmeshlay)
    if(allocated(nmeshplay)) deallocate(nmeshplay)
    allocate(winlay(mlayer,2))
    allocate(nmeshlay(mlayer))
    allocate(nmeshplay(mlayer))
  end subroutine alloc_winlay

  subroutine get_adjustfactor(winlay,m,nmeshl,factor,nwrite)
    integer, intent(in) :: m,nwrite
    real(kind=DP),intent(in), dimension(m,2) :: winlay
    integer, intent(in), dimension(m) :: nmeshl
    real(kind=DP),intent(out), dimension(m) :: factor
    integer :: i,nall
    nall = 0
    do i = 1, m
       nall = nall + nmeshl(i)
    end do
    do i = 1, m
       if(nmeshl(i)>0) then
          factor(i) = nall*(winlay(i,2)-winlay(i,1))/height/nmeshl(i)
       else
          factor(i) = 1.d0
       end if
    end do
    if(nwrite==0) then
       if(ipridos>=1) then
          write(nfout,'("!!ldos  no  delta_window   nmesh  factor  adjustednmesh")')
          do i = 1, m
             write(nfout,'("!!ldos ",i4, f12.5,  i8, f8.4, f16.5)') &
                  & i,winlay(i,2)-winlay(i,1),nmeshl(i),factor(i),nmeshl(i)*factor(i)
          end do
       end if
    end if
  end subroutine get_adjustfactor

  subroutine dealloc_winlay()
    deallocate(nmeshplay,nmeshlay,winlay)
  end subroutine dealloc_winlay

  subroutine fillup_mesh(fillmode)
    integer, intent(in) :: fillmode
    integer :: id, nd2, nl, nm, nn
    real(kind=DP),allocatable,dimension(:,:) :: wk_catoms, cps_full ! d(natm2,3)
    real(kind=DP),allocatable,dimension(:) :: wk_dstnc  ! d(natm2)
    integer, allocatable, dimension(:) ::     wk_ioddst, ityp_full ! d((natm2+1)*4)
    real(kind=DP),allocatable,dimension(:,:) :: winlay_internal    ! d(mlayer,2)
    integer :: nmmax,nmmax0,nsize
    integer :: ia,i,ii,ierr
    integer :: id_sname = -1

    call tstatc0_begin('fillup_mesh ',id_sname)
    call m_Parallel_init_mpi_atm2(nfout,ipriparallel,printable,natm2)
    call get_fftbox_size(SOFTPART,id,nd2,nl,nm,nn) ! fft_box_size_WF -> id,nd2,nl,nm,nn
    if(fillmode==ALDOS) then
       allocate(wk_catoms(natm2,3))
       allocate(wk_dstnc(natm2))
       allocate(wk_ioddst((natm2+1)*4))
       allocate(cps_full(natm2,3))
       allocate(ityp_full(natm2))
       cps_full(1:natm,1:3) = cps(1:natm,1:3)
       ityp_full(1:natm) = ityp(1:natm)
       call rplcps(cps_full,ityp_full,1,natm2,natm,iwei)
       if(sw_ac_mesh==ON)then
         allocate(ac_mesh(natm2,1,3));ac_mesh=0
         allocate(nac_mesh(natm2));nac_mesh=0
         call anlmes_ac(ista_atm2,iend_atm2,nfout,ipridos,1,ac_mesh,nac_mesh,id,nd2 &
              &        ,nl*acmesh_factor,nm*acmesh_factor,nn*acmesh_factor &
              &        ,altv,rltv,cps_full,natm2,crtdst_aldos,DELTA10,wk_catoms &
              &        ,wk_ioddst,wk_dstnc,.true.,nmmax)
         call mpi_allreduce(mpi_in_place,nmmax,1,mpi_integer,mpi_max,MPI_CommGroup,ierr)
         nsize = natm2*nmmax*3
         deallocate(ac_mesh);allocate(ac_mesh(natm2,nmmax,3));ac_mesh=0
         call anlmes_ac(ista_atm2,iend_atm2,nfout,ipridos,nmmax,ac_mesh,nac_mesh,id,nd2 &
              &        ,nl*acmesh_factor,nm*acmesh_factor,nn*acmesh_factor &
              &        ,altv,rltv,cps_full,natm2,crtdst_aldos,DELTA10,wk_catoms &
              &        ,wk_ioddst,wk_dstnc,.false.,nmmax0)
         call mpi_allreduce(mpi_in_place,nac_mesh,natm2,mpi_integer,mpi_sum,MPI_CommGroup,ierr)
         call mpi_allreduce(mpi_in_place,ac_mesh,nsize,mpi_integer,mpi_sum,MPI_CommGroup,ierr)
         write(nfout,'(a)') ' !** AC MESH number of mesh allocated to each atom '
         do ia=1,natm2
            write(nfout,'(2i8)') ia,nac_mesh(ia)
         enddo
         if(ipridos>1)then
           write(nfout,'(a)') ' !** AC MESH AC MESH for each atom '
           do ia=1,natm2
              do i=1,nac_mesh(ia)
                 write(nfout,'(5i8,f20.10)') ia,i,ac_mesh(ia,i,1),ac_mesh(ia,i,2),ac_mesh(ia,i,3)
              enddo
           enddo
         endif
       else
         call anlmes(nfout,ipridos,mesh,id,nd2,nl,nm,nn,altv,rltv,cps_full &
              &           ,natm2,crtdst_aldos,DELTA10,wk_catoms,wk_ioddst,wk_dstnc)
       endif
    else if(fillmode==LAYERDOS) then
       allocate(winlay_internal(mlayer,2))
       winlay_internal = winlay/height
       call anlmesl(nfout,ipridos,meshl,id,nd2,nl,nm,nn,winlay_internal,mlayer,normal_axis_winlay,nmeshlay)
    end if

    call get_fftbox_size(HARDPART,id,nd2,nl,nm,nn) ! fft_box_size_CD -> id,nd2,nl,nm,nn
    if(fillmode==ALDOS) then
       if(sw_ac_mesh==ON) then
          allocate(ac_mesh_cd(natm2,1,3));ac_mesh_cd=0
          allocate(nac_mesh_cd(natm2));nac_mesh_cd=0
          call anlmes_ac(ista_atm2,iend_atm2,nfout,ipridos,1,ac_mesh_cd,nac_mesh_cd,id,nd2 &
               &        ,nl*acmesh_factor,nm*acmesh_factor,nn*acmesh_factor &
               &        ,altv,rltv,cps_full,natm2,crtdst_aldos,DELTA10 &
               &        ,wk_catoms,wk_ioddst,wk_dstnc,.true.,nmmax)
          call mpi_allreduce(mpi_in_place,nmmax,1,mpi_integer,mpi_max,MPI_CommGroup,ierr)
          deallocate(ac_mesh_cd);allocate(ac_mesh_cd(natm2,nmmax,3));ac_mesh_cd=0
          nsize = natm2*nmmax*3
          call anlmes_ac(ista_atm2,iend_atm2,nfout,ipridos,nmmax,ac_mesh_cd,nac_mesh_cd,id,nd2 &
               &        ,nl*acmesh_factor,nm*acmesh_factor,nn*acmesh_factor &
               &        ,altv,rltv,cps_full,natm2,crtdst_aldos,DELTA10 &
               &        ,wk_catoms,wk_ioddst,wk_dstnc,.false.,nmmax0)
          call mpi_allreduce(mpi_in_place,nac_mesh_cd,natm2,mpi_integer,mpi_sum,MPI_CommGroup,ierr)
          call mpi_allreduce(mpi_in_place,ac_mesh_cd,nsize,mpi_integer,mpi_sum,MPI_CommGroup,ierr)
          write(nfout,'(a)') ' !** AC MESH number of mesh allocated to each atom (hardpart)'
          do ia=1,natm2
             write(nfout,'(2i8)') ia,nac_mesh_cd(ia)
          enddo
          if(ipridos>1)then
            write(nfout,'(a)') ' !** AC MESH AC MESH for each atom '
            do ia=1,natm2
               do i=1,nac_mesh_cd(ia)
                  write(nfout,'(5i8,f20.10)') ia,i,ac_mesh_cd(ia,i,1),ac_mesh_cd(ia,i,2),ac_mesh_cd(ia,i,3)
               enddo
            enddo
          endif
       else
          call anlmes(nfout,ipridos-1,meshp,id,nd2,nl,nm,nn,altv,rltv,cps_full &
         &           ,natm2,crtdst_aldos,DELTA10,wk_catoms,wk_ioddst,wk_dstnc)
       endif
       deallocate(ityp_full,cps_full,wk_ioddst,wk_dstnc,wk_catoms)
    else if(fillmode==LAYERDOS) then
       call anlmesl(nfout,ipridos,meshpl,id,nd2,nl,nm,nn,winlay_internal,mlayer,normal_axis_winlay,nmeshplay)
       deallocate(winlay_internal)
    end if
    call tstatc0_end(id_sname)

  end subroutine fillup_mesh

  subroutine get_height(h)
    real(kind=DP), intent(out) :: h
    integer :: iaxis,i2,i3
    real(kind=DP) :: a,b,S

    iaxis = normal_axis_winlay
    i2 = iaxis + 1
    if(i2 >= 4) i2 = i2 - 3
    i3 = i2 + 1
    if(i3 >= 4) i3 = i3 - 3
    a = dsqrt(altv(1,i2)**2 + altv(2,i2)**2 + altv(3,i2)**2)
    b = dsqrt(altv(1,i3)**2 + altv(2,i3)**2 + altv(3,i3)**2)
! ----------------------------
!        29 Jan. 2008
!     Revised by M. Usami
!!$    S = a*a * b*b - (altv(1,i2)*altv(1,i3)+altv(2,i2)*altv(2,i3)+altv(3,i2)*altv(3,i3))
    S = a*a * b*b - (altv(1,i2)*altv(1,i3)+altv(2,i2)*altv(2,i3)+altv(3,i2)*altv(3,i3))**2
! --------------------------
    S = dsqrt(S)
    h = univol/S
    if(ipridos>=1) write(nfout,'(" !!ldos S, h, V = ",3f16.8," <<m_Ldos.set_nlayer>>")') S,h,univol
  end subroutine get_height

  subroutine set_winlay()
! ----------------------------------------------
!       Nov. 17th 1992 by T.Yamasaki
!
!    the original subroutine name was "dwinly2"
!
!       modified by T. Yamasaki, 09th Feb 2004
! ----------------------------------------------
    real(kind=DP) :: zmin, tmp, h
    integer ::  i, nl
    i = normal_axis_winlay
!!$    call get_height(h)
    h = height
    if(kimg == 1) then
       maxhv = h*0.5d0
       maxv  = maxhv*2
    else
       maxhv = h
       maxv  = maxhv
    end if

    winlay = 0.d0
    if(slicing_way_winlay == REGULAR_INTERVALS) then
       if(ipridos>=1) then
          write(nfout,'(" !!ldos slicing_way = REGULAR_INTERVALS")')
          write(nfout,'(" !!ldos nlayer, mlayer = ",2i6)') nlayer, mlayer
          write(nfout,'(" !!ldos maxhv, maxv = ",2f8.4)') maxhv, maxv
          write(nfout,'(" !!deltaz_winlay = ",f8.4)') deltaz_winlay
       end if
       if(kimg == 1) then
          winlay(1,1) = 0.d0
          winlay(mlayer,2) = maxhv
          do i = 1, nlayer
             winlay(i,2) = deltaz_winlay * i
             if(winlay(i,2) > maxhv) winlay(i,2) = maxhv
             if(i+1 <= mlayer) winlay(i+1,1) = winlay(i,2)
          end do
       else if(kimg == 2) then
          zmin = minval(cps(1:natm,normal_axis_winlay))
          winlay(1,1) = zmin - crtdst_winlay
          winlay(mlayer,2) = winlay(1,1) + maxv
          do i = 1, nlayer
             winlay(i,2) = winlay(i,1) + deltaz_winlay
             if(winlay(i,2) > winlay(mlayer,2)) winlay(i,2) = winlay(mlayer,2)
             if(i+1 <= mlayer) winlay(i+1,1) = winlay(i,2)
          end do
       end if
    else if(slicing_way_winlay == BY_ATOMIC_POSITIONS) then
       if(ipridos>=1) &
            & write(nfout,'(" !!ldos slicing_selection = BY_ATOMIC_POSITIONS")')

       do i = 1, mlayer
          winlay(i,1) = +maxv*30
          winlay(i,2) = -maxv*30
       end do

       do i = 1, natm
          nl = numlay(i)
          if(winlay(nl,1) > cps(i,normal_axis_winlay)) winlay(nl,1) = cps(i,normal_axis_winlay)
          if(winlay(nl,2) < cps(i,normal_axis_winlay)) winlay(nl,2) = cps(i,normal_axis_winlay)
       end do

       if(ipridos>=1) then
          write(nfout,'(" !!ldos a range of atomic positions of each layer ")')
          do i = 1, nlayer
             write(nfout,'(2x,i4,2f20.8)') i, winlay(i,1), winlay(i,2)
          end do
       end if

       tmp = winlay(1,1) - crtdst_winlay
       if(kimg == 2) then
          if( tmp < winlay(nlayer,2) - maxhv + crtdst_winlay ) then
!     ( a unit cell has no vacuum region )
             winlay(1,1) = (winlay(1,1) + winlay(nlayer,2) - maxv)*0.5
             winlay(nlayer,2) = winlay(1,1)+maxv
             winlay(mlayer,1) = winlay(1,1)
             winlay(mlayer,2) = winlay(1,1)
          else
!     ( a unit cell has a vacuum region)
             winlay(1,1) = winlay(1,1) - crtdst_winlay
             winlay(nlayer,2) = winlay(nlayer,2) + crtdst_winlay
             winlay(mlayer,1) = winlay(nlayer,2)
             winlay(mlayer,2) = winlay(1,1) + maxv
          endif
       else if(kimg == 1) then
          winlay(1,1) = 0.d0
          winlay(mlayer,2) = maxhv
          winlay(nlayer,2) = winlay(nlayer,2) + crtdst_winlay
          winlay(mlayer,1) = winlay(nlayer,2)
       end if

       do i = 1, nlayer-1
          winlay(i,2) = (winlay(i,2) + winlay(i+1,1))*0.5
          winlay(i+1,1) = winlay(i,2)
       end do
    else
       if(ipridos>=1) write(nfout,'(" !!ldos slicing_selection is illegal")')
    end if

    if(ipridos>=1) then
       write(nfout,'(" !!ldos",3x,a2,5x,a9,11x,a9,8x,a11,5x,a11)') &
            &   "no","min(Bohr)","max(Bohr)","min(Angst.)","max(Angst.)"
!!$       write(nfout,'(" !!ldos    no,        min(Bohr),       max(Bohr),     min(Angst.),    max(Angst.)")')
       do i = 1, nlayer
          write(6,'(" !!ldos ",i4,2f20.8,2f16.4)') i, winlay(i,1), winlay(i,2), winlay(i,1)*BOHR,winlay(i,2)*BOHR
       enddo
       if(mlayer > nlayer) then
          write(6,'(" !!ldos ",i4,2f20.8,2f16.4)') mlayer, winlay(mlayer,1), winlay(mlayer,2) &
               &    ,  winlay(mlayer,1)*BOHR, winlay(mlayer,2)*BOHR
       endif
    end if

  end subroutine set_winlay

  subroutine m_Ldos_alloc_weiwsc_etc()
    integer :: n_ldos_allocated
    integer :: n_lband_allocated

    if ( hardpart_subroutine /= 2 ) then
       if ( sw_rspace_ldos == OFF .or. sw_rspace_lband == OFF ) then
          call m_CD_cp_chgq_to_chgqo() !ASMS
       endif
    endif

    n_total_ldoscal = 0

    if ( sw_aldos == ON .or. sw_calc_wf_atom_decomposition == ON ) then
       if ( noncol ) then
          call alloc_weiwsc_noncl(n_ldos_allocated)
       else
          call alloc_weiwsc(n_ldos_allocated)
       endif
       n_total_ldoscal = n_total_ldoscal + n_ldos_allocated
    end if
    if ( sw_layerdos == ON .or. sw_calc_wf_layer_decomposition == ON ) then
       if ( noncol ) then
          call alloc_weilay_noncl(n_ldos_allocated)
       else
          call alloc_weilay(n_ldos_allocated)
       endif
       n_total_ldoscal = n_total_ldoscal + n_ldos_allocated
    end if
  end subroutine m_Ldos_alloc_weiwsc_etc

  subroutine m_Ldos_dealloc_weiwsc_etc()
    if ( sw_aldos == ON .or. sw_calc_wf_atom_decomposition == ON ) then
       if ( noncol ) then
          call dealloc_weiwsc_noncl()
       else
          call dealloc_weiwsc()
       endif
    endif
    if ( sw_layerdos == ON .or. sw_calc_wf_layer_decomposition == ON ) then
       if ( noncol ) then
          call dealloc_weilay_noncl()
       else
          call dealloc_weilay()
       endif
    endif
  end subroutine m_Ldos_dealloc_weiwsc_etc

  subroutine alloc_weiwsc(n_ldos_allocated)
    integer, intent(out) :: n_ldos_allocated
    n_ldos_allocated = 0
    if(ekmode == ON) then
       if ( sw_lband == ON ) then
          if (hardpart_subroutine==2 .or. sw_rspace_lband==ON )then
             allocate(weiwsc(maldos,np_e,ista_k:iend_k)); weiwsc = 0.d0
          else
             allocate(weiwsc(maldos,neg,kv3)); weiwsc = 0.d0
          endif
       else
          allocate(weiwsc(maldos,neg,kv3_ek)); weiwsc = 0.d0
       endif
    else
       if (hardpart_subroutine==2 .or. sw_rspace_ldos==ON &
            &                     .or. sw_rspace_lband==ON )then
          allocate(weiwsc(maldos,np_e,ista_k:iend_k)); weiwsc = 0.d0
       else
          allocate(weiwsc(maldos,neg,kv3)); weiwsc = 0.d0
       endif
    end if
    n_ldos_allocated = maldos
  end subroutine alloc_weiwsc

! ================================= added by K. Tagami ============= 11.0
  subroutine alloc_weiwsc_noncl(n_ldos_allocated)
    integer, intent(out) :: n_ldos_allocated
    n_ldos_allocated = 0
    if(ekmode == ON) then
       if ( sw_lband == ON ) then
          if(hardpart_subroutine==2 .or. sw_rspace_lband==ON )then
             allocate(weiwsc_noncl(maldos,np_e,ista_k:iend_k,ndim_magmom))
          else
             allocate(weiwsc_noncl(maldos,neg,kv3,ndim_magmom))
          endif
       else
          allocate(weiwsc_noncl(maldos,neg,kv3_ek,ndim_magmom))
       endif
    else if(ekmode == OFF) then
       if(hardpart_subroutine==2 .or. sw_rspace_ldos==ON &
            &                    .or. sw_rspace_lband==ON )then
          allocate(weiwsc_noncl(maldos,np_e,ista_k:iend_k,ndim_magmom))
       else
          allocate(weiwsc_noncl(maldos,neg,kv3,ndim_magmom))
       endif
    end if

    weiwsc_noncl = 0.0d0
    n_ldos_allocated = maldos

  end subroutine alloc_weiwsc_noncl
! =============================================================== 11.0

  subroutine dealloc_weiwsc()
    deallocate(weiwsc)
  end subroutine dealloc_weiwsc

! ================================== added by K. Tagami ========== 11.0
  subroutine dealloc_weiwsc_noncl()
    deallocate(weiwsc_noncl)
  end subroutine dealloc_weiwsc_noncl
! =============================================================== 11.0

  subroutine alloc_weilay(n_ldos_allocated)
    integer, intent(out) :: n_ldos_allocated
!!$    allocate(weilay(mlayer,np_e)); weilay = 0.d0
    n_ldos_allocated = 0
    if(ekmode == ON) then
       if ( sw_lband == ON ) then
          if(hardpart_subroutine==2.or.sw_rspace_lband==ON)then
             allocate(weilay(mlayer,np_e,ista_k:iend_k)); weilay = 0.d0
          else
             allocate(weilay(mlayer,neg,kv3)); weilay = 0.d0
          endif
       else
          allocate(weilay(mlayer,neg,kv3_ek)); weilay = 0.d0
       endif
    else if(ekmode == OFF) then
       if(hardpart_subroutine==2.or.sw_rspace_ldos==ON)then
          allocate(weilay(mlayer,np_e,ista_k:iend_k)); weilay = 0.d0
       else
          allocate(weilay(mlayer,neg,kv3)); weilay = 0.d0
       endif
    end if
    n_ldos_allocated = mlayer
  end subroutine alloc_weilay

! ============================ added by K. Tagami ================= 11.0
  subroutine alloc_weilay_noncl(n_ldos_allocated)
    integer, intent(out) :: n_ldos_allocated

    n_ldos_allocated = 0
    if(ekmode == ON) then
       if ( sw_lband == ON ) then
          if(hardpart_subroutine==2.or.sw_rspace_lband==ON)then
             allocate(weilay_noncl(mlayer,np_e,ista_k:iend_k,ndim_magmom))
          else
             allocate(weilay_noncl(mlayer,neg,kv3,ndim_magmom))
          endif
       else
          allocate(weilay_noncl(mlayer,neg,kv3_ek,ndim_magmom))
       endif
    else if(ekmode == OFF) then
       if(hardpart_subroutine==2.or.sw_rspace_ldos==ON)then
          allocate(weilay_noncl(mlayer,np_e,ista_k:iend_k,ndim_magmom))
       else
          allocate(weilay_noncl(mlayer,neg,kv3,ndim_magmom))
       endif
    end if
    weilay_noncl = 0.d0
    n_ldos_allocated = mlayer
  end subroutine alloc_weilay_noncl
! =============================================================== 11.0

  subroutine dealloc_weilay()
    deallocate(weilay)
  end subroutine dealloc_weilay

! ============================== added by K. Tagami ========== 11.0
  subroutine dealloc_weilay_noncl()
    deallocate(weilay_noncl)
  end subroutine dealloc_weilay_noncl
! ============================================================= 11.0

  subroutine m_Ldos_wd_natm2_and_totch( nfldos )
    integer, intent(in) :: nfldos
    if(mype == 0 .and. ipridos>=0 ) then
!       write(nfldos,'(" natm2 = ",i8)') naldos ! natm2
!       write(nfldos,'(" totch = ",f20.8)') totch
       write(nfldos,'(" num_atoms  = ",i8)') naldos ! natm2
       write(nfldos,'(" num_layers = ",i8)') nlayer
       write(nfldos,'(" num_bands  = ",i8)') neg -num_extra_bands
       write(nfldos,*)
    end if
  end subroutine m_Ldos_wd_natm2_and_totch

  subroutine m_Ldos_cal_ek( nfldos )
    use m_Files, only : nfzaj_kall
    use m_IterationNumbers, only : m_Iter_reset_nk_in_the_process
    use m_ES_IO, only : m_ESIO_rd_next_wfs_ek
    use m_Files, only : m_Files_open_nfzaj_kall
    use m_ES_ortho, only : m_ES_modified_gram_schmidt
    use m_ES_occup, only : m_ESoc_check_if_metalic
    use m_ES_nonlocal, only : m_ES_betar_dot_WFs_4_each_k

    integer, intent(in) :: nfldos

    integer :: id_sname = -1
    integer :: i, is, ikt, nk, ierr, ik, ib
    integer :: sw_atom_decomp, sw_layer_decomp
    real(kind=DP), allocatable, dimension(:,:,:) :: wei
    integer, allocatable, dimension(:) :: kmap

    call m_Files_open_nfzaj_kall()
    rewind(nfzaj_kall)
    allocate(kmap(kv3))

    call tstatc0_begin('m_Ldos_cal_ek ',id_sname,1)
    call m_Iter_reset_nk_in_the_process()

    sw_atom_decomp = OFF;  sw_layer_decomp = OFF
    if ( sw_aldos == ON .or. sw_calc_wf_atom_decomposition == ON ) then
       sw_atom_decomp = ON
    endif
    if ( sw_layerdos == ON .or. sw_calc_wf_layer_decomposition == ON ) then
       sw_layer_decomp = ON
    endif

    do nk=1, kv3_ek, kv3
       if ( sw_atom_decomp==ON ) weiwsc = 0.d0
       if ( sw_layer_decomp==ON ) weilay = 0.d0
       call KpointNumber_Setting2()
       call Preparation_ek()
       call Preparation_for_mpi_ek()
       call PseudoPotential_ek()

       kmap = 0
       if(sw_modified_kpoint_increment == ON)then
          do i=0,nrank_k-1
             do is=1,nspin
                ikt = nspin*(nis_kv3_ek(i)-1)+(nk-1)*nspin+is
                kmap(i+1) = ikt
                if(ikt>kv3_ek) exit
                call m_ESIO_rd_next_wfs_ek(i+1, nfout, nfzaj_kall)
                if(map_k(i+1) /= myrank_k) cycle
                call m_ES_betar_dot_WFs_4_each_k(nfout, i+1)
             enddo
          enddo
       else
          do is=1, nspin, af+1
             do ik=is, kv3+is-nspin, nspin
                ikt = nk + ik - 1
                if(ikt>kv3_ek) exit
                kmap(ik) = ikt
                call m_ESIO_rd_next_wfs_ek(ik, nfout, nfzaj_kall)
                if(map_k(ik) /= myrank_k) cycle
                call m_ES_betar_dot_WFs_4_each_k(nfout, ik)
             enddo
          enddo
       endif
       if ( sw_ldos == ON .or. sw_lband ==ON ) call m_Ldos_cal( nfldos, .false. )

       if (sw_atom_decomp==ON) then
          call wd_weight_ek(maldos,kmap,weiwsc,ALDOS,nfldos)
       endif
       if (sw_layer_decomp==ON) then
          call wd_weight_ek(mlayer,kmap,weilay,LAYERDOS,nfldos)
       endif

    enddo
    deallocate(kmap)

    if(sw_atom_decomp==ON .and. sw_ac_mesh==ON) then
      call check_sum_full('total ',.true.)
      if(sw_checksum == ON) call check_sum_full('total ')
    endif

    call tstatc0_end(id_sname)
  end subroutine m_Ldos_cal_ek


  subroutine wd_weight_ek(m,kmap,wei,mode,nfldos)
    real(kind=DP), intent(in), dimension(m,neg,kv3_ek) :: wei
    integer, intent(in), dimension(kv3) :: kmap
    integer, intent(in) :: mode, nfldos
    integer :: m, ib, ib_ordr, kend, kstep, iws, ik, ii

    if ( nrank_e > 1 ) then
       Do ii=1, kv3_ek
          call mpi_allreduce( mpi_in_place, wei(:,:,ii), m*neg, mpi_double_precision, &
               &              mpi_sum, MPI_CommGroup, ierr )
       End Do
    endif

    if(mype == 0) then
       do ik=1,kv3
         if (kmap(ik)>kv3_ek .or. kmap(ik)==0 ) cycle
         write(nfout,'(a,3i8)') '!** in wd_weight ',ik,kmap(ik),mode
!         write(nfldos,'(2i8)') mode,kmap(ik)
         if ( mode == 1 ) then
            write(nfldos,'(2i8,10x,a)') mode, kmap(ik), " # Atom decomposed,  ik"
         else
            write(nfldos,'(2i8,10x,a)') mode, kmap(ik), " # Layer decomposed, ik"
         endif

         do ib = 1, neg -num_extra_bands
            if(ilen >= m) then
               write(nfldos,'(i4,")",5f15.10)') ib,(wei(iws,ib,ik),iws=1,m)
            else
               write(nfldos,'(i4,")",5f15.10)') ib,(wei(iws,ib,ik),iws=1,ilen)
               write(nfldos,'(5x,5f15.10)') (wei(iws,ib,ik),iws=ilen+1,m)
            end if
         end do
      enddo
    endif
  end subroutine wd_weight_ek

  subroutine check_sum(ik,aword)
    integer,intent(in) :: ik
    character*(*),intent(in) :: aword

    integer :: sw_atom_decomp, sw_layer_decomp
    integer :: ib, ia, ilay, ib_ordr, ik_t
    real(kind=DP), allocatable, dimension(:) :: sum_aldos, sum_layer

    sw_atom_decomp = OFF;  sw_layer_decomp = OFF
    if ( sw_aldos == ON .or. sw_calc_wf_atom_decomposition == ON ) then
       sw_atom_decomp = ON
    endif
    if ( sw_layerdos == ON .or. sw_calc_wf_layer_decomposition == ON ) then
       sw_layer_decomp = ON
    endif

    if(sw_atom_decomp == ON) allocate(sum_aldos(neg))
    if(sw_layer_decomp == ON) allocate(sum_layer(neg))

    if ( sw_lband == ON ) then
       ik_t = ik
    else
       if(ekmode /= ON) then
          ik_t = ik
       else if(ekmode == ON) then
          ik_t = 1
       end if
    endif

    if(sw_atom_decomp == ON) then
       sum_aldos = 0.d0
       do ib = ista_e, iend_e, istep_e
!!$               ibt = neg_g(ib)
          do ia = 1, natm2 + 1
             sum_aldos(ib) = sum_aldos(ib) + weiwsc(ia,map_z(ib),ik_t)
!!$               sum_aldos(neg_g(ib)) = sum_aldos(neg_g(ib)) + weiwsc(ia,ibt,ik_t)
          end do
       end do
       call mpi_allreduce(MPI_IN_PLACE,sum_aldos,neg,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr)
    end if

    if(sw_layer_decomp == ON) then
       sum_layer = 0.d0
       do ib = ista_e, iend_e, istep_e
          do ilay = 1, mlayer
             sum_layer(ib) = sum_layer(ib) + weilay(ilay,map_z(ib),ik_t)
          end do
       end do
       call mpi_allreduce(MPI_IN_PLACE,sum_layer,neg,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr)
    end if

!!$      if(mype == 0) then
    if(ipridos>=1) then
       if(sw_atom_decomp == ON .and. sw_layer_decomp == ON) then
          write(nfout,9001) aword,"atomic, layer",neg,natm2,ik
          do ib = 1, neg
             write(nfout,'(" !!ldos ", i4,2f16.8)') ib, sum_aldos(ib), sum_layer(ib)
          end do
       else if(sw_atom_decomp == ON) then
          write(nfout,9002) aword,"atomic",neg,natm2,ik
          do ib = 1, neg
             write(nfout,'(" !!ldos ", i4,f16.8)') ib, sum_aldos(ib)
          end do
       else if(sw_layer_decomp == ON) then
          write(nfout,9002) aword,"layer ",neg,natm2,ik
          do ib = 1, neg
             write(nfout,'(" !!ldos ", i4,f16.8)') ib, sum_layer(ib)
          end do
       end if
9001     format(' !!ldos ',a13,' -- iban , sum(',a6,') --  neg = ',i6,' natm2 = ',i6,' ik = ',i6)
9002     format(' !!ldos ',a8, ' -- iban , sum(',a6,') --  neg = ',i6,' natm2 = ',i6,' ik = ',i6)
    end if

    if(sw_atom_decomp == ON) deallocate(sum_aldos)
    if(sw_layer_decomp == ON) deallocate(sum_layer)

  end subroutine check_sum

  subroutine check_sum_full(aword,normalize)
    character*(*),intent(in) :: aword
    logical, intent(in), optional :: normalize
    integer :: ik, ib, ia, ilay, ike,j,ierr
    integer, parameter :: num_kset = 4
    integer :: sw_atom_decomp, sw_layer_decomp
    real(kind=DP), allocatable, dimension(:,:) :: sum_aldos, sum_layer, sum_mpi
    logical :: norma
    norma = .false.
    if(present(normalize)) norma = .true.

    sw_atom_decomp = OFF;  sw_layer_decomp = OFF
    if ( sw_aldos == ON .or. sw_calc_wf_atom_decomposition == ON ) then
       sw_atom_decomp = ON
    endif
    if ( sw_layerdos == ON .or. sw_calc_wf_layer_decomposition == ON ) then
       sw_layer_decomp = ON
    endif

    if(sw_atom_decomp == ON) allocate(sum_aldos(neg,kv3))
    if(sw_layer_decomp == ON) allocate(sum_layer(neg,kv3))

    if(sw_atom_decomp==ON) then
       sum_aldos = 0.d0
       do ik = ista_k, iend_k
          do ib = ista_e, iend_e, istep_e
             do ia = 1, natm2 + 1
                sum_aldos(ib,ik) = sum_aldos(ib,ik) + weiwsc(ia,map_z(ib),ik)
             end do
          end do
       end do
       call mpi_allreduce(MPI_IN_PLACE,sum_aldos,neg*kv3,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
    end if
    if(sw_layer_decomp == ON) then
       sum_layer = 0.d0
       do ik = ista_k, iend_k
          do ib = ista_e, iend_e, istep_e
             do ilay = 1, mlayer
                sum_layer(ib,ik) = sum_layer(ib,ik) + weilay(ilay,map_z(ib),ik)
             end do
          end do
       end do
       call mpi_allreduce(MPI_IN_PLACE,sum_layer,neg*kv3,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
    end if

    if(ipridos>=1) then
       do ik = 1, kv3, num_kset
          ike = min(kv3,ik+num_kset-1)
          if(sw_atom_decomp == ON .and. sw_layer_decomp == ON) then
             write(nfout,9001) aword,"atomic,layer",neg,natm2,ik,ike
             write(nfout,'(" !!ldos iban ",4a24)')   (("   aldos      layerdos  "),j=1,4)
             do ib = 1, neg
                write(nfout,'(" !!ldos ",i4,8f12.8)') ib,(sum_aldos(ib,j),sum_layer(ib,j),j=ik,ike)
             end do
          else if(sw_atom_decomp == ON) then
             write(nfout,9002) aword,"atomic",neg,natm2,ik,ike
             do ib = 1, neg
                write(nfout,'(" !!ldos ",i4,4f14.8)') ib,(sum_aldos(ib,j),j=ik,ike)
             end do
          else if(sw_layer_decomp == ON) then
             write(nfout,9002) aword,"layer ",neg,natm2,ik,ike
             do ib = 1, neg
                write(nfout,'(" !!ldos ",i4,4f14.8)') ib,(sum_layer(ib,j),j=ik,ike)
             end do
          end if
       end do
9001     format(' !!ldos ',a13,' -- iban , sum(',a6,') --  neg = ',i6,' natm2 = ',i6,' ik = [',i4,':',i4,']')
9002     format(' !!ldos ',a8, ' -- iban , sum(',a6,') --  neg = ',i6,' natm2 = ',i6,' ik = [',i4,':',i4,']')
    end if
    if(norma .and. sw_atom_decomp==ON .and. sw_ac_mesh==ON)then
       do ik = ista_k, iend_k
          do ib = ista_e, iend_e, istep_e
             weiwsc(:,ib,ik) = weiwsc(:,ib,ik)/sum_aldos(ib,ik)
          enddo
       enddo
    endif
  end subroutine check_sum_full



  subroutine cal_weiwsc(ik,is)
    integer, intent(in) :: ik, is
  end subroutine cal_weiwsc

  subroutine m_Ldos_cal( nfldos, print_flag )
    use m_IterationNumbers,   only : nk_in_the_process
    use m_Ionic_System,       only : speciesname

    integer, intent(in) :: nfldos
    logical, intent(in) :: print_flag
    integer,allocatable,dimension(:,:) :: meshwk, meshwk_layer
    integer :: ik,is,ip, ib,j,i,ia
    integer :: ekmode_t
    integer :: sw_full_fftmesh
    integer :: sw_atom_decomp, sw_layer_decomp

    real(kind=DP) :: denom, chgq0
    real(kind=DP), allocatable, dimension(:) :: zrhoik ! d(kg*kimg)
    real(kind=DP), allocatable, dimension(:) :: bfft   ! d(nfft)

    real(kind=DP), allocatable, dimension(:,:,:) :: wsctmp,laytmp
    real(kind=DP), allocatable, dimension(:,:,:) :: wscsoft,wschard
    real(kind=DP), allocatable, dimension(:) :: zfcos,zfsin
    real(kind=DP), allocatable, dimension(:,:) :: psi_tmp

    sw_atom_decomp = OFF;  sw_layer_decomp = OFF
    if ( sw_aldos == ON .or. sw_calc_wf_atom_decomposition == ON ) then
       sw_atom_decomp = ON
    endif
    if ( sw_layerdos == ON .or. sw_calc_wf_layer_decomposition == ON ) then
       sw_layer_decomp = ON
    endif

    sw_full_fftmesh=OFF
    if(integration_dimension_winlay==3) sw_full_fftmesh=ON
    if(sw_layer_decomp == ON) allocate(zrhoik(kg*kimg))

    denom = 1.d0/product(fft_box_size_WF(1:3,1))

    call m_CD_keep_retrieve_hsr(.true.)
    if(sw_atom_decomp==ON .and. sw_ac_mesh==ON)then
       allocate(wscsoft(maldos,neg,kv3));wscsoft=0.d0
       allocate(wschard(maldos,neg,kv3));wschard=0.d0
    endif

     if(modnrm == EXECUT) call m_CD_set_ylm_enl_etc()

     if ( mype == 0 ) call wd_header

!==================================== added by K. Tagami ============= 11.0
     if ( sw_atom_decomp == ON ) weiwsc = 0.d0
     if ( sw_layer_decomp == ON ) weilay = 0.d0
!===================================================================== 11.0

     do is = 1, nspin, af+1
        if(hardpart_subroutine==2.or.sw_rspace_ldos==ON &
             & .or.  sw_rspace_lband ==ON ) then

           do ik = is, kv3+is-nspin, nspin
              if(map_k(ik) /= myrank_k) cycle
              if ( ekmode==ON ) then
                 if ( nk_in_the_process +ik -1 > kv3_ek ) cycle
              endif

              if(ipridos >= 2) write(nfout,'(" !!ldos ik = ",i6)') ik
              call m_FFT_alloc_WF_work()
              if(sw_atom_decomp==ON .and. sw_ac_mesh==OFF) then
                 call fillup_meshwk(mode=SOFTPART)
              endif
              if(sw_layer_decomp==ON .and. sw_full_fftmesh==ON) then
                 call fillup_meshwk_layer(mode=SOFTPART)
              endif

              allocate(bfft(nfft))
              if(ipridos >= 2 .and. printable) write(nfout,'(" !!ldos nfft = ",i6," <<m_Ldos_cal>>")') nfft
              if(ipridos >= 2 .and. printable) &
                   & call m_ES_wd_zaj_small_portion(nfout,ik," -- before m_ES_WF_in_Rspace <<m_Ldos_cal>> --",46)
              do ib = ista_e, iend_e, istep_e
                 bfft = 0.d0
                 if ( sw_band_unfolding == ON .and. band_unfolding_active ) then
                    call m_ES_WF_in_Rspace_kt( ik, ik, ik, &
                         &                     zaj_l(:,map_z(ib),ik,:), bfft )
                 else
                    call m_ES_WF_in_Rspace(ik,ib,bfft)
                 endif

                 if(sw_atom_decomp == ON.and.sw_ac_mesh==OFF) call substitute_weiwsc(ik,ib,.false.)
                 if(sw_atom_decomp == ON.and.sw_ac_mesh==ON)  call substitute_weiwsc_acmesh(ik,ib,.false.,.false.,wscsoft)
                 if(sw_layer_decomp==ON) then
                    if(sw_full_fftmesh==ON) call substitute_weilay_fftmesh(ik,ib,.false.)
                    if(sw_full_fftmesh==OFF) call substitute_weilay(ik,ib,.false.)
                 end if
              end do
              call check_sum(ik,"softpart",.false.,modnrm /= EXECUT)
              if(sw_ac_mesh==ON .and. modnrm /= EXECUT) &
              &  call check_sum(ik,"softpart",.false.,.false.)
              deallocate(bfft)
              call m_FFT_dealloc_WF_work()

              if(modnrm == EXECUT) then
                 call m_FFT_alloc_CD_box()
                 allocate(bfft(nfftp_nonpara))

                 if(sw_atom_decomp==ON .and. sw_ac_mesh==OFF) call fillup_meshwk(mode=HARDPART)
                 if(sw_layer_decomp==ON .and. sw_full_fftmesh==ON) call fillup_meshwk_layer(mode=HARDPART)

                 do ib = ista_e, iend_e, istep_e
                    if(ipridos>=2) write(nfout,'(" !ldos   ib = ", i5)') ib
                    if(sw_rspace_ldos==ON .or. sw_rspace_lband==ON)then
                       if(sw_layer_decomp == ON .and. sw_full_fftmesh==OFF) then
                          call m_CD_hardpart_sub2_rs(nfout,is,ik,ib,bfft,chgq0)
                       else
                          call m_CD_hardpart_sub2_rs(nfout,is,ik,ib,bfft)
                       endif
                    else
                       call m_CD_hardpart_sub2(nfout,is,ik,ib,chgq0)  ! -> chgq_enl
                       call m_CD_map_chgqenl_to_fft_box(is,nfftp_nonpara,bfft) ! -> bfft
                       call m_FFT_CD_inverse0(nfout,bfft) ! bfft(G_sp.) -> bfft(R_sp.)
                    endif

                    if(sw_atom_decomp == ON .and. sw_ac_mesh==OFF) &
                   &   call substitute_weiwsc_cd(ik,ib,.false.)

                    if(sw_atom_decomp == ON.and.sw_ac_mesh==ON)  &
                   &   call substitute_weiwsc_acmesh(ik,ib,.false.,.true.,wschard)

                    if(sw_layer_decomp==ON .and.sw_full_fftmesh==ON) call substitute_weilay_cd_fftmesh(1,ib,.false.)
                    if(sw_layer_decomp==ON .and. sw_full_fftmesh==OFF) call substitute_weilay_cd(ik,ib,chgq0,.false.)
                 end do
                 if(sw_atom_decomp==ON .and. sw_ac_mesh==OFF) deallocate(meshwk)
                 if(sw_layer_decomp==ON .and. sw_full_fftmesh==ON) deallocate(meshwk_layer)
                 deallocate(bfft)
                 call m_FFT_dealloc_CD_box()
                 call check_sum(ik,"total   ",.false.,.true.)
                 if(sw_ac_mesh==ON) call check_sum(ik,"total   ",.false.,.false.)
              end if

           end do

        else
!========================= modiifed by K. Tagami ============= 11.0
            weiwsc = 0.d0
            weilay = 0.d0
!============================================================= 11.0

           do ik = is, kv3+is-nspin, nspin
              if(map_k(ik) /= myrank_k) cycle
              if ( ekmode==ON ) then
                 if ( nk_in_the_process +ik -1 > kv3_ek ) cycle
              endif

              if(ipridos >= 2) write(nfout,'(" !!ldos ik = ",i6)') ik
              call m_FFT_alloc_WF_work()
              if(sw_atom_decomp==ON .and. sw_ac_mesh==OFF) call fillup_meshwk(mode=SOFTPART)
              if(sw_layer_decomp==ON .and. sw_full_fftmesh==ON) call fillup_meshwk_layer(mode=SOFTPART)

              allocate(bfft(nfft));bfft=0.d0
              if(ipridos >= 2 .and. printable) write(nfout,'(" !!ldos nfft = ",i6," <<m_Ldos_cal>>")') nfft
              if(ipridos >= 2 .and. printable) &
                   & call m_ES_wd_zaj_small_portion(nfout,ik," -- before m_ES_WF_in_Rspace <<m_Ldos_cal>> --",46)
              do ib = 1,neg
                 if(map_e(ib)/=myrank_e) cycle
                 bfft = 0.d0
                 if ( sw_band_unfolding == ON .and. band_unfolding_active ) then
                    call m_ES_WF_in_Rspace_kt( ik, ik, ik, &
                         &                     zaj_l(:,map_z(ib),ik,:), bfft )
                 else
                    call m_ES_WF_in_Rspace(ik,ib,bfft)
                 endif

                 if(sw_atom_decomp==ON.and.sw_ac_mesh==OFF) call substitute_weiwsc(ik,ib,.true.)
                 if(sw_atom_decomp==ON.and.sw_ac_mesh==ON)  call substitute_weiwsc_acmesh(ik,ib,.true.,.false.,wscsoft)
                 if(sw_layer_decomp==ON .and. sw_full_fftmesh==ON) call substitute_weilay_fftmesh(ik,ib,.true.)
                 if(sw_layer_decomp==ON .and. sw_full_fftmesh==OFF) call substitute_weilay(ik,ib,.true.)
              end do
              deallocate(bfft)
              call m_FFT_dealloc_WF_work()
           enddo
           if(npes>1)then
              if(sw_atom_decomp==ON.and.sw_ac_mesh==OFF) then
                 allocate(wsctmp(maldos,neg,kv3));wsctmp=0.d0
                 call mpi_allreduce( weiwsc, wsctmp, kv3*maldos*neg, &
                      &              mpi_double_precision, &
                      &              mpi_sum,MPI_CommGroup,ierr )
!============================== modified by K. Tagami =============== 11.0
                  weiwsc = wsctmp

                 Do ik=is, kv3+is-nspin,nspin
                    weiwsc(:,:,ik) = wsctmp(:,:,ik)
                 End do
!==================================================================== 11.0
                 deallocate(wsctmp)
              endif
              if(sw_layer_decomp==ON)then
                 allocate(laytmp(mlayer,neg,kv3));laytmp=0.d0
                 call mpi_allreduce( weilay, laytmp, kv3*mlayer*neg, &
                      &              mpi_double_precision, &
                      &              mpi_sum,MPI_CommGroup,ierr )
!============================== modified by K. Tagami =============== 11.0
                  weilay = laytmp

                 Do ik=is, kv3+is-nspin,nspin
                    weilay(:,:,ik) = laytmp(:,:,ik)
                 End do
!==================================================================== 11.0
                 deallocate(laytmp)
              endif
           endif
           if(sw_ac_mesh==OFF)then
              do ik=is,kv3+is-nspin,nspin
                 call check_sum(ik,"softpart",.true.,modnrm /= EXECUT)
              enddo
           endif

           if(modnrm == EXECUT) then
              allocate(bfft(nfftp_nonpara))
              do ik = is, kv3+is-nspin, nspin
                 call m_FFT_alloc_CD_box()
                 if(sw_atom_decomp==ON .and. sw_ac_mesh==OFF) call fillup_meshwk(mode=HARDPART)
                 if(sw_layer_decomp==ON .and. sw_full_fftmesh==ON) call fillup_meshwk_layer(mode=HARDPART)

                 do ib = 1,neg
                    if(ipridos>=2) write(nfout,'(" !ldos   ib = ", i5)') ib
                    call m_CD_hardpart_sub(nfout,is,ik,ib,chgq0)  ! -> chgq_l
                    call m_CD_map_valence_charge_to_fft_box(is,nfftp_nonpara,bfft)
                    call m_FFT_CD_inverse0(nfout,bfft) ! bfft(G_sp.) -> bfft(R_sp.)
                    if(sw_atom_decomp==ON .and.sw_ac_mesh==OFF) call substitute_weiwsc_cd(ik,ib,.true.)
                    if(sw_atom_decomp==ON.and.sw_ac_mesh==ON)   call substitute_weiwsc_acmesh(ik,ib,.true.,.true.,wschard)
                    if(sw_layer_decomp==ON .and. sw_full_fftmesh==ON) call substitute_weilay_cd_fftmesh(ik,ib,.true.)
                    if(sw_layer_decomp==ON .and. sw_full_fftmesh==OFF) call substitute_weilay_cd(ik,ib,chgq0,.true.)
                 end do
                 if(sw_atom_decomp==ON .and. sw_ac_mesh==OFF) deallocate(meshwk)
                 if(sw_layer_decomp==ON .and. sw_full_fftmesh==ON) deallocate(meshwk_layer)
                 call m_FFT_dealloc_CD_box()
                 if(sw_ac_mesh==OFF) call check_sum(ik,"total   ",.true.,.true.)
              enddo
              deallocate(bfft)
           endif

        endif
     enddo

     if(sw_atom_decomp==ON .and. sw_ac_mesh==ON) then
        call mpi_allreduce(mpi_in_place,wscsoft,kv3*neg*maldos,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
        call mpi_allreduce(mpi_in_place,wschard,kv3*neg*maldos,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
        weiwsc = wscsoft+wschard
        do is=1,nspin,af+1
           do ik = is, kv3+is-nspin, nspin
              call check_sum(ik,"total   ",.true.,.true.)
              call check_sum(ik,"total   ",.true.,.false.)
           enddo
        enddo
     endif

     if(npes > 1) call mpi_barrier(MPI_CommGroup,ierr)
!     if(sw_save_ldos_weight == ON .and. ekmode == OFF) then
     if(sw_save_ldos_weight == ON .and. print_flag ) then
        if ( sw_lband == ON ) then
           if(sw_atom_decomp==ON) call wd_weight2(ALDOS,nfldos)
           if(sw_layer_decomp==ON) call wd_weight2(LAYERDOS,nfldos)
        else
           if(hardpart_subroutine==2 .or. sw_rspace_ldos==ON ) then
              if(sw_atom_decomp==ON) call wd_weight(ALDOS,nfldos)
              if(sw_layer_decomp==ON) call wd_weight(LAYERDOS,nfldos)
           else
              if(sw_atom_decomp==ON) call wd_weight_serial(ALDOS,nfldos)
              if(sw_layer_decomp==ON) call wd_weight_serial(LAYERDOS,nfldos)
           endif
        endif
     end if

    call m_CD_keep_retrieve_hsr(.false.)

    if(sw_atom_decomp==ON .and. sw_ac_mesh==ON)then
       deallocate(wscsoft)
       deallocate(wschard)
    endif
    if(modnrm == EXECUT) call m_CD_dealloc_ylm_enl_etc()
    if(hardpart_subroutine/=2 .and. &
         &  (sw_rspace_ldos==OFF .or. sw_rspace_lband==OFF) ) call m_CD_restore_chgq()
    if(allocated(meshwk)) deallocate(meshwk)

    if(sw_layer_decomp==ON) deallocate(zrhoik)

  contains

    subroutine wd_header
      integer :: i, j

      if ( sw_lband == OFF ) then
         if(sw_save_ldos_weight == ON .and. ekmode==OFF) then
            call m_Ldos_wd_natm2_and_totch( nfldos )
         endif
         if ( ekmode == ON .and. nk_in_the_process == 1 ) then
            call m_Ldos_wd_natm2_and_totch( nfldos )
         endif
      endif
      if ( sw_lband == ON ) then
         if ( ( ekmode == ON .and. nk_in_the_process == 1 ) &
              &  .or. ( ekmode == OFF) ) then
            write(nfldos,'(A)') "# Local Decomposition for bands"
            write(nfldos,*) ""
            write(nfldos,'(A,I8)') ' num_kpoints = ', max(kv3,kv3_ek) /ndim_spinor
            write(nfldos,'(A,I8)') ' num_bands   = ', neg -num_extra_bands
            write(nfldos,'(A,I8)') ' nspin       = ', nspin /ndim_spinor
            write(nfldos,*) ""
            write(nfldos,'(A)') "# Decomposition Info."

            if ( sw_atom_decomp == ON ) then
               write(nfldos,'(A,I8,A)') ' num_atoms   = ', natom_decomp +1, &
                    &        "  ( last 1 corresponds to external region )"
            endif
            if ( sw_layer_decomp == ON ) then
               write(nfldos,'(A,I8,A)') ' num_layers  = ', nlayer +1, &
                    &        "  ( last 1 corresponds to external region )"
            endif
            write(nfldos,*) ""

            if ( sw_atom_decomp == ON ) then
               write(nfldos,'(A)') "# Atom Info."
               write(nfldos,'(A)') "   no.    target_atom    species"
               Do i=1, natom_decomp
                  j = atom_decomp_map(i)
                  write(nfldos,'(I5,5X,I5,12X,A)') i, j, speciesname(ityp(j))
               end Do
               write(nfldos,'(I5)') natom_decomp +1
               write(nfldos,*)
            endif
            if ( sw_layer_decomp == ON ) then
               write(nfldos,'(A)') "# Layer Info."
               write(nfldos,'(A)') "   no.       min.         max.  (Angstrom)"
               Do i=1, nlayer +1
                  write(nfldos,'(I5,2(5X,F8.4))') i, winlay(i,1)*BOHR, winlay(i,2)*BOHR
               end Do
               write(nfldos,*)
            endif
         endif
      endif
    end subroutine wd_header

    subroutine wd_weight2(mode,nfldos)
      integer, intent(in) :: mode, nfldos
      integer :: m, ib, ib_ordr, kend, kstep, iws, ik, iks
      integer :: itmp, count
      real(kind=DP),allocatable,dimension(:,:,:):: wei

!      write(nfout,'(a,2i8)') '!** in wd2_weight ',ik,mode

      if(mode==ALDOS) then
         m = maldos
      else if(mode==LAYERDOS) then
         m = mlayer
      end if
      allocate(wei(m,neg,kv3)); wei = 0.d0

      if(mode==ALDOS)then
         do ik=ista_k, iend_k
            do ib = ista_e, iend_e, istep_e
               wei(1:m,ib,ik) = weiwsc(1:m,map_z(ib),ik)
            end do
         end do
      else if(mode==LAYERDOS)then
         do ik=ista_k, iend_k
            do ib = ista_e, iend_e, istep_e
               wei(1:m,ib,ik) = weilay(1:m,map_z(ib),ik)
            end do
         end do
      endif

      if (npes>1) then
         call mpi_allreduce( MPI_IN_PLACE, wei, m*neg*kv3, mpi_double_precision, &
              &              mpi_sum, MPI_CommGroup, ierr )
      end if

      if(mype == 0) then
         Do ik=1, kv3
            iks = 0
            if ( sw_lband == ON .and. ekmode == ON ) then
               iks = nk_in_the_process -1
               if ( ik +iks > kv3_ek ) cycle
            endif
            if(mode==ALDOS) then
               write(nfldos,'(" weight for each atomic cell    nk = ",i8)') ik+iks
            else if(mode==LAYERDOS) then
               write(nfldos,'(" weight for each layer          nk = ",i8)') ik+iks
            end if

            if ( mode==ALDOS .and. sw_wd_only_specified_atoms == ON ) then
               do ib = 1, neg -num_extra_bands
                  write(nfldos,'(i4,A)',advance="no") ib, ")"
                  count = 0
                  do iws=1, m-1
                     if ( if_aldos_full(iws)==1 ) then
                        count = count +1
                        if ( count > 1 .and. mod(count,ilen)==1 ) then
                           write(nfldos,'(5x)',advance="no")
                        endif
                        write(nfldos,'(f15.10)',advance="no") wei(iws,ib,ik)
                        if ( mod(count,ilen)==0 ) write(nfldos,*)
                     endif
                  end do
                  count = count +1
                  if ( count > 1 .and. mod(count,ilen)==1 ) then
                     write(nfldos,'(5x)',advance="no")
                  endif
                  write(nfldos,'(f15.10)') wei(maldos,ib,ik)
               end do
            else
               do ib = 1, neg -num_extra_bands
                  if(ilen >= m) then
                     write(nfldos,'(i4,")",5f15.10)') ib,(wei(iws,ib,ik),iws=1,m)
                  else
                     write(nfldos,'(i4,")",5f15.10)') ib,(wei(iws,ib,ik),iws=1,ilen)
                     write(nfldos,'(5x,5f15.10)') (wei(iws,ib,ik),iws=ilen+1,m)
                  end if
               end do
            endif
         end Do
      end if
      deallocate(wei)
    end subroutine wd_weight2

    subroutine wd_weight(mode,nfldos)
      integer, intent(in) :: mode, nfldos
      integer :: m, ik, ib, ib_ordr, kend, kstep, iws, ip, iks
      real(kind=DP),allocatable,dimension(:):: wei

      if(mode==ALDOS) then
         m = maldos
      else if(mode==LAYERDOS) then
         m = mlayer
      end if

      allocate(wei(m))

      if(ekmode==OFF) then
         kend = kv3; kstep = af+1
      else if(ekmode==ON) then
!         kend = 1;   kstep = 1
         kend = kv3;   kstep = 1
      end if

      if(ipridos >= 3) then
         do ik = 1, kv3, af+1
            if(map_k(ik) /= myrank_k) cycle
            write(nfout,'(" !ldos: ik = ",i8)') ik
            do ib = 1, neg
               write(nfout,'(" !ldos: ib, ib_ordr = ",2i20)') ib, neordr(ib,ik)
            end do
         end do
      end if

      do ik = 1, kend, kstep
         iks = 0
         if ( sw_lband == ON .and. ekmode == ON ) then
            iks = nk_in_the_process -1
            if ( ik +iks > kv3_ek ) cycle
         endif
         if(mype == 0) then
            if(mode==ALDOS) then
               write(nfldos,'(" weight for each atomic cell    nk = ",i8)') ik+iks
            else if(mode==LAYERDOS) then
               write(nfldos,'(" weight for each layer          nk = ",i8)') ik+iks
            end if
         end if

         do ib = 1, neg -num_extra_bands

            if(nrank_k >= 2) then
               if(map_ek(ib,ik) == mype) ib_ordr = neordr(ib,ik)
               call mpi_bcast(ib_ordr,1,mpi_integer,map_ek(ib,ik),MPI_CommGroup,ierr)
            else
               ib_ordr = neordr(ib,ik)
            end if

            if(map_ek(ib_ordr,ik) == mype) then
               if(mode == ALDOS) then
                  wei(1:m) = weiwsc(1:m,map_z(ib_ordr),ik)
               else if(mode == LAYERDOS) then
                  wei(1:m) = weilay(1:m,map_z(ib_ordr),ik)
               end if
               if(map_ek(ib_ordr,ik) /= 0) &
                  call mpi_send(wei,m,mpi_double_precision,0,1,MPI_CommGroup,ierr)
             else if(mype == 0 .and. map_ek(ib_ordr,ik) /= 0) then
                call mpi_recv(wei,m,mpi_double_precision,map_ek(ib_ordr,ik),1,MPI_CommGroup,istatus,ierr)
             end if
             if(mype == 0)  then
               if(ilen >= m) then
                  write(nfldos,'(i4,")",5f15.10)') ib,(wei(iws),iws=1,m)
               else
                  write(nfldos,'(i4,")",5f15.10)') ib,(wei(iws),iws=1,ilen)
                  write(nfldos,'(5x,5f15.10)') (wei(iws),iws=ilen+1,m)
               end if
             end if
          end do
       end do

       deallocate(wei)
     end subroutine wd_weight

    subroutine wd_weight_serial(mode,nfldos)
      integer, intent(in) :: mode, nfldos
      integer :: m, ik, ib, ib_ordr, kend, kstep, iws, ip
      real(kind=DP),allocatable,dimension(:):: wei

      if(mode==ALDOS) then
         m = maldos
      else if(mode==LAYERDOS) then
         m = mlayer
      end if

      allocate(wei(m))

      if(ekmode==OFF) then
         kend = kv3; kstep = af+1
      else if(ekmode==ON) then
!         kend = 1;   kstep = 1
         kend = kv3;   kstep = 1
      end if

      if(ipridos >= 3) then
         do ik = 1, kv3, af+1
            if(map_k(ik) /= myrank_k) cycle
            write(nfout,'(" !ldos: ik = ",i8)') ik
            do ib = 1, neg
               write(nfout,'(" !ldos: ib, ib_ordr = ",2i20)') ib, neordr(ib,ik)
            end do
         end do
      end if

      do ik = 1, kend, kstep
         if(mype == 0) then
            if(mode==ALDOS) then
               write(nfldos,'(" weight for each atomic cell    nk = ",i8)') ik
            else if(mode==LAYERDOS) then
               write(nfldos,'(" weight for each layer          nk = ",i8)') ik
            end if
         end if

         do ib = 1, neg -num_extra_bands

             if(mode == ALDOS) then
                wei(1:m) = weiwsc(1:m,ib,ik)
             else if(mode == LAYERDOS) then
                wei(1:m) = weilay(1:m,ib,ik)
             end if
             if(mype == 0)  then
               if(ilen >= m) then
                  write(nfldos,'(i4,")",5f15.10)') ib,(wei(iws),iws=1,m)
               else
                  write(nfldos,'(i4,")",5f15.10)') ib,(wei(iws),iws=1,ilen)
                  write(nfldos,'(5x,5f15.10)') (wei(iws),iws=ilen+1,m)
               end if
             end if
         end do
      end do

       deallocate(wei)
     end subroutine wd_weight_serial

    subroutine check_sum(ik,aword,serial,norma)
      integer,intent(in) :: ik
      character*(*),intent(in) :: aword
      logical, intent(in) :: serial
      logical, intent(in) :: norma

      integer :: ib, ia, ilay, ib_ordr, ik_t
      real(kind=DP), allocatable, dimension(:) :: sum_aldos, sum_layer, sum_mpi


      if(sw_atom_decomp==ON) allocate(sum_aldos(neg))
      if(sw_layer_decomp==ON) allocate(sum_layer(neg))

      if ( sw_lband == ON ) then
         ik_t = ik
      else
         if(ekmode /= ON) then
            ik_t = ik
         else if(ekmode == ON) then
            ik_t = 1
         end if
      endif

      if(sw_atom_decomp==ON .or. sw_layer_decomp==ON) then
         if(nrank_e >= 2) then
            allocate(sum_mpi(neg)); sum_mpi = 0.d0
         end if
      end if
      if(sw_atom_decomp==ON) then
         sum_aldos = 0.d0
         do ib = 1, neg
            if(.not.serial)then
               if(map_ek(ib,ik_t) == mype) then
                  do ia = 1, natm2 + 1
                     sum_aldos(ib) = sum_aldos(ib) + weiwsc(ia,map_z(ib),ik_t)
                  end do
               endif
            else
               do ia = 1, natm2 + 1
                  sum_aldos(ib) = sum_aldos(ib) + weiwsc(ia,ib,ik_t)
               end do
            end if
         end do
         if(.not.serial)then
            if(nrank_e >= 2) then
               sum_mpi = 0.d0
               call mpi_allreduce(sum_aldos,sum_mpi,neg,mpi_double_precision &
                    & , mpi_sum, mpi_k_world(myrank_k),ierr)
               sum_aldos = sum_mpi
            end if
         endif
      end if

      if(sw_layer_decomp==ON) then
         sum_layer = 0.d0
         do ib = 1, neg
            if(.not.serial)then
               if(map_e(ib) == myrank_e) then
                  do ilay = 1, mlayer
                     sum_layer(ib) = sum_layer(ib) + weilay(ilay,map_z(ib),ik_t)
                  end do
               endif
            else
               do ilay = 1, mlayer
                  sum_layer(ib) = sum_layer(ib) + weilay(ilay,ib,ik_t)
               end do
            end if
         end do
         if(.not.serial)then
            if(nrank_e >= 2) then
               sum_mpi = 0.d0
               call mpi_allreduce(sum_layer,sum_mpi,neg,mpi_double_precision &
                    & , mpi_sum, mpi_k_world(myrank_k),ierr)
               sum_layer = sum_mpi
            end if
         endif
      end if

!!$      if(mype == 0) then
      if(ipridos>=1) then
         if(sw_atom_decomp==ON .and. sw_layer_decomp==ON) then
            write(nfout,'(" !!ldos ",a8," -- iban , sum(atomic, layer) --  neg = ",i6 &
                 & ," natm2 = ",i6, " ik = ",i6)') aword,neg,natm2,ik
            do ib = 1, neg
               if(serial) then
                  ib_ordr = ib
               else
                  ib_ordr = neordr(ib,ik)
               endif
               write(nfout,'(" !!ldos ", i4,2f16.8)') ib, sum_aldos(ib_ordr), sum_layer(ib_ordr)
            end do
         else if(sw_atom_decomp==ON) then
            write(nfout,'(" !!ldos ",a8," -- iban , sum(atomic) --  neg = ",i6 &
                 & ," natm2 = ",i6," ik = ",i6)') aword,neg,natm2,ik
            do ib = 1, neg
               if(serial) then
                  ib_ordr = ib
                  write(nfout,'(" !!ldos ", i4,f16.8)') ib, sum_aldos(ib_ordr)
               else
                  if(map_ek(ib,ik_t) == mype) then
                    ib_ordr = neordr(ib,ik)
                    write(nfout,'(" !!ldos ", i4,f16.8)') ib, sum_aldos(ib_ordr)
                  endif
               endif
            end do
         else if(sw_layer_decomp == ON) then
            write(nfout,'(" !!ldos ",a8," -- iban , sum(layer ) --  neg = ",i6 &
                 & ," natm2 = ",i6," ik = ",i6)') aword,neg,natm2,ik
            do ib = 1, neg
               if(serial) then
                  ib_ordr = ib
               else
                  ib_ordr = neordr(ib,ik)
               endif
               write(nfout,'(" !!ldos ",i4,f16.8)') ib, sum_layer(ib_ordr)
            end do
         end if
      end if
      if (norma .and. sw_atom_decomp==ON .and. sw_ac_mesh==ON) then
        if(.not.serial)then
           do ib=ista_e,iend_e,istep_e
              weiwsc(:,map_z(ib),ik_t) = weiwsc(:,map_z(ib),ik_t)/sum_aldos(ib)
           enddo
        else
           do ib=1,neg
              weiwsc(:,ib,ik_t) = weiwsc(:,ib,ik_t)/sum_aldos(ib)
           enddo
        endif
      endif

      if(sw_atom_decomp==ON) deallocate(sum_aldos)
      if(sw_layer_decomp==ON) deallocate(sum_layer)
      if(allocated(sum_mpi)) deallocate(sum_mpi)
!!$      if(ekmode /= ON .and. (sw_aldos == ON .or. sw_layerdos == ON) .and. npes >= 2) deallocate(sum_mpi)

    end subroutine check_sum

    subroutine substitute_weilay(ik,ib,serial)
!   *********************************
!     by T.Yamasaki
!           10th Jun 1992
!   *********************************
!     Revised by T. Yamasaki, Feb. 2004
!
      integer, intent(in) :: ik,ib
      logical, intent(in) :: serial

      integer :: i, ilay, iax1, iax2, iax3, ib_t
      real(kind=DP) :: a1,a2,width, g,ss
      complex(kind=CMPLDP) :: zi,zsum,zchg
      real(kind=DP),allocatable,dimension(:) :: zsum_layer, zsum_layer_mpi

      allocate(zsum_layer(mlayer)); zsum_layer = 0.d0

!!$      zi = cmplx(0.d0,1.d0)
      zi = cmplx(0.d0,1.d0)*PAI2

      iax1 = normal_axis_winlay
      iax2 = iax1 + 1
      if(iax2 >= 4) iax2 = iax2 - 3
      iax3 = iax2 + 1
      if(iax3 >= 4) iax3 = iax3 - 3


!!$      rl11 = rltv(iax1,iax1)
!!$      ss = 1.d0/maxhv
      ss = 3.d0 - kimg

!!$      call get_height(height)

!!$      n = fft_box_size_WF(1,0)*fft_box_size_WF(2,0)*fft_box_size_WF(3,1)
!!$      n = nfft*kimg/kimg/2
      do i = 1, nfft/2
         bfft(2*i-1) = bfft(2*i-1)**2 + bfft(2*i)**2
         bfft(2*i  ) = 0.d0
      end do
      call m_FFT_WF(ELECTRON,nfout,bfft,DIRECT,ON)
      if(kimg == 1) then
         do i = 1, kg
            zrhoik(i) = bfft(igf(i))*denom
         end do
      else
         do i = 1, kg
            zrhoik(i*2-1) = bfft(igf(i)*2-1)*denom
            zrhoik(i*2  ) = bfft(igf(i)*2  )*denom
         end do
      end if

      do ilay = 1, mlayer
!!$         a2 = winlay(ilay,2)
!!$         a1 = winlay(ilay,1)
         a2 = winlay(ilay,2)/height
         a1 = winlay(ilay,1)/height
         width = a2 - a1
         zsum = cmplx(zrhoik(1)*width,0.d0)
!!$         zsum = univol*cmplx(zrhoik(1)*width,0.d0)
         if(kimg == 1) then
            do i = 2, kg
               if(ngabc(i,iax2) == 0 .and. ngabc(i,iax3) == 0) then
!!$                  g = rl11*ngabc(i,iax1)
                  g = ngabc(i,iax1)
                  zsum = zsum  &
                       & + zrhoik(i)*((cdexp(zi*g*a2)-cdexp(zi*g*a1))/( zi*g ))
               end if
            end do
         else if(kimg == 2) then
            do i = 2, kg
               if(ngabc(i,iax2) ==  0 .and. ngabc(i,iax3) == 0) then
!!$                  g = rl11*ngabc(i,iax1)
                  g = ngabc(i,iax1)
                  zchg = dcmplx(zrhoik(i*2-1),zrhoik(i*2))
                  zsum = zsum  &
                       & + zchg*((cdexp(zi*g*a2)-cdexp(zi*g*a1))/( zi*g ))
               endif
            enddo
         endif
         zsum_layer(ilay) = real(zsum)*ss
      end do

      if ( sw_lband == ON ) then
         ib_t = map_z(ib)
      else
         if(ekmode == ON.or.serial) then
            ib_t = ib
         else
            ib_t = map_z(ib)
         end if
      endif

      do ilay = 1, mlayer
         weilay(ilay,ib_t,ik) = zsum_layer(ilay)
      end do

      deallocate(zsum_layer)

    end subroutine substitute_weilay

    subroutine substitute_weilay_fftmesh(ik,ib,serial)
!   *********************************
!     by T.Yamasaki
!           7th Apr 2017
!   *********************************
!
      integer, intent(in) :: ik,ib
      logical, intent(in) :: serial
      integer :: ijk, nwsc,n, nel, ib_t
      real(kind=DP) :: xx, denom
      real(kind=DP),allocatable,dimension(:) :: adjustfactor ! d(mlayer)
      integer, save :: nwritecount = 0

      nel = product(fft_box_size_WF(1:3,1))
      denom = 1.d0/dble(nel)
      if(ipridos >= 2) write(nfout,'(" !!ldos: ib, nfft, natm2 = ",i6,i12,i6,"  <<substitute_weiwsc>>")') &
           & ib, nfft, natm2

      if ( sw_lband == ON ) then
         ib_t = map_z(ib)
      else
         if(ekmode == ON.or.serial) then
            ib_t = ib
         else if(ekmode == OFF) then
            ib_t = map_z(ib)
         end if
      endif

      allocate(adjustfactor(mlayer))
      call get_adjustfactor(winlay,mlayer,nmeshlay,adjustfactor,nwritecount) ! winlay,height,nmeshl -> adustfactor
      if(nwritecount==0) nwritecount=1
      do ijk = 1, nel
         nwsc = meshwk_layer(ijk,1)
         n    = meshwk_layer(ijk,2)
         if((n < 1 .or. 2*n > nfft) .or. (nwsc < 1 .or. nwsc >mlayer+1) ) then
            if(ipridos>=1) then
               write(nfout,'(" !!ldos: n = ",i6," ijk = ",i6)') n, ijk
               write(nfout,'(" !!ldos: nwsc = ",i9)') nwsc
            end if
         else
!!$            xx = denom*(bfft(2*n-1)**2 + bfft(2*n)**2)
            xx = denom*(bfft(2*n-1)**2 + bfft(2*n)**2)*adjustfactor(nwsc)
            weilay(nwsc,ib_t,ik) = weilay(nwsc,ib_t,ik) + xx
         end if
      end do
      deallocate(adjustfactor)
    end subroutine substitute_weilay_fftmesh

    subroutine substitute_weiwsc(ik,ib,serial)
      integer, intent(in) :: ik,ib
      logical, intent(in) :: serial
!   *********************************
!     by T.Yamasaki
!           10th Jun 1992
!   *********************************
!     Revised by T. Yamasaki, Feb. 2004
!
      integer :: ijk, nwsc,n, nel, ib_t
      real(kind=DP) :: xx, denom

      nel = product(fft_box_size_WF(1:3,1))
      denom = 1.d0/dble(nel)
      if(ipridos >= 2) write(nfout,'(" !!ldos: ib, nfft, natm2 = ",i6,i12,i6,"  <<substitute_weiwsc>>")') &
           & ib, nfft, natm2

      if ( sw_lband == ON ) then
         ib_t = map_z(ib)
      else
         if(ekmode == ON.or.serial) then
            ib_t = ib
         else if(ekmode == OFF) then
            ib_t = map_z(ib)
         end if
      endif

      do ijk = 1, nel
         nwsc = meshwk(ijk,1)
         n    = meshwk(ijk,2)
         if((n < 1 .or. 2*n > nfft) .or. (nwsc < 1 .or. nwsc >natm2+1) ) then
            if(ipridos>=1) then
               write(nfout,'(" !!ldos: n = ",i6," ijk = ",i6)') n, ijk
               write(nfout,'(" !!ldos: nwsc = ",i9)') nwsc
            end if
         else
            xx = denom*(bfft(2*n-1)**2 + bfft(2*n)**2)
            weiwsc(nwsc,ib_t,ik) = weiwsc(nwsc,ib_t,ik) + xx
         end if
      end do

      if(ipridos >= 2) then
         write(nfout,'(" !ldos: ik, ib = ",2i8," <<substitute_weiwsc>>")') ik, ib
         write(nfout,'(" !ldos: ",8f10.5)') (weiwsc(nwsc,ib_t,ik),nwsc=1,maldos)
      end if

    end subroutine substitute_weiwsc

    subroutine substitute_weiwsc_acmesh(ik,ib,serial,cd,wsc)
      integer, intent(in) :: ik,ib
      logical, intent(in) :: serial,cd
      real(kind=DP), intent(inout), dimension(maldos,neg,kv3) :: wsc
      integer :: ijk, nwsc,n, ib_t
      real(kind=DP) :: xx, denom
      integer :: i,j,k, id, nl, nm, nn, nlhf,inew,jnew,knew,ip,mm,is,nl2,nm2,nn2
      integer :: ii,jj,kk,iatm
      real(kind=DP) :: ri,rj,rk
      real(kind=DP), allocatable, dimension(:) :: xmat,ymat,zmat
      real(kind=DP), allocatable, dimension(:,:,:) :: cmatr
      real(kind=DP) :: rho
      real(kind=DP) :: trilinear_interpolation,closest
      integer, pointer, dimension(:,:) :: amesh
      integer, pointer, dimension(:) :: namesh
      integer :: ista,iend
      integer :: id_sname = -1
      call tstatc0_begin('substitute_weiwsc_acmesh ',id_sname,level=1)
      if(printable .and. ib==1)then
        if(.not.cd) then
          write(nfout,'(a,i10,a)') ' !** AC MESH processing kpoint no. ',ik,' (softpart)'
        else
          write(nfout,'(a,i10,a)') ' !** AC MESH processing kpoint no. ',ik,' (hardpart)'
        endif
      endif
      if(cd)then
        id = fft_box_size_CD_nonpara(1,0)
        mm = fft_box_size_CD_nonpara(2,0)
        nl = fft_box_size_CD(1,1)
        nm = fft_box_size_CD(2,1)
        nn = fft_box_size_CD(3,1)
        ista = ista_atm2
        iend = iend_atm2
      else
        id = fft_box_size_WF(1,0)
        mm = fft_box_size_WF(2,0)
        nl = fft_box_size_WF(1,1)
        nm = fft_box_size_WF(2,1)
        nn = fft_box_size_WF(3,1)
        ista = 1
        iend = natm2
      endif
      nl2=nl+2;nm2=nm+2;nn2=nn+2
      allocate(xmat(nl2));allocate(ymat(nm2));allocate(zmat(nn2))
      allocate(cmatr(nl2,nm2,nn2))
      if(kimg == 1) then
         nlhf = id/2
      else
         nlhf = id
      end if
      if(cd) then
        do i = 0, nm+1
          do j = 0, nn+1
            do k = 0, nl+1
              ii=i;jj=j;kk=k
              if(i.eq.0) ii = nm
              if(j.eq.0) jj = nn
              if(k.eq.0) kk = nl
              if(i.eq.nm+1) ii = 1
              if(j.eq.nn+1) jj = 1
              if(k.eq.nl+1) kk = 1
              if(kimg == 1 .and. kk > nlhf) then
                 knew = id - kk
                 jnew = nn+2 - jj
                 inew = nm+2 - ii
                 if(jnew > nn) then
                    jnew = jnew - nn
                 end if
                 if(inew > nm) then
                    inew = inew - nm
                 end if
              else
                 knew = kk; jnew = jj; inew = ii
              end if
              ip = nlhf*mm*(jnew-1) + nlhf*(inew-1) + knew
              xmat(k+1) = dble(k)/dble(nl)
              ymat(i+1) = dble(i)/dble(nm)
              zmat(j+1) = dble(j)/dble(nn)
              cmatr(k+1,i+1,j+1) = bfft(ip*2-1)
             enddo
          enddo
        enddo
      else
        do i = 0, nm+1
          do j = 0, nn+1
            do k = 0, nl+1
              ii=i;jj=j;kk=k
              if(i.eq.0) ii = nm
              if(j.eq.0) jj = nn
              if(k.eq.0) kk = nl
              if(i.eq.nm+1) ii = 1
              if(j.eq.nn+1) jj = 1
              if(k.eq.nl+1) kk = 1
              if(kimg == 1 .and. kk > nlhf) then
                 knew = id - kk
                 jnew = nn+2 - jj
                 inew = nm+2 - ii
                 if(jnew > nn) then
                    jnew = jnew - nn
                 end if
                 if(inew > nm) then
                    inew = inew - nm
                 end if
              else
                 knew = kk; jnew = jj; inew = ii
              end if
              ip = nlhf*mm*(jnew-1) + nlhf*(inew-1) + knew
              xmat(k+1) = dble(k)/dble(nl)
              ymat(i+1) = dble(i)/dble(nm)
              zmat(j+1) = dble(j)/dble(nn)
              cmatr(k+1,i+1,j+1) = bfft(ip*2-1)**2+bfft(ip*2)**2
            enddo
          enddo
        enddo
      endif

      if ( sw_lband == ON ) then
         ib_t = map_z(ib)
      else
         ib_t = map_z(ib)
         if(ekmode == ON.or.serial) ib_t = ib
      endif

      if(cd) then
        namesh => nac_mesh_cd
      else
        namesh => nac_mesh
      endif

      if(cd) then
        denom = 0.5d0*univol/dble(product(fft_box_size_CD(1:3,1)))/dble(acmesh_factor*acmesh_factor*acmesh_factor)
      else
        denom = 1.d0/dble(product(fft_box_size_WF(1:3,1)))/dble(acmesh_factor*acmesh_factor*acmesh_factor)
      endif
      do iatm=ista,iend
         if(cd)then
           amesh => ac_mesh_cd(iatm,:,:)
         else
           amesh => ac_mesh(iatm,:,:)
         endif
         do ijk=1,namesh(iatm)
            ii = amesh(ijk,1)
            jj = amesh(ijk,2)
            kk = amesh(ijk,3)
            ri = dble(ii+acmesh_factor)/dble(nl*acmesh_factor)+pos(iatm,1)
            rj = dble(jj+acmesh_factor)/dble(nm*acmesh_factor)+pos(iatm,2)
            rk = dble(kk+acmesh_factor)/dble(nn*acmesh_factor)+pos(iatm,3)
            if(ri < 0.d0) ri = (ceiling(-ri)) + ri
            if(ri > 1.d0) ri = ri - floor(ri)
            if(rj < 0.d0) rj = (ceiling(-rj)) + rj
            if(rj > 1.d0) rj = rj - floor(rj)
            if(rk < 0.d0) rk = (ceiling(-rk)) + rk
            if(rk > 1.d0) rk = rk - floor(rk)
            rho = trilinear_interpolation(nl2,nm2,nn2,cmatr,xmat,ymat,zmat,ri,rj,rk)
            xx = denom*rho
            wsc(iatm,ib_t,ik) = wsc(iatm,ib_t,ik)+xx
         enddo
      enddo

      deallocate(xmat)
      deallocate(ymat)
      deallocate(zmat)
      deallocate(cmatr)
      call tstatc0_end(id_sname)
    end subroutine substitute_weiwsc_acmesh

    subroutine substitute_weiwsc_cd(ik,ib,serial)
      integer, intent(in) :: ik,ib
      logical, intent(in) :: serial
!   *********************************
!     by T.Yamasaki
!           10th Jun 1992
!   *********************************
!     Revised by T. Yamasaki, Feb. 2004
!
      integer :: ijk, nwsc,n, ib_t
      real(kind=DP) :: denom

      denom = 0.5d0*univol/product(fft_box_size_CD(1:3,1))

      if ( sw_lband == ON ) then
         ib_t = map_z(ib)
      else
         if(ekmode == ON.or.serial) then
            ib_t = ib
         else if(ekmode == OFF) then
            ib_t = map_z(ib)
         end if
      endif

      do ijk = 1, product(fft_box_size_CD(1:3,1))
         nwsc = meshwk(ijk,1)
         n    = meshwk(ijk,2)
         if((n < 1 .or. 2*n > nfftp_nonpara) .or. (nwsc < 1 .or. nwsc >natm2+1) ) then
            if(ipridos>=1) then
               write(nfout,'(" !!ldos: n = ",i6," ijk = ",i6)') n, ijk
               write(nfout,'(" !!ldos: nwsc = ",i9)') nwsc
            end if
         else
            weiwsc(nwsc,ib_t,ik) = weiwsc(nwsc,ib_t,ik) + denom*bfft(2*n-1)
! =============================== modified by K. Tagami ============ 11.0
!            weiwsc(nwsc,ib_t,ik) = weiwsc(nwsc,ib_t,ik) + denom*bfft(2*n-1)
!
!            if (map_ek(ib,ik) == mype) then
!               weiwsc(nwsc,ib_t,ik) = weiwsc(nwsc,ib_t,ik) + denom*bfft(2*n-1)
!            endif
! ================================================================== 11.0
         end if
      end do

      if(ipridos >=2) then
         write(nfout,'(" !ldos: ik ,ib = ",2i8)') ik, ib
         write(nfout,'(" !ldos: ",8f10.5)') (weiwsc(nwsc,ib_t,ik),nwsc=1,maldos)
      end if

    end subroutine substitute_weiwsc_cd

    subroutine substitute_weilay_cd(ik,ib,chgq0,serial)
      integer, intent(in) :: ik, ib
      real(kind=DP), intent(in) :: chgq0
      logical, intent(in) :: serial
!   *********************************
!     by T.Yamasaki
!           10th Jun 1992
!   *********************************
!     Revised by T. Yamasaki, Feb. 2004
!
      integer :: i, ista, iax1, iax2, iax3, ilay,  ib_t
      real(kind=DP) :: g, z_work, ss, a1, a2, width
      complex(kind=CMPLDP) :: zi,zsum,zchg
      real(kind=DP),allocatable,dimension(:,:) :: zsum_layer, zsum_layer_mpi

      allocate(zsum_layer(mlayer,2)); zsum_layer = 0.d0

      zi = cmplx(0.d0,1.d0)*PAI2

      iax1 = normal_axis_winlay
      iax2 = iax1 + 1
      if(iax2 >= 4) iax2 = iax2 - 3
      iax3 = iax2 + 1
      if(iax3 >= 4) iax3 = iax3 - 3

!!$      rl11 = rltv(iax1,iax1)
!!$      ss = 0.5d0*univol/maxhv
      ss = 0.5d0*univol*(3.d0-kimg)


      ista = ista_kngp
      if(ista == 1) ista = 2
      do ilay = 1, mlayer
         if(ipridos >= 2) write(nfout,'(" !!ldos: ib = ",i6," ilay = ",i6 &
              & ," <<substitute_weilay_cd>>")') ib,ilay
         a2 = winlay(ilay,2)/height
         a1 = winlay(ilay,1)/height
         width = a2 - a1
         z_work = chgq0*width
         if(serial)then
            if(ista/=2) z_work=0.d0
         endif
         zsum = 0.d0

         if(kimg == 1) then
            if ( sw_lband == ON ) then
               do i = 2, kgp
                  if(ngabc(i,iax2) == 0 .and. ngabc(i,iax3) == 0) then
                     g = ngabc(i,iax1)
                     zsum = zsum  &
                          & + chgq_enl(i,1)*( ( cdexp(zi*g*a2) - cdexp(zi*g*a1))/(zi*g))
                  end if
               end do
            else
               if(ekmode == ON.or.serial) then
                  do i = ista, iend_kngp
                     if(ngabc(i,iax2) == 0 .and. ngabc(i,iax3) == 0) then
                        g = ngabc(i,iax1)
                        zsum = zsum  &
                             & + chgq_l(i,1,1)*((cdexp(zi*g*a2) -cdexp(zi*g*a1))/(zi*g))
                     end if
                  end do
               else if(ekmode == OFF) then
                  do i = 2, kgp
                     if(ngabc(i,iax2) == 0 .and. ngabc(i,iax3) == 0) then
                        g = ngabc(i,iax1)
                        zsum = zsum  &
                          & + chgq_enl(i,1)*( ( cdexp(zi*g*a2) - cdexp(zi*g*a1))/(zi*g))
                     end if
                  end do
               end if
            endif
         else if(kimg == 2) then
            if ( sw_lband == ON ) then
               do i = 2, kgp
                  if(ngabc(i,iax2) == 0 .and. ngabc(i,iax3) == 0) then
                     g = ngabc(i,iax1)
                     zchg = cmplx(chgq_enl(i,1),chgq_enl(i,2))
                     zsum = zsum &
                          & + zchg * ((cdexp(zi*g*a2) - cdexp(zi*g*a1))/(zi*g))
                  end if
               end do
            else
               if(ekmode == ON.or.serial) then
                  do i = ista, iend_kngp
                     if(ngabc(i,iax2) == 0 .and. ngabc(i,iax3) == 0) then
                        g = ngabc(i,iax1)
                        zchg = cmplx(chgq_l(i,1,1),chgq_l(i,2,1))
                        zsum = zsum &
                             & + zchg * ((cdexp(zi*g*a2) - cdexp(zi*g*a1))/(zi*g))
                     end if
                  end do
               else if(ekmode == OFF) then
                  do i = 2, kgp
                     if(ngabc(i,iax2) == 0 .and. ngabc(i,iax3) == 0) then
                        g = ngabc(i,iax1)
                        zchg = cmplx(chgq_enl(i,1),chgq_enl(i,2))
                        zsum = zsum &
                             & + zchg * ((cdexp(zi*g*a2) - cdexp(zi*g*a1))/(zi*g))
                     end if
                  end do
               end if
            endif
         end if
         zsum = zsum + cmplx(z_work,0.d0)

         zsum_layer(ilay,1) = real(zsum)
         zsum_layer(ilay,2) = dimag(zsum)

      end do

#if 0
      if(ekmode == ON.or.serial) then
         if(npes > 1) then
            allocate(zsum_layer_mpi(mlayer,2))
            call mpi_allreduce(zsum_layer,zsum_layer_mpi,mlayer*2 &
                 & ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
            zsum_layer = zsum_layer_mpi
            deallocate(zsum_layer_mpi)
         end if
      end if

      if(mype == 0) then
         do ilay = 1, mlayer
            if(zsum_layer(ilay,2) > 1.d-6) then
               write(6,'(" !!ldos: zsum,ne.real imag = ",d16.8," ilay = ",i6)') &
                    & zsum_layer(ilay,2), ilay
            endif
!!$            write(6,'(" !!ldos zsum_layer = ",d16.8," ilay = ",i6)') &
!!$                 & zsum_layer(ilay,1), ilay
         end do
      end if
#endif

      if ( sw_lband == ON ) then
         ib_t = map_z(ib)
      else
         if(ekmode == ON.or.serial) then
            ib_t = ib
         else
            ib_t = map_z(ib)
         end if
      endif

      do ilay = 1, mlayer
         weilay(ilay,ib_t,ik) = weilay(ilay,ib_t,ik) + zsum_layer(ilay,1)*ss
      end do
! ========================= modified by K. Tagami =============== 11.0
!      do ilay = 1, mlayer
!         weilay(ilay,ib_t,ik) = weilay(ilay,ib_t,ik) + zsum_layer(ilay,1)*ss
!      end do
!
!      if (map_ek(ib,ik) == mype) then
!         do ilay = 1, mlayer
!            weilay(ilay,ib_t,ik) = weilay(ilay,ib_t,ik) + zsum_layer(ilay,1)*ss
!         end do
!      endif
! ============================================================== 11.0

      deallocate(zsum_layer)

    end subroutine substitute_weilay_cd

    subroutine substitute_weilay_cd_fftmesh(ik,ib,serial)
      integer, intent(in) :: ik, ib
      logical, intent(in) :: serial
!   *********************************
!     by T.Yamasaki
!           7th Apr 2017
!   *********************************
!
      integer :: ijk, nwsc,n, ib_t
      real(kind=DP) :: denom
      real(kind=DP),allocatable,dimension(:) :: adjustfactor ! d(mlayer)
      integer, save :: nwritecount = 0

      denom = 0.5d0*univol/product(fft_box_size_CD(1:3,1))

      if ( sw_lband == ON ) then
         ib_t = map_z(ib)
      else
         if(ekmode == ON.or.serial) then
            ib_t = ib
         else if(ekmode == OFF) then
            ib_t = map_z(ib)
         end if
      endif

      allocate(adjustfactor(mlayer))
      call get_adjustfactor(winlay,mlayer,nmeshplay,adjustfactor,nwritecount) ! winlay,height,nmeshl -> adustfactor
      if(nwritecount==0) nwritecount=1

      do ijk = 1, product(fft_box_size_CD(1:3,1))
         nwsc = meshwk_layer(ijk,1)
         n    = meshwk_layer(ijk,2)
         if((n < 1 .or. 2*n > nfftp_nonpara) .or. (nwsc < 1 .or. nwsc >mlayer+1) ) then
            if(ipridos>=1) then
               write(nfout,'(" !!ldos: n = ",i6," ijk = ",i6)') n, ijk
               write(nfout,'(" !!ldos: nwsc = ",i9)') nwsc
            end if
         else
!!$            weilay(nwsc,ib_t,ik) = weilay(nwsc,ib_t,ik) + denom*bfft(2*n-1)
            weilay(nwsc,ib_t,ik) = weilay(nwsc,ib_t,ik) + denom*bfft(2*n-1)*adjustfactor(nwsc)
         end if
      end do
      deallocate(adjustfactor)
    end subroutine substitute_weilay_cd_fftmesh

    subroutine fillup_meshwk(mode)
!   *********************************
!     by T.Yamasaki
!           10th Jun 1992
!   *********************************
!     Revised by T. Yamasaki, Feb. 2004
!
      integer, intent(in) :: mode
!!$ASASASASAS
!!$      integer, pointer, dimension(:,:,:) :: mesh_t
      integer, allocatable :: mesh_t(:,:,:)
!!$ASASASASAS
      integer :: i,j,k, icount, nwsc, max_nwsc, min_nwsc, n, nl, nm,nn, nd2,id, nfftwk &
           &    , inew, jnew, knew, idr2, idh

      if(mode == SOFTPART) then

         nl = fft_box_size_WF(1,1)
         nm = fft_box_size_WF(2,1)
         nn = fft_box_size_WF(3,1)

         id  = fft_box_size_WF(1,0)
         nd2 = fft_box_size_WF(2,0)
         allocate(mesh_t(id,nd2,nn))
         mesh_t = mesh

         nfftwk = product(fft_box_size_WF(1:3,1))
      else if(mode == HARDPART) then

         nl = fft_box_size_CD(1,1)
         nm = fft_box_size_CD(2,1)
         nn = fft_box_size_CD(3,1)

         id  = fft_box_size_CD_nonpara(1,0)
         nd2 = fft_box_size_CD_nonpara(2,0)
         allocate(mesh_t(id,nd2,nn))
         mesh_t = meshp
         nfftwk = product(fft_box_size_CD(1:3,1))
      end if
      idr2 = (nl+2)/2
      idh  = id/2

      if(allocated(meshwk)) deallocate(meshwk)
      allocate(meshwk(nfftwk,2)); meshwk = 0

      if(ipridos>=2) write(nfout,'(" !! nfftwk = ",i8)') nfftwk

      max_nwsc = 0
      min_nwsc = natm2+2
      icount = 1
      do k = 1, nn
         do j = 1, nm
            do i = 1, nl
!!$               if(icount > nfftwk) then
!!$                  write(nfout,'(" ! icount = ",i7," > nfftwk")') icount
!!$                  stop  ' large icount <<fillup_meshwk.m_Ldos_cal>>'
!!$               end if
               nwsc = mesh_t(i,j,k)
               if(nwsc > max_nwsc) max_nwsc = nwsc
               if(nwsc < min_nwsc) min_nwsc = nwsc
               meshwk(icount,1) = nwsc
               if(kimg == 1) then
                  if( i <= idr2 ) then
                     n = idh*nd2*(k-1) + idr2*(j-1) + i
                  else
                     inew = nl + 2 - i
                     jnew = nm + 2 - j
                     knew = nn + 2 - k
                     if(jnew > nm) jnew = jnew - nm
                     if(knew > nn) knew = knew - nn
                     n = idh*nd2*(knew-1) + idh*(jnew-1) + inew
                  end if
               else
                  n = id*nd2*(k-1) + id*(j-1) + i
               end if
               meshwk(icount,2) = n
               icount = icount + 1
            end do
         end do
      end do
      if(ipridos>=2) then
         write(nfout,*) ' !!ldos:  max_nwsc =  ', max_nwsc
         write(nfout,*) ' !!ldos:  min_nwsc =  ', min_nwsc
#ifdef DEBUG_MESH
         n = 1000
         if(n*2 > nfftwk) n = nfftwk/2
         write(nfout,'(" !ldos nfftwk = ",i8)') nfftwk
         write(nfout,'(" !ldos: meshwk(*,1) first ",i8," elements")') n
         write(nfout,'(14i5)') (meshwk(i,1),i=1,n)
         write(nfout,'(" !ldos: meshwk(*,1) last  ",i8," elements")') n
         write(nfout,'(14i6)') (meshwk(i,1),i=nfftwk-n+1, nfftwk)

         write(nfout,'(" !ldos: meshwk(*,2) first ",i8," elements")') n
         write(nfout,'(14i5)') (meshwk(i,2),i=1,n)
         write(nfout,'(" !ldos: meshwk(*,2) last  ",i8," elements")') n
         write(nfout,'(14i6)') (meshwk(i,2),i=nfftwk-n+1, nfftwk)
#endif
      end if
      deallocate(mesh_t)
    end subroutine fillup_meshwk

    subroutine fillup_meshwk_layer(mode)
!   *********************************
!     by T.Yamasaki
!           7th Apr 2017
!   *********************************

      integer, intent(in) :: mode
      integer, allocatable :: mesh_t(:,:,:)
      integer :: i,j,k, icount, nwsc, max_nwsc, min_nwsc, n, nl, nm,nn, nd2,id, nfftwk &
           &    , inew, jnew, knew, idr2, idh

      if(mode == SOFTPART) then

         nl = fft_box_size_WF(1,1)
         nm = fft_box_size_WF(2,1)
         nn = fft_box_size_WF(3,1)

         id  = fft_box_size_WF(1,0)
         nd2 = fft_box_size_WF(2,0)
         allocate(mesh_t(id,nd2,nn))
         mesh_t = meshl

         nfftwk = product(fft_box_size_WF(1:3,1))
      else if(mode == HARDPART) then

         nl = fft_box_size_CD(1,1)
         nm = fft_box_size_CD(2,1)
         nn = fft_box_size_CD(3,1)

         id  = fft_box_size_CD_nonpara(1,0)
         nd2 = fft_box_size_CD_nonpara(2,0)
         allocate(mesh_t(id,nd2,nn))
         mesh_t = meshpl
         nfftwk = product(fft_box_size_CD(1:3,1))
      end if
      idr2 = (nl+2)/2
      idh  = id/2

      if(allocated(meshwk_layer)) deallocate(meshwk_layer)
      allocate(meshwk_layer(nfftwk,2)); meshwk_layer = 0

      if(ipridos>=2) write(nfout,'(" !! nfftwk = ",i8)') nfftwk

      max_nwsc = 0
      min_nwsc = mlayer+1
      icount = 1
      do k = 1, nn
         do j = 1, nm
            do i = 1, nl
!!$               if(icount > nfftwk) then
!!$                  write(nfout,'(" ! icount = ",i7," > nfftwk")') icount
!!$                  stop  ' large icount <<fillup_meshwk.m_Ldos_cal>>'
!!$               end if
               nwsc = mesh_t(i,j,k)
               if(nwsc > max_nwsc) max_nwsc = nwsc
               if(nwsc < min_nwsc) min_nwsc = nwsc
               meshwk_layer(icount,1) = nwsc
               if(kimg == 1) then
                  if( i <= idr2 ) then
                     n = idh*nd2*(k-1) + idr2*(j-1) + i
                  else
                     inew = nl + 2 - i
                     jnew = nm + 2 - j
                     knew = nn + 2 - k
                     if(jnew > nm) jnew = jnew - nm
                     if(knew > nn) knew = knew - nn
                     n = idh*nd2*(knew-1) + idh*(jnew-1) + inew
                  end if
               else
                  n = id*nd2*(k-1) + id*(j-1) + i
               end if
               meshwk_layer(icount,2) = n
               icount = icount + 1
            end do
         end do
      end do
      if(ipridos>=2) then
         write(nfout,*) ' !!ldos:  max_nwsc =  ', max_nwsc
         write(nfout,*) ' !!ldos:  min_nwsc =  ', min_nwsc
#ifdef DEBUG_MESH
         n = 1000
         if(n*2 > nfftwk) n = nfftwk/2
         write(nfout,'(" !ldos nfftwk = ",i8)') nfftwk
         write(nfout,'(" !ldos: meshwk_layer(*,1) first ",i8," elements")') n
         write(nfout,'(14i5)') (meshwk_layer(i,1),i=1,n)
         write(nfout,'(" !ldos: meshwk_layer(*,1) last  ",i8," elements")') n
         write(nfout,'(14i6)') (meshwk_layer(i,1),i=nfftwk-n+1, nfftwk)

         write(nfout,'(" !ldos: meshwk_layer(*,2) first ",i8," elements")') n
         write(nfout,'(14i5)') (meshwk_layer(i,2),i=1,n)
         write(nfout,'(" !ldos: meshwk_layer(*,2) last  ",i8," elements")') n
         write(nfout,'(14i6)') (meshwk_layer(i,2),i=nfftwk-n+1, nfftwk)
#endif
      end if
      deallocate(mesh_t)
    end subroutine fillup_meshwk_layer

  end subroutine m_Ldos_cal

  subroutine m_Ldos_cal_partial()
  end subroutine m_Ldos_cal_partial

! ======================= added by K. Tagami =================== 11.0
  subroutine m_Ldos_cal_noncl( nfldos, print_flag )
    use m_IterationNumbers,   only : nk_in_the_process
    use m_Ionic_System,       only : speciesname

    integer, intent(in) :: nfldos
    logical, intent(in) :: print_flag

    integer,allocatable,dimension(:,:) :: meshwk
    integer :: is, ik, ip, ib
    integer :: ekmode_t
    integer :: istmp
    integer :: sw_atom_decomp, sw_layer_decomp

    real(kind=DP) :: chgq0( ndim_magmom )
    real(kind=DP), allocatable, dimension(:,:) :: zrhoik_kt
    real(kind=DP), allocatable, dimension(:,:) :: bfft_kt
    real(kind=DP), allocatable, dimension(:,:) :: chgq_enl_kt(:,:,:)

    real(kind=DP), allocatable, dimension(:,:,:,:) :: wsctmp,laytmp

    sw_atom_decomp = OFF;  sw_layer_decomp = OFF
    if ( sw_aldos == ON .or. sw_calc_wf_atom_decomposition == ON ) then
       sw_atom_decomp = ON
    endif
    if ( sw_layerdos == ON .or. sw_calc_wf_layer_decomposition == ON ) then
       sw_layer_decomp = ON
    endif

    if ( mype == 0 ) call wd_header

    if (sw_layer_decomp == ON) then
       allocate( zrhoik_kt(kg*kimg,ndim_magmom) )
       zrhoik_kt = 0.0d0
    endif

    allocate( chgq_enl_kt(kgp,kimg,ndim_magmom) )
    chgq_enl_kt = 0.0d0

    call m_CD_keep_retrieve_hsr(.true.)

    if (ekmode == ON) then
       do ik = 1, kv3, ndim_spinor
          if ( nk_in_the_process +ik -1 > kv3_ek ) cycle

          call m_FFT_alloc_WF_work()

          if(ipridos >= 3) write(nfout,'(" !!ldos ik = ",i6)') ik
          if(sw_atom_decomp == ON) call fillup_meshwk(mode=SOFTPART)

          if(ipridos >= 3) write(nfout,'(" !!ldos nfft = ",i6," <<m_Ldos_cal>>")') nfft

          allocate(bfft_kt(nfft,ndim_spinor))

          do ib = 1, neg
             bfft_kt = 0.d0
             if (map_ek(ib,ik) == mype ) then
                Do is=1, ndim_spinor
                   if ( sw_band_unfolding == ON .and. band_unfolding_active ) then
                      call m_ES_WF_in_Rspace_kt( ik+is-1, ik+is-1, ik+is-1, &
                           &                     zaj_l(:,map_z(ib),ik+is-1,:), &
                           &                     bfft_kt(:,is) )
                   else
                      call m_ES_WF_in_Rspace(ik+is-1,ib,bfft_kt(:,is))
                   endif
                End do
!!                ip = mype
!!             end if
!!             if (npes > 1) then
!!                call mpi_bcast( bfft_kt, nfft*ndim_spinor, mpi_double_precision,&
!!                     &          ip, MPI_CommGroup, ierr )
!!             endif

!!             if(sw_atom_decomp == ON) call substitute_weiwsc_noncl(1,ib,.false.)
!!             if(sw_layer_decomp == ON) call substitute_weilay_noncl(1,ib,.false.)
                if(sw_atom_decomp == ON) call substitute_weiwsc_noncl(ik,ib,.false.)
                if(sw_layer_decomp == ON) call substitute_weilay_noncl(ik,ib,.false.)
             endif
          end do
          if ( map_k(ik) == 0 ) call check_sum_noncl(ik,"softpart",.false.)

          deallocate(bfft_kt)
          call m_FFT_dealloc_WF_work()

          if (modnrm == EXECUT) then
             call m_FFT_alloc_CD_box()
             allocate( bfft_kt(nfftp_nonpara,ndim_magmom) )

             if (sw_atom_decomp==ON) call fillup_meshwk(mode=HARDPART)
             do ib = 1, neg
                call m_CD_hardpart_sub_noncl( nfout, ik, ib, chgq0 )  ! -> chgq_l
                Do istmp=1, ndim_magmom
                   call m_CD_map_chgq_to_fft_box( istmp, nfftp_nonpara, &
                        &                         bfft_kt(:,istmp) )   ! -> bfft

                   call m_FFT_CD_inverse0( nfout, bfft_kt(:,istmp) )
                                                   ! bfft(G_sp.) -> bfft(R_sp.)
                End do

                if (sw_atom_decomp==ON) call substitute_weiwsc_cd_noncl(ik,ib,.false.)
                if (sw_layer_decomp==ON) then
                   call substitute_weilay_cd_noncl(ik,ib,chgq0,.false., chgq_enl_kt)
                endif
             end do

             if (sw_atom_decomp==ON) deallocate(meshwk)
             deallocate(bfft_kt)
             call m_FFT_dealloc_CD_box()

             if ( map_k(ik) == 0 ) call check_sum_noncl(ik,"total   ",.false.)
          end if

          if (sw_save_ldos_weight == ON .and. sw_lband == OFF ) then
             if (sw_atom_decomp==ON) call wd_weight_noncl(ALDOS,nfldos)
             if (sw_layer_decomp==ON) call wd_weight_noncl(LAYERDOS,nfldos)
          end if

       end do

       if ( sw_lband == ON ) then
          if(sw_atom_decomp==ON) call wd_weight2_noncl(ALDOS,nfldos)
          if(sw_layer_decomp==ON) call wd_weight2_noncl(LAYERDOS,nfldos)
       endif

    else if (ekmode == OFF) then
       if (modnrm == EXECUT) call m_CD_set_ylm_enl_etc()
!       if (sw_save_ldos_weight == ON) call m_Ldos_wd_natm2_and_totch( nfldos )

       if (hardpart_subroutine==2 .or. sw_rspace_ldos==ON &
            &                     .or. sw_rspace_lband==ON )then

          do ik = 1, kv3, ndim_spinor
             if (map_k(ik) /= myrank_k) cycle

             if (ipridos >= 2) write(nfout,'(" !!ldos ik = ",i6)') ik
             call m_FFT_alloc_WF_work()
             if (sw_atom_decomp==ON) call fillup_meshwk(mode=SOFTPART)

             allocate(bfft_kt(nfft,ndim_spinor))

             if (ipridos >= 2 .and. printable) then
                write(nfout,'(" !!ldos nfft = ",i6," <<m_Ldos_cal>>")') nfft
             endif
             if (ipridos >= 2 .and. printable) then
                call m_ES_wd_zaj_small_portion(nfout,ik, &
                     &    " -- before m_ES_WF_in_Rspace <<m_Ldos_cal>> --",46)
             endif

             do ib = ista_e, iend_e, istep_e
                bfft_kt = 0.d0
                Do is=1, ndim_spinor
                   if ( sw_band_unfolding == ON .and. band_unfolding_active ) then
                      call m_ES_WF_in_Rspace_kt( ik+is-1, ik+is-1, ik+is-1, &
                           &                     zaj_l(:,map_z(ib),ik+is-1,:), &
                           &                     bfft_kt(:,is) )
                   else
                      call m_ES_WF_in_Rspace(ik+is-1,ib,bfft_kt(:,is))
                   endif

                End do

                if (sw_atom_decomp==ON) call substitute_weiwsc_noncl(ik,ib,.false.)
                if (sw_layer_decomp==ON) call substitute_weilay_noncl(ik,ib,.false.)

             end do

             call check_sum_noncl(ik,"softpart",.false.)

             deallocate(bfft_kt)
             call m_FFT_dealloc_WF_work()

             if (modnrm == EXECUT) then
                call m_FFT_alloc_CD_box()

                allocate( bfft_kt(nfftp_nonpara,ndim_spinor) )

                if (sw_atom_decomp==ON) call fillup_meshwk(mode=HARDPART)
                do ib = ista_e, iend_e, istep_e
                   if (ipridos>=2) write(nfout,'(" !ldos   ib = ", i5)') ib

                   call m_CD_hardpart_sub2_noncl( nfout, ik, ib, chgq0, &
                        &                         chgq_enl_kt )     ! -> chgq_enl

                   Do istmp=1, ndim_magmom
                      call m_CD_map_chgqenl_to_fft_box_kt( istmp, nfftp_nonpara, &
                           &                               bfft_kt(:,istmp), &
                           &                               chgq_enl_kt )
                                                                 ! -> bfft
                      call m_FFT_CD_inverse0( nfout, bfft_kt(:,istmp) )
                                                    ! bfft(G_sp.) -> bfft(R_sp.)
                   End do

                   if (sw_atom_decomp==ON) call substitute_weiwsc_cd_noncl(ik,ib,.false.)
                   if (sw_layer_decomp==ON) then
                      call substitute_weilay_cd_noncl( ik, ib, chgq0, .false., &
                           &                          chgq_enl_kt )
                   endif

                end do

                if (sw_atom_decomp==ON) deallocate(meshwk)
                deallocate(bfft_kt)
                call m_FFT_dealloc_CD_box()

                call check_sum_noncl(ik,"total   ",.false.)
             end if

          end do

       else
          if ( sw_atom_decomp==ON ) weiwsc_noncl = 0.d0
          if ( sw_layer_decomp==ON ) weilay_noncl = 0.d0

          do ik = 1, kv3, ndim_spinor
             if (map_k(ik) /= myrank_k) cycle
             if (ipridos >= 2) write(nfout,'(" !!ldos ik = ",i6)') ik
             call m_FFT_alloc_WF_work()

             if (sw_atom_decomp==ON) call fillup_meshwk(mode=SOFTPART)

             allocate( bfft_kt(nfft,ndim_spinor) ); bfft_kt = 0.d0

             if (ipridos >= 2 .and. printable) then
                write(nfout,'(" !!ldos nfft = ",i6," <<m_Ldos_cal>>")') nfft
             endif
             if (ipridos >= 2 .and. printable) then
                call m_ES_wd_zaj_small_portion(nfout,ik, &
                     &  " -- before m_ES_WF_in_Rspace <<m_Ldos_cal>> --",46)
             endif

             do ib = 1,neg
                if(map_e(ib)/=myrank_e) cycle
                bfft_kt = 0.d0
                Do is=1, ndim_spinor
                   if ( sw_band_unfolding == ON .and. band_unfolding_active ) then
                      call m_ES_WF_in_Rspace_kt( ik+is-1, ik+is-1, ik+is-1, &
                           &                     zaj_l(:,map_z(ib),ik+is-1,:), &
                           &                     bfft_kt(:,is) )
                   else
                      call m_ES_WF_in_Rspace(ik+is-1,ib,bfft_kt(:,is))
                   endif
                End do

                if ( sw_lband == ON ) then
                   if (sw_atom_decomp==ON) call substitute_weiwsc_noncl(ik,ib,.false.)
                   if (sw_layer_decomp==ON) call substitute_weilay_noncl(ik,ib,.false.)
                else
                   if (sw_atom_decomp==ON) call substitute_weiwsc_noncl(ik,ib,.true.)
                   if (sw_layer_decomp==ON) call substitute_weilay_noncl(ik,ib,.true.)
                endif
             end do
             deallocate(bfft_kt)
             call m_FFT_dealloc_WF_work()
          enddo

          if ( sw_lband == OFF .and. npes > 1 ) then
             if (sw_atom_decomp==ON)then
                allocate( wsctmp(maldos,neg,kv3,ndim_magmom) );wsctmp = 0.d0
                call mpi_allreduce( weiwsc_noncl, wsctmp, kv3*maldos*neg*ndim_magmom, &
                     &              mpi_double_precision, &
                     &              mpi_sum,MPI_CommGroup,ierr )
                weiwsc_noncl = wsctmp
                deallocate(wsctmp)
             endif
             if (sw_layer_decomp==ON)then
                allocate( laytmp(mlayer,neg,kv3,ndim_magmom) );laytmp=0.d0
                call mpi_allreduce( weilay_noncl, laytmp, kv3*mlayer*neg*ndim_magmom, &
                     &              mpi_double_precision, &
                     &              mpi_sum,MPI_CommGroup, ierr )
                weilay_noncl = laytmp
                deallocate(laytmp)
             endif
          endif

          do ik=1, kv3, ndim_spinor
             if ( sw_lband == ON ) then
                if ( map_k(ik) == 0 ) call check_sum_noncl(ik,"softpart",.false.)
             else
                call check_sum_noncl(ik,"softpart",.true.)
             endif
          enddo

          if (modnrm == EXECUT) then
             allocate( bfft_kt(nfftp_nonpara,ndim_magmom) )
             do ik = 1, kv3, ndim_spinor
                call m_FFT_alloc_CD_box()
                if (sw_atom_decomp==ON) call fillup_meshwk(mode=HARDPART)
                do ib = 1,neg
                   if (ipridos>=2) write(nfout,'(" !ldos   ib = ", i5)') ib
                   call m_CD_hardpart_sub_noncl( nfout, ik, ib, chgq0 )
                                                ! -> chgq_enl

                   Do istmp=1, ndim_magmom
                      call m_CD_map_valence_charge_to_fft_box( istmp, nfftp_nonpara,&
                           &                                   bfft_kt(:,istmp) )
                      call m_FFT_CD_inverse0( nfout, bfft_kt(:,istmp) )
                                                    ! bfft(G_sp.) -> bfft(R_sp.)
                   End do

                   if ( sw_lband == ON ) then
                      if (sw_atom_decomp==ON) &
                           &   call substitute_weiwsc_cd_noncl(ik,ib,.false.)
                      if (sw_layer_decomp==ON) then
                         call substitute_weilay_cd_noncl( ik, ib, chgq0,.false., &
                              &                           chgq_enl_kt )
                      endif
                   else
                      if (sw_atom_decomp==ON) &
                           &   call substitute_weiwsc_cd_noncl(ik,ib,.true.)
                      if (sw_layer_decomp==ON) then
                         call substitute_weilay_cd_noncl( ik, ib, chgq0,.true., &
                              &                           chgq_enl_kt )
                      endif
                   endif
                end do

                if (sw_atom_decomp==ON) deallocate(meshwk)
                call m_FFT_dealloc_CD_box()

                if ( sw_lband == ON ) then
                   if ( map_k(ik) == 0 ) call check_sum_noncl(ik,"total   ",.false.)
                else
                   call check_sum_noncl(ik,"total   ",.true.)
                endif
             enddo
             deallocate(bfft_kt)
          endif

       endif

       if (npes > 1) call mpi_barrier(MPI_CommGroup,ierr)

       if (sw_save_ldos_weight == ON .and. print_flag ) then
          if ( sw_lband == ON ) then
             if(sw_atom_decomp==ON) call wd_weight2_noncl(ALDOS,nfldos)
             if(sw_layer_decomp==ON) call wd_weight2_noncl(LAYERDOS,nfldos)
          else
             if (hardpart_subroutine==2 .or. sw_rspace_ldos==ON &
                  &                     .or. sw_rspace_lband==ON )then
                if (sw_atom_decomp==ON) call wd_weight_noncl(ALDOS,nfldos)
                if (sw_layer_decomp==ON) call wd_weight_noncl(LAYERDOS,nfldos)
             else
                if(sw_atom_decomp==ON) call wd_weight_serial_noncl(ALDOS,nfldos)
                if(sw_layer_decomp==ON) call wd_weight_serial_noncl(LAYERDOS,nfldos)
             endif
          end if
       endif

    end if

    call m_CD_keep_retrieve_hsr(.false.)

    if (modnrm == EXECUT) call m_CD_dealloc_ylm_enl_etc()

    if (hardpart_subroutine/=2 .and. &
         &  (sw_rspace_ldos==OFF .or. sw_rspace_lband==OFF) ) call m_CD_restore_chgq()

    if (allocated(meshwk)) deallocate(meshwk)
    if (sw_layer_decomp==ON) deallocate(zrhoik_kt)
    deallocate( chgq_enl_kt )

  contains

    subroutine wd_header
      integer :: i, j

      if ( sw_lband == OFF ) then
         if(sw_save_ldos_weight == ON .and. ekmode==OFF) then
            call m_Ldos_wd_natm2_and_totch( nfldos )
         endif
         if ( ekmode == ON .and. nk_in_the_process == 1 ) then
            call m_Ldos_wd_natm2_and_totch( nfldos )
         endif
      endif
      if ( sw_lband == ON ) then
         if ( ( ekmode == ON .and. nk_in_the_process == 1 ) &
              &  .or. ( ekmode == OFF) ) then
            write(nfldos,'(A)') "# Local Decomposition for bands"
            write(nfldos,*) ""
            write(nfldos,'(A,I8)') ' num_kpoints = ', max(kv3,kv3_ek) /ndim_spinor
            write(nfldos,'(A,I8)') ' num_bands   = ', neg -num_extra_bands
            write(nfldos,'(A,I8)') ' nspin       = ', nspin /ndim_spinor
            write(nfldos,'(A,I8)') ' ndim_magmom = ', ndim_magmom
            write(nfldos,*) ""
            write(nfldos,'(A)') "# Decomposition Info."

            if ( sw_atom_decomp == ON ) then
               write(nfldos,'(A,I8,A)') ' num_atoms   = ', natom_decomp +1, &
                    &        "  ( last 1 corresponds to external region )"
            endif
            if ( sw_layer_decomp == ON ) then
               write(nfldos,'(A,I8,A)') ' num_layers  = ', nlayer +1, &
                    &        "  ( last 1 corresponds to external region )"
            endif
            write(nfldos,*) ""

            if ( sw_atom_decomp == ON ) then
               write(nfldos,'(A)') "# Atom Info."
               write(nfldos,'(A)') "   no.    target_atom    species"
               Do i=1, natom_decomp
                  j = atom_decomp_map(i)
                  write(nfldos,'(I5,5X,I5,12X,A)') i, j, speciesname(ityp(j))
               end Do
               write(nfldos,'(I5)') natom_decomp +1
               write(nfldos,*)
            endif
            if ( sw_layer_decomp == ON ) then
               write(nfldos,'(A)') "# Layer Info."
               write(nfldos,'(A)') "   no.       min.         max.  (Angstrom)"
               Do i=1, nlayer +1
                  write(nfldos,'(I5,2(5X,F8.4))') i, winlay(i,1)*BOHR, winlay(i,2)*BOHR
               end Do
               write(nfldos,*)
            endif
         endif
      endif
    end subroutine wd_header

    subroutine wd_weight2_noncl(mode,nfldos)
      integer, intent(in) :: mode, nfldos
      integer :: m, ib, ib_ordr, kend, kstep, iws, ik, iks
      integer :: iktmp, istmp
      integer :: count
      real(kind=DP),allocatable,dimension(:,:,:,:) :: wei

      if(mode==ALDOS) then
         m = maldos
      else if(mode==LAYERDOS) then
         m = mlayer
      end if
      allocate(wei(m,neg,kv3/ndim_spinor,ndim_magmom)); wei = 0.d0

      if(mode==ALDOS)then
         do ik=ista_k, iend_k, ndim_spinor
            iktmp = (ik-1)/ndim_spinor +1
            do ib = ista_e, iend_e, istep_e
               wei(1:m,ib,iktmp,1:ndim_magmom) &
!                    &   = weiwsc_noncl(1:m,ib,ik,1:ndim_magmom)
                    &   = weiwsc_noncl(1:m,map_z(ib),ik,1:ndim_magmom)
            end do
         end do
      else if(mode==LAYERDOS)then
         do ik=ista_k, iend_k, ndim_spinor
            iktmp = (ik-1)/ndim_spinor +1
            do ib = ista_e, iend_e, istep_e
               wei(1:m,ib,iktmp,1:ndim_magmom) &
                    &   = weilay_noncl(1:m,map_z(ib),ik,1:ndim_magmom)
            end do
         end do
      endif

      if (npes>1) then
         call mpi_allreduce( MPI_IN_PLACE, wei, m*neg*kv3/ndim_spinor*ndim_magmom, &
              &              mpi_double_precision, &
              &              mpi_sum, MPI_CommGroup, ierr )
      end if

      if(mype == 0) then
         Do ik=1, kv3, ndim_spinor
            iks = 0
            if ( sw_lband == ON .and. ekmode == ON ) then
               iks = nk_in_the_process -1
               if ( ik +iks > kv3_ek ) cycle
            endif
            if(mode==ALDOS) then
               write(nfldos,'(" weight for each atomic cell    nk = ",i8)') ik+iks
            else if(mode==LAYERDOS) then
               write(nfldos,'(" weight for each layer          nk = ",i8)') ik+iks
            end if

            iktmp = ( ik-1 )/ndim_spinor +1
            Do istmp=1, ndim_magmom
               select case( istmp )
               case (1)
                  write(nfldos,*) '! -------- for total charge ----- '
               case (2)
                  write(nfldos,*) '! -------- for magnetic moment X ----- '
               case (3)
                  write(nfldos,*) '! -------- for magnetic moment Y ----- '
               case (4)
                  write(nfldos,*) '! -------- for magnetic moment Z ----- '
               end select

               if ( mode==ALDOS .and. sw_wd_only_specified_atoms == ON ) then
                  do ib = 1, neg -num_extra_bands
                     write(nfldos,'(i4,A)',advance="no") ib, ")"
                     count = 0
                     do iws=1, m-1
                        if ( if_aldos_full(iws)==1 ) then
                           count = count +1
                           if ( count > 1 .and. mod(count,ilen)==1 ) then
                              write(nfldos,'(5x)',advance="no")
                           endif
                           write(nfldos,'(f15.10)',advance="no") wei(iws,ib,iktmp,istmp)
                           if ( mod(count,ilen)==0 ) write(nfldos,*)
                        endif
                     end do
                     count = count +1
                     if ( count > 1 .and. mod(count,ilen)==1 ) then
                        write(nfldos,'(5x)',advance="no")
                     endif
                     write(nfldos,'(f15.10)') wei(maldos,ib,iktmp,istmp)
                  end do
               else
                  do ib = 1, neg -num_extra_bands
                     if(ilen >= m) then
                        write(nfldos,'(i4,")",5f15.10)') &
                             &         ib,(wei(iws,ib,iktmp,istmp),iws=1,m)
                     else
                        write(nfldos,'(i4,")",5f15.10)') &
                             &         ib,(wei(iws,ib,iktmp,istmp),iws=1,ilen)
                        write(nfldos,'(5x,5f15.10)') &
                             &         (wei(iws,ib,iktmp,istmp),iws=ilen+1,m)
                     end if
                  end do
               endif
            end Do
         end Do
      end if

      deallocate(wei)
    end subroutine wd_weight2_noncl

    subroutine wd_weight_noncl(mode,nfldos)
      integer, intent(in) :: mode, nfldos
      integer :: m, ik, ib, ib_ordr, kend, kstep, iws, ip, count
      real(kind=DP),allocatable,dimension(:) :: wei

      integer :: istmp

      if(mode==ALDOS) then
         m = maldos
      else if(mode==LAYERDOS) then
         m = mlayer
      end if

      allocate(wei(m))

      if (ekmode==OFF) then
         kend = kv3;
      else if (ekmode==ON) then
!         kend = 1;
         kend = kv3;
      end if
      kstep = ndim_spinor

      if (ipridos >= 3) then
         do ik = 1, kv3, kstep
            if(map_k(ik) /= myrank_k) cycle
            write(nfout,'(" !ldos: ik = ",i8)') ik
            do ib = 1, neg
               write(nfout,'(" !ldos: ib, ib_ordr = ",2i20)') ib, neordr(ib,ik)
            end do
         end do
      end if

      do ik = 1, kend, kstep
         if (mype == 0) then
            if (mode==ALDOS) then
               write(nfldos,'(" weight for each atomic cell    nk = ",i8)') ik
            else if (mode==LAYERDOS) then
               write(nfldos,'(" weight for each layer          nk = ",i8)') ik
            end if
         end if

         Do istmp=1, ndim_magmom

            if ( mype == 0 ) then
               select case( istmp )
               case (1)
                  write(nfldos,*) '! -------- for total charge ----- '
               case (2)
                  write(nfldos,*) '! -------- for magnetic moment X ----- '
               case (3)
                  write(nfldos,*) '! -------- for magnetic moment Y ----- '
               case (4)
                  write(nfldos,*) '! -------- for magnetic moment Z ----- '
               end select
            endif

            do ib = 1, neg -num_extra_bands

               if (nrank_k >= 2) then
                  if (map_ek(ib,ik) == mype) ib_ordr = neordr(ib,ik)
                  call mpi_bcast( ib_ordr, 1, mpi_integer, map_ek(ib,ik), &
                       &          MPI_CommGroup, ierr )
               else
                  ib_ordr = neordr(ib,ik)
               end if

               if (map_ek(ib_ordr,ik) == mype) then
                  if (mode == ALDOS) then
                     wei(1:m) = weiwsc_noncl(1:m,map_z(ib_ordr),ik,istmp)
                  else if (mode == LAYERDOS) then
                     wei(1:m) = weilay_noncl(1:m,map_z(ib_ordr),ik,istmp)
                  end if
                  if (map_ek(ib_ordr,ik) /= 0) &
                       call mpi_send( wei, m, mpi_double_precision, &
                       &              0, 1, MPI_CommGroup, ierr )
               else if (mype == 0 .and. map_ek(ib_ordr,ik) /= 0) then
                  call mpi_recv( wei, m, mpi_double_precision, &
                       &         map_ek(ib_ordr,ik), 1, MPI_CommGroup, istatus, ierr )
               end if

               if ( mode==ALDOS .and. sw_wd_only_specified_atoms == ON ) then
                  write(nfldos,'(i4,A)',advance="no") ib, ")"
                  count = 0
                  do iws=1, m-1
                     if ( if_aldos_full(iws)==1 ) then
                        count = count +1
                        if ( count > 1 .and. mod(count,ilen)==1 ) then
                           write(nfldos,'(5x)',advance="no")
                        endif
                        write(nfldos,'(f15.10)',advance="no") wei(iws)
                        if ( mod(count,ilen)==0 ) write(nfldos,*)
                     endif
                  end do
                  count = count +1
                  if ( count > 1 .and. mod(count,ilen)==1 ) then
                     write(nfldos,'(5x)',advance="no")
                  endif
                  write(nfldos,'(f15.10)') wei(maldos)
               else
                  if(ilen >= m) then
                     write(nfldos,'(i4,")",5f15.10)') ib,(wei(iws),iws=1,m)
                  else
                     write(nfldos,'(i4,")",5f15.10)') ib,(wei(iws),iws=1,ilen)
                     write(nfldos,'(5x,5f15.10)') (wei(iws),iws=ilen+1,m)
                  end if
               endif

            end do
         End do
      end do

      deallocate(wei)
    end subroutine wd_weight_noncl

    subroutine wd_weight_serial_noncl(mode,nfldos)
      integer, intent(in) :: mode, nfldos
      integer :: m, ik, ib, ib_ordr, kend, kstep, iws, ip, count
      real(kind=DP),allocatable,dimension(:):: wei

      integer :: istmp

      if (mode==ALDOS) then
         m = maldos
      else if (mode==LAYERDOS) then
         m = mlayer
      end if

      allocate(wei(m))

      if(ekmode==OFF) then
         kend = kv3
      else if(ekmode==ON) then
!         kend = 1
         kend = kv3
      end if
      kstep = ndim_spinor

      if (ipridos >= 3) then
         do ik = 1, kv3, kstep
            if(map_k(ik) /= myrank_k) cycle
            write(nfout,'(" !ldos: ik = ",i8)') ik
            do ib = 1, neg
               write(nfout,'(" !ldos: ib, ib_ordr = ",2i20)') ib, neordr(ib,ik)
            end do
         end do
      end if

      do ik = 1, kend, kstep
         if (mype == 0) then
            if (mode==ALDOS) then
               write(nfldos,'(" weight for each atomic cell    nk = ",i8)') ik
            else if (mode==LAYERDOS) then
               write(nfldos,'(" weight for each layer          nk = ",i8)') ik
            end if
         end if

         Do istmp=1, ndim_magmom

            if ( mype == 0 ) then
               select case( istmp )
               case (1)
                  write(nfldos,*) '! -------- for total charge ----- '
               case (2)
                  write(nfldos,*) '! -------- for magnetic moment X ----- '
               case (3)
                  write(nfldos,*) '! -------- for magnetic moment Y ----- '
               case (4)
                  write(nfldos,*) '! -------- for magnetic moment Z ----- '
               end select
            endif

            do ib = 1, neg -num_extra_bands

               if (mode == ALDOS) then
                  wei(1:m) = weiwsc_noncl(1:m,ib,ik,istmp)
               else if(mode == LAYERDOS) then
                  wei(1:m) = weilay_noncl(1:m,ib,ik,istmp)
               end if

               if ( mode==ALDOS .and. sw_wd_only_specified_atoms == ON ) then
                  write(nfldos,'(i4,A)',advance="no") ib, ")"
                  count = 0
                  do iws=1, m-1
                     if ( if_aldos_full(iws)==1 ) then
                        count = count +1
                        if ( count > 1 .and. mod(count,ilen)==1 ) then
                           write(nfldos,'(5x)',advance="no")
                        endif
                        write(nfldos,'(f15.10)',advance="no") wei(iws)
                        if ( mod(count,ilen)==0 ) write(nfldos,*)
                     endif
                  end do
                  count = count +1
                  if ( count > 1 .and. mod(count,ilen)==1 ) then
                     write(nfldos,'(5x)',advance="no")
                  endif
                  write(nfldos,'(f15.10)') wei(maldos)
               else
                  if(ilen >= m) then
                     write(nfldos,'(i4,")",5f15.10)') ib,(wei(iws),iws=1,m)
                  else
                     write(nfldos,'(i4,")",5f15.10)') ib,(wei(iws),iws=1,ilen)
                     write(nfldos,'(5x,5f15.10)') (wei(iws),iws=ilen+1,m)
                  end if
               endif
            end do
         End Do
      end do

      deallocate(wei)
    end subroutine wd_weight_serial_noncl

    subroutine check_sum_noncl(ik,aword,serial)
      integer,intent(in) :: ik
      character*(*),intent(in) :: aword
      logical, intent(in) :: serial

      integer :: ib, ia, ilay, ib_ordr, ik_t
      real(kind=DP), allocatable, dimension(:,:) :: sum_aldos, sum_layer, sum_mpi

      if (sw_atom_decomp==ON) allocate(sum_aldos(neg,ndim_magmom))
      if (sw_layer_decomp==ON) allocate(sum_layer(neg,ndim_magmom))

      if ( sw_lband == ON ) then
         ik_t = ik
      else
         if (ekmode /= ON) then
            ik_t = ik
         else if (ekmode == ON) then
            ik_t = 1
         end if
      endif

      if (sw_atom_decomp==ON .or. sw_layer_decomp==ON) then
         if (nrank_e >= 2) then
            allocate(sum_mpi(neg,ndim_magmom)); sum_mpi = 0.d0
         end if
      end if

      if (sw_atom_decomp==ON) then
         sum_aldos = 0.d0
         do ib = 1, neg
            if (.not.serial)then
               if (map_ek(ib,ik_t) == mype) then
                  do ia = 1, natm2 + 1
                     sum_aldos(ib,:) = sum_aldos(ib,:) &
                          &          + weiwsc_noncl(ia,map_z(ib),ik_t,:)
                  end do
               endif
            else
               do ia = 1, natm2 + 1
                  sum_aldos(ib,:) = sum_aldos(ib,:) &
                       &          + weiwsc_noncl(ia,ib,ik_t,:)
               end do
            end if
         end do
         if (.not.serial)then
            if (nrank_e >= 2) then
               sum_mpi = 0.d0
               call mpi_allreduce( sum_aldos, sum_mpi, neg*ndim_magmom, &
                    &              mpi_double_precision, mpi_sum, &
                    &              mpi_k_world(myrank_k), ierr )
               sum_aldos = sum_mpi
            end if
         endif
      end if

      if (sw_layer_decomp==ON) then
         sum_layer = 0.d0
         do ib = 1, neg
            if (.not.serial)then
               if (map_e(ib) == myrank_e) then
                  do ilay = 1, mlayer
                     sum_layer(ib,:) = sum_layer(ib,:) &
                          &          + weilay_noncl(ilay,map_z(ib),ik_t,:)
                  end do
               endif
            else
               do ilay = 1, mlayer
                  sum_layer(ib,:) = sum_layer(ib,:) &
                       &          + weilay_noncl(ilay,ib,ik_t,:)
               end do
            end if
         end do
         if (.not.serial)then
            if (nrank_e >= 2) then
               sum_mpi = 0.d0
               call mpi_allreduce( sum_layer, sum_mpi, neg*ndim_magmom, &
                    &              mpi_double_precision, mpi_sum, &
                    &              mpi_k_world(myrank_k), ierr )
               sum_layer = sum_mpi
            end if
         endif
      end if

!!$      if(mype == 0) then
      if (ipridos>=1) then
         if (sw_atom_decomp==ON .and. sw_layer_decomp==ON) then
            write(nfout,'(" !!ldos ",a8," -- iban , sum(atomic, layer) --  neg = ",i6 &
                 & ," natm2 = ",i6, " ik = ",i6)') aword,neg,natm2,ik
            do ib = 1, neg
               if ( sw_lband == ON ) then
                  ib_ordr = ib
               else
                  if(serial) then
                     ib_ordr = ib
                  else
                     ib_ordr = neordr(ib,ik)
                  endif
               endif
               write(nfout,'(" !!ldos ", i4,2f16.8)') ib, sum_aldos(ib_ordr,1),&
                    &                                     sum_layer(ib_ordr,1)
            end do

         else if( sw_atom_decomp==ON ) then
            write(nfout,'(" !!ldos ",a8," -- iban , sum(atomic) --  neg = ",i6 &
                 & ," natm2 = ",i6," ik = ",i6)') aword,neg,natm2,ik
            do ib = 1, neg
               if ( sw_lband == ON ) then
                  ib_ordr = ib
               else
                  if(serial) then
                     ib_ordr = ib
                  else
                     ib_ordr = neordr(ib,ik)
                  endif
               endif
               write(nfout,'(" !!ldos ", i4,f16.8)') ib, sum_aldos(ib_ordr,1)
            end do

         else if ( sw_layer_decomp==ON ) then
            write(nfout,'(" !!ldos ",a8," -- iban , sum(layer ) --  neg = ",i6 &
                 & ," natm2 = ",i6," ik = ",i6)') aword,neg,natm2,ik
            do ib = 1, neg
               if ( sw_lband == ON ) then
                  ib_ordr = ib
               else
                  if(serial) then
                     ib_ordr = ib
                  else
                     ib_ordr = neordr(ib,ik)
                  endif
               endif
               write(nfout,'(" !!ldos ",i4,f16.8)') ib, sum_layer(ib_ordr,1)
            end do
         end if
      end if

      if (sw_atom_decomp==ON) deallocate(sum_aldos)
      if (sw_layer_decomp==ON) deallocate(sum_layer)
      if (allocated(sum_mpi)) deallocate(sum_mpi)
!!$      if(ekmode /= ON .and. (sw_aldos == ON .or. sw_layerdos == ON) .and. npes >= 2) deallocate(sum_mpi)

    end subroutine check_sum_noncl

    subroutine substitute_weilay_noncl(ik,ib,serial)
!   *********************************
!     by T.Yamasaki
!           10th Jun 1992
!   *********************************
!     Revised by T. Yamasaki, Feb. 2004
!
      integer, intent(in) :: ik,ib
      logical, intent(in) :: serial

      integer :: i, ilay, iax1, iax2, iax3, ib_t
      real(kind=DP) :: a1,a2,width, g,ss

      complex(kind=CMPLDP) :: zi, zsum, zchg
      real(kind=DP),allocatable,dimension(:,:) :: zsum_layer, zsum_layer_mpi

      integer :: is1, is2, istmp
      real(kind=DP) :: rho_magmom(ndim_magmom)
      real(kind=DP) :: c_nr(ndim_spinor,ndim_spinor)
      real(kind=DP) :: c_ni(ndim_spinor,ndim_spinor)
      real(kind=DP) :: cr, ci, denom
      real(kind=DP), allocatable :: afft_kt(:,:)

      denom = 1.d0/product(fft_box_size_WF(1:3,1))

      allocate( zsum_layer(mlayer,ndim_magmom) ); zsum_layer = 0.d0

!!$      zi = cmplx(0.d0,1.d0)
      zi = cmplx(0.d0,1.d0)*PAI2

      iax1 = normal_axis_winlay
      iax2 = iax1 + 1
      if(iax2 >= 4) iax2 = iax2 - 3
      iax3 = iax2 + 1
      if(iax3 >= 4) iax3 = iax3 - 3


!!$      rl11 = rltv(iax1,iax1)
!!$      ss = 1.d0/maxhv
      ss = 3.d0 - kimg

!!$      call get_height(height)

!!$      n = fft_box_size_WF(1,0)*fft_box_size_WF(2,0)*fft_box_size_WF(3,1)
!!$      n = nfft*kimg/kimg/2


      allocate( afft_kt(nfft,ndim_magmom) );  afft_kt = 0.0d0

      do i = 1, nfft-1, 2

         Do is1 = 1, ndim_spinor
            Do is2 = 1, ndim_spinor
               istmp = ( is1 -1 )*ndim_spinor + is2

               cr =    bfft_kt(i,  is1) *bfft_kt(i,  is2) &
                    & +bfft_kt(i+1,is1) *bfft_kt(i+1,is2)

               ci =   -bfft_kt(i,  is1) *bfft_kt(i+1,is2) &
                    & +bfft_kt(i+1,is1) *bfft_kt(i,  is2)

               c_nr(is1,is2) = cr;   c_ni(is1,is2) = ci
            End do
         End do
! --------------------------- convert charge from ss-rep to magmom-rep --
         rho_magmom(1) = c_nr(1,1) + c_nr(2,2);         ! ctot
         rho_magmom(2) = c_nr(1,2) + c_nr(2,1);         ! mx
         rho_magmom(3) = c_ni(2,1) - c_ni(1,2);         ! my
         rho_magmom(4) = c_nr(1,1) - c_nr(2,2);          ! mz
! -----------------------------------------------------------------------
         afft_kt(i,:) = rho_magmom(:)
      End do
!
      Do istmp=1, ndim_magmom
         call m_FFT_WF( ELECTRON, nfout, afft_kt(:,istmp), DIRECT, ON )
         if(kimg == 1) then
            do i = 1, kg
               zrhoik_kt(i,istmp) = afft_kt(igf(i),istmp)*denom
            end do
         else
            do i = 1, kg
               zrhoik_kt(i*2-1,istmp) = afft_kt(igf(i)*2-1,istmp)*denom
               zrhoik_kt(i*2,  istmp) = afft_kt(igf(i)*2,  istmp)*denom
            end do
         end if
      End do
! --
      do ilay = 1, mlayer
!!$         a2 = winlay(ilay,2)
!!$         a1 = winlay(ilay,1)
         a2 = winlay(ilay,2)/height
         a1 = winlay(ilay,1)/height
         width = a2 - a1

         Do istmp=1, ndim_magmom
            zsum = cmplx( zrhoik_kt(1,istmp)*width, 0.d0 )
!!$         zsum = univol*cmplx(zrhoik(1)*width,0.d0)

            if (kimg == 1) then
               do i = 2, kg
                  if(ngabc(i,iax2) == 0 .and. ngabc(i,iax3) == 0) then
!!$                  g = rl11*ngabc(i,iax1)
                     g = ngabc(i,iax1)
                     zsum = zsum + zrhoik_kt(i,istmp) &
                          &        *((cdexp(zi*g*a2)-cdexp(zi*g*a1))/( zi*g ))
                  end if
               end do
            else if(kimg == 2) then
               do i = 2, kg
                  if(ngabc(i,iax2) ==  0 .and. ngabc(i,iax3) == 0) then
!!$                  g = rl11*ngabc(i,iax1)
                     g = ngabc(i,iax1)
                     zchg = dcmplx( zrhoik_kt(i*2-1,istmp),zrhoik_kt(i*2,istmp) )
                     zsum = zsum  + zchg*((cdexp(zi*g*a2)-cdexp(zi*g*a1))/( zi*g ))
                  endif
               enddo
            endif
            zsum_layer(ilay,istmp) = real(zsum) *ss
         End do
      end do

      if ( sw_lband == ON ) then
         ib_t = map_z(ib)
      else
         if(ekmode == ON.or.serial) then
            ib_t = ib
         else
            ib_t = map_z(ib)
         end if
      endif

      do ilay = 1, mlayer
         weilay_noncl(ilay,ib_t,ik,:) = zsum_layer(ilay,:)
      end do

      deallocate(zsum_layer)
      deallocate(afft_kt)

    end subroutine substitute_weilay_noncl

    subroutine substitute_weiwsc_noncl(ik,ib,serial)
      integer, intent(in) :: ik,ib
      logical, intent(in) :: serial
!   *********************************
!     by T.Yamasaki
!           10th Jun 1992
!   *********************************
!     Revised by T. Yamasaki, Feb. 2004
!
      integer :: ijk, nwsc,n, nel, ib_t
      real(kind=DP) :: denom

      integer :: is1, is2
      real(kind=DP) :: rho_magmom(ndim_magmom)
      real(kind=DP) :: c_nr(ndim_spinor,ndim_spinor)
      real(kind=DP) :: c_ni(ndim_spinor,ndim_spinor)
      real(kind=DP) :: cr, ci

      nel = product(fft_box_size_WF(1:3,1))

      denom = 1.d0/dble(nel)
!      denom = 2.d0/dble(nel)

      if(ipridos >= 2) write(nfout,'(" !!ldos: ib, nfft, natm2 = ",i6,i12,i6,"  <<substitute_weiwsc_noncl>>")') &
           & ib, nfft, natm2

      if ( sw_lband == ON ) then
         ib_t = map_z(ib)
      else
         if(ekmode == ON.or.serial) then
            ib_t = ib
         else if(ekmode == OFF) then
            ib_t = map_z(ib)
         end if
      endif

      do ijk = 1, nel
         nwsc = meshwk(ijk,1)
         n    = meshwk(ijk,2)
         if((n < 1 .or. 2*n > nfft) .or. (nwsc < 1 .or. nwsc >natm2+1) ) then
            if(ipridos>=1) then
               write(nfout,'(" !!ldos: n = ",i6," ijk = ",i6)') n, ijk
               write(nfout,'(" !!ldos: nwsc = ",i9)') nwsc
            end if
         else

            Do is1 = 1, ndim_spinor
               Do is2 = 1, ndim_spinor
                  istmp = ( is1 -1 )*ndim_spinor + is2

                  cr =    bfft_kt(2*n-1,is1) *bfft_kt(2*n-1,is2) &
                       & +bfft_kt(2*n,  is1) *bfft_kt(2*n,  is2)

                  ci =   -bfft_kt(2*n-1,is1) *bfft_kt(2*n,  is2) &
                       & +bfft_kt(2*n,  is1) *bfft_kt(2*n-1,is2)

                  c_nr(is1,is2) = cr;   c_ni(is1,is2) = ci
               End do
            End do
! --------------------------- convert charge from ss-rep to magmom-rep --
            rho_magmom(1) = c_nr(1,1) + c_nr(2,2);         ! ctot
            rho_magmom(2) = c_nr(1,2) + c_nr(2,1);         ! mx
            rho_magmom(3) = c_ni(2,1) - c_ni(1,2);         ! my
            rho_magmom(4) = c_nr(1,1) - c_nr(2,2);          ! mz
! -----------------------------------------------------------------------

            weiwsc_noncl(nwsc,ib_t,ik,:) = weiwsc_noncl(nwsc,ib_t,ik,:) &
                 &                       + rho_magmom(:)*denom

         end if
      end do

      if(ipridos >= 2) then
         write(nfout,'(" !ldos: ik, ib = ",2i8," <<substitute_weiwsc_noncl>>")') ik, ib
         write(nfout,'(" !ldos: ",8f10.5)') (weiwsc_noncl(nwsc,ib_t,ik,1),nwsc=1,maldos)
      end if

    end subroutine substitute_weiwsc_noncl

    subroutine substitute_weiwsc_cd_noncl(ik,ib,serial)
      integer, intent(in) :: ik,ib
      logical, intent(in) :: serial
!   *********************************
!     by T.Yamasaki
!           10th Jun 1992
!   *********************************
!     Revised by T. Yamasaki, Feb. 2004
!
      integer :: ijk, nwsc,n, ib_t
      real(kind=DP) :: denom

!      denom = 0.5d0*univol/product(fft_box_size_CD(1:3,1))
      denom = univol/product(fft_box_size_CD(1:3,1))

      if ( sw_lband == ON ) then
         ib_t = map_z(ib)
      else
         if(ekmode == ON.or.serial) then
            ib_t = ib
         else if(ekmode == OFF) then
            ib_t = map_z(ib)
         end if
      endif

      do ijk = 1, product(fft_box_size_CD(1:3,1))
         nwsc = meshwk(ijk,1)
         n    = meshwk(ijk,2)
         if((n < 1 .or. 2*n > nfftp_nonpara) .or. (nwsc < 1 .or. nwsc >natm2+1) ) then
            if(ipridos>=1) then
               write(nfout,'(" !!ldos: n = ",i6," ijk = ",i6)') n, ijk
               write(nfout,'(" !!ldos: nwsc = ",i9)') nwsc
            end if
         else
            if ( map_ek(ib,ik) == mype) then
               weiwsc_noncl(nwsc,ib_t,ik,:) = weiwsc_noncl(nwsc,ib_t,ik,:) &
                    &                        + denom *bfft_kt(2*n-1,:)
            endif
         end if
      end do

      if(ipridos >=2) then
         write(nfout,'(" !ldos: ik ,ib = ",2i8)') ik, ib
         write(nfout,'(" !ldos: ",8f10.5)') (weiwsc_noncl(nwsc,ib_t,ik,1),nwsc=1,maldos)
      end if

    end subroutine substitute_weiwsc_cd_noncl

    subroutine substitute_weilay_cd_noncl( ik, ib, chgq0, serial, chgq_enl_kt )
      integer, intent(in) :: ik, ib
      real(kind=DP), intent(in) :: chgq0( ndim_magmom )
      real(kind=DP), intent(in) :: chgq_enl_kt( kgp,kimg,ndim_magmom )
      logical, intent(in) :: serial
!   *********************************
!     by T.Yamasaki
!           10th Jun 1992
!   *********************************
!     Revised by T. Yamasaki, Feb. 2004
!
      integer :: i, ista, iax1, iax2, iax3, ilay,  ib_t
      real(kind=DP) :: g, z_work, ss, a1, a2, width
      complex(kind=CMPLDP) :: zi,zsum,zchg
      real(kind=DP),allocatable,dimension(:,:,:) :: zsum_layer, zsum_layer_mpi

      allocate(zsum_layer(mlayer,2,ndim_magmom)); zsum_layer = 0.d0

      zi = cmplx(0.d0,1.d0)*PAI2

      iax1 = normal_axis_winlay
      iax2 = iax1 + 1
      if(iax2 >= 4) iax2 = iax2 - 3
      iax3 = iax2 + 1
      if(iax3 >= 4) iax3 = iax3 - 3

!!$      rl11 = rltv(iax1,iax1)
!!$      ss = 0.5d0*univol/maxhv

!!!!!!!!      ss = 0.5d0*univol*(3.d0-kimg)
      ss = univol*(3.d0-kimg)


      ista = ista_kngp
      if(ista == 1) ista = 2
      do ilay = 1, mlayer
         if(ipridos >= 2) write(nfout,'(" !!ldos: ib = ",i6," ilay = ",i6 &
              & ," <<substitute_weilay_cd_noncl>>")') ib,ilay
         a2 = winlay(ilay,2)/height
         a1 = winlay(ilay,1)/height
         width = a2 - a1

         Do istmp=1, ndim_magmom
            z_work = chgq0(istmp) *width
            if(serial)then
               if(ista/=2) z_work=0.d0
            endif
            zsum = 0.d0

            if(kimg == 1) then
               if ( sw_lband == ON ) then
                  do i = 2, kgp
                     if(ngabc(i,iax2) == 0 .and. ngabc(i,iax3) == 0) then
                        g = ngabc(i,iax1)
                        zsum = zsum + chgq_enl_kt(i,1,istmp) &
                             &     *( ( cdexp(zi*g*a2) - cdexp(zi*g*a1))/(zi*g))
                     end if
                  end do
               else
                  if(ekmode == ON.or.serial) then
                     do i = ista, iend_kngp
                        if(ngabc(i,iax2) == 0 .and. ngabc(i,iax3) == 0) then
                           g = ngabc(i,iax1)
                           zsum = zsum + chgq_l(i,1,istmp) &
                                &     *( ( cdexp(zi*g*a2) - cdexp(zi*g*a1))/(zi*g))

                        end if
                     end do
                  else if(ekmode == OFF) then
                     do i = 2, kgp
                        if(ngabc(i,iax2) == 0 .and. ngabc(i,iax3) == 0) then
                           g = ngabc(i,iax1)
                           zsum = zsum + chgq_enl_kt(i,1,istmp) &
                                &     *( ( cdexp(zi*g*a2) - cdexp(zi*g*a1))/(zi*g))
                        end if
                     end do
                  end if
               endif
            else if(kimg == 2) then
               if ( sw_lband == ON ) then
                  do i = 2, kgp
                     if(ngabc(i,iax2) == 0 .and. ngabc(i,iax3) == 0) then
                        g = ngabc(i,iax1)
                        zchg = cmplx(chgq_enl_kt(i,1,istmp),chgq_enl_kt(i,2,istmp))
                        zsum = zsum &
                             & + zchg * ((cdexp(zi*g*a2) - cdexp(zi*g*a1))/(zi*g))
                     end if
                  end do
               else
                  if(ekmode == ON.or.serial) then
                     do i = ista, iend_kngp
                        if(ngabc(i,iax2) == 0 .and. ngabc(i,iax3) == 0) then
                           g = ngabc(i,iax1)
                           zchg = cmplx(chgq_l(i,1,istmp),chgq_l(i,2,istmp))
                           zsum = zsum &
                                & + zchg * ((cdexp(zi*g*a2) - cdexp(zi*g*a1))/(zi*g))
                        end if
                     end do
                  else if(ekmode == OFF) then
                     do i = 2, kgp
                        if(ngabc(i,iax2) == 0 .and. ngabc(i,iax3) == 0) then
                           g = ngabc(i,iax1)
                           zchg = cmplx(chgq_enl_kt(i,1,istmp),chgq_enl_kt(i,2,istmp))
                           zsum = zsum &
                                & + zchg * ((cdexp(zi*g*a2) - cdexp(zi*g*a1))/(zi*g))
                        end if
                     end do
                  end if
               endif
            end if
            zsum = zsum + cmplx(z_work,0.d0)

            zsum_layer(ilay,1,istmp) = real(zsum)
            zsum_layer(ilay,2,istmp) = dimag(zsum)

         End do
      end do

#if 0
      if(ekmode == ON.or.serial) then
         if(npes > 1) then
            allocate( zsum_layer_mpi(mlayer,2,ndim_magmom) )
            call mpi_allreduce( zsum_layer, zsum_layer_mpi, mlayer*2*ndim_magmom, &
                 &              mpi_double_precision, mpi_sum, MPI_CommGroup,ierr )
            zsum_layer = zsum_layer_mpi
            deallocate(zsum_layer_mpi)
         end if
      end if

      if(mype == 0) then
         do ilay = 1, mlayer
            if(zsum_layer(ilay,2,1) > 1.d-6) then
               write(6,'(" !!ldos: zsum,ne.real imag = ",d16.8," ilay = ",i6)') &
                    & zsum_layer(ilay,2,1), ilay
            endif
!!$            write(6,'(" !!ldos zsum_layer = ",d16.8," ilay = ",i6)') &
!!$                 & zsum_layer(ilay,1), ilay
         end do
      end if
#endif

      if ( sw_lband == ON ) then
         ib_t = map_z(ib)
      else
         if(ekmode == ON.or.serial) then
            ib_t = ib
         else
            ib_t = map_z(ib)
         end if
      endif

      if ( map_ek(ib,ik) == mype) then
         do ilay = 1, mlayer
            weilay_noncl(ilay,ib_t,ik,:) = weilay_noncl(ilay,ib_t,ik,:) &
                 &                           + zsum_layer(ilay,1,:)*ss
         end do
      endif

      deallocate(zsum_layer)

    end subroutine substitute_weilay_cd_noncl

    subroutine fillup_meshwk(mode)
!   *********************************
!     by T.Yamasaki
!           10th Jun 1992
!   *********************************
!     Revised by T. Yamasaki, Feb. 2004
!
      integer, intent(in) :: mode
!!$ASASASASAS
!!$      integer, pointer, dimension(:,:,:) :: mesh_t
      integer, allocatable :: mesh_t(:,:,:)
!!$ASASASASAS
      integer :: i,j,k, icount, nwsc, max_nwsc, min_nwsc, n, nl, nm,nn, nd2,id, nfftwk &
           &    , inew, jnew, knew, idr2, idh

      if(mode == SOFTPART) then

         nl = fft_box_size_WF(1,1)
         nm = fft_box_size_WF(2,1)
         nn = fft_box_size_WF(3,1)

         id  = fft_box_size_WF(1,0)
         nd2 = fft_box_size_WF(2,0)
         allocate(mesh_t(id,nd2,nn))
         mesh_t = mesh

         nfftwk = product(fft_box_size_WF(1:3,1))
      else if(mode == HARDPART) then

         nl = fft_box_size_CD(1,1)
         nm = fft_box_size_CD(2,1)
         nn = fft_box_size_CD(3,1)

         id  = fft_box_size_CD_nonpara(1,0)
         nd2 = fft_box_size_CD_nonpara(2,0)
         allocate(mesh_t(id,nd2,nn))
         mesh_t = meshp
         nfftwk = product(fft_box_size_CD(1:3,1))
      end if
      idr2 = (nl+2)/2
      idh  = id/2

      if(allocated(meshwk)) deallocate(meshwk)
      allocate(meshwk(nfftwk,2)); meshwk = 0

      if(ipridos>=2) write(nfout,'(" !! nfftwk = ",i8)') nfftwk

      max_nwsc = 0
      min_nwsc = natm2+2
      icount = 1
      do k = 1, nn
         do j = 1, nm
            do i = 1, nl
!!$               if(icount > nfftwk) then
!!$                  write(nfout,'(" ! icount = ",i7," > nfftwk")') icount
!!$                  stop  ' large icount <<fillup_meshwk.m_Ldos_cal>>'
!!$               end if
               nwsc = mesh_t(i,j,k)
               if(nwsc > max_nwsc) max_nwsc = nwsc
               if(nwsc < min_nwsc) min_nwsc = nwsc
               meshwk(icount,1) = nwsc
               if(kimg == 1) then
                  if( i <= idr2 ) then
                     n = idh*nd2*(k-1) + idr2*(j-1) + i
                  else
                     inew = nl + 2 - i
                     jnew = nm + 2 - j
                     knew = nn + 2 - k
                     if(jnew > nm) jnew = jnew - nm
                     if(knew > nn) knew = knew - nn
                     n = idh*nd2*(knew-1) + idh*(jnew-1) + inew
                  end if
               else
                  n = id*nd2*(k-1) + id*(j-1) + i
               end if
               meshwk(icount,2) = n
               icount = icount + 1
            end do
         end do
      end do
      if(ipridos>=2) then
         write(nfout,*) ' !!ldos:  max_nwsc =  ', max_nwsc
         write(nfout,*) ' !!ldos:  min_nwsc =  ', min_nwsc
#ifdef DEBUG_MESH
         n = 1000
         if(n*2 > nfftwk) n = nfftwk/2
         write(nfout,'(" !ldos nfftwk = ",i8)') nfftwk
         write(nfout,'(" !ldos: meshwk(*,1) first ",i8," elements")') n
         write(nfout,'(14i5)') (meshwk(i,1),i=1,n)
         write(nfout,'(" !ldos: meshwk(*,1) last  ",i8," elements")') n
         write(nfout,'(14i6)') (meshwk(i,1),i=nfftwk-n+1, nfftwk)

         write(nfout,'(" !ldos: meshwk(*,2) first ",i8," elements")') n
         write(nfout,'(14i5)') (meshwk(i,2),i=1,n)
         write(nfout,'(" !ldos: meshwk(*,2) last  ",i8," elements")') n
         write(nfout,'(14i6)') (meshwk(i,2),i=nfftwk-n+1, nfftwk)
#endif
      end if
      deallocate(mesh_t)
    end subroutine fillup_meshwk

  end subroutine m_Ldos_cal_noncl
! ================================================================= 11.0
!$$#endif

  subroutine m_Ldos_get_dos_weight(aldos_or_layerdos,ip,nfldos,ne,nk,dos_weight)
    integer, intent(in) :: aldos_or_layerdos, ip, nfldos
    integer, intent(in) :: ne, nk
    real(kind=DP), intent(out), dimension(ne,nk) :: dos_weight

    real(kind=DP), allocatable, dimension(:,:) :: dos_weight_mpi

    integer :: i, j
    real(kind=DP) :: w
    logical :: serial
    serial = .false.
    if(hardpart_subroutine/=2.and.sw_rspace_ldos==OFF) serial=.true.

    if(ne < neg) then
       if(ipridos>=1) write(nfout,'(" !ldos: ne (= ",i6," ) < neg (= ",i6 &
            & ," ) <<m_Ldos_get_dos_weight>>")') ne, neg
       call phase_error_with_msg(nfout,' ne < neg <<m_Ldos_get_dos_weight>>',__LINE__,__FILE__)
    end if

    if(ekmode == OFF) then
       if(nk < kv3_ek) then
          if(ipridos>=1) write(nfout,'(" !ldos: nk (= ",i6," ) < kv3_ek (= ",i6 &
               & ," ) <<m_Ldos_get_dos_weight>>")') nk, kv3_ek
          call phase_error_with_msg(nfout,' nk < kv3_ek <<m_Ldos_get_dos_weight>>',__LINE__,__FILE__)
       end if

!       if(ne <= iend_e .and. nk <= iend_k) then
        if(ne <= neg .and. nk <= kv3) then
          dos_weight = 0.d0
          if(aldos_or_layerdos == ALDOS) then
             if(.not.serial)then
               do j = ista_k, iend_k
                  do i = ista_e, iend_e, istep_e
                     dos_weight(i,j) = weiwsc(ip,map_z(i),j)
                  end do
               end do
             else
               do j = 1,kv3
                  do i = 1,neg
                     dos_weight(i,j) = weiwsc(ip,i,j)
                  end do
               end do
             endif
          else if(aldos_or_layerdos == LAYERDOS) then
             if(.not.serial)then
               do j = ista_k, iend_k
                  do i = ista_e, iend_e, istep_e
                     dos_weight(i,j) = weilay(ip,map_z(i),j)
                  end do
               end do
             else
               do j = 1, kv3
                  do i = 1,neg
                     dos_weight(i,j) = weilay(ip,i,j)
                  end do
               end do
             endif
          end if
          if(npes > 1.and..not.serial) then
            allocate(dos_weight_mpi(ne,nk))
            call mpi_allreduce(dos_weight, dos_weight_mpi,ne*nk &
                 & ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
            dos_weight = dos_weight_mpi
            deallocate(dos_weight_mpi)
          end if
        end if
    else if(ekmode == ON) then
       if(nk < kv3) then
          if(ipridos>=1) write(nfout,'(" !ldos: nk (= ",i6," ) < kv3 (= ",i6 &
               & ," ) <<m_Ldos_get_dos_weight>>")') nk, kv3_ek
          call phase_error_with_msg(nfout,' nk < kv3 <<m_Ldos_get_dos_weight>>',__LINE__,__FILE__)
       end if
       dos_weight(:,:) = dos_weight_ek(aldos_or_layerdos,ip,:,:)

!       call m_Files_open_nfldos()
!       call skip_headerpart()
!       if(aldos_or_layerdos == ALDOS) then
!          do j = 1, nk
!             do i = 1, ne
!                call read_aldos_weight(ip,w)
!                dos_weight(i,j) = w
!             end do
!          end do
!       else if(aldos_or_layerdos == LAYERDOS) then
!          do j = 1, nk
!             do i = 1, ne
!                call read_layerdos_weight(ip,w)
!                dos_weight(i,j) = w
!             end do
!          end do
!       end if
!       call m_Files_close_nfldos()
    end if
!  contains
!    subroutine skip_headerpart()
!      read(nfldos,*)
!    end subroutine skip_headerpart
!
!    subroutine read_aldos_weight(ip,w)
!      integer, intent(in) :: ip
!      real(kind=DP), intent(out) :: w
!      w = 0.d0
!    end subroutine read_aldos_weight
!
!    subroutine read_layerdos_weight(ip,w)
!      integer, intent(in) :: ip
!      real(kind=DP), intent(out) :: w
!      w = 0.d0
!    end subroutine read_layerdos_weight

  end subroutine m_Ldos_get_dos_weight

! ============================== added by K. Tagami ================== 11.0
  subroutine m_Ldos_get_dos_weight_noncl( aldos_or_layerdos, ip, nfldos, &
       &                                  ne, nk, dos_weight_noncl )
    integer, intent(in) :: aldos_or_layerdos, ip, nfldos
    integer, intent(in) :: ne, nk
    real(kind=DP), intent(out), dimension(ne,nk,ndim_magmom) :: dos_weight_noncl

    real(kind=DP), allocatable, dimension(:,:,:) :: dos_weight_mpi

    integer :: i, j
    real(kind=DP) :: w
    logical :: serial

    serial = .false.
    if(hardpart_subroutine/=2.and.sw_rspace_ldos==OFF) serial=.true.

    if(ne < neg) then
       if(ipridos>=1) write(nfout,'(" !ldos: ne (= ",i6," ) < neg (= ",i6 &
            & ," ) <<m_Ldos_get_dos_weight>>")') ne, neg
       call phase_error_with_msg(nfout,' ne < neg <<m_Ldos_get_dos_weight>>',__LINE__,__FILE__)
    end if

    if(ekmode == OFF) then
       if(nk < kv3_ek) then
          if(ipridos>=1) write(nfout,'(" !ldos: nk (= ",i6," ) < kv3_ek (= ",i6 &
               & ," ) <<m_Ldos_get_dos_weight>>")') nk, kv3_ek
          call phase_error_with_msg(nfout,' nk < kv3_ek <<m_Ldos_get_dos_weight>>',__LINE__,__FILE__)
       end if

!       if(ne <= iend_e .and. nk <= iend_k) then
        if(ne <= neg .and. nk <= kv3) then
          dos_weight_noncl = 0.d0
          if(aldos_or_layerdos == ALDOS) then
             if(.not.serial)then
               do j = ista_k, iend_k
                  do i = ista_e, iend_e, istep_e
                     dos_weight_noncl(i,j,:) = weiwsc_noncl(ip,map_z(i),j,:)
                  end do
               end do
             else
               do j = 1,kv3
                  do i = 1,neg
                     dos_weight_noncl(i,j,:) = weiwsc_noncl(ip,i,j,:)
                  end do
               end do
             endif
          else if(aldos_or_layerdos == LAYERDOS) then
             if(.not.serial)then
               do j = ista_k, iend_k
                  do i = ista_e, iend_e, istep_e
                     dos_weight_noncl(i,j,:) = weilay_noncl(ip,map_z(i),j,:)
                  end do
               end do
             else
               do j = 1, kv3
                  do i = 1,neg
                     dos_weight_noncl(i,j,:) = weilay_noncl(ip,i,j,:)
                  end do
               end do
             endif
          end if

          if(npes > 1.and..not.serial) then
            allocate(dos_weight_mpi(ne,nk,ndim_magmom))
            call mpi_allreduce( dos_weight_noncl, dos_weight_mpi, ne*nk*ndim_magmom, &
                 &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
            dos_weight_noncl = dos_weight_mpi
            deallocate(dos_weight_mpi)
          end if
        end if

    else if(ekmode == ON) then
       if(nk < kv3) then
          if(ipridos>=1) write(nfout,'(" !ldos: nk (= ",i6," ) < kv3 (= ",i6 &
               & ," ) <<m_Ldos_get_dos_weight>>")') nk, kv3_ek
          call phase_error_with_msg(nfout,' nk < kv3 <<m_Ldos_get_dos_weight>>',__LINE__,__FILE__)
       end if

       call m_Files_open_nfldos()
       call skip_headerpart()
       if(aldos_or_layerdos == ALDOS) then
          do j = 1, nk
             do i = 1, ne
                call read_aldos_weight(ip,w)
                dos_weight_noncl(i,j,:) = w
             end do
          end do
       else if(aldos_or_layerdos == LAYERDOS) then
          do j = 1, nk
             do i = 1, ne
                call read_layerdos_weight(ip,w)
                dos_weight_noncl(i,j,:) = w
             end do
          end do
       end if
       call m_Files_close_nfldos()
    end if
  contains
    subroutine skip_headerpart()
      read(nfldos,*)
    end subroutine skip_headerpart

    subroutine read_aldos_weight(ip,w)
      integer, intent(in) :: ip
      real(kind=DP), intent(out) :: w
      w = 0.d0
    end subroutine read_aldos_weight

    subroutine read_layerdos_weight(ip,w)
      integer, intent(in) :: ip
      real(kind=DP), intent(out) :: w
      w = 0.d0
    end subroutine read_layerdos_weight

  end subroutine m_Ldos_get_dos_weight_noncl
! ==================================================================== 11.0

  subroutine m_Ldos_get_ldos_index(ip,aldos_or_layerdos,ipdos,tagwords)
    integer, intent(in) ::  ip
    integer, intent(out) :: aldos_or_layerdos, ipdos
    character(len=40), intent(out) :: tagwords
    integer :: i, na

    aldos_or_layerdos = NO
    ipdos = 0
    tagwords = ''
    if(sw_aldos == ON) then
!!$       if(1 <= ip .and. ip <= maldos ) then
       if(1 <= ip .and. ip <= naldos_write ) then
          aldos_or_layerdos = ALDOS
!!$          ipdos = naldos_from + ip - 1
          na = 0
          do i = naldos_from, naldos_to
             if(if_aldos_full(i) == ON) na = na + 1
             if(na == ip) then
                ipdos = i
                exit
             end if
          end do
          if(ipdos > maldos .or. ipdos < 1 ) then
             write(nfout,'(" !! ipdos ( = ",i8," ) is illegal << m_Ldos_get_ldos_index>>")') ipdos
             call phase_error_with_msg(nfout,' ipdos is illegal << m_Ldos_get_ldos_index>>',__LINE__,__FILE__)
          end if
       else if( ip > naldos_write ) then
          if(sw_layerdos == ON) then
             aldos_or_layerdos = LAYERDOS
             ipdos = ip - naldos_write
             write(tagwords,999) winlay(ipdos,1)*BOHR,winlay(ipdos,2)*BOHR
          end if
       end if
    else if(sw_layerdos == ON) then
       if( 1 <= ip .and. ip <= mlayer) then
          aldos_or_layerdos = LAYERDOS
          ipdos = ip
          write(tagwords,999) winlay(ipdos,1)*BOHR,winlay(ipdos,2)*BOHR
       end if
    end if
999 format(' range = [',f8.4,' : ',f8.4,']  (Angst.)')

  end subroutine m_Ldos_get_ldos_index

  subroutine m_Ldos_load_dos_weight()
    use m_Const_Parameters, only :FMAXTAGLEN
    integer :: nldos,imode,ik,kread,m,i,ierr
    real(kind=DP),allocatable,dimension(:) :: w
    integer :: ios
    logical :: eof

    call m_Files_open_nfldos()
    nldos = maldos
    if(mlayer>nldos) nldos = mlayer
    if(.not.allocated(dos_weight_ek)) allocate(dos_weight_ek(2,nldos,neg,kv3_ek))
    dos_weight_ek = 0.d0

    if (mype==0) then
      rewind(nfldos)
      kread = 0
      do i=1, 4
         read(nfldos,*,iostat=ios)
      end do

      do
         read(nfldos,*,iostat=ios) imode,ik
         if(ios==-1) then
           write(nfout,*) '!** EOF reached'
           exit
         endif
         if(imode==ALDOS) then
           m = maldos
         else if (imode==LAYERDOS) then
           m = mlayer
         endif
         allocate(w(m));w=0.d0

         do i = 1, neg -num_extra_bands
            call read_dos_weight(m,w,eof)
            dos_weight_ek(imode,1:m,i,ik) = w(1:m)
            if(eof) then
               write(nfout,*) '!** EOF reached'
               exit
            endif
         end do
         deallocate(w)
         kread = kread+1
!         if (kread == kv3_ek) exit
      enddo
    endif

    call m_Files_close_nfldos()
    if(npes>1) call mpi_bcast(dos_weight_ek,neg*kv3_ek*2*nldos, &
                    mpi_double_precision,0,MPI_CommGroup,ierr)

    contains

    subroutine read_dos_weight(m,w,eof)
      integer, intent(in) :: m
      real(kind=DP), dimension(m), intent(out) :: w
      logical, intent(out) :: eof
      character(len=FMAXTAGLEN) :: chr
      integer :: j,ios
      eof = .false.
      read(nfldos,*,iostat=ios) chr,(w(j),j=1,m)
      if(ios == -1) then
         eof = .true.
         return
      endif
    end subroutine read_dos_weight

  end subroutine m_Ldos_load_dos_weight


  subroutine m_Ldos_dealloc()
     if(allocated(if_aldos_full)) deallocate(if_aldos_full)
     if(allocated(mesh)) deallocate(mesh)
     if(allocated(meshp)) deallocate(meshp)
     if(allocated(winlay)) deallocate(winlay)
     if(allocated(ac_mesh)) deallocate(ac_mesh)
     if(allocated(nac_mesh)) deallocate(nac_mesh)
     if(allocated(ac_mesh_cd)) deallocate(ac_mesh_cd)
     if(allocated(nac_mesh_cd)) deallocate(nac_mesh_cd)
  end subroutine m_Ldos_dealloc

  subroutine get_fftbox_size(mode,nl0,nm0,nl,nm,nn)
    integer, intent(in)  :: mode
    integer, intent(out) :: nl0,nm0,nl,nm,nn
    if(mode==SOFTPART) then
       nl0 = fft_box_size_WF(1,0)
       nm0 = fft_box_size_WF(2,0)
       nl  = fft_box_size_WF(1,1)
       nm  = fft_box_size_WF(2,1)
       nn  = fft_box_size_WF(3,1)
    else if(mode==HARDPART) then
       nl0 = fft_box_size_CD_nonpara(1,0)
       nm0 = fft_box_size_CD_nonpara(2,0)
       nl  = fft_box_size_CD(1,1)
       nm  = fft_box_size_CD(2,1)
       nn  = fft_box_size_CD(3,1)
    end if
  end subroutine get_fftbox_size

end module m_Ldos

