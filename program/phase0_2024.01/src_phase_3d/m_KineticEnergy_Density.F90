#ifdef __TIMER_SUB__
#   define __TIMER_SUB_START(a)  call timer_sta(a)
#   define __TIMER_SUB_STOP(a)   call timer_end(a)
#else
#   define __TIMER_SUB_START(a)
#   define __TIMER_SUB_STOP(a)
#endif
#ifdef __TIMER_DO__
#   define __TIMER_DO_START(a)   call timer_sta(a)
#   define __TIMER_DO_STOP(a)    call timer_end(a)
#else
#   define __TIMER_DO_START(a)
#   define __TIMER_DO_STOP(a)
#endif
#ifdef __TIMER_IODO__
#   define __TIMER_IODO_START(a)   call timer_sta(a)
#   define __TIMER_IODO_STOP(a)    call timer_end(a)
#else
#   define __TIMER_IODO_START(a)
#   define __TIMER_IODO_STOP(a)
#endif
#ifdef __TIMER_COMM__
#   define __TIMER_COMM_START_w_BARRIER(str,a)   call timer_barrier(str) ;   call timer_sta(a)
#   define __TIMER_COMM_START(a)       call timer_sta(a)
#   define __TIMER_COMM_STOP(a)        call timer_end(a)
#else
#   define __TIMER_COMM_START_w_BARRIER(str,a)
#   define __TIMER_COMM_START(a)
#   define __TIMER_COMM_STOP(a)
#endif
#ifdef __TIMER_IOCOMM__
#   define __TIMER_IOCOMM_START_w_BARRIER(str,a)   call timer_barrier(str) ;   call timer_sta(a)
#   define __TIMER_IOCOMM_START(a)       call timer_sta(a)
#   define __TIMER_IOCOMM_STOP(a)        call timer_end(a)
#else
#   define __TIMER_IOCOMM_START_w_BARRIER(str,a)
#   define __TIMER_IOCOMM_START(a)
#   define __TIMER_IOCOMM_STOP(a)
#endif

module m_KineticEnergy_Density
  use m_Control_Parameters,  only : nspin, kimg, noncol, &
       &                            use_symm_ekin_density, use_asymm_ekin_density, &
       &                            af

  use m_KPoints,             only : kv3, vkxyz
  use m_Const_Parameters,    only : DP, DIRECT, OFF, ELECTRON, DELTA, BUCS, PAI2, NO, &
       &                            ANTIFERRO
  use m_PlaneWaveBasisSet,    only : kgp, kg1, ngabc, iba, nbase, igf, kg, gr_l, &
       &                             ngpt_l, ista_k, iend_k
  use m_Crystal_Structure,    only : rltv, univol, tau, nopr
  use m_Electronic_Structure, only : zaj_l, occup_l, m_ES_WF_in_Rspace_3D1
  use m_FFT,                  only : m_FFT_alloc_WF_work, m_FFT_dealloc_WF_work, &
       &                             nfft, m_FFT_WF, fft_box_size_WF
  use m_Parallelization,      only : np_e, map_z, map_k, myrank_k, MPI_CommGroup, &
       &                             npes, ista_kngp, iend_kngp, mype
  use m_Files,                only : nfout

! === KT_add === 2014/09/19
  use m_Const_Parameters,   only : YES, NO, EXECUT, SKIP, zi, SIMPLE_CUBIC, PAI, ON
  use m_Control_Parameters,  only : m_CtrlP_cachesize, ndim_magmom, nel_ylm, &
       &                            ekin_density_type, sw_rspace_ekin_density, &
       &                            sw_calc_ekin_density_hardpart, &
       &                            sw_add_ekin_hardpart_on_Gspace, &
       &                            sw_fft_xzy, sw_communicator_for_chg, &
       &                            nblocksize_fftw_is_given, nblocksize_fftw
  use m_Parallelization,     only : nel_fft_z , nel_fft_y, nel_fft_x &
       &                          , np_fft_x, np_fft_y, np_fft_z  &
       &                          , fft_X_x_nel, fft_X_y_nel, fft_X_z_nel  &
       &                          , xyz_fft_x, mp_fft_x, ista_e, iend_e, istep_e &
       &                          , nrank_e, myrank_g, nrank_g, ista_k, iend_k   &
       &                          , nrank_chg, myrank_chg, np_g1k, np_kngp_gw   &
       &                          , ista_kngp_gw, iend_kngp_gw &
       &                          , mpi_ke_world, mpi_chg_world, mpi_kg_world &
       &                          , mpi_k_world, mpi_ge_world, ierr &
       &                          , ista_g1k , iend_g1k, np_g1k, is_kngp, ie_kngp
  use m_Ionic_System,           only : ntyp, ityp, natm, iwei, pos
  use m_Crystal_Structure,    only : nbztyp
  use m_PseudoPotential,    only :  kina_qitg_l, kins_qitg_l, modnrm, nqitg, &
       &                            m_PP_include_vanderbilt_pot, nlmt, dl2p, &
       &                            m_PP_find_maximum_l, &
       &                            m_PP_set_index_arrays1, m_PP_set_index_arrays2, &
       &                            mmesh, nmesh, xh, rmax, ilmt, ltp, mtp, taup, &
       &                            psirpw, phirpw
  use m_PlaneWaveBasisSet,    only : ylm_l, m_pwBS_sphrp2, fp_l
  use m_Charge_Density,       only : map_fft_to_chgq_3D

  use m_FFT,                 only : m_FFT_Direct_3D, m_FFT_Direct_XYZ_3D
#ifdef FFT_3D_DIVISION
  use m_FFT,                 only : m_FFT_Direct_3DIV_3D
#endif
  use m_Realspace,            only : nmesh_rs_aug_max, nmesh_rs_aug, &
       &                             meshx_rs_aug, meshy_rs_aug, meshz_rs_aug, &
       &                             meshxyz_rs_aug, m_RS_R_minus_pos
! ============== 2014/09/19

  use m_Files,  only : nf_ekindens
  use m_Const_Parameters,   only : CMPLDP
  use mpi

  implicit none
!  include 'mpif.h'

  real(kind=DP), allocatable, target :: ekins_l(:,:,:)  ! symmetric, positive definite
  real(kind=DP), allocatable, target :: ekina_l(:,:,:)  ! asymmetric

  real(kind=DP), allocatable, target :: ekins_old(:,:,:)
  real(kind=DP), allocatable, target :: ekina_old(:,:,:)

  real(kind=DP), allocatable :: kins_qrs_rspace(:,:,:,:)
  real(kind=DP), allocatable :: kina_qrs_rspace(:,:,:,:)

  real(kind=DP), public, pointer, dimension(:,:)        :: work

contains

  subroutine m_KE_alloc_ekin_density
    if ( use_symm_ekin_density ) then
       if ( noncol ) then
       else
          allocate( ekins_l(ista_kngp:iend_kngp,kimg,nspin) )
          allocate( ekins_old(ista_kngp:iend_kngp,kimg,nspin) )
       endif
       ekins_l = 0.0d0;   ekins_old = 0.0d0
    endif

    if ( use_asymm_ekin_density ) then
       if ( noncol ) then
       else
          allocate( ekina_l(ista_kngp:iend_kngp,kimg,nspin) )
          allocate( ekina_old(ista_kngp:iend_kngp,kimg,nspin) )
       endif
       ekina_l = 0.0d0;   ekina_old = 0.0d0
    endif
  end subroutine m_KE_alloc_ekin_density

  subroutine m_KE_dealloc_ekin_density
    if ( allocated( ekins_l ) ) deallocate( ekins_l )
    if ( allocated( ekina_l ) ) deallocate( ekina_l )
    if ( allocated( ekins_old ) ) deallocate( ekins_old )
    if ( allocated( ekina_old ) ) deallocate( ekina_old )

    if ( allocated( kins_qrs_rspace ) ) deallocate( kins_qrs_rspace )
  end subroutine m_KE_dealloc_ekin_density

  subroutine m_KE_cp_ekin_density_to_old
    if ( use_symm_ekin_density )  ekins_old = ekins_l
    if ( use_asymm_ekin_density ) ekina_old = ekina_l
  end subroutine m_KE_cp_ekin_density_to_old


  subroutine m_KE_symm_softpart_3D
    integer :: ispin, ik, ib, ixyz, ierr
    integer :: i, ig, ri, iend, i1
    integer :: isrsize, fft_l_size
    integer :: ib1, ib2, ibsize, lsize, ibesize

    real(kind=DP) :: fac, occupation
    real(kind=DP), allocatable :: qxyz(:,:), psi_l(:,:,:,:)
    real(kind=DP), allocatable,dimension(:,:) :: afft_l, afft_l_mpi
    real(kind=DP), allocatable,dimension(:,:) :: bfft_l, map_afft_l

    integer       :: ist,ien
    real(kind=DP), pointer, dimension(:,:,:) :: ekinq_p
    real(kind=DP), allocatable, target, dimension(:,:,:) :: ekinqtmp

    allocate( qxyz(maxval(np_g1k),3 ) ); qxyz = 0.0d0
    ekins_l = 0.0d0

#ifdef FFT_3D_DIVISION
    lsize = fft_X_x_nel*fft_X_y_nel*fft_X_z_nel
#else
    lsize = max(maxval(nel_fft_x(:)),maxval(nel_fft_y(:)),maxval(nel_fft_z(:)))
#endif

    ibsize = 1
    if (nblocksize_fftw_is_given) then
       ibsize = nblocksize_fftw
       if (ibsize < 1) ibsize = 1
    endif

#ifdef FFT_3D_DIVISION
    allocate(afft_l(lsize*2,1) ,stat=ierr)
    allocate(bfft_l(lsize*2,ibsize) ,stat=ierr)
    afft_l = 0.0d0
    bfft_l = 0.0d0
#else
    allocate(afft_l(lsize*kimg,1) ,stat=ierr)
    allocate(bfft_l(lsize*kimg,ibsize) ,stat=ierr)
#endif

    fac = 2.d0/(univol*kv3*product(fft_box_size_WF(1:3,1)))

    Do ispin=1, nspin
       afft_l = 0.0d0
       do ik = ispin, kv3+ispin-nspin, nspin
          if ( map_k(ik) /= myrank_k ) cycle

          call k_plus_G_vectors_m_3D( ik, kg, kg1, kv3, iba, nbase, vkxyz, ngabc, rltv,&
               &                      qxyz(:,1), qxyz(:,2), qxyz(:,3) )

          allocate( psi_l( maxval(np_g1k), np_e, ik:ik, kimg ) ); psi_l = 0.0d0

          Do ib1 = 1, np_e, ibsize
             ib2 = min(ib1+ibsize-1,np_e)
             ibesize = ib2 - ib1 + 1

             occupation = occup_l( ib1,ik )
             if ( occupation < DELTA ) cycle

             Do ixyz=1, 3
                psi_l = 0.0d0
                if ( kimg == 1 ) then
                   call phase_error_with_msg(nfout,"Not supported",__LINE__,__FILE__)
                else
                   Do ig=1, np_g1k(ik)
                      psi_l(ig,ib1:ib2,ik,2) =  zaj_l(ig,ib1:ib2,ik,1) &
                           &                       *qxyz(ig,ixyz)
                      psi_l(ig,ib1:ib2,ik,1) = -zaj_l(ig,ib1:ib2,ik,2) &
                           &                       *qxyz(ig,ixyz)
                   End do
                endif

                bfft_l = 0.0d0
                call m_ES_WF_in_Rspace_3D1( ik, ik, ik, ib1, ib2, ibsize, lsize, &
                     &                      psi_l, bfft_l )
                call add_occupied_densities()
             End Do
          End do
          deallocate( psi_l )
       End do

#ifdef FFT_3D_DIVISION
       allocate(afft_l_mpi(lsize*2,1) ,stat=ierr)
#else
       allocate(afft_l_mpi(lsize*kimg,1) ,stat=ierr)
#endif

#ifdef FFT_3D_DIVISION
       call mpi_allreduce(afft_l,afft_l_mpi,lsize*2,mpi_double_precision,   &
            &                   mpi_sum,mpi_kg_world,ierr)
       afft_l = afft_l_mpi
       call mpi_allreduce(afft_l,afft_l_mpi,lsize*2,mpi_double_precision,   &
            &                   mpi_sum,mpi_ge_world,ierr)
#else
       if (sw_fft_xzy > 0) then
          call mpi_allreduce(afft_l,afft_l_mpi,np_fft_y*kimg,mpi_double_precision,   &
               &                   mpi_sum,mpi_kg_world,ierr)
          afft_l = afft_l_mpi
          call mpi_allreduce(afft_l,afft_l_mpi,np_fft_y*kimg,mpi_double_precision,   &
               &                   mpi_sum,mpi_ge_world,ierr)
       else
          call mpi_allreduce(afft_l,afft_l_mpi,lsize*kimg,mpi_double_precision,   &
               &                   mpi_sum,mpi_kg_world,ierr)
          afft_l = afft_l_mpi
          call mpi_allreduce(afft_l,afft_l_mpi,lsize*kimg,mpi_double_precision,   &
               &                   mpi_sum,mpi_ge_world,ierr)
       end if
#endif
       __TIMER_COMM_STOP(825)
       afft_l = afft_l_mpi

       deallocate(afft_l_mpi)

#ifdef FFT_3D_DIVISION
       call m_FFT_Direct_3DIV_3D (nfout, afft_l, lsize, 1)
#else
       if (sw_fft_xzy > 0) then
          call m_FFT_Direct_3D (nfout, afft_l, lsize, 1)
       else
          call m_FFT_Direct_XYZ_3D (nfout, afft_l, lsize, 1)
       end if
#endif
#ifdef FFT_3D_DIVISION
       allocate(map_afft_l(np_kngp_gw*2,1) ,stat=ierr)
#else
       allocate(map_afft_l(np_kngp_gw*kimg,1) ,stat=ierr)
#endif
       call map_fft_to_chgq_3D(lsize, 1, afft_l, map_afft_l, nfout)
       call substitute_KE_for_ekinq
       deallocate(map_afft_l)
    End Do

    ekins_l = ekins_l /2.0d0

    call symmetrize_ekin_3D( NO, ekins_l )

    deallocate( afft_l ); deallocate(bfft_l); deallocate( qxyz )

  contains

    subroutine add_occupied_densities
      integer  :: i, ib
      real(kind=DP) :: occupation
                                                 __TIMER_SUB_START(717)
                                                 __TIMER_DO_START(826)
      do ib = ib1, ib2
         occupation = occup_l(ib,ik)
         if(occupation < DELTA) cycle
#ifdef FFT_3D_DIVISION
         do i = 1, lsize*2-1, 2
            afft_l(i,1) = afft_l(i,1) + occupation*(bfft_l(i,ib-ib1+1)**2+bfft_l(i+1,ib-ib1+1)**2) ! MPI
         end do
#else
         if (sw_fft_xzy > 0) then
            do i = 1, np_fft_y*kimg-1, 2
               afft_l(i,1) = afft_l(i,1) + occupation*(bfft_l(i,ib-ib1+1)**2+bfft_l(i+1,ib-ib1+1)**2) ! MPI
            end do
         else
            do i = 1, np_fft_z*kimg-1, 2
               afft_l(i,1) = afft_l(i,1) + occupation*(bfft_l(i,ib-ib1+1)**2+bfft_l(i+1,ib-ib1+1)**2) ! MPI
            end do
         end if
#endif
      end do
                                                 __TIMER_DO_STOP(826)
                                                 __TIMER_SUB_STOP(717)
    end subroutine add_occupied_densities

    subroutine substitute_KE_for_ekinq

      real(kind=DP) :: fac
      integer       :: i, iend !mpi
      integer       :: ist,ien
      real(kind=DP), pointer, dimension(:,:,:) :: ekinq_p
      real(kind=DP), allocatable, target, dimension(:,:,:) :: ekinqtmp

      if(sw_communicator_for_chg == ON)then
          allocate(ekinqtmp(ista_kngp_gw:iend_kngp_gw,kimg,nspin));ekinqtmp = 0.d0
          ekinq_p => ekinqtmp
          ist = ista_kngp_gw
          ien = iend_kngp_gw
      else
          ekinq_p => ekins_l
          ist = ista_kngp
          ien = iend_kngp
      endif
                                                 __TIMER_SUB_START(719)
      fac = 2.d0/(univol*kv3*product(fft_box_size_WF(1:3,1)))
      iend = ien
      if( iend > kg ) iend = kg
                                                 __TIMER_DO_START(830)
      if( ist <= iend ) then
         if (kimg == 1) then
            do i = ist, iend  !for mpi
               ekinq_p(i,1,ispin) = map_afft_l(i-ist+1,1)*fac
            end do
         else
            do i = ist, iend  !for mpi
               ekinq_p(i,1,ispin) = map_afft_l((i-ist+1)*2-1,1)*fac
               ekinq_p(i,2,ispin) = map_afft_l((i-ist+1)*2  ,1)*fac
            end do
         end if
      endif
                                                 __TIMER_DO_STOP(830)
                                                 __TIMER_SUB_STOP(719)
      if(sw_communicator_for_chg == ON)then
        do i=ista_kngp,iend_kngp
           ekins_l(i,:,:) = ekinq_p(i,:,:)
        enddo
        deallocate(ekinqtmp)
      endif
    end subroutine substitute_KE_for_ekinq

  end subroutine m_KE_symm_softpart_3D

  subroutine symmetrize_ekin_3D(mode,ekin)
    integer,intent(in)           :: mode
    real(kind=DP), intent(inout) :: ekin(ista_kngp:iend_kngp,kimg,nspin)
    integer ::       ispin, ng, no, ngp, no1, no2
    real(kind=DP) :: fi, tx,ty,tz, fp, fc, fs, zcr, zci
    real(kind=DP), pointer, dimension(:,:) :: work2
    logical, save                           :: firstcall = .true.
    integer, save                           :: sendmax, recvmax, sendranks, recvranks
    integer                                 :: mymin, mymax, lrank
    integer,                    dimension(2)   :: my_from_to
    integer,       allocatable, dimension(:,:) :: all_from_to
    integer,       allocatable, dimension(:)   :: is_ngpt, ie_ngpt
    integer, save, allocatable, dimension(:,:) :: sendinfo, recvinfo
                                                 __TIMER_SUB_START(731)
    allocate(work2(ista_kngp:iend_kngp,kimg)); work2 = 0.d0

    if(mode == ANTIFERRO) then
       fi = 1.d0/af
       no1 = nopr + 1; no2 = nopr + af
    else
       fi = 1.d0/nopr
       no1 = 1; no2 = nopr
    end if

!XX!if (firstcall) then
    allocate(all_from_to(2,0:nrank_chg-1))
    allocate(is_ngpt(0:nrank_chg-1))
    allocate(ie_ngpt(0:nrank_chg-1))
    allocate(sendinfo(2,0:nrank_chg-1))
    allocate(recvinfo(2,0:nrank_chg-1))
    sendinfo = 0
    recvinfo = 0
    sendmax = 0
    recvmax = 0
    sendranks = 0
    recvranks = 0

    mymin = kgp
    mymax = 1
                                                 __TIMER_DO_START(848)
    do no = no1, no2
       do ng = ista_kngp, iend_kngp !for mp
          ngp = ngpt_l(ng,no)
          if (mymin > ngp) mymin = ngp
          if (mymax < ngp) mymax = ngp
       end do
    end do
                                                 __TIMER_DO_STOP(848)

    my_from_to(1) = mymin
    my_from_to(2) = mymax

                                                 __TIMER_COMM_START_w_BARRIER(mpi_chg_world,849)
    call mpi_allgather(my_from_to , 2, MPI_INTEGER, &
   &                   all_from_to, 2, MPI_INTEGER, mpi_chg_world, ierr)
                                                 __TIMER_COMM_STOP(849)
    is_ngpt(0:nrank_chg-1) = all_from_to(1,0:nrank_chg-1)
    ie_ngpt(0:nrank_chg-1) = all_from_to(2,0:nrank_chg-1)

                                                 __TIMER_DO_START(850)
    do lrank = 0, nrank_chg-1
      if(lrank == myrank_chg) cycle
      if(ie_ngpt(lrank) < ista_kngp) cycle
      if(iend_kngp < is_ngpt(lrank)) cycle
      if(ista_kngp <= ie_ngpt(lrank)) then
         sendinfo(1,lrank) = max(ista_kngp,is_ngpt(lrank))
         sendinfo(2,lrank) = min(iend_kngp,ie_ngpt(lrank))
         sendranks = sendranks + 1
      else if(is_ngpt(lrank) <= iend_kngp) then
         sendinfo(1,lrank) = max(ista_kngp,is_ngpt(lrank))
         sendinfo(2,lrank) = min(iend_kngp,ie_ngpt(lrank))
         sendranks = sendranks + 1
      endif
      sendmax = max(sendmax,sendinfo(2,lrank)-sendinfo(1,lrank)+1)
    end do
                                                 __TIMER_DO_STOP(850)
                                                 __TIMER_DO_START(851)
    do lrank = 0, nrank_chg-1
      if(lrank == myrank_chg) cycle
      if(ie_kngp(lrank) < mymin) cycle
      if(mymax < is_kngp(lrank)) cycle
      if(mymin <= ie_kngp(lrank)) then
         recvinfo(1,lrank) = max(mymin,is_kngp(lrank))
         recvinfo(2,lrank) = min(mymax,ie_kngp(lrank))
         recvranks = recvranks + 1
      else if(is_kngp(lrank) <= mymax) then
         recvinfo(1,lrank) = max(mymin,is_kngp(lrank))
         recvinfo(2,lrank) = min(mymax,ie_kngp(lrank))
         recvranks = recvranks + 1
      endif
      recvmax = max(recvmax,recvinfo(2,lrank)-recvinfo(1,lrank)+1)
    end do
                                                 __TIMER_DO_STOP(851)
    deallocate(all_from_to)
    deallocate(is_ngpt)
    deallocate(ie_ngpt)

!XX!firstcall = .false.
!XX!end if

    allocate(work(mymin:mymax,kimg))

    do ispin = 1, nspin, af+1
       call cp_chgq_by_ngpt() ! chg -> work
       work2 = 0.d0                   ! initialization
                                                 __TIMER_DO_START(852)
       do no = no1, no2
!!$          tx = tau(1,no,BUCS)*PAI2
!!$          ty = tau(2,no,BUCS)*PAI2
!!$          tz = tau(3,no,BUCS)*PAI2
          if(kimg == 1) then
             do ng = ista_kngp, iend_kngp !for mpi
                ngp = ngpt_l(ng,no)
!                fp = ngabc(ngp,1)*tx + ngabc(ngp,2)*ty + ngabc(ngp,3)*tz
!                fp = ngabc_kngp_l(ngp,1)*tx + ngabc_kngp_l(ngp,2)*ty + ngabc_kngp_l(ngp,3)*tz
                fp = fp_l(ng,no)
                work2(ng,1)        = work2(ng,1) + dcos(fp)*work(ngp,1)
             end do
          else if(kimg == 2) then
             do ng = ista_kngp, iend_kngp !for mpi
                ngp= ngpt_l(ng,no)
!                fp = ngabc(ngp,1)*tx + ngabc(ngp,2)*ty + ngabc(ngp,3)*tz
!                fp = ngabc_kngp_l(ngp,1)*tx + ngabc_kngp_l(ngp,2)*ty + ngabc_kngp_l(ngp,3)*tz
                fp = fp_l(ng,no)
                fc = dcos(fp);     fs = dsin(fp)
                zcr= work(ngp,1);  zci= work(ngp,kimg)
                work2(ng,1)        = work2(ng,1) + fc*zcr - fs*zci
                work2(ng,2)        = work2(ng,2) + fc*zci + fs*zcr
             end do
          end if
       end do
                                                 __TIMER_DO_STOP(852)
       if(mode /= ANTIFERRO) ekin(:,:,ispin) = work2(:,:)*fi
    end do

    if(mode == ANTIFERRO) ekin(:,:,nspin) = work2(:,:)*fi

    deallocate(work2)
    deallocate(sendinfo)
    deallocate(recvinfo)
    deallocate(work)
                                                 __TIMER_SUB_STOP(731)
   contains

    subroutine cp_chgq_by_ngpt()
      real(DP), allocatable, dimension(:,:) :: sendbuf
      real(DP), allocatable, dimension(:,:) :: recvbuf
      integer, dimension(sendranks) :: req_s
      integer, dimension(recvranks) :: req_r, src
      integer, dimension(MPI_STATUS_SIZE,sendranks) :: sta_s
      integer, dimension(MPI_STATUS_SIZE,recvranks) :: sta_r

      integer :: i, j, ri, ista, iend, nel
      integer :: icnt_s, icnt_r, ierr, itag=30
! === DEBUG by tkato 2012/06/04 ================================================
#ifndef USE_NONBLK_COMM
#ifndef USE_ALLTOALLV
      real(kind=DP), allocatable, dimension(:,:) :: sbuf, rbuf
      integer :: maxbuf
#else
      integer, allocatable, dimension(:) :: sdsp, rdsp
      integer, allocatable, dimension(:) :: scnt, rcnt
      real(kind=DP), allocatable, dimension(:,:) :: sbuf, rbuf
#endif
#endif
! ==============================================================================

                                                 __TIMER_SUB_START(732)
      allocate(recvbuf(recvmax*kimg,recvranks))
      allocate(sendbuf(sendmax*kimg,sendranks))

      icnt_r = 0
#ifdef USE_NONBLK_COMM
                                                 __TIMER_COMM_START_w_BARRIER(mpi_chg_world,853)
#endif
      do lrank = 0, nrank_chg - 1
         if(recvinfo(1,lrank) == 0) cycle
         ista = recvinfo(1,lrank)
         iend = recvinfo(2,lrank)
         nel  = iend-ista+1
         icnt_r = icnt_r + 1
         src(icnt_r) = lrank
#ifdef USE_NONBLK_COMM
         call mpi_irecv(recvbuf(1,icnt_r), nel*kimg, mpi_double_precision, &
        &               lrank, itag, mpi_chg_world, req_r(icnt_r), ierr)
          if (ierr /= 0) then
             call mpi_abort(mpi_comm_world, 171, ierr)
          endif
#endif
      end do

      icnt_s = 0
      do lrank = 0, nrank_chg-1
         if(sendinfo(1,lrank) == 0) cycle
         ista = sendinfo(1,lrank)
         iend = sendinfo(2,lrank)
         nel  = iend-ista+1
         icnt_s = icnt_s + 1
                                                 __TIMER_DO_START(854)
         do ri = 1, kimg
            do i = ista, iend
               sendbuf((i-ista+1+nel*(ri-1)),icnt_s) = ekin(i,ri,ispin)
            end do
         end do
                                                 __TIMER_DO_STOP(854)
#ifdef USE_NONBLK_COMM
         call mpi_isend(sendbuf(1,icnt_s), nel*kimg, mpi_double_precision, &
        &               lrank, itag, mpi_chg_world, req_s(icnt_s), ierr)
          if (ierr /= 0) then
             call mpi_abort(mpi_comm_world, 172, ierr)
          endif
#endif
      end do
                                                 __TIMER_DO_START(855)
      do ri = 1, kimg
         do i = ista_kngp, iend_kngp
            if (i < mymin) cycle
            if (mymax < i) cycle
            work(i,ri) = ekin(i,ri,ispin)
         end do
      end do
                                                 __TIMER_DO_STOP(855)
#ifdef USE_NONBLK_COMM
      call mpi_waitall(icnt_s, req_s, sta_s, ierr)
       if (ierr /= 0) then
          call mpi_abort(mpi_comm_world, 173, ierr)
       endif
      call mpi_waitall(icnt_r, req_r, sta_r, ierr)
       if (ierr /= 0) then
          call mpi_abort(mpi_comm_world, 174, ierr)
       endif
                                                 __TIMER_COMM_STOP(853)
#else
                                                 __TIMER_COMM_START_w_BARRIER(mpi_chg_world,694)
#ifndef USE_ALLTOALLV
! === DEBUG by tkato 2012/06/05 ================================================
!      if(sendmax/=0)then
! ==============================================================================
! === DEBUG by tkato 2012/06/04 ================================================
!      real(kind=DP), allocatable, dimension(:,:) :: sbuf, rbuf
! ==============================================================================
! === DEBUG by tkato 2012/06/04 ================================================
!      allocate(sbuf(sendmax*kimg,0:sendranks-1), stat=ierr)
!      allocate(rbuf(recvmax*kimg,0:recvranks-1), stat=ierr)
!      do i = 0, nrank_g - 1
!         sbuf(:,i)=sendbuf(:,i+1)
!      enddo
!      call MPI_ALLTOALL( sbuf, nel*kimg, mpi_double_precision, &
!     &                   rbuf, nel*kimg, mpi_double_precision, &
!     &                                    mpi_ke_world, ierr )
!      if (ierr /= 0) then
!         call mpi_abort(mpi_comm_world, 175, ierr)
!      endif
!      do i = 0, nrank_g - 1
!         recvbuf(:,i+1)=rbuf(:,i)
!      enddo
       maxbuf = max(sendmax,recvmax)
       call mpi_allreduce(MPI_IN_PLACE,maxbuf,1,mpi_integer,mpi_max,mpi_chg_world,ierr)
       allocate(sbuf(maxbuf*kimg,0:nrank_chg-1), stat=ierr)
       allocate(rbuf(maxbuf*kimg,0:nrank_chg-1), stat=ierr)
       icnt_s = 0
       do i = 0, nrank_chg - 1
          if(sendinfo(1,i) == 0) cycle
          icnt_s = icnt_s + 1
          sbuf(1:sendmax*kimg,i)=sendbuf(1:sendmax*kimg,icnt_s)
       enddo
       call MPI_ALLTOALL( sbuf, maxbuf*kimg, mpi_double_precision, &
      &                   rbuf, maxbuf*kimg, mpi_double_precision, &
      &                                    mpi_chg_world, ierr )
       if (ierr /= 0) then
          call mpi_abort(mpi_comm_world, 175, ierr)
       endif
       icnt_r = 0
       do i = 0, nrank_chg - 1
          if(recvinfo(1,i) == 0) cycle
          icnt_r = icnt_r + 1
          recvbuf(1:recvmax*kimg,icnt_r)=rbuf(1:recvmax*kimg,i)
       enddo
! ==============================================================================
       deallocate(sbuf)
       deallocate(rbuf)
! === DEBUG by tkato 2012/06/05 ================================================
!      endif
! ==============================================================================
#else
! === DEBUG by tkato 2012/06/05 ================================================
!      if(sendmax/=0)then
! ==============================================================================
! === DEBUG by tkato 2012/06/05 ================================================
!      integer, allocatable, dimension(:) :: sdsp, rdsp
!      integer, allocatable, dimension(:) :: scnt, rcnt
!      real(kind=DP), allocatable, dimension(:,:) :: sbuf, rbuf
!      allocate(sdsp(0:sendranks-1), stat=ierr)
!      allocate(rdsp(0:recvranks-1), stat=ierr)
!      allocate(scnt(0:sendranks-1), stat=ierr)
!      allocate(rcnt(0:recvranks-1), stat=ierr)
!      allocate(sbuf(sendmax*kimg,0:sendranks-1), stat=ierr)
!      allocate(rbuf(recvmax*kimg,0:recvranks-1), stat=ierr)
!      write(6,*) sendmax, recvmax, kimg
!      do i = 0, nrank_g - 1
!         sdsp(i)=sendmax*kimg*i
!         rdsp(i)=recvmax*kimg*i
!         scnt(i)=(sendinfo(2,i)-sendinfo(1,i)+1)*kimg
!         rcnt(i)=(recvinfo(2,i)-recvinfo(1,i)+1)*kimg
!         sbuf(:,i)=sendbuf(:,i+1)
!      enddo
!      call MPI_ALLTOALLV(      sendb, scnt, sdsp, &
!     &   mpi_double_precision, recvb, rcnt, rdsp, &
!     &   mpi_double_precision, mpi_ke_world, ierr )
!      if (ierr /= 0) then
!         call mpi_abort(mpi_comm_world, 175, ierr)
!      endif
!      do i = 0, nrank_g - 1
!         recvbuf(:,i+1)=rbuf(:,i)
!      enddo
       allocate(sdsp(0:nrank_chg-1), stat=ierr)
       allocate(rdsp(0:nrank_chg-1), stat=ierr)
       allocate(scnt(0:nrank_chg-1), stat=ierr)
       allocate(rcnt(0:nrank_chg-1), stat=ierr)
! === DEBUG by tkato 2012/06/05 ================================================
!      allocate(sbuf(sendmax*kimg,0:nrank_g-1), stat=ierr)
!      allocate(rbuf(recvmax*kimg,0:nrank_g-1), stat=ierr)
       if(sendmax /= 0) then
          allocate(sbuf(sendmax*kimg,0:nrank_chg-1), stat=ierr)
       else
          allocate(sbuf(1,0:nrank_chg-1), stat=ierr)
       endif
       if(recvmax /= 0) then
          allocate(rbuf(recvmax*kimg,0:nrank_chg-1), stat=ierr)
       else
          allocate(rbuf(1,0:nrank_chg-1), stat=ierr)
       endif
! ==============================================================================
       write(6,*) sendmax, recvmax, kimg
       icnt_s = 0
       do i = 0, nrank_chg - 1
          sdsp(i)=sendmax*kimg*i
          rdsp(i)=recvmax*kimg*i
          scnt(i)=(sendinfo(2,i)-sendinfo(1,i)+1)*kimg
          rcnt(i)=(recvinfo(2,i)-recvinfo(1,i)+1)*kimg
          if(sendinfo(1,i) == 0) cycle
          icnt_s = icnt_s + 1
          sbuf(:,i)=sendbuf(:,icnt_s)
       enddo
       call MPI_ALLTOALLV(      sbuf, scnt, sdsp, &
      &   mpi_double_precision, rbuf, rcnt, rdsp, &
      &   mpi_double_precision, mpi_chg_world, ierr )
       if (ierr /= 0) then
          call mpi_abort(mpi_comm_world, 175, ierr)
       endif
       icnt_r = 0
       do i = 0, nrank_chg - 1
          if(recvinfo(1,i) == 0) cycle
          icnt_r = icnt_r + 1
          recvbuf(:,icnt_r)=rbuf(:,i)
       enddo
! ==============================================================================
       deallocate(sdsp)
       deallocate(rdsp)
       deallocate(scnt)
       deallocate(rcnt)
       deallocate(sbuf)
       deallocate(rbuf)
! === DEBUG by tkato 2012/06/05 ================================================
!      endif
! ==============================================================================
#endif
                                                 __TIMER_COMM_STOP(694)
#endif
                                                 __TIMER_DO_START(856)
      do j = 1, icnt_r
         lrank = src(j)
         ista = recvinfo(1,lrank)
         iend = recvinfo(2,lrank)
         nel  = iend-ista+1
         do ri = 1, kimg
            do i = ista, iend
               work(i,ri) = recvbuf((i-ista+1+nel*(ri-1)),j)
            end do
         end do
      end do
                                                 __TIMER_DO_STOP(856)
      deallocate(sendbuf)
      deallocate(recvbuf)
                                                 __TIMER_SUB_STOP(732)
    end subroutine cp_chgq_by_ngpt
  end subroutine symmetrize_ekin_3D

  subroutine m_KE_set_modeled_ekin_density( ista_rho, iend_rho, &
       &                                    rho, grad_rho, grad2_rho, ekin_dens )
    integer,intent(in) :: ista_rho, iend_rho
    real(kind=DP), intent(in) :: rho( ista_rho:iend_rho, nspin )
    real(kind=DP), intent(in) :: grad_rho( ista_rho:iend_rho, nspin )
    real(kind=DP), intent(in) :: grad2_rho( ista_rho:iend_rho, nspin )
    real(kind=DP), intent(out) :: ekin_dens( ista_rho:iend_rho, nspin )

    select case ( ekin_density_type )

    case (1)
       call calc_Abramov_ekin_density( rho, grad_rho, grad2_rho, ekin_dens )

    case (2)                        ! Thomas-Fermi von Weizsacker
       call calc_TF_vWeiz_ekin_density( rho, grad_rho, grad2_rho, &
         &                              1.0d0, 1.0d0, ekin_dens )

    case (3)                        ! extended Thomas-Fermi
       call calc_Extended_TF_ekin_density( rho, grad_rho, grad2_rho, ekin_dens )

    end select

  contains

    subroutine calc_extended_TF_ekin_density( rho, grad_rho, grad2_rho, &
         &                                    ekin_dens )
      real(kind=DP), intent(in) :: rho( ista_rho:iend_rho, nspin )
      real(kind=DP), intent(in) :: grad_rho( ista_rho:iend_rho, nspin )
      real(kind=DP), intent(in) :: grad2_rho( ista_rho:iend_rho, nspin )
      real(kind=DP), intent(out) :: ekin_dens( ista_rho:iend_rho, nspin )

      integer :: is, i
      real(kind=DP) :: c1, c2, c3, kappa0, kappa2, kappa4, term1, term2, term3
      real(kind=DP), parameter :: delta1 = 1.0E-12

      c1 = (3.0d0*PAI**2)**(2.0d0/3.0d0)
      kappa0 = 3.0d0/5.0d0 *c1
      kappa2 = 1.0d0/36.0d0
      kappa4 = 1.0d0/6480.0d0/c1

      ekin_dens = 0.0d0

      Do is=1, nspin
         Do i=ista_rho, iend_rho
            c1 = rho(i,is)

            if ( c1 > delta1 ) then
               c2 =  grad_rho(i,is)/c1
               c3 = grad2_rho(i,is)/c1

               term1 = kappa0 *c1**(5.0d0/3.0d0)
               term2 = kappa2 *grad_rho(i,is)*c2
               term3 = kappa4 *c1**(1.0d0/3.0d0) &
                    &         *( 8.0d0 *c2**4 -27.0d0 *c2**2 *c3 +24.0d0 *c3**2 )

               ekin_dens(i,is) = term1 +term2 +term3
            endif
         End do
      End Do
      ekin_dens = ekin_dens /2.0d0

    end subroutine calc_extended_TF_ekin_density

    subroutine calc_TF_vWeiz_ekin_density( rho, grad_rho, grad2_rho, &
         &                                 weight_TF, weight_vWeiz, ekin_dens )
      real(kind=DP), intent(in) :: weight_TF,  weight_vWeiz
      real(kind=DP), intent(in) :: rho( ista_rho:iend_rho, nspin )
      real(kind=DP), intent(in) :: grad_rho( ista_rho:iend_rho, nspin )
      real(kind=DP), intent(in) :: grad2_rho( ista_rho:iend_rho, nspin )
      real(kind=DP), intent(out) :: ekin_dens( ista_rho:iend_rho, nspin )

      integer :: is, i
      real(kind=DP) :: c1, coeff, term1, term2, term3
      real(kind=DP), parameter :: delta1 = 1.0E-12

      coeff = 3.0d0/10.0d0 *(3.0D0*PAI**2)**( 2.0d0/3.0d0 )

      ekin_dens = 0.0d0

      Do is=1, nspin
         Do i=ista_rho, iend_rho
            c1 = rho(i,is)

            if ( c1 > delta1 ) then
               term1 = coeff *c1**(5.0d0/3.0d0)
               term2 = grad_rho(i,is)**2 /c1 /8.0d0
               ekin_dens(i,is) = weight_TF *term1 + weight_vWeiz *term2
            endif

         End do
      End Do
    end subroutine calc_TF_vWeiz_ekin_density

    subroutine calc_Abramov_ekin_density( rho, grad_rho, grad2_rho, ekin_dens )
      real(kind=DP), intent(in) :: rho( ista_rho:iend_rho, nspin )
      real(kind=DP), intent(in) :: grad_rho( ista_rho:iend_rho, nspin )
      real(kind=DP), intent(in) :: grad2_rho( ista_rho:iend_rho, nspin )
      real(kind=DP), intent(out) :: ekin_dens( ista_rho:iend_rho, nspin )

      integer :: is, i
      real(kind=DP) :: c1, c2, coeff, term1, term2, term3
      real(kind=DP), parameter :: delta1 = 1.0E-15

      coeff = 3.0d0/10.0d0 *(3.0D0*PAI**2)**( 2.0d0/3.0d0 )

      ekin_dens = 0.0d0

      Do is=1, nspin
         Do i=ista_rho, iend_rho
            c1 = rho(i,is) *nspin
            term1 = coeff *c1**(5.0d0/3.0d0) /nspin

            if ( c1 > delta1 ) then
               c2 = grad_rho(i,is) *nspin
               term2 = c2**2 /c1 /8.0 /nspin

               term2 = term2 /9.0d0
            else
               term2 = 0.0d0
            endif

            term3 = grad2_rho(i,is) /6.0d0
            term3 = term3 *1.05

            ekin_dens(i,is) = term1 + term2 + term3
         End do
      End Do
    end subroutine calc_Abramov_ekin_density

  end subroutine m_KE_set_modeled_ekin_density



  subroutine m_KE_add_ekindens_hard_rspace( ista_rho, iend_rho, ekin_dens )
    integer, intent(in) :: ista_rho, iend_rho
    real(kind=DP), intent(inout) :: ekin_dens( ista_rho,iend_rho, nspin )
    return
  end subroutine m_KE_add_ekindens_hard_rspace

  subroutine m_KE_wd_ekindens
    real(kind=DP), allocatable :: wk_mpi1(:,:,:), wk_mpi2(:,:,:)

    if(mype==0) rewind nf_ekindens
    if(npes > 1) then
       allocate(wk_mpi1(kgp,kimg,nspin)); wk_mpi1 = 0.d0
       allocate(wk_mpi2(kgp,kimg,nspin)); wk_mpi2 = 0.d0

       if ( use_symm_ekin_density ) then
          wk_mpi1(ista_kngp:iend_kngp,:,:) = ekins_l(ista_kngp:iend_kngp,:,:)
       else if ( use_asymm_ekin_density ) then
          wk_mpi1(ista_kngp:iend_kngp,:,:) = ekina_l(ista_kngp:iend_kngp,:,:)
       endif
       call mpi_allreduce( wk_mpi1, wk_mpi2, kgp*kimg*nspin, &
            &              mpi_double_precision, mpi_sum, mpi_chg_world, ierr)

       if(mype==0) write(nf_ekindens) wk_mpi2

       deallocate(wk_mpi1); deallocate(wk_mpi2)
    else
       if ( use_symm_ekin_density ) then
          write(nf_ekindens) ekins_l
       else if ( use_asymm_ekin_density ) then
          write(nf_ekindens) ekina_l
       end if
    endif

  end subroutine m_KE_wd_ekindens

  subroutine m_KE_rd_ekindens
    real(kind=DP), allocatable :: wk_mpi1(:,:,:)

    if(mype==0) rewind nf_ekindens
    if(npes > 1) then
       allocate(wk_mpi1(kgp,kimg,ndim_magmom)); wk_mpi1 = 0.d0

       if ( mype==0 ) read(nf_ekindens) wk_mpi1
       call mpi_bcast( wk_mpi1, kgp*kimg*ndim_magmom, mpi_double_precision, &
            &          0, MPI_CommGroup, ierr )

       if ( use_symm_ekin_density ) then
          ekins_l(ista_kngp:iend_kngp,:,:) = wk_mpi1(ista_kngp:iend_kngp,:,:)
       else if ( use_asymm_ekin_density ) then
          ekina_l(ista_kngp:iend_kngp,:,:) = wk_mpi1(ista_kngp:iend_kngp,:,:)
       endif
       deallocate(wk_mpi1)
    else
       if ( use_symm_ekin_density ) then
          read(nf_ekindens) ekins_l
       else if ( use_asymm_ekin_density ) then
          read(nf_ekindens) ekina_l
       end if
    endif
  end subroutine m_KE_rd_ekindens

end module m_KineticEnergy_Density
