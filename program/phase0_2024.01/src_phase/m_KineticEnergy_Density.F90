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
  use m_Electronic_Structure, only : zaj_l, occup_l, m_ES_WF_in_Rspace1
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
       &                            sw_add_ekin_hardpart_on_Gspace
  use m_Parallelization,      only : ista_fftph, iend_fftph, ista_fftp, iend_fftp, &
       &                             nel_fftp, idisp_fftp, ierr, &
       &                             ista_e, iend_e, istep_e
  use m_Ionic_System,           only : ntyp, ityp, natm, iwei, pos
  use m_Crystal_Structure,    only : nbztyp
  use m_PseudoPotential,    only :  kina_qitg_l, kins_qitg_l, modnrm, nqitg, &
       &                            m_PP_include_vanderbilt_pot, nlmt, dl2p, &
       &                            m_PP_find_maximum_l, &
       &                            m_PP_set_index_arrays1, m_PP_set_index_arrays2, &
       &                            mmesh, nmesh, xh, rmax, ilmt, ltp, mtp, taup, &
       &                            psirpw, phirpw
  use m_PlaneWaveBasisSet,    only : ylm_l,  m_pwBS_sphrp2, igfp_l
  use m_Charge_Density,       only : hsr, chgq_l

  use m_FFT,                  only : fft_box_size_CD, nfftp_nonpara, nfftp, &
       &                             m_FFT_CD_direct_c, m_FFT_alloc_CD_box, &
       &                             m_FFT_dealloc_CD_box, m_FFT_CD_inverse_c
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

  subroutine m_KE_calc_asymm_ekin_density
    ekina_l = 0.0d0

    call calc_local_part
    if ( sw_calc_ekin_density_hardpart == OFF ) return
    if ( sw_rspace_ekin_density == ON ) then
       if ( sw_add_ekin_hardpart_on_Gspace == OFF ) return
    endif

    if ( sw_rspace_ekin_density == ON ) then
       call add_hardpart_evaluated_rspace
    else
       call add_hardpart_to_ekinq_l( nfout, ndim_magmom, hsr, NO, .false. , &
            &                        kina_qitg_l, ekina_l )
    endif

  contains

    subroutine add_hardpart_evaluated_rspace
      integer :: is, ia, it, nma, ind, i
      integer :: lmt1, lmt2
      real(kind=DP) :: fac, fac0
      real(kind=DP), allocatable :: kina_tmp(:), ekina_hard(:,:)
      real(kind=DP), allocatable :: afft(:)

      call m_FFT_alloc_CD_box

      allocate( afft(ista_fftp:iend_fftp) ); afft = 0.0d0
      allocate( kina_tmp(nmesh_rs_aug_max)); kina_tmp = 0.0d0
      allocate( ekina_hard(ista_kngp:iend_kngp,kimg) ); ekina_hard = 0.0d0

      Do is=1, nspin
         afft = 0.0d0
         ekina_hard = 0.0d0

         Do ia=1, natm
            it = ityp(ia);   nma = nmesh_rs_aug(ia)
            kina_tmp = 0.0d0

            Do lmt1=1, ilmt(it)
               Do lmt2=1, ilmt(it)

                  fac0=1.0d0
                  fac = hsr(ia,lmt1,lmt2,is) *fac0
                  Do i=1, nma
                     kina_tmp(i) = kina_tmp(i) &
                          &            +fac *kina_qrs_rspace(i,ia,lmt1,lmt2)
                  End do
               End do
            End do

            Do i=1, nma
               ind = meshxyz_rs_aug(i,ia)
               if ( ind >= ista_fftph .and. ind <= iend_fftph ) then
                  afft(2*ind-1) = afft(2*ind-1) +kina_tmp(i)
               endif
            End do
         End do

         afft = afft / univol

         call m_FFT_CD_direct_c( nfout, afft )         ! afft(R) -> afft(G)
         call extract_matrix_from_FFTcd_mesh( afft, ekina_hard(:,:) )
!
         ekina_l(:,:,is) = ekina_l(:,:,is) +ekina_hard(:,:)
      End Do

      call m_FFT_dealloc_CD_box
      deallocate( kina_tmp );  deallocate( afft )
      deallocate( ekina_hard )

    end subroutine add_hardpart_evaluated_rspace

    subroutine calc_local_part
      integer :: ispin, ik, ib, ixyz, ierr
      integer :: i, ig, ri, iend, i1
      real(kind=DP) :: fac, occupation, c1, c2
      real(kind=DP), allocatable :: qxyz(:,:), psi_l(:,:,:,:)
      real(kind=DP), allocatable :: afft(:), bfft(:), cfft(:)

      allocate(afft(nfft)); afft =0.0d0;
      allocate(bfft(nfft)); bfft= 0.0d0
      allocate(cfft(nfft)); cfft= 0.0d0
      allocate( qxyz(kg1,3 ) ); qxyz = 0.0d0

      call m_FFT_alloc_WF_work()

      fac = 2.d0/(univol*kv3*product(fft_box_size_WF(1:3,1)))

      Do ispin=1, nspin
         afft = 0.0d0

         Do ik=ispin, kv3+ispin-nspin, nspin
            if ( map_k(ik) /= myrank_k ) cycle

            call k_plus_G_vectors_m( ik, kgp, kg1, kv3, iba, nbase, ngabc, rltv, &
                 &                   qxyz(:,1), qxyz(:,2), qxyz(:,3) )

            allocate( psi_l( kg1, np_e, ik:ik, kimg ) ); psi_l = 0.0d0

            Do ib = ista_e, iend_e, istep_e
               occupation = occup_l( map_z(ib),ik )
               if ( occupation < DELTA ) cycle

               call m_ES_WF_in_Rspace1( ista_k, iend_k, ik, ib, zaj_l, cfft )

               Do ig=1, iba(ik)
                  c1 = qxyz(ig,1)**2 + qxyz(ig,2)**2 + qxyz(ig,3)**2
                  psi_l(ig,map_z(ib),ik,1:2) = zaj_l(ig,map_z(ib),ik,1:2) *c1
               End Do
               call m_ES_WF_in_Rspace1( ik, ik, ik, ib, psi_l, bfft )

               Do i = 1, nfft-1, 2
                  c1 = bfft(i)*cfft(i)   +bfft(i+1)*cfft(i+1)
                  c2 = bfft(i)*cfft(i+1) -bfft(i+1)*cfft(i)
                  afft(i)   = afft(i)   + occupation *c1
                  afft(i+1) = afft(i+1) + occupation *c2
               End do

            End do
            deallocate( psi_l )
         End do
         !
         if ( npes > 1 ) then
            call mpi_allreduce( afft, bfft, nfft, mpi_double_precision, &
                 &              mpi_sum, MPI_CommGroup, ierr )
            afft = bfft
         end if
         call m_FFT_WF(ELECTRON,nfout,afft,DIRECT,OFF)
         !
         do ri = 1, kimg
            iend = iend_kngp
            if( iend_kngp > kg ) iend = kg
            if( ista_kngp <= iend ) then
               do i = ista_kngp, iend  !for mpi
                  i1 = kimg*igf(i) + (ri - kimg)
                  ekina_l(i,ri,ispin) = afft(i1)*fac
               end do
            endif
         end do

      End Do

      ekina_l = ekina_l /2.0d0

      call symmetrize_ekin( NO, ekina_l )

      deallocate( afft ); deallocate(bfft ); deallocate( cfft ); deallocate( qxyz )
      call m_FFT_dealloc_WF_work()

    end subroutine calc_local_part

  end subroutine m_KE_calc_asymm_ekin_density

  subroutine m_KE_calc_symm_ekin_density
    ekins_l = 0.0d0

    call calc_local_part

    if ( sw_calc_ekin_density_hardpart == OFF ) return

    if ( sw_rspace_ekin_density == ON ) then
       if ( sw_add_ekin_hardpart_on_Gspace == OFF ) return
    endif

    if ( sw_rspace_ekin_density == ON ) then
       call add_hardpart_evaluated_rspace
    else
       call add_hardpart_to_ekinq_l( nfout, ndim_magmom, hsr, NO, .false. , &
            &                        kins_qitg_l, ekins_l )
    endif

  contains

    subroutine add_hardpart_evaluated_rspace
      integer :: is, ia, it, nma, ind, i
      integer :: lmt1, lmt2
      real(kind=DP) :: fac, fac0
      real(kind=DP), allocatable :: kins_tmp(:), ekins_hard(:,:)
      real(kind=DP), allocatable :: afft(:)

      call m_FFT_alloc_CD_box

      allocate( afft(ista_fftp:iend_fftp) ); afft = 0.0d0
      allocate( kins_tmp(nmesh_rs_aug_max)); kins_tmp = 0.0d0
      allocate( ekins_hard(ista_kngp:iend_kngp,kimg) ); ekins_hard = 0.0d0

      Do is=1, nspin
         afft = 0.0d0
         ekins_hard = 0.0d0

         Do ia=1, natm
            it = ityp(ia);   nma = nmesh_rs_aug(ia)
            kins_tmp = 0.0d0

            Do lmt1=1, ilmt(it)
               Do lmt2=lmt1, ilmt(it)

                  fac0=2.0d0
                  if (lmt1.eq.lmt2) fac0=1.d0

                  fac = hsr(ia,lmt1,lmt2,is) *fac0
                  Do i=1, nma
                     kins_tmp(i) = kins_tmp(i) &
                          &            +fac *kins_qrs_rspace(i,ia,lmt1,lmt2)
                  End do
               End do
            End do

            Do i=1, nma
               ind = meshxyz_rs_aug(i,ia)
               if ( ind >= ista_fftph .and. ind <= iend_fftph ) then
                  afft(2*ind-1) = afft(2*ind-1) +kins_tmp(i)
               endif
            End do
         End do

         afft = afft / univol

         call m_FFT_CD_direct_c( nfout, afft )         ! afft(R) -> afft(G)
         call extract_matrix_from_FFTcd_mesh( afft, ekins_hard(:,:) )
!
         ekins_l(:,:,is) = ekins_l(:,:,is) +ekins_hard(:,:)
      End Do

      call m_FFT_dealloc_CD_box
      deallocate( kins_tmp );  deallocate( afft )
      deallocate( ekins_hard )

    end subroutine add_hardpart_evaluated_rspace

    subroutine calc_local_part
      integer :: ispin, ik, ib, ixyz, ierr
      integer :: i, ig, ri, iend, i1
      real(kind=DP) :: fac, occupation
      real(kind=DP), allocatable :: qxyz(:,:), psi_l(:,:,:,:)
      real(kind=DP), allocatable :: afft(:), bfft(:)

      allocate(afft(nfft)); afft = 0.0d0
      allocate(bfft(nfft)); bfft= 0.0d0
      allocate( qxyz(kg1,3 ) ); qxyz = 0.0d0

      call m_FFT_alloc_WF_work()

      fac = 2.d0/(univol*kv3*product(fft_box_size_WF(1:3,1)))

      Do ispin=1, nspin
         afft = 0.0d0

         Do ik=ispin, kv3+ispin-nspin, nspin
!!         Do ik=ispin, kv3, nspin
            if ( map_k(ik) /= myrank_k ) cycle

            call k_plus_G_vectors_m( ik, kgp, kg1, kv3, iba, nbase, ngabc, rltv, &
                 &                   qxyz(:,1), qxyz(:,2), qxyz(:,3) )

            allocate( psi_l( kg1, np_e, ik:ik, kimg ) ); psi_l = 0.0d0

            Do ib = ista_e, iend_e, istep_e
               occupation = occup_l( map_z(ib),ik )
               if ( occupation < DELTA ) cycle

               Do ixyz=1, 3
                  psi_l = 0.0d0
                  if ( kimg == 1 ) then
                     call phase_error_with_msg(nfout,"kimg=1 unsupported",__LINE__,__FILE__)
                  else
                     Do ig=1, iba(ik)
                        psi_l(ig,map_z(ib),ik,2) =  zaj_l(ig,map_z(ib),ik,1) &
                             &                       *qxyz(ig,ixyz)
                        psi_l(ig,map_z(ib),ik,1) = -zaj_l(ig,map_z(ib),ik,2) &
                             &                       *qxyz(ig,ixyz)
                     End do
                  endif

                  bfft = 0.0d0
                  call m_ES_WF_in_Rspace1( ik, ik, ik, ib, psi_l, bfft )

                  Do i = 1, nfft-1, 2
                     afft(i) = afft(i) + occupation*(bfft(i)**2+bfft(i+1)**2) ! MPI
                  End do

               End Do
            End do
            deallocate( psi_l )
         End do
       !
         if ( npes > 1 ) then
            call mpi_allreduce( afft, bfft, nfft, mpi_double_precision, &
                 &              mpi_sum, MPI_CommGroup, ierr )
            afft = bfft
         end if
         call m_FFT_WF(ELECTRON,nfout,afft,DIRECT,OFF)
         !
         do ri = 1, kimg
            iend = iend_kngp
            if( iend_kngp > kg ) iend = kg
            if( ista_kngp <= iend ) then
               do i = ista_kngp, iend  !for mpi
                  i1 = kimg*igf(i) + (ri - kimg)
                  ekins_l(i,ri,ispin) = afft(i1)*fac
               end do
            endif
         end do

      End Do

      ekins_l = ekins_l /2.0d0

      call symmetrize_ekin( NO, ekins_l )

      deallocate( afft ); deallocate(bfft ); deallocate( qxyz )
      call m_FFT_dealloc_WF_work()

    end subroutine calc_local_part

  end subroutine m_KE_calc_symm_ekin_density


  subroutine extract_matrix_from_FFTcd_mesh( afft, matrix )
    real(kind=DP), intent(in) :: afft(ista_fftp:iend_fftp)
    real(kind=DP), intent(out) :: matrix(ista_kngp:iend_kngp,kimg)

    integer :: ik, i, ip
    real(kind=DP) :: rinplw
    real(kind=DP), allocatable :: afft_mpi1(:)

    matrix = 0.0d0

    allocate(afft_mpi1(nfftp))
    if ( npes >1 ) then
       call mpi_allgatherv( afft, nel_fftp(mype), mpi_double_precision, &  ! MPI
            &               afft_mpi1, nel_fftp, idisp_fftp, mpi_double_precision, &
            &               MPI_CommGroup, ierr )
    else
       afft_mpi1 = afft
    end if

    rinplw = 1.d0 /product(fft_box_size_CD(1:3,1))
    do ik=1, kimg
       do i=ista_kngp, iend_kngp          !for mpi
          ip = (igfp_l(i)-1)*kimg + ik
          matrix(i,ik) = afft_mpi1(ip)*rinplw
       end do
    end do

  end subroutine extract_matrix_from_FFTcd_mesh

  subroutine k_plus_G_vectors_m(ik,kgp,kg1,knv3,iba,nbase,ngabc,rltv, &
       &                        qx,qy,qz)
    integer, intent(in)        :: ik, kgp,kg1,knv3,iba(knv3),nbase(kg1,knv3)
    integer, intent(in)        :: ngabc(kgp,3)
    real(kind=DP), intent(in)  :: rltv(3,3)
    real(kind=DP), intent(out) :: qx(kg1),qy(kg1),qz(kg1)

    integer :: i, ip
    real(kind=DP) :: ga, gb, gc

    do i = 1, iba(ik)
       ip = nbase(i,ik)
       ga = vkxyz(ik,1,BUCS) + real(ngabc(ip,1),kind=DP)
       gb = vkxyz(ik,2,BUCS) + real(ngabc(ip,2),kind=DP)
       gc = vkxyz(ik,3,BUCS) + real(ngabc(ip,3),kind=DP)
       qx(i)  = rltv(1,1)*ga + rltv(1,2)*gb + rltv(1,3)*gc
       qy(i)  = rltv(2,1)*ga + rltv(2,2)*gb + rltv(2,3)*gc
       qz(i)  = rltv(3,1)*ga + rltv(3,2)*gb + rltv(3,3)*gc
    end do
  end subroutine k_plus_G_vectors_m


  subroutine symmetrize_ekin(mode,ekin)
    integer,intent(in)           :: mode
    real(kind=DP), intent(inout) :: ekin(ista_kngp:iend_kngp,kimg,nspin)

    integer ::       ispin, ng, no, ngp, no1, no2
    real(kind=DP) :: fi, tx,ty,tz, fp, fc, fs, zcr, zci
    real(kind=DP), allocatable :: work(:,:), work2(:,:)

    allocate(work(kgp,kimg)); work = 0.d0
    allocate(work2(ista_kngp:iend_kngp,kimg)); work2 = 0.d0

    if(mode == ANTIFERRO) then
       fi = 1.d0/af
       no1 = nopr + 1; no2 = nopr + af
    else
       fi = 1.d0/nopr
       no1 = 1; no2 = nopr
    end if

    do ispin = 1, nspin, af+1
       work = 0.0d0
       call cp_ekin_to_work(ispin,ekin)     ! ekin -> work
       work2 = 0.d0                          ! initialization

       do no = no1, no2
          tx = tau(1,no,BUCS)*PAI2
          ty = tau(2,no,BUCS)*PAI2
          tz = tau(3,no,BUCS)*PAI2
          if(kimg == 1) then
             do ng = ista_kngp, iend_kngp !for mpi
                ngp = ngpt_l(ng,no)
                fp = ngabc(ngp,1)*tx + ngabc(ngp,2)*ty + ngabc(ngp,3)*tz
                work2(ng,1)        = work2(ng,1) + dcos(fp)*work(ngp,1)
             end do
          else if(kimg == 2) then
             do ng = ista_kngp, iend_kngp !for mpi
                ngp= ngpt_l(ng,no)
                fp = ngabc(ngp,1)*tx + ngabc(ngp,2)*ty + ngabc(ngp,3)*tz
                fc = dcos(fp);     fs = dsin(fp)
                zcr= work(ngp,1);  zci= work(ngp,kimg)
                work2(ng,1)        = work2(ng,1) + fc*zcr - fs*zci
                work2(ng,2)        = work2(ng,2) + fc*zci + fs*zcr
             end do
          end if
       end do
       if(mode /= ANTIFERRO) ekin(:,:,ispin) = work2(:,:)*fi
    end do

    if(mode == ANTIFERRO) ekin(:,:,nspin) = work2(:,:)*fi

    deallocate( work )
    deallocate(work2)

  contains

    subroutine cp_ekin_to_work(ispin,ekin)
      integer, intent(in) :: ispin
      real(DP),intent(in),dimension(ista_kngp:iend_kngp,kimg,nspin) :: ekin

      integer :: ng,ri, ierr
      real(kind=DP), allocatable, dimension(:,:) :: work_mpi

      do ri = 1, kimg
         do ng = ista_kngp, iend_kngp  !for mpi
            work(ng,ri) = ekin(ng,ri,ispin)
         end do
      end do

      if(npes >= 2) then
         allocate(work_mpi(kgp,kimg)); work_mpi = 0.d0
         call mpi_allreduce( work, work_mpi, kgp*kimg, mpi_double_precision, &
              &              mpi_sum, MPI_CommGroup, ierr )
         work = work_mpi
         deallocate(work_mpi)
      end if

    end subroutine cp_ekin_to_work

  end subroutine symmetrize_ekin


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

! =================================== modified by K. Tagami ============ 11.0
!!  subroutine add_hardpart_to_ekinq_l(nfout,kspin,hsr)
  subroutine add_hardpart_to_ekinq_l( nfout, kspin, hsr, singlemode, &
       &                              do_symmetrization, &
       &                              kin_qitg_l, ekinq_l )
! ===================================================================== 11.0
    !  The total operation number has been reduced not only for the gamma-point
    ! but also for other k-points by T. Yamasaki in April 2006.
    !  ----
    ! (Rev) T. Yamaskai, 31, Aug, 2007
    !     1. 'call set_index_arrays1' that included a bug is replaced
    !       by 'call m_PP_set_index_arrays1', whose bug is fixed.
    !     2. 'call set_index_arrays2' is also replaced by 'call
    !       m_PP_set_index_arrays2' that can be referred from other modules.
    !     3. contained subroutines, set_index_arrays1 and set_index_arrays2 were
    !       deleted.

    integer, intent(in) :: nfout, kspin, singlemode
    logical, intent(in) :: do_symmetrization

! ================================ modified by K. Tagami ================= 11.0
!    real(kind=DP), intent(in), dimension(natm,nlmt,nlmt,nspin):: hsr
    real(kind=DP), intent(in), dimension(natm,nlmt,nlmt,ndim_magmom):: hsr
! ======================================================================== 11.0

    real(kind=DP), intent(in) :: kin_qitg_l( ista_kngp:iend_kngp,nqitg )
    real(kind=DP), intent(out) :: ekinq_l(ista_kngp:iend_kngp,kimg,ndim_magmom)

    real(kind=DP), pointer, dimension(:)               :: ylm
    real(kind=DP), allocatable, target, dimension(:)   :: ylm_t
    real(kind=DP), allocatable, target, dimension(:,:) :: ylm_ext

    integer :: is,it,lmt1,lmt2,n,ia,mdvdb,il1,tau1,il2,tau2,ilm3,l3,iiqitg
    real(kind=DP) :: fac !, tpos(3)

    integer :: kngp_adj, n_ialist, n_ialist0, ia_start, ia_end, n_iagroup, n_ia, ia_g

    integer, allocatable, dimension(:)                :: il3

    real(kind=DP), allocatable, dimension(:,:) :: zfcos_x, zfsin_x
    real(kind=DP), allocatable, target, dimension(:) :: zfcos, zfsin
    real(kind=DP), allocatable, dimension(:) :: work(:,:)

    integer, allocatable, dimension(:) :: ia_list
#ifdef _VECTOR_TUNING_
    real(kind=DP), allocatable, dimension(:,:,:) :: shdg_x ! d(n_ialist0,maxm,nqitg)
#else
    real(kind=DP), allocatable, dimension(:,:) :: ylm_red, qitg_red
    real(kind=DP), allocatable, dimension(:) :: ylm_sum
    real(kind=DP), allocatable, dimension(:,:,:) :: ekinq_red
    real(kind=DP), allocatable, dimension(:,:) :: shdg  ! d(max,nqitg,nspin)
    integer ::          iqm, iqmmax
    real(kind=DP) :: zdga
#endif

    integer :: m, maxm, ip, np, iq, sw_spin
    integer, parameter :: mcritical = 4*2+1
    integer, allocatable, dimension(:) :: nqitg_sp, nqitg_sp0 !d(ntyp)
    integer, allocatable, dimension(:) :: iq2l3 ! d(nqitg)
    integer, allocatable, dimension(:,:) :: nc  ! d(maxm,nqitg)
    integer :: mc ! maxval(nc)
    integer, allocatable, dimension(:,:,:) :: nc2lmt1, nc2lmt2, nc2n ! d(mc,maxm,nqitg)
! NEC tune
    integer :: ibl1,ibl2,ibsize,ncache,iwidth

! =============================== added by K. Tagami ================ 11.0
    integer :: nspin_kt, ispi_start
! =================================================================== 11.0

! =============================== added by K. Tagami =============== 11.0
    if ( singlemode == YES ) then
       ispi_start = kspin
    else
       ispi_start = 1
    endif
! ================================================================== 11.0

!!!!    if(modnrm == EXECUT) then
       call m_PP_find_maximum_l(n)   !  n-1: maximum l
       n = (n-1) + (n-1) + 1
! =============================== Modified by K. Tagami =========
!       allocate(il3(n**2)); call substitute_il3(n**2,il3) ! -(b_Elec..)
       allocate(il3(n**2)); il3=0; call substitute_il3(n**2,il3) ! -(b_Elec..)
! ==============================================================

       allocate(nqitg_sp(ntyp)); allocate(nqitg_sp0(ntyp))
       allocate(iq2l3(nqitg))
       allocate(nc(mcritical,nqitg));nc=0
! =================================== Added by K. Tagami =======
	nqitg_sp = 0; nqitg_sp0 = 0; iq2l3 = 0
! ==============================================================
       call m_PP_set_index_arrays1(nfout,ntyp,nqitg,mcritical,n**2,il3 &
            & ,maxm,mc,nqitg_sp,nqitg_sp0,iq2l3,nc)
       allocate(nc2lmt1(mc,maxm,nqitg))
       allocate(nc2lmt2(mc,maxm,nqitg))
       allocate(nc2n(mc,maxm,nqitg))
! ==================================== Added by K. Tagami ======
       nc2lmt1 = 0;  nc2lmt2 = 0; nc2n = 0
! =============================================================
       call m_PP_set_index_arrays2(nfout,mc,maxm,nqitg,mcritical,n**2,il3,iq2l3 &
            & ,nc2lmt1,nc2lmt2,nc2n,nc) ! -> nc2lmt1, nc2lmt2, nc2n, nc

#ifdef _VECTOR_TUNING_

       n_ialist = 1
#ifdef HIUX
       n_ialist = 4
#endif
#ifdef VPP
       n_ialist = 8
#endif
#ifdef SX
       n_ialist = 8
#endif
       kngp_adj = iend_kngp - ista_kngp + 1
       if(mod(kngp_adj,2) == 0) kngp_adj = kngp_adj + 1
       allocate(zfcos_x(kngp_adj,n_ialist)); zfcos_x = 0.d0
       allocate(zfsin_x(kngp_adj,n_ialist)); zfsin_x = 0.d0
       allocate(ia_list(n_ialist)); ia_list = 0
! NEC tune
!       allocate(shdg_x(n_ialist))
       allocate(shdg_x(n_ialist,maxm,nqitg))
! ================================================ by K. Tagami =======
        shdg_x = 0.0d0
! ===================================================================

       if(n**2 > nel_Ylm) then
          allocate(ylm_ext(ista_kngp:iend_kngp,nel_Ylm+1:n**2)); ylm_ext = 0.d0
       end if
       allocate(ylm_t(ista_kngp:iend_kngp)); ylm_t = 0.d0
       do ilm3 = nel_Ylm+1, n**2
          call m_pwBS_sphrp2(ilm3,rltv,ista_kngp,iend_kngp,ylm_t)
          ylm_ext(:,ilm3) = ylm_t(:)
       end do

!!$       do is = 1, kspin, af+1
          do it = 1, ntyp
             mdvdb = m_PP_include_vanderbilt_pot(it)
!!!!!             if(mdvdb == SKIP) cycle

             n_ia = 0
             do ia = 1, natm
                if(ityp(ia) == it) n_ia = n_ia + 1
             end do

             if(n_ialist <=0) call phase_error_with_msg(nfout,'n_ialist is illegal <<m_Charge_Density.add_hardpart_to_ekinq_l>>'&
                                                       ,__LINE__,__FILE__)
             n_iagroup = n_ia/n_ialist + 1
             ia_start = 1

             do ia_g = 1, n_iagroup
                n_ialist0 = 0
                ia_list = 0
                AtomcountLoop: do ia = ia_start, natm
                   if(ityp(ia) == it) then
                      n_ialist0 = n_ialist0 + 1
                      ia_list(n_ialist0) = ia
                   end if
                   if(n_ialist0 >= n_ialist) exit AtomcountLoop
                end do AtomcountLoop
                ia_start = ia+1
                if(n_ialist0 >= 1 )then


! NEC tune ------------------------------------------------------------------->
!!$      ncache = (cachesize(3)*1024)*3/4
      ncache = (m_CtrlP_cachesize()*1024)*3/4
      if(ncache == 0) then
         ibsize = iend_kngp - ista_kngp + 1
      else

      iwidth = nqitg_sp(it) - nqitg_sp0(it)
      if(n_ialist0 == 1) then
         if(kimg == 1) then ! kin_qitg_l(i,iq)*ylm(iy)*zfsc1_1(iy)
           ibsize=ncache/(8*(2+iwidth))
         else ! kin_qitg_l(i,iq),ylm(iy),zfsc1_1(iy),zfsc2_1(iy)
           ibsize=ncache/(8*(3+iwidth))
         endif
      else if(n_ialist0 == 2) then
         if(kimg == 1) then ! kin_qitg_l(i,iq),ylm(iy),zfsc1_1(iy),zfsc1_2(iy)
           ibsize=ncache/(8*(3+iwidth))
         else ! kin_qitg_l(i,iq)
              ! ylm(iy),zfsc1_1(iy),zfsc1_2(iy),zfsc2_1(iy),zfsc2_2(iy)
           ibsize=ncache/(8*(5+iwidth))
         endif
      else if(n_ialist0 == 3) then
         if(kimg == 1) then ! kin_qitg_l(i,iq)
                            ! ylm(iy),zfsc1_1(iy),zfsc1_2(iy),zfsc1_3(iy)
           ibsize=ncache/(8*(4+iwidth))
         else ! kin_qitg_l(i,iq),ylm(iy),zfsc1_1(iy),zfsc1_2(iy),zfsc1_3(iy)
              !                      zfsc2_1(iy),zfsc2_2(iy),zfsc2_3(iy)
           ibsize=ncache/(8*(7+iwidth))
         endif
      else if(n_ialist0 == 4) then
         if(kimg == 1) then ! kin_qitg_l(i,iq)
                   ! ylm(iy),zfsc1_1(iy),zfsc1_2(iy),zfsc1_3(iy),zfsc1_4(iy)
           ibsize=ncache/(8*(5+iwidth))
         else ! kin_qitg_l(i,iq),ylm(iy)
              ! zfsc1_1(iy),zfsc1_2(iy),zfsc1_3(iy),zfsc1_4(iy)
              ! zfsc2_1(iy),zfsc2_2(iy),zfsc2_3(iy),zfsc2_4(iy)
           ibsize=ncache/(8*(9+iwidth))
         endif
      else if(n_ialist0 == 5) then
         if(kimg == 1) then ! kin_qitg_l(i,iq),ylm(iy)
               ! zfsc1_1(iy),zfsc1_2(iy),zfsc1_3(iy),zfsc1_4(iy),zfsc1_5(iy))
           ibsize=ncache/(8*(6+iwidth))
         else ! kin_qitg_l(i,iq),ylm(iy)
              ! zfsc1_1(iy),zfsc1_2(iy),zfsc1_3(iy),zfsc1_4(iy),zfsc1_5(iy)
              ! zfsc2_1(iy),zfsc2_2(iy),zfsc2_3(iy),zfsc2_4(iy),zfsc2_5(iy)
           ibsize=ncache/(8*(11+iwidth))
         endif
      else if(n_ialist0 == 6) then
         if(kimg == 1) then ! kin_qitg_l(i,iq),ylm(iy)
                 ! zfsc1_1(iy),zfsc1_2(iy),zfsc1_3(iy),zfsc1_4(iy)
                 ! zfsc1_5(iy),zfsc1_6(iy)
           ibsize=ncache/(8*(7+iwidth))
         else ! kin_qitg_l(i,iq),ylm(iy)
              ! zfsc1_1(iy),zfsc1_2(iy),zfsc1_3(iy),zfsc1_4(iy)
              ! zfsc1_5(iy),zfsc1_6(iy))
              ! zfsc2_1(iy),zfsc2_2(iy),zfsc2_3(iy),zfsc2_4(iy)
              ! zfsc2_5(iy),zfsc2_6(iy))
           ibsize=ncache/(8*(13+iwidth))
         endif
      else if(n_ialist0 == 7) then
         if(kimg == 1) then
             !  kin_qitg_l(i,iq),ylm(iy)
             !  zfsc1_1(iy),zfsc1_2(iy),zfsc1_3(iy),zfsc1_4(iy)
             !  zfsc1_5(iy),zfsc1_6(iy),zfsc1_7(iy)
           ibsize=ncache/(8*(8+iwidth))
         else
             !  kin_qitg_l(i,iq),ylm(iy)
             !  zfsc1_1(iy),zfsc1_2(iy),zfsc1_3(iy),zfsc1_4(iy)
             !  zfsc1_5(iy),zfsc1_6(iy),zfsc1_7(iy)
             !  zfsc2_1(iy),zfsc2_2(iy),zfsc2_3(iy),zfsc2_4(iy)
             !  zfsc2_5(iy),zfsc2_6(iy),zfsc2_7(iy)
           ibsize=ncache/(8*(15+iwidth))
         endif
      else if(n_ialist0 >= 8) then
         if(kimg == 1) then
             !  kin_qitg_l(i,iq),ylm(iy)
             !  zfsc1_1(iy),zfsc1_2(iy),zfsc1_3(iy),zfsc1_4(iy)
             !  zfsc1_5(iy),zfsc1_6(iy),zfsc1_7(iy),zfsc1_8(iy)
           ibsize=ncache/(8*(9+iwidth))
         else
            !   kin_qitg_l(i,iq),ylm(iy)
            !   zfsc1_1(iy),zfsc1_2(iy),zfsc1_3(iy),zfsc1_4(iy)
            !   zfsc1_5(iy),zfsc1_6(iy),zfsc1_7(iy),zfsc1_8(iy)
            !   zfsc2_1(iy),zfsc2_2(iy),zfsc2_3(iy),zfsc2_4(iy)
            !   zfsc2_5(iy),zfsc2_6(iy),zfsc2_7(iy),zfsc2_8(iy)
           ibsize=ncache/(8*(17+iwidth))
         endif
      end if
      endif
! debug
!write(6,990) 'n_ialist0,kimg,ibsize,ista_kngp,iend_kngp,iwidth=',&
!n_ialist0,kimg,ibsize,ista_kngp,iend_kngp,iwidth
!990 format(a,i2,i2,4i8)
      call calc_phase_b(natm,pos,ia_list,n_ialist0,kgp,ngabc,ista_kngp,iend_kngp,1,kngp_adj,zfcos_x,zfsin_x)
      do ibl1=ista_kngp,iend_kngp,ibsize
        ibl2=ibl1+ibsize-1
        if(ibl2.gt.iend_kngp) ibl2=iend_kngp
! NEC tune <-------------------------------------------------------------------

!!$          do ia = 1, natm
!!$             it = ityp(ia)
!!$             mdvdb = m_PP_include_vanderbilt_pot(it)
!!$             if(mdvdb == SKIP) cycle
! NEC tune (move 1 line to up)
!                   call calc_phase_b(natm,pos,ia_list,n_ialist0,kgp,ngabc,ista_kngp,iend_kngp,1,kngp_adj,zfcos_x,zfsin_x)
!!$             call calc_phase2(natm,pos,ia,kgp,ngabc,ista_kngp,iend_kngp,zfcos,zfsin)
                   ! -(b_Elec.)  -> zfcos, zfsin

! =========================== modified by K. Tagami ============= 11.0
!                do is = 1, kspin, af+1
                do is = ispi_start, kspin, af+1
! =============================================================== 11.0

! NEC tune --------------------------------------------------->
                   do iq = nqitg_sp0(it), nqitg_sp(it)
                      l3 = iq2l3(iq)
                      do m = 1, 2*l3+1
                         call sum_hsr_dot_gauntc(is,it,iq,m) ! hsr, dl2p -> shdg_x(n_ialist0)
                      end do
                   end do
! NEC tune <---------------------------------------------------

                   do iq = nqitg_sp0(it), nqitg_sp(it)
                      l3 = iq2l3(iq)
                      do m = 1, 2*l3+1
                         ilm3 = l3*l3+m
! NEC tune
!                         call sum_hsr_dot_gauntc(is,it,iq,m) ! hsr, dl2p -> shdg_x(n_ialist0)
                         if(ilm3 <= nel_Ylm) then
                            ylm => ylm_l(ista_kngp:iend_kngp,ilm3)
                         else
                            ylm => ylm_ext(ista_kngp:iend_kngp,ilm3)
                         end if
                         call add_hardpart_to_ekinq_l_core4(iq) ! iq, shdg_x, exp(-iGR), kin_qitg_l, ylm -> ekinq_l
                         !    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                      end do
                   end do
! NEC tune
                end do
             end do ! ibl1
                end if
             end do! ia_g
          end do! it
!!$       end do! is
       deallocate(ylm_t)
       if(allocated(ylm_ext)) deallocate(ylm_ext)
       deallocate(shdg_x)
       deallocate(ia_list,zfsin_x,zfcos_x,il3)
#else

       ncache = (m_CtrlP_cachesize()*1024)*3/4
       if(ncache == 0) then
          ibsize = iend_kngp - ista_kngp + 1
       else
          iwidth = 0.d0
          do it = 1, ntyp
             iwidth = max(iwidth,nqitg_sp(it)-nqitg_sp0(it)+1)
          end do
          if(kimg == 1) then ! kin_qitg_l(i,iq)*ylm(iy)*zfsc1_1(iy)
             ibsize=ncache/(8*(iwidth + 1))
          else ! kin_qitg_l(i,iq),ylm(iy),zfsc1_1(iy),zfsc2_1(iy)
             ibsize=ncache/(8*(iwidth + 2))
          endif
       end if

       allocate(zfcos(ibsize)); zfcos = 0.d0
       allocate(zfsin(ibsize)); zfsin = 0.d0
       allocate(qitg_red(ibsize,nqitg))
       allocate(ylm_red(ibsize,n**2))
       allocate(ylm_sum(ibsize))

! =========================== modified by K. Tagami =================== 11.0
!       if(kspin == 2 .and. af == 0) then
!          sw_spin = ON
!       else
!          sw_spin = OFF
!       end if
! ===================================================================== 11.0

! ============================== added by K. Tagami ================== 11.0
       if ( noncol ) then
          nspin_kt = kspin
       else
          nspin_kt = nspin /(af+1)
       endif
! ================================================================ 11.0

       iqmmax = 0
       do it = 1, ntyp
          iqm = 0
          do iq = nqitg_sp0(it), nqitg_sp(it)
             l3 = iq2l3(iq)
             do m = 1, 2*l3+1
                iqm = iqm+1
             end do
          end do
          if(iqmmax < iqm) iqmmax = iqm
       end do

! ===================== modified by K. Tagami ============== 11.0
!       if(sw_spin == ON) then
!          allocate(ekinq_red(ibsize,kimg,2))
!          allocate(shdg(iqmmax,2))
!       else if(sw_spin == OFF) then
!          allocate(ekinq_red(ibsize,kimg,1))
!          allocate(shdg(iqmmax,1))
!       end if
! ========================================================== 11.0
! ========================= added by K. Tagami =========== 11.0
       allocate( ekinq_red( ibsize, kimg, nspin_kt ))
       allocate( shdg( iqmmax,nspin_kt ))
! ========================================================= 11.0

       do ibl1=ista_kngp,iend_kngp,ibsize
          ibl2=ibl1+ibsize-1
          if(ibl2.gt.iend_kngp) ibl2=iend_kngp
          if(ibl2.gt.kgp) ibl2 = kgp

          ekinq_red = 0.d0
          call substitute_qitgred()  ! kin_qitg_l -> qitg_red
          call substitute_ylmred() ! ylm_l, ylm_ext -> ylm_red
          do ia = 1, natm
             it = ityp(ia)
             mdvdb = m_PP_include_vanderbilt_pot(it)
!!!!             if(mdvdb == SKIP) cycle

             call calc_phase_div(ia) ! -> zfsin, zfcos

!!$             iqm = 0
!!$             do iq = nqitg_sp0(it), nqitg_sp(it)
!!$                l3 = iq2l3(iq)
!!$                do m = 1, 2*l3+1
!!$                   iqm = iqm+1
!!$                   call sum_hsr_dot_gauntc0(it,ia,iq,m,iqm) ! hsr, dl2p -> shdg
!!$                end do
!!$             end do
!!$             iqmmax = iqm

! ================================= modified by K. Tagami ========== 11.0
!             do is = 1, kspin, af+1
             do is = ispi_start, kspin, af+1
! ================================================================== 11.0

                iqm = 0
                do iq = nqitg_sp0(it), nqitg_sp(it)
                   l3 = iq2l3(iq)
                   ylm_sum = 0.d0
                   do m = 1, 2*l3+1
                      ilm3 = l3*l3+m
                      iqm = iqm+1
                      call sum_hsr_dot_gauntc0(it,ia,iq,m,iqm) ! hsr, dl2p -> shdg
                      ylm_sum(:) = ylm_sum(:) + shdg(iqm,is)*ylm_red(:,ilm3)
                   end do
                   if(mod(l3,2) == 0) then
                      zdga = real(zi**(-l3))
                      call add_hardpart_to_ekinq_l_div0(zdga,iq) ! iq, shdg_x, exp(-iGR), kin_qitg_l, ylm -> ekinq_l
                   else
                      zdga = aimag(zi**(-l3))
                      call add_hardpart_to_ekinq_l_div1(zdga,iq) ! iq, shdg_x, exp(-iGR), kin_qitg_l, ylm -> ekinq_l
                   end if
!!$                   call add_hardpart_to_ekinq_l_div(iq) ! iq, shdg_x, exp(-iGR), kin_qitg_l, ylm -> ekinq_l
                   !    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                end do
             end do
          end do
          call cp_ekinqred2ekinq()
       end do
       deallocate(ylm_sum,ylm_red,qitg_red,shdg)
       deallocate(ekinq_red)
       deallocate(zfsin,zfcos)
#endif

       deallocate(nc2n,nc2lmt2,nc2lmt1,nc,iq2l3,nqitg_sp,nqitg_sp0)
!!!    endif

    if(nbztyp >= SIMPLE_CUBIC .or. af /=0 ) then
       allocate(work(kgp,kimg))
! =========================================== Added by K. Tagami ====
        work = 0.0d0
! ===================================================================

! ==================================== modified by K. Tagami ============== 11.0
      if(nbztyp >= SIMPLE_CUBIC .and. do_symmetrization ) then
         if ( noncol ) then
!!!            call symmetrize_ekin_noncl( NO,ekinq_l )
         else
            call symmetrize_ekin( NO,ekinq_l )
         endif
      endif
! ========================================================================== 11.0

       deallocate(work)
    end if

  contains

#ifdef _VECTOR_TUNING_
    subroutine sum_hsr_dot_gauntc(is,it,iq,m)
      integer, intent(in) :: is,it,iq,m
      integer :: ip, lmt1, lmt2, np
      real(kind=DP) :: fac
! NEC tune
!      shdg_x(1:n_ialist0) = 0.d0
      shdg_x(1:n_ialist0,m,iq) = 0.d0
      do ip = 1, nc(m,iq)
         lmt1 = nc2lmt1(ip,m,iq)
         lmt2 = nc2lmt2(ip,m,iq)
         np = nc2n(ip,m,iq)
         fac = 2.d0; if(lmt1 == lmt2) fac = 1.d0
         do ia = 1, n_ialist0
! NEC tune
!            shdg_x(ia) = shdg_x(ia) + &
            shdg_x(ia,m,iq) = shdg_x(ia,m,iq) + &
                 & fac*iwei(ia_list(ia))*hsr(ia_list(ia),lmt1,lmt2,is)*dl2p(lmt1,lmt2,np,it)
         end do
      end do
    end subroutine sum_hsr_dot_gauntc
    subroutine add_hardpart_to_ekinq_l_core4(iq)
      integer, intent(in) :: iq
      integer       :: i,iy, ia
!!$      real(kind=DP) :: dga, flekinq, f,f2
      real(kind=DP) :: flekinq, f,f2, flekinq_1, flekinq_2, flekinq_3, flekinq_4 &
           &                            , flekinq_5, flekinq_6, flekinq_7, flekinq_8, zdga
      real(kind=DP) :: qf1,qf2, qy
!!$      real(kind=DP), pointer, dimension(:) :: zfsc1,zfsc2
      real(kind=DP), pointer, dimension(:) :: zfsc1_1,zfsc2_1, zfsc1_2, zfsc2_2 &
           &         , zfsc1_3,zfsc2_3,zfsc1_4,zfsc2_4, zfsc1_5,zfsc2_5,zfsc1_6,zfsc2_6 &
           &         , zfsc1_7,zfsc2_7,zfsc1_8,zfsc2_8
!!$      real(kind=DP), allocatable, dimension(:) :: w_f ! d(ista_kngp:iend_kngp)
!!$      if(n_ialist == 1) allocate(w_f(ista_kngp:iend_kngp))
      if(mod(l3,2) == 0) then
         if(n_ialist0 >= 1) zdga = real(zi**(-l3))
         if(kimg == 1) then
!!$            zfsc1 => zfcos
            if(n_ialist0 >= 1) zfsc1_1 => zfcos_x(:,1)
            if(n_ialist0 >= 2) zfsc1_2 => zfcos_x(:,2)
            if(n_ialist0 >= 3) zfsc1_3 => zfcos_x(:,3)
            if(n_ialist0 >= 4) zfsc1_4 => zfcos_x(:,4)
            if(n_ialist0 >= 5) zfsc1_5 => zfcos_x(:,5)
            if(n_ialist0 >= 6) zfsc1_6 => zfcos_x(:,6)
            if(n_ialist0 >= 7) zfsc1_7 => zfcos_x(:,7)
            if(n_ialist0 >= 8) zfsc1_8 => zfcos_x(:,8)
         else
!!$            f2 = -1; zfsc1 => zfcos; zfsc2 => zfsin
            f2 = -1
            if(n_ialist0 >= 1) then
               zfsc1_1 => zfcos_x(:,1); zfsc2_1 => zfsin_x(:,1)
            end if
            if(n_ialist0 >= 2) then
               zfsc1_2 => zfcos_x(:,2); zfsc2_2 => zfsin_x(:,2)
            end if
            if(n_ialist0 >= 3) then
               zfsc1_3 => zfcos_x(:,3); zfsc2_3 => zfsin_x(:,3)
            end if
            if(n_ialist0 >= 4) then
               zfsc1_4 => zfcos_x(:,4); zfsc2_4 => zfsin_x(:,4)
            end if
            if(n_ialist0 >= 5) then
               zfsc1_5 => zfcos_x(:,5); zfsc2_5 => zfsin_x(:,5)
            end if
            if(n_ialist0 >= 6) then
               zfsc1_6 => zfcos_x(:,6); zfsc2_6 => zfsin_x(:,6)
            end if
            if(n_ialist0 >= 7) then
               zfsc1_7 => zfcos_x(:,7); zfsc2_7 => zfsin_x(:,7)
            end if
            if(n_ialist0 >= 8) then
               zfsc1_8 => zfcos_x(:,8); zfsc2_8 => zfsin_x(:,8)
            end if
         end if
      else
!!$         flekinq = fac*aimag(zi**(-l3))*dga*hsr(ia,lmt1,lmt2,is)
         if(n_ialist0 >= 1) zdga = aimag(zi**(-l3))
         if(kimg == 1) then
!!$            zfsc1 => zfsin
            if(n_ialist0 >= 1) zfsc1_1 => zfsin_x(:,1)
            if(n_ialist0 >= 2) zfsc1_2 => zfsin_x(:,2)
            if(n_ialist0 >= 3) zfsc1_3 => zfsin_x(:,3)
            if(n_ialist0 >= 4) zfsc1_4 => zfsin_x(:,4)
            if(n_ialist0 >= 5) zfsc1_5 => zfsin_x(:,5)
            if(n_ialist0 >= 6) zfsc1_6 => zfsin_x(:,6)
            if(n_ialist0 >= 7) zfsc1_7 => zfsin_x(:,7)
            if(n_ialist0 >= 8) zfsc1_8 => zfsin_x(:,8)
         else
!!$            f2 = 1; zfsc1 => zfsin; zfsc2 => zfcos
            f2 = 1
            if(n_ialist0 >= 1) then
               zfsc1_1 => zfsin_x(:,1); zfsc2_1 => zfcos_x(:,1)
            end if
            if(n_ialist0 >= 2) then
               zfsc1_2 => zfsin_x(:,2); zfsc2_2 => zfcos_x(:,2)
            end if
            if(n_ialist0 >= 3) then
               zfsc1_3 => zfsin_x(:,3); zfsc2_3 => zfcos_x(:,3)
            end if
            if(n_ialist0 >= 4) then
               zfsc1_4 => zfsin_x(:,4); zfsc2_4 => zfcos_x(:,4)
            end if
            if(n_ialist0 >= 5) then
               zfsc1_5 => zfsin_x(:,5); zfsc2_5 => zfcos_x(:,5)
            end if
            if(n_ialist0 >= 6) then
               zfsc1_6 => zfsin_x(:,6); zfsc2_6 => zfcos_x(:,6)
            end if
            if(n_ialist0 >= 7) then
               zfsc1_7 => zfsin_x(:,7); zfsc2_7 => zfcos_x(:,7)
            end if
            if(n_ialist0 >= 8) then
               zfsc1_8 => zfsin_x(:,8); zfsc2_8 => zfcos_x(:,8)
            end if
         end if
      end if
! NEC tune ---------------------------------------->
!      if(n_ialist0 >= 1) flekinq_1 = zdga*shdg_x(1)
!      if(n_ialist0 >= 2) flekinq_2 = zdga*shdg_x(2)
!      if(n_ialist0 >= 3) flekinq_3 = zdga*shdg_x(3)
!      if(n_ialist0 >= 4) flekinq_4 = zdga*shdg_x(4)
!      if(n_ialist0 >= 5) flekinq_5 = zdga*shdg_x(5)
!      if(n_ialist0 >= 6) flekinq_6 = zdga*shdg_x(6)
!      if(n_ialist0 >= 7) flekinq_7 = zdga*shdg_x(7)
!      if(n_ialist0 >= 8) flekinq_8 = zdga*shdg_x(8)
! NEC tune -----------------------------------------
      if(n_ialist0 >= 1) flekinq_1 = zdga*shdg_x(1,m,iq)
      if(n_ialist0 >= 2) flekinq_2 = zdga*shdg_x(2,m,iq)
      if(n_ialist0 >= 3) flekinq_3 = zdga*shdg_x(3,m,iq)
      if(n_ialist0 >= 4) flekinq_4 = zdga*shdg_x(4,m,iq)
      if(n_ialist0 >= 5) flekinq_5 = zdga*shdg_x(5,m,iq)
      if(n_ialist0 >= 6) flekinq_6 = zdga*shdg_x(6,m,iq)
      if(n_ialist0 >= 7) flekinq_7 = zdga*shdg_x(7,m,iq)
      if(n_ialist0 >= 8) flekinq_8 = zdga*shdg_x(8,m,iq)
! NEC tune <----------------------------------------

      if(n_ialist0 == 1) then
         if(kimg == 1) then
#ifdef HIUX
*POPTION PARALLEL
#endif
#ifdef NEC_TUNE_HARD
!CDIR PARALLEL DO PRIVATE(i,iy)
#endif
! NEC tune
!            do i = ista_kngp, iend_kngp     !for mpi
            do i = ibl1, ibl2
               iy = i - ista_kngp+1
               ekinq_l(i,1,is) = ekinq_l(i,1,is) &
                    & +flekinq_1*kin_qitg_l(i,iq)*ylm(iy)*zfsc1_1(iy)
            end do
         else
#ifdef HIUX
*POPTION PARALLEL
#endif
#ifdef NEC_TUNE_HARD
!CDIR PARALLEL DO PRIVATE(i,iy)
#endif
! NEC tune
!            do i = ista_kngp, iend_kngp   !for mpi
            do i = ibl1, ibl2
               iy = i - ista_kngp+1
               f = flekinq_1*kin_qitg_l(i,iq)*ylm(iy)
               ekinq_l(i,1,is) = ekinq_l(i,1,is) + f * zfsc1_1(iy)
               ekinq_l(i,2,is) = ekinq_l(i,2,is) + f2 * f * zfsc2_1(iy)
            end do
         end if
      else if(n_ialist0 == 2) then
         if(kimg == 1) then
#ifdef HIUX
*POPTION PARALLEL
#endif
#ifdef NEC_TUNE_HARD
!CDIR PARALLEL DO PRIVATE(i,iy)
#endif
! NEC tune
!            do i = ista_kngp, iend_kngp     !for mpi
            do i = ibl1, ibl2
               iy = i - ista_kngp+1
               ekinq_l(i,1,is) = ekinq_l(i,1,is) + kin_qitg_l(i,iq)*ylm(iy) &
                    &         *( flekinq_1*zfsc1_1(iy) + flekinq_2*zfsc1_2(iy) )
            end do
         else
#ifdef HIUX
*POPTION PARALLEL
#endif
#ifdef NEC_TUNE_HARD
!CDIR PARALLEL DO PRIVATE(i,iy)
#endif
! NEC tune
!            do i = ista_kngp, iend_kngp   !for mpi
            do i = ibl1, ibl2
               iy = i - ista_kngp+1
               ekinq_l(i,1,is) = ekinq_l(i,1,is) + kin_qitg_l(i,iq)*ylm(iy) &
                    &         *( flekinq_1*zfsc1_1(iy) + flekinq_2*zfsc1_2(iy) )
               ekinq_l(i,2,is) = ekinq_l(i,2,is) + f2 * kin_qitg_l(i,iq)*ylm(iy) &
                    &         *( flekinq_1*zfsc2_1(iy) + flekinq_2*zfsc2_2(iy) )
            end do
         end if
      else if(n_ialist0 == 3) then
         if(kimg == 1) then
#ifdef HIUX
*POPTION PARALLEL
#endif
#ifdef NEC_TUNE_HARD
!CDIR PARALLEL DO PRIVATE(i,iy)
#endif
! NEC tune
!            do i = ista_kngp, iend_kngp     !for mpi
            do i = ibl1, ibl2
               iy = i - ista_kngp+1
               ekinq_l(i,1,is) = ekinq_l(i,1,is) + kin_qitg_l(i,iq)*ylm(iy) &
                    &     *( flekinq_1*zfsc1_1(iy) + flekinq_2*zfsc1_2(iy) + flekinq_3*zfsc1_3(iy) )
            end do
         else
#ifdef HIUX
*POPTION PARALLEL
#endif
#ifdef NEC_TUNE_HARD
!CDIR PARALLEL DO PRIVATE(i,iy)
#endif
! NEC tune
!            do i = ista_kngp, iend_kngp   !for mpi
            do i = ibl1, ibl2
               iy = i - ista_kngp+1
               ekinq_l(i,1,is) = ekinq_l(i,1,is) + kin_qitg_l(i,iq)*ylm(iy) &
                    &     *( flekinq_1*zfsc1_1(iy) + flekinq_2*zfsc1_2(iy) + flekinq_3*zfsc1_3(iy))
               ekinq_l(i,2,is) = ekinq_l(i,2,is) + f2 * kin_qitg_l(i,iq)*ylm(iy) &
                    &     *( flekinq_1*zfsc2_1(iy) + flekinq_2*zfsc2_2(iy) + flekinq_3*zfsc2_3(iy))
            end do
         end if
      else if(n_ialist0 == 4) then
         if(kimg == 1) then
#ifdef HIUX
*POPTION PARALLEL
#endif
#ifdef NEC_TUNE_HARD
!CDIR PARALLEL DO PRIVATE(i,iy)
#endif
! NEC tune
!            do i = ista_kngp, iend_kngp     !for mpi
            do i = ibl1, ibl2
               iy = i - ista_kngp+1
               ekinq_l(i,1,is) = ekinq_l(i,1,is) + kin_qitg_l(i,iq)*ylm(iy) &
                    &   *(flekinq_1*zfsc1_1(iy) + flekinq_2*zfsc1_2(iy) + flekinq_3*zfsc1_3(iy)+ flekinq_4*zfsc1_4(iy))
            end do
         else
#ifdef HIUX
*POPTION PARALLEL
#endif
#ifdef NEC_TUNE_HARD
!CDIR PARALLEL DO PRIVATE(i,iy)
#endif
! NEC tune
!            do i = ista_kngp, iend_kngp   !for mpi
            do i = ibl1, ibl2
               iy = i - ista_kngp+1
               ekinq_l(i,1,is) = ekinq_l(i,1,is) + kin_qitg_l(i,iq)*ylm(iy) &
                    &   *( flekinq_1*zfsc1_1(iy) + flekinq_2*zfsc1_2(iy) + flekinq_3*zfsc1_3(iy) + flekinq_4*zfsc1_4(iy))
               ekinq_l(i,2,is) = ekinq_l(i,2,is) + f2 * kin_qitg_l(i,iq)*ylm(iy) &
                    &   *( flekinq_1*zfsc2_1(iy) + flekinq_2*zfsc2_2(iy) + flekinq_3*zfsc2_3(iy) + flekinq_4*zfsc2_4(iy))
            end do
         end if
      else if(n_ialist0 == 5) then
         if(kimg == 1) then
#ifdef HIUX
*POPTION PARALLEL
#endif
#ifdef NEC_TUNE_HARD
!CDIR PARALLEL DO PRIVATE(i,iy)
#endif
! NEC tune
!            do i = ista_kngp, iend_kngp     !for mpi
            do i = ibl1, ibl2
               iy = i - ista_kngp+1
               ekinq_l(i,1,is) = ekinq_l(i,1,is) + kin_qitg_l(i,iq)*ylm(iy) &
                    &   *(flekinq_1*zfsc1_1(iy) + flekinq_2*zfsc1_2(iy) + flekinq_3*zfsc1_3(iy)+ flekinq_4*zfsc1_4(iy) &
                    &    +flekinq_5*zfsc1_5(iy))
            end do
         else
#ifdef HIUX
*POPTION PARALLEL
#endif
#ifdef NEC_TUNE_HARD
!CDIR PARALLEL DO PRIVATE(i,iy)
#endif
! NEC tune
!            do i = ista_kngp, iend_kngp   !for mpi
            do i = ibl1, ibl2
               iy = i - ista_kngp+1
               ekinq_l(i,1,is) = ekinq_l(i,1,is) + kin_qitg_l(i,iq)*ylm(iy) &
                    &   *( flekinq_1*zfsc1_1(iy) + flekinq_2*zfsc1_2(iy) + flekinq_3*zfsc1_3(iy) + flekinq_4*zfsc1_4(iy)&
                    &     +flekinq_5*zfsc1_5(iy))
               ekinq_l(i,2,is) = ekinq_l(i,2,is) + f2 * kin_qitg_l(i,iq)*ylm(iy) &
                    &   *( flekinq_1*zfsc2_1(iy) + flekinq_2*zfsc2_2(iy) + flekinq_3*zfsc2_3(iy) + flekinq_4*zfsc2_4(iy)&
                    &     +flekinq_5*zfsc2_5(iy))
            end do
         end if
      else if(n_ialist0 == 6) then
         if(kimg == 1) then
#ifdef HIUX
*POPTION PARALLEL
#endif
#ifdef NEC_TUNE_HARD
!CDIR PARALLEL DO PRIVATE(i,iy)
#endif
! NEC tune
!            do i = ista_kngp, iend_kngp     !for mpi
            do i = ibl1, ibl2
               iy = i - ista_kngp+1
               ekinq_l(i,1,is) = ekinq_l(i,1,is) + kin_qitg_l(i,iq)*ylm(iy) &
                    &   *(flekinq_1*zfsc1_1(iy) + flekinq_2*zfsc1_2(iy) + flekinq_3*zfsc1_3(iy)+ flekinq_4*zfsc1_4(iy) &
                    &    +flekinq_5*zfsc1_5(iy) + flekinq_6*zfsc1_6(iy))
            end do
         else
#ifdef HIUX
*POPTION PARALLEL
#endif
#ifdef NEC_TUNE_HARD
!CDIR PARALLEL DO PRIVATE(i,iy)
#endif
! NEC tune
!            do i = ista_kngp, iend_kngp   !for mpi
            do i = ibl1, ibl2
               iy = i - ista_kngp+1
               ekinq_l(i,1,is) = ekinq_l(i,1,is) + kin_qitg_l(i,iq)*ylm(iy) &
                    &   *( flekinq_1*zfsc1_1(iy) + flekinq_2*zfsc1_2(iy) + flekinq_3*zfsc1_3(iy) + flekinq_4*zfsc1_4(iy)&
                    &     +flekinq_5*zfsc1_5(iy) + flekinq_6*zfsc1_6(iy))
               ekinq_l(i,2,is) = ekinq_l(i,2,is) + f2 * kin_qitg_l(i,iq)*ylm(iy) &
                    &   *( flekinq_1*zfsc2_1(iy) + flekinq_2*zfsc2_2(iy) + flekinq_3*zfsc2_3(iy) + flekinq_4*zfsc2_4(iy)&
                    &     +flekinq_5*zfsc2_5(iy) + flekinq_6*zfsc2_6(iy))
            end do
         end if
      else if(n_ialist0 == 7) then
         if(kimg == 1) then
#ifdef HIUX
*POPTION PARALLEL
#endif
#ifdef NEC_TUNE_HARD
!CDIR PARALLEL DO PRIVATE(i,iy)
#endif
! NEC tune
!            do i = ista_kngp, iend_kngp     !for mpi
            do i = ibl1, ibl2
               iy = i - ista_kngp+1
               ekinq_l(i,1,is) = ekinq_l(i,1,is) + kin_qitg_l(i,iq)*ylm(iy) &
                    &   *(flekinq_1*zfsc1_1(iy) + flekinq_2*zfsc1_2(iy) + flekinq_3*zfsc1_3(iy) + flekinq_4*zfsc1_4(iy) &
                    &    +flekinq_5*zfsc1_5(iy) + flekinq_6*zfsc1_6(iy) + flekinq_7*zfsc1_7(iy))
            end do
         else
#ifdef HIUX
*POPTION PARALLEL
#endif
#ifdef NEC_TUNE_HARD
!CDIR PARALLEL DO PRIVATE(i,iy)
#endif
! NEC tune
!            do i = ista_kngp, iend_kngp   !for mpi
            do i = ibl1, ibl2
               iy = i - ista_kngp+1
               ekinq_l(i,1,is) = ekinq_l(i,1,is) + kin_qitg_l(i,iq)*ylm(iy) &
                    &   *( flekinq_1*zfsc1_1(iy) + flekinq_2*zfsc1_2(iy) + flekinq_3*zfsc1_3(iy) + flekinq_4*zfsc1_4(iy)&
                    &     +flekinq_5*zfsc1_5(iy) + flekinq_6*zfsc1_6(iy) + flekinq_7*zfsc1_7(iy))
               ekinq_l(i,2,is) = ekinq_l(i,2,is) + f2 * kin_qitg_l(i,iq)*ylm(iy) &
                    &   *( flekinq_1*zfsc2_1(iy) + flekinq_2*zfsc2_2(iy) + flekinq_3*zfsc2_3(iy) + flekinq_4*zfsc2_4(iy)&
                    &     +flekinq_5*zfsc2_5(iy) + flekinq_6*zfsc2_6(iy) + flekinq_7*zfsc2_7(iy))
            end do
         end if
      else if(n_ialist0 >= 8) then
         if(kimg == 1) then
#ifdef HIUX
*POPTION PARALLEL
#endif
#ifdef NEC_TUNE_HARD
!CDIR PARALLEL DO PRIVATE(i,iy)
#endif
! NEC tune
!            do i = ista_kngp, iend_kngp     !for mpi
            do i = ibl1, ibl2
               iy = i - ista_kngp+1
               ekinq_l(i,1,is) = ekinq_l(i,1,is) + kin_qitg_l(i,iq)*ylm(iy) &
                    &   *(flekinq_1*zfsc1_1(iy) + flekinq_2*zfsc1_2(iy) + flekinq_3*zfsc1_3(iy) + flekinq_4*zfsc1_4(iy) &
                    &    +flekinq_5*zfsc1_5(iy) + flekinq_6*zfsc1_6(iy) + flekinq_7*zfsc1_7(iy) + flekinq_8*zfsc1_8(iy))
            end do
         else
#ifdef HIUX
*POPTION PARALLEL
#endif
#ifdef NEC_TUNE_HARD
!CDIR PARALLEL DO PRIVATE(i,iy)
#endif
! NEC tune
!            do i = ista_kngp, iend_kngp   !for mpi
            do i = ibl1, ibl2
               iy = i - ista_kngp+1
               ekinq_l(i,1,is) = ekinq_l(i,1,is) + kin_qitg_l(i,iq)*ylm(iy) &
                    &   *( flekinq_1*zfsc1_1(iy) + flekinq_2*zfsc1_2(iy) + flekinq_3*zfsc1_3(iy) + flekinq_4*zfsc1_4(iy)&
                    &     +flekinq_5*zfsc1_5(iy) + flekinq_6*zfsc1_6(iy) + flekinq_7*zfsc1_7(iy) + flekinq_8*zfsc1_8(iy))
               ekinq_l(i,2,is) = ekinq_l(i,2,is) + f2 * kin_qitg_l(i,iq)*ylm(iy) &
                    &   *( flekinq_1*zfsc2_1(iy) + flekinq_2*zfsc2_2(iy) + flekinq_3*zfsc2_3(iy) + flekinq_4*zfsc2_4(iy)&
                    &     +flekinq_5*zfsc2_5(iy) + flekinq_6*zfsc2_6(iy) + flekinq_7*zfsc2_7(iy) + flekinq_8*zfsc2_8(iy))
            end do
         end if
      end if
!!$      if(n_ialist == 1) deallocate(w_f)
    end subroutine add_hardpart_to_ekinq_l_core4
#else
    subroutine substitute_qitgred()
      integer :: iq, i
      do iq = 1, nqitg
         do i = 1, ibl2-ibl1+1
            qitg_red(i, iq) = kin_qitg_l(i+ibl1-1,iq)
         end do
      end do
    end subroutine substitute_qitgred

    subroutine substitute_ylmred
      integer :: ilm, i

      do ilm = 1, nel_Ylm
         do i = 1, ibl2-ibl1+1
            ylm_red(i,ilm) = ylm_l(i+ibl1-1,ilm)
         end do
      end do

      if(n**2 > nel_Ylm) then
         allocate(ylm_t(ibl1:ibl2)); ylm_t = 0.d0
         do ilm = nel_ylm+1, n**2
            call m_pwBS_sphrp2(ilm,rltv,ibl1,ibl2,ylm_t)
            do i = 1, ibl2-ibl1+1
               ylm_red(i,ilm) = ylm_t(i+ibl1-1)
            end do
         end do
         deallocate(ylm_t)
      end if
    end subroutine substitute_ylmred

    subroutine sum_hsr_dot_gauntc0(it,ia,iq,m,iqm)
      integer, intent(in) :: it,ia,iq,m,iqm
      integer :: ip, lmt1, lmt2, np
      real(kind=DP) :: fac

! ========================= modified by K. Tagami ====================== 11.0
!      if(sw_spin == OFF) then
!         shdg(iqm,1) = 0.d0
!         do ip = 1, nc(m,iq)
!            lmt1 = nc2lmt1(ip,m,iq)
!            lmt2 = nc2lmt2(ip,m,iq)
!            np = nc2n(ip,m,iq)
!            fac = 2.d0; if(lmt1 == lmt2) fac = 1.d0
!            shdg(iqm,1) = shdg(iqm,1) + fac*iwei(ia)*hsr(ia,lmt1,lmt2,1)*dl2p(lmt1,lmt2,np,it)
!         end do
!      else if(sw_spin == ON) then
!         shdg(iqm,1) = 0.d0; shdg(iqm,2) = 0.d0
!         do ip = 1, nc(m,iq)
!            lmt1 = nc2lmt1(ip,m,iq)
!            lmt2 = nc2lmt2(ip,m,iq)
!            np = nc2n(ip,m,iq)
!            fac = 2.d0; if(lmt1 == lmt2) fac = 1.d0
!            shdg(iqm,1) = shdg(iqm,1) + fac*iwei(ia)*hsr(ia,lmt1,lmt2,1) &
!                           &           *dl2p(lmt1,lmt2,np,it)
!            shdg(iqm,2) = shdg(iqm,2) + fac*iwei(ia)*hsr(ia,lmt1,lmt2,2) &
!                           &          *dl2p(lmt1,lmt2,np,it)
!         end do
!      end if

      shdg(iqm,:) = 0.d0

      do ip = 1, nc(m,iq)
         lmt1 = nc2lmt1(ip,m,iq);    lmt2 = nc2lmt2(ip,m,iq)
         np = nc2n(ip,m,iq)
         fac = 2.d0; if(lmt1 == lmt2) fac = 1.d0
         shdg(iqm,:) = shdg(iqm,:) &
              &  + fac *iwei(ia) *hsr(ia,lmt1,lmt2,:) *dl2p(lmt1,lmt2,np,it)
      end do
! ============================================================== 11.0

    end subroutine sum_hsr_dot_gauntc0

    subroutine calc_phase_div(ia)
      integer, intent(in) :: ia
      real(kind=DP) :: fx, fy, fz, ph
      integer :: i, iy
      fx = pos(ia,1)*PAI2
      fy = pos(ia,2)*PAI2
      fz = pos(ia,3)*PAI2
      do i = 1, ibl2-ibl1+1
         iy = i + ibl1 - 1
         ph = ngabc(iy,1)*fx+ngabc(iy,2)*fy+ngabc(iy,3)*fz
         zfcos(i) = dcos(ph)
         zfsin(i) = dsin(ph)
      end do
    end subroutine calc_phase_div

    subroutine add_hardpart_to_ekinq_l_div(iq)
      integer, intent(in) :: iq
      integer       :: i,iy, iy2
      real(kind=DP) :: flekinq, f,f2,  zdga, flekinq_up, flekinq_dw, f_up, f_dw
      real(kind=DP), pointer, dimension(:) :: zfsc1,zfsc2

      if(mod(l3,2) == 0) then
         zdga = real(zi**(-l3))
         if(kimg == 1) then
            zfsc1 => zfcos
         else
            f2 = -1
            zfsc1 => zfcos; zfsc2 => zfsin
         end if
      else
         zdga = aimag(zi**(-l3))
         if(kimg == 1) then
            zfsc1 => zfsin
         else
            f2 = 1
            zfsc1 => zfsin; zfsc2 => zfcos
         end if
      end if

!!$      if(sw_spin == ON) then
!!$         flekinq_up = zdga*shdg(iqm,1)
!!$         flekinq_dw = zdga*shdg(iqm,2)
!!$      else if(sw_spin == OFF) then
!!$         flekinq = zdga*shdg(iqm,1)
!!$      end if

      if(kimg == 1) then
!!$         if(sw_spin == ON) then
!!$            do i = 1, ibl2-ibl1+1
!!$               f = qitg_red(i,iq)*ylm_sum(i)*zfsc1(i)
!!$               ekinq_red(i,1,1) = ekinq_red(i,1,1) +zdga*f
!!$               ekinq_red(i,1,2) = ekinq_red(i,1,2) +zdga*f
!!$            end do
!!$         else if(sw_spin == OFF) then
         do i = 1, ibl2-ibl1+1
            ekinq_red(i,1,is) = ekinq_red(i,1,is) +zdga*qitg_red(i,iq)*ylm_sum(i)*zfsc1(i)
         end do
!!$         end if
      else
!!$         if(sw_spin == ON) then
!!$            do i = 1, ibl2-ibl1+1
!!$               f = qitg_red(i,iq)*ylm_red(i,ilm3)
!!$               ekinq_red(i,1,1) = ekinq_red(i,1,1) +      flekinq_up*f * zfsc1(i)
!!$               ekinq_red(i,2,1) = ekinq_red(i,2,1) + f2 * flekinq_up*f * zfsc2(i)
!!$               ekinq_red(i,1,2) = ekinq_red(i,1,2) +      flekinq_dw*f * zfsc1(i)
!!$               ekinq_red(i,2,2) = ekinq_red(i,2,2) + f2 * flekinq_dw*f * zfsc2(i)
!!$            end do
!!$         else if(sw_spin == OFF) then
         do i = 1, ibl2-ibl1+1
            f = zdga*qitg_red(i,iq)*ylm_sum(i)
            ekinq_red(i,1,is) = ekinq_red(i,1,is) +      f * zfsc1(i)
            ekinq_red(i,2,is) = ekinq_red(i,2,is) + f2 * f * zfsc2(i)
         end do
!!$         end if
      end if
    end subroutine add_hardpart_to_ekinq_l_div

    subroutine add_hardpart_to_ekinq_l_div0(zdga,iq)
      real(kind=DP), intent(in) :: zdga
      integer, intent(in) :: iq
      integer       :: i
      real(kind=DP) :: f

      if(kimg == 1) then
         do i = 1, ibl2-ibl1+1
            ekinq_red(i,1,is) = ekinq_red(i,1,is) +zdga*qitg_red(i,iq)*ylm_sum(i)*zfcos(i)
         end do
      else
         do i = 1, ibl2-ibl1+1
            f = zdga*qitg_red(i,iq)*ylm_sum(i)
            ekinq_red(i,1,is) = ekinq_red(i,1,is) + f * zfcos(i)
            ekinq_red(i,2,is) = ekinq_red(i,2,is) - f * zfsin(i)
         end do
      end if
    end subroutine add_hardpart_to_ekinq_l_div0

    subroutine add_hardpart_to_ekinq_l_div1(zdga,iq)
      real(kind=DP), intent(in) :: zdga
      integer, intent(in) :: iq
      integer       :: i
      real(kind=DP) :: f

      if(kimg == 1) then
         do i = 1, ibl2-ibl1+1
            ekinq_red(i,1,is) = ekinq_red(i,1,is) +zdga*qitg_red(i,iq)*ylm_sum(i)*zfsin(i)
         end do
      else
         do i = 1, ibl2-ibl1+1
            f = zdga*qitg_red(i,iq)*ylm_sum(i)
            ekinq_red(i,1,is) = ekinq_red(i,1,is) +  f * zfsin(i)
            ekinq_red(i,2,is) = ekinq_red(i,2,is) +  f * zfcos(i)
         end do
      end if
    end subroutine add_hardpart_to_ekinq_l_div1

    subroutine cp_ekinqred2ekinq()
      integer :: i, ir

! =============================== modified by K. Tagami ============= 11.0
!      if(sw_spin == ON) then
!         do ir = 1, kimg
!            do i = ibl1, ibl2
!               ekinq_l(i,ir,1) = ekinq_l(i,ir,1) + ekinq_red(i-ibl1+1,ir,1)
!               ekinq_l(i,ir,2) = ekinq_l(i,ir,2) + ekinq_red(i-ibl1+1,ir,2)
!            end do
!         end do
!      else if(sw_spin == OFF) then
!         do ir = 1, kimg
!            do i = ibl1, ibl2
!               ekinq_l(i,ir,1) = ekinq_l(i,ir,1) + ekinq_red(i-ibl1+1,ir,1)
!            end do
!         end do
!      end if

      do ir = 1, kimg
         do i = ibl1, ibl2
            ekinq_l(i,ir,:) = ekinq_l(i,ir,:) + ekinq_red(i-ibl1+1,ir,:)
         end do
      end do
! ==================================================================== 11.0

    end subroutine cp_ekinqred2ekinq
#endif
  end subroutine add_hardpart_to_ekinq_l

  subroutine m_KE_calc_kins_qrs_rspace
    integer :: imesh
    integer :: ia, it, nma
    integer :: il1, il2, im1, im2, tau1, tau2, lmt1, lmt2
    integer :: nspher1, nspher2
    integer :: ier, ir

    real(kind=DP) :: ps1, ps2, dps1, dps2
    real(kind=DP) :: ph1, ph2, dph1, dph2
    real(kind=DP) :: cosx, cosy, cosz
    real(kind=DP) :: ctmp_psi1(3), ctmp_psi2(3)
    real(kind=DP) :: ctmp_phi1(3), ctmp_phi2(3)
    real(kind=DP) :: csum_psi, csum_phi
    real(kind=DP) :: hn, inl, inm, inn

    real(kind=DP), allocatable :: cx(:), cy(:), cz(:), rdiff(:)
    real(kind=DP), allocatable :: ylm1(:), ylm2(:), dylm1(:,:), dylm2(:,:)
    real(kind=DP), allocatable :: psi1(:), psi2(:), wk_psi1(:), wk_psi2(:)
    real(kind=DP), allocatable :: phi1(:), phi2(:), wk_phi1(:), wk_phi2(:)
    real(kind=DP), allocatable :: radr(:)

    inl = 1.d0/dble(fft_box_size_CD(1,1))
    inm = 1.d0/dble(fft_box_size_CD(2,1))
    inn = 1.d0/dble(fft_box_size_CD(3,1))

    allocate(cx(nmesh_rs_aug_max)); cx=0.d0
    allocate(cy(nmesh_rs_aug_max)); cy=0.d0
    allocate(cz(nmesh_rs_aug_max)); cz=0.d0
    allocate(rdiff(nmesh_rs_aug_max)); rdiff=0.d0

    allocate(ylm1(nmesh_rs_aug_max)); ylm1 = 0.d0
    allocate(ylm2(nmesh_rs_aug_max)); ylm2 = 0.d0
    allocate(dylm1(nmesh_rs_aug_max,3)); dylm1 = 0.d0
    allocate(dylm2(nmesh_rs_aug_max,3)); dylm2 = 0.d0

    allocate( psi1(mmesh) ); psi1 = 0.0d0
    allocate( psi2(mmesh) ); psi2 = 0.0d0
    allocate( phi1(mmesh) ); phi1 = 0.0d0
    allocate( phi2(mmesh) ); phi2 = 0.0d0
    allocate( wk_psi1(mmesh) ); wk_psi1 = 0.0d0
    allocate( wk_psi2(mmesh) ); wk_psi2 = 0.0d0
    allocate( wk_phi1(mmesh) ); wk_phi1 = 0.0d0
    allocate( wk_phi2(mmesh) ); wk_phi2 = 0.0d0
    allocate( radr(mmesh) ); radr = 0.0d0

!    write(*,*) 'nmesh_rs_aug_max = ', nmesh_rs_aug_max
    if ( nmesh_rs_aug_max == 0 ) return

    allocate( kins_qrs_rspace( nmesh_rs_aug_max,natm,nlmt,nlmt) )
    kins_qrs_rspace = 0.0d0

    Do ia=1, natm
       it = ityp(ia)
       nma = nmesh_rs_aug(ia)

       if ( nma == 0 ) cycle

       call rmeshs(nmesh(it),nmesh(it),xh(it),rmax(it),radr,hn) ! -(b_PP)
       call m_RS_R_minus_pos( pos, ia, nmesh_rs_aug_max, inl, inm, inn, &
            &                 cx, cy, cz, rdiff, &
            &                 meshx_rs_aug, meshy_rs_aug, meshz_rs_aug )

       Do lmt1 = 1, ilmt(it)
          il1 = ltp(lmt1,it); im1 = mtp(lmt1,it);  tau1 = taup(lmt1,it)
          nspher1 = ( il1 -1 )**2 +im1

          psi1(1:nmesh(it)) = psirpw(1:nmesh(it),il1,tau1,it) /radr(1:nmesh(it))
          phi1(1:nmesh(it)) = phirpw(1:nmesh(it),il1,tau1,it) /radr(1:nmesh(it))

          call sphr( nma, nspher1, cx(1:nma), cy(1:nma), cz(1:nma), ylm1(1:nma) )
          call sphr_diff( nma, nma, nspher1, cx(1:nma), cy(1:nma), cz(1:nma), &
               &          dylm1(1:nma,1:3) )
          call init_cubic_spline( nmesh(it), radr(1:nmesh(it)),&
               &                  psi1(1:nmesh(it)), wk_psi1(1:nmesh(it) ) )
          call init_cubic_spline( nmesh(it), radr(1:nmesh(it)),&
               &                  phi1(1:nmesh(it)), wk_phi1(1:nmesh(it) ) )

          Do lmt2=lmt1, ilmt(it)
             il2 = ltp(lmt2,it); im2 = mtp(lmt2,it);  tau2 = taup(lmt2,it)
             nspher2 = ( il2 -1 )**2 +im2

             psi2(1:nmesh(it)) = psirpw(1:nmesh(it),il2,tau2,it) /radr(1:nmesh(it))
             phi2(1:nmesh(it)) = phirpw(1:nmesh(it),il2,tau2,it) /radr(1:nmesh(it))

             call sphr( nma, nspher2, cx(1:nma), cy(1:nma), cz(1:nma), ylm2(1:nma) )
             call sphr_diff( nma, nma, nspher2, cx(1:nma), cy(1:nma), cz(1:nma), &
                  &          dylm2(1:nma,1:3) )
             call init_cubic_spline( nmesh(it), radr(1:nmesh(it)),&
                  &                  psi2(1:nmesh(it)), wk_psi2(1:nmesh(it) ) )
             call init_cubic_spline( nmesh(it), radr(1:nmesh(it)),&
                  &                  phi2(1:nmesh(it)), wk_phi2(1:nmesh(it) ) )

             do imesh=1,nma
                call cubic_spline( nmesh(it), radr, psi1(1:nmesh(it)),&
                     &             wk_psi1(1:nmesh(it)), rdiff(imesh), ps1, dps1 )
                call cubic_spline( nmesh(it), radr, psi2(1:nmesh(it)),&
                     &             wk_psi2(1:nmesh(it)), rdiff(imesh), ps2, dps2 )
                call cubic_spline( nmesh(it), radr, phi1(1:nmesh(it)),&
                     &             wk_phi1(1:nmesh(it)), rdiff(imesh), ph1, dph1 )
                call cubic_spline( nmesh(it), radr, phi2(1:nmesh(it)),&
                     &             wk_phi2(1:nmesh(it)), rdiff(imesh), ph2, dph2 )

                if (rdiff(imesh) < radr(2) ) cycle

                cosx = cx(imesh)/rdiff(imesh)
                cosy = cy(imesh)/rdiff(imesh)
                cosz = cz(imesh)/rdiff(imesh)

                ctmp_psi1(1) = cosx *dps1 *ylm1(imesh) +ps1 *dylm1(imesh,1)
                ctmp_psi1(2) = cosy *dps1 *ylm1(imesh) +ps1 *dylm1(imesh,2)
                ctmp_psi1(3) = cosz *dps1 *ylm1(imesh) +ps1 *dylm1(imesh,3)

                ctmp_psi2(1) = cosx *dps2 *ylm2(imesh) +ps2 *dylm2(imesh,1)
                ctmp_psi2(2) = cosy *dps2 *ylm2(imesh) +ps2 *dylm2(imesh,2)
                ctmp_psi2(3) = cosz *dps2 *ylm2(imesh) +ps2 *dylm2(imesh,3)

                ctmp_phi1(1) = cosx *dph1 *ylm1(imesh) +ph1 *dylm1(imesh,1)
                ctmp_phi1(2) = cosy *dph1 *ylm1(imesh) +ph1 *dylm1(imesh,2)
                ctmp_phi1(3) = cosz *dph1 *ylm1(imesh) +ph1 *dylm1(imesh,3)

                ctmp_phi2(1) = cosx *dph2 *ylm2(imesh) +ph2 *dylm2(imesh,1)
                ctmp_phi2(2) = cosy *dph2 *ylm2(imesh) +ph2 *dylm2(imesh,2)
                ctmp_phi2(3) = cosz *dph2 *ylm2(imesh) +ph2 *dylm2(imesh,3)

                csum_psi = ctmp_psi1(1) *ctmp_psi2(1) &
                     &      + ctmp_psi1(2) *ctmp_psi2(2) &
                     &      + ctmp_psi1(3) *ctmp_psi2(3)
                csum_phi = ctmp_phi1(1) *ctmp_phi2(1) &
                     &      + ctmp_phi1(2) *ctmp_phi2(2) &
                     &      + ctmp_phi1(3) *ctmp_phi2(3)
                kins_qrs_rspace(imesh,ia,lmt1,lmt2) = ( csum_psi -csum_phi )/2.0d0
             End do
          End do
       End Do
    End Do

    deallocate( cx ); deallocate( cy ); deallocate( cz ); deallocate( rdiff )
    deallocate( ylm1 ); deallocate( ylm2 ); deallocate( dylm1 ); deallocate( dylm2 )
    deallocate( psi1 ); deallocate( psi2 ); deallocate( wk_psi1 ); deallocate( wk_psi2 )
    deallocate( phi1 ); deallocate( phi2 ); deallocate( wk_phi1 ); deallocate( wk_phi2 )
    deallocate( radr )

  end subroutine m_KE_calc_kins_qrs_rspace

  subroutine m_KE_calc_kina_qrs_rspace
    integer :: imesh
    integer :: ia, it, nma
    integer :: il1, il2, im1, im2, tau1, tau2, lmt1, lmt2
    integer :: nspher1, nspher2
    integer :: ier, ir

    real(kind=DP) :: ps1, ps2, dps1, dps2, ddps2
    real(kind=DP) :: ph1, ph2, dph1, dph2, ddph2
    real(kind=DP) :: csum_psi, csum_phi
    real(kind=DP) :: hn, inl, inm, inn

    real(kind=DP), allocatable :: cx(:), cy(:), cz(:), rdiff(:)
    real(kind=DP), allocatable :: ylm1(:), ylm2(:), dylm1(:,:), dylm2(:,:)
    real(kind=DP), allocatable :: psi1(:), psi2(:), wk_psi1(:), wk_psi2(:)
    real(kind=DP), allocatable :: phi1(:), phi2(:), wk_phi1(:), wk_phi2(:)
    real(kind=DP), allocatable :: dpsi2(:), dphi2(:), wk_dpsi2(:), wk_dphi2(:)
    real(kind=DP), allocatable :: radr(:)

    inl = 1.d0/dble(fft_box_size_CD(1,1))
    inm = 1.d0/dble(fft_box_size_CD(2,1))
    inn = 1.d0/dble(fft_box_size_CD(3,1))

    allocate(cx(nmesh_rs_aug_max)); cx=0.d0
    allocate(cy(nmesh_rs_aug_max)); cy=0.d0
    allocate(cz(nmesh_rs_aug_max)); cz=0.d0
    allocate(rdiff(nmesh_rs_aug_max)); rdiff=0.d0

    allocate(ylm1(nmesh_rs_aug_max)); ylm1 = 0.d0
    allocate(ylm2(nmesh_rs_aug_max)); ylm2 = 0.d0
    allocate(dylm1(nmesh_rs_aug_max,3)); dylm1 = 0.d0
    allocate(dylm2(nmesh_rs_aug_max,3)); dylm2 = 0.d0

    allocate( psi1(mmesh) ); psi1 = 0.0d0
    allocate( psi2(mmesh) ); psi2 = 0.0d0
    allocate( phi1(mmesh) ); phi1 = 0.0d0
    allocate( phi2(mmesh) ); phi2 = 0.0d0
    allocate( dpsi2(mmesh) ); dpsi2 = 0.0d0
    allocate( dphi2(mmesh) ); dphi2 = 0.0d0

    allocate( wk_psi1(mmesh) ); wk_psi1 = 0.0d0
    allocate( wk_psi2(mmesh) ); wk_psi2 = 0.0d0
    allocate( wk_phi1(mmesh) ); wk_phi1 = 0.0d0
    allocate( wk_phi2(mmesh) ); wk_phi2 = 0.0d0
    allocate( wk_dpsi2(mmesh) ); wk_dpsi2 = 0.0d0
    allocate( wk_dphi2(mmesh) ); wk_dphi2 = 0.0d0

    allocate( radr(mmesh) ); radr = 0.0d0

    allocate( kina_qrs_rspace( nmesh_rs_aug_max,natm,nlmt,nlmt) )
    kina_qrs_rspace = 0.0d0

    Do ia=1, natm
       it = ityp(ia)
       nma = nmesh_rs_aug(ia)

       call rmeshs(nmesh(it),nmesh(it),xh(it),rmax(it),radr,hn) ! -(b_PP)
       call m_RS_R_minus_pos( pos, ia, nmesh_rs_aug_max, inl, inm, inn, &
            &                 cx, cy, cz, rdiff, &
            &                 meshx_rs_aug, meshy_rs_aug, meshz_rs_aug )

       Do lmt1 = 1, ilmt(it)
          il1 = ltp(lmt1,it); im1 = mtp(lmt1,it);  tau1 = taup(lmt1,it)
          nspher1 = ( il1 -1 )**2 +im1

          psi1(1:nmesh(it)) = psirpw(1:nmesh(it),il1,tau1,it) /radr(1:nmesh(it))
          phi1(1:nmesh(it)) = phirpw(1:nmesh(it),il1,tau1,it) /radr(1:nmesh(it))

          call sphr( nma, nspher1, cx(1:nma), cy(1:nma), cz(1:nma), ylm1(1:nma) )

          call init_cubic_spline( nmesh(it), radr(1:nmesh(it)),&
               &                  psi1(1:nmesh(it)), wk_psi1(1:nmesh(it) ) )
          call init_cubic_spline( nmesh(it), radr(1:nmesh(it)),&
               &                  phi1(1:nmesh(it)), wk_phi1(1:nmesh(it) ) )

          Do lmt2=1, ilmt(it)
             il2 = ltp(lmt2,it); im2 = mtp(lmt2,it);  tau2 = taup(lmt2,it)
             nspher2 = ( il2 -1 )**2 +im2

             psi2(1:nmesh(it)) = psirpw(1:nmesh(it),il2,tau2,it) /radr(1:nmesh(it))
             phi2(1:nmesh(it)) = phirpw(1:nmesh(it),il2,tau2,it) /radr(1:nmesh(it))

             call calc_diff_exp( ier,5, nmesh(it), radr, psi2, dpsi2 )
             call calc_diff_exp( ier,5, nmesh(it), radr, phi2, dphi2 )

             call sphr( nma, nspher2, cx(1:nma), cy(1:nma), cz(1:nma), ylm2(1:nma) )

             call init_cubic_spline( nmesh(it), radr(1:nmesh(it)),&
                  &                  psi2(1:nmesh(it)), wk_psi2(1:nmesh(it) ) )
             call init_cubic_spline( nmesh(it), radr(1:nmesh(it)),&
                  &                  phi2(1:nmesh(it)), wk_phi2(1:nmesh(it) ) )
             call init_cubic_spline( nmesh(it), radr(1:nmesh(it)),&
                  &                  dpsi2(1:nmesh(it)), wk_dpsi2(1:nmesh(it) ) )
             call init_cubic_spline( nmesh(it), radr(1:nmesh(it)),&
                  &                  dphi2(1:nmesh(it)), wk_dphi2(1:nmesh(it) ) )

             do imesh=1,nma
                call cubic_spline( nmesh(it), radr, psi1(1:nmesh(it)),&
                     &             wk_psi1(1:nmesh(it)), rdiff(imesh), ps1, dps1 )
                call cubic_spline( nmesh(it), radr, psi2(1:nmesh(it)),&
                     &             wk_psi2(1:nmesh(it)), rdiff(imesh), ps2, dps2 )
                call cubic_spline( nmesh(it), radr, dpsi2(1:nmesh(it)),&
                     &             wk_dpsi2(1:nmesh(it)), rdiff(imesh), dps2, ddps2 )

                call cubic_spline( nmesh(it), radr, phi1(1:nmesh(it)),&
                     &             wk_phi1(1:nmesh(it)), rdiff(imesh), ph1, dph1 )
                call cubic_spline( nmesh(it), radr, phi2(1:nmesh(it)),&
                     &             wk_phi2(1:nmesh(it)), rdiff(imesh), ph2, dph2 )
                call cubic_spline( nmesh(it), radr, dphi2(1:nmesh(it)),&
                     &             wk_dphi2(1:nmesh(it)), rdiff(imesh), dph2, ddph2 )

                if (rdiff(imesh) < radr(2) ) cycle

                csum_psi = ps1 *( -ddps2 -2.0d0/rdiff(imesh) *dps2 &
                     &            +il2*(il2+1.0d0)/rdiff(imesh)**2 *ps2  )&
                     &         *ylm1(imesh) *ylm2(imesh)
                csum_phi = ph1 *( -ddph2 -2.0d0/rdiff(imesh) *dph2 &
                     &            +il2*(il2+1.0d0)/rdiff(imesh)**2 *ph2  )&
                     &         *ylm1(imesh) *ylm2(imesh)
                kina_qrs_rspace(imesh,ia,lmt1,lmt2) = ( csum_psi -csum_phi )/2.0d0
             End do
          End do
       End Do
    End Do

    deallocate( cx ); deallocate( cy ); deallocate( cz ); deallocate( rdiff )
    deallocate( ylm1 ); deallocate( ylm2 ); deallocate( dylm1 ); deallocate( dylm2 )
    deallocate( psi1 ); deallocate( psi2 ); deallocate( wk_psi1 ); deallocate( wk_psi2 )
    deallocate( phi1 ); deallocate( phi2 ); deallocate( wk_phi1 ); deallocate( wk_phi2 )
    deallocate( wk_dpsi2 ); deallocate( wk_dphi2 )

    deallocate( radr )

  end subroutine m_KE_calc_kina_qrs_rspace

  subroutine m_KE_add_ekindens_hard_rspace( ista_rho, iend_rho, ekin_dens )
    integer, intent(in) :: ista_rho, iend_rho
    real(kind=DP), intent(inout) :: ekin_dens( ista_rho:iend_rho, nspin )

    integer :: ia, it, nma
    integer :: lmt1, lmt2, is, i, ind

    real(kind=DP) :: fac0, fac

    real(kind=DP), allocatable :: kins_tmp(:,:)

!    return

!    write(*,*) 'nffp_nonp ', nfftp_nonpara
!    write(*,*) iend_fftph
!    stop

    allocate( kins_tmp(nmesh_rs_aug_max, nspin)); kins_tmp = 0.0d0

    Do ia=1, natm
       it = ityp(ia)
       nma = nmesh_rs_aug(ia)

       kins_tmp = 0.0d0
       Do lmt1=1, ilmt(it)
          Do lmt2=lmt1, ilmt(it)

             fac0=2.0d0
             if (lmt1.eq.lmt2) fac0=1.d0

             Do is=1, nspin
                fac = hsr(ia,lmt1,lmt2,is) *fac0

                if ( use_symm_ekin_density ) then
                   Do i=1, nma
                      kins_tmp(i,is) = kins_tmp(i,is) &
                           &            +fac *kins_qrs_rspace(i,ia,lmt1,lmt2)

                   End do
                else if ( use_asymm_ekin_density ) then
                   Do i=1, nma
                      kins_tmp(i,is) = kins_tmp(i,is) &
                           &            +fac *kina_qrs_rspace(i,ia,lmt1,lmt2)

                   End do
                endif
             End do
           End do
        End do

!        kins_tmp = kins_tmp /univol**2

        Do i=1, nma
           ind = meshxyz_rs_aug(i,ia)
           if ( ind >= ista_rho .and. ind <= iend_rho ) then
              ekin_dens(ind,1:nspin) = ekin_dens(ind,1:nspin) +kins_tmp(i,1:nspin)

!              write(900+mype,*) ind, kins_tmp(i,1)
           endif
        End do
     End Do

     deallocate( kins_tmp )

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
            &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )

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

  subroutine m_KE_init_ekindens
    use m_Parallelization,  only :  nel_fftp, idisp_fftp, ista_fftp, iend_fftp, &
         &                          mp_fftp, nis_fftp, nie_fftp

    integer :: is, i
    real(kind=DP) :: factor
    real(kind=DP), allocatable :: afft(:), wk_l(:,:)

!    return

    call m_FFT_alloc_CD_box()
    allocate(afft(ista_fftp:iend_fftp))
    allocate(wk_l(ista_kngp:iend_kngp,kimg) )

    factor = 0.3d0 *( 6.0d0 *PAI**2 )**(2./3.)

    call set_mat_on_fftmesh_rspace( chgq_l, afft )
    Do i=ista_fftph, iend_fftph
       afft(2*i-1) =  abs( afft(2*i-1) )**(5./3.) *factor
    End Do
    call extract_matrix_from_FFTcd_mesh( afft, wk_l )

    if ( ndim_magmom == 1 .or. ndim_magmom == 4 ) then
       ekins_l(:,:,1) = wk_l(:,:)
    else
       ekins_l(:,:,1) = wk_l(:,:) /2.0d0
       ekins_l(:,:,2) = ekins_l(:,:,1)
    endif

    deallocate( afft )
    deallocate( wk_l )
    call m_FFT_dealloc_CD_box()

  contains

    subroutine set_mat_on_fftmesh_rspace( wkq_l, afft )
      real(kind=DP), intent(in) :: wkq_l(ista_kngp:iend_kngp,kimg,ndim_magmom)
      real(kind=DP), intent(out) :: afft(ista_fftp:iend_fftp)

      integer :: j, i, ip
      real(kind=DP) :: rinplw
      real(kind=DP), allocatable :: afft_mpi1(:), afft_mpi2(:), afft_mpi3(:)

      rinplw = 1.d0 /product(fft_box_size_CD(1:3,1))

      afft = 0.0d0;

      allocate(afft_mpi1(nfftp)); afft_mpi1 = 0.0d0
      if(npes >= 2) then
         allocate(afft_mpi2(mp_fftp)); afft_mpi2 = 0.0d0
         allocate(afft_mpi3(mp_fftp)); afft_mpi3 = 0.0d0
      end if
      if ( ndim_magmom == 1 .or. ndim_magmom == 4 ) then
         do j = 1, kimg
            do i = ista_kngp, iend_kngp  !for mpi
               ip = (igfp_l(i)-1)*kimg + j
               afft_mpi1(ip) = wkq_l(i,j,1)
            end do
         end do
      else
         do j = 1, kimg
            do i = ista_kngp, iend_kngp  !for mpi
               ip = (igfp_l(i)-1)*kimg + j
               afft_mpi1(ip) = wkq_l(i,j,1) +wkq_l(i,j,2)
            end do
         end do
      endif

      if (npes >= 2) then
         call mpi_barrier(MPI_CommGroup,ierr)
         do j = 0, npes-1
            do i = nis_fftp(j),nie_fftp(j)
               afft_mpi2(i-nis_fftp(j)+1) = afft_mpi1(i)
            end do
            call mpi_allreduce(afft_mpi2,afft_mpi3,mp_fftp &
                 &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
            if(j == mype) then
               do i = ista_fftp, iend_fftp
                  afft(i) = afft_mpi3(i - ista_fftp + 1)
               end do
            end if
         end do
      else
         afft = afft_mpi1
      end if

      call m_FFT_CD_inverse_c(nfout,afft)        ! G-->R space
      deallocate(afft_mpi1)
      if (allocated(afft_mpi2) ) deallocate(afft_mpi2)
      if (allocated(afft_mpi3) ) deallocate(afft_mpi3)

    end subroutine set_mat_on_fftmesh_rspace

  end subroutine m_KE_init_ekindens

end module m_KineticEnergy_Density
