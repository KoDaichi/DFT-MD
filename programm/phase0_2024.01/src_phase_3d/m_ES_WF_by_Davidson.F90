#undef NEC_TIMER
#ifdef NEC_TIMER
#  define START_TIMER(a) call start_timer(a)
#  define STOP_TIMER(a)  call stop_timer(a)
#else
#  define START_TIMER(a)
#  define STOP_TIMER(a)
#endif
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  MODULE:  m_ES_WF_by_Davidson
!
!  AUTHOR(S): T. Yamamoto   June/12/2005
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
!
module m_ES_WF_by_Davidson
! $Id: m_ES_WF_by_Davidson.F90 570 2017-04-21 20:34:50Z yamasaki $
#ifdef TRANSPOSE
#ifdef VPP
#define _ODD_BOUNDARY_
#endif
#ifdef SX
#define _ODD_BOUNDARY_
#endif
#ifdef NEC_TUNE1
#define _ODD_BOUNDARY_
#endif
  use m_Const_Parameters,    only : DP,SP,DIRECT,ON,OFF,SCF,GAMMA &
        &                         , SmallestPositiveNumber,SKIP,EXECUT,zi, ELECTRON, OLD
  use m_Parallelization,     only : MPI_CommGroup &
       &                          , myrank_e,myrank_k,map_e,map_k,ista_e,iend_e,istep_e &
       &                          , ista_k,iend_k,np_g1k,ista_g1,mpi_k_world,ierr,map_z &
       &                          , np_e,npes,nrank_e,mype
  use m_Control_Parameters,  only : nspin,ipridavidson,kimg,neg,af,ndavid &
       &                          , max_iter_david, delta_eig_occup, delta_eig_empty &
#ifdef SAVE_FFT_TIMES
       &                          , eps_david, max_subspace_size, sw_save_fft
#else
       &                          , eps_david, max_subspace_size
#endif
  use m_Control_Parameters,  only : sw_dav_scalapack, dav_nprow, dav_npcol, dav_divide_square &
                                  , dav_block_size, sw_fft_xzy
  use m_Files,               only : nfout
  use m_Timing,              only : tstatc0_begin, tstatc0_end
!fj$$F  use m_FFT,                 only : nfft,fft_box_size_WF, m_FFT_Vlocal_W, m_FFT_WF
  use m_FFT,                 only : nfft,fft_box_size_WF
  use m_Kpoints,             only : kv3,vkxyz, k_symmetry
  use m_PlaneWaveBasisSet,   only : kg1, iba, igf, nbase, m_pwBS_kinetic_energies
  use m_Electronic_Structure,only : zaj_l, neordr, nrvf_ordr, eko_l, vlhxcQ &
#ifdef SAVE_FFT_TIMES
 &                                , status_saved_phifftr &
#endif
       &                          , occup_l&
       &                          , fsr_l,fsi_l, vnlph_l, vlhxc_l
  use m_ES_ortho            ,only : np_g1k_x
  use m_ES_nonlocal         ,only : m_ES_alloc_scss_etc_3D                            &
 &                                , m_ES_dealloc_scss_etc
  use m_Ionic_System,       only : natm, iwei, ityp, ntyp
  use m_PseudoPotential,    only : ilmt,nlmta,lmta,q,dion &
       &                         , lmtt,ltp,mtp &
       &                         , m_PP_include_vanderbilt_pot &
       &                         , ipaw,dion_paw
  use m_NonLocal_Potential, only : snl
!fj$$F use m_Electronic_Structure,only : zaj_l_3D, eko_l_3D, vnlph_l_3D, snl_l_3D          &
!fj$$F &                                , fsr_l_3D, fsi_l_3D                                &
 use m_Electronic_Structure,only :  m_ES_wd_zaj_small_portion_3D                      &
 &                                , m_ES_wd_eko_3D
  use m_ES_nonlocal         ,only : m_ES_Vnonlocal_W_3D                               &
 &                                , m_ES_betar_dot_WFs_4_each_k_3D                    &
!fj$$F &                                , m_ES_Vnonlocal_W_3D_new                           &
!fj$$F &                                , m_ES_betar_dot_WFs_4_each_k_3D_new
 &                                , m_ES_Vnonlocal_W_3D                           &
 &                                , m_ES_betar_dot_WFs_4_each_k_3D
  use m_Parallelization,     only : ista_g1k, iend_g1k, neg_g, mpi_ke_world &
       &                          , ista_k, iend_k, mpi_kg_world &
       &                          , np_fs, myrank_g, nis_fs
!fj$$F  use z_interface_3D, only        : decomp_zaj_l_3D_ik, decomp_zaj_l_r_3D_ik &
!fj$$F       &                          , decomp_eko_l_3D_new, decomp_vnlph_l_3D &
!fj$$F       &                          , decomp_snl_l_3D_2, decomp_fsr_l_r_3D_ik, decomp_fsr_l_3D_ik &
!fj$$F       &                          , decomp_eko_l_r_3D_new
#ifdef NEC_TIMER
  use nec_timer
#endif
! === FFT Marge. by T.Kato ===============================================================
  use m_Parallelization,     only : nel_fft_x , nel_fft_y, nel_fft_z &
                                  , fft_X_x_nel, fft_X_y_nel, fft_X_z_nel &
                                  , mp_g1k, myrank_g
  use m_Electronic_Structure,only : m_ES_Vlocal_in_Rspace_3D, m_ES_WF_in_Rspace_3D &
!fj$$F                                  , vlhxc_l, vlhxc_l_3D
                                  , vlhxc_l
  use m_FFT,                 only : m_FFT_Vlocal_W_3D, m_FFT_Direct_3D             &
#ifdef FFT_3D_DIVISION
 &                                , m_FFT_Vlocal_W_3DIV_3D, m_FFT_Direct_3DIV_3D   &
#endif
 &                                , m_FFT_Direct_XYZ_3D
  use m_ES_WF_by_SDorCG,     only : map_fft_to_WF_3D
!fj$$F  use z_interface_3D,        only : decomp_vlhxc_l_3D
! ========================================================================================
! === Marge to Riken Source by T.Kato ====================================================
  use m_Electronic_Structure,only : zaj_ball, nblocksize_mgs_default
!fj$$F  use z_interface_3D,        only : replacement_zaj_ball_sequence
! ========================================================================================
! === DEBUG by tkato 2012/04/06 ================================================
  use m_Parallelization,     only : nrank_g
! ==============================================================================


! ==================================== added by K. Tagami ============== 11.0
  use m_Const_Parameters,    only : CMPLDP, Neglected, BuiltIn
  use m_Control_Parameters,  only : noncol, ndim_spinor, ndim_chgpot, ndim_magmom, &
       &                            SpinOrbit_mode, sw_hubbard
  use m_PseudoPotential,     only : q_noncl, dion_scr_noncl
!  use m_Electronic_Structure,  only : m_ES_Vlocal_in_Rspace_noncl, &
!       &                              m_ES_sort_eigen_vals_noncl
!
!  use m_FFT,                 only : m_FFT_Vlocal_W_noncl
! ====================================================================== 11.0
  use mpi

  implicit none

  integer :: nsize_subspace, nsize_matrix
  integer :: nsize_sb_now, nsize_mt_now, nsize_mt_old
  real(kind=DP), allocatable, dimension(:) :: w1hw2
  real(kind=DP), allocatable, dimension(:) :: w1sw2
  real(kind=DP), allocatable, dimension(:) :: w1_mpi
  real(kind=DP), allocatable, dimension(:,:,:,:) :: zat_l
  real(kind=DP), allocatable, dimension(:,:,:) :: zah_l
  real(kind=DP), allocatable, dimension(:,:,:) :: fsr,fsi
  real(kind=DP), allocatable, dimension(:,:,:) :: zaj_l_backup
  logical, allocatable, dimension(:) :: feigconv
  integer, allocatable, dimension(:) :: ibover
#ifdef _USE_SCALAPACK_
  integer :: nprow, npcol, nrowa, ncola
  integer, dimension(9) :: desca
  integer :: lwork1, liwork1
  integer :: lwork2, liwork2, lrwork2
  integer :: ctxt
! === DEBUG by tkato 2012/12/20 ================================================
  integer :: myprow, mypcol
! ==============================================================================
#endif

! ================================= added by K. Tagami ================= 11.0
  real(kind=DP), allocatable, dimension(:,:,:,:,:) :: zat_l_noncl
  real(kind=DP), allocatable, dimension(:,:,:,:) :: zah_l_noncl
  real(kind=DP), allocatable, dimension(:,:,:,:) :: fsr_noncl, fsi_noncl
  real(kind=DP), allocatable, dimension(:,:,:,:) :: zaj_l_backup_noncl
! ===================================================================== 11.0

!  include 'mpif.h'
contains

  subroutine m_ESdavidson_Renew_WF(nfout,precon)
    integer, intent(in) :: nfout,precon

    integer             :: ispin, ik, iksnl, switch_of_eko_part
    integer :: idavid,itot
    real(kind=DP), allocatable, dimension(:) ::  afft, bfft
    real(kind=DP), dimension(maxval(np_g1k)) :: ekin_l,vnldi_l,hdiag_l,sdiag_l
    real(kind=DP) :: vlhxc0
    logical :: frestart
    integer :: max_idavid, max_itot
    integer :: idavid_now, itot_now, ipri0
! === FFT Marge. by T.Kato ===============================================================
    integer :: lsize, ibsize, isrsize, fft_l_size
    real(kind=DP), allocatable, dimension(:) :: afft_l
    real(kind=DP), allocatable, dimension(:,:) :: wk_bfft_l
    real(kind=DP), allocatable, dimension(:,:) :: bfft_l
! ========================================================================================
    integer :: n_unconv
! ========================================================================================

    max_idavid = 1
    max_itot = 1

    call m_ES_alloc_scss_etc_3D()
    allocate(afft(nfft)); allocate(bfft(nfft))

    call allocate_matrix
! ==============================================================================
START_TIMER('Davidson_Initialize')
!fj$$F    do ispin = 1, nspin, (af+1)
!fj$$F       call decomp_vlhxc_l_3D(vlhxc_l,vlhxc_l_3D,ispin)
!fj$$F       do ik = ispin, kv3-nspin+ispin, nspin
!fj$$F          if(map_k(ik) /= myrank_k) cycle
!fj$$F          call decomp_zaj_l_3D_ik(zaj_l,zaj_l_3D,ik,nrvf_ordr,"    ")
!fj$$F          call decomp_snl_l_3D_2(snl,snl_l_3D,ik)
!fj$$F          call decomp_fsr_l_3D_ik(fsr_l,fsr_l_3D,ik,nrvf_ordr,"    ")
!fj$$F          if(.not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) &
!fj$$F             call decomp_fsr_l_3D_ik(fsi_l,fsi_l_3D,ik,nrvf_ordr,"    ")
!fj$$F          call decomp_eko_l_3D_new(eko_l,eko_l_3D,ik,nrvf_ordr,"    ")
!fj$$F       enddo
!fj$$F    enddo
STOP_TIMER('Davidson_Initialize')
! ==============================================================================
! === FFT Marge. by T.Kato ===============================================================

#ifdef FFT_3D_DIVISION
    lsize = fft_X_x_nel*fft_X_y_nel*fft_X_z_nel
    allocate(afft_l(lsize*2), stat=ierr)
#else
    lsize = max(maxval(nel_fft_x(:)),maxval(nel_fft_y(:)),maxval(nel_fft_z(:)))
    allocate(afft_l(lsize*kimg), stat=ierr)
#endif
    if(ierr /= 0) then
       write(nfout,*)' m_ESdavidson_Renew_WF : Not allocated afft_l array'
       call flush(nfout)
       call mpi_abort(mpi_comm_world, 201, ierr)
    endif
    ibsize = 1
!   if (nblocksize_fftw_is_given) then
!      ibsize = nblocksize_fftw
!      if (ibsize < 1) ibsize = 1
!   endif
! ========================================================================================
#ifdef _USE_SCALAPACK_
! === DEBUG by tkato 2012/12/20 ================================================
    if(sw_dav_scalapack == ON) call create_scalapack_context()
! ==============================================================================
#endif
    do ispin = 1, nspin, (af+1)
! === FFT Marge. by T.Kato ===============================================================
       call m_ES_Vlocal_in_Rspace_3D(ispin,afft_l,lsize,1,OFF) ! (ptfft1) ->afft
! ========================================================================================
       call vlhxc_l_zero_term(vlhxc0,ispin)        ! vlhxc_l -> vlhxc0
       do ik = ispin, kv3-nspin+ispin, nspin
          if(map_k(ik) /= myrank_k) cycle          ! MPI
! === FFT Marge. by T.Kato ===============================================================
          isrsize = min(lsize,mp_g1k(ik))
          fft_l_size  = nel_fft_x(myrank_g)
#ifdef FFT_3D_DIVISION
          allocate(wk_bfft_l(lsize*2,ibsize) ,stat=ierr)
          allocate(bfft_l(lsize*2,ibsize) ,stat=ierr)
#else
          allocate(wk_bfft_l(lsize*kimg,ibsize) ,stat=ierr)
          allocate(bfft_l(lsize*kimg,ibsize) ,stat=ierr)
#endif
          if (ierr /= 0) then
             write(nfout,*)' m_ESdavidson_Renew_WF :  Not allocate '
             call flush(nfout)
             call mpi_abort(mpi_comm_world, 205, ierr)
          endif
! ========================================================================================
          iksnl = (ik-1)/nspin + 1
          call allocate_t_matrix_3D(ik) ! -> np_g1k_x
          call m_pwBS_kinetic_energies(ik,vkxyz,ekin_l) ! (diakin) ->ekin
          if(precon==ON) then
             call Vnonlocal_Diagonal_part_3D(ispin,ik,iksnl,vnldi_l)
             hdiag_l(1:maxval(np_g1k)) = ekin_l(1:maxval(np_g1k)) + vlhxc0 + vnldi_l(1:maxval(np_g1k))
             call S_Diagonal_part_3D(ik,iksnl,sdiag_l)
          end if
          !!$write(nfout,*) 'debug loop ik=',ik
          max_itot=0
          max_idavid=0
          Loop: do itot=1,max_iter_david
             itot_now = itot
             if(itot>max_itot) max_itot=itot
             feigconv = .false.
             David_Loop: do idavid=1,ndavid
                idavid_now = idavid
                if(idavid>max_idavid) max_idavid=idavid
                !!$write(nfout,*) 'debug loop itot,idavid=',itot,idavid
                if(precon==ON) then
!fj$$F                   call m_ES_Vnonlocal_W_3D_new(ik,iksnl,ispin,switch_of_eko_part=OFF) ! -> vnlph_l
                   call m_ES_Vnonlocal_W_3D(ik,iksnl,ispin,switch_of_eko_part=OFF) ! -> vnlph_l
                else
!fj$$F                   call m_ES_Vnonlocal_W_3D_new(ik,iksnl,ispin,switch_of_eko_part=ON) ! -> vnlph_l
                   call m_ES_Vnonlocal_W_3D(ik,iksnl,ispin,switch_of_eko_part=ON) ! -> vnlph_l
                end if
START_TIMER('Davidson_build_subspace_3D')
!               call build_subspace_3D(ik,ekin,hdiag,sdiag,afft,bfft,idavid,precon) ! -> zajold_l
                call build_subspace_3D(ik,ekin_l,hdiag_l,sdiag_l,afft_l,bfft_l,wk_bfft_l, &
                                       lsize,ibsize,isrsize,fft_l_size,idavid,precon) ! -> zajold_l
STOP_TIMER('Davidson_build_subspace_3D')
START_TIMER('Davidson_allreduce_fs_3D')
                call allreduce_fs_3D(ik,idavid) ! -> fsr,fsi
STOP_TIMER('Davidson_allreduce_fs_3D')
START_TIMER('Davidson_evolve_WFs_in_subspace_3D')
!               call evolve_WFs_in_subspace_3D&     !-(m_ES_WF_by_Davidson)
!                 &(ik,ispin,ekin,afft,bfft,idavid,frestart) !-> zaj_l
                call evolve_WFs_in_subspace_3D&     !-(m_ES_WF_by_Davidson)
                  &(ik,ispin,ekin_l,afft_l,bfft_l,wk_bfft_l, &
                    lsize,ibsize,isrsize,fft_l_size,idavid,frestart) !-> zaj_l
STOP_TIMER('Davidson_evolve_WFs_in_subspace_3D')
                if(ik==1.and.ipridavidson>= 2) then
! ==============================================================================
!fj$$F                   call decomp_zaj_l_r_3D_ik(zaj_l,zaj_l_3D,ik,neordr,"    ")
! ==============================================================================
                   call m_ES_wd_zaj_small_portion_3D(nfout,ik," -- after davidson --",21)
                endif
!fj$$F                call m_ES_betar_dot_WFs_4_each_k_3D_new(nfout,ik)   ! -> fsr_l,fsi_l
                call m_ES_betar_dot_WFs_4_each_k_3D(nfout,ik)   ! -> fsr_l,fsi_l
!               if(nsize_sb_now > max_subspace_size) exit David_Loop
                if(frestart) exit David_Loop
                if(eigenvalues_are_converged_3D(n_unconv)) exit Loop
                if(nsize_sb_now + n_unconv > max_subspace_size) exit David_Loop
             end do David_Loop
          end do Loop
          if(ipridavidson>=2) then
             write(nfout,'("Davidson: ik=",i5," itot=",i5," idavid=",i5," subspace=",i5)') &
                  & ik, itot_now, idavid_now, nsize_sb_now
          end if
          call deallocate_t_matrix_3D
! === FFT Marge. by T.Kato ===============================================================
          deallocate(wk_bfft_l)
          deallocate(bfft_l)
! ========================================================================================
       enddo      ! k-point loop
    enddo      ! spin loop
#ifdef _USE_SCALAPACK_
! === DEBUG by tkato 2012/12/20 ================================================
    if(sw_dav_scalapack == ON) call destroy_scalapack_context()
! ==============================================================================
#endif
! === FFT Marge. by T.Kato ===============================================================
    deallocate(afft_l)
! ========================================================================================
! ==============================================================================
START_TIMER('Davidson_Finalize')
!fj$$F    do ispin = 1, nspin, (af+1)
!fj$$F       do ik = ispin, kv3-nspin+ispin, nspin
!fj$$F          if(map_k(ik) /= myrank_k) cycle
!fj$$F          call decomp_zaj_l_r_3D_ik(zaj_l,zaj_l_3D,ik,neordr,"    ")
!fj$$F          call decomp_fsr_l_r_3D_ik(fsr_l,fsr_l_3D,ik,neordr,"    ",0)
!fj$$F          if(.not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) &
!fj$$F             call decomp_fsr_l_r_3D_ik(fsi_l,fsi_l_3D,ik,neordr,"    ",0)
!fj$$F          call decomp_eko_l_r_3D_new(eko_l,eko_l_3D,ik,neordr,"    ")
!fj$$F          call replacement_zaj_ball_sequence(zaj_ball,ik,neordr,nblocksize_mgs_default)
!fj$$F       enddo
!fj$$F    enddo
STOP_TIMER('Davidson_Finalize')
! ==============================================================================
    call deallocate_matrix
    if(ipridavidson>=2) then
       write(nfout,'("Davidson: max_itot=",i5," max_idavid=",i5)') max_itot,max_idavid
    end if

!!  ( in case of af=1 )
    if(af /= 0) then
       call cp_eigen_values_for_af       !-(contained here)
       call expand_neordr_and_nrvf_ordr  !-(contained here)
    end if

    call get_ipri0(ipridavidson,ipri0)
    if(ipri0 >= 2) call m_ES_wd_eko_3D(nfout,mode=SCF)
    deallocate(bfft);   deallocate(afft)
    call m_ES_dealloc_scss_etc()


  contains

    subroutine get_ipri0(ipri_in, ipri_out)
      integer, intent(in)  :: ipri_in
      integer, intent(out) :: ipri_out
      if(npes > 1) then
         if(mype == 0) ipri_out = ipri_in
         call mpi_bcast(ipri_out,1,mpi_integer,0,MPI_CommGroup,ierr)
      else
         ipri_out = ipri_in
      end if
    end subroutine get_ipri0

    logical function eigenvalues_are_converged_3D(n_unconv)
       integer, intent(out) :: n_unconv
       integer :: ib
       n_unconv = 0
       eigenvalues_are_converged_3D = .true.
       do ib=1,neg
          if(.not.feigconv(ib)) then
             eigenvalues_are_converged_3D = .false.
             n_unconv = n_unconv + 1
          end if
       end do
    end function eigenvalues_are_converged_3D

    subroutine cp_eigen_values_for_af
      integer :: ik,ib
      do ik = 1, kv3, af+1
         if(map_k(ik) /= myrank_k) cycle    ! MPI
         do ib = 1, np_e                    ! MPI
            eko_l(ib,ik+af) = eko_l(ib,ik)
         enddo
      enddo
    end subroutine cp_eigen_values_for_af

    subroutine expand_neordr_and_nrvf_ordr
      integer :: ik
      do ik = 1, kv3, af+1
         if(map_k(ik) /= myrank_k) cycle     ! MPI
         neordr(1:neg,ik+af) = neordr(1:neg,ik)
         nrvf_ordr(1:neg,ik+af) = nrvf_ordr(1:neg,ik)
      end do
    end subroutine expand_neordr_and_nrvf_ordr

  end subroutine m_ESdavidson_Renew_WF


  subroutine allocate_matrix
    nsize_subspace = neg*ndavid
    nsize_matrix =  max_subspace_size*(max_subspace_size+1)/2
    allocate(feigconv(neg))
    allocate(ibover(nsize_subspace))
    allocate(fsr(neg,nlmta,ndavid))
    if(.not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
       allocate(fsi(neg,nlmta,ndavid))
    end if
! ==============================================================================
    allocate(zaj_l_backup(maxval(np_g1k),np_e,kimg)) ! MPI
! ==============================================================================
  end subroutine allocate_matrix


  subroutine deallocate_matrix
    deallocate(feigconv)
    deallocate(ibover)
    deallocate(fsr)
    if(.not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
       deallocate(fsi)
    end if
    deallocate(zaj_l_backup)
  end subroutine deallocate_matrix


  subroutine allocate_t_matrix_3D(ik)
    integer, intent(in) :: ik
    integer :: kimg_t
    if(k_symmetry(ik) == GAMMA) then
       kimg_t = 1
    else
       kimg_t = kimg
    end if
#ifdef _ODD_BOUNDARY_
    if(mod(np_g1k(ik),2) == 0) then
       np_g1k_x = np_g1k(ik) + 1
    else
       np_g1k_x = np_g1k(ik)
    end if
#else
    np_g1k_x = np_g1k(ik)
#endif
    allocate(zat_l(maxval(np_g1k),neg,kimg,ndavid)) ! MPI
    allocate(zah_l(maxval(np_g1k),np_e,kimg)) ! MPI
    allocate(w1hw2(nsize_matrix*kimg_t))
    allocate(w1sw2(nsize_matrix*kimg_t))
  end subroutine allocate_t_matrix_3D

  subroutine deallocate_t_matrix_3D
    deallocate(zat_l) ! MPI
    deallocate(zah_l) ! MPI
    deallocate(w1hw2)
    deallocate(w1sw2)
  end subroutine deallocate_t_matrix_3D

  subroutine allreduce_fs_3D(ik,idavid)
    integer, intent(in) :: ik,idavid
    integer :: ib,ib1,iter,kimg_t,is,is1
    real(kind=DP), allocatable, dimension(:,:) :: fs_mpi,fs_mpi2

    allocate(fs_mpi(neg,nlmta))
    allocate(fs_mpi2(neg,nlmta))
    if(idavid>1) then
!fj$$F       call m_ES_betar_dot_WFs_4_each_k_3D_new(nfout,ik)   ! -> fsr_l,fsi_l
       call m_ES_betar_dot_WFs_4_each_k_3D(nfout,ik)   ! -> fsr_l,fsi_l
    end if
    fs_mpi=0.d0
    do is = 1, np_fs ! MPI
       is1=nis_fs(myrank_g)+is-1
       do ib = 1, np_e ! MPI
          ib1=neg_g(ib)
          fs_mpi(ib1,is1) = fsr_l(ib,is,ik)
       end do
    end do
    fs_mpi2=0.d0
    call mpi_allreduce(fs_mpi,fs_mpi2,neg*nlmta &
      & ,mpi_double_precision,mpi_sum &
      & ,mpi_k_world(myrank_k),ierr) ! MPI
    fsr(1:neg,1:nlmta,idavid) = fs_mpi2(1:neg,1:nlmta)
    if(.not. k_symmetry(ik) == GAMMA) then
       fs_mpi=0.d0
       do is = 1, np_fs ! MPI
          is1=nis_fs(myrank_g)+is-1
          do ib = 1, np_e ! MPI
             ib1=neg_g(ib)
             fs_mpi(ib1,is1) = fsi_l(ib,is,ik)
          end do
       end do
       fs_mpi2=0.d0
       call mpi_allreduce(fs_mpi,fs_mpi2,neg*nlmta &
         & ,mpi_double_precision,mpi_sum &
         & ,mpi_k_world(myrank_k),ierr) ! MPI
       fsi(1:neg,1:nlmta,idavid) = fs_mpi2(1:neg,1:nlmta)
    end if
    deallocate(fs_mpi)
    deallocate(fs_mpi2)
  end subroutine allreduce_fs_3D

  subroutine build_subspace_3D(ik,ekin_l,hdiag_l,sdiag_l,afft_l,bfft_l,wk_bfft_l, &
                               lsize,ibsize,isrsize,fft_l_size,idavid,precon)
    integer, intent(in) :: ik,idavid,precon,lsize,ibsize,isrsize,fft_l_size
    real(kind=DP), intent(in)  :: ekin_l(maxval(np_g1k)),hdiag_l(maxval(np_g1k)), &
                                  sdiag_l(maxval(np_g1k))
#ifdef FFT_3D_DIVISION
    real(kind=DP), intent(in)  :: afft_l(lsize*2   )
    real(kind=DP), intent(out) :: bfft_l(lsize*2   ,1)
    real(kind=DP), intent(inout) :: wk_bfft_l(lsize*2   ,ibsize)
#else
    real(kind=DP), intent(in)  :: afft_l(lsize*kimg)
!   real(kind=DP), intent(out) :: bfft_l(lsize*kimg,ibesize)
    real(kind=DP), intent(out) :: bfft_l(lsize*kimg,1)
    real(kind=DP), intent(inout) :: wk_bfft_l(lsize*kimg,ibsize)
#endif
    integer :: iter,ib1,ib,ii,i1,iadd
    real(kind=DP) :: dr1,dr2,di1
    real(kind=DP) :: p(maxval(np_g1k))
    real(kind=DP) :: denom, norm, zajr, zaji

    denom = 1.d0/product(fft_box_size_WF(1:3,1))

    if(idavid==1) then
       return
    end if
    ! if precon == ON
    ! |R> = P*(H-Hdiag)|Psi>
    ! P = (Hdiag-e*Sdiag)^(-1)
    ! else
    ! |R> = (H-eS)|Psi>
    do ib1 = 1, np_e ! MPI
       if(feigconv(neg_g(ib1))) cycle
       ib = ib1 ! MPI
START_TIMER('Davidson2_temprary_FFT_1')
#ifdef __TIMER_COMM__
       call m_ES_WF_in_Rspace_3D(ik,ib,ib,ibsize,lsize,wk_bfft_l,0)
#else
       call m_ES_WF_in_Rspace_3D(ik,ib,ib,ibsize,lsize,wk_bfft_l)
#endif
#ifdef FFT_3D_DIVISION
       call m_FFT_Vlocal_W_3DIV_3D(afft_l,wk_bfft_l,lsize,ibsize,nel_fft_x(myrank_g))
       call m_FFT_Direct_3DIV_3D(nfout,wk_bfft_l,lsize,ibsize)
#else
       if (sw_fft_xzy > 0) then
          call m_FFT_Vlocal_W_3D(afft_l,wk_bfft_l,lsize,ibsize,nel_fft_y(myrank_g))
          call m_FFT_Direct_3D(nfout,wk_bfft_l,lsize,ibsize)
       else
          call m_FFT_Vlocal_W_3D(afft_l,wk_bfft_l,lsize,ibsize,nel_fft_z(myrank_g))
          call m_FFT_Direct_XYZ_3D(nfout,wk_bfft_l,lsize,ibsize)
       end if
#endif
       call map_fft_to_WF_3D(ik,lsize,ibsize,wk_bfft_l,bfft_l,isrsize,fft_l_size)
#ifdef SAVE_FFT_TIMES
       if(sw_save_fft == ON) status_saved_phifftr(ib,ik) = OLD
#endif
STOP_TIMER('Davidson2_temprary_FFT_1')
       if(precon==ON) then
          call decide_precon_factor_david(ik,hdiag_l,sdiag_l,eko_l(ib,ik),p) ! -> p
          if(kimg == 1) then
             do ii=ista_g1k(ik),iend_g1k(ik)
                iadd = ii - ista_g1k(ik) + 1
                i1  = igf(nbase(ii,ik))
                dr1 = zaj_l(iadd,ib,ik,kimg)
                dr2 = bfft_l(iadd,1)*denom
                zajr= (ekin_l(iadd)-hdiag_l(iadd))*dr1+dr2+vnlph_l(iadd,ib,kimg)
                zaj_l(iadd,ib,ik,kimg) = p(iadd)*zajr
             enddo
          else
             do ii=ista_g1k(ik),iend_g1k(ik)
                iadd = ii - ista_g1k(ik) + 1
                i1  = igf(nbase(ii,ik))
                dr1  = zaj_l(iadd,ib,ik,1)
                di1  = zaj_l(iadd,ib,ik,kimg)
                dr2  = ekin_l(iadd) - hdiag_l(iadd)
                zajr= dr2*dr1+bfft_l(2*iadd-1,1)*denom + vnlph_l(iadd,ib,1)
                zaji= dr2*di1+bfft_l(2*iadd,  1)*denom + vnlph_l(iadd,ib,2)
                zaj_l(iadd,ib,ik,1   )=p(iadd)*zajr
                zaj_l(iadd,ib,ik,kimg)=p(iadd)*zaji
             enddo
          endif
       else
          if(kimg == 1) then
             do ii=ista_g1k(ik),iend_g1k(ik)
                iadd = ii - ista_g1k(ik) + 1
                i1  = igf(nbase(ii,ik))
                dr1 = zaj_l(iadd,ib,ik,kimg)
                dr2 = bfft_l(iadd,1)*denom
                zajr= ekin_l(iadd)*dr1+dr2+vnlph_l(iadd,ib,kimg)
                zaj_l(iadd,ib,ik,kimg) = zajr
             enddo
          else
             do ii=ista_g1k(ik),iend_g1k(ik)
                iadd = ii - ista_g1k(ik) + 1
                i1  = igf(nbase(ii,ik))
                dr1  = zaj_l(iadd,ib,ik,1)
                di1  = zaj_l(iadd,ib,ik,kimg)
                zajr= ekin_l(iadd)*dr1+bfft_l(2*iadd-1,1)*denom + vnlph_l(iadd,ib,1)
                zaji= ekin_l(iadd)*di1+bfft_l(2*iadd,  1)*denom + vnlph_l(iadd,ib,2)
                zaj_l(iadd,ib,ik,1   )=zajr
                zaj_l(iadd,ib,ik,kimg)=zaji
             enddo
          endif
       end if

       ! Normalize
       norm = 0.d0
       if(kimg==1) then
          do ii=ista_g1k(ik),iend_g1k(ik)
             iadd = ii - ista_g1k(ik) + 1
             zajr = zaj_l(iadd,ib,ik,kimg)
             norm = norm + zajr*zajr
          end do
          call mpi_allreduce(MPI_IN_PLACE,norm,1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
          norm = 1.d0/sqrt(norm)
          do ii=ista_g1k(ik),iend_g1k(ik)
             iadd = ii - ista_g1k(ik) + 1
             zaj_l(iadd,ib,ik,kimg) = zaj_l(iadd,ib,ik,kimg)*norm
          end do
       else
          do ii=ista_g1k(ik),iend_g1k(ik)
             iadd = ii - ista_g1k(ik) + 1
             zajr = zaj_l(iadd,ib,ik,1   )
             zaji = zaj_l(iadd,ib,ik,kimg)
             norm = norm + zajr*zajr+zaji*zaji
          end do
          call mpi_allreduce(MPI_IN_PLACE,norm,1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
          norm = 1.d0/sqrt(norm)
          do ii=ista_g1k(ik),iend_g1k(ik)
             iadd = ii - ista_g1k(ik) + 1
             zaj_l(iadd,ib,ik,1)   = zaj_l(iadd,ib,ik,1)*norm
             zaj_l(iadd,ib,ik,kimg)= zaj_l(iadd,ib,ik,kimg)*norm
          enddo
       end if
    end do
  contains
    subroutine decide_precon_factor_david(ik,hdiag_l,sdiag_l,eig,p)
      integer, intent(in) :: ik
      real(kind=DP), intent(in) :: hdiag_l(maxval(np_g1k)), &
                                   sdiag_l(maxval(np_g1k)),eig
      real(kind=DP), intent(out) :: p(maxval(np_g1k))

      integer :: ii, iadd
      real(kind=DP) :: denom

      do ii=ista_g1k(ik),iend_g1k(ik)
         iadd = ii - ista_g1k(ik) + 1
         denom = hdiag_l(iadd)-eig*sdiag_l(iadd)
         if(abs(denom) < eps_david) then
            denom = sign(eps_david,denom)
         end if
         p(iadd) = 1.d0/denom
      end do

    end subroutine decide_precon_factor_david
  end subroutine build_subspace_3D

  subroutine evolve_WFs_in_subspace_3D&
       &(ik,ispin,ekin_l,afft_l,bfft_l,wk_bfft_l, &
         lsize,ibsize,isrsize,fft_l_size,idavid,frestart)
    integer, intent(in) :: ik,ispin,lsize,ibsize,isrsize,fft_l_size
    integer, intent(in) :: idavid
    logical, intent(out) :: frestart
    real(kind=DP), intent(in)  :: ekin_l(maxval(np_g1k))
#ifdef FFT_3D_DIVISION
    real(kind=DP), intent(in)  :: afft_l(lsize*2   )
    real(kind=DP), intent(out) :: bfft_l(lsize*2   ,1)
    real(kind=DP), intent(inout) :: wk_bfft_l(lsize*2   ,ibsize)
#else
    real(kind=DP), intent(in)  :: afft_l(lsize*kimg)
!   real(kind=DP), intent(out) :: bfft_l(lsize*kimg,ibesize)
    real(kind=DP), intent(out) :: bfft_l(lsize*kimg,1)
    real(kind=DP), intent(inout) :: wk_bfft_l(lsize*kimg,ibsize)
#endif
! (allocatable variables)
    real(kind=DP), allocatable,dimension(:) ::     eig
    real(kind=DP), allocatable,dimension(:,:) ::   vec
    real(kind=DP), allocatable,dimension(:) ::     eko_d
    real(kind=DP), allocatable,dimension(:) ::     eko_d_mpi
    integer, allocatable,dimension(:) ::     occup
    integer, allocatable,dimension(:) ::     occup_mpi

    integer       :: ib1,ib2,ib1to,ib2to,i1,ii,ri,ib
    integer       :: ibb1,ibb2
    integer       :: ii1,ii2,iter,iter1,iter2
    real(kind=DP) :: denom, eko1, eko2, ekod
    real(kind=DP) :: hr2,hi2,dr1,dr2,di1,di2,dd
    integer :: ip0,ip0b,ip1,ip1b,ib1n,ib2n,ndata,nshift,kimg_t,ig1
    integer :: noffset
    integer :: nsize_max_sb_now
    integer :: ierr_diag
    integer :: id_sname = -1, ipri0
! ========================================================================================
    integer :: ib1_, ib2_, ibmax, iadd
! ========================================================================================
    call tstatc0_begin('evolve_WFs_in_subspace_3D (davidson) ', id_sname,1)

    call get_ipri0(ipridavidson,ipri0)

    denom = 1.d0/product(fft_box_size_WF(1:3,1))
    if(k_symmetry(ik) == GAMMA) then
       kimg_t = 1
    else
       kimg_t = kimg
    end if

    if(idavid==1) then
       ip0=neg
       do ib=1,neg
          ibover(ib) = ib
       end do
       nsize_mt_now=0
    else
       noffset = neg*(idavid-1)
       ip0=nsize_sb_now
       do ib=1,neg
          if(.not.feigconv(ib)) then
             ip0=ip0+1
             ibover(noffset+ib) = ip0
          else
             ibover(noffset+ib) = -1
          end if
       end do
    end if
    nsize_sb_now = ip0
    nsize_mt_old = nsize_mt_now
    nsize_mt_now = nsize_sb_now*(nsize_sb_now+1)/2
    nsize_max_sb_now = neg*idavid
    if(ipridavidson >=2) then
       write(nfout,*) 'ibover=',ibover(1:nsize_max_sb_now)
    end if

    allocate(eig(nsize_sb_now)); eig=0.d0
    allocate(vec(nsize_sb_now*kimg_t,nsize_sb_now))
    allocate(eko_d(neg));     eko_d = 0.d0
    allocate(eko_d_mpi(neg))
    allocate(occup(neg)); occup=0
    allocate(occup_mpi(neg))

    if(ipridavidson >=2) then
       write(nfout,*) 'Davidson:ik,idavid,nsize_sb_now=', ik,idavid,nsize_sb_now
    end if


    do ib1 = 1, np_e
       eko_d(neg_g(ib1)) = eko_l(ib1,ik)  ! MPI
    end do
    call mpi_allreduce(eko_d,eko_d_mpi,neg,mpi_double_precision,mpi_sum &
         & ,mpi_kg_world,ierr) ! MPI
    eko_d = eko_d_mpi                         ! MPI
    eko1 = sum(eko_d(1:neg))

    do ib1=1,neg
       if(map_e(ib1) == myrank_e) then !MPI
          if( occup_l(map_z(ib1),ik) > 0.d0 ) occup(ib1) = 1  ! MPI
       end if
    end do
    call mpi_allreduce(occup,occup_mpi,neg,mpi_integer,mpi_sum &
         & ,mpi_k_world(myrank_k),ierr)       ! MPI
    occup = occup_mpi                         ! MPI

! (zaj_l <- (T+Vloc)|phi> )
!( tenchi ) (zat_l <- zaj_l)
!   call m_ES_W_transpose(ista_k,iend_k,ik,zaj_l,zat_l(1,1,1,idavid))
! ==============================================================================
    zat_l(:,:,:,idavid) = 0.0d0
    do ib1 = 1, np_e
       zat_l(:,neg_g(ib1),:,idavid) = zaj_l(:,ib1,ik,:)
    enddo
    call mpi_allreduce(MPI_IN_PLACE,zat_l(:,:,:,idavid),maxval(np_g1k)*neg*kimg,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_kg_world,ierr)
! ==============================================================================
    do ib1 = 1, np_e ! MPI
       ib = ib1 ! MPI
START_TIMER('Davidson2_temprary_FFT_2')
#ifdef __TIMER_COMM__
       call m_ES_WF_in_Rspace_3D(ik,ib,ib,ibsize,lsize,wk_bfft_l,0)
#else
       call m_ES_WF_in_Rspace_3D(ik,ib,ib,ibsize,lsize,wk_bfft_l)
#endif
#ifdef FFT_3D_DIVISION
       call m_FFT_Vlocal_W_3DIV_3D(afft_l,wk_bfft_l,lsize,ibsize,nel_fft_x(myrank_g))
       call m_FFT_Direct_3DIV_3D(nfout,wk_bfft_l,lsize,ibsize)
#else
       if (sw_fft_xzy > 0) then
          call m_FFT_Vlocal_W_3D(afft_l,wk_bfft_l,lsize,ibsize,nel_fft_y(myrank_g))
          call m_FFT_Direct_3D(nfout,wk_bfft_l,lsize,ibsize)
       else
          call m_FFT_Vlocal_W_3D(afft_l,wk_bfft_l,lsize,ibsize,nel_fft_z(myrank_g))
          call m_FFT_Direct_XYZ_3D(nfout,wk_bfft_l,lsize,ibsize)
       end if
#endif
       call map_fft_to_WF_3D(ik,lsize,ibsize,wk_bfft_l,bfft_l,isrsize,fft_l_size)
STOP_TIMER('Davidson2_temprary_FFT_2')
       if(kimg == 1) then
          do ii=ista_g1k(ik),iend_g1k(ik)
             iadd = ii - ista_g1k(ik) + 1
             i1  = igf(nbase(ii,ik))
             dr1 = zaj_l(iadd,ib,ik,1)
             dr2 = bfft_l(iadd,1)*denom
             zaj_l(iadd,ib,ik,1)= ekin_l(iadd)*dr1+dr2
          enddo
       else
          do ii=ista_g1k(ik), iend_g1k(ik)
             iadd = ii - ista_g1k(ik) + 1
             i1  = igf(nbase(ii,ik))
             dr1  = zaj_l(iadd,ib,ik,1)
             di1  = zaj_l(iadd,ib,ik,kimg)
             zaj_l(iadd,ib,ik,1)= ekin_l(iadd)*dr1+bfft_l(2*iadd-1,1)*denom
             zaj_l(iadd,ib,ik,kimg)= ekin_l(iadd)*di1+bfft_l(2*iadd,1)*denom
          enddo
       endif
#ifdef SAVE_FFT_TIMES
       if(sw_save_fft == ON) status_saved_phifftr(ib,ik) = OLD
#endif
    enddo
!( tenchi ) (zah_l <- zaj_l)
!   call m_ES_W_transpose(ista_k,iend_k,ik,zaj_l,zah_l(1,1,1,idavid))
! ==============================================================================
    zah_l(:,:,:) = zaj_l(:,:,ik,:)
! ==============================================================================

! (make matrix elements )
! parallel loop
    ! <n|T+Vloc|m> G-wise parallel
! === Need Initialization to Change Decomposition!!! by T.Kato =================
    w1hw2(nsize_mt_old*kimg_t+1:nsize_mt_now*kimg_t) = 0.0d0
    w1sw2(nsize_mt_old*kimg_t+1:nsize_mt_now*kimg_t) = 0.0d0
! ==============================================================================
    do ib2_ = 1, np_e
       ibb2 = neg_g(ib2_) + neg*(idavid - 1)
       if(ibover(ibb2)<0) cycle
       ib2 = ibover(ibb2)
       iter2 = idavid
       ip0b = ib2*(ib2-1)/2
! ==============================================================================
       do iter = 1, idavid
! ==============================================================================
       ibmax = neg
       if(iter == idavid) ibmax = neg_g(ib2_)
       do ib1_ = 1, ibmax
          ibb1 = ib1_ + neg*(iter - 1)
          if(ibover(ibb1)<0) cycle
          ib1 = ibover(ibb1)
          iter1 = iter
          ip0 = ip0b + ib1
          if(kimg == 1) then
             do ii = ista_g1k(ik), iend_g1k(ik) ! MPI
                iadd = ii - ista_g1k(ik) + 1
                hr2 = zah_l(iadd,ib2_,1)
                dr2 = zat_l(iadd,neg_g(ib2_),1,iter2)
                dr1 = zat_l(iadd,ib1_,1,iter1)
                w1hw2(ip0) = w1hw2(ip0) + dr1*hr2
                w1sw2(ip0) = w1sw2(ip0) + dr1*dr2
             end do
          else
             if(k_symmetry(ik) == GAMMA) then
                do ii = max(ista_g1k(ik),2), iend_g1k(ik) ! MPI
                   iadd = ii - ista_g1k(ik) + 1
                   hr2 = zah_l(iadd,ib2_,1) ! MPI
                   hi2 = zah_l(iadd,ib2_,2) ! MPI
                   dr2 = zat_l(iadd,neg_g(ib2_),1,iter2) ! MPI
                   di2 = zat_l(iadd,neg_g(ib2_),2,iter2) ! MPI
                   dr1 = zat_l(iadd,ib1_,1,iter1) ! MPI
                   di1 = zat_l(iadd,ib1_,2,iter1) ! MPI
                   w1hw2(ip0) =w1hw2(ip0)+(dr1*hr2+di1*hi2)*2.d0
                   w1sw2(ip0) =w1sw2(ip0)+(dr1*dr2+di1*di2)*2.d0
                end do
                if(ista_g1k(ik) == 1) then
                   hr2 = zah_l(1,ib2_,1) ! MPI
                   hi2 = zah_l(1,ib2_,2) ! MPI
                   dr2 = zat_l(1,neg_g(ib2_),1,iter2) ! MPI
                   di2 = zat_l(1,neg_g(ib2_),2,iter2) ! MPI
                   dr1 = zat_l(1,ib1_,1,iter1) ! MPI
                   di1 = zat_l(1,ib1_,2,iter1) ! MPI
                   w1hw2(ip0) =w1hw2(ip0)+dr1*hr2+di1*hi2
                   w1sw2(ip0) =w1sw2(ip0)+dr1*dr2+di1*di2
                endif
             else
                do ii = ista_g1k(ik), iend_g1k(ik) ! MPI
                   iadd = ii - ista_g1k(ik) + 1
                   hr2 = zah_l(iadd,ib2_,1) ! MPI
                   hi2 = zah_l(iadd,ib2_,2) ! MPI
                   dr2 = zat_l(iadd,neg_g(ib2_),1,iter2) ! MPI
                   di2 = zat_l(iadd,neg_g(ib2_),2,iter2) ! MPI
                   dr1 = zat_l(iadd,ib1_,1,iter1) ! MPI
                   di1 = zat_l(iadd,ib1_,2,iter1) ! MPI
                   w1hw2(2*ip0-1) =w1hw2(2*ip0-1)+dr1*hr2+di1*hi2
                   w1hw2(2*ip0  ) =w1hw2(2*ip0  )+dr1*hi2-di1*hr2
                   w1sw2(2*ip0-1) =w1sw2(2*ip0-1)+dr1*dr2+di1*di2
                   w1sw2(2*ip0  ) =w1sw2(2*ip0  )+dr1*di2-di1*dr2
                end do
             end if
          end if
       end do
! ==============================================================================
       enddo
! ==============================================================================
    end do

    if(ipridavidson >= 3) call wd_w1hw2(" -- w1hw2 without nl part--")
    ! <n|Vnl|m> G-wise parallel
START_TIMER('Davidson2_add_nonlocal_part')
! === DEBUG by tkato 2011/07/14 ================================================
    if(myrank_g == 0) then
! ==============================================================================
    call add_nonlocal_part ! w1hw2 = w1hw2 + w1Vnlw2
                           ! w1sw2 = w1sw2 + w1qw2
! === DEBUG by tkato 2011/07/14 ================================================
    endif
! ==============================================================================
STOP_TIMER('Davidson2_add_nonlocal_part')
    if(ipridavidson >= 3) call wd_w1hw2(" -- w1hw2 with nl part--")

!! (spread sum of w1hw2 and w1sw2)
    if(npes > 1) then
       allocate(w1_mpi(nsize_matrix*kimg_t))
       nshift = nsize_mt_old*kimg_t
       ndata = (nsize_mt_now-nsize_mt_old)*kimg_t
       w1_mpi = 0.0d0
       call mpi_allreduce(w1hw2(nshift+1),w1_mpi,ndata,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr) ! MPI
       w1hw2(nshift+1:nshift+ndata) = w1_mpi(1:ndata) ! MPI
       call mpi_allreduce(w1sw2(nshift+1),w1_mpi,ndata,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr) ! MPI
       w1sw2(nshift+1:nshift+ndata) = w1_mpi(1:ndata) ! MPI
       deallocate(w1_mpi)
    end if

    if(ipridavidson >= 2) call wd_w1hw2(" -- just after making w1hw2 --")
    if(ipridavidson >= 2) then
       write(nfout,*) 'neordr for ik = ',ik
       write(nfout,9002) (neordr(ib1,ik),ib1=1,neg)
       write(nfout,*) 'nrvf_ordr for ik = ',ik
       write(nfout,9002) (nrvf_ordr(ib1,ik),ib1=1,neg)
       write(nfout,*) 'eig'
       write(nfout,'(5x,10f8.4)') (eko_l(ib1,ik),ib1=1,np_e)
    endif
9002 format(5x,10i8)

!! (Diagonalization )  !!

#ifdef _USE_LAPACK_
#ifdef _USE_SCALAPACK_
    if(sw_dav_scalapack == ON) then
       call initialize_scalapack()
       if(kimg_t == 1) then
START_TIMER('Davidson_pdsygvx')
          call pdsygvx_driver(eig,vec,w1hw2,w1sw2,ierr_diag)
STOP_TIMER('Davidson_pdsygvx')
       else
START_TIMER('Davidson_pzhegvx')
          call pzhegvx_driver(eig,vec,w1hw2,w1sw2,ierr_diag)
STOP_TIMER('pzhegvx_zhpgvx')
       endif
! === DEBUG by tkato 2012/12/20 ================================================
!      call finalize_scalapack()
! ==============================================================================
    else
       if(kimg_t == 1) then
START_TIMER('Davidson_dspgvx')
          call dspgvx_driver(eig,vec,w1hw2,w1sw2,ierr_diag)
STOP_TIMER('Davidson_dspgvx')
       else
START_TIMER('Davidson_zhpgvx')
          call zhpgvx_driver(eig,vec,w1hw2,w1sw2,ierr_diag)
STOP_TIMER('Davidson_zhpgvx')
       endif
    endif
#else
    if(kimg_t == 1) then
START_TIMER('Davidson_dspgvx')
       call dspgvx_driver(eig,vec,w1hw2,w1sw2,ierr_diag)
STOP_TIMER('Davidson_dspgvx')
    else
START_TIMER('Davidson_zhpgvx')
       call zhpgvx_driver(eig,vec,w1hw2,w1sw2,ierr_diag)
STOP_TIMER('Davidson_zhpgvx')
    endif
#endif
#else
    call phase_error_with_msg(nfout,'Complie me with -D_USE_LAPACK_ to use the Davidoson solver',__LINE__,__FILE__)
    !!if(kimg_t == 1) then
    !!   call hobsvw_driver(neg,eig,w1hw2)
    !!else
    !!   call chobsd_driver(neg,eig,w1hw2)
    !!end if
#endif

    frestart = .false.
    if(ierr_diag /= 0) then
       zaj_l(:,:,ik,:) = zaj_l_backup(:,:,:)
       do ib1 = 1, np_e
          eko_l(ib1,ik)=eko_d(neg_g(ib1))
       end do
       frestart = .true.
       if(ipridavidson >= 2) then
          write(nfout,*) '** restart Davidson iteration **'
       end if
    else

       feigconv = .false.
       do ib=1,neg
          if(occup(ib) == 1) then
             if(abs(eko_d(ib)-eig(ib)) < delta_eig_occup) feigconv(ib) = .true.
          else
             if(abs(eko_d(ib)-eig(ib)) < delta_eig_empty) feigconv(ib) = .true.
          end if
       end do

       if(ipridavidson >= 2) then
          write(nfout,*) 'eko_d for ik = ',ik
          write(nfout,9001) (eko_d(ib),ib=1,neg)
          write(nfout,*) 'eig for ik = ',ik
          write(nfout,9001) (eig(ib),ib=1,neg)
          call wd_w1hw2(" -- after diagonalization --")
!sum eko
          dr1=0.d0;dr2=0.d0
          do ib1=1,np_e
             ib1to = neordr(neg_g(ib1),ik)
             if(map_e(ib1to) == myrank_e) dr1=dr1+eko_l(ib1,ik) ! MPI
             dr2=dr2+eig(neg_g(ib1))
          enddo
          call mpi_allreduce(dr1,di1,1,mpi_double_precision,mpi_sum,mpi_kg_world,ierr) ! MPI
          dr1 = di1  ! MPI
          write(nfout,'(" sum of eko_l, eig, abs diff =",3e25.10)') dr1,dr2,abs(dr2-dr1)
       endif
!! (subspace rotation) !!
START_TIMER('Davidson2_subspace_rotation')
       call subspace_rotation ! vec,zat_l -> zat_l
STOP_TIMER('Davidson2_subspace_rotation')
!( tenchi ) (zaj_l <- zat_l)
       iter = min(idavid+1,ndavid)
!      call m_ES_W_transpose_back(ista_k,iend_k,ik,zaj_l,zat_l(1,1,1,iter))
! ==============================================================================
       do ib1 = 1, np_e
          zaj_l(:,ib1,ik,:) = zat_l(:,neg_g(ib1),:,iter)
       enddo
! ==============================================================================
       zaj_l_backup(:,:,:) = zaj_l(:,:,ik,:)
!! (eko_l)
       do ib1 = 1, np_e
          eko_l(ib1,ik)=eig(neg_g(ib1))
       end do
       if(ipridavidson >= 2) then
          eko2 = sum(eig(1:neg))
          write(nfout,1201) ik,eko1,ekod,eko2

          write(nfout,*) 'eko_l'
          write(nfout,9001) (eko_l(ib1,ik),ib1=1,np_e)
       endif
1201 format(' %% for ik = ',i4,4x,' eko1&ekod&eko2 = ',3f14.7)
9001 format(5x,6f12.5)
!! (neordr & nrvf_ordr)

    end if

    neordr(1:neg,ik) = (/(ib1,ib1=1,neg)/)
    nrvf_ordr(1:neg,ik) = (/(ib1,ib1=1,neg)/)

! (deallocate)
    deallocate(eko_d)
    deallocate(eko_d_mpi)
    deallocate(eig)
    deallocate(vec)
    deallocate(occup)
    deallocate(occup_mpi)

    call tstatc0_end(id_sname)

  contains
    subroutine get_ipri0(ipri_in, ipri_out)
      integer, intent(in)  :: ipri_in
      integer, intent(out) :: ipri_out
      if(npes > 1) then
         if(mype == 0) ipri_out = ipri_in
         call mpi_bcast(ipri_out,1,mpi_integer,0,mpi_k_world(myrank_k),ierr)
      else
         ipri_out = ipri_in
      end if
    end subroutine get_ipri0

    subroutine wd_w1hw2(somecomment)
      character(len=*), intent(in) :: somecomment
      integer :: ib1, ib2, neg_wd, nsb_wd
      write(nfout,'(a35)') somecomment
      write(nfout,*) 'w1hw2 for ik = ',ik
      neg_wd = 8
      nsb_wd = 8
      if(neg_wd > neg) neg_wd = neg
      if(nsb_wd > nsize_sb_now) nsb_wd = nsize_sb_now
      if(kimg_t==1) then
         do ib2=1,nsb_wd
            ip0b = ib2*(ib2-1)/2
            write(nfout,9001) (w1hw2(ip0),ip0=ip0b+1,ip0b+ib2)
         end do
      else
         do ib2=1,nsb_wd
            ip0b = ib2*(ib2-1)/2
            write(nfout,9001) (w1hw2(2*ip0-1),ip0=ip0b+1,ip0b+ib2)
         end do
         do ib2=1,nsb_wd
            ip0b = ib2*(ib2-1)/2
            write(nfout,9001) (w1hw2(2*ip0),ip0=ip0b+1,ip0b+ib2)
         end do
      end if
      write(nfout,*) 'w1sw2 for ik = ',ik
      if(kimg_t==1) then
         do ib2=1,nsb_wd
            ip0b = ib2*(ib2-1)/2
            write(nfout,9001) (w1sw2(ip0),ip0=ip0b+1,ip0b+ib2)
         end do
      else
         do ib2=1,nsb_wd
            ip0b = ib2*(ib2-1)/2
            write(nfout,9001) (w1sw2(2*ip0-1),ip0=ip0b+1,ip0b+ib2)
         end do
         do ib2=1,nsb_wd
            ip0b = ib2*(ib2-1)/2
            write(nfout,9001) (w1sw2(2*ip0),ip0=ip0b+1,ip0b+ib2)
         end do
      end if
9001  format(5x,9f12.5)
      write(nfout,*) 'eko_l for ik = ',ik
      write(nfout,9001) (eko_d(neordr(ib1,ik)),ib1=1,neg_wd)
    end subroutine wd_w1hw2

    subroutine add_nonlocal_part
      integer :: ip,ib1,ib2,ibb1,ibb2
      integer       :: ia, lmt1, lmt2, it, p, s, ib
      real(kind=DP) :: facv,facq,vr,vi,qr,qi
      real(kind=DP) :: tmpr,tmpi
      do ibb2 = nsize_max_sb_now-neg+1,nsize_max_sb_now
         if(ibover(ibb2)<0) cycle
         ib2 = ibover(ibb2)
         iter2= (ibb2-1)/neg+1
         ii2  = ibb2-neg*(iter2-1)
         ip0b = ib2*(ib2-1)/2
         do ibb1 = 1,ibb2
            if(ibover(ibb1)<0) cycle
            ib1 = ibover(ibb1)
            iter1=(ibb1-1)/neg+1
            ii1  = ibb1-neg*(iter1-1)
            ip0 = ip0b + ib1
            if(mod(ip0-1,nrank_e)/=myrank_e) cycle
            if(kimg_t==1) then
               vr=0.d0
               qr=0.d0
            else
               vr=0.d0
               vi=0.d0
               qr=0.d0
               qi=0.d0
            end if
            do ia = 1, natm
               it = ityp(ia)
               do lmt1 = 1, ilmt(it)
                  p = lmta(lmt1,ia)
                  do lmt2 = 1, ilmt(it)
                     s = lmta(lmt2,ia)
!!$                     facv   = iwei(ia)*(dion(lmt1,lmt2,it) + vlhxcQ(lmt1,lmt2,ia,ispin))
                     if(ipaw(it)==0)then
                        facv   = iwei(ia)*(dion(lmt1,lmt2,it) + vlhxcQ(lmt1,lmt2,ia,ispin))
                     else
                        facv   = iwei(ia)*(dion_paw(lmt1,lmt2,ispin,ia) + vlhxcQ(lmt1,lmt2,ia,ispin))
                     endif
                     facq   = iwei(ia)*q(lmt1,lmt2,it)
                     if(kimg==1) then
                        tmpr = fsr(ii1,p,iter1)*fsr(ii2,s,iter2)&
                    &        + fsi(ii1,p,iter1)*fsi(ii2,s,iter2)
                        vr = vr + facv*tmpr
                        qr = qr + facq*tmpr
                     else
                        if(k_symmetry(ik) == GAMMA) then
                           tmpr = fsr(ii1,p,iter1)*fsr(ii2,s,iter2)
                           vr = vr + facv*tmpr
                           qr = qr + facq*tmpr
                        else
                           tmpr = fsr(ii1,p,iter1)*fsr(ii2,s,iter2)&
                    &        + fsi(ii1,p,iter1)*fsi(ii2,s,iter2)
                           tmpi = fsr(ii1,p,iter1)*fsi(ii2,s,iter2)&
                    &        - fsi(ii1,p,iter1)*fsr(ii2,s,iter2)
                           vr = vr + facv*tmpr
                           vi = vi + facv*tmpi
                           qr = qr + facq*tmpr
                           qi = qi + facq*tmpi
                        end if
                     end if
                  end do
               end do
            end do
            if(kimg_t==1) then
               w1hw2(ip0) = w1hw2(ip0) + vr
               w1sw2(ip0) = w1sw2(ip0) + qr
            else
               w1hw2(2*ip0-1) = w1hw2(2*ip0-1) + vr
               w1hw2(2*ip0  ) = w1hw2(2*ip0  ) + vi
               w1sw2(2*ip0-1) = w1sw2(2*ip0-1) + qr
               w1sw2(2*ip0  ) = w1sw2(2*ip0  ) + qi
            end if
         end do
      end do
    end subroutine add_nonlocal_part

    subroutine subspace_rotation
      integer :: ib1,ib2,ibb2,iadd
      real(kind=DP), dimension(maxval(np_g1k),np_e,kimg) :: zaj_wk

      zaj_wk(:,:,:) = 0.d0
      if(kimg==1) then
         do ib1=1,np_e
! ==============================================================================
            do iter = 1, idavid
! ==============================================================================
            do ib2_=1,neg
               ibb2 = ib2_ + neg*(iter - 1)
               if(ibover(ibb2)<0) cycle
               ib2 = ibover(ibb2)
               do ii=ista_g1k(ik),iend_g1k(ik)
                  iadd = ii - ista_g1k(ik) + 1
                  zaj_wk(iadd,ib1,kimg) = zaj_wk(iadd,ib1,kimg) + zat_l(iadd,ib2_,kimg,iter)*vec(ib2,neg_g(ib1))
               end do
            end do
! ==============================================================================
            end do
! ==============================================================================
         end do
      else
         if(k_symmetry(ik) == GAMMA) then
            do ib1=1,np_e
! ==============================================================================
               do iter = 1, idavid
! ==============================================================================
               do ib2_=1,neg
                  ibb2 = ib2_ + neg*(iter - 1)
                  if(ibover(ibb2)<0) cycle
                  ib2 = ibover(ibb2)
                  hr2=vec(ib2,neg_g(ib1))
                  do ii=ista_g1k(ik),iend_g1k(ik)
                     iadd = ii - ista_g1k(ik) + 1
                     dr1=zat_l(iadd,ib2_,1   ,iter)
                     di1=zat_l(iadd,ib2_,kimg,iter)
                     zaj_wk(iadd,ib1,1   ) = zaj_wk(iadd,ib1,1   ) + dr1*hr2
                     zaj_wk(iadd,ib1,kimg) = zaj_wk(iadd,ib1,kimg) + di1*hr2
                  end do
               end do
! ==============================================================================
               end do
! ==============================================================================
            end do
         else
            do ib1=1,np_e
! ==============================================================================
            do iter = 1, idavid
! ==============================================================================
               do ib2_=1,neg
                  ibb2 = ib2_ +  neg*(iter - 1)
                  if(ibover(ibb2)<0) cycle
                  ib2 = ibover(ibb2)
                  hr2=vec(2*ib2-1,neg_g(ib1))
                  hi2=vec(2*ib2  ,neg_g(ib1))
                  do ii=ista_g1k(ik),iend_g1k(ik)
                     iadd = ii - ista_g1k(ik) + 1
                     dr1=zat_l(iadd,ib2_,1   ,iter)
                     di1=zat_l(iadd,ib2_,kimg,iter)
                     zaj_wk(iadd,ib1,1   ) = zaj_wk(iadd,ib1,1   ) + dr1*hr2 - di1*hi2
                     zaj_wk(iadd,ib1,kimg) = zaj_wk(iadd,ib1,kimg) + dr1*hi2 + di1*hr2
                  end do
               end do
! ==============================================================================
               end do
! ==============================================================================
            end do
         end if
      end if
      iter = min(idavid+1,ndavid)
!     zat_l(:,:,:,iter) = zaj_wk(:,:,:)
! ==============================================================================
      zat_l(:,:,:,iter) = 0.0d0
      do ib1 = 1, np_e
         zat_l(:,neg_g(ib1),:,iter) = zaj_wk(:,ib1,:)
      enddo
      call mpi_allreduce(MPI_IN_PLACE,zat_l(:,:,:,iter),maxval(np_g1k)*neg*kimg,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_kg_world,ierr)
! ==============================================================================
    end subroutine subspace_rotation

  end subroutine evolve_WFs_in_subspace_3D

  subroutine vlhxc_l_zero_term(vlhxc0,ispin)
    real(kind=DP), intent(out) :: vlhxc0
    integer, intent(in)        :: ispin

    if(mype == 0) vlhxc0 = vlhxc_l(1,1,ispin)
    call mpi_bcast(vlhxc0,1,mpi_double_precision,0,MPI_CommGroup,ierr)
  end subroutine vlhxc_l_zero_term

! ===================================== added by K. Tagami ============== 11.0
  subroutine vlhxc_l_zero_term_noncl( vlhxc0 )
    real(kind=DP), intent(out) :: vlhxc0(ndim_chgpot)

    vlhxc0 = 0.0d0
    if (mype == 0) then
!	vlhxc0(1)           = vlhxc_l(1,1,1) + vlhxc_l(1,1,ndim_magmom)
!	vlhxc0(ndim_chgpot) = vlhxc_l(1,1,1) - vlhxc_l(1,1,ndim_magmom)
!
       vlhxc0 = vlhxc_l(1,1,1)
    endif

    call mpi_bcast( vlhxc0, ndim_chgpot, mpi_double_precision,0,MPI_CommGroup,ierr)

  end subroutine vlhxc_l_zero_term_noncl
! ========================================================================= 11.0

  subroutine Vnonlocal_Diagonal_part_3D(ispin,ik,iksnl,vnldi_l)
    integer, intent(in)                        :: ispin, ik, iksnl
    real(kind=DP), intent(out), dimension(maxval(np_g1k)) :: vnldi_l

    integer :: it,mdvdb

    vnldi_l = 0.d0
    do it = 1, ntyp
       mdvdb = m_PP_include_vanderbilt_pot(it)
       if(mdvdb == SKIP) then
          call Vnonlocal_D_norm_conserve_case
       else if(mdvdb == EXECUT) then
          call Vnonlocal_D_vanderbilt_case
       end if
    end do
  contains
    subroutine Vnonlocal_D_vanderbilt_case
      integer       :: ia, p1,p2,lmtt1,il1,lmtt2,il2,i,iadd
      real(kind=DP) :: ph,fac

      do p1 = 1, ilmt(it)
         lmtt1 = lmtt(p1,it); il1 = ltp(p1,it)
         do p2 = p1, ilmt(it)
            lmtt2 = lmtt(p2,it); il2 = ltp(p2,it)
            if( p1 /= p2) then
               ph = 2.d0*real(zi**(il2-il1))
            else
               ph = 1.d0
            endif
            if(mod(il1+il2,2) == 1) cycle
            do ia = 1, natm
               if(ityp(ia) /= it) cycle
               if(ipaw(it)==0) then
                  fac = ph*iwei(ia) * (dion(p1,p2,it)+vlhxcQ(p1,p2,ia,ispin))
               else
                  fac = ph*iwei(ia) * (dion_paw(p1,p2,ispin,ia)+vlhxcQ(p1,p2,ia,ispin))
               endif
               do i = ista_g1k(ik), iend_g1k(ik)
                  iadd = i - ista_g1k(ik) + 1
                  vnldi_l(iadd) = vnldi_l(iadd)+fac*snl(iadd,lmtt1,iksnl)*snl(iadd,lmtt2,iksnl)
               end do
            end do
         end do
      end do
    end subroutine Vnonlocal_D_vanderbilt_case

    subroutine Vnonlocal_D_norm_conserve_case
      integer       :: ia, lmt1,lmt2,lmtt1,il1,im1,il2,im2,i,iadd
      real(kind=DP) :: ph,fac

      ph = 0.d0
      do ia = 1, natm
         if(ityp(ia) /= it) cycle
         ph = ph + iwei(ia)
      end do
      do lmt1 = 1, ilmt(it)
         lmtt1 = lmtt(lmt1,it); il1 = ltp(lmt1,it); im1 = mtp(lmt1,it)
         do lmt2 = lmt1, ilmt(it)
            il2 = ltp(lmt2,it); im2 = mtp(lmt2,it)
            if(il1 /= il2 .or. im1 /= im2) cycle
            if(mod(il1+il2,2) == 1) cycle
            if(ipaw(it)==0)then
               fac = ph * dion(lmt1,lmt2,it)
            else
               fac = ph * dion_paw(lmt1,lmt2,ispin,ia)
            endif
            do i = ista_g1k(ik), iend_g1k(ik)
               iadd = i - ista_g1k(ik) + 1
               vnldi_l(iadd)  = vnldi_l(iadd) + fac * snl(iadd,lmtt1,iksnl)*snl(iadd,lmtt1,iksnl)
            end do
         end do
      end do
    end subroutine Vnonlocal_D_norm_conserve_case
  end subroutine Vnonlocal_Diagonal_part_3D

  subroutine S_Diagonal_part_3D(ik,iksnl,sdiag_l)
    integer, intent(in)                        :: ik, iksnl
    real(kind=DP), intent(out), dimension(maxval(np_g1k)) :: sdiag_l

    integer :: it,mdvdb

    sdiag_l = 1.d0
    do it = 1, ntyp
       mdvdb = m_PP_include_vanderbilt_pot(it)
       if(mdvdb /= SKIP) call Vanderbilt_case
    end do
  contains
    subroutine Vanderbilt_case
      integer       :: ia, p1,p2,lmtt1,il1,lmtt2,il2,i,iadd
      real(kind=DP) :: ph,fac

      do p1 = 1, ilmt(it)
         lmtt1 = lmtt(p1,it); il1 = ltp(p1,it)
         do p2 = p1, ilmt(it)
            lmtt2 = lmtt(p2,it); il2 = ltp(p2,it)
            if( p1 /= p2) then
               ph = 2.d0*real(zi**(il2-il1))
            else
               ph = 1.d0
            endif
            if(mod(il1+il2,2) == 1) cycle
            do ia = 1, natm
               if(ityp(ia) /= it) cycle
               fac = ph*iwei(ia) * q(p1,p2,it)
               do i = ista_g1k(ik), iend_g1k(ik)
                  iadd = i - ista_g1k(ik) + 1
                  sdiag_l(iadd) = sdiag_l(iadd)+fac*snl(iadd,lmtt1,iksnl)*snl(iadd,lmtt2,iksnl)
               end do
            end do
         end do
      end do
    end subroutine Vanderbilt_case
  end subroutine S_Diagonal_part_3D

  subroutine dspgvx_driver(eig,vec,w1hw2,w1sw2,ierr)
    real(kind=DP), intent(out) ,dimension(nsize_sb_now) :: eig
    real(kind=DP), intent(out) ,dimension(nsize_sb_now*neg) :: vec
    real(kind=DP), intent(inout) ,dimension(nsize_mt_now) :: w1hw2,w1sw2
    integer, intent(out) :: ierr
    integer :: ITYPE
    character(len=1) :: JOBZ,RANGE,UPLO
    integer :: il,iu
    real(kind=DP),allocatable,dimension(:) :: work_lapack
    integer, allocatable, dimension(:) :: iwork_lapack, ifail_lapack
    real(kind=DP) :: vl,vu,abstol
    integer :: info,m
    real(kind=DP), external :: dlamch
!!$    real(kind=DP), dimension(nsize_mt_now) :: ap,bp
    real(kind=DP), allocatable, dimension(:) :: ap,bp

    allocate(ap(nsize_mt_now),bp(nsize_mt_now))
    abstol = 2*dlamch('S')

    il=1; iu=neg
!(LAPACK)  ITYPE = 1:  A*x = (lambda)*B*x, 2:  A*B*x = (lambda)*x, 3:  B*A*x = (lambda)*x
    ITYPE = 1
!(LAPACK)  JOBZ = N : eigenvalue, V : eigenvalue + eigenvector
    JOBZ = 'V'
!(LAPACK)  RANGE = A : all eigenvalues, V: all eigenvalues in (VL, VU], I: the IL-th through IU-th eigenvalues
    RANGE = 'I'
!(LAPACK)  UPLO = U : upper triangle matrix,  L : lower triangle matrix
    UPLO = 'U'
    allocate(work_lapack(8*nsize_sb_now))
    allocate(iwork_lapack(5*nsize_sb_now))
    allocate(ifail_lapack(nsize_sb_now))

    ap = w1hw2
    bp = w1sw2

    call dspgvx(ITYPE,JOBZ,RANGE,UPLO,nsize_sb_now,ap,bp &
    &          ,vl,vu,il,iu,abstol,m,eig,vec,nsize_sb_now &
    &          ,work_lapack,iwork_lapack,ifail_lapack,info)

    if(ipridavidson >=2 .and. info/=0) then
       write(nfout,*) "debug(dspgvx) info=",info
       write(nfout,*) "debug(dspgvx) ifail"
       write(nfout,'(8(1x,i3))') ifail_lapack
       write(nfout,*) "debug(dspgvx) eig"
       write(nfout,'(8(1x,f10.5))') eig
    end if

    deallocate(work_lapack)
    deallocate(iwork_lapack)
    deallocate(ifail_lapack)

    deallocate(ap,bp)

    if(info/=0) then
       !!write(nfout,*) "dspgvx: info=",info
       !!stop 'error in dspgvx_driver'
       ierr = 1
    else
       ierr = 0
    end if

  end subroutine dspgvx_driver

#ifdef _USE_SCALAPACK_
  subroutine pdsygvx_driver(eig, vec, w1hw2, w1sw2, ierr)
    real(kind=DP), intent(out),  dimension(nsize_sb_now)      :: eig
    real(kind=DP), intent(out),  dimension(nsize_sb_now, neg) :: vec
    real(kind=DP), intent(inout),dimension(nsize_mt_now)      :: w1hw2, w1sw2
    integer, intent(out) :: ierr
    integer :: itype
    integer :: il,iu
    integer, allocatable, dimension(:) :: ifail
    real(kind=DP) :: vl, vu, abstol
    integer :: info, m
    real(kind=DP), external :: dlamch
! ==============================================================================
    real(kind=DP), allocatable, dimension(:,:) :: ap, bp, z
    integer :: nz, count
    real(kind=DP) :: orfac = -1.0d0
    real(kind=DP), allocatable, dimension(:) :: work, gap
    integer,       allocatable, dimension(:) :: iwork, iclustr
    integer :: i, j
! ==============================================================================

    abstol = 2*dlamch('S')
    il=1; iu=neg
    itype = 1

    allocate(ap(nrowa,ncola))
    allocate(bp(nrowa,ncola))
    allocate(z (nrowa,ncola))

    count = 1
    do j = 1, nsize_sb_now
      do i = 1, nsize_sb_now
        if(i .le. j) then
          call pdelset(ap,i,j,desca,w1hw2(count))
          call pdelset(bp,i,j,desca,w1sw2(count))
          count = count + 1
        endif
      enddo
    enddo

    allocate(work   (lwork1))
    allocate(iwork  (liwork1))
    allocate(ifail  (nsize_sb_now))
    allocate(iclustr(2*nprow*npcol))
    allocate(gap    (nprow*npcol))

    call pdsygvx(itype,'V','I','U',nsize_sb_now, &
                 ap,1,1,desca, &
                 bp,1,1,desca, &
                 vl,vu,il,iu,abstol,m,nz,eig,orfac, &
                 z,1,1,desca,&
                 work,lwork1,iwork,liwork1, &
                 ifail,iclustr,gap,info)

    do j = 1, neg
      do i = 1, nsize_sb_now
        call pdelget('A',' ',vec(i,j),z,i,j,desca)
      enddo
    enddo

    if(ipridavidson >=2 .and. info/=0) then
       write(nfout,*) "debug(dspgvx) info=",info
       write(nfout,*) "debug(dspgvx) ifail"
       write(nfout,'(8(1x,i3))') ifail
       write(nfout,*) "debug(dspgvx) eig"
       write(nfout,'(8(1x,f10.5))') eig
    end if

    deallocate(ap)
    deallocate(bp)
    deallocate(z)
    deallocate(work)
    deallocate(iwork)
    deallocate(iclustr)
    deallocate(gap)
    deallocate(ifail)

    if(info/=0) then
       ierr = 1
    else
       ierr = 0
    end if

  end subroutine pdsygvx_driver
#endif

  subroutine zhpgvx_driver(eig,vec,w1hw2,w1sw2,ierr)
    real(kind=DP), intent(out) ,dimension(nsize_sb_now) :: eig
    real(kind=DP), intent(out) ,dimension(nsize_sb_now*kimg,nsize_sb_now) :: vec
    real(kind=DP), intent(in) ,dimension(nsize_mt_now*kimg) :: w1hw2,w1sw2
    integer, intent(out) :: ierr
    integer :: ITYPE
    character(len=1) :: JOBZ,RANGE,UPLO
    integer :: il,iu
    real(kind=DP),allocatable,dimension(:) :: work_lapack
    real(kind=DP),allocatable,dimension(:) :: rwork_lapack
    integer, allocatable, dimension(:) :: iwork_lapack, ifail_lapack
    real(kind=DP) :: vl,vu,abstol
    integer :: info,m
    real(kind=DP), external :: dlamch
!!$    real(kind=DP), dimension(nsize_mt_now*kimg) :: ap,bp
    real(kind=DP), allocatable, dimension(:) :: ap,bp
    integer :: ib,i

    abstol = 2*dlamch('S')

    il=1; iu=neg
!(LAPACK)  ITYPE = 1:  A*x = (lambda)*B*x, 2:  A*B*x = (lambda)*x, 3:  B*A*x = (lambda)*x
    ITYPE = 1
!(LAPACK)  JOBZ = N : eigenvalue, V : eigenvalue + eigenvector
    JOBZ = 'V'
!(LAPACK)  RANGE = A : all eigenvalues, V: all eigenvalues in (VL, VU], I: the IL-th through IU-th eigenvalues
    RANGE = 'I'
!(LAPACK)  UPLO = U : upper triangle matrix,  L : lower triangle matrix
    UPLO = 'U'
    allocate(work_lapack(2*nsize_sb_now*kimg)); work_lapack = 0.d0
    allocate(rwork_lapack(7*nsize_sb_now)); rwork_lapack = 0.d0
    allocate(iwork_lapack(5*nsize_sb_now)); iwork_lapack = 0
    allocate(ifail_lapack(nsize_sb_now)); ifail_lapack=0
    allocate(ap(nsize_mt_now*kimg));ap=0.d0
    allocate(bp(nsize_mt_now*kimg));bp=0.d0

    if(ipridavidson >=3 ) then
       write(nfout,*) "debug(zhpgvx) i,w1hw2,w1sw2"
       do ib=1,nsize_sb_now
          i = ib*(ib-1)/2 + ib
          write(nfout,*) ib,w1hw2(2*i-1),w1sw2(2*i-1)
       end do
       ap = w1sw2
       write(nfout,*) "debug(zhpgvx) i,ap,w1sw2"
       do i=1,nsize_mt_now
          write(nfout,*) i,ap(2*i-1),w1sw2(2*i-1)
       end do
       call zhpevx('N','A','U',nsize_sb_now,ap,vl,vu,il,iu,abstol,m &
    &          ,eig,vec,nsize_sb_now &
    &          ,work_lapack,rwork_lapack,iwork_lapack,ifail_lapack,info)
       write(nfout,*) "debug(zhpgvx) eig of w1sw2"
       do ib=1,nsize_sb_now
          write(nfout,*) ib,eig(ib)
       end do
    end if

    ap = w1hw2
    bp = w1sw2

    call zhpgvx(ITYPE,JOBZ,RANGE,UPLO,nsize_sb_now,ap,bp &
    &          ,vl,vu,il,iu,abstol,m,eig,vec,nsize_sb_now &
    &          ,work_lapack,rwork_lapack,iwork_lapack,ifail_lapack,info)

    if(ipridavidson >=2 .and. info/=0) then
       write(nfout,*) "debug(zhpgvx) info=",info
       write(nfout,*) "debug(zhpgvx) ifail"
       write(nfout,'(8(1x,i3))') ifail_lapack
       write(nfout,*) "debug(zhpgvx) eig"
       write(nfout,'(8(1x,f10.5))') eig
    end if

    deallocate(work_lapack)
    deallocate(rwork_lapack)
    deallocate(iwork_lapack)
    deallocate(ifail_lapack)
    deallocate(ap,bp)

    if(info/=0) then
       !!write(nfout,*) "zhpgvx: info=",info
       !!stop 'error in zhpgvx_driver'
       ierr = 1
    else
       ierr = 0
    end if

  end subroutine zhpgvx_driver

#ifdef _USE_SCALAPACK_
  subroutine pzhegvx_driver(eig, vec, w1hw2, w1sw2, ierr)
    real(kind=DP), intent(out),  dimension(nsize_sb_now)           :: eig
    real(kind=DP), intent(out),  dimension(nsize_sb_now*kimg, neg) :: vec
    real(kind=DP), intent(inout),dimension(nsize_mt_now*kimg)      :: w1hw2, w1sw2
    integer, intent(out) :: ierr
    integer :: itype
    integer :: il,iu
    integer, allocatable, dimension(:) :: ifail
    real(kind=DP) :: vl, vu, abstol
    integer :: info, m
    real(kind=DP), external :: dlamch
! ==============================================================================
!!$    complex(kind=16), allocatable, dimension(:,:) :: ap, bp, z
    complex(kind=CMPLDP), allocatable, dimension(:,:) :: ap, bp, z
    integer :: nz, count
    real(kind=DP) :: orfac = -1.0d0
!!$    complex(kind=16), allocatable, dimension(:) :: work, rwork
    complex(kind=CMPLDP), allocatable, dimension(:) :: work, rwork
    integer,          allocatable, dimension(:) :: iwork, iclustr
    real(kind=DP),    allocatable, dimension(:) :: gap
    integer :: i, j
! ==============================================================================

    abstol = 2*dlamch('S')
    il=1; iu=neg
    itype = 1

    allocate(ap(nrowa,ncola))
    allocate(bp(nrowa,ncola))
    allocate(z (nrowa,ncola))

    count = 1
    do j = 1, nsize_sb_now
      do i = 1, nsize_sb_now
        if(i .le. j) then
          call pzelset(ap,i,j,desca,w1hw2(count))
          call pzelset(bp,i,j,desca,w1sw2(count))
          count = count + 2
        endif
      enddo
    enddo

    allocate(work   (lwork2))
    allocate(rwork   (lrwork2))
    allocate(iwork  (liwork2))
    allocate(ifail  (nsize_sb_now))
    allocate(iclustr(2*nprow*npcol))
    allocate(gap    (nprow*npcol))

    call pzhegvx(itype,'V','I','U',nsize_sb_now, &
                 ap,1,1,desca, &
                 bp,1,1,desca, &
                 vl,vu,il,iu,abstol,m,nz,eig,orfac, &
                 z,1,1,desca,&
                 work,lwork2,rwork,lrwork2,iwork,liwork2, &
                 ifail,iclustr,gap,info)

    do j = 1, neg
      do i = 1, nsize_sb_now
        call pzelget('A',' ',vec(2*i-1,j),z,i,j,desca)
      enddo
    enddo

    if(ipridavidson >=2 .and. info/=0) then
       write(nfout,*) "debug(dspgvx) info=",info
       write(nfout,*) "debug(dspgvx) ifail"
       write(nfout,'(8(1x,i3))') ifail
       write(nfout,*) "debug(dspgvx) eig"
       write(nfout,'(8(1x,f10.5))') eig
    end if

    deallocate(ap)
    deallocate(bp)
    deallocate(z)
    deallocate(work)
    deallocate(rwork)
    deallocate(iwork)
    deallocate(iclustr)
    deallocate(gap)
    deallocate(ifail)

    if(info/=0) then
       ierr = 1
    else
       ierr = 0
    end if

  end subroutine pzhegvx_driver

! === DEBUG by tkato 2012/12/20 ================================================
  subroutine create_scalapack_context()
    integer :: myrank, nprocs
    integer, allocatable :: usermap(:,:)
    integer :: i, j

    call blacs_setup(myrank, nprocs)

    nprow = dav_nprow
    npcol = dav_npcol
    if(dav_divide_square > 0) then
       if(nprow == 0 .or. npcol == 0) then
          nprow = int(sqrt(real(nrank_e*nrank_g)))
          npcol = nprow
       end if
    else
       if(nprow == 0 .or. npcol == 0) then
          nprow = nrank_g
          npcol = nrank_e
       end if
    end if
    if(nprow*npcol /= nrank_g*nrank_e) then
       nprow = nrank_g
       npcol = nrank_e
    end if

    write(nfout, '("set nprow, npcol = ", 4i5)') nprow, npcol

    allocate(usermap(nprow,npcol))
    do j = 1, npcol
       do i = 1, nprow
          usermap(i,j) = myrank_k*nrank_e*nrank_g + (i - 1)*npcol + j - 1
       end do
    end do

    call blacs_get(-1, 0, ctxt)
    call blacs_gridmap(ctxt, usermap, nprow, nprow, npcol)
    call blacs_gridinfo(ctxt, nprow, npcol, myprow, mypcol)

    deallocate(usermap)
  end subroutine create_scalapack_context

  subroutine initialize_scalapack()
    integer :: nb, nb_max
    integer, external :: numroc, iceil ! External Functions
    integer :: nn, np0, mq0, nq0
    integer :: info

    nb = dav_block_size
    nb_max = min((nsize_sb_now - 1)/nprow + 1, (nsize_sb_now - 1)/npcol + 1)
    if(nb > nb_max) nb = nb_max

    write(nfout, '("set nb_max, nb = ", 4i5)') nb_max, nb

    nrowa = numroc(nsize_sb_now,nb,myprow,0,nprow) ! number of row    elements per proc
    ncola = numroc(nsize_sb_now,nb,mypcol,0,npcol) ! number of column elements per proc
    call descinit(desca,nsize_sb_now,nsize_sb_now,nb,nb,0,0,ctxt,nrowa,info)

    ! For pdsygvx
    nn      = max(max(nsize_sb_now,nb),2)
    np0     = numroc(nn,nb,0,0,nprow)
    mq0     = numroc(max(max(neg,nb),2),nb,0,0,npcol)
    lwork1  = 5*nsize_sb_now + max(5*nn,np0*mq0+2*nb*nb) &
            + iceil(neg,nprow*npcol)*nn + (neg-1)*nsize_sb_now
    liwork1 = 6*max(max(nsize_sb_now,nprow*npcol+1),4)

    ! For pzhegvx
    np0     = numroc(nsize_sb_now,nb,0,0,nprow)
    nq0     = numroc(nsize_sb_now,nb,0,0,npcol)
    lwork2  = nsize_sb_now + (np0+nq0+nb)*nb
    nn      = max(max(nsize_sb_now,nb),2)
    np0     = numroc(nn,nb,0,0,nprow)
    mq0     = numroc(max(max(neg,nb),2),nb,0,0,npcol)
    lrwork2 = 4*nsize_sb_now + max(5*nn,np0*mq0+2*nb*nb) &
            + iceil(neg,nprow*npcol)*nn + (neg-1)*nsize_sb_now
    liwork2 = 6*max(max(nsize_sb_now,nprow*npcol+1),4)
  end subroutine initialize_scalapack

  subroutine destroy_scalapack_context()
     call blacs_gridexit(ctxt)
  end subroutine destroy_scalapack_context
! ==============================================================================

#endif
#endif
end module m_ES_WF_by_Davidson

