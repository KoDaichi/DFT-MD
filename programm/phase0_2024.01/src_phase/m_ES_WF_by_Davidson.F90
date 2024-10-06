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
                                  , dav_block_size
  use m_Files,               only : nfout
  use m_Timing,              only : tstatc0_begin, tstatc0_end
!fj$$F  use m_FFT,                 only : nfft,fft_box_size_WF, m_FFT_Vlocal_W, m_FFT_WF
  use m_FFT,                 only : nfft,fft_box_size_WF
  use m_FFT,                 only : m_FFT_Vlocal_W, m_FFT_WF
  use m_Kpoints,             only : kv3,vkxyz, k_symmetry
  use m_PlaneWaveBasisSet,   only : kg1, iba, igf, nbase, m_pwBS_kinetic_energies
  use m_Electronic_Structure,only : zaj_l, neordr, nrvf_ordr, eko_l, vlhxcQ &
#ifdef SAVE_FFT_TIMES
 &                                , status_saved_phifftr &
#endif
       &                          , occup_l&
       &                          , fsr_l,fsi_l, vnlph_l, vlhxc_l
  use m_Electronic_Structure,only : m_ES_Vlocal_in_Rspace &
 &                                , m_ES_WF_in_Rspace &
 &                                , m_ES_wd_zaj_small_portion &
 &                                , m_ES_decide_precon_factor &
 &                                , m_ES_wd_eko
  use m_ES_ortho            ,only : np_g1k_x
  use m_ES_ortho            ,only : m_ES_W_transpose_r                              &
 &                                , m_ES_W_transpose_back_r
  use m_ES_nonlocal         ,only : m_ES_Vnonlocal_W                                  &
 &                                , m_ES_betar_dot_WFs_4_each_k                       &
 &                                , m_ES_alloc_scss_etc                               &
 &                                , m_ES_dealloc_scss_etc
#endif
  use m_Ionic_System,       only : natm, iwei, ityp, ntyp
  use m_PseudoPotential,    only : ilmt,nlmta,lmta,q,dion &
       &                         , lmtt,ltp,mtp &
       &                         , m_PP_include_vanderbilt_pot &
       &                         , ipaw,dion_paw
  use m_NonLocal_Potential, only : snl
#ifdef NEC_TIMER
  use nec_timer
#endif

! ==================================== added by K. Tagami ============== 11.0
  use m_Const_Parameters,    only : CMPLDP, Neglected, BuiltIn
  use m_Control_Parameters,  only : noncol, ndim_spinor, ndim_chgpot, ndim_magmom, &
       &                            SpinOrbit_mode, sw_hubbard
  use m_PseudoPotential,     only : q_noncl, dion_scr_noncl
  use m_Electronic_Structure,  only : m_ES_Vlocal_in_Rspace_noncl, &
       &                              m_ES_sort_eigen_vals_noncl

  use m_FFT,                 only : m_FFT_Vlocal_W_noncl
! ====================================================================== 11.0
  use mpi


  implicit none

  integer :: nsize_subspace, nsize_matrix
  integer :: nsize_sb_now, nsize_mt_now, nsize_mt_old
  real(kind=DP), allocatable, dimension(:) :: w1hw2,w1hw2_mpi
  real(kind=DP), allocatable, dimension(:) :: w1sw2,w1sw2_mpi
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
!!$    real(kind=DP), dimension(kg1) :: ekin,vnldi,hdiag,sdiag
    real(kind=DP), allocatable, dimension(:) :: ekin,vnldi,hdiag,sdiag
    real(kind=DP) :: vlhxc0
    logical :: frestart
    integer :: max_idavid, max_itot
    integer :: idavid_now, itot_now, ipri0
! ========================================================================================
    integer :: n_unconv
! ========================================================================================

    max_idavid = 1
    max_itot = 1

    allocate(ekin(kg1),vnldi(kg1),hdiag(kg1),sdiag(kg1))
    call m_ES_alloc_scss_etc()
    allocate(afft(nfft)); allocate(bfft(nfft))

    call allocate_matrix
#ifdef _USE_SCALAPACK_
! === DEBUG by tkato 2012/12/20 ================================================
    if(sw_dav_scalapack == ON) call create_scalapack_context()
! ==============================================================================
#endif
    do ispin = 1, nspin, (af+1)
       call m_ES_Vlocal_in_Rspace(ispin,afft)      ! (ptfft1) vlhxc_l->afft
       call vlhxc_l_zero_term(vlhxc0,ispin)        ! vlhxc_l -> vlhxc0
       do ik = ispin, kv3-nspin+ispin, nspin
          if(map_k(ik) /= myrank_k) cycle          ! MPI
          iksnl = (ik-1)/nspin + 1
          call allocate_t_matrix(ik) ! -> np_g1k_x
          call m_pwBS_kinetic_energies(ik,vkxyz,ekin) ! (diakin) ->ekin
          if(precon==ON) then
             call Vnonlocal_Diagonal_part(ispin,ik,iksnl,vnldi)
             hdiag(1:iba(ik)) = ekin(1:iba(ik)) + vlhxc0 + vnldi(1:iba(ik))

! ======================================= modified by K. Tagami ====== 11.0
!             call S_Diagonal_part(ik,iksnl,sdiag)
             call S_Diagonal_part(ispin,ik,iksnl,sdiag)
! ==================================================================== 11.0

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
                   call m_ES_Vnonlocal_W(ik,iksnl,ispin,switch_of_eko_part=OFF) ! -> vnlph_l
                else
                   call m_ES_Vnonlocal_W(ik,iksnl,ispin,switch_of_eko_part=ON) ! -> vnlph_l
                end if
                call build_subspace(ik,ekin,hdiag,sdiag,afft,bfft,idavid,precon) ! -> zajold_l
                call allreduce_fs(ik,idavid) ! -> fsr,fsi
                call evolve_WFs_in_subspace&     !-(m_ES_WF_by_Davidson)
                  &(ik,ispin,ekin,afft,bfft,idavid,frestart) !-> zaj_l
                if(ik==1.and.ipridavidson>= 2) &
                  & call m_ES_wd_zaj_small_portion(nfout,ik," -- after davidson --",21)
                call m_ES_betar_dot_WFs_4_each_k(nfout,ik)   ! -> fsr_l,fsi_l
                if(frestart) exit David_Loop
                if(eigenvalues_are_converged(n_unconv)) exit Loop
                if(nsize_sb_now + n_unconv > max_subspace_size) exit David_Loop
             end do David_Loop
          end do Loop
          if(ipridavidson>=2) then
             write(nfout,'("Davidson: ik=",i5," itot=",i5," idavid=",i5," subspace=",i5)') ik, itot_now, idavid_now, nsize_sb_now
          end if
          call deallocate_t_matrix
       enddo      ! k-point loop
    enddo      ! spin loop
!fj$$#endif
#ifdef _USE_SCALAPACK_
! === DEBUG by tkato 2012/12/20 ================================================
    if(sw_dav_scalapack == ON) call destroy_scalapack_context()
! ==============================================================================
#endif
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
    if(ipri0 >= 2) call m_ES_wd_eko(nfout,mode=SCF)
    deallocate(bfft);   deallocate(afft)
    call m_ES_dealloc_scss_etc()

    deallocate(ekin,vnldi,hdiag,sdiag)

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

    logical function eigenvalues_are_converged(n_unconv)
       integer, intent(out) :: n_unconv
       integer :: ib

       n_unconv = 0
       eigenvalues_are_converged = .true.
       do ib=1,neg
          if(.not.feigconv(ib)) then
             eigenvalues_are_converged = .false.
             n_unconv = n_unconv + 1
          end if
       end do
    end function eigenvalues_are_converged

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

! ================================ added by K. Tagami =================== 11.0
  subroutine m_ESdavidson_Renew_WF_noncl(nfout,precon)
    integer, intent(in) :: nfout,precon

    integer             :: ik, iksnl, switch_of_eko_part
    integer :: idavid,itot
    real(kind=DP), allocatable ::  afft_kt(:,:)
    real(kind=DP), allocatable ::  bfft_kt(:,:)
    real(kind=DP), allocatable :: vnlph_noncl(:,:,:,:)

    real(kind=DP), allocatable :: ekin(:)

    real(kind=DP) :: vlhxc0( ndim_chgpot )

#ifdef USE_COMPLEX_VNLDI
    complex(kind=CMPLDP), allocatable :: hdiag(:,:), sdiag(:,:)
    complex(kind=CMPLDP), allocatable :: vnldi_noncl(:,:)
#else
    real(kind=DP), allocatable :: hdiag(:,:), sdiag(:,:)
    real(kind=DP), allocatable :: vnldi_noncl(:,:)
#endif

    logical :: frestart
    integer :: max_idavid, max_itot
    integer :: idavid_now, itot_now, ipri0
    integer :: n_unconv

    integer :: is1, is2, istmp, k2

    integer :: iitmp, jjtmp

    max_idavid = 1
    max_itot = 1

!$$#ifndef PARA3D
    allocate( ekin(kg1) )
    allocate( vnldi_noncl(kg1,ndim_spinor) )
    allocate( hdiag(kg1,ndim_spinor) ); hdiag = 0.0d0
    allocate( sdiag(kg1,ndim_spinor) ); sdiag = 0.0d0
    call m_ES_alloc_scss_etc()
!$$#endif
    allocate(afft_kt(nfft,ndim_chgpot)); afft_kt = 0.0d0
    allocate(bfft_kt(nfft,ndim_spinor)); bfft_kt = 0.0d0
    allocate(vnlph_noncl(kg1,np_e,kimg,ndim_spinor)); vnlph_noncl = 0.0d0

    call allocate_matrix_noncl
!$$#ifndef PARA3D

    call m_ES_Vlocal_in_Rspace_noncl( afft_kt )    ! (ptfft1) vlhxc_l->afft
    call vlhxc_l_zero_term_noncl( vlhxc0 )        ! vlhxc_l -> vlhxc0
#ifdef _USE_SCALAPACK_
! === DEBUG by tkato 2012/12/20 ================================================
    if(sw_dav_scalapack == ON) call create_scalapack_context()
! ==============================================================================
#endif

    Do ik=1, kv3, ndim_spinor
       if ( map_k(ik) /= myrank_k )  cycle          ! MPI
       iksnl = (ik-1)/ndim_spinor + 1

       call allocate_t_matrix(ik) ! -> np_g1k_x
       call m_pwBS_kinetic_energies( ik, vkxyz, ekin ) ! (diakin) ->ekin

       if ( precon==ON ) then
          vnldi_noncl = 0.0d0;   sdiag = 1.0d0

          Do is1=1, ndim_spinor
             Do is2=1, ndim_spinor
                istmp = ( is1-1 )*ndim_spinor +is2

                call Vnonlocal_Diagonal_part_noncl( istmp, ik, iksnl, &
                     &                              vnldi_noncl(:,is1) )
             End do
          End do

          Do is1=1, ndim_spinor
             hdiag(1:iba(ik),is1) = ekin(1:iba(ik)) &
                  &                + vlhxc0(1) + vnldi_noncl(1:iba(ik),is1)

          End do
          Do is1=1, ndim_spinor
             Do is2=1, ndim_spinor
                istmp = ( is1-1 )*ndim_spinor +is2

                call S_Diagonal_part_noncl( istmp, ik, iksnl, sdiag(:,is1) )
             End do
          End do

       end if

       max_itot=0;     max_idavid=0
       Loop: do itot=1,max_iter_david
         itot_now = itot
         if ( itot>max_itot ) max_itot=itot
         feigconv = .false.

         David_Loop: do idavid=1,ndavid
           idavid_now = idavid
           if ( idavid>max_idavid ) max_idavid=idavid

           vnlph_noncl = 0.0d0
           Do is1=1, ndim_spinor
             Do is2=1, ndim_spinor
               istmp = ( is1-1 )*ndim_spinor +is2
               k2 = ik +is2 -1

               if ( precon==ON ) then
                  call m_ES_Vnonlocal_W( k2, iksnl, istmp, switch_of_eko_part=OFF )
                                                          ! -> vnlph_l
               else
                  call m_ES_Vnonlocal_W( k2, iksnl, istmp, switch_of_eko_part=ON )
                                                          ! -> vnlph_l
               endif
               vnlph_noncl(:,:,:,is1) = vnlph_noncl(:,:,:,is1) &
                  &                     + vnlph_l(:,:,:)
             End do
           End do

           call build_subspace_noncl( ik, ekin, hdiag, sdiag, afft_kt, bfft_kt, &
	&                             vnlph_noncl, idavid, precon )    ! -> zajold_l

           call allreduce_fs_noncl( ik, idavid )     ! -> fsr,fsi

           call evolve_WFs_in_subspace_noncl( ik, ekin, afft_kt, bfft_kt,&
	&                                     idavid, frestart ) !-> zaj_l

           if ( ik==1 .and. ipridavidson>= 2 ) then
             Do is1=1, ndim_spinor
               call m_ES_wd_zaj_small_portion( nfout, ik+is1-1, &
	&                                      " -- after davidson --",21 )
             End do
           endif
           Do is1=1, ndim_spinor
              call m_ES_betar_dot_WFs_4_each_k( nfout,ik+is1-1 )
                                             ! -> fsr_l,fsi_l
           End do

           if (frestart) exit David_Loop
           if (eigenvalues_are_converged(n_unconv)) exit Loop
           if (nsize_sb_now + n_unconv > max_subspace_size) exit David_Loop
         end do David_Loop
       end do Loop

!       stop

      if (ipridavidson>=2) then
          write(nfout,'("Davidson: ik=",i5," itot=",i5," idavid=",i5, &
	&               " subspace=",i5)') ik, itot_now, idavid_now, nsize_sb_now
       end if
       call deallocate_t_matrix

    enddo      ! k-point loop

!$$#endif
#ifdef _USE_SCALAPACK_
! === DEBUG by tkato 2012/12/20 ================================================
    if(sw_dav_scalapack == ON) call destroy_scalapack_context()
! ==============================================================================
#endif

    call deallocate_matrix_noncl
    if(ipridavidson>=2) then
       write(nfout,'("Davidson: max_itot=",i5," max_idavid=",i5)') max_itot,max_idavid
    end if


! ****************************************** This is important ******
    call m_ES_sort_eigen_vals_noncl()
! *******************************************************************


!$$#ifndef PARA3D
    call get_ipri0(ipridavidson,ipri0)
    if(ipri0 >= 2) call m_ES_wd_eko(nfout,mode=SCF)
!$$#endif
    deallocate(bfft_kt);   deallocate(afft_kt)
    call m_ES_dealloc_scss_etc()

    deallocate(ekin,vnldi_noncl,hdiag,sdiag)

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

!$$#ifndef PARA3D
    logical function eigenvalues_are_converged(n_unconv)
       integer, intent(out) :: n_unconv
       integer :: ib

       n_unconv = 0
       eigenvalues_are_converged = .true.
       do ib=1,neg
          if(.not.feigconv(ib)) then
             eigenvalues_are_converged = .false.
             n_unconv = n_unconv + 1
          end if
       end do
    end function eigenvalues_are_converged
!$$#endif

  end subroutine m_ESdavidson_Renew_WF_noncl
! =========================================================== 11.0

  subroutine allocate_matrix
    nsize_subspace = neg*ndavid
    nsize_matrix =  max_subspace_size*(max_subspace_size+1)/2
    allocate(feigconv(neg))
    allocate(ibover(nsize_subspace))
    allocate(fsr(neg,nlmta,ndavid))
    if(.not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
       allocate(fsi(neg,nlmta,ndavid))
    end if
    allocate(zaj_l_backup(kg1,np_e,kimg)) ! MPI
  end subroutine allocate_matrix

! ================================ added by K. Tagami =============== 11.0
  subroutine allocate_matrix_noncl
    nsize_subspace = neg*ndavid
    nsize_matrix =  max_subspace_size*(max_subspace_size+1)/2
    allocate(feigconv(neg))
    allocate(ibover(nsize_subspace))

    allocate( fsr_noncl( neg, nlmta, ndavid, ndim_spinor ) )

    if(.not.(kv3/ndim_spinor==1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
       allocate( fsi_noncl( neg, nlmta, ndavid, ndim_spinor ) )
    end if

    allocate( zaj_l_backup_noncl( kg1, np_e, kimg, ndim_spinor ))

  end subroutine allocate_matrix_noncl
! ======================================================================= 11.0

  subroutine deallocate_matrix
    deallocate(feigconv)
    deallocate(ibover)
    deallocate(fsr)
    if(.not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
       deallocate(fsi)
    end if
    deallocate(zaj_l_backup)
  end subroutine deallocate_matrix

! ================================ added by K. Tagami =============== 11.0
  subroutine deallocate_matrix_noncl
    deallocate(feigconv)
    deallocate(ibover)
    deallocate(fsr_noncl)
    if (allocated(fsi_noncl)) deallocate(fsi_noncl)
    deallocate(zaj_l_backup_noncl)
  end subroutine deallocate_matrix_noncl
! ======================================================================= 11.0

  subroutine allocate_t_matrix(ik)
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

! ================================ modified by K. Tagami ============== 11.0
!    allocate(zat_l(np_g1k_x,neg,kimg,ndavid)) ! MPI
!    allocate(zah_l(np_g1k_x,neg,kimg)) ! MPI

    if ( noncol ) then
       allocate(zat_l_noncl(np_g1k_x,neg,kimg,ndavid,ndim_spinor)) ! MPI
       allocate(zah_l_noncl(np_g1k_x,neg,kimg,ndim_spinor)) ! MPI
    else
       allocate(zat_l(np_g1k_x,neg,kimg,ndavid)) ! MPI
       allocate(zah_l(np_g1k_x,neg,kimg)) ! MPI
    endif
! ======================================================================== 11.0

    allocate(w1hw2(nsize_matrix*kimg_t))
    allocate(w1sw2(nsize_matrix*kimg_t))
    if(npes>1) then
       allocate(w1hw2_mpi(nsize_matrix*kimg_t))
       allocate(w1sw2_mpi(nsize_matrix*kimg_t))
    end if
  end subroutine allocate_t_matrix

  subroutine deallocate_t_matrix
! =============================== modified by K. Tagami ============ 11.0
!    deallocate(zat_l) ! MPI
!    deallocate(zah_l) ! MPI

    if ( noncol ) then
       deallocate(zat_l_noncl) ! MPI
       deallocate(zah_l_noncl) ! MPI
    else
       deallocate(zat_l) ! MPI
       deallocate(zah_l) ! MPI
    endif
! ====================================================================== 11.0
    deallocate(w1hw2)
    deallocate(w1sw2)
    if(npes>1) then
       deallocate(w1hw2_mpi)
       deallocate(w1sw2_mpi)
    end if
  end subroutine deallocate_t_matrix

  subroutine allreduce_fs(ik,idavid)
    integer, intent(in) :: ik,idavid
    integer :: ib,ib1,iter,kimg_t
    real(kind=DP), allocatable, dimension(:,:) :: fs_mpi,fs_mpi2

    allocate(fs_mpi(neg,nlmta))
    allocate(fs_mpi2(neg,nlmta))
    if(idavid>1) then
       call m_ES_betar_dot_WFs_4_each_k(nfout,ik)   ! -> fsr_l,fsi_l
    end if
    fs_mpi=0.d0
    do ib1 = ista_e, iend_e, istep_e     ! MPI
       ib=map_z(ib1)
       fs_mpi(ib1,1:nlmta) = fsr_l(ib,1:nlmta,ik)
    end do
    fs_mpi2=0.d0
    call mpi_allreduce(fs_mpi,fs_mpi2,neg*nlmta &
      & ,mpi_double_precision,mpi_sum &
      & ,mpi_k_world(myrank_k),ierr)       ! MPI
    fsr(1:neg,1:nlmta,idavid) = fs_mpi2(1:neg,1:nlmta)
    if(.not. k_symmetry(ik) == GAMMA) then
       fs_mpi=0.d0
       do ib1 = ista_e, iend_e, istep_e     ! MPI
          ib=map_z(ib1)
          fs_mpi(ib1,1:nlmta) = fsi_l(ib,1:nlmta,ik)
       end do
       fs_mpi2=0.d0
       call mpi_allreduce(fs_mpi,fs_mpi2,neg*nlmta &
         & ,mpi_double_precision,mpi_sum &
         & ,mpi_k_world(myrank_k),ierr)       ! MPI
       fsi(1:neg,1:nlmta,idavid) = fs_mpi2(1:neg,1:nlmta)
    end if
    deallocate(fs_mpi)
    deallocate(fs_mpi2)
  end subroutine allreduce_fs

! =================================== added by K. Tagami ================ 11.0
  subroutine allreduce_fs_noncl( ik, idavid )
    integer, intent(in) :: ik,idavid

    integer :: ib, ib1, iter, is, kimg_t
    real(kind=DP), allocatable, dimension(:,:,:) :: fs_mpi,fs_mpi2

    allocate(fs_mpi( neg,nlmta,ndim_spinor));  fs_mpi  = 0.0d0
    allocate(fs_mpi2(neg,nlmta,ndim_spinor));  fs_mpi2 = 0.0d0

    if ( idavid>1 ) then
       Do is=1, ndim_spinor
         call m_ES_betar_dot_WFs_4_each_k( nfout,ik+is-1 )   ! -> fsr_l,fsi_l
       End do
    end if

    do ib1 = ista_e, iend_e, istep_e     ! MPI
       ib = map_z(ib1)
       Do is=1, ndim_spinor
         fs_mpi( ib1, 1:nlmta, is ) = fsr_l( ib, 1:nlmta, ik+is-1 )
       End do
    end do
    call mpi_allreduce( fs_mpi, fs_mpi2, neg*nlmta*ndim_spinor, &
      &                 mpi_double_precision, mpi_sum, mpi_k_world(myrank_k), ierr)

!
    Do is=1, ndim_spinor
      fsr_noncl( 1:neg, 1:nlmta, idavid, is ) = fs_mpi2( 1:neg, 1:nlmta, is )
    End do
!
    if ( .not. k_symmetry(ik) == GAMMA ) then
       fs_mpi=0.d0;  fs_mpi2 = 0.0d0
       do ib1 = ista_e, iend_e, istep_e     ! MPI
          ib=map_z(ib1)
          Do is=1, ndim_spinor
             fs_mpi( ib1, 1:nlmta, is ) = fsi_l( ib, 1:nlmta, ik+is-1 )
          End do
       end do

       call mpi_allreduce( fs_mpi, fs_mpi2, neg*nlmta*ndim_spinor, &
         &                 mpi_double_precision, mpi_sum, mpi_k_world(myrank_k), ierr )
       Do is=1, ndim_spinor
         fsi_noncl( 1:neg, 1:nlmta, idavid, is ) = fs_mpi2( 1:neg, 1:nlmta, is )
       End do
    end if

    deallocate(fs_mpi);   deallocate(fs_mpi2)
  end subroutine allreduce_fs_noncl
! ===================================================================== 11.0

  subroutine build_subspace(ik,ekin,hdiag,sdiag,afft,bfft,idavid,precon)
    integer, intent(in) :: ik,idavid,precon
    real(kind=DP), intent(in)  :: ekin(kg1),hdiag(kg1),sdiag(kg1)
    real(kind=DP), intent(in)  :: afft(nfft)
    real(kind=DP), intent(out) :: bfft(nfft)

    integer :: iter,ib1,ib,ii,i1
    real(kind=DP) :: dr1,dr2,di1
    real(kind=DP) :: p(kg1)
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
    do ib1 = ista_e, iend_e, istep_e     ! MPI
       if(feigconv(ib1)) cycle
       call m_ES_WF_in_Rspace(ik,ib1,bfft)!(swffft)
       call m_FFT_Vlocal_W(afft,bfft)    ! (afft, bfft)-> (bfft)
       call m_FFT_WF(ELECTRON,nfout,bfft,DIRECT,ON)
       ib = map_z(ib1)                 ! MPI
#ifdef SAVE_FFT_TIMES
       if(sw_save_fft == ON) status_saved_phifftr(ib,ik) = OLD
#endif
       if(precon==ON) then
          call decide_precon_factor_david(ik,hdiag,sdiag,eko_l(ib,ik),p) ! -> p
          if(kimg == 1) then
             do ii=1,iba(ik)
                i1  = igf(nbase(ii,ik))
                dr1 = zaj_l(ii,ib,ik,kimg)
                dr2 = bfft(i1)*denom
                zajr= (ekin(ii)-hdiag(ii))*dr1+dr2+vnlph_l(ii,ib,kimg)
                zaj_l(ii,ib,ik,kimg) = p(ii)*zajr
             enddo
          else
             do ii=1,iba(ik)
                i1  = igf(nbase(ii,ik))
                dr1  = zaj_l(ii,ib,ik,1)
                di1  = zaj_l(ii,ib,ik,kimg)
                dr2  = ekin(ii) - hdiag(ii)
                zajr= dr2*dr1+bfft(2*i1-1)*denom + vnlph_l(ii,ib,1)
                zaji= dr2*di1+bfft(2*i1)*denom + vnlph_l(ii,ib,2)
                zaj_l(ii,ib,ik,1   )=p(ii)*zajr
                zaj_l(ii,ib,ik,kimg)=p(ii)*zaji
             enddo
          endif
       else
          if(kimg == 1) then
             do ii=1,iba(ik)
                i1  = igf(nbase(ii,ik))
                dr1 = zaj_l(ii,ib,ik,kimg)
                dr2 = bfft(i1)*denom
                zajr= ekin(ii)*dr1+dr2+vnlph_l(ii,ib,kimg)
                zaj_l(ii,ib,ik,kimg) = zajr
             enddo
          else
             do ii=1,iba(ik)
                i1  = igf(nbase(ii,ik))
                dr1  = zaj_l(ii,ib,ik,1)
                di1  = zaj_l(ii,ib,ik,kimg)
                zajr= ekin(ii)*dr1+bfft(2*i1-1)*denom + vnlph_l(ii,ib,1)
                zaji= ekin(ii)*di1+bfft(2*i1)*denom + vnlph_l(ii,ib,2)
                zaj_l(ii,ib,ik,1   )=zajr
                zaj_l(ii,ib,ik,kimg)=zaji
             enddo
          endif
       end if

       ! Normalize
       norm = 0.d0
       if(kimg==1) then
          do ii=1,iba(ik)
             zajr = zaj_l(ii,ib,ik,kimg)
             norm = norm + zajr*zajr
          end do
          norm = 1.d0/sqrt(norm)
          do ii=1,iba(ik)
             zaj_l(ii,ib,ik,1) = zaj_l(ii,ib,ik,1)*norm
          end do
       else
          do ii=1,iba(ik)
             zajr = zaj_l(ii,ib,ik,1   )
             zaji = zaj_l(ii,ib,ik,kimg)
             norm = norm + zajr*zajr+zaji*zaji
          end do
          norm = 1.d0/sqrt(norm)
          do ii=1,iba(ik)
             zaj_l(ii,ib,ik,1)= zaj_l(ii,ib,ik,1)*norm
! === DEBUG by tkato 2011/07/13 ================================================
!            zaj_l(ii,ib,ik,2)= zaj_l(ii,ib,ik,1)*norm
             zaj_l(ii,ib,ik,2)= zaj_l(ii,ib,ik,2)*norm
! ==============================================================================
          enddo
       end if
    end do
  contains
    subroutine decide_precon_factor_david(ik,hdiag,sdiag,eig,p)
      integer, intent(in) :: ik
      real(kind=DP), intent(in) :: hdiag(kg1),sdiag(kg1),eig
      real(kind=DP), intent(out) :: p(kg1)

      integer :: ii
      real(kind=DP) :: denom

      do ii=1,iba(ik)
         denom = hdiag(ii)-eig*sdiag(ii)
         if(abs(denom) < eps_david) then
            denom = sign(eps_david,denom)
         end if
         p(ii) = 1.d0/denom
      end do

    end subroutine decide_precon_factor_david
  end subroutine build_subspace

! ================================ added by K. Tagami ==================== 11.0
  subroutine build_subspace_noncl( ik, ekin, hdiag, sdiag, afft_kt, bfft_kt, &
	&                          vnlph_noncl, idavid, precon )
    integer, intent(in) :: ik,idavid,precon
    real(kind=DP), intent(in)  :: ekin(kg1)


#ifdef USE_COMPLEX_VNLDI
    complex(kind=CMPLDP), intent(in)  :: hdiag(kg1,ndim_spinor), sdiag(kg1,ndim_spinor)
#else
    real(kind=DP), intent(in)  :: hdiag(kg1,ndim_spinor), sdiag(kg1,ndim_spinor)
#endif

    real(kind=DP), intent(inout)  :: afft_kt(nfft,ndim_chgpot)
    real(kind=DP), intent(out) :: bfft_kt(nfft,ndim_spinor)
    real(kind=DP), intent(in)  :: vnlph_noncl(kg1,np_e,kimg,ndim_spinor)

    integer :: iter,ib1,ib,ii,i1
    integer :: is, k1

    real(kind=DP) :: dr1,dr2,di1
    real(kind=DP) :: denom, norm, zajr, zaji

    complex(kind=CMPLDP) :: zdr2
    real(kind=DP) :: ctmp1, stmp1, ctmp2, stmp2

#ifdef USE_COMPLEX_VNLDI
    complex(kind=CMPLDP), allocatable :: p(:,:)
#else
    real(kind=DP), allocatable :: p(:,:)
#endif

    denom = 1.d0/product(fft_box_size_WF(1:3,1))
    allocate( p( kg1,ndim_spinor ) ); p = 0.0d0

    if(idavid==1) then
       return
    end if

    ! if precon == ON
    ! |R> = P*(H-Hdiag)|Psi>
    ! P = (Hdiag-e*Sdiag)^(-1)
    ! else
    ! |R> = (H-eS)|Psi>
    do ib1 = ista_e, iend_e, istep_e     ! MPI
       if ( feigconv(ib1) ) cycle

       Do is=1, ndim_spinor
         call m_ES_WF_in_Rspace( ik+is-1, ib1, bfft_kt(:,is) )!(swffft)
       End do
       call m_FFT_Vlocal_W_noncl( afft_kt, bfft_kt, ndim_chgpot, ndim_spinor )
                                           ! (afft, bfft)-> (bfft)
       Do is=1, ndim_spinor
         call m_FFT_WF( ELECTRON, nfout, bfft_kt(:,is), DIRECT, ON )
       End do

       ib = map_z(ib1)                 ! MPI
#ifdef SAVE_FFT_TIMES
       if(sw_save_fft == ON) then
          do is = 1, ndim_spinor
             k1 = ik+is-1
             status_saved_phifftr(ib,k1) = OLD
          end do
       end if
#endif
       if ( precon==ON ) then
          call decide_precon_fac_david_noncl( ik, hdiag, sdiag, eko_l(ib,ik), p )
                                                              ! -> p
#ifdef USE_COMPLEX_VNLDI
          if ( kimg == 1 ) then
          else
             Do is=1, ndim_spinor
                k1 = ik +is -1
                do ii=1,iba(ik)
                   i1  = igf( nbase(ii,ik) )

                   dr1  = zaj_l( ii, ib, k1, 1 )
                   di1  = zaj_l( ii, ib, k1, kimg )
                   zdr2  = ekin(ii) - hdiag(ii,is)

                   ctmp1 = real(zdr2)*dr1 - aimag(zdr2)*di1
                   stmp1 = aimag(zdr2)*dr1 + real(zdr2)*di1

                   zajr= ctmp1 +bfft_kt( 2*i1-1,is )*denom &
        &                        +vnlph_noncl( ii, ib, 1, is )
                   zaji= stmp1 +bfft_kt( 2*i1,  is )*denom &
        &                        +vnlph_noncl( ii, ib, 2, is )

                   ctmp2 = real(p(ii,is)) *zajr - aimag(p(ii,is))*zaji
                   stmp2 = aimag(p(ii,is)) *zajr + real(p(ii,is))*zaji

                   zaj_l( ii, ib, k1, 1    ) = ctmp2
                   zaj_l( ii, ib, k1, kimg ) = stmp2
                End do
             enddo
          endif

#else
          if ( kimg == 1 ) then
             Do is=1, ndim_spinor
                k1 = ik +is -1
                do ii=1, iba(ik)
                   i1 = igf( nbase(ii,ik) )

	           dr1 = zaj_l( ii, ib, k1, kimg )
                   dr2 = bfft_kt( i1,is )*denom
                   zajr= ( ekin(ii) - hdiag(ii,is) )*dr1 &
	&                   + dr2 +vnlph_noncl( ii, ib, kimg, is  )
                   zaj_l( ii, ib, k1, kimg ) = p(ii,is) *zajr
                End do
             enddo
          else
             Do is=1, ndim_spinor
                k1 = ik +is -1
                do ii=1,iba(ik)
                   i1  = igf( nbase(ii,ik) )

                   dr1  = zaj_l( ii, ib, k1, 1 )
                   di1  = zaj_l( ii, ib, k1, kimg )
                   dr2  = ekin(ii) - hdiag(ii,is)
                   zajr= dr2*dr1 +bfft_kt( 2*i1-1,is )*denom &
	&                        +vnlph_noncl( ii, ib, 1, is )
                   zaji= dr2*di1 +bfft_kt( 2*i1,  is )*denom &
	&                        +vnlph_noncl( ii, ib, 2, is )
                   zaj_l( ii, ib, k1, 1    ) = p(ii,is) *zajr
                   zaj_l( ii, ib, k1, kimg ) = p(ii,is) *zaji
                End do
             enddo
          endif

#endif
       else
          if ( kimg == 1 ) then
             Do is=1, ndim_spinor
                k1 = ik +is -1
                do ii=1, iba(ik)
                   i1  = igf(nbase(ii,ik))

                   dr1 = zaj_l( ii, ib, k1, kimg )
                   dr2 = bfft_kt( i1,is )*denom
                   zajr = ekin(ii)*dr1 +dr2 +vnlph_noncl( ii, ib, kimg, is )
                   zaj_l( ii, ib, k1, kimg ) = zajr
                enddo
             End do
          else
             Do is=1, ndim_spinor
                k1 = ik +is -1
                do ii=1, iba(ik)
                   i1  = igf(nbase(ii,ik))

                   dr1  = zaj_l( ii, ib, k1, 1 )
                   di1  = zaj_l( ii, ib, k1, kimg )
                   zajr = ekin(ii)*dr1 +bfft_kt( 2*i1-1,is )*denom &
	&                              +vnlph_noncl( ii, ib, 1,is )
                   zaji = ekin(ii)*di1 +bfft_kt( 2*i1,  is )*denom &
	&                              +vnlph_noncl( ii, ib, 2,is )
                   zaj_l( ii, ib, k1, 1    ) = zajr
                   zaj_l( ii, ib, k1, kimg ) = zaji
                enddo
             End do
          endif
       end if

       ! Normalize
       norm = 0.d0
       if ( kimg==1 ) then
          Do is=1, ndim_spinor
             k1 = ik +is -1
             do ii=1, iba(ik)
                zajr = zaj_l( ii, ib, k1, kimg )
                norm = norm + zajr*zajr
             End do
          end do
          norm = 1.d0 /sqrt(norm)
          Do is=1, ndim_spinor
             k1 = ik +is -1
             do ii=1, iba(ik)
                zaj_l( ii, ib, k1, 1 ) = zaj_l( ii, ib, k1, 1 )*norm
             End do
          end do
       else
          Do is=1, ndim_spinor
             k1 = ik +is -1
             do ii=1, iba(ik)
                zajr = zaj_l( ii, ib, k1, 1    )
                zaji = zaj_l( ii, ib, k1, kimg )
                norm = norm + zajr*zajr +zaji*zaji
             end do
          end do
          norm = 1.d0 /sqrt(norm)
          Do is=1, ndim_spinor
             k1 = ik +is -1
             do ii=1,iba(ik)
                zaj_l( ii, ib, k1, 1 )    = zaj_l( ii, ib, k1, 1    )*norm
                zaj_l( ii, ib, k1, kimg ) = zaj_l( ii, ib, k1, kimg )*norm
             end do
          enddo
       end if
    end do
!
    deallocate(p)

  contains

    subroutine decide_precon_fac_david_noncl( ik, hdiag, sdiag, eig, p )
      integer, intent(in) :: ik

#ifdef USE_COMPLEX_VNLDI
      complex(kind=CMPLDP), intent(in) :: hdiag(kg1,ndim_spinor), &
           &                              sdiag(kg1,ndim_spinor)
      complex(kind=CMPLDP), intent(out) :: p(kg1,ndim_spinor)
#else
      real(kind=DP), intent(in) :: hdiag(kg1,ndim_spinor), sdiag(kg1,ndim_spinor)
      real(kind=DP), intent(out) :: p(kg1,ndim_spinor)
#endif

      real(kind=DP), intent(in) :: eig

      integer :: ii, is
      real(kind=DP) :: denom
      complex(kind=CMPLDP) :: ztmp

#ifdef USE_COMPLEX_VNLDI
      Do is=1, ndim_spinor
         do ii=1,iba(ik)
            ztmp = hdiag(ii,is)-eig*sdiag(ii,is)
            denom = real(ztmp)
            if ( abs(denom) < eps_david ) then
               if ( denom > 0.0d0 ) then
                  ztmp = eps_david
               else
                  ztmp = -eps_david
               endif
            end if
            p(ii,is) = 1.d0 / ztmp
         End do
      end do
#else
      Do is=1, ndim_spinor
         do ii=1,iba(ik)
           denom = hdiag(ii,is)-eig*sdiag(ii,is)
           if(abs(denom) < eps_david) then
              denom = sign(eps_david,denom)
           end if
           p(ii,is) = 1.d0/denom
         End do
      end do
#endif

    end subroutine decide_precon_fac_david_noncl

  end subroutine build_subspace_noncl
! ===================================================================== 11.0


  subroutine evolve_WFs_in_subspace&
       &(ik,ispin,ekin,afft,bfft,idavid,frestart)
    integer, intent(in) :: ik,ispin
    integer, intent(in) :: idavid
    logical, intent(out) :: frestart
    real(kind=DP), intent(in)  :: ekin(kg1)
    real(kind=DP), intent(in)  :: afft(nfft)
    real(kind=DP), intent(out) :: bfft(nfft)
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
    call tstatc0_begin('evolve_WFs_in_subspace (davidson) ', id_sname,1)

    call get_ipri0(ipridavidson,ipri0)

! ==================== added by K. Tagami ========== 11.0
#ifdef forsafe
    ekod = 0.0d0
#endif
! ================================================== 11.0

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


    do ib1 = 1, neg
       if(map_e(ib1) == myrank_e) eko_d(ib1) = eko_l(map_z(ib1),ik)  ! MPI
    end do
    call mpi_allreduce(eko_d,eko_d_mpi,neg,mpi_double_precision,mpi_sum &
         & ,mpi_k_world(myrank_k),ierr)       ! MPI
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
    call m_ES_W_transpose_r(.false.,ista_k,iend_k,ik,zaj_l,zat_l(1,1,1,idavid))
    do ib1 = ista_e, iend_e, istep_e     ! MPI
       call m_ES_WF_in_Rspace(ik,ib1,bfft)!(swffft)
       call m_FFT_Vlocal_W(afft,bfft)    ! (afft, bfft)-> (bfft)
       call m_FFT_WF(ELECTRON,nfout,bfft,DIRECT,ON)
       ib = map_z(ib1)                 ! MPI
       if(kimg == 1) then
          do ii=1,iba(ik)
             i1  = igf(nbase(ii,ik))
             dr1 = zaj_l(ii,ib,ik,1)
             dr2 = bfft(i1)*denom
             zaj_l(ii,ib,ik,1)= ekin(ii)*dr1+dr2
          enddo
       else
          do ii=1,iba(ik)
             i1  = igf(nbase(ii,ik))
             dr1  = zaj_l(ii,ib,ik,1)
             di1  = zaj_l(ii,ib,ik,kimg)
             zaj_l(ii,ib,ik,1)= ekin(ii)*dr1+bfft(2*i1-1)*denom
             zaj_l(ii,ib,ik,kimg)= ekin(ii)*di1+bfft(2*i1)*denom
          enddo
       endif
#ifdef SAVE_FFT_TIMES
       if(sw_save_fft == ON) status_saved_phifftr(ib,ik) = OLD
#endif
    enddo
    if(ipridavidson>=2 .and. ik==1) write(nfout,'(" !### zaj_l is new,  bfft is old")')
!( tenchi ) (zah_l <- zaj_l)
    call m_ES_W_transpose_r(.false.,ista_k,iend_k,ik,zaj_l,zah_l(1,1,1))

! (make matrix elements )
! parallel loop
    ! <n|T+Vloc|m> G-wise parallel
    do ibb2 = nsize_max_sb_now-neg+1,nsize_max_sb_now
       if(ibover(ibb2)<0) cycle
       ib2 = ibover(ibb2)
       iter2 = (ibb2-1)/neg+1
       ii2  = ibb2-neg*(iter2-1)
       ip0b = ib2*(ib2-1)/2
       do ibb1 = 1,ibb2
          if(ibover(ibb1)<0) cycle
          ib1 = ibover(ibb1)
          iter1 = (ibb1-1)/neg+1
          ii1 = ibb1-neg*(iter1-1)
          ip0 = ip0b + ib1
          if(kimg == 1) then
             w1hw2(ip0) = 0.d0
             w1sw2(ip0) = 0.d0
             do ii = 1, np_g1k(ik)            ! MPI
                hr2 = zah_l(ii,ii2,1)
                dr2 = zat_l(ii,ii2,1,iter2)
                dr1 = zat_l(ii,ii1,1,iter1)
                w1hw2(ip0) = w1hw2(ip0) + dr1*hr2
                w1sw2(ip0) = w1sw2(ip0) + dr1*dr2
             end do
          else
             if(k_symmetry(ik) == GAMMA) then
                w1hw2(ip0) = 0.d0
                w1sw2(ip0) = 0.d0
!!$                ig1 = 1;  if(mype == 0) ig1 = 2
                ig1 = 1;  if(myrank_e == 0) ig1 = 2
                do ii = ig1, np_g1k(ik)            ! MPI
                   hr2 = zah_l(ii,ii2,1) ! MPI
                   hi2 = zah_l(ii,ii2,2) ! MPI
                   dr2 = zat_l(ii,ii2,1,iter2) ! MPI
                   di2 = zat_l(ii,ii2,2,iter2) ! MPI
                   dr1 = zat_l(ii,ii1,1,iter1) ! MPI
                   di1 = zat_l(ii,ii1,2,iter1) ! MPI
                   w1hw2(ip0) =w1hw2(ip0)+(dr1*hr2+di1*hi2)*2.d0
                   w1sw2(ip0) =w1sw2(ip0)+(dr1*dr2+di1*di2)*2.d0
                end do
!!$                if(mype == 0) then
                if(myrank_e == 0) then
                   hr2 = zah_l(1,ii2,1) ! MPI
                   hi2 = zah_l(1,ii2,2) ! MPI
                   dr2 = zat_l(1,ii2,1,iter2) ! MPI
                   di2 = zat_l(1,ii2,2,iter2) ! MPI
                   dr1 = zat_l(1,ii1,1,iter1) ! MPI
                   di1 = zat_l(1,ii1,2,iter1) ! MPI
                   w1hw2(ip0) =w1hw2(ip0)+dr1*hr2+di1*hi2
                   w1sw2(ip0) =w1sw2(ip0)+dr1*dr2+di1*di2
                end if
             else
                w1hw2(2*ip0-1:2*ip0) = 0.d0
                w1sw2(2*ip0-1:2*ip0) = 0.d0
                do ii = 1, np_g1k(ik)           ! MPI
                   hr2 = zah_l(ii,ii2,1) ! MPI
                   hi2 = zah_l(ii,ii2,2) ! MPI
                   dr2 = zat_l(ii,ii2,1,iter2) ! MPI
                   di2 = zat_l(ii,ii2,2,iter2) ! MPI
                   dr1 = zat_l(ii,ii1,1,iter1) ! MPI
                   di1 = zat_l(ii,ii1,2,iter1) ! MPI
                   w1hw2(2*ip0-1) =w1hw2(2*ip0-1)+dr1*hr2+di1*hi2
                   w1hw2(2*ip0  ) =w1hw2(2*ip0  )+dr1*hi2-di1*hr2
                   w1sw2(2*ip0-1) =w1sw2(2*ip0-1)+dr1*dr2+di1*di2
                   w1sw2(2*ip0  ) =w1sw2(2*ip0  )+dr1*di2-di1*dr2
                end do
             end if
          end if
       end do
    end do
    if(ipridavidson >= 3) call wd_w1hw2(" -- w1hw2 without nl part--")
    ! <n|Vnl|m> G-wise parallel
    call add_nonlocal_part ! w1hw2 = w1hw2 + w1Vnlw2
                           ! w1sw2 = w1sw2 + w1qw2
    if(ipridavidson >= 3) call wd_w1hw2(" -- w1hw2 with nl part--")

!! (spread sum of w1hw2 and w1sw2)
    if(npes > 1) then
       w1hw2_mpi = 0.d0
       w1sw2_mpi = 0.d0
       nshift = nsize_mt_old*kimg_t
       ndata = (nsize_mt_now-nsize_mt_old)*kimg_t
       call mpi_allreduce(w1hw2(nshift+1),w1hw2_mpi,ndata,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr) ! MPI
       w1hw2(nshift+1:nshift+ndata) = w1hw2_mpi(1:ndata) ! MPI
       call mpi_allreduce(w1sw2(nshift+1),w1sw2_mpi,ndata,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr) ! MPI
       w1sw2(nshift+1:nshift+ndata) = w1sw2_mpi(1:ndata) ! MPI
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

#ifdef _USE_SCALAPACK_
    if(sw_dav_scalapack == ON) then
       call initialize_scalapack()
       if(kimg_t == 1) then
          call pdsygvx_driver(eig,vec,w1hw2,w1sw2,ierr_diag)
       else
          call pzhegvx_driver(eig,vec,w1hw2,w1sw2,ierr_diag)
       endif
! === DEBUG by tkato 2012/12/20 ================================================
!      call finalize_scalapack()
! ==============================================================================
    else
       if(kimg_t == 1) then
          call dspgvx_driver(eig,vec,w1hw2,w1sw2,ierr_diag)
       else
          call zhpgvx_driver(eig,vec,w1hw2,w1sw2,ierr_diag)
       endif
    endif
#else
    if(kimg_t == 1) then
       call dspgvx_driver(eig,vec,w1hw2,w1sw2,ierr_diag)
    else
       call zhpgvx_driver(eig,vec,w1hw2,w1sw2,ierr_diag)
    endif
#endif

    frestart = .false.
    if(ierr_diag /= 0) then
       zaj_l(:,:,ik,:) = zaj_l_backup(:,:,:)
       if(ipridavidson>=2 .and. ik==1) write(nfout,'(" !### zaj_l is new,  bfft is old")')
       do ib1 = 1, neg
          if(map_e(ib1) == myrank_e) then         ! MPI
             eko_l(map_z(ib1),ik)=eko_d(ib1)
          end if
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

!!$       if(ipridavidson >= 2) then
       if(ipri0 >= 2) then
          write(nfout,*) 'eko_d for ik = ',ik
          write(nfout,9001) (eko_d(ib),ib=1,neg)
          write(nfout,*) 'eig for ik = ',ik
          write(nfout,9001) (eig(ib),ib=1,neg)
          call wd_w1hw2(" -- after diagonalization --")
!sum eko
          dr1=0.d0;dr2=0.d0
          do ib1=1,neg
             ib1to = neordr(ib1,ik)
             if(map_e(ib1to) == myrank_e) dr1=dr1+eko_l(map_z(ib1to),ik) ! MPI
             dr2=dr2+eig(ib1)
          enddo
          call mpi_allreduce(dr1,di1,1,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr) ! MPI
          dr1 = di1  ! MPI
          write(nfout,'(" sum of eko_l, eig, abs diff =",3e25.10)') dr1,dr2,abs(dr2-dr1)
       endif
!! (subspace rotation) !!
       call subspace_rotation ! vec,zat_l -> zat_l
!( tenchi ) (zaj_l <- zat_l)
       iter = min(idavid+1,ndavid)
       call m_ES_W_transpose_back_r(.false.,ista_k,iend_k,ik,zaj_l,zat_l(1,1,1,iter))
       zaj_l_backup(:,:,:) = zaj_l(:,:,ik,:)
!! (eko_l)
       do ib1 = 1, neg
          if(map_e(ib1) == myrank_e) then         ! MPI
             eko_l(map_z(ib1),ik)=eig(ib1)
          end if
       end do
       if(ipridavidson >= 2) then
          eko2 = sum(eig(1:neg))
          write(nfout,1201) ik,eko1,ekod,eko2

          write(nfout,*) 'eko_l'
          write(nfout,9001) (eko_l(ib1,ik),ib1=1,neg)
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

! ========================== added by K. Tagami ========== 11.0
#ifdef forsafe
      integer :: ipaw_tmp
#endif
! ======================================================== 11.0

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

! ========================== added by K. Tagami =================== 11.0
#ifdef forsafe
               ipaw_tmp = ipaw(it)
#endif
! ================================================================= 11.0

               do lmt1 = 1, ilmt(it)
                  p = lmta(lmt1,ia)
                  do lmt2 = 1, ilmt(it)
                     s = lmta(lmt2,ia)
!!$                     facv   = iwei(ia)*(dion(lmt1,lmt2,it) + vlhxcQ(lmt1,lmt2,ia,ispin))
! ========================== modified by K. Tagami =================== 11.0
!                     if(ipaw(it)==0)then
#ifdef forsafe
                     if ( ipaw_tmp == 0 ) then
#else
                     if (ipaw(it)==0 ) then
#endif
! ================================================================= 11.0

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
      integer :: ib1,ib2,ibb2
!!$      real(kind=DP), dimension(np_g1k_x,neg,kimg) :: zaj_wk
      real(kind=DP), allocatable, dimension(:,:,:) :: zaj_wk
      allocate(zaj_wk(np_g1k_x,neg,kimg))

      zaj_wk(:,:,:) = 0.d0
      if(kimg==1) then
         do ib1=1,neg
            do ibb2=1,nsize_max_sb_now
               if(ibover(ibb2)<0) cycle
               ib2 = ibover(ibb2)
               iter2=(ibb2-1)/neg+1
               ii2=ibb2-neg*(iter2-1)
               do ii=1,np_g1k(ik)
                  zaj_wk(ii,ib1,kimg) = zaj_wk(ii,ib1,kimg) + zat_l(ii,ii2,kimg,iter2)*vec(ib2,ib1)
               end do
            end do
         end do
      else
         if(k_symmetry(ik) == GAMMA) then
            do ib1=1,neg
               do ibb2=1,nsize_max_sb_now
                  if(ibover(ibb2)<0) cycle
                  ib2 = ibover(ibb2)
                  iter2=(ibb2-1)/neg+1
                  ii2=ibb2-neg*(iter2-1)
                  hr2=vec(ib2,ib1)
                  do ii=1,np_g1k(ik)
                     dr1=zat_l(ii,ii2,1   ,iter2)
                     di1=zat_l(ii,ii2,kimg,iter2)
                     zaj_wk(ii,ib1,1   ) = zaj_wk(ii,ib1,1   ) + dr1*hr2
                     zaj_wk(ii,ib1,kimg) = zaj_wk(ii,ib1,kimg) + di1*hr2
                  end do
               end do
            end do
         else
            do ib1=1,neg
               do ibb2=1,nsize_max_sb_now
                  if(ibover(ibb2)<0) cycle
                  ib2 = ibover(ibb2)
                  iter2=(ibb2-1)/neg+1
                  ii2=ibb2-neg*(iter2-1)
                  hr2=vec(2*ib2-1,ib1)
                  hi2=vec(2*ib2  ,ib1)
                  do ii=1,np_g1k(ik)
                     dr1=zat_l(ii,ii2,1   ,iter2)
                     di1=zat_l(ii,ii2,kimg,iter2)
                     zaj_wk(ii,ib1,1   ) = zaj_wk(ii,ib1,1   ) + dr1*hr2 - di1*hi2
                     zaj_wk(ii,ib1,kimg) = zaj_wk(ii,ib1,kimg) + dr1*hi2 + di1*hr2
                  end do
               end do
            end do
         end if
      end if
      iter = min(idavid+1,ndavid)
      zat_l(:,:,:,iter) = zaj_wk(:,:,:)
      deallocate(zaj_wk)
    end subroutine subspace_rotation

  end subroutine evolve_WFs_in_subspace

! ==================================== added by K. Tagami ================= 11.0
  subroutine evolve_WFs_in_subspace_noncl( ik, ekin, afft_kt, bfft_kt, &
	&                                  idavid, frestart )
    integer, intent(in) :: ik
    integer, intent(in) :: idavid
    logical, intent(out) :: frestart
    real(kind=DP), intent(in)  :: ekin(kg1)
    real(kind=DP), intent(inout)  :: afft_kt(nfft,ndim_chgpot)
    real(kind=DP), intent(out) :: bfft_kt(nfft,ndim_spinor)

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
    integer :: is, k1

    integer :: ip0,ip0b,ip1,ip1b,ib1n,ib2n,ndata,nshift,kimg_t,ig1
    integer :: noffset
    integer :: nsize_max_sb_now
    integer :: ierr_diag
    integer :: id_sname = -1, ipri0

    call tstatc0_begin('evolve_WFs_in_subspace (davidson) ', id_sname,1)
    call get_ipri0(ipridavidson,ipri0)

! ==================== added by K. Tagami ========== 11.0
#ifdef forsafe
    ekod = 0.0d0
#endif
! ================================================== 11.0

    denom = 1.d0 /product(fft_box_size_WF(1:3,1))
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

    do ib1 = 1, neg
       if(map_e(ib1) == myrank_e) eko_d(ib1) = eko_l(map_z(ib1),ik)  ! MPI
    end do
    call mpi_allreduce(eko_d,eko_d_mpi,neg,mpi_double_precision,mpi_sum &
         & ,mpi_k_world(myrank_k),ierr)       ! MPI
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

    Do is=1, ndim_spinor
      k1 = ik +is -1

      call m_ES_W_transpose_r(.false., ista_k, iend_k, k1, zaj_l, &
	&                    zat_l_noncl(1,1,1,idavid,is) )
    End do

    do ib1 = ista_e, iend_e, istep_e     ! MPI
       Do is=1, ndim_spinor
	  k1 = ik +is- 1
         call m_ES_WF_in_Rspace( k1, ib1, bfft_kt(:,is) )!(swffft)
       End do
       call m_FFT_Vlocal_W_noncl( afft_kt, bfft_kt, ndim_chgpot, ndim_spinor )
                                       ! (afft, bfft)-> (bfft)
       Do is=1, ndim_spinor
         call m_FFT_WF( ELECTRON, nfout, bfft_kt(:,is), DIRECT, ON )
       End do

       ib = map_z(ib1)                 ! MPI
       if ( kimg == 1 ) then
          Do is=1, ndim_spinor
            k1 = ik + is -1
            do ii=1,iba(ik)
               i1  = igf(nbase(ii,ik))
               dr1 = zaj_l( ii, ib, k1, 1 )
               dr2 = bfft_kt( i1,is )*denom
               zaj_l( ii, ib, k1, 1 )= ekin(ii)*dr1 +dr2
            end do
#ifdef SAVE_FFT_TIMES
            if(sw_save_fft == ON) status_saved_phifftr(ib,k1) = OLD
#endif
          Enddo
       else
          Do is=1, ndim_spinor
            k1 = ik + is -1
            do ii=1,iba(ik)
               i1  = igf(nbase(ii,ik))
               dr1  = zaj_l( ii, ib, k1,   1 )
               di1  = zaj_l( ii, ib, k1,kimg )
               zaj_l( ii, ib, k1,    1 ) = ekin(ii)*dr1 &
	&                                + bfft_kt( 2*i1-1,is )*denom
               zaj_l( ii, ib, k1, kimg ) = ekin(ii)*di1 &
	&                                + bfft_kt( 2*i1,  is )*denom
            enddo
#ifdef SAVE_FFT_TIMES
            if(sw_save_fft == ON) status_saved_phifftr(ib,k1) = OLD
#endif
          End do
       endif
    enddo
    if(ipridavidson>=2 .and. ik==1) write(nfout,'(" !### zaj_l is new,  bfft is old")')
!( tenchi ) (zah_l <- zaj_l)
    Do is=1, ndim_spinor
      k1 = ik +is -1
      call m_ES_W_transpose_r(.false., ista_k, iend_k, k1, zaj_l, zah_l_noncl(1,1,1,is) )
    End do

! (make matrix elements )
! parallel loop
    ! <n|T+Vloc|m> G-wise parallel
    do ibb2 = nsize_max_sb_now-neg+1,nsize_max_sb_now
       if(ibover(ibb2)<0) cycle
       ib2 = ibover(ibb2)
       iter2 = (ibb2-1)/neg+1
       ii2  = ibb2-neg*(iter2-1)
       ip0b = ib2*(ib2-1)/2

       do ibb1 = 1,ibb2
          if(ibover(ibb1)<0) cycle
          ib1 = ibover(ibb1)
          iter1 = (ibb1-1)/neg+1
          ii1 = ibb1-neg*(iter1-1)
          ip0 = ip0b + ib1

          if ( kimg == 1 ) then
             w1hw2(ip0) = 0.d0
             w1sw2(ip0) = 0.d0

             Do is=1, ndim_spinor
               do ii = 1, np_g1k(ik)            ! MPI
                  hr2 = zah_l_noncl( ii, ii2, 1, is )
                  dr2 = zat_l_noncl( ii, ii2, 1, iter2, is )
                  dr1 = zat_l_noncl( ii, ii1, 1, iter1, is )
                  w1hw2(ip0) = w1hw2(ip0) + dr1*hr2
                  w1sw2(ip0) = w1sw2(ip0) + dr1*dr2
               end do
             End do
          else
             if ( k_symmetry(ik) == GAMMA ) then
                w1hw2(ip0) = 0.d0
                w1sw2(ip0) = 0.d0
!!$                ig1 = 1;  if(mype == 0) ig1 = 2
                ig1 = 1;  if(myrank_e == 0) ig1 = 2
!
                Do is=1, ndim_spinor
                  do ii = ig1, np_g1k(ik)
                     hr2 = zah_l_noncl( ii, ii2, 1, is )
                     hi2 = zah_l_noncl( ii, ii2, 2, is )
                     dr2 = zat_l_noncl( ii, ii2, 1, iter2, is )
                     di2 = zat_l_noncl( ii, ii2, 2, iter2, is )
                     dr1 = zat_l_noncl( ii, ii1, 1, iter1, is )
                     di1 = zat_l_noncl( ii, ii1, 2, iter1, is )
                     w1hw2(ip0) = w1hw2(ip0) +( dr1*hr2 +di1*hi2 )*2.d0
                     w1sw2(ip0) = w1sw2(ip0) +( dr1*dr2 +di1*di2 )*2.d0
                  end do
                End do
!!$                if(mype == 0) then
                if ( myrank_e == 0 ) then
                   Do is=1, ndim_spinor
                      hr2 = zah_l_noncl( 1, ii2, 1, is )
                      hi2 = zah_l_noncl( 1, ii2, 2, is )
                      dr2 = zat_l_noncl( 1, ii2, 1, iter2, is )
                      di2 = zat_l_noncl( 1, ii2, 2, iter2, is )
                      dr1 = zat_l_noncl( 1, ii1, 1, iter1, is )
                      di1 = zat_l_noncl( 1, ii1, 2, iter1, is )
                      w1hw2(ip0) = w1hw2(ip0) +dr1*hr2 +di1*hi2
                      w1sw2(ip0) = w1sw2(ip0) +dr1*dr2 +di1*di2
                   End do
                end if
             else
                w1hw2( 2*ip0-1:2*ip0 ) = 0.d0
                w1sw2( 2*ip0-1:2*ip0 ) = 0.d0
                Do is=1, ndim_spinor
                   do ii = 1, np_g1k(ik)
                      hr2 = zah_l_noncl( ii, ii2, 1, is )
                      hi2 = zah_l_noncl( ii, ii2, 2, is )
                      dr2 = zat_l_noncl( ii, ii2, 1, iter2, is )
                      di2 = zat_l_noncl( ii, ii2, 2, iter2, is )
                      dr1 = zat_l_noncl( ii, ii1, 1, iter1, is )
                      di1 = zat_l_noncl( ii, ii1, 2, iter1, is )
                      w1hw2( 2*ip0-1 ) = w1hw2( 2*ip0-1 ) +dr1*hr2 +di1*hi2
                      w1hw2( 2*ip0   ) = w1hw2( 2*ip0   ) +dr1*hi2 -di1*hr2
                      w1sw2( 2*ip0-1 ) = w1sw2( 2*ip0-1 ) +dr1*dr2 +di1*di2
                      w1sw2( 2*ip0   ) = w1sw2( 2*ip0   ) +dr1*di2 -di1*dr2
                  end do
                End do
             end if
          end if
       end do
    end do
! --
    if ( ipridavidson >= 3 ) call wd_w1hw2(" -- w1hw2 without nl part--")
    ! <n|Vnl|m> G-wise parallel
    call add_nonlocal_part_noncl      ! w1hw2 = w1hw2 + w1Vnlw2
                                      ! w1sw2 = w1sw2 + w1qw2
    if(ipridavidson >= 3) call wd_w1hw2(" -- w1hw2 with nl part--")

!! (spread sum of w1hw2 and w1sw2)

    if (npes > 1) then
       w1hw2_mpi = 0.d0
       w1sw2_mpi = 0.d0
       nshift = nsize_mt_old*kimg_t
       ndata = (nsize_mt_now-nsize_mt_old)*kimg_t
       call mpi_allreduce(w1hw2(nshift+1),w1hw2_mpi,ndata,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr) ! MPI

       w1hw2(nshift+1:nshift+ndata) = w1hw2_mpi(1:ndata) ! MPI
       call mpi_allreduce(w1sw2(nshift+1),w1sw2_mpi,ndata,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr) ! MPI
       w1sw2(nshift+1:nshift+ndata) = w1sw2_mpi(1:ndata) ! MPI
    end if

    if (ipridavidson >= 2) call wd_w1hw2(" -- just after making w1hw2 --")
    if (ipridavidson >= 2) then
       write(nfout,*) 'neordr for ik = ',ik
       write(nfout,9002) (neordr(ib1,ik),ib1=1,neg)
       write(nfout,*) 'nrvf_ordr for ik = ',ik
       write(nfout,9002) (nrvf_ordr(ib1,ik),ib1=1,neg)
       write(nfout,*) 'eig'
       write(nfout,'(5x,10f8.4)') (eko_l(ib1,ik),ib1=1,np_e)
    endif
9002 format(5x,10i8)

!! (Diagonalization )  !!

#ifdef _USE_SCALAPACK_
    if(sw_dav_scalapack == ON) then
       call initialize_scalapack()
       if(kimg_t == 1) then
          call pdsygvx_driver(eig,vec,w1hw2,w1sw2,ierr_diag)
       else
          call pzhegvx_driver(eig,vec,w1hw2,w1sw2,ierr_diag)
       endif
! === DEBUG by tkato 2012/12/20 ================================================
!      call finalize_scalapack()
! ==============================================================================
    else
       if(kimg_t == 1) then
          call dspgvx_driver(eig,vec,w1hw2,w1sw2,ierr_diag)
       else
          call zhpgvx_driver(eig,vec,w1hw2,w1sw2,ierr_diag)
       endif
    endif
#else
    if(kimg_t == 1) then
       call dspgvx_driver(eig,vec,w1hw2,w1sw2,ierr_diag)
    else
       call zhpgvx_driver(eig,vec,w1hw2,w1sw2,ierr_diag)
    endif
#endif

    frestart = .false.
    if ( ierr_diag /= 0 ) then
       if(ipridavidson>=2 .and. ik==1) write(nfout,'(" !### zaj_l is new,  bfft is old")')
       Do is=1, ndim_spinor
         zaj_l(:,:,ik+is-1,:) = zaj_l_backup_noncl(:,:,:,is)   ! ????
       End do

       do ib1 = 1, neg
          if(map_e(ib1) == myrank_e) then         ! MPI
             eko_l(map_z(ib1),ik)=eko_d(ib1)
          end if
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

!!$       if(ipridavidson >= 2) then
       if(ipri0 >= 2) then
          write(nfout,*) 'eko_d for ik = ',ik
          write(nfout,9001) (eko_d(ib),ib=1,neg)
          write(nfout,*) 'eig for ik = ',ik
          write(nfout,9001) (eig(ib),ib=1,neg)
          call wd_w1hw2(" -- after diagonalization --")
!sum eko
          dr1=0.d0;dr2=0.d0
          do ib1=1,neg
             ib1to = neordr(ib1,ik)
             if(map_e(ib1to) == myrank_e) dr1=dr1+eko_l(map_z(ib1to),ik) ! MPI
             dr2=dr2+eig(ib1)
          enddo
          call mpi_allreduce(dr1,di1,1,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr) ! MPI
          dr1 = di1  ! MPI
          write(nfout,'(" sum of eko_l, eig, abs diff =",3e25.10)') dr1,dr2,abs(dr2-dr1)
       endif

!! (subspace rotation) !!
       call subspace_rotation_noncl        ! vec,zat_l -> zat_l

!( tenchi ) (zaj_l <- zat_l)
       iter = min(idavid+1,ndavid)
       Do is=1, ndim_spinor
         call m_ES_W_transpose_back_r(.false., ista_k, iend_k, ik+is-1, &
	&                            zaj_l, zat_l_noncl(1,1,1,iter,is) )
         zaj_l_backup_noncl(:,:,:,is) = zaj_l(:,:,ik+is-1,:)
       End do

!! (eko_l)
       do ib1 = 1, neg
          if(map_e(ib1) == myrank_e) then         ! MPI
             eko_l(map_z(ib1),ik)=eig(ib1)
          end if
       end do
       if(ipridavidson >= 2) then
          eko2 = sum(eig(1:neg))
          write(nfout,1201) ik,eko1,ekod,eko2

          write(nfout,*) 'eko_l'
          write(nfout,9001) (eko_l(ib1,ik),ib1=1,neg)
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

    subroutine add_nonlocal_part_noncl
      integer :: ip,ib1,ib2,ibb1,ibb2
      integer       :: ia, lmt1, lmt2, it, p, s, ib

      integer :: is1, is2, is_tmp
      integer :: il1, il2, im1, im2
      integer :: mdvdb

      complex(kind=CMPLDP) :: facv,facq
      real(kind=DP) :: vr,vi,qr,qi, cv1, cv2, cq1, cq2

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

! ---------------------------------------------------- 11.0S
#ifdef SKIP_TEST
               mdvdb = m_PP_include_vanderbilt_pot(it)
#endif
! ---------------------------------------------------- 11.0S

               do lmt1 = 1, ilmt(it)
                  p = lmta(lmt1,ia)
                  il1 = ltp(lmt1,it); im1 = mtp(lmt1,it)

                  do lmt2 = 1, ilmt(it)
                     s = lmta(lmt2,ia)
                     il2 = ltp(lmt2,it); im2 = mtp(lmt2,it)

! ---------------------------------------------------- 11.0S
#ifdef SKIP_TEST
                     if ( mdvdb == SKIP ) then
                        if ( il1 /= il2 ) cycle
                        if ( SpinOrbit_mode == Neglected .and. sw_hubbard == OFF ) then
                           if ( im1 /= im2 ) cycle
                        endif
                     endif
#endif
! ---------------------------------------------------- 11.0S

                     Do is1=1, ndim_spinor
                        Do is2=1, ndim_spinor
                           is_tmp = ndim_spinor*( is1 -1 ) +is2
                           facv = iwei(ia) *dion_scr_noncl( lmt1, lmt2, is_tmp, ia )
                           facq = iwei(ia) *q_noncl( lmt1, lmt2, is_tmp, it )

                           if ( kimg==1 ) then
                              tmpr = fsr_noncl( ii1, p, iter1, is1 ) &
	&                           *fsr_noncl( ii2, s, iter2, is2 ) &
        &                          + fsi_noncl( ii1, p, iter1, is1 ) &
	&                           *fsi_noncl( ii2, s, iter2, is2 )
                              tmpi = fsr_noncl( ii1, p, iter1, is1 )&
        &                           *fsi_noncl( ii2, s, iter2, is2 )&
        &                          - fsi_noncl( ii1, p, iter1, is1 )&
	&                           *fsr_noncl( ii2, s, iter2, is2 )
                              vr = vr + real( facv )*tmpr -aimag( facv )*tmpi
                              qr = qr + real( facq )*tmpr -aimag( facq )*tmpi
                           else
                              if ( k_symmetry(ik) == GAMMA ) then
                                 tmpr = fsr_noncl( ii1, p, iter1, is1 )&
	&                              *fsr_noncl( ii2, s, iter2, is2 )
                                 vr = vr + real( facv )*tmpr
                                 qr = qr + real( facq )*tmpr
                              else
                                 tmpr = fsr_noncl( ii1, p, iter1, is1 )&
	&                              *fsr_noncl( ii2, s, iter2, is2 )&
                    &                 + fsi_noncl( ii1, p, iter1, is1 )&
	&                              *fsi_noncl( ii2, s, iter2, is2 )
	                         tmpi = fsr_noncl( ii1, p, iter1, is1 )&
	&                              *fsi_noncl( ii2, s, iter2, is2 )&
        &                             - fsi_noncl( ii1, p, iter1, is1 )&
	&                              *fsr_noncl( ii2, s, iter2, is2 )

                                 cv1 = real(facv);  cv2 = aimag(facv)
                                 cq1 = real(facq);  cq2 = aimag(facq)

                                 vr = vr + cv1 *tmpr - cv2 *tmpi
                                 vi = vi + cv1 *tmpi + cv2 *tmpr
                                 qr = qr + cq1 *tmpr - cq2 *tmpi
                                 qi = qi + cq1 *tmpi + cq2 *tmpr

                              endif
                           end if
                        End do
                     End do
                  end do
               end do
            end do

            if ( kimg_t==1 ) then
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
    end subroutine add_nonlocal_part_noncl

    subroutine subspace_rotation_noncl
      integer :: ib1,ib2,ibb2
      real(kind=DP), allocatable, dimension(:,:,:,:) :: zaj_wk

      allocate(zaj_wk(np_g1k_x,neg,kimg,ndim_spinor))
      zaj_wk = 0.d0

      if ( kimg==1 ) then
         do ib1=1,neg
            do ibb2=1,nsize_max_sb_now
               if(ibover(ibb2)<0) cycle
               ib2 = ibover(ibb2)
               iter2=(ibb2-1)/neg+1
               ii2=ibb2-neg*(iter2-1)

               Do is=1, ndim_spinor
                 k1 = ik +is -1
                 do ii=1,np_g1k(ik)
                    zaj_wk( ii, ib1, kimg, is ) = zaj_wk( ii, ib1, kimg,is ) &
	&                  + zat_l_noncl( ii, ii2, kimg, iter2, is )*vec( ib2,ib1 )
                 end do
               End do
            end do
         end do
      else
         if ( k_symmetry(ik) == GAMMA ) then
            do ib1=1,neg
               do ibb2=1,nsize_max_sb_now
                  if(ibover(ibb2)<0) cycle
                  ib2 = ibover(ibb2)
                  iter2=(ibb2-1)/neg+1
                  ii2=ibb2-neg*(iter2-1)
                  hr2=vec(ib2,ib1)

                  Do is=1, ndim_spinor
                     k1 = ik +is -1
                     do ii=1,np_g1k(ik)
                        dr1 = zat_l_noncl( ii, ii2, 1   , iter2, is )
                        di1 = zat_l_noncl( ii, ii2, kimg, iter2, is )
                        zaj_wk( ii, ib1, 1, is ) &
	&                   = zaj_wk( ii, ib1, 1, is ) + dr1*hr2
                        zaj_wk( ii, ib1, kimg, is ) &
	&                   = zaj_wk( ii, ib1, kimg, is ) + di1*hr2
                     end do
                  End do
               end do
            end do
         else
            do ib1=1,neg
               do ibb2=1,nsize_max_sb_now
                  if(ibover(ibb2)<0) cycle
                  ib2 = ibover(ibb2)
                  iter2=(ibb2-1)/neg+1
                  ii2=ibb2-neg*(iter2-1)
                  hr2=vec(2*ib2-1,ib1)
                  hi2=vec(2*ib2  ,ib1)

                  Do is=1, ndim_spinor
                     k1 = ik +is -1
                     do ii=1,np_g1k(ik)
                        dr1 = zat_l_noncl( ii, ii2, 1   , iter2, is )
                        di1 = zat_l_noncl( ii, ii2, kimg, iter2, is )
                        zaj_wk( ii, ib1, 1, is ) &
	&                   = zaj_wk( ii, ib1, 1, is ) + dr1*hr2 - di1*hi2
                        zaj_wk( ii, ib1, kimg, is ) &
	&                   = zaj_wk( ii, ib1, kimg,is ) + dr1*hi2 + di1*hr2
                     end do
                  End do
               end do
            end do
         end if
      end if
      iter = min(idavid+1,ndavid)

      Do is=1, ndim_spinor
        zat_l_noncl(:,:,:,iter,is) = zaj_wk(:,:,:,is)
      End do

      deallocate(zaj_wk)
    end subroutine subspace_rotation_noncl

  end subroutine evolve_WFs_in_subspace_noncl
!======================================================================= 11.0


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

#ifndef PARA3D
  subroutine Vnonlocal_Diagonal_part(ispin,ik,iksnl,vnldi)
    integer, intent(in)                        :: ispin, ik, iksnl
    real(kind=DP), intent(out), dimension(kg1) :: vnldi

    integer :: it,mdvdb

    vnldi = 0.d0

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
      integer       :: ia, p1,p2,lmtt1,il1,lmtt2,il2,i
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
               do i = 1, iba(ik)
                  vnldi(i) = vnldi(i)+fac*snl(i,lmtt1,iksnl)*snl(i,lmtt2,iksnl)
               end do
            end do
         end do
      end do
    end subroutine Vnonlocal_D_vanderbilt_case

    subroutine Vnonlocal_D_norm_conserve_case
      integer       :: ia, lmt1,lmt2,lmtt1,il1,im1,il2,im2,i
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

! ======================================== modified by K. Tagami =========== 11.0
!            if(ipaw(it)==0)then
!               fac = ph * dion(lmt1,lmt2,it)
!            else
!               fac = ph * dion_paw(lmt1,lmt2,ispin,ia)
!            endif
!
            fac = 0.d0
            Do ia=1, natm
              if(ityp(ia) /= it) cycle
              if(ipaw(it)==0)then
                 fac = fac + iwei(ia) * dion(lmt1,lmt2,it)
              else
                 fac = fac + iwei(ia) * dion_paw(lmt1,lmt2,ispin,ia)
              endif
            End do
! ========================================================================= 11.0
            do i = 1, iba(ik)
               vnldi(i)  = vnldi(i) + fac * snl(i,lmtt1,iksnl)*snl(i,lmtt1,iksnl)
            end do
         end do
      end do
    end subroutine Vnonlocal_D_norm_conserve_case

  end subroutine Vnonlocal_Diagonal_part

! ====================================== added by K. Tagami ================ 11.0
  subroutine Vnonlocal_Diagonal_part_noncl(ispin,ik,iksnl,vnldi)
    integer, intent(in)                        :: ispin, ik, iksnl

#ifdef USE_COMPLEX_VNLDI
    complex(kind=CMPLDP), intent(out), dimension(kg1) :: vnldi
#else
    real(kind=DP), intent(out), dimension(kg1) :: vnldi
#endif

    integer :: it,mdvdb,ib

!!$    integer :: id_sname = -1
!!$    call tstatc0_begin('Vnonlocal_Diagonal_part ', id_sname)

!    vnldi = 0.d0

    do it = 1, ntyp
       mdvdb = m_PP_include_vanderbilt_pot(it)
       if(mdvdb == SKIP) then
          call Vnonlocal_D_norm_consv_noncl
       else if(mdvdb == EXECUT) then
          call Vnonlocal_D_vanderbilt_noncl
       end if
    end do

!!$    call tstatc0_end(id_sname)
  contains

    subroutine Vnonlocal_D_vanderbilt_noncl
     integer       :: ia, p1,p2,lmtt1,il1,lmtt2,il2,i
      integer :: im1, im2
      complex(kind=CMPLDP) :: ph, fac

      do p1 = 1, ilmt(it)
         lmtt1 = lmtt(p1,it); il1 = ltp(p1,it);  im1 = mtp(p1,it)

         do p2 = 1, ilmt(it)
            lmtt2 = lmtt(p2,it); il2 = ltp(p2,it);  im2 = mtp(p2,it)

            if( p1 /= p2 ) then
               ph = zi**(il2-il1)
            else
               ph = 1.d0
            endif

!            if(mod(il1+il2,2) == 1) cycle

            fac = 0.d0

            do ia = 1, natm
               if(ityp(ia) /= it) cycle
               fac = fac + ph*iwei(ia)* dion_scr_noncl( p1,p2,ispin,ia )
            end do

            do i = 1, iba(ik)
               vnldi(i) = vnldi(i)+fac*snl(i,lmtt1,iksnl)*snl(i,lmtt2,iksnl)
            end do
         end do
      end do

    end subroutine Vnonlocal_D_vanderbilt_noncl

    subroutine Vnonlocal_D_norm_consv_noncl
      integer       :: ia, p1,p2,lmtt1,il1,lmtt2,il2,i
      integer :: im1, im2
      complex(kind=CMPLDP) :: ph, fac

      do p1 = 1, ilmt(it)
         lmtt1 = lmtt(p1,it); il1 = ltp(p1,it);  im1 = mtp(p1,it)

         do p2 = 1, ilmt(it)
            lmtt2 = lmtt(p2,it); il2 = ltp(p2,it);  im2 = mtp(p2,it)

            if ( il1 /= il2 ) cycle
            if ( SpinOrbit_mode == Neglected .and. sw_hubbard == OFF ) then
               if ( im1 /= im2 ) cycle
            endif

            if( p1 /= p2 ) then
               ph = zi**(il2-il1)
            else
               ph = 1.d0
            endif

!            if(mod(il1+il2,2) == 1) cycle

            fac = 0.d0

            do ia = 1, natm
               if(ityp(ia) /= it) cycle
               fac = fac + ph*iwei(ia)* dion_scr_noncl( p1,p2,ispin,ia )

            end do

            do i = 1, iba(ik)
               vnldi(i) = vnldi(i)+fac*snl(i,lmtt1,iksnl)*snl(i,lmtt2,iksnl)
            end do
         end do
      end do

    end subroutine Vnonlocal_D_norm_consv_noncl

  end subroutine Vnonlocal_Diagonal_part_noncl
! ==================================================================== 11.0



! ====================================== modified by K. Tagami ============ 11.0
!!!  subroutine S_Diagonal_part(ik,iksnl,sdiag)
!!    integer, intent(in)                        :: ik, iksnl
!
  subroutine S_Diagonal_part(ispin,ik,iksnl,sdiag)
    integer, intent(in)                        :: ik, iksnl, ispin
! ========================================================================= 11.0

    real(kind=DP), intent(out), dimension(kg1) :: sdiag

    integer :: it,mdvdb

    sdiag = 1.d0
    do it = 1, ntyp
       mdvdb = m_PP_include_vanderbilt_pot(it)
       if(mdvdb /= SKIP) call Vanderbilt_case
    end do

  contains

    subroutine Vanderbilt_case
      integer       :: ia, p1,p2,lmtt1,il1,lmtt2,il2,i
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
               do i = 1, iba(ik)
                  sdiag(i) = sdiag(i)+fac*snl(i,lmtt1,iksnl)*snl(i,lmtt2,iksnl)
               end do
            end do
         end do
      end do
    end subroutine Vanderbilt_case

  end subroutine S_Diagonal_part

! ================================== added by K. Tagami ====================== 11.0
  subroutine S_Diagonal_part_noncl(ispin,ik,iksnl,sdiag)
    integer, intent(in)                        :: ik, iksnl, ispin

#ifdef USE_COMPLEX_VNLDI
    complex(kind=CMPLDP), intent(out), dimension(kg1) :: sdiag
#else
    real(kind=DP), intent(out), dimension(kg1) :: sdiag
#endif

    integer :: it,mdvdb

!    sdiag = 1.d0

    do it = 1, ntyp
       mdvdb = m_PP_include_vanderbilt_pot(it)
       if (mdvdb /= SKIP) then
          Call Vanderbilt_case_noncl2
       endif
    end do

  contains

    subroutine Vanderbilt_case_noncl2
      integer       :: ia, p1,p2,lmtt1,il1,lmtt2,il2,i
      integer :: im1, im2
      complex(kind=CMPLDP) :: ph,fac


      do p1 = 1, ilmt(it)
         lmtt1 = lmtt(p1,it); il1 = ltp(p1,it);  im1 = mtp(p1,it)

         do p2 = 1, ilmt(it)
            lmtt2 = lmtt(p2,it); il2 = ltp(p2,it);  im2 = mtp(p2,it)

            if ( il1 /= il2 ) cycle
            if ( SpinOrbit_mode /= BuiltIn ) then
               if ( im1 /= im2 ) cycle
            endif

            if( p1 /= p2) then
               ph = zi**(il2-il1)
            else
               ph = 1.d0
            endif

!            if(mod(il1+il2,2) == 1) cycle

            fac = 0.0d0
            do ia = 1, natm
               if(ityp(ia) /= it) cycle
               fac = fac + ph*iwei(ia) * q_noncl(p1,p2,ispin,it)
            End do

            do i = 1, iba(ik)
               sdiag(i) = sdiag(i)+ fac *snl(i,lmtt1,iksnl)*snl(i,lmtt2,iksnl)
            end do
         end do
      end do
    end subroutine Vanderbilt_case_noncl2

  end subroutine S_Diagonal_part_noncl
! ============================================================================ 11.0


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
    call set_nprow_npcol(nprow, npcol)

    write(nfout, '("set nprow, npcol = ", 4i5)') nprow, npcol

    allocate(usermap(nprow,npcol))
    do j = 1, npcol
       do i = 1, nprow
          usermap(i,j) = myrank_k*nrank_e + (i - 1)*npcol + j - 1
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

  subroutine set_nprow_npcol(npr, npc)
    integer, intent(inout) :: npr, npc
    logical :: cflag
    integer :: j, n

    cflag = .true.

    if(npr*npc /= nrank_e) then
       npr = 1
       npc = 1
       j = nrank_e
 120   continue
       if(mod(j, 2) == 0) then
          j = j/2
          n = 2
          call prod(cflag, n, npr, npc)
          goto 120
       else if(mod(j, 3) == 0) then
         j = j/3
         n = 3
         call prod(cflag, n, npr, npc)
         goto 120
       else if(mod(j, 5) == 0) then
         j = j/5
         n = 5
         call prod(cflag, n, npr, npc)
         goto 120
       else if(j == 1) then
         n = 1
       else
         n = j
       end if
       call prod(cflag, n, npr, npc)
    end if
  end subroutine set_nprow_npcol

  subroutine prod(cflag, n, nr, nc)
    logical, intent(inout) :: cflag
    integer, intent(in)    :: n
    integer, intent(inout) :: nr,nc
    if(cflag) then
       nr = nr * n
    else
       nc = nc * n
    endif
    cflag = .not. cflag
  end subroutine prod
#endif
#endif
end module m_ES_WF_by_Davidson

