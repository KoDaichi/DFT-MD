!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  MODULE: m_epc_potential
!
!  AUTHOR(S): M. Saito, T. Yamasaki   November/04/2003
!
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)
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
!
!
module m_epc_potential
! $Id: m_epc_potential.F90 570 2017-04-21 20:34:50Z yamasaki $
  use m_Electronic_Structure, only : vlhxc_l
!$$#ifndef PARA3D
  use m_PlaneWaveBasisSet,    only : ngabc,gr_l,igfp_l,igfp_nonpara,kg,kgp,kgp_reduced,ylm_l&
       &                           , m_pwBS_sphrp2,m_pwBS_sphrp2_diff
!$$#endif
#ifdef _MPIFFTTEST_
  use m_PlaneWaveBasisSet,    only : igfp_l_c
#endif
  use m_PseudoPotential,      only : itpcc,rhcg_l,rhceg_l

  use m_Crystal_Structure,    only : rltv,univol
  use m_Ionic_System,         only : ntyp,natm,iwei,ityp,pos,zfm3_l,iatomn
  use m_FFT,                  only : fft_box_size_CD,fft_box_size_CD_nonpara,nfftp,nfftp_nonpara &
       &                           , m_FFT_CD0 &
       &                           , m_FFT_check_of_negative_CD
#ifdef _MPIFFTTEST_
  use m_FFT,                  only : m_FFT_set_cdata  &
       &                           , lx,ly,lz,ly_d,lz_d, ny_d,nz_d
#else
  use m_FFT,                  only : m_FFT_alloc_CD_box &
       &                           , m_FFT_dealloc_CD_box
#endif
  use m_Timing,               only : tstatc0_begin, tstatc0_end
  use m_Control_Parameters,   only : xctype,len_xctype,ipripositron,nspin,kimg,af,nel_Ylm &
       &                           , sw_gga_p
  use m_Const_Parameters,     only : DP,DOWN,UP,Partial_Core_Charge,Valence_plus_PC_Charge,PAI4,ON &
       &                           , LDA,GGA,PAI,PAI2,VXC_AND_EXC,STRESS_,zi,SKIP, INVERSE, DIRECT
  use m_Parallelization,      only : MPI_CommGroup,ista_kngp,iend_kngp,npes,mype, ierr &
       &                           , ista_sfftp,iend_sfftp,ista_sfftph,iend_sfftph &
       &                           , nis_sfftp,nie_sfftp,nel_sfftp,idisp_sfftp,np_sfftp,mp_sfftp &
       &                           , npes_cdfft

! === Postitron SCF === 2015/11/28
  use m_Control_Parameters,  only : positron_method
  use m_Const_Parameters,   only : positron_CONV
  use m_Charge_Density,     only : chgq_l
  use m_PseudoPotential,      only : psc_l
! ===================== 2015/11/28
  use mpi

  implicit none

!  include 'mpif.h'

  real(kind=DP), pointer :: vlhxc_p( :,:,: )

  real(kind=DP), allocatable, dimension(:,:,:) :: tchgq_l  ! d(ista_kngp:iend_kngp,kimg)
  real(kind=DP), allocatable, dimension(:,:) ::   tchgr_l  ! d(ista_sfftph:iend_sfftph,nspin)
!!$  real(kind=DP), allocatable, dimension(:) ::     p_potential_l ! d(ista_kngp:iend_kngp)
  real(kind=DP), allocatable, dimension(:,:) ::   grad_tchgr_l !d(ista_sfftph:iend_sfftph,nspin)
  real(kind=DP), allocatable, dimension(:,:) ::   chden_l !d(ista_sfftp:iend_sfftp,nspin)
  integer, private,allocatable,dimension(:)  ::   inx,jnx,knx ! MPI d(ista_sfftp:iend_sfftp)
  real(kind=DP), allocatable, dimension(:,:,:) :: vepc_l   ! d(ista_kngp:iend_kngp,kimg,nspin)  MPI

  real(kind=DP)                                :: epc

  real(kind=DP),private,allocatable,dimension(:):: afft        ! d(ista_sfftp:iend_sfftp)
  real(kind=DP),private,pointer,dimension(:)    :: afft_mpi1,afft_CD   ! MPI d(nfftp)
  real(kind=DP),private,pointer,dimension(:)    :: afft_mpi2,afft_mpi3   ! MPI d(mp_sfftp)

  real(kind=DP),private, allocatable, dimension(:)  :: f2or1 ! MPI d(ista_sfftph,iend_sfftph)
  logical,      private, dimension(3)               :: lmn_even

!   1. m_epc_alloc
!   2. m_epc_dealloc
!   3. epc_allocate
!   4. epc_deallocate
!   5. m_epc_cal_potential
!       map_charge_onto_a_fft_box
!       cp_afft_to_tchgrhr
!       check_lmn_even
!       cpafft
!       cpafft_to_vepc
!   6. check_of_pcc
!   7. afft_allgatherv


contains
  subroutine m_epc_alloc()
    integer :: idp, i, n, nn, nlp, nmp, nnp, nd2p, j, k
#ifdef _MPIFFTTEST_
    integer                  :: np0, j0, i2, j2, k2, n0
#endif

    if ( allocated(tchgq_l) ) deallocate( tchgq_l )
    if ( allocated(tchgr_l) ) deallocate( tchgr_l )
    if ( allocated( grad_tchgr_l ) ) deallocate( grad_tchgr_l )
    if ( allocated( chden_l ) ) deallocate( chden_l )
    if ( allocated( inx ) ) deallocate( inx )
    if ( allocated( jnx ) ) deallocate( jnx )
    if ( allocated( knx ) ) deallocate( knx )
    if ( allocated( vepc_l ) ) deallocate( vepc_l )
    if ( allocated( f2or1) ) deallocate( f2or1 )

    allocate(tchgq_l(ista_kngp:iend_kngp,kimg,nspin));  tchgq_l = 0.d0
    allocate(tchgr_l(ista_sfftph:iend_sfftph,nspin));     tchgr_l = 0.d0
    if(sw_gga_p == ON) then
       allocate(grad_tchgr_l(ista_sfftph:iend_sfftph,nspin));       grad_tchgr_l = 0.d0
       allocate(chden_l(ista_sfftp:iend_sfftp,nspin)); chden_l = 0.d0
       allocate(inx(ista_sfftp:iend_sfftp)); inx = 0
       allocate(jnx(ista_sfftp:iend_sfftp)); jnx = 0
       allocate(knx(ista_sfftp:iend_sfftp)); knx = 0
       nlp  = fft_box_size_CD(1,1)
       nmp  = fft_box_size_CD(2,1)
       nnp  = fft_box_size_CD(3,1)
#ifdef _MPIFFT_
       idp  = fft_box_size_CD_nonpara(1,0)
       nd2p = fft_box_size_CD_nonpara(2,0)
#elif _MPIFFTTEST_
       idp  = fft_box_size_CD_c(1,0)
       nd2p = fft_box_size_CD_c(2,0)
#else
       idp  = fft_box_size_CD(1,0)
       nd2p = fft_box_size_CD(2,0)
#endif
       do n = ista_sfftp, iend_sfftp
          nn = (n+kimg-1)/kimg
#ifdef _MPIFFTTEST_
          np0 = nn - idp*ly_d*lz*myrank_cdfft
          i  = mod(np0-1,idp)+1
          j0 = mod((np0-i)/idp,ly_d) + 1
          j  = j0 + ny_d*myrank_cdfft
          k  = (np0-i-(j0-1)*idp)/(idp*ly_d) + 1
          if(i >= idp) i = i - idp
#else
          i  = mod(nn,idp)
          j  = mod((nn-1)/idp,nd2p) + 1
          k  = (nn - (j-1)*idp - i)/(idp*nd2p) + 1
#endif
          inx(n) = i - 1;  if(2*inx(n) > nlp) inx(n) = inx(n) - nlp
          jnx(n) = j - 1;  if(2*jnx(n) > nmp) jnx(n) = jnx(n) - nmp
          knx(n) = k - 1;  if(2*knx(n) > nnp) knx(n) = knx(n) - nnp
       end do
    end if
!!$    allocate(p_potential_l(ista_kngp:iend_kngp)); p_potential_l = 0.d0
    allocate(vepc_l(ista_kngp:iend_kngp,kimg,nspin)); vepc_l = 0.d0

    allocate(f2or1(ista_sfftph:iend_sfftph))
    if(kimg == 1) then
       f2or1 = 2.d0
#ifdef _MPIFFT_
       idp = fft_box_size_CD_nonpara(1,0)
#else
       idp = fft_box_size_CD(1,0)
#endif
       do i = ista_sfftph, iend_sfftph
          if(mod(i*2,idp) == 2 .or. mod(i*2,idp) == 0) f2or1(i) = 1.d0
       end do
    else
       f2or1 = 1.d0
    end if

  end subroutine m_epc_alloc

  subroutine m_epc_dealloc()
!!$    if(allocated(tchgq_l)) deallocate(tchgq_l)
    if(allocated(tchgr_l)) deallocate(tchgr_l)
    if(allocated(grad_tchgr_l)) deallocate(grad_tchgr_l)
    if(allocated(chden_l)) deallocate(chden_l)
    if(allocated(inx)) deallocate(inx)
    if(allocated(jnx)) deallocate(jnx)
    if(allocated(knx)) deallocate(knx)
!!$    if(allocated(p_potential_l)) deallocate(p_potential_l)
    if(allocated(f2or1)) deallocate(f2or1)
!Saito
!    if(allocated(vepc_l)) deallocate(vepc_l)
  end subroutine m_epc_dealloc

  subroutine epc_allocate()
#ifndef _MPIFFTTEST_
    call m_FFT_alloc_CD_box()
#endif
    allocate(afft(ista_sfftp:iend_sfftp))
    allocate(afft_mpi1(nfftp_nonpara))                         ! MPI
    allocate(afft_CD(nfftp_nonpara))
    if(npes >= 2) then
       allocate(afft_mpi2(mp_sfftp))                       ! MPI
       allocate(afft_mpi3(mp_sfftp))                       ! MPI
    end if

  end subroutine epc_allocate

  subroutine epc_deallocate()
    call m_FFT_dealloc_CD_box()
    deallocate(afft)
    deallocate(afft_mpi1)
    deallocate(afft_CD)
    if(npes >= 2) then
       deallocate(afft_mpi2); deallocate(afft_mpi3)
    end if
  end subroutine epc_deallocate

  subroutine m_epc_cal_potential(nfout,chgq_l)
    integer, intent(in) :: nfout
    real(kind=DP), intent(in) :: chgq_l(ista_kngp:iend_kngp,kimg,nspin)

    integer :: i, n, j, ispin, iloop, in
    real(kind=DP) :: rinplw
    real(kind=DP), allocatable, dimension(:,:) :: vepcr_l

    real(kind=DP) :: epc = 0.d0
    real(kind=DP) :: epc_mpi     ! MPI
    real(kind=DP) :: DELTA
    data DELTA/1.d-40/
#ifdef _MPIFFTTEST_
    rinplw = 1.d0/product(fft_box_size_CD_c(1:3,1))
#else
    rinplw = 1.d0/product(fft_box_size_CD(1:3,1))
#endif

!!$    do j = 1, kimg
!!$       do n = 1, ntyp
!!$          do i = ista_kngp, iend_kngp
!!$             tchgq_l(i,j,1) = tchgq_l(i,j,1) + (rhcg_l(i,n)-iatomn(n))*zfm3_l(i,n,j)
!!$          end do
!!$       end do
!!$    end do
!!$
!!$    if(nspin == 2) then
!!$       tchgq_l(:,:,2) = tchgq_l(:,:,1)*0.5d0
!!$       tchgq_l(:,:,1) = tchgq_l(:,:,2)
!!$    end if

    tchgq_l = 0.0d0
    do ispin = 1, nspin
       do j = 1, kimg
          do i = ista_kngp, iend_kngp
             tchgq_l(i,j,ispin) = tchgq_l(i,j,ispin) + chgq_l(i,j,ispin)
          end do
       end do
    end do

    if(ipripositron >= 2) then
       write(nfout,'(" -- tchgq_l, chgq_l <<m_epc_cal_potential>>")')
       do ispin = 1, nspin
          write(nfout,'(" ispin = ",i8)') ispin
          do j = 1, kimg
             write(nfout,'(" j = ",i8)') j
             do i = ista_kngp, min(ista_kngp+20,iend_kngp)
                write(nfout,'(3f10.5)') tchgq_l(i,j,ispin),chgq_l(i,j,ispin)
             end do
          end do
       end do
    end if

    call epc_allocate()

    do iloop = 1, nspin
       call map_charge_onto_a_fft_box ! tchgq_l -> afft_CD (G-space -> G-space)
       if(sw_gga_p == ON) &
            & chden_l(ista_sfftp:iend_sfftp,iloop) = afft_CD(ista_sfftp:iend_sfftp) ! G-space
       call m_FFT_CD0(nfout,afft_CD,INVERSE) ! afft -> afft (G-space representation -> r-space representation)
       afft(ista_sfftp:iend_sfftp) = afft_CD(ista_sfftp:iend_sfftp)
       call m_FFT_check_of_negative_CD(npes,ista_sfftp,iend_sfftp,ista_sfftph,iend_sfftph &
            &                         ,afft,nfout,nspin,iloop)
       call cp_afft_to_tchgrhr       ! -> tchgr_l (R-space)
    end do

    call check_lmn_even

    if(sw_gga_p == ON) then
       do iloop = 1, nspin
          grad_tchgr_l(:,iloop) = 0.d0
          do in = 1, 3
             call g_xyz_chden_l(in,iloop) ! -> G_xyz * chden_l --> afft_CD
             call m_FFT_CD0(nfout,afft_CD,INVERSE)
             call add_sq_afft_CD_to_grad_tchgr(iloop)
          end do
          grad_tchgr_l(:,iloop) = dsqrt(grad_tchgr_l(:,iloop))
       end do
    end if

    ! ******************************************
    ! epcor_0 is the subroutine coded by Dr. Mineo Saito.
    allocate(vepcr_l(ista_sfftph:iend_sfftph,nspin))
    vepcr_l = tchgr_l
    call epcor_0(nspin,DELTA,vepcr_l,f2or1,epc)
    ! -> b_XC_Potential,  Electron-Positron CORrelation , 0-density limit

    epc = epc*univol*rinplw
    if(npes >= 2) then               ! MPI
       call mpi_allreduce(epc,epc_mpi,1,mpi_double_precision,mpi_sum &
            & ,MPI_CommGroup,ierr)  ! MPI
       epc = epc_mpi                 ! MPI
    end if

    if(ipripositron >= 3) then
       write(nfout,'(" epc = ",d20.8,"<<m_epc_cal_potential>>")') epc
       write(nfout,'(" -- vepcr_l -- <<m_epc_cal_potential>>")')
       do ispin = 1, nspin
          write(nfout,'(" ispin = ",i8)') ispin
          do i = ista_sfftph, min(ista_sfftph+50,iend_sfftph)
             write(nfout,'(3f12.6)') vepcr_l(i,ispin)
          end do
       end do
    end if


    do iloop = 1, nspin
       call cpafft_CD(iloop)           ! vepcr_l -> afft_CD
       call m_FFT_CD0(nfout,afft_CD,DIRECT)
       call cpafft_CD_to_vepc(iloop)   ! afft_CD -> vepc_l
    end do
    deallocate(vepcr_l)

    if(ipripositron >= 3) then
       write(nfout,'(" -- vepc_l, rhceg_l, rhcg_l <<m_epc_cal_potential>>")')
       do n = 1, ntyp
          write(nfout,'(" ntyp = ",i8)') ntyp
          do i = ista_kngp, min(ista_kngp+20,iend_kngp)
             write(nfout,'(3f10.5)') vepc_l(i,1,1),rhceg_l(i,n), rhcg_l(i,n)
          end do
       end do
    end if

    do j = 1, kimg
       do n = 1, ntyp
          do i = ista_kngp, iend_kngp
             vepc_l(i,j,1) = vepc_l(i,j,1) + rhceg_l(i,n)*zfm3_l(i,n,j)
          end do
       end do
    end do

    if(nspin == 2) then
!!$       vepc_l(:,:,1) = tchgq_l(:,:,1)+tchgq_l(:,:,2)
       vepc_l(:,:,2) = vepc_l(:,:,1)
!!$       vepc_l(:,:,2) = tchgq_l(:,:,1)*0.5d0
!!$       vepc_l(:,:,1) = tchgq_l(:,:,2)
    end if

    tchgq_l=0.d0
    do j = 1, kimg
       do n = 1, ntyp
          do i = ista_kngp, iend_kngp
             tchgq_l(i,j,1) = tchgq_l(i,j,1) + (rhcg_l(i,n)-iatomn(n)/univol)*zfm3_l(i,n,j)
          end do
       end do
    end do
    if(nspin == 2) then
       tchgq_l(:,:,2) = tchgq_l(:,:,1)*0.5d0
       tchgq_l(:,:,1) = tchgq_l(:,:,2)
    end if
    do ispin = 1, nspin
       do j = 1, kimg
          do i = ista_kngp, iend_kngp
             tchgq_l(i,j,ispin) = tchgq_l(i,j,ispin) + chgq_l(i,j,ispin)
             !             write(671,*) i,tchgq_l(i,j,ispin)*univol
          end do
       end do
    end do
    call epc_deallocate()
  contains
    subroutine g_xyz_chden_l(in,is)
      integer, intent(in) :: in, is
      integer       :: n
      real(kind=DP) :: gxyz
      afft_CD = 0.d0
      do n = ista_sfftp, iend_sfftp
         gxyz = rltv(in,1)*inx(n)+rltv(in,2)*jnx(n)+rltv(in,3)*knx(n)
         afft_CD(n) = gxyz*chden_l(n,is)
      end do
      call boundary_zero_into_afft(in)
    end subroutine g_xyz_chden_l

    subroutine boundary_zero_into_afft(in)
      integer, intent(in) :: in
      integer             :: i,j,k,nn,n,idp,nlp,nmp,nnp,nd2p, j0

      nlp = fft_box_size_CD(1,1)
      nmp = fft_box_size_CD(2,1)
      nnp = fft_box_size_CD(3,1)
#ifdef _MPIFFT_
      idp  = fft_box_size_CD_nonpara(1,0)
      nd2p = fft_box_size_CD_nonpara(2,0)
#elif _MPIFFTTEST_
      idp = fft_box_size_CD_c(1,0)
      nd2p = fft_box_size_CD_c(2,0)
#else
      idp = fft_box_size_CD(1,0)
      nd2p = fft_box_size_CD(2,0)
#endif

      if(kimg == 1) then
         if(lmn_even(in)) then
            if( in == 1) then
#ifdef  MPIFFTTEST_
               do j = 1, ly_d*lz
                  nn = nlp/2 + 1 + idp*(j-1) + idp*ly_d*lz*myrank_cdfft
                  afft(nn) = 0.d0
               end do
#else
               do j = 1, nmp      ! y
                  do k = 1, nnp   ! z
                     nn = nlp/2 + 1 + idp*(j-1) + idp*nd2p*(k-1)
                     if(nn >= ista_sfftp .and. nn <= iend_sfftp) afft(nn) = 0.d0
                  end do
               end do
#endif
            else if( in == 2) then
#ifdef _MPIFFTTEST_
               j0 = nmp/2 + 1 - ny_d*myrank_cdfft
               if(j0 >= 1 .and. j0 <= ny_d) then
                  do i = 1, idp
                     do k = 1, nnp
                        nn = i + idp*(j0-1) + idp*ly_d*(k-1) + idp*ly_d*lz*myrank_cdfft
                        afft(nn) = 0.d0
                     end do
                  end do
               end if
#else
               do i = 1, idp      ! x
                  do k = 1, nnp   ! z
                     nn = i + idp*(nmp/2) + idp*nd2p*(k-1)
                     if(nn >= ista_sfftp .and. nn <= iend_sfftp) afft(nn) = 0.d0
                  end do
               end do
#endif
            else if(in == 3) then
#ifdef _MPIFFTTEST_
               do i = 1, idp
                  do j = 1+ny_d*myrank_cdfft, min(nmp,ny_d*(myrank_cdfft+1))
                     nn = i + idp*(j-1-ny_d*myrank_cdfft) + idp*ly_d*(nnp/2) &
                          &                          + idp*ly_d*lz*myrank_cdfft
                     afft(nn) = 0.d0
                  end do
               end do
#else
               do i = 1, idp      ! x
                  do j = 1, nmp   ! y
                     nn = i + idp*(j-1) + idp*nd2p*(nnp/2)
                     if(nn >= ista_sfftp .and. nn <= iend_sfftp) afft(nn) = 0.d0
                  end do
               end do
#endif
            end if
         end if
      else if(kimg == 2) then
         if(lmn_even(in)) then
            if( in == 1) then
#ifdef _MPIFFTTEST_
               do j = 1, ly_d*lz
                  nn = nlp/2 + 1 + idp*(j-1) + idp*ly_d*lz*myrank_cdfft
                  n = nn*2 - 1
                  afft(n) = 0.d0
                  afft(n+1) = 0.d0
               end do
#else
               do j = 1, nmp
                  do k = 1, nnp
                     nn = nlp/2 + 1 + idp*(j-1) + idp*nd2p*(k-1)
                     n  = nn*2 - 1
                     if(n >= ista_sfftp .and. n <= iend_sfftp) afft(n) = 0.d0
                     n  = nn*2
                     if(n >= ista_sfftp .and. n <= iend_sfftp) afft(n) = 0.d0
                  end do
               end do
#endif
            else if( in == 2) then
#ifdef _MPIFFTTEST_
               j0 = nmp/2 + 1 - ny_d*myrank_cdfft
               if(j0 >= 1 .and. j0 <= ny_d) then
                  do i = 1, idp
                     do k = 1, nnp
                        nn = i + idp*(j0-1) + idp*ly_d*(k-1) + idp*ly_d*lz*myrank_cdfft
                        n  = nn*2 - 1
                        afft(n) = 0.d0; afft(n+1) = 0.d0
                     end do
                  end do
               end if
#else
               do i = 1, idp
                  do k = 1, nnp
                     nn = i + idp*(nmp/2) + idp*nd2p*(k-1)
                     n  = nn*2 - 1
                     if(n >= ista_sfftp .and. n <= iend_sfftp) afft(n) = 0.d0
                     n  = nn*2
                     if(n >= ista_sfftp .and. n <= iend_sfftp) afft(n) = 0.d0
                  end do
               end do
#endif
            else if(in == 3) then
#ifdef _MPIFFTTEST_
               do i = 1, idp
                  do j = 1+ny_d*myrank_cdfft, min(nmp,ny_d*(myrank_cdfft+1))
                     nn = i + idp*(j-1-ny_d*myrank_cdfft) + idp*ly_d*(nnp/2) &
                          &                          + idp*ly_d*lz*myrank_cdfft
                     n  = nn*2 - 1
                     afft(n) = 0.d0
                     afft(n+1) = 0.d0
                  end do
               end do
#else
               do i = 1, idp
                  do j = 1, nmp
!!$                  nn = i + idp*(j-1) + idp*nmp*(nnp/2)
                     nn = i + idp*(j-1) + idp*nd2p*(nnp/2)
                     n  = nn*2 - 1
                     if(n >= ista_sfftp .and. n <= iend_sfftp) afft(n) = 0.d0
                     n  = nn*2
                     if(n >= ista_sfftp .and. n <= iend_sfftp) afft(n) = 0.d0
                  end do
               end do
#endif
            end if
         end if
      end if
    end subroutine boundary_zero_into_afft

    subroutine add_sq_afft_CD_to_grad_tchgr(is)
      integer, intent(in) :: is
      integer :: i
      do i = ista_sfftph, iend_sfftph     ! MPI
         grad_tchgr_l(i,is) = grad_tchgr_l(i,is) + afft_CD(i*2)**2
      end do
    end subroutine add_sq_afft_CD_to_grad_tchgr

    subroutine map_charge_onto_a_fft_box
      real(kind=DP)      :: fac
      integer            :: mm, it, i, j, ip

      fac = 1.d0/nspin

      mm = 0
      afft_CD = 0.d0

      do j = 1, kimg
         do i = ista_kngp, min(kgp_reduced,iend_kngp)  !for mpi
#ifdef _MPIFFT_
            ip = (igfp_nonpara(i)-1)*kimg + j
#elif _MPIFFTTEST_
            ip = (igfp_l_c(i)-1)*kimg + j
#else
            ip = (igfp_l(i)-1)*kimg + j
#endif
            afft_CD(ip) = afft_CD(ip) + tchgq_l(i,j,iloop) !mpi
         end do
      end do

      if(npes >= 2) then
         call mpi_allreduce(afft_CD,afft_mpi1,nfftp_nonpara ,mpi_double_precision&
              &,mpi_sum,MPI_CommGroup,ierr)
         afft_CD = afft_mpi1
!!$         call mpi_barrier(MPI_CommGroup,ierr)
!!$         do j = 0, npes-1
!!$            do i = nis_sfftp(j),nie_sfftp(j)
!!$               afft_mpi2(i-nis_sfftp(j)+1) = afft_mpi1(i)  !
         !! afft_mpi2(mp_sfftp), afft_mpi1(1:nfftp)
!!$            end do
!!$            call mpi_allreduce(afft_mpi2,afft_mpi3,mp_sfftp !!$
         !!      &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
!!$
!!$            if(j == mype) then
!!$!!$#ifdef _MPIFFTTEST_
!!$!!$               call m_FFT_set_scdata(afft_mpi3,afft)
!!$!!$#else
!!$               do i = ista_sfftp, iend_sfftp
!!$                  afft(i) = afft_mpi3(i - ista_sfftp + 1)
!!$               end do
!!$!!$#endif
!!$            end if
!!$         end do
!!$      else
!!$         afft = afft_mpi1
      end if

    end subroutine map_charge_onto_a_fft_box

    subroutine cp_afft_to_tchgrhr
      integer  :: i
      do i = ista_sfftph, iend_sfftph    ! MPI
         tchgr_l(i,iloop) = afft(i*2-1)
      end do

      if(ipripositron >= 2) then
         write(nfout,'(" tchgr_l")')
         do i = ista_sfftph, min(ista_sfftph+10,iend_sfftph)
            write(nfout,'(i3,f20.8)') i, tchgr_l(i,iloop)
         end do
      end if
    end subroutine cp_afft_to_tchgrhr

    subroutine check_lmn_even
      integer             :: nlmn, i

      do i = 1, 3
         nlmn = fft_box_size_CD(i,1)/2
         if(2*nlmn == fft_box_size_CD(i,1)) then
            lmn_even(i) = .true.
         else
            lmn_even(i) = .false.
         end if
      end do
    end subroutine check_lmn_even

    subroutine cpafft_CD(is)
      integer, intent(in) :: is
      integer             :: i

      afft_CD = 0.d0
      do i = ista_sfftph, iend_sfftph                  ! MPI
!!$         afft(2*i-1) = tchgr_l(i,is)
         afft_CD(2*i-1) = vepcr_l(i,is)
      end do
      if(npes > 1) then
         call mpi_allreduce(afft_CD, afft_mpi1, nfftp_nonpara , mpi_double_precision,&
              & mpi_sum, MPI_CommGroup,ierr)
         afft_CD = afft_mpi1
      end if

    end subroutine cpafft_CD

    subroutine cpafft_CD_to_vepc(is)
      integer, intent(in) :: is
      integer :: i, ik, ip


!!$      if(npes >= 2) then
!!$         call afft_allgatherv(afft,afft_mpi1)
!!$      else
!!$         afft_mpi1 = afft
!!$      end if

      do ik = 1, kimg
         do i = ista_kngp, min(kgp_reduced,iend_kngp)  !for mpi
#ifdef _MPIFFT_
            ip = (igfp_nonpara(i)-1)*kimg + ik
#elif _MPIFFTTEST_
            ip = (igfp_l_c(i)-1)*kimg + ik
#else
            ip = (igfp_l(i)-1)*kimg + ik
#endif
            vepc_l(i,ik,is)   = afft_CD(ip)*rinplw
         end do
      end do
    end subroutine cpafft_CD_to_vepc

  end subroutine m_epc_cal_potential

  subroutine check_of_pcc(nopcc)
    logical, intent(out) :: nopcc
    integer              :: i

    nopcc = .true.

    do i = 1, ntyp
       if(itpcc(i) /= 0) then
          nopcc = .false.
          goto 2201
       endif
    end do

    epc = 0.d0
2201 continue
  end subroutine check_of_pcc

  subroutine afft_allgatherv(afft_mpi0,afft)
    real(kind=DP),intent(in),dimension(:) :: afft_mpi0(ista_sfftp:iend_sfftp)
    real(kind=DP),intent(out),dimension(:):: afft(nfftp_nonpara)
    integer  :: id_sname = -1
    if(npes >= 2) then
       call tstatc0_begin('afft_allgatherv(in m_epc_potential) ',id_sname)
       call mpi_allgatherv(afft_mpi0,nel_sfftp(mype),mpi_double_precision &  ! MPI
            &       ,afft,nel_sfftp,idisp_sfftp,mpi_double_precision&
            &       ,MPI_CommGroup,ierr)
    else     ! npes == 1
       afft = afft_mpi0
    end if
    if(npes >= 2) call tstatc0_end(id_sname)
  end subroutine afft_allgatherv

! ==== POSITRON SCF ===== 2015/11/28
  subroutine m_epc_alloc_vlhxc_p
    if ( positron_method == positron_CONV ) then
       vlhxc_p => vlhxc_l
    else
       if ( .not. associated( vlhxc_p ) ) then
          allocate( vlhxc_p(ista_kngp:iend_kngp,kimg,nspin) )
          vlhxc_p = 0.0d0
       endif
    endif
  end subroutine m_epc_alloc_vlhxc_p

  subroutine m_epc_dealloc_vlhxc_p
    if ( positron_method /= positron_CONV ) then
       if ( associated( vlhxc_p ) ) deallocate( vlhxc_p )
    endif
  end subroutine m_epc_dealloc_vlhxc_p
! ======================= 2015/11/28

  subroutine m_epc_ESlhxc_potential(nfout)
    integer, intent(in)       :: nfout
!    real(kind=DP), intent(in) :: chg(ista_kngp:iend_kngp,kimg,nspin)
!    real(kind=DP), intent(in) :: vxc(ista_kngp:iend_kngp,kimg,nspin)

    integer :: is,ik,i,it
    integer :: ist
    integer :: id_sname = -1

!    call tstatc0_begin('m_ESlhxc_potential ',id_sname)

    vlhxc_p = 0.d0
    ist = ista_kngp
    if(ist == 1) ist = 2
    do is = 1, nspin
       do ik = 1, kimg
          if(mype==0) vlhxc_p(1,ik,is)   = vepc_l(1,ik,is)
          if(nspin == 1) then
             do i = ist, iend_kngp  !for mpi
                vlhxc_p(i,ik,is) = vepc_l(i,ik,is) -PAI4*tchgq_l(i,ik,is)&
                     &/gr_l(i)**2
             end do
          else if(nspin == 2) then
             do i = ist, iend_kngp  !for mpi
                vlhxc_p(i,ik,is) = vepc_l(i,ik,is) -PAI4*(tchgq_l(i,ik,UP)&
                     &+tchgq_l(i,ik,DOWN))/gr_l(i)**2
             end do
          endif
!          do it    = 1,ntyp
!             do i = ista_kngp, iend_kngp  !for mpi
!                vlhxc_p(i,ik,is) !                     & = vlhxc_p(i,ik,is)
          !+psc_l(i,it)*zfm3_l(i,it,ik)
!             end do
!          end do
       end do
    end do
!    call tstatc0_end(id_sname)
    if(ipripositron >= 3) then
       write(nfout,'(" -- vepc_l tchgq_l vlhxc_p << m_epc_ESlhxc_potential>&
            &>")')
       do is = 1, nspin
          if(nspin == 2) write(nfout,'("  ispin = ",i5)') is
          do ik = 1, kimg
             do i = ist, min(ist+20,iend_kngp)
                write(nfout,'(3f10.4)') vepc_l(i,ik,is),tchgq_l(i,ik,is)&
                     &,vlhxc_p(i,ik,is)
             end do
          end do
       end do
    end if
  end subroutine m_epc_ESlhxc_potential

! ==== POSITRON SCF ===== 2015/11/28
  subroutine m_epc_ESlhxc_potential_mod(nfout)    ! for scf
    integer, intent(in)       :: nfout

    integer :: is,ik,i,it
    integer :: ist
    integer :: id_sname = -1

    vlhxc_p = 0.d0
    ist = ista_kngp

    if(ist == 1) ist = 2
    do is = 1, nspin
       do ik = 1, kimg
          if(mype==0) vlhxc_p(1,ik,is)   = vepc_l(1,ik,is)
          if(nspin == 1) then
             do i = ist, iend_kngp  !for mpi
                vlhxc_p(i,ik,is) = vepc_l(i,ik,is) &
                     &            -PAI4 *chgq_l(i,ik,is) /gr_l(i)**2
             end do
          else if(nspin == 2) then
             do i = ist, iend_kngp  !for mpi
                vlhxc_p(i,ik,is) = vepc_l(i,ik,is) &
                     &           -PAI4*( chgq_l(i,ik,UP)&
                     &                  +chgq_l(i,ik,DOWN) )/gr_l(i)**2
             end do
          endif
          do it = 1,ntyp
             do i = ista_kngp, iend_kngp  !for mpi
                vlhxc_p(i,ik,is) = vlhxc_p(i,ik,is) &
                     &            -psc_l(i,it) *zfm3_l(i,it,ik)
             end do
          end do
       end do
    end do

  end subroutine m_epc_ESlhxc_potential_mod
! ============== 2015/11/28

end module m_epc_potential
