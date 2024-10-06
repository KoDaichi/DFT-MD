!=======================================================================
!
!  SOFTWARE NAME : PHASE (ver. 6.01)
!
!  MODULE: m_ELF
!
!  AUTHOR(S): T. Yamamoto   April/30/2004
!
!  Contact address :  Phase System Consortium
!                     E-mail: phase_system@nims.go.jp URL https://azuma.nims.go.jp
!
!
!
!=======================================================================
!
!
module m_ELF
!$$#ifndef PARA3D
! $Id: m_ELF.f90 376 2014-06-17 07:48:31Z jkoga $
  use m_PlaneWaveBasisSet,    only : ngabc,gr_l,igfp_l,kg,kgp,kgp_reduced, igf,iba,nbase,ngpt_l
  use m_PseudoPotential,      only : iatomn,ival
  use m_Crystal_Structure,    only : nopr,tau,altv,rltv,univol
  use m_Kpoints,              only : kv3
  use m_Ionic_System,         only : natm2, m_IS_pack_all_ions_in_uc
  use m_Charge_Density,       only : chgq_l
  use m_Electronic_Structure, only : occup_l,zaj_l
  use m_FFT,                  only : fft_box_size_CD,nfftp &
       &                           , m_FFT_CD_inverse &
       &                           , m_FFT_CD_direct &
       &                           , m_FFT_alloc_CD_box &
       &                           , m_FFT_dealloc_CD_box
  use m_FFT,                  only : fft_box_size_WF,nfft &
       &                           , m_FFT_WF &
       &                           , m_FFT_alloc_WF_work &
       &                           , m_FFT_dealloc_WF_work
  use m_Timing,               only : tstatc0_begin, tstatc0_end
  use m_Control_Parameters,   only : ipri,nspin,kimg,af,nel_Ylm,elf_filetype,elf_title
  use m_Const_Parameters,     only : BUCS,DP,DOWN,UP &
       &                           , PAI,PAI2,SKIP,DENSITY_ONLY,CUBE &
       &                           , DIRECT,INVERSE,OFF,ON &
       &                           , unit_conv_byname, ELECTRON
  use m_Parallelization,      only : MPI_CommGroup,ista_kngp,iend_kngp,npes,mype,ierr &
       &                           , ista_sfftp,iend_sfftp,ista_sfftph,iend_sfftph &
       &                           , nis_sfftp,nie_sfftp,nel_sfftp,idisp_sfftp,np_sfftp,mp_sfftp &
       &                           , map_k,myrank_k,np_e,ista_e,iend_e,istep_e,map_z
  use mpi
  implicit none

!  include 'mpif.h'

  real(kind=DP),private,allocatable, dimension(:)   :: chden_l ! d(ista_sfftp:iend_sfftp)
  real(kind=DP),private,allocatable, dimension(:)   :: keden_l ! d(ista_sfftp:iend_sfftp)
  real(kind=DP),private,allocatable, dimension(:)   :: grho2_l ! d(ista_sfftp:iend_sfftp)
  real(kind=DP),private,allocatable, dimension(:)   :: ita_l ! d(ista_sfftp:iend_sfftp)
  real(kind=DP),private,allocatable, dimension(:)   :: ita   ! d(nfftp)
  logical,      private, dimension(3)               :: lmn_even
  integer,private, allocatable, dimension(:) :: inx,jnx,knx ! d(ista_sfftp:iend_sfftp)

contains

  subroutine m_ELF_alloc()
    integer :: idp,nlp,nmp,nnp,nd2p,n,nn,i,j,k

    allocate(inx(ista_sfftp:iend_sfftp))
    allocate(jnx(ista_sfftp:iend_sfftp))
    allocate(knx(ista_sfftp:iend_sfftp))

    idp  = fft_box_size_CD(1,0)
    nlp  = fft_box_size_CD(1,1)
    nmp  = fft_box_size_CD(2,1)
    nnp  = fft_box_size_CD(3,1)
    nd2p = fft_box_size_CD(2,0)

    do n = ista_sfftp, iend_sfftp
       nn = (n+kimg-1)/kimg
       i  = mod(nn,idp)
       j  = mod((nn-1)/idp,nd2p) + 1
       k  = (nn - (j-1)*idp - i)/(idp*nd2p) + 1
       inx(n) = i - 1;  if(2*inx(n) > nlp) inx(n) = inx(n) - nlp
       jnx(n) = j - 1;  if(2*jnx(n) > nmp) jnx(n) = jnx(n) - nmp
       knx(n) = k - 1;  if(2*knx(n) > nnp) knx(n) = knx(n) - nnp
    end do

    allocate(chden_l(ista_sfftp:iend_sfftp))
    allocate(keden_l(ista_sfftp:iend_sfftp))
    allocate(grho2_l(ista_sfftp:iend_sfftp))
    allocate(ita_l(ista_sfftp:iend_sfftp))
    allocate(ita(nfftp))
  end subroutine m_ELF_alloc

  subroutine m_ELF_dealloc()
    deallocate(inx)
    deallocate(jnx)
    deallocate(knx)
    deallocate(chden_l)
    deallocate(keden_l)
    deallocate(grho2_l)
    deallocate(ita_l)
    deallocate(ita)
  end subroutine m_ELF_dealloc

  subroutine m_ELF_cnstrct_elf(nfout)
    integer, intent(in) :: nfout

    real(kind=DP), parameter :: cf = 7.15294982885971d0 ! 3.d-1*(3.d0*PAI**2)**(2.d0/3.d0)
    real(kind=DP), parameter :: third5 = 5.d0/3.d0
    real(kind=DP) :: d,tp,th,xs
    integer :: i

    call map_charge_onto_a_fft_box(nspin,chgq_l,chden_l) ! chgq_l -> chden_l
    call m_ELF_kinetic_energy_density(nfout,keden_l)
    call m_ELF_gradRho2(grho2_l)

    call m_FFT_alloc_CD_box
    call m_FFT_CD_inverse(nfout,chden_l)  ! G space --> R space
    call m_FFT_dealloc_CD_box

    ita_l = 0.d0
    do i=ista_sfftp,iend_sfftp
       d = chden_l(i)
       if(d >= 1.d-30) then
          tp = keden_l(i) - 0.125d0*grho2_l(i)/d
          th = cf*d**third5
          xs = tp/th
          ita_l(i) = 1.d0/(1.d0+xs*xs)
       end if
    end do

  end subroutine m_ELF_cnstrct_elf

  subroutine m_ELF_wd_elf(nfout,nfelf)
    integer, intent(in) :: nfout,nfelf

    integer :: i,j,k, idp, nlp, nmp, nnp, nlphf,inew,jnew,knew,ip,mmp
    real(kind=DP),allocatable,dimension(:,:,:) :: wkelf
    real(kind=DP) :: x,y,z
    real(kind=DP),allocatable,dimension(:,:) :: cps_full
    integer, allocatable,dimension(:) :: ityp_full
    real(kind=DP), dimension(3) :: r_wk
    integer :: ucret, m

    real(kind=DP) :: ita_mpi(nfftp)

    if(npes >= 2) then
       ita = 0.d0
       ita_mpi = 0.d0
       do i=ista_sfftp,iend_sfftp
          ita_mpi(i) = ita_l(i)
       end do
       call mpi_allreduce(ita_mpi,ita,nfftp,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr) ! MPI
    else
       ita = ita_l
    end if

    if(mype /= 0) return

    idp = fft_box_size_CD(1,0)
    mmp = fft_box_size_CD(2,0)
    nlp = fft_box_size_CD(1,1)
    nmp = fft_box_size_CD(2,1)
    nnp = fft_box_size_CD(3,1)

    if(kimg == 1) then
       nlphf = idp/2
    else
       nlphf = idp
    end if

    if(elf_filetype == DENSITY_ONLY) then
       allocate(wkelf(nlp,nmp,nnp)); wkelf = 0.d0
    else if(elf_filetype == CUBE) then
       allocate(wkelf(nnp,nmp,nlp)); wkelf = 0.d0
    end if

    write(nfout,9001) nlp*nmp*nnp, nlp, nmp, nnp
9001  format(' Electron Localization Function NE = ',i8,'(',3i5,')')

    write(nfout,*) ' !D FFT cube mapping start'
    do i = 1, nmp
       do j = 1, nnp
          do k = 1, nlp
             if(kimg == 1 .and. k > nlphf) then
                knew = idp - k
                jnew = nnp+2 - j
                inew = nmp+2 - i
                if(jnew > nnp) then
                   jnew = jnew - nnp
                end if
                if(inew > nmp) then
                   inew = inew - nmp
                end if
             else
                knew = k; jnew = j; inew = i
             end if
             ip = nlphf*mmp*(jnew-1) + nlphf*(inew-1) + knew
             if(elf_filetype == DENSITY_ONLY) then
                wkelf(k,i,j) = ita(ip*2-1)
             else if(elf_filetype == CUBE) then
                wkelf(j,i,k) = ita(ip*2-1)
             end if
          end do
       end do
    end do
    if(elf_filetype == DENSITY_ONLY) then
       write(nfelf,9001) nlp*nmp*nnp, nlp, nmp, nnp
       write(nfelf,'(6e13.5)') wkelf
    else if(elf_filetype == CUBE) then
       if(len_trim(elf_title) >= 1) then
          write(nfelf,*) trim(elf_title)
       else
          write(nfelf,'(" Calculated by phase")')
       end if
       write(nfelf,'(" SCF Total Density")')
       x = 0.d0; y = 0.d0; z = 0.d0
       write(nfelf,'(i6,3f10.4)') natm2, x,y,z
       do i = 1, 3
       !!   do m = 1, 3
       !!      ucret = unit_conv_byname( altv(m,i), r_wk(m), 'bohr', 'angstrom' )
       !!   end do
       !!   write(nfelf,'(i6,3f10.6)') fft_box_size_CD(i,1), r_wk(1:3)/dble(fft_box_size_CD(i,1))
          write(nfelf,'(i6,3f10.6)') fft_box_size_CD(i,1), altv(1:3,i)/dble(fft_box_size_CD(i,1))
       end do

       allocate(cps_full(natm2,3))
       allocate(ityp_full(natm2))
       call m_IS_pack_all_ions_in_uc(ityp_full,cps_full)
       do i = 1, natm2
          m = ityp_full(i)
          write(nfelf,'(f8.4,4f10.6)') iatomn(m), ival(m), cps_full(i,1:3)
       end do
       deallocate(ityp_full,cps_full)

       write(nfelf,'(6e13.5)') wkelf
       deallocate(wkelf)

    end if

  end subroutine m_ELF_wd_elf


  subroutine m_ELF_kinetic_energy_density(nfout,keden)
    integer, intent(in) :: nfout
    real(kind=DP), intent(out), dimension(ista_sfftp:iend_sfftp) :: keden

    real(kind=DP), dimension(nfft) :: afft,bfft
    real(kind=DP), dimension(ista_sfftp:iend_sfftp) :: afft_mpi
    real(kind=DP), dimension(ista_kngp:iend_kngp,kimg) :: kedq_l
    integer :: ik,ib1,in,i,i1,i2,ri,iend,ip
    real(kind=DP) :: occupation,fac

    afft = 0.d0

    ! debug
    !  write(*,*) 'debug: (1)'
    !  write(*,*) 'ista_e,iend_e,istep_e=',ista_e, iend_e, istep_e
    ! end debug

    call m_FFT_alloc_WF_work

    do ib1 = ista_e, iend_e, istep_e     ! MPI
       do ik=1,kv3
          if(map_k(ik) /= myrank_k) cycle ! MPI
          do in=1,3
            ! debug
            !  write(*,*) 'ib1,ik,in=',ib1,ik,in
            ! end debug
             call m_ELF_gradWF_in_Rspace(in,ik,ib1,bfft)
             occupation = occup_l(map_z(ib1),ik)
             do i = 1, nfft-1, 2
                afft(i) = afft(i) + occupation*(bfft(i)**2+bfft(i+1)**2) ! MPI
             end do
          end do
       end do
    end do
    if(npes >= 2) then
       call mpi_allreduce(afft,bfft,nfft,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr) ! MPI
       afft = bfft                          ! MPI
    end if

    call m_FFT_WF(ELECTRON,nfout,afft,DIRECT,ON)

    call m_FFT_dealloc_WF_work

! if nopr==1
!    keden = 0.d0
!    fac = 0.5d0*2.d0/(univol*kv3*product(fft_box_size_WF(1:3,1)))
!    do ri = 1, kimg
!       iend = iend_kngp
!       if( iend_kngp > kg ) iend = kg
!       if( ista_kngp <= iend ) then
!          do i = ista_kngp, iend  !for mpi
!             i1 = kimg*igf(i) + (ri - kimg)
!             ip = (igfp_l(i)-1)*kimg + ri
!             keden(ip) = afft(i1)*fac
!          end do
!       endif
!    end do
! end if

    fac = 0.5d0*2.d0/(univol*kv3*product(fft_box_size_WF(1:3,1)))
    do ri = 1, kimg
       iend = iend_kngp
       if( iend_kngp > kg ) iend = kg
       if( ista_kngp <= iend ) then
          do i = ista_kngp, iend  !for mpi
             i1 = kimg*igf(i) + (ri - kimg)
             kedq_l(i,ri) = afft(i1)*fac
          end do
       endif
    end do

    call ked_average(kedq_l)
    call map_charge_onto_a_fft_box(1,kedq_l,keden)

    call m_FFT_alloc_CD_box
    call m_FFT_CD_inverse(nfout,keden)  ! G space --> R space
    call m_FFT_dealloc_CD_box

  end subroutine m_ELF_kinetic_energy_density


  subroutine m_ELF_gradRho2(grho2)
    real(kind=DP), intent(out), dimension(ista_sfftp:iend_sfftp) :: grho2

    real(kind=DP) :: afft(ista_sfftp:iend_sfftp)
    integer :: in,i

    grho2 = 0.d0
    do in=1,3
       call m_ELF_gradRho_in_Rspace(in,afft)
       do i= ista_sfftp, iend_sfftp, 2
          grho2(i) = grho2(i) + afft(i)**2 + afft(i+1)**2
       end do
    end do

  end subroutine m_ELF_gradRho2


  subroutine map_charge_onto_a_fft_box(nspin,chgq_l,chden_l)
    integer, intent(in) :: nspin
    real(kind=DP), dimension(ista_kngp:iend_kngp,kimg,nspin) :: chgq_l
    real(kind=DP), dimension(ista_sfftp:iend_sfftp) :: chden_l

    real(kind=DP), dimension(nfftp) :: afft_mpi1
    real(kind=DP), dimension(mp_sfftp) :: afft_mpi2,afft_mpi3

    integer :: j, i, ip, is

    afft_mpi1 = 0.d0
    do is = 1, nspin
    do j = 1, kimg
       do i = ista_kngp, iend_kngp
          if(kgp_reduced <i ) cycle
          ip = (igfp_l(i)-1)*kimg + j
          afft_mpi1(ip) = afft_mpi1(ip) + chgq_l(i,j,is) !MPI
       end do
    end do
    end do

    if(npes >= 2) then
       call mpi_barrier(MPI_CommGroup,ierr)
       do j = 0, npes-1
          do i = nis_sfftp(j),nie_sfftp(j)
             afft_mpi2(i-nis_sfftp(j)+1) = afft_mpi1(i)
          end do
          call mpi_allreduce(afft_mpi2,afft_mpi3,mp_sfftp &
               &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)

          if(j == mype) then
             do i = ista_sfftp, iend_sfftp
                chden_l(i) = afft_mpi3(i - ista_sfftp + 1)
             end do
          end if
       end do
    else
       chden_l = afft_mpi1
    end if

  end subroutine map_charge_onto_a_fft_box


  subroutine m_ELF_gradWF_in_Rspace(in,ik,ib,bfft)
    use m_Files, only : nfout
    integer, intent(in)                           :: in, ik, ib
    real(kind=DP), intent(inout), dimension(nfft) :: bfft

    integer :: i,i1,ri,ip
    real(kind=DP) :: g

    bfft = 0.d0
    do ri = 1, kimg
       do i = 1, iba(ik)
          ip = nbase(i,ik)
          i1 = kimg*igf(ip) + (ri - kimg)
          g = rltv(in,1)*ngabc(ip,1)+rltv(in,2)*ngabc(ip,2)+rltv(in,3)*ngabc(ip,3)
          bfft(i1) = g*zaj_l(i,map_z(ib),ik,ri)   ! MPI
       end do
    end do
    if(ipri >= 2) then
       if(ik <= 2 .and. ib <= 3) then
          write(6,'(" ! bfft G-space in = ",i1," ik = ",i3," ib = ",i3," <<m_ELF_in_Rspace>>")') in, ik, ib
          write(6,'(8f8.4)') (bfft(i),i=1,120)
       end if
    end if
    call m_FFT_WF(ELECTRON,nfout,bfft,INVERSE,OFF)
    if(ipri >= 2) then
       if(ik <=  2 .and. ib <= 3) then
          write(6,'(" ! bfft R-space in = ",i1," ik = ",i3," ib = ",i3," <<m_ELF_in_Rspace>>")') in, ik, ib
          write(6,'(8f8.4)') (bfft(i),i=1,120)
       end if
    end if

  end subroutine m_ELF_gradWF_in_Rspace


  subroutine m_ELF_gradRho_in_Rspace(in,afft)
    use m_Files, only : nfout
    integer, intent(in) :: in
    real(kind=DP), intent(out) :: afft(ista_sfftp:iend_sfftp)

    integer       :: n
    real(kind=DP) :: gxyz

    afft = 0.d0
    do n = ista_sfftp, iend_sfftp  ! MPI
       gxyz = rltv(in,1)*inx(n)+rltv(in,2)*jnx(n)+rltv(in,3)*knx(n)
       afft(n) = afft(n) + gxyz*chden_l(n)
    enddo
    call check_lmn_even
    call boundary_zero(afft,in)

    call m_FFT_alloc_CD_box
    call m_FFT_CD_inverse(nfout,afft)  ! G space --> R space
    call m_FFT_dealloc_CD_box

  end subroutine m_ELF_gradRho_in_Rspace

  subroutine boundary_zero(afft,in)
    real(kind=DP), intent(inout) :: afft(ista_sfftp:iend_sfftp)
    integer, intent(in) :: in
    integer             :: i,j,k,nn,n,idp,nlp,nmp,nnp,nd2p

    idp = fft_box_size_CD(1,0)
    nlp = fft_box_size_CD(1,1)
    nmp = fft_box_size_CD(2,1)
    nd2p = fft_box_size_CD(2,0)
    nnp = fft_box_size_CD(3,1)

    if(lmn_even(in)) then
       if( in == 1) then
          do j = 1, nmp
             do k = 1, nnp
!!$                  nn = nlp/2 + 1 + idp*(j-1) + idp*nmp*(k-1)
                nn = nlp/2 + 1 + idp*(j-1) + idp*nd2p*(k-1)
                n  = nn*kimg - (kimg-1)
                if(n >= ista_sfftp .and. n <= iend_sfftp) afft(n) = 0.d0
                n  = nn*kimg
                if(n >= ista_sfftp .and. n <= iend_sfftp) afft(n) = 0.d0
             end do
          end do
       else if( in == 2) then
          do i = 1, idp
             do k = 1, nnp
!!$                  nn = i + idp*(nmp/2) + idp*nmp*(k-1)
                nn = i + idp*(nmp/2) + idp*nd2p*(k-1)
                n  = nn*kimg - (kimg-1)
                if(n >= ista_sfftp .and. n <= iend_sfftp) afft(n) = 0.d0
                n  = nn*kimg
                if(n >= ista_sfftp .and. n <= iend_sfftp) afft(n) = 0.d0
             end do
          end do
       else if(in == 3) then
          do i = 1, idp
             do j = 1, nmp
!!$                  nn = i + idp*(j-1) + idp*nmp*(nnp/2)
                nn = i + idp*(j-1) + idp*nd2p*(nnp/2)
                n  = nn*kimg - (kimg-1)
                if(n >= ista_sfftp .and. n <= iend_sfftp) afft(n) = 0.d0
                n  = nn*kimg
                if(n >= ista_sfftp .and. n <= iend_sfftp) afft(n) = 0.d0
             end do
          end do
       end if
    end if
  end subroutine boundary_zero

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

  subroutine ked_average(ked)
    real(kind=DP), intent(inout) :: ked(ista_kngp:iend_kngp,kimg)
    integer ::       ng, no, ngp, no1, no2
    real(kind=DP) :: fi, tx,ty,tz, fp, fc, fs, zcr, zci
    real(kind=DP), dimension(kgp,kimg) :: work
    real(kind=DP), pointer, dimension(:,:) :: work2

    allocate(work2(ista_kngp:iend_kngp,kimg)); work2 = 0.d0

    fi = 1.d0/dble(nopr)
    no1 = 1; no2 = nopr

    call cp_ked_to_work(ked,work) ! ked -> work
    work2 = 0.d0                  ! initialization
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
    ked(:,:) = work2(:,:)*fi

    deallocate(work2)

  contains

    subroutine cp_ked_to_work(ked,work)
      real(DP),intent(in),dimension(ista_kngp:iend_kngp,kimg) :: ked
      real(DP),intent(out),dimension(kgp,kimg) :: work
      integer :: ng,ri
      real(kind=DP), pointer, dimension(:,:) :: work_mpi
      allocate(work_mpi(kgp,kimg)); work_mpi = 0.d0
      work = 0.d0
      do ri = 1, kimg
         do ng = ista_kngp, iend_kngp  !for mpi
            work_mpi(ng,ri) = ked(ng,ri)
         end do
      end do
      call mpi_allreduce(work_mpi,work,kgp*kimg &
            &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
      deallocate(work_mpi)
     end subroutine cp_ked_to_work

  end subroutine ked_average

!$$#endif

end module m_ELF
