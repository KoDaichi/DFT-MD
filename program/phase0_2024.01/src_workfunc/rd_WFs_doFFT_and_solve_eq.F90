!================================================
!  Software name : STM
!  Subroutine(s) : rd_WFs_doFFT_and_solve_eq, alloc_afft_etc, dealloc_afft_etc,
!                  print_ck_etc,print_VLC_Rspace, sumup_afft_to, print_dz_etc,
!                  print_summative_charges
!  Author(s)     : Koichi Kato and Takahiro Yamasaki (June 7, 2004)
!
!  FURTHER MODIFICATION: Junichiro  Koga (June 24, 2004)
!
!  Contact address :  IIS,The University of Tokyo RSS21 project
!  
!  "Multiscale Simulation System for Function Analysis of Nanomaterials"  
!
!================================================
!
!     The original version of this program "STM" was developed in the
!  period 1998-2001 at JRCAT by Koichi Kato (Toshiba Corporate
!  Research and Development Center) and Takahiro Yamasaki (Fujitsu 
!  Laboratories Ltd.) through a joint research between JRCAT and
!  Toshiba Corporate Research and Development Center and another joint
!  research between JRCAT and FUJITSU Laboratories Ltd.
!     Since 2003, this program has been tuned and revised as it
!  matches to the first-principles molecular dynamics program "PHASE"
!  as a part of the national project "Frontier Simulation Software for
!  Industrial Science (FSIS)", which is supported by the IT program of
!  the Ministry of Education, Culture, Sports, Science and Technology
!  (MEXT) of Japan.
!     Since 2006, this program set has been developed as a part of the
!  national project "Revolutionary Simulation Software (RSS21)", which
!  is supported by the next-generation IT program of MEXT of Japan.

subroutine rd_WFs_doFFT_and_solve_eq
! $Id: rd_WFs_doFFT_and_solve_eq.F90,v 1.1.1.1 2004/06/26 11:21:37 yamasaki Exp $
  use m_Const_Parameters,     only : DP, BUCS, CARTS &
       &                           , HIGHER, LOWER, IN_BETWEEN, INVERSE, DIRECT &
       &                           , Hartree, NEGATIVE, BOHR
  use m_Control_Parameters,   only : e1, e2, izi, izf,rini,rfin,nfin,erlmt &
       &                           , n_fc,n_fc_i
  use m_Charge_File,          only : natm2,x_origin,y_origin,z_origin &
       &                           , fft_param,cell1,cell2,cell3 &
       &                           , iatomtype,ival,atom_pos
  use m_Kpoints, only :              kv3, nspin
  use m_Electronic_Structure, only : neordr,eko,efermi &
       &                           , m_ES_rd_WFs, m_ES_rd_VLCs &
       &                           , m_ES_WF_in_Rspace_fine_mesh &
       &                           , m_ES_VLC_in_Rspace_fine_mesh
  use m_FFT,                  only : nfft, nfftp, nfftpf, nlpf, fft_box_size_fine &
       &                           , m_FFT_WF_2dFFT_fine, m_FFT_Vlocal_W_fine
  use m_Crystal_Structure,    only : neg, zl
  use m_Files,                only : nfzaj, nfinp, nfout, nfchgu, nfchgd, nfvlc, nfchgu_p &
       &                           , nfchgd_p, nfvlcr,nfvlcr_av
  use m_Kpoints,              only : vkxyz
  use m_PlaneWaveBasisSet,    only : m_pwBS_kinetic_energies
  use m_ArraySize_Parameters, only : kspin, kimg

  implicit none
  integer :: i, kg2d, nls, it, ig, nlfid, n,n0, nmp, nnp
  real(kind=DP) :: ertmp, ermx, rmix, dz, ck, rdz, dz2, errsum
  real(kind=DP), allocatable, dimension(:,:,:) :: sfft, tfft
  real(kind=DP), allocatable, dimension(:,:,:) :: c, y
  real(kind=DP), allocatable, dimension(:,:)   :: am, w
  real(kind=DP), allocatable, dimension(:) :: ekin
  real(kind=DP), allocatable, dimension(:) :: bfft, cfft
  real(kind=DP), allocatable, dimension(:) :: ertmp_t, domin
  real(kind=DP), allocatable, dimension(:)  :: afft
  real(kind=DP) :: fac
  integer                                   :: id_sname = -1
  integer                                   :: ik, is, ib, ibto, ikp
  integer, dimension(3) :: tmpind

  nmp = fft_box_size_fine(2)
  nnp = fft_box_size_fine(3)

  allocate(bfft(nfftpf));bfft=0.d0
  call m_ES_rd_VLCs(nfvlc)               !-(m_E.S.)
  do is = 1, nspin
     call m_ES_VLC_in_Rspace_fine_mesh(is,bfft) !-(m_E.S.) \Vlc(G) -> \Vlc(r)
     call print_VLC_Rspace(is)             !-(contained here)
  enddo
  deallocate(bfft)

contains

  subroutine print_VLC_Rspace(is)
    integer, intent(in) :: is
    integer :: i, n, idp,nlp,nmp,nnp,nlphf, j,k, ip(4),knew,inew,jnew,nlph
    real(kind=DP), allocatable, dimension(:,:,:) :: wk,wk2
    real(kind=DP) :: vav

    idp = fft_box_size_fine(0)
    nlp = fft_box_size_fine(1); nmp = fft_box_size_fine(2); nnp = fft_box_size_fine(3)
    if(kimg.eq.1) then
       nlphf = idp/2
    else
       nlphf = idp
    end if
    
    allocate(wk(nlp,nmp,nnp))
    nlph = nlp
    if(kimg.eq.1) nlph = nlp/2
    do i = 1, nmp
       do j = 1, nnp
          do k = 1, nlp
             if(kimg.eq.1.and.k.gt.nlphf) then
                knew = idp -k
                jnew = nnp+2 -j
                inew = nmp+2 -i
                if(jnew.gt.nnp) then
                   jnew = jnew - nnp
                end if
                if(inew.gt.nmp) then
                   inew = inew - nmp
                end if
             else
                knew = k
                jnew = j
                inew = i
             end if
             ip(1) = nlphf*nmp*(jnew-1) + nlphf*(inew-1) + knew
             wk(k,i,j) = bfft(ip(1)*2-1) * Hartree ! to eV units
          end do
       end do
    end do

    write(nfvlcr_av,'(a,f10.5)') '# Fermi energy (eV) ',efermi*Hartree
    write(nfvlcr_av,'(a)')       '# distance along the z-axis(Angstrom)  averaged local potential (eV) '
    
    ! output the local potential in cube format
    do i=1,nlpf
       vav=0.d0
       do j=1,nnp
          do k=1,nmp
             vav = vav+wk(i,k,j)
          enddo
       enddo
       if ( n_fc(1).eq.1 ) then
           dz = cell1(1) * BOHR
       else if ( n_fc(1).eq.2 ) then
           dz = cell2(2) * BOHR
       else if ( n_fc(1).eq.3 ) then
           dz = cell3(3) * BOHR
       endif

       write(nfvlcr_av,'(f14.6,3e20.6)') dz*i, vav/dble(nnp*nmp)
    enddo

    write(nfvlcr,9001) nlpf*nmp*nnp,nmp,nnp,nlpf
    write(nfvlcr,*) ' -- local potential -- '

    if ( kimg .eq. 1 ) then
      fac = 2
    else
      fac = 1
    endif

    write(nfvlcr,'(i8,3f20.10)') nint(natm2/fac), x_origin,y_origin,z_origin
    tmpind(n_fc(1)) = nlpf
    tmpind(n_fc(2)) = nmp
    tmpind(n_fc(3)) = nnp

    write(nfvlcr,'(i8,3f20.10)') tmpind(1),cell1(1:3)
    write(nfvlcr,'(i8,3f20.10)') tmpind(2),cell2(1:3)
    write(nfvlcr,'(i8,3f20.10)') tmpind(3),cell3(1:3)

    do i=1,natm2/fac
      write(nfvlcr,'(i8,4f20.10)') iatomtype(i),ival(i),atom_pos(i,1:3)
    enddo

    allocate(wk2(tmpind(3),tmpind(2),tmpind(1)));wk2=0.d0

    do j=1,nnp
       do k=1,nmp
          do i=1,nlpf
             tmpind(n_fc_i(1)) = j
             tmpind(n_fc_i(2)) = k
             tmpind(n_fc_i(3)) = i
             wk2(tmpind(1),tmpind(2),tmpind(3)) = wk(i,k,j)
          enddo
       enddo
    enddo

    write(nfvlcr,9002) wk2

    deallocate(wk2)
    deallocate(wk)
9001 format(' CHARGE DENSITY NE = ',i8,'(',3i5,')')
9002 format(6e12.4)

  end subroutine print_VLC_Rspace

end subroutine rd_WFs_doFFT_and_solve_eq
