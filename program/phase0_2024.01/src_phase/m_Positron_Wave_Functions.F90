!=======================================================================
!
!  SOFTWARE NAME : PHASE (ver. 1200)
!
!  MODULE: m_Positron_Wave_Functions
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!
!    Further modification by T. Yamasaki   April/10/2007
!
!  Contact address :  Phase System Consortium
!                     E-mail: phase_system@nims.go.jp URL https://azuma.nims.go.jp
!
!
!
!=======================================================================
module m_Positron_Wave_Functions
! $Id: m_Positron_Wave_Functions.F90 472 2015-11-28 09:01:17Z ktagami $
  use m_Const_Parameters, only   : DP,POSITRON,OFF,ON,DIRECT,INVERSE, DELTAevdff, SD, MSD &
       &                         , DENSITY_ONLY, VTK, CUBE, SmallestPositiveNumber, BOHR, PAI
  use m_Control_Parameters, only : af, nspin, npeg,kimg,ipripositron &
       &                         , delta_pev, evaluation_pev_diff &
       &                         , num_extra_pev, sw_gga_p, sw_epsilon_ele, epsilon_ele &
       &                         , sw_positron_file,positron_filetype, positron_title &
       &                         , m_CtrlP_pstrn_ntcnvg_incre &
       &                         , m_CtrlP_pstrn_ntcnvg_reset &
       &                         , m_CtrlP_pstrn_ntcnvg_clear
  use m_IterationNumbers, only :   iteration_positron_wf
  use m_Timing, only :             tstatc0_begin, tstatc0_end
  use m_Files, only :              nfout, nfpstrn, nfvelec, nfeppair, nfeppair2, nfvelec_grad &
       &                         , m_Files_open_nfpstrn, m_Files_open_nfvelec &
       &                         , m_Files_open_nfeppair,m_Files_open_nfeppair2 &
       &                         , m_Files_open_nfvelec_grad &
       &                         , m_Files_close_nfpstrn,m_Files_close_nfvelec &
       &                         , m_Files_close_nfeppair,m_Files_close_nfeppair2 &
       &                         , m_Files_close_nfvelec_grad
  use m_Crystal_Structure, only  : univol
  use m_PlaneWaveBasisSet, only  : kg_pwf, kg1_pwf, igf_pstrn, igfp_l, igfp_nonpara &
       &                         , m_pwBS_pstrn_kinetic_energies
!!  use m_Electronic_Structure,only: vlhxc_l
  use m_Parallelization, only :    mype, npes, ista_kngp, iend_kngp &
       &                         , ista_sfftph, iend_sfftph, ierr, MPI_CommGroup
  use m_FFT, only   :              nfft, nfftp, nfftp_nonpara, nfft_pstrn, fft_box_size_pWF &
       &                         , fft_box_size_CD, fft_box_size_CD_nonpara &
       &                         , m_FFT_alloc_pWF_work, m_FFT_dealloc_WF_work &
       &                         , m_FFT_alloc_CD_box, m_FFT_dealloc_CD_box &
       &                         , m_FFT_Vlocal_pW &
       &                         , m_FFT_W_Vlocal_W &
       &                         , m_FFT_WF, m_FFT_CD_inverse0 &
       &                         , m_FFT_coef_CD_integration
  use m_Crystal_Structure,  only : altv
  use m_Ionic_System,         only : ntyp,ityp,zfm3_l,natm2, iatomn &
       &                         , m_IS_pack_all_ions_in_uc
  use m_epc_potential, only      : tchgr_l, grad_tchgr_l, vlhxc_p
  use m_PseudoPotential, only    : rhchg_l, ival
  use mpi

  implicit none
!  include 'mpif.h'

  real(kind=DP),allocatable,dimension(:,:,:) :: pzaj   ! positron wave functions d(kg1_pwf,npeg,kimg)
  real(kind=DP),allocatable,dimension(:,:,:) :: pzaj_old ! d(kg1_pwf,npeg,kimg)
  real(kind=DP),allocatable,dimension(:,:) ::   pchg_l ! positron charge in g-space, d(ista_kngp:iend_kngp,kimg)
  real(kind=DP),allocatable,dimension(:,:) ::   pchgo_l ! positron charge in g-space, d(ista_kngp:iend_kngp,kimg)
  real(kind=DP),allocatable,dimension(:) ::     pchr_l ! positron charge in r-space, d(ista_sfftph:iend_sfftph)
  integer, allocatable, dimension(:) ::         npeordr, nprvf_ordr !d(npeg)
  real(kind=DP),allocatable,dimension(:) ::     pev, pev1   ! d(npeg)
  real(kind=DP),private,allocatable,dimension(:) :: pevdff !d(3)

!!$  real(kind=DP),private,allocatable,dimension(:,:,:) :: vlhepc_l !d(ista_kngp:iend_kngp,kimg,nspin)
  real(kind=DP),private,allocatable, dimension(:) ::  afft, bfft, ekin, afft_mpi
  real(kind=DP),private,allocatable, dimension(:) ::  p
  real(kind=DP),private::valence_annihilation_rate,core_annihilation_rate
  real(kind=DP),private::valence_annihilation_rate_down,valence_annihilation_rate_total
  real(kind=DP),private::p_new_lifetime,p_new_lifetime_down,p_new_lifetime_total &
       &                ,p_old_lifetime,p_core_rate,p_core_rate_down

contains
  subroutine m_pWF_construct_pcharge()
    allocate(afft(nfft_pstrn)); allocate(bfft(nfft_pstrn))
    call m_FFT_alloc_pWF_work()
    pchg_l = 0.d0
    call m_pWF_WF_in_Rspace(npeordr(1),bfft)
    call pdensity() ! bfft -> afft
!!$    call m_FFT_pWF(nfout,afft,DIRECT,OFF) !  R-space -> G-space
    call m_FFT_WF(POSITRON, nfout,afft,DIRECT,OFF) !  R-space -> G-space
    call substitute_pCD_for_pchg()
    deallocate(afft); deallocate(bfft)
    call m_FFT_dealloc_WF_work()
  contains
    subroutine pdensity()
      integer :: i
      do i = 1, nfft_pstrn, 2
         afft(i) = bfft(i)**2 + bfft(i+1)**2
      end do
    end subroutine pdensity

    subroutine substitute_pCD_for_pchg()
      integer :: i, ri, i1, iend
      real(kind=DP) :: fac
      fac = 2.d0/(univol*product(fft_box_size_pWF(1:3,1)))
      do ri = 1, kimg
         iend = iend_kngp
         if( iend_kngp > kg_pwf) iend = kg_pwf
         if( ista_kngp <= iend ) then
            do i = ista_kngp, iend
               i1 = kimg*igf_pstrn(i) + (ri-kimg)
               pchg_l(i,ri) = afft(i1)*fac
            end do
         end if
      end do
#ifdef SINGLE_POSITRON
      if ( nspin == 1 ) pchg_l = pchg_l /2.0d0
#endif
    end subroutine substitute_pCD_for_pchg

  end subroutine m_pWF_construct_pcharge

  subroutine m_pWF_charge_rspace()
    integer :: i, j, ip

    call m_FFT_alloc_CD_box();    allocate(afft(nfftp_nonpara))
    if(npes >= 2) allocate(afft_mpi(nfftp_nonpara))

    afft = 0.d0
    do j = 1, kimg
       do i = ista_kngp, iend_kngp ! for mpi
#ifdef _MPIFFT_
          ip = (igfp_nonpara(i)-1)*kimg + j
#else
          ip = (igfp_l(i)-1)*kimg + j
#endif
          afft(ip) = afft(ip) + pchg_l(i,j)
       end do
    end do

    if(npes >= 2) then
       call mpi_allreduce(afft,afft_mpi,nfftp_nonpara,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
       afft = afft_mpi
    endif
    call m_FFT_CD_inverse0(nfout,afft)
    call cp_afft_to_pchr_l()             ! afft -> pchr_l

    call m_pWF_core_annihilation()
    call m_pWF_valence_annihilation()
    if(npes >=2) deallocate(afft_mpi)
    deallocate(afft);  call m_FFT_dealloc_CD_box()
  contains
    subroutine cp_afft_to_pchr_l()
      integer :: i
      do i = ista_sfftph, iend_sfftph
         pchr_l(i) = afft(i*2-1)
      end do
    end subroutine cp_afft_to_pchr_l
  end subroutine m_pWF_charge_rspace

  subroutine m_pWF_copy_pzaj_to_pzaj_old()
    pzaj_old = pzaj
  end subroutine m_pWF_copy_pzaj_to_pzaj_old

  subroutine m_pWF_allocate_pzaj_etc()
    integer :: ib
    allocate(pzaj(kg1_pwf,npeg,kimg))
    allocate(pzaj_old(kg1_pwf,npeg,kimg))
    allocate(npeordr(npeg))
    allocate(nprvf_ordr(npeg))
    npeordr(1:npeg) = (/(ib,ib=1,npeg)/)
    nprvf_ordr(1:npeg) = (/(ib,ib=1,npeg)/)
    allocate(pev(1:npeg))
    allocate(pev1(1:npeg)); pev1 = 0.d0
    allocate(pevdff(3))
    allocate(pchg_l(ista_kngp:iend_kngp,kimg))
    allocate(pchgo_l(ista_kngp:iend_kngp,kimg))
    allocate(pchr_l(ista_sfftph:iend_sfftph))
  end subroutine m_pWF_allocate_pzaj_etc

  subroutine m_pWF_deallocate_pzaj_etc()
    if(allocated(pchr_l)) deallocate(pchr_l)
    if(allocated(pchg_l)) deallocate(pchg_l)
    if(allocated(pchgo_l)) deallocate(pchgo_l)
    if(allocated(pevdff)) deallocate(pevdff)
    if(allocated(pev)) deallocate(pev)
    if(allocated(pev1)) deallocate(pev1)
    if(allocated(nprvf_ordr)) deallocate(nprvf_ordr)
    if(allocated(npeordr)) deallocate(npeordr)
    if(allocated(pzaj_old)) deallocate(pzaj_old)
    if(allocated(pzaj)) deallocate(pzaj)
  end subroutine m_pWF_deallocate_pzaj_etc

  subroutine m_pWF_alloc_afft_etc()
    allocate(afft(nfft_pstrn))
    allocate(bfft(nfft_pstrn))
    allocate(ekin(kg1_pwf))
    call m_FFT_alloc_pWF_work()
  end subroutine m_pWF_alloc_afft_etc

  subroutine m_pWF_dealloc_afft_etc()
    deallocate(afft,bfft,ekin)
    call m_FFT_dealloc_WF_work()
  end subroutine m_pWF_dealloc_afft_etc

  subroutine m_pWF_IW_by_randomnumbers()
    real(kind=DP)   :: a,b,p,xn
    data a,b,p/32771.d0,1234567891.d0,2147483648.d0/

    integer :: iimg, ieg, i
    integer :: id_sname = -1
    call tstatc0_begin('m_pWF_IW_by_randomnumbers ',id_sname)

    if(ipripositron >= 2) then
       write(nfout,*)
       write(nfout,'(" <<< m_pWF_IW_by_randomnumbers >>>")')
    end if

    pzaj = 0.d0

    xn = 0.d0
!!$    do iimg = 1, kimg
    do iimg = 1, 1
       do ieg = 1, npeg
          do i = 1,  kg1_pwf
             xn = mod(xn*a+b,p)
             pzaj(i,ieg,iimg) = xn/p
          enddo
       enddo
    enddo

    call tstatc0_end(id_sname)
  end subroutine m_pWF_IW_by_randomnumbers

  subroutine m_pWF_renew_WF_by_SDorCG(nfout,isolver,precon,dtim)
    integer, intent(in) ::       nfout,isolver,precon
    real(kind=DP), intent(in) :: dtim
    integer :: is, ib
    real(kind=DP) :: vlhxc0

    call m_pWF_alloc_afft_etc()
    allocate(p(kg1_pwf))

    call m_pwBS_pstrn_kinetic_energies(ekin)
!!$    do is = 1, nspin, af+1
    is = 1
       if(isolver == MSD) call vlhxc_p_zero_term(vlhxc0,is)
       call Vlocal_in_Rspace(is,afft)
       if(isolver == SD) then
          do ib = 1, npeg
             call m_pWF_WF_in_Rspace(ib,bfft)
             call m_FFT_Vlocal_pW(afft,bfft)
!!$             call m_FFT_pWF(nfout,bfft,DIRECT,ON)
             call m_FFT_WF(POSITRON,nfout,bfft,DIRECT,ON)
             call steepest_descent_p(precon,ib,dtim,ekin,bfft,p)
          end do
       else if(isolver == MSD) then
          do ib = 1, npeg
             call m_pWF_WF_in_Rspace(ib,bfft)
             call m_FFT_Vlocal_pW(afft,bfft)
!!$             call m_FFT_pWF(nfout,bfft,DIRECT,ON)
             call m_FFT_WF(POSITRON,nfout,bfft,DIRECT,ON)
             call modified_sd_p(precon,ib,dtim,vlhxc0,ekin,bfft,p)
          end do
       end if
       call m_pWF_modified_gram_schmidt()
!!$    end do

    call m_pWF_energy_eigen_values()
    call m_pWF_sort_eigen_values()
    call m_pWF_wd_pev(nfout)

    deallocate(p)
    call m_pWF_dealloc_afft_etc()
  contains
    subroutine vlhxc_p_zero_term(vlhxc0,ispin)
      real(kind=DP), intent(out) :: vlhxc0
      integer, intent(in)        :: ispin

      if(mype == 0) vlhxc0 = vlhxc_p(1,1,ispin)
      call mpi_bcast(vlhxc0,1,mpi_double_precision,0,MPI_CommGroup,ierr)
    end subroutine vlhxc_p_zero_term
  end subroutine m_pWF_renew_WF_by_SDorCG

  subroutine m_pWF_evolve_WFs_again(nfout,mode,dtim_old,dtim_new)
    integer, intent(in) :: nfout, mode ! ! mode = {ORTHONORMALIZATION | NORMALIZATION}
    real(kind=DP), intent(in) :: dtim_old, dtim_new
    integer :: is

    call m_pWF_alloc_afft_etc()
!!$    do is = 1, nspin, af+1
    call evolve_each_pWF_again(dtim_new,dtim_old)
    call m_pWF_modified_gram_schmidt()
!!$    end do
    call m_pWF_energy_eigen_values()
    call m_pWF_sort_eigen_values()
    call m_pWF_dealloc_afft_etc()
  end subroutine m_pWF_evolve_WFs_again

  subroutine evolve_each_pWF_again(dt_new,dt_old)
    real(kind=DP), intent(in) :: dt_new, dt_old
    integer :: ir, ib, ig
    real(kind=DP) :: dtt

    dtt = dt_new/dt_old
    do ir=1,kimg
       do ib=1, npeg
          do ig=1, kg1_pwf
             pzaj(ig,ib,ir) = (1-dtt)*pzaj_old(ig,ib,ir) + dtt*pzaj(ig,ib,ir)
          end do
       end do
    end do
  end subroutine evolve_each_pWF_again

  subroutine steepest_descent_p(precon,ib,dtim,ekin,VlocalpW,p)
    integer, intent(in) ::                             precon,ib
    real(kind=DP), intent(in) ::                       dtim
    real(kind=DP), intent(in), dimension(kg1_pwf) ::   ekin
    real(kind=DP), intent(in), dimension(nfft_pstrn) ::VlocalpW
    real(kind=DP),             dimension(kg1_pwf) ::   p

    integer :: i, i1
    real(kind=DP) :: evr,devr,denom,evi,e1,devi

    denom = 1.d0/product(fft_box_size_pWF(1:3,1))
    call decide_precon_factor_p(precon,ib,ekin,p) ! -> p(1:kg1_pwf)

    if(kimg==1) then
       do i = 1, kg1_pwf
          i1 = igf_pstrn(i)
          evr = pzaj(i,ib,1)
          devr = (ekin(i) - pev(ib))*evr + VlocalpW(i1)*denom
          pzaj(i,ib,1) = evr - p(i)*dtim*devr
       end do
    else if(kimg==2) then
       do i = 1, kg1_pwf
          i1 = igf_pstrn(i)
          evr = pzaj(i,ib,1); evi = pzaj(i,ib,kimg)
          e1 = ekin(i) - pev(ib)
          devr = e1*evr+VlocalpW(2*i1-1)*denom
          devi = e1*evi+VlocalpW(2*i1  )*denom
          pzaj(i,ib,1)    = evr - p(i)*dtim*devr
          pzaj(i,ib,kimg) = evi - p(i)*dtim*devi
       end do
    end if
  end subroutine steepest_descent_p

  subroutine modified_sd_p(precon,ib,dtim,vlhxc0,ekin,VlocalpW,p)
    integer, intent(in) ::                             precon,ib
    real(kind=DP), intent(in) ::                       dtim, vlhxc0
    real(kind=DP), intent(in), dimension(kg1_pwf) ::   ekin
    real(kind=DP), intent(in), dimension(nfft_pstrn) ::VlocalpW
    real(kind=DP),             dimension(kg1_pwf) ::   p

    integer :: i, i1
    real(kind=DP) :: evr,devr,denom,evi,e1,devi, wdi, fdexp

    denom = 1.d0/product(fft_box_size_pWF(1:3,1))
    call decide_precon_factor_p(precon,ib,ekin,p) ! -> p(1:kg1_pwf)

    if(kimg==1) then
       do i = 1, kg1_pwf
          i1 = igf_pstrn(i)
          evr = pzaj(i,ib,1)
          devr = (ekin(i) - pev(ib))*evr + VlocalpW(i1)*denom
          wdi = ekin(i) + vlhxc0 - pev(ib)
          if (dabs(wdi) < SmallestPositiveNumber) then
             pzaj(i,ib,1) = -p(i)*devr*dtim + evr
          else
             fdexp = dexp( -p(i) * wdi * dtim)
             pzaj(i,ib,1) = (fdexp - 1) * devr/wdi + evr
          endif
       end do
    else if(kimg==2) then
       do i = 1, kg1_pwf
          i1 = igf_pstrn(i)
          evr = pzaj(i,ib,1); evi = pzaj(i,ib,kimg)
          e1 = ekin(i) - pev(ib)
          devr = e1*evr+VlocalpW(2*i1-1)*denom
          devi = e1*evi+VlocalpW(2*i1  )*denom
          wdi  = ekin(i) + vlhxc0 - pev(ib)
          if (dabs(wdi) < SmallestPositiveNumber) then
             pzaj(i,ib,1) = -p(i)*devr*dtim + evr
             pzaj(i,ib,2) = -p(i)*devi*dtim + evi
          else
             fdexp = dexp( -p(i) * wdi * dtim)
             pzaj(i,ib,1)    = (fdexp -1)*devr/wdi + evr
             pzaj(i,ib,kimg) = (fdexp -1)*devi/wdi + evi
          endif
       end do
    end if
  end subroutine modified_sd_p

  subroutine m_pWF_submat(nfout)
    integer,intent(in) :: nfout
    integer :: is
    real(kind=DP) :: damp = 1.d0

    call m_pWF_alloc_afft_etc()

    call m_pwBS_pstrn_kinetic_energies(ekin)

!!$    do is = 1, nspin, af+1
    is = 1
    call Vlocal_in_Rspace(is,afft)  ! vlhxc_p -> afft
    call evolve_pWFs_in_subspace ! (is,npeg,damp,ekin,afft,bfft)
!!$    end do
    call m_pWF_dealloc_afft_etc()
  contains
    subroutine evolve_pWFs_in_subspace ! (is,npeg,damp,ekin,afft,bfft)
      integer :: ib, i, i1, ib1, ib2, ib2to, ibto, ib1to, m
      real(DP) :: denom, dr1, di1, dr2, di2
      real(DP), allocatable, dimension(:,:,:) :: pzah, pzaj_wk ! d(kg1_pwf,npeg,kimg)
      real(DP), allocatable, dimension(:,:)   :: zmat ! d(npeg*kimg,npeg)
      real(DP), allocatable, dimension(:)     :: w1hw2, eig

      allocate(pzah(kg1_pwf,npeg,kimg)); pzah = 0.d0
      denom = 1.d0/product(fft_box_size_pWF(1:3,1))
      do ib = 1, npeg
         ibto = nprvf_ordr(ib)
         if(ibto > npeg) cycle
         call m_pWF_WF_in_Rspace(ib,bfft)
         call m_FFT_Vlocal_pW(afft,bfft)
         call m_FFT_WF(POSITRON,nfout,bfft,DIRECT,ON)
         if(kimg == 1) then
            do i = 1, kg1_pwf
               i1 = igf_pstrn(i)
               pzah(i,ib,1) = ekin(i)*pzaj(i,ib,1)+bfft(i1)*denom
            end do
         else
            do i = 1, kg1_pwf
               i1 = igf_pstrn(i)
               pzah(i,ib,1) = ekin(i)*pzaj(i,ib,1) + bfft(2*i1-1)*denom
               pzah(i,ib,2) = ekin(i)*pzaj(i,ib,2) + bfft(2*i1  )*denom
            end do
         end if
      end do

      allocate(zmat(npeg*kimg,npeg))
      allocate(w1hw2(npeg*(npeg+1)/2*kimg)); w1hw2 = 0.d0
      allocate(eig(npeg)); eig = 0.d0

      do ib1 = 1, npeg
         if(ib1 > npeg) cycle
         ib1to = npeordr(ib1)
         do ib2 = 1, npeg
            ib2to = nprvf_ordr(ib2)
            if(ib2to > ib1) cycle
            m = ib1*(ib1-1)/2 + ib2to
            if(kimg == 1) then
               do i = 1, kg1_pwf
                  w1hw2(m) = w1hw2(m) + pzah(i,ib1to,1)*pzaj(i,ib2,1)
               end do
            else
               do i = 1, kg1_pwf
                  dr1 = pzah(i,ib1to,1)
                  di1 = pzah(i,ib1to,2)
                  dr2 = pzaj(i,ib2,1)
                  di2 = pzaj(i,ib2,2)
                  w1hw2(2*m-1)=w1hw2(2*m-1)+dr1*dr2+di1*di2
                  w1hw2(2*m  )=w1hw2(2*m  )+dr1*di2-di1*dr2
               end do
            end if
         end do
      end do

      call set_hmat(npeg,w1hw2,zmat)
      deallocate(w1hw2)

      if(kimg==1) then
         call dsyev_driver(npeg,eig,zmat)
      else
         call zheev_driver(npeg,eig,zmat)
      end if

      allocate(pzaj_wk(kg1_pwf,npeg,kimg)); pzaj_wk = 0.d0
      if(kimg == 1) then
         call subspace_rotation_real(1,npeg,pzaj_wk,zmat)
      else
         call subspace_rotation_imag(1,npeg,pzaj_wk,zmat)
      end if
      pzaj = pzaj_wk
      pev = eig
      deallocate(eig)
      deallocate(pzaj_wk)
      deallocate(zmat)
      deallocate(pzah)

     if(ipripositron >= 2) then
        write(nfout,*) ' !pstrn eigen values (pev) <<evolve_pWFs_in_subspace>>'
        write(nfout,'(5x,6f12.5)') (pev(ib1),ib1=1,npeg)
     endif

    end subroutine evolve_pWFs_in_subspace

    subroutine set_hmat(npeg,w1hw2,zmat)
      integer, intent(in) :: npeg
      real(kind=DP), intent(inout) :: w1hw2(*)
      real(kind=DP), intent(out) ::   zmat(npeg*kimg,npeg)
      real(kind=DP), allocatable, dimension(:) :: w1hw2_mpi

      integer :: i,j, m, msize
      msize = npeg*(npeg+1)/2
      if(kimg == 1) then
         do j = 1, npeg
            do i = 1, j
               m = j*(j-1)/2+i
               zmat(i,j)= w1hw2(m)
            end do
         end do
      else
         do j = 1, npeg
            do i = 1, j
               m = j*(j-1)/2+i
               zmat(2*i-1,j) =  w1hw2(2*m-1)
               zmat(2*i,  j) = -w1hw2(2*m)
            end do
         end do
      end if
    end subroutine set_hmat

    subroutine dsyev_driver(ndim,eig,w1hw2)
      integer, intent(in) :: ndim
      real(kind=DP), intent(out), dimension(ndim) :: eig
      real(kind=DP), intent(inout), dimension(ndim,ndim) :: w1hw2
      character(len=1) :: JOBZ,UPLO

      integer :: lda
      integer :: lwork,liwork
      complex(kind=DP),allocatable,dimension(:) :: work
      real(kind=DP),allocatable,dimension(:) :: rwork
      real(kind=DP),allocatable,dimension(:) :: iwork
      integer :: info

      lda=ndim
!(LAPACK)  JOBZ = N : eigenvalue, V : eigenvalue + eigenvector
      JOBZ = 'V'
!(LAPACK)  UPLO = U : upper triangle matrix,  L : lower triangle matrix
      UPLO = 'U'
! (define lwork)

      lwork  = 1 + ndim*(6+2*ndim)
      liwork = 3 + 5*ndim

      allocate(work(lwork)); work = 0.d0
      allocate(iwork(liwork)); iwork = 0.d0

      call dsyevd(JOBZ,UPLO,ndim,w1hw2,lda,eig,work,lwork,iwork,liwork,info)

      deallocate(work)
      deallocate(iwork)

    end subroutine dsyev_driver

    subroutine zheev_driver(ndim,eig,w1hw2)
      integer, intent(in):: ndim
      real(kind=DP), intent(out) ,dimension(ndim) :: eig
      real(kind=DP), intent(inout) ,dimension(ndim*2,ndim) :: w1hw2
      character(len=1) :: JOBZ,UPLO

      integer :: lda
      integer :: lwork,lrwork,liwork
      complex(kind=DP),allocatable,dimension(:) :: work
      real(kind=DP),allocatable,dimension(:) :: rwork
      real(kind=DP),allocatable,dimension(:) :: iwork
      integer :: info

      lda=ndim
!(LAPACK)  JOBZ = N : eigenvalue, V : eigenvalue + eigenvector
      JOBZ = 'V'
!(LAPACK)  UPLO = U : upper triangle matrix,  L : lower triangle matrix
      UPLO = 'U'
      lwork=ndim*(2+ndim)
      lrwork=1+ndim*(5+2*ndim)
      liwork=3+5*ndim

      allocate(work(lwork))
      allocate(rwork(lrwork))
      allocate(iwork(liwork))

      call zheevd(JOBZ,UPLO,ndim,w1hw2,lda,eig,work,lwork,rwork,lrwork,iwork,liwork,info)

      deallocate(work)
      deallocate(rwork)
      deallocate(iwork)

    end subroutine zheev_driver

    subroutine subspace_rotation_real(is,ie,zat_wk,zz)
      integer, intent(in) :: is,ie
      real(kind=DP), intent(out) :: zat_wk(kg1_pwf,npeg,kimg)
      real(kind=DP), intent(in) :: zz(npeg,is:ie)

      integer :: ib1to, ib1, ib2, ii
      do ib2=is,ie
         do ib1=1,npeg
            ib1to=npeordr(ib1)
            do ii=1,kg1_pwf
               zat_wk(ii,ib2,1)=zat_wk(ii,ib2,1)+zz(ib1,ib2)*pzaj(ii,ib1to,1)
            enddo
         enddo
      enddo
    end subroutine subspace_rotation_real

    subroutine subspace_rotation_imag(is,ie,zat_wk,zz)
      integer, intent(in) :: is,ie
      real(kind=DP), intent(out) :: zat_wk(kg1_pwf,npeg,kimg)
      real(kind=DP), intent(in) :: zz(npeg*kimg,is:ie)
      integer :: ib1to, ib1, ib2, ii
      real(kind=DP) :: dr1, di1, dr2, di2

      do ib2=is,ie
         do ib1=1,npeg
            ib1to = npeordr(ib1)
            dr1=zz(2*ib1-1,ib2)
            di1=zz(2*ib1  ,ib2)
            do ii=1,kg1_pwf
               dr2=pzaj(ii,ib1to,1)
               di2=pzaj(ii,ib1to,2)
               zat_wk(ii,ib2,1)=zat_wk(ii,ib2,1)+dr1*dr2-di1*di2
               zat_wk(ii,ib2,2)=zat_wk(ii,ib2,2)+dr1*di2+di1*dr2
            enddo
         enddo! ib1 loop
      enddo! ib2 loop
    end subroutine subspace_rotation_imag
  end subroutine m_pWF_submat

  subroutine decide_precon_factor_p(precon,ib,ekin,p)
    integer,intent(in) ::                           precon,ib
    real(kind=DP),intent(in),dimension(kg1_pwf)  :: ekin
    real(kind=DP),intent(out),dimension(kg1_pwf) :: p

    integer  :: i
    real(kind=DP) :: ektot, x, x1, x2, d_ektot

    if(precon == ON) then
       call kinetic_energy_p(ib,ekin,ektot)
       d_ektot = 1.d0/ektot
       do i = 1, kg1_pwf
          x = ekin(i)*d_ektot
          x1 = 27 + (18 + (12 + 8*x)*x)*x
          x2 = 16*(x*x)*(x*x)
          p(i) = x1/(x1 + x2)
       end do
    else
       p = 1.d0
    end if
  end subroutine decide_precon_factor_p

  subroutine kinetic_energy_p(ib,dekin,ektot)
    integer, intent(in) :: ib
    real(kind=DP), intent(in), dimension(kg1_pwf) :: dekin
    real(kind=DP), intent(out) ::                    ektot
    integer :: i, ri
    ektot = 0.d0
    do ri = 1, kimg
       do i = 1, kg1_pwf
          ektot = ektot + dekin(i)*pzaj(i,ib,ri)*pzaj(i,ib,ri)
       end do
    end do
  end subroutine kinetic_energy_p


  subroutine m_pWF_modified_gram_schmidt()
    integer :: i, ito, nmax
    real(kind=DP) :: fr
    real(kind=DP),pointer, dimension(:,:) :: p1p2
    integer,      pointer, dimension(:)   :: ib2to_a, ib2back_a

    allocate(ib2to_a(npeg), ib2back_a(npeg))
    allocate(p1p2(npeg,kimg))
    do i = 1, npeg
       ito = npeordr(i)
       call WW(ito,fr)
       call normalize_pzaj(ito,fr)
       call substitute_jto_ib2back(i,nmax)
       call W1W2(ito)
       call modify_pzaj(i,ito)
    end do
    deallocate(p1p2,ib2back_a,ib2to_a)

  contains
    subroutine WW(j,fr)
      integer, intent(in) ::        j
      real(kind=DP), intent(out) :: fr

      integer ::       ri, ig
      real(kind=DP) :: fr1, fi1
      if(kimg == 1) then
         fr1 = 0.d0
         do ig = 1, kg1_pwf
            fr1 = fr1 + pzaj(ig,j,1)*pzaj(ig,j,1)
         end do
         fr = fr1
      else
         fr1 = 0.d0; fi1 = 0.d0
         do ig = 1, kg1_pwf
!!$            fr1 = fr1 + pzaj(ig,j,1)**2+pzaj(ig,j,2)**2
            fr1 = fr1 + pzaj(ig,j,1)*pzaj(ig,j,1)
            fi1 = fi1 + pzaj(ig,j,2)*pzaj(ig,j,2)
         end do
         fr = fr1 + fi1
      end if

      if(fr < 0.d0) then
         if(ipripositron >= 1) then
            write(nfout,'(" i, ito, j = ",3i8)') i, ito, j
            write(nfout,'(" fr = ",d23.10)') fr
            if(kimg == 2) write(nfout,'(" fr1, fi1 = ",2d23.10)') fr1, fi1
         end if
         call phase_error_with_msg(nfout,' fr < 0 <<m_pWF_modified_gram_schmidt.WW>>',__LINE__,__FILE__)
      else
         fr = 1.d0/dsqrt(fr)
      end if
    end subroutine WW

    subroutine normalize_pzaj(ibo,fr)
      integer, intent(in) :: ibo
      real(kind=DP), intent(in) :: fr
      integer :: ri, i
      do ri = 1, kimg
         do i = 1, kg1_pwf
            pzaj(i,ibo,ri) = fr*pzaj(i,ibo,ri)
         end do
      end do
    end subroutine normalize_pzaj

    subroutine substitute_jto_ib2back(i,nmax)
      integer, intent(in)  :: i
      integer, intent(out) :: nmax

      integer              ::jto, j
      jto = 0
      do j = 1, npeg
         if(nprvf_ordr(j) <= i) cycle
         jto = jto + 1
         ib2to_a(jto)  = j
         ib2back_a(j)  = jto
      end do
      nmax = jto
    end subroutine substitute_jto_ib2back

    subroutine W1W2(ito)
      integer, intent(in) :: ito

      real(kind=DP) :: ar, ai
      integer :: jto, j, ia

      p1p2 = 0.d0
      if(kimg == 1) then
         do jto = 1, nmax
            j = ib2to_a(jto)
            do ia = 1, kg1_pwf
               p1p2(jto,1) = p1p2(jto,1) + pzaj(ia,ito,1)*pzaj(ia,j,1)
            end do
         end do
      else if(kimg == 2) then
         do jto = 1, nmax
            j = ib2to_a(jto)
            do ia = 1, kg1_pwf
               ar = pzaj(ia,ito,1)
               ai = pzaj(ia,ito,2)
               p1p2(jto,1) = p1p2(jto,1) + ar*pzaj(ia,j,1)+ai*pzaj(ia,j,2)
               p1p2(jto,2) = p1p2(jto,2) + ar*pzaj(ia,j,2)-ai*pzaj(ia,j,1)
            end do
         end do
      end if
    end subroutine W1W2

    subroutine modify_pzaj(i,ito)

      integer, intent(in) :: i, ito
      integer :: j, ia, jto
      real(kind=DP) :: sr, si
      do j = 1, npeg
         if(nprvf_ordr(j) <= i) cycle
         jto = ib2back_a(j)
         if(kimg == 1) then
            do ia = 1, kg1_pwf
               pzaj(ia,j,1) = pzaj(ia,j,1) - p1p2(jto,1)*pzaj(ia,ito,1)
            end do
         else if(kimg == 2) then
            do ia = 1, kg1_pwf
               sr = pzaj(ia,ito,1);         si = pzaj(ia,ito,2)
               pzaj(ia,j,1) = pzaj(ia,j,1) - p1p2(jto,1)*sr+p1p2(jto,2)*si
               pzaj(ia,j,2) = pzaj(ia,j,2) - p1p2(jto,1)*si-p1p2(jto,2)*sr
            end do
         end if
      end do
    end subroutine modify_pzaj

  end subroutine m_pWF_modified_gram_schmidt

  logical function m_pWF_pevdff()
    integer :: ib
    real(kind=DP) :: fac, tmp
    if(evaluation_pev_diff == OFF) then
       m_pWF_pevdff = .false.
       return
    end if

    pevdff = 0.d0
    fac = 1.d0/npeg
    do ib = 1, npeg
       if(nprvf_ordr(ib) > npeg- num_extra_pev) cycle
       tmp = pev(ib) - pev1(ib)
       pevdff(1) = pevdff(1) + tmp*tmp
       pevdff(2) = pevdff(2) + dabs(tmp)
       if(dabs(tmp) > DELTAevdff) &
            & pevdff(3) = pevdff(3) + dsqrt(dabs(pev1(ib)**2 - pev(ib)**2))
    end do

    pevdff(1) = dsqrt(fac*pevdff(1))
    pevdff(2) = fac*pevdff(2)
    pevdff(3) = fac*pevdff(3)

    pev1 = pev

    if(ipripositron >= 1) write(nfout,*)'lifetime: ',p_old_lifetime,p_new_lifetime
    if(max(pevdff(1),pevdff(2)) < delta_pev.and. &
         & dabs(p_old_lifetime-p_new_lifetime) < 1.d0  ) then
       ib = m_CtrlP_pstrn_ntcnvg_incre()
       if(ipripositron >= 1) &
            & write(nfout,'(" !iter = ",i7," ntcnvg = ",i7)') iteration_positron_wf,ib
    else
       call m_CtrlP_pstrn_ntcnvg_reset()
    end if
    if(ipripositron >= 1) &
         & write(nfout,'(" >> (",i6,") <pev_old-pev_new>:(",3d13.5,")")') &
         &    iteration_positron_wf, pevdff(1), pevdff(2), pevdff(3)

    p_old_lifetime=p_new_lifetime
    m_pWF_pevdff = m_CtrlP_pstrn_ntcnvg_clear()

  end function m_pWF_pevdff

  subroutine m_pWF_wd_pev(nf)
    integer, intent(in) :: nf
    integer :: ib

    if(ipripositron >= 1) then
       write(nf,*) '=== positron eigen values ==='
       write(nf,'(99(4f18.10,/))') (pev(npeordr(ib)),ib=1,npeg-num_extra_pev)
       write(nfout,'(" -- extra_bands --",/99(4f18.10,/))') &
            & (pev(npeordr(ib)),ib=npeg-num_extra_pev+1, npeg)
    end if
  end subroutine m_pWF_wd_pev

!!$  subroutine m_pWF_wd_pev(nfout)
!!$    integer, intent(in) :: nfout
!!$
!!$    integer :: ib
!!$    write(nfout,'(" -- Positron Energy Eigen Values --")')
!!$    write(nfout,'(" ",99(10f8.4,/))') (pev(npeordr(ib)),ib=1,npeg)
!!$  end subroutine m_pWF_wd_pev

  subroutine m_pWF_energy_eigen_values()
    integer :: is, ib
    real(kind=DP) :: eg

!!$    allocate(afft(nfft_pstrn))
!!$    allocate(bfft(nfft_pstrn))
!!$    allocate(ekin(kg1_pwf))
!!$
!!$    call m_FFT_alloc_pWF_work()
    call m_pwBS_pstrn_kinetic_energies(ekin)
!!$    do is = 1, nspin, af+1
    is = 1
       call Vlocal_in_Rspace(is,afft) ! vlhxc_p -> afft
       if(ipripositron >= 3) then
          write(nfout,'(" afft <<m_pWF_energy_eigen_values>>")')
          write(nfout,'( 8f8.4)') (afft(ib),ib=1,100)
       end if
       do ib = 1, npeg
          call m_pWF_W_T_W(ib)  ! -> pev
          call m_pWF_WF_in_Rspace(ib,bfft)
          call m_FFT_W_Vlocal_W(POSITRON,nfft_pstrn,afft,bfft,eg)
          if(ipripositron >= 2) then
             if(ib==1) write(nfout,'(" ---- ib,  pev, eg, pev+eg <<m_pWF_energy_eigen_values>>")')
             write(nfout,'(i5,3f15.7)') ib, pev(ib),eg,pev(ib)+eg
          end if
          pev(ib) = pev(ib) + eg
       end do
!!$    end do

!!$    call m_FFT_dealloc_WF_work()
!!$    deallocate(ekin)
!!$    deallocate(bfft)
!!$    deallocate(afft)
  end subroutine m_pWF_energy_eigen_values

  subroutine m_pWF_sort_eigen_values()
    integer :: ib, jb, ibo, jbo
    real(kind=DP), parameter :: delta = 1.d-12

    npeordr(1:npeg) = (/(ib,ib=1,npeg)/)
    do ib = 1, npeg-1
       do jb = ib+1, npeg
          ibo = npeordr(ib)
          jbo = npeordr(jb)
          if(pev(jbo) < pev(ibo) - delta) then
             npeordr(jb) = ibo
             npeordr(ib) = jbo
          end if
       end do
    end do
    do ib = 1, npeg
       do jb = 1, npeg
          if(ib == npeordr(jb)) then
             nprvf_ordr(ib) = jb
             exit
          end if
       end do
    end do

  end subroutine m_pWF_sort_eigen_values

  subroutine Vlocal_in_Rspace(is,afft)
    integer, intent(in) :: is
    real(kind=DP), intent(out), dimension(nfft_pstrn) :: afft
    integer :: i
    call map_vlhxc_p_to_afft()
    if(ipripositron >= 3) then
       write(nfout,'(" afft <<Vlocal_in_Rspace>>")')
       write(nfout,'( 8f8.4)') (afft(i),i=1,100)
    end if
!!$    call m_FFT_pWF(nfout,afft,INVERSE,OFF)  ! afft -> afft
    call m_FFT_WF(POSITRON,nfout,afft,INVERSE,OFF)  ! afft -> afft
  contains
    subroutine map_vlhxc_p_to_afft()
      integer :: i,i1,ri, iend
      if(npes >= 2)  allocate(afft_mpi(nfft_pstrn))
      afft = 0.d0
      iend = iend_kngp
      if( iend > kg_pwf) iend = kg_pwf
      do ri = 1, kimg
         do i = ista_kngp, iend
            i1 = kimg*igf_pstrn(i) + (ri-kimg)
            afft(i1) = vlhxc_p(i,ri,is)
         end do
      end do
      if(ipripositron >= 3) then
         write(nfout,'(" -- vlhxc_p <<map_vlhxc_p_to_afft>>")')
         do ri = 1, kimg
            write(nfout,'(" kimg = ", i8)') kimg
            write(nfout,'(10f8.4)') (vlhxc_p(i,ri,is),i=ista_kngp,iend)
         end do
      end if
      if(npes >= 2) then
         call mpi_allreduce(afft,afft_mpi,nfft_pstrn &
              & , mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
         afft = afft_mpi
      end if
      if(npes >= 2) deallocate(afft_mpi)
    end subroutine map_vlhxc_p_to_afft

  end subroutine Vlocal_in_Rspace

  subroutine m_pWF_W_T_W(ib)
    integer, intent(in) :: ib
    integer :: i, ri
    pev(ib) = 0.d0
    do ri = 1, kimg
       do i = 1, kg1_pwf
          pev(ib) = pev(ib) + ekin(i)*pzaj(i,ib,ri)**2
       end do
    end do
  end subroutine m_pWF_W_T_W

  subroutine m_pWF_WF_in_Rspace(ib,bfft)
    integer, intent(in) :: ib
    real(kind=DP), intent(inout), dimension(nfft_pstrn) :: bfft

    integer :: i,i1,ri
    bfft = 0.d0
    do ri = 1, kimg
       do i = 1, kg1_pwf
          i1 = kimg*igf_pstrn(i) + (ri-kimg)
          bfft(i1) = pzaj(i,ib,ri)
       end do
    end do
!!$    call m_FFT_pWF(nfout,bfft,INVERSE,ON)
    call m_FFT_WF(POSITRON,nfout,bfft,INVERSE,ON)
    if(ipripositron >= 3) then
       write(nfout,'(" bfft <<m_pWF_WF_in_Rspace>>")')
       write(nfout,'( 8f8.4)') (bfft(i),i=1,100)
    end if
  end subroutine m_pWF_WF_in_Rspace

  real(kind=DP) function m_pWF_tell_band_energy()
    integer :: ib
    real(kind=DP) :: eband
    eband = 0.d0
    do ib = 1, npeg
       if(nprvf_ordr(ib) > npeg - num_extra_pev ) cycle
       eband = eband + pev(ib)
    end do
    m_pWF_tell_band_energy = eband
  end function m_pWF_tell_band_energy

  subroutine m_pWF_core_annihilation()
    integer :: it, ri, i
    real(kind=DP) :: sss, sss_mpi
!!$    integer :: i, j, ip,idp,mmp,mnp,ipp,k,nlp,nmp,nnp,up_down &
!!$         ,inew,jnew,knew, nlphf,nmesh,n,ri,it
!!$    real(kind=DP) ::fac,fac2
!!$    real(kind=DP),allocatable,dimension(:,:,:) ::wkchr
!!$    real(kind=DP),allocatable,dimension(:,:,:) ::wkchr2
!!$    real(kind=DP),allocatable,dimension(:) ::gr_l
!!$    real(kind=DP),allocatable,dimension(:) ::radr,wos,rhcr,wky,wkx
!!$    real(kind=DP) ::sss,ssk,enhanc,sss1,sss2,gabs
!!$
!!$    rewind 501

    sss=0.d0
    do it=1,ntyp
!!$       read(501,*)nmesh
!!$       allocate(radr(1:nmesh),wos(1:nmesh),rhcr(1:nmesh),gr_l(ista_kngp:iend_kngp))
!!$       allocate(wky(1:nmesh),wkx(1:nmesh))
!!$       read(501,*)radr(1:nmesh)
!!$       read(501,*)wos(1:nmesh)
!!$       read(501,*)rhcr(1:nmesh)
!!$       read(501,*)wky(1:nmesh)
!!$       read(501,*)gr_l(ista_kngp:iend_kngp)
       do ri = 1, kimg
          do i = ista_kngp, iend_kngp
             sss = sss + rhchg_l(i,it)*pchg_l(i,ri)*zfm3_l(i,it,ri)
          end do
       end do

       if(ipripositron >= 3) then
          if(it == 1) write(nfout,'(" -- rhchg_l, pchgq_l, zfm3_l <<m_pWF_core_annihilation>>")')
          write(nfout,'(" -- it = ",i8)') it
          write(nfout,'(" --- i, rhchg_l pchg_l zfm3_l --")')
          do i = ista_kngp, min(ista_kngp+max(10,20*ipripositron), iend_kngp)
             write(nfout,'(i8,3f20.8)') i, rhchg_l(i,it), pchg_l(i,1), zfm3_l(i,it,1)
          end do
          write(nfout,'(" sss = ",f20.8)') sss
       end if
!!$          gabs = gr_l(i)
!!$          wkx(1:nmesh) = gabs*radr(1:nmesh)
!!$          call dsjnv(0,nmesh,wkx,wky)
!!$          do ri=1,kimg
!!$             do n = 1, nmesh
!!$                fac2=rhcr(n)/radr(n)/4.d0/PAI
!!$                call enhance_0(fac2,enhanc)
!!$                fac = wos(n)*rhcr(n)/univol*enhanc
!!$                sss=sss+fac*wky(n)* pchg_l(i,ri)*zfm3_l(i,it,ri)
!!$!                sss=sss+fac*wky(n)* pchg_l(i,ri)*2.d0
!!$             end do
!!$          enddo
!!$       end do
    enddo

    if(npes >= 2) then
       call mpi_allreduce(sss,sss_mpi,1,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
       sss = sss_mpi
    end if
!                       valence_annihilation_rate=sss1
!                       ssk=((BOHR)**3)*1.d4/(2.8d0**2)/3.d0/PAI
!    core_annihilation_rate=sss/dsqrt(4.d0*PAI)
     core_annihilation_rate=sss/2.d0*univol

#ifdef SINGLE_POSITRON
     core_annihilation_rate = sss *univol
#endif

  end subroutine m_pWF_core_annihilation

  subroutine m_pWF_valence_annihilation()
    integer :: i, j, ip,idp,mmp,mnp,ipp,k,nlp,nmp,nnp,up_down &
         & ,inew,jnew,knew, nlphf, nd1,nd2,nd3

    real(kind=DP),allocatable,dimension(:)     :: f2or1
    real(kind=DP),allocatable,dimension(:,:,:) :: wkchr
    real(kind=DP),allocatable,dimension(:,:,:) :: wkchr2,wkchr3, wkchr4,wkchr22,wkchr32
    real(kind=DP) ::sss,sss_down,sss_total,ssk,enhanc,sss1,sss2, sss_mpi
    allocate(f2or1(ista_sfftph:iend_sfftph))
    call m_FFT_coef_CD_integration(f2or1)

#ifdef FFTW3
    idp = fft_box_size_CD(1,1)
    mmp = fft_box_size_CD(2,1)
    mnp = fft_box_size_CD(3,1)
#else
    idp = fft_box_size_CD(1,0)
    mmp = fft_box_size_CD(2,0)
    mnp = fft_box_size_CD(3,0)
#endif
    nlp = fft_box_size_CD(1,1)
    nmp = fft_box_size_CD(2,1)
    nnp = fft_box_size_CD(3,1)

    if(sw_positron_file == ON) then

    if(kimg == 1) then
       nlphf = idp/2
    else
       nlphf = idp
    end if

    if(positron_filetype == DENSITY_ONLY .or. positron_filetype == VTK) then
       nd1 = nlp; nd2 = nmp; nd3 = nnp
    else if(positron_filetype == CUBE) then
       nd1 = nnp; nd2 = nmp; nd3 = nlp
    end if
    allocate(wkchr(nd1,nd2,nd3)); wkchr = 0.d0
    allocate(wkchr2(nd1,nd2,nd3)); wkchr2 = 0.d0
    allocate(wkchr3(nd1,nd2,nd3)); wkchr3 = 0.d0
    if(nspin==2) then
       allocate(wkchr22(nd1,nd2,nd3)); wkchr22 = 0.d0
       allocate(wkchr32(nd1,nd2,nd3)); wkchr32 = 0.d0
    end if
    if(sw_gga_p == ON) then
       allocate(wkchr4(nd1,nd2,nd3)); wkchr4 = 0.d0
    end if

    up_down=1
    do i = 1, nmp
       do j = 1, nnp
          do k = 1,nlp
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
             ip=nlphf*mmp*(jnew-1) + nlphf*(inew-1) + knew
             if(positron_filetype == DENSITY_ONLY .or. positron_filetype == VTK) then
!                wkchr(k,i,j) = afft(ip*2-2+up_down)
                wkchr(k,i,j) = afft(ip*2-2+up_down)/2
                if(ip >= ista_sfftph .and. ip <= iend_sfftph) then
                   wkchr2(k,i,j) = tchgr_l(ip,1)
                   if(nspin==2) wkchr22(k,i,j) = tchgr_l(ip,2)
                   if(sw_gga_p == ON) wkchr4(k,i,j) = grad_tchgr_l(ip,1)
                end if
             else if(positron_filetype == CUBE) then
! lin
!             if(positron_filetype == DENSITY_ONLY .or. positron_filetype == VTK) then

!                wkchr(k,i,j) = afft(ip*2-2+up_down)
                if(ip*2-2+up_down<1) cycle
                wkchr(j,i,k) = afft(ip*2-2+up_down)/2
                if(ip >= ista_sfftph .and. ip <= iend_sfftph) then
                   wkchr2(j,i,k) = tchgr_l(ip,1)
                   if(nspin==2) wkchr22(j,i,k) = tchgr_l(ip,2)
                   if(sw_gga_p == ON) wkchr4(j,i,k) = grad_tchgr_l(ip,1)
                end if
! lin
!             end if

             end if
          end do
       end do
    end do
    if(npes>=2)then
       call mpi_allreduce(MPI_IN_PLACE,wkchr2,nnp*nmp*nlp,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
       if(nspin==2) &
         & call mpi_allreduce(MPI_IN_PLACE,wkchr22,nnp*nmp*nlp,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
       if(sw_gga_p==ON) then
          call mpi_allreduce(MPI_IN_PLACE,wkchr4,nnp*nmp*nlp,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
       end if
    end if
    wkchr3 = 0.d0

    up_down=1

    if(nspin==1) then

    if(sw_epsilon_ele == OFF) then
       do i = 1, nmp
          do j = 1, nnp
             do k = 1, nlp
!!$             if(kimg == 1 .and. k > nlphf) then
!!$                knew = idp - k
!!$                jnew = nnp+2 - j
!!$                inew = nmp+2 - i
!!$                if(jnew > nnp) then
!!$                   jnew = jnew - nnp
!!$                end if
!!$                if(inew > nmp) then
!!$                   inew = inew - nmp
!!$                end if
!!$             else
!!$                knew = k; jnew = j; inew = i
!!$             end if
!!$             ip=nlphf*mmp*(jnew-1) + nlphf*(inew-1) + knew
!!$             wkchr(j,i,k) = afft(ip*2-2+up_down)
!!$             wkchr2(j,i,k) = tchgr_l(ip,1)
!c --------------->
!!$             if(sw_gga_p == ON) then
!!$                call enhance_gga_0(wkchr2(j,i,k)*nspin,wkchr4(j,i,k)*nspin,enhanc)
!!$             else
                call enhance_0(wkchr2(j,i,k)*nspin,enhanc)
!!$             end if
!c <--------------
                if(positron_filetype == DENSITY_ONLY .or. positron_filetype == VTK) then
                   wkchr3(k,i,j) = wkchr(k,i,j)*wkchr2(k,i,j)/(nnp*nmp*nlp)*univol*enhanc
                else
                   wkchr3(j,i,k) = wkchr(j,i,k)*wkchr2(j,i,k)/(nnp*nmp*nlp)*univol*enhanc
                end if
             end do
          end do
       end do
    else if(sw_epsilon_ele == ON) then
       do i = 1, nmp
          do j = 1, nnp
             do k = 1, nlp
                call enhance_01(wkchr2(j,i,k)*nspin,enhanc,epsilon_ele)
                if(positron_filetype == DENSITY_ONLY .or. positron_filetype == VTK) then
                   wkchr3(k,i,j) = wkchr(k,i,j)*wkchr2(k,i,j)/(nnp*nmp*nlp)*univol*enhanc
                else if(positron_filetype == CUBE) then
                   wkchr3(j,i,k) = wkchr(j,i,k)*wkchr2(j,i,k)/(nnp*nmp*nlp)*univol*enhanc
                end if
             end do
          end do
       end do
    end if

    else if(nspin==2) then
    if(sw_epsilon_ele == OFF) then
       do i = 1, nmp
          do j = 1, nnp
             do k = 1, nlp
                call enhance_0(wkchr2(j,i,k)+wkchr22(j,i,k),enhanc)
                if(positron_filetype == DENSITY_ONLY .or. positron_filetype == VTK) then
                   wkchr3(k,i,j) = wkchr(k,i,j)*wkchr2(k,i,j)/(nnp*nmp*nlp)*univol*enhanc
                   wkchr32(k,i,j) = wkchr(k,i,j)*wkchr22(k,i,j)/(nnp*nmp*nlp)*univol*enhanc
                else
                   wkchr3(j,i,k) = wkchr(j,i,k)*wkchr2(j,i,k)/(nnp*nmp*nlp)*univol*enhanc
                   wkchr32(j,i,k) = wkchr(j,i,k)*wkchr22(j,i,k)/(nnp*nmp*nlp)*univol*enhanc
                end if
             end do
          end do
       end do
    else if(sw_epsilon_ele == ON) then
       do i = 1, nmp
          do j = 1, nnp
             do k = 1, nlp
                call enhance_01(wkchr2(j,i,k)+wkchr22(j,i,k),enhanc,epsilon_ele)
                if(positron_filetype == DENSITY_ONLY .or. positron_filetype == VTK) then
                   wkchr3(k,i,j) = wkchr(k,i,j)*wkchr2(k,i,j)/(nnp*nmp*nlp)*univol*enhanc
                   wkchr32(k,i,j) = wkchr(k,i,j)*wkchr22(k,i,j)/(nnp*nmp*nlp)*univol*enhanc
                else if(positron_filetype == CUBE) then
                   wkchr3(j,i,k) = wkchr(j,i,k)*wkchr2(j,i,k)/(nnp*nmp*nlp)*univol*enhanc
                   wkchr32(j,i,k) = wkchr(j,i,k)*wkchr22(j,i,k)/(nnp*nmp*nlp)*univol*enhanc
                end if
             end do
          end do
       end do
    end if
    end if
!end of seperate the cases of nonmagnetic and magnetic

    end if
!!$      rewind 500
!!$      read(500,*) wkchr2
!!$      sss=0.d0
!!$      sss1=0.d0
!!$      sss2=0.d0
!!$      do j=1,nnp
!!$         do i=1,nmp
!!$            do k=1,nlp
!!$               call enhance_0(wkchr2(j,i,k),enhanc)
!!$               wkchr3(j,i,k)=wkchr(j,i,k)*wkchr2(j,i,k)/(nnp*nmp*nlp)*univol*enhanc
!!$               sss=sss+wkchr3(j,i,k)
!!$               sss1=sss1+wkchr(j,i,k)/(nnp*nmp*nlp)*univol
!!$               sss2=sss2+wkchr2(j,i,k)/(nnp*nmp*nlp)*univol
!!$            enddo
!!$         enddo
!!$      enddo
!!$      write(nfout,'(" (1) sss, sss1, sss2 = ",3f12.4)') sss, sss1, sss2

!!$      write(6,'(" tchgr_l")')
!!$    do i = ista_sfftph, ista_sfftph+10
!!$         write(6,'(i3,f20.8)') i, tchgr_l(i,1)
!!$    end do

   ! wkchr  <- afft
   ! wkchr2 <- tchgr_l
      sss=0.d0
      sss1=0.d0
      sss2=0.d0
      sss_down = 0.d0
      sss_total = 0.d0

    if(nspin == 1) then
      if(sw_epsilon_ele == OFF) then
         do i = ista_sfftph, iend_sfftph
!c --------------->
!!$         if(sw_gga_p == ON) then
!!$            call enhance_gga_0(tchgr_l(i,1)*nspin,grad_tchgr_l(i,1)*nspin,enhanc)
!!$         else
            call enhance_0(tchgr_l(i,1)*nspin,enhanc)
!!$         end if
!c <--------------
            sss =  sss + afft(i*2-1)*tchgr_l(i,1)*nspin*f2or1(i)*enhanc
            sss1 = sss1 + afft(i*2-1)*f2or1(i)
            sss2 = sss2 + tchgr_l(i,1)*nspin*f2or1(i)
         end do
      else if(sw_epsilon_ele == ON) then
         do i = ista_sfftph, iend_sfftph
            call enhance_01(tchgr_l(i,1)*nspin,enhanc,epsilon_ele)
            sss =  sss + afft(i*2-1)*tchgr_l(i,1)*nspin*f2or1(i)*enhanc
            sss1 = sss1 + afft(i*2-1)*f2or1(i)
            sss2 = sss2 + tchgr_l(i,1)*nspin*f2or1(i)
         end do
      end if
!magnetic system,nspin=2
    else if(nspin==2) then
       if(sw_epsilon_ele == OFF) then
         do i = ista_sfftph, iend_sfftph
!c --------------->
!!$         if(sw_gga_p == ON) then
!!$            call enhance_gga_0(tchgr_l(i,1)*nspin,grad_tchgr_l(i,1)*nspin,enhanc)
!!$         else
            call enhance_0((tchgr_l(i,1)+tchgr_l(i,2)),enhanc)
!!$         end if
!c <--------------
            sss =  sss + afft(i*2-1)*tchgr_l(i,1)*nspin*f2or1(i)*enhanc     ! up_spin
            sss_down =  sss_down + afft(i*2-1)*tchgr_l(i,2)*nspin*f2or1(i)*enhanc    ! down_spin
            sss_total =  sss_total + afft(i*2-1)*(tchgr_l(i,1)+tchgr_l(i,2))*f2or1(i)*enhanc      ! total_spin
            sss1 = sss1 + afft(i*2-1)*f2or1(i)
            sss2 = sss2 + (tchgr_l(i,1)+tchgr_l(i,2))*f2or1(i)
         end do
      else if(sw_epsilon_ele == ON) then
         do i = ista_sfftph, iend_sfftph
            call enhance_01((tchgr_l(i,1)+tchgr_l(i,2)),enhanc,epsilon_ele)
            sss =  sss + afft(i*2-1)*tchgr_l(i,1)*nspin*f2or1(i)*enhanc     ! up_spin
            sss_down =  sss_down + afft(i*2-1)*tchgr_l(i,2)*nspin*f2or1(i)*enhanc    ! down_spin
            sss_total =  sss_total + afft(i*2-1)*(tchgr_l(i,1)+tchgr_l(i,2))*f2or1(i)*enhanc      ! total_spin
            sss1 = sss1 + afft(i*2-1)*f2or1(i)
            sss2 = sss2 + (tchgr_l(i,1)+tchgr_l(i,2))*f2or1(i)
         end do
      end if
!end of magnetic system
    end if

      if(npes >= 2) then
         call mpi_allreduce(MPI_IN_PLACE,sss,1,mpi_double_precision,mpi_sum, MPI_CommGroup,ierr)
         call mpi_allreduce(MPI_IN_PLACE,sss1,1,mpi_double_precision,mpi_sum, MPI_CommGroup,ierr)
         call mpi_allreduce(MPI_IN_PLACE,sss2,1,mpi_double_precision,mpi_sum, MPI_CommGroup,ierr)
         if(nspin==2) then
            call mpi_allreduce(MPI_IN_PLACE,sss_down,1,mpi_double_precision,mpi_sum, MPI_CommGroup,ierr)
            call mpi_allreduce(MPI_IN_PLACE,sss_total,1,mpi_double_precision,mpi_sum, MPI_CommGroup,ierr)
         end if
      end if
      if(nspin == 1) then
         sss_down = sss
         sss_total = sss
      end if

      sss  = sss/(nnp*nmp*nlp)*univol
      sss1 = sss1/(nnp*nmp*nlp)*univol
      sss2 = sss2/(nnp*nmp*nlp)*univol
      sss_down  = sss_down/(nnp*nmp*nlp)*univol
      sss_total  = sss_total/(nnp*nmp*nlp)*univol

      if(ipripositron >= 1) write(nfout,'(" (2) sss, sss1, sss2 = ",&
           & 3f12.4)') sss, sss1, sss2

      if(sw_positron_file == ON) then
         call m_Files_open_nfpstrn()
         call m_Files_open_nfvelec()
         call m_Files_open_nfeppair()
         if(nspin==2) call m_Files_open_nfeppair2()
         if(sw_gga_p == ON) call m_Files_open_nfvelec_grad()

!!$rewind 300
!!$rewind 310
!!$rewind 320
         if(mype == 0) call wdchgr(nfout,nfpstrn,1,positron_title(1))
         if(mype == 0) call wdchgr(nfout,nfvelec,2,positron_title(2))
         if(mype == 0) call wdchgr(nfout,nfeppair,3,positron_title(3))
         if(nspin==2) then
            if(mype == 0) call wdchgr(nfout,nfeppair2,5,positron_title(5))
         end if
!!$write(300,*) 'cube file of positron density'
!!$write(300,'(6e13.5)') wkchr
!!$write(310,*) 'cube file of valence electron density'
!!$write(310,'(6e13.5)') wkchr2
!!$write(320,*) 'cube file of e-p pair'
!!$write(320,'(6e13.5)') wkchr3
         if(sw_gga_p == ON) then
!!$         rewind 330
            if(mype == 0) call wdchgr(nfout,nfvelec_grad,4,positron_title(4))
!!$         write(330,*)  'cube file of gradient of valence electron density'
!!$         write(330,'(6e13.5)') wkchr4
         end if

         call m_Files_close_nfpstrn()
         call m_Files_close_nfvelec()
         call m_Files_close_nfeppair()
         if(nspin==2) call m_Files_close_nfeppair2()
         if(sw_gga_p == ON) call m_Files_close_nfvelec_grad()
         if(sw_gga_p == ON) deallocate(wkchr4)
         deallocate(wkchr)
         deallocate(wkchr2)
         deallocate(wkchr3)
         if(nspin==2) deallocate(wkchr22)
         if(nspin==2) deallocate(wkchr32)
      end if

     if(nspin==1) then
       valence_annihilation_rate=sss/2.d0

#ifdef SINGLE_POSITRON
       valence_annihilation_rate = sss
#endif

       ssk=((BOHR)**3)*1.d4/(2.8d0**2)/3.d0/PAI
       ssk=ssk
       sss=core_annihilation_rate+valence_annihilation_rate
       p_new_lifetime=ssk/sss
       p_core_rate=core_annihilation_rate/sss*100.d0
     else if(nspin==2) then
        valence_annihilation_rate=sss/2.d0
        valence_annihilation_rate_down=sss_down/2.d0
        valence_annihilation_rate_total=sss_total/2.d0
        ssk=((BOHR)**3)*1.d4/(2.8d0**2)/3.d0/PAI
        ssk=ssk
        sss=core_annihilation_rate+valence_annihilation_rate
        sss_down=core_annihilation_rate+valence_annihilation_rate_down
        sss_total=core_annihilation_rate+valence_annihilation_rate_total
        p_new_lifetime=ssk/sss
        p_new_lifetime_down=ssk/sss_down
        p_new_lifetime_total=ssk/sss_total
        p_core_rate=core_annihilation_rate/sss*100.d0
        p_core_rate_down=core_annihilation_rate/sss_down*100.d0
     end if
!!$     if(ipripositron >= 1) write(nfout,*)core_annihilation_rate, valence_annihilation_rate

     deallocate(f2or1)
   contains
     subroutine wdchgr(nfout,nfchr,np,str)
       integer, intent(in) :: nfout, nfchr,np
       character(len=*), intent(in) :: str
       integer :: n1,n2,n3,i,j,k, m, iloop
       real(kind=DP) :: dn1,dn2,dn3,x,y,z
       real(kind=DP),allocatable,dimension(:,:) :: cps_full
       integer, allocatable,dimension(:) :: ityp_full

       if(positron_filetype == DENSITY_ONLY) then
          write(nfchr,9001) nlp*nmp*nnp, nlp, nmp, nnp
9001      format(' CHARGE DENSITY NE = ',i8,'(',3i5,')')
          if(np == 1) then
             write(nfchr,'(6e13.5)') wkchr
          else if(np == 2) then
             write(nfchr,'(6e13.5)') wkchr2
          else if(np == 3) then
             write(nfchr,'(6e13.5)') wkchr3
          else if(np == 4) then
             write(nfchr,'(6e13.5)') wkchr4
          else if(np == 5) then
             if(allocated(wkchr32)) &
                  & write(nfchr,'(6e13.5)') wkchr32
          end if
       else if(positron_filetype == VTK) then
          write(nfchr,'("# vtk DataFile Version 2.0")')
          if(len_trim(str) >= 1) then
             write(nfchr,*) trim(str)
          else
             write(nfchr,'(" Calculated by PHASE")')
          end if

          write(nfchr,'("ASCII")')
          write(nfchr,'("DATASET STRUCTURED_GRID")')
          write(nfchr,'("DIMENSIONS",3(1x,i5))') nlp+1,nmp+1,nnp+1
          write(nfchr,'("POINTS",1x,i7,1x,"float")') (nlp+1)*(nmp+1)*(nnp+1)
          do n1=0,nlp
             do n2=0,nmp
                do n3=0,nnp
                   dn1 = n1/dble(nlp)
                   dn2 = n2/dble(nmp)
                   dn3 = n3/dble(nnp)
                   x = altv(1,1)*dn1 + altv(1,2)*dn2 + altv(1,3)*dn3
                   y = altv(2,1)*dn1 + altv(2,2)*dn2 + altv(2,3)*dn3
                   z = altv(3,1)*dn1 + altv(3,2)*dn2 + altv(3,3)*dn3
                   write(nfchr,'(3(1x,e13.5))') x,y,z
                end do
             end do
          end do
          if(np.le.4.or.(np==5.and.allocated(wkchr32))) then
             write(nfchr,'("")')
             write(nfchr,'("POINT_DATA",1x,i7)') (nlp+1)*(nmp+1)*(nnp+1)
             write(nfchr,'("SCALARS scalars float")')
             write(nfchr,'("LOOKUP_TABLE default")')
             do n1=0,nlp
                i=n1+1
                if(n1==nlp) i=1
                do n2=0,nmp
                   j=n2+1
                   if(n2==nmp) j=1
                   do n3=0,nnp
                      k=n3+1
                      if(n3==nnp) k=1
                      if(np==1) then
                         write(nfchr,'(e13.5)') wkchr(i,j,k)
                      else if(np==2) then
                         write(nfchr,'(e13.5)') wkchr2(i,j,k)
                      else if(np==3) then
                         write(nfchr,'(e13.5)') wkchr3(i,j,k)
                      else if(np==4) then
                         write(nfchr,'(e13.5)') wkchr4(i,j,k)
                      else if(np==5) then
                         write(nfchr,'(e13.5)') wkchr32(i,j,k)
                      end if
                   end do
                end do
             end do
          end if
       else if(positron_filetype == CUBE) then
          if(len_trim(str) >= 1) then
             write(nfchr,*) trim(str)
          else
             write(nfchr,'(" Calculated by phase")')
          end if
          if(nspin == 2) then
             if(iloop == 1) then
                write(nfchr,'(" SCF Total Density  UP")')
             else
                write(nfchr,'(" SCF Total Density  DOWN")')
             end if
          else
             write(nfchr,'(" SCF Total Density")')
          end if
          x = 0.d0; y = 0.d0; z = 0.d0
          write(nfchr,'(i6,3f10.4)') natm2, x,y,z
          do i = 1, 3
             write(nfchr,'(i6,3f10.6)') fft_box_size_CD(i,1), altv(1:3,i)/dble(fft_box_size_CD(i,1))
          end do

          allocate(cps_full(natm2,3))
          allocate(ityp_full(natm2))
          call m_IS_pack_all_ions_in_uc(ityp_full,cps_full)
          do i = 1, natm2
             m = ityp_full(i)
             write(nfchr,'(i6,4f10.6)') nint(iatomn(m)), ival(m), cps_full(i,1:3)
          end do
          deallocate(ityp_full,cps_full)

          if(np==1) then
             write(nfchr,'(6e13.5)') wkchr
          else if(np==2) then
             write(nfchr,'(6e13.5)') wkchr2
          else if(np==3) then
             write(nfchr,'(6e13.5)') wkchr3
          else if(np==4) then
             write(nfchr,'(6e13.5)') wkchr4
          else if(np==5) then
             if(allocated(wkchr32)) &
                  & write(nfchr,'(6e13.5)') wkchr32
          end if
       end if
     end subroutine wdchgr
end subroutine m_pWF_valence_annihilation

subroutine m_pWF_wlifetime()

  if(ipripositron >= 1) then
  if(nspin==1)then
     write(nfout,*) '***************************************************'
     write(nfout,*)'positron lifetime(ps)', p_new_lifetime
     write(nfout,*)'core rate',p_core_rate,'%'
     write(nfout,*)'core and valence',core_annihilation_rate, valence_annihilation_rate
     write(nfout,*) '*************************************************'
  elseif(nspin==2)then
     write(nfout,*) '******************** Majority &      Minority ****'
     write(nfout,*)'positron lifetime(ps)', p_new_lifetime,p_new_lifetime_down
     write(nfout,*)'core rate            ',p_core_rate,p_core_rate_down,'%'
     write(nfout,*)'core and valence(up,down)',core_annihilation_rate, valence_annihilation_rate,valence_annihilation_rate_down
     write(nfout,*) '*************************************************'
  end if
  end if
end subroutine m_pWF_wlifetime

  subroutine m_pWF_update_lifetime
    if(ipripositron >= 1) write(nfout,*)'lifetime: ',p_old_lifetime,p_new_lifetime
    p_old_lifetime=p_new_lifetime
  end subroutine m_pWF_update_lifetime

  subroutine m_pWF_wd_pzaj(nfout,comment,nc)
    integer,        intent(in) :: nfout, nc
    character(len=nc), intent(in) :: comment

    character(len=5) :: a
    integer :: i, ib, ri, j
    integer, parameter :: NZAJSIZE = 20
    integer :: nelm

    if(ipripositron <= 1) return
    write(nfout,*) comment
    a = "     "
    do ib = 1, npeg
       do ri = 1, kimg
          if(ri == 1 .and. kimg == 2) a = "(Re) "
          if(ri == 2) a = "(Im) "

          if(ri == 1) write(nfout,'(" pev(",i3,")= ",e14.6," ",a4,5e14.6)') &
               & ib,pev(ib),a,(pzaj(i,ib,ri),i=1,5)
          if(ri == 2) write(nfout,'(31x,a4,5e14.6)') a,(pzaj(i,ib,ri),i=1,5)
       end do
       if(ipripositron >= 3) then
          nelm = min(kg1_pwf,NZAJSIZE)
          do ri = 1, kimg
             if(ri == 1 .and. kimg == 2) a = "(Re) "
             if(ri == 2) a = "(Im) "

             if(ri == 1) write(nfout,'(" pev(",i3,")= ",e14.6," ",a4,5e14.6)') &
                  & ib,pev(ib),a,(pzaj(i,ib,ri),i=1,5)
             if(ri == 2) write(nfout,'(31x,a4,5e14.6)') a,(pzaj(i,ib,ri),i=1,5)
             write(nfout,'(35x,5e14.6)') (pzaj(i,ib,ri),i=6,nelm)

          end do
       end if
    end do
  end subroutine m_pWF_wd_pzaj
end module m_Positron_Wave_Functions
