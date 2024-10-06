!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 552 $)
!
!  "First-principles Electronic Structure Calculation Program"
!
!  PROGRAM: vdW-Soler
!
!  AUTHOR(S): Y. Ono
!
!  Contact address :  Phase System Consortium
!                     E-mail: phase_system@nims.go.jp URL https://azuma.nims.go.jp
!
!
!
!=======================================================================
!
!   The original version of this set of the computer programs "PHASE" was developed by
!  the members of the Theory Group of Joint Research Center for Atom Technology
!  (JRCAT), based in Tsukuba, in the period 1993-2001.
!   Since 2002, this program set had been intensively developed as a part of the following
!  national projects supported by the Ministry of Education, Culture, Sports, Science and
!  Technology (MEXT) of Japan; "Frontier Simulation Software for Industrial Science
!  (FSIS)" from 2002 to 2005, "Revolutionary Simulation Software (RSS21)" from 2006 to
!  2008. "Research and Development of Innovative Simulation Software (RISS)" from 2008
!  to 2013. These projects is lead by the Center for Research on Innovative Simulation
!  Software (CISS), the Institute of Industrial Science (IIS), the University of Tokyo.
!   Since 2013, this program set has been further developed centering on PHASE System
!  Consortium.
!   The activity of development of this program set has been supervised by Takahisa Ohno.
!
!
!********************************** Note ******************************************
! This program calculates the non-local correlation energy (Ecnl) and
! the local correlation energy (EcLDA) as the post calculation by utilizing
! the output files from the PHASE. By adding these correlation terms as
!    E_total = E_total(GGAx) + EcLDA + Ecnl,
! the van der Waals interaction will be included into the total energy.
!
! This program follows the Dion et al's 1-shot method (so called vdW-DF), and
! the Roman-Perez et al's convolution algorithm. Because of this convolution,
! the CPU cost in this program is reduced to O(NlogN) from the original O(N**2).
!                                 (Reference; M.Dion et al. PRL 92 (2004) 246401.)
!                         (Reference; G.Roman-Perez et al. PRL 103 (2009) 096102.)
!
!
! Periodic systems are assumed.
! The atomic units (Hartree) are used.
!
! ======= modification =====
!
!  2016/06/06 : by asms
!         The FFT normalization parameters are changed in order to be consistent
!         with the other subroutines.
!
!         rho(G) = 1/V *Int rho(r) exp(-iGr);   rho(r) = sum_{G} rho(G) exp(iGr)
!
!         In the dicretized grids, rho(G) = FFT[rho(r)] /(na*nb*nc);
!                                  rho(r) = FFT[rho(G)]
!
!         By this, the nonlocal vdW energy is written as
!             Encl = V**2 /2 sum_{ij} sum_{G} conjg(theta(i,G))*phi(i,j,G)*theta(i,G)
!         and
!             phi(i,j,G) = 1/V Int phi(i,j,r) exp(-iGr)
!
! =========================
!
! ++++++++ List of subroutines ++++++
! All subroutine files listed below are included in this file.
! In addition to the below list, it is necessary to link to a FFT library, for example FFTW3.
!
! CPU_TIME             : Check the calculation time.
! derivation           : Calculate the derivations of the electron density distribution.
! d_q0                 : Obtain q(rho(r)) from the electron density rho(r) and its derivations.
! spline               : Determine the function 'p(q)' for bi-cubic spline interpolation.
!                        Here p(r) is defined as, phi(q1,q2,r12) = SUM_ab(p(qa)*p(qb)*phi(qa,qb,r12)).
! theta_ab             : Multiply theta(r) = p(r)*rho(r) by use of the given parameter qa or qb.
! RtoG                 : Calculate the Fourier transform of theta(r) within the FFT algorithm.
! phiab                : Calculate the Fourier transform of the kernel phi(r) by use of the given
!                        parameter qa and qb. phi(r) does not depends on angular components,
!                        so that it can be done within the radial mesh.
! convolution_3d       : Calculate the integral SUM(dk*theta_a(k)*theta_b(k)*phi(k)).
! piDphi               : Calculate the local part of dr**2*rho(r1)*rho(r2)*phi(r1,r2,r12) directory.
!                        The local part includes the singular point, thus it is necessary to execute
!                        the 3D space integral directory.
! cLDA                 : Obtain the LDA version of correlation term EcLDA.
! outputs              : Output the results.
! kernel_phi           : Calculate the value of kernel phi(d1,d2).
! gauleg               : Prepare the grid points for the Gauss-Legendre integral.
!
!
!
! ++++++++ Input & Output +++++++++++
! Input files  (Both of the input files will be output by PHASE.)
!   nfchr.cube        : An electron density distribution file given by
!                       GGA(exchange only) calculation
!   nfefn.data        : The total energy given by GGA(exchange only) calculation
!
! Output   E_total(vdW-DF) = E_total(GGAx) + EcLDA + Ecnl
!   E_total(GGAx)     : Total energy given by GGA(exchange only) calculation
!   Ecnl              : Non-local Correlation energy
!   EcLDA             : LDA Correlation energy
!
!
!
! ++++++ Internal parameters ++++++++
! dq      (       )   : The smallest distance between the grid points for bi-cubic
!                       spline interpolation. The grid points are plotted in logarithmic manner.
! lambda   (       )   : The power of dq
! nr12    (       )   : The number of the grid points for the kernel phi(r12)
! phi0    (       )   : One of the coefficients in the 'soft' function,
!                       phi_s(r12) = phi0 + phi2*d**2 + phi4*d**4. phi2 and phi4 will be
!                       eventually determined to match the value and the slope of phi_s and
!                       those of phi at d=d_s, in the program.
! ds      (       )   : The connection point of phi_s and phi
! r12max  (       )   : The cutoff for r12
! rhomin  (       )   : Minimum value of electron density.
!                       If the electron density was smaller than rhomin, then the program
!                       will read as this minimum value in order to avoid a divergence.
! nk      (       )   : The number of the grid points plotted for the Fourier
!                       transformed kernel phi(k)
! maxk    (       )   : The cutoff for k
!
!
!                                              Written by Youky Ono in June/2013
!**********************************************************************************
#ifndef DISABLE_VDWDF
module progress_bar
  implicit none
  logical,private :: printable = .false.
  integer,private :: iend=10
  integer,private :: uni=6
  integer,private :: j=0
  logical, dimension(0:9),private :: done

contains

  subroutine set_printable(pri)
    logical, intent(in) :: pri
    printable = pri
  end subroutine set_printable

  subroutine reset_progress()
    j = 0
    done = .false.
    if(printable) write(unit=uni,fmt="(a10)") "0%    100%"
  end subroutine reset_progress

  subroutine set_unit(un)
    integer, intent(in) :: un
    uni = un
  end subroutine set_unit

  subroutine set_end(en)
    integer, intent(in) :: en
    iend = en
  end subroutine set_end

  subroutine progress()
    use m_Const_Parameters, only : DP
    implicit none
    integer(kind=4)::k
    character(len=17)::bar="???% |          |"
    integer :: jj
    real(kind=DP) :: jjj
    j = j+1
    if(j.gt.iend) j=iend
    jj = int(10*(dble(j)/dble(iend))*0.999999d0)
    if(.not.done(jj).and.j<=iend) then
       if(printable) write(unit=uni,fmt="(a1,$)") '*'
       call flush(uni)
       done(jj) = .true.
    endif
    if(j==iend.and.printable) write(unit=uni,fmt=*)
    !     jjj = (dble(j)/dble(iend))
    !     write(unit=bar(1:3),fmt="(i3)") int(100*jjj)
    !     do k=1, jj
    !       bar(6+k:6+k)="*"
    !     enddo
    !     ! print the progress bar.
    !     write(unit=uni,fmt="(a1,a17,$)") char(13), bar
    !!     write(unit=uni,fmt="(a1,a17)") char(13), bar
    !     if (j/=iend) then
    !       flush(unit=6)
    !     else
    !       write(unit=6,fmt=*)
    !       do k=1,jj
    !          bar(6+k:6+k) = ""
    !       enddo
    !     endif
    return
  end subroutine progress

end module progress_bar


module m_vdWDF
  use m_Const_Parameters, only : DP, ON, PAI, FMAXVALLEN, LOWER, CMPLDP, DIRECT, INVERSE
  use m_Control_Parameters, only : printable,nspin,eval_kernel_by_interpolation, &
       &                           na_gl,a1,a2,dq_vdw,lambda,q0cut,ds,ndel,nphiD, &
       &                           nr12,nk,maxk, r12max,oneshot, ipri, &
       &                           sw_save_memory_vdw, kimg, iprixc, &
       &                           q0min, sw_use_WuGygi_method, sw_fft_xzy
  use m_Files, only : nfout
  use m_Charge_Density, only : m_CD_get_rspace_charge

  use m_FFT, only : fft_box_size_CD,fft_box_size_CD_nonpara, m_FFT_CD0,nfftp_nonpara, fft_box_size_CD_3D
  use m_Crystal_Structure, only : altv,univol, rltv

  use m_Parallelization, only : npes,mype,MPI_CommGroup,m_Parallel_init_mpi_nq,  &
       &                        ista_nq,iend_nq,np_nq,mp_nq,is_nq,ie_nq,nel_nq, &
       &                        map_z_nq, np_fftcd_y, np_fftcd_z, mp_fftcd_y, mp_fftcd_z, &
       &                        mpi_g_world, nrank_ke,myrank_ke, nel_fftcd_x, nel_fftcd_y, nel_fftcd_z, &
       &                        mpi_ke_world, nrank_g, myrank_g
  use m_Timing, only : tstatc0_begin,tstatc0_end
  use mpi

  implicit none
!  include 'mpif.h'

  real(kind=DP),parameter :: pi=PAI

  ! Physical values
  Real(kind=DP) ExGGA,Ecnl,Ecnl_12,Ecnl_12_ab,Ecnl_3,Ecnl_3s,EcLDA

  real(kind=DP) :: phi0=2.77d0

  Integer na,nb,nc,nabc
  Real(kind=DP) :: aa(3,3),dv
  Real(kind=DP), Allocatable :: rho(:),grad(:), cgrad(:,:)
  real(kind=DP),allocatable,dimension(:) :: phi_ab

  ! Grid points
!!!! Spline curves
  Integer nq0
  real(kind=DP) :: qa,qb,q0max

!!!! The table of phidD
  real(kind=DP), allocatable, dimension(:,:) :: phidD

  ! Internal parameters
  Real(kind=DP),parameter :: rhomin=1.d-9

  ! Real-space and reciprocal-space Functions

  Complex*16, Allocatable :: theta_G(:),     &
       &                            theta_G_ab(:,:),theta_G_a(:),theta_G_b(:)
  Real(kind=DP), Allocatable :: theta_R(:),dtheta_R(:,:),ddtheta_R(:,:)

  complex(kind=DP), allocatable, dimension(:,:) :: ualpha_g1,ualpha_r1

  real(kind=DP) :: etai = 1.3d0,eta1 = 8.d0

  real(kind=DP), allocatable, dimension(:)   :: qar
  real(kind=DP), allocatable, dimension(:,:) :: q2ar

  real(kind=DP) :: rinplw
  logical :: firstcall = .true.

  integer, parameter :: FFTW_ESTIMATE = 64
  integer, parameter :: FFTW_FORWARD  = -1
  integer, parameter :: FFTW_BACKWARD = +1

  real(kind=DP), allocatable, dimension(:) :: dFdrho, dFddrho
  real(kind=DP), allocatable, dimension(:) :: rkar

  real(kind=DP) :: s_cnl1(3,3), s_cnl2(3,3), ecnl_vdwdf
  real(kind=DP) :: s_cnl1_pc(3,3), s_cnl2_pc(3,3)

  !  logical :: grad_rho_eq_0 = .true.
  logical :: grad_rho_eq_0 = .false.

! ---- extension to vdWDF2 ----
  integer :: version_no_vdwdf = 1

  real(kind=DP) :: Zab
  real(kind=DP), parameter :: Zab_0 = -0.8491d0
! ---------------------------

  integer :: nfft_div_size
  integer, allocatable, dimension(:) :: mp_fftcd, imp_fftcd

  integer, allocatable, dimension(:) :: ica,icb,icc

  logical :: lprog = .false.
contains

  integer function get_index_1D_kimg(cix,ciy,ciz,ignore_kimg)
    integer, intent(in) :: cix,ciy,ciz
    logical, intent(in), optional :: ignore_kimg
    integer :: cix2,ciy2,ciz2
    integer :: idp,mmp,nlphf
    logical :: ig
    ig = .false.
    if(present(ignore_kimg)) ig = ignore_kimg

    idp = fft_box_size_CD_nonpara(1,0)
    mmp = fft_box_size_CD_nonpara(2,0)
    if(kimg == 1 .and. .not.ig) then
       nlphf = idp/2
    else
       nlphf = idp
    end if
    if ( cix > nlphf ) then
       cix2 = idp -cix
       ciy2 = nb +2 -ciy
       ciz2 = nc +2 -ciz
       if ( ciy2 > nb ) ciy2 = ciy2 -nb
       if ( ciz2 > nc ) ciz2 = ciz2 -nc
    else
       cix2 = cix;  ciy2 = ciy; ciz2 = ciz
    endif
    get_index_1D_kimg = (ciz2-1)*mmp*nlphf +(ciy2-1)*nlphf +cix2
  end function get_index_1D_kimg

  integer function get_index_1D(cix,ciy,ciz)
    integer, intent(in) :: cix,ciy,ciz
    integer :: cix2,ciy2,ciz2
    integer :: idp,mmp,nlphf
    get_index_1D = (ciz-1)*nb*na + (ciy-1)*na + cix
  end function get_index_1D

  function get_index_3D(cir) result(ret)
    integer, intent(in) :: cir
    integer :: ret(3)
    integer :: cix,ciy,ciz
    cix = 1+(cir-1-MOD(cir-1+nb*nc,nb*nc))/(nb*nc)
    ciy = 1+(cir-nb*nc*(cix-1)-1-MOD(cir-1+nc,nc))/nc
    ciz = cir-nc*(nb*(cix-1)-1+ciy)
    ret(1) = cix
    ret(2) = ciy
    ret(3) = ciz
  end function get_index_3D

  subroutine print_vdw_parameters()
    if(printable)then
       write(nfout,'(a)')            '-- parameters for the vdW-DF calculations --'
       write(nfout,'(a,2f15.10,i8)') '   rmax, kmax, nmesh     : ',r12max,maxk,nr12
       write(nfout,'(a,3f15.10)')    '   dq, lambda, and q0cut : ',dq_vdw,lambda,q0cut
       write(nfout,'(a,2f15.10,i5)') '   q0min, q0max, nq0     : ',q0min,q0max,nq0
       write(nfout,'(a,f15.10)' )    '   ds                    : ',ds
       write(nfout,'(a,i5,2f15.10)') '   na_gl, a1, a2         : ',na_gl,a1,a2
       if(eval_kernel_by_interpolation)then
          write(nfout,'(a)')            '   kernel evaluation     : by interpolation'
          write(nfout,'(a,2i8)')        '   ndel, nphiD           : ',ndel,nphiD
       else
          write(nfout,'(a)')            '   kernel evaluation     : direct'
       endif
       write(nfout,'(A,I3)') '   sw_use_WuGygi_method : ', sw_use_WuGygi_method
    endif
  end subroutine print_vdw_parameters

  subroutine initialize_vdwdf_scf(nspin,ispin,nfft_div_size, na, nb, nc, chgr,grad_rho,version_no)
    use progress_bar, only : set_printable

    integer, intent(in) :: nspin,ispin,nfft_div_size,na,nb,nc, version_no
    real(kind=DP), dimension(nfft_div_size), intent(in) :: chgr,grad_rho

    integer :: i,cix,ciy,ciz,nrxyz,cir,cir2
    integer :: idp, mmp, nlphf, cix2, ciy2, ciz2
    real(kind=DP) :: q

    lprog = .false.
    if ( version_no == 1 ) then
       Zab = Zab_0
    else if ( version_no == 2 ) then
       Zab = Zab_0 *2.222d0
    endif

    call set_printable(printable)
    call do_cell_params()

    rinplw = 1.d0/(dble(na*nb*nc))
    q0max = q0cut*1.01d0
!!!    q0min = 0.09d0
    nq0 = DINT(dLOG((q0max-q0min)*(lambda-1.d0)/dq_vdw+1.d0)/dLOG(lambda))+1
    maxk = dble(nr12)/r12max

    if(firstcall) call m_Parallel_init_mpi_nq(nfout,ipri,printable,nq0)
    if(firstcall) call print_vdw_parameters()
    call alloc_vdw()

    rho(:) = chgr(:)
    grad(:) = grad_rho(:)
!    do cix = 1,na
!       do ciy = 1,nb
!          do ciz = 1,nc
!             cir = get_index_1D_kimg(cix,ciy,ciz)
!             cir2 = get_index_1D(cix,ciy,ciz)
!             rho(cir2) = chgr(cir)
!             grad(cir2) = grad_rho(cir)
!          enddo
!       enddo
!    enddo

    do i=1,nq0
       q = q0min + dq_vdw*(lambda**DBLE(i-1)-1.d0)/(lambda-1.d0)
       qar(i) = q
    enddo
    call spline0(nq0,qar,q2ar)

    if (eval_kernel_by_interpolation.and.firstcall) then
       call build_lookup_table(ndel,nphiD,phidD)
    endif
    firstcall = .false.
  end subroutine initialize_vdwdf_scf

  subroutine initialize_vdwdf_oneshot(is)
    use progress_bar, only : set_printable
    integer, intent(in) :: is
    integer :: i,cir,cix,ciy,ciz,nrxyz
    real(kind=DP) :: q
    lprog = .true.
    call set_printable(printable)
    if(printable) write(nfout,'(a)') 'initialization ...'

    call do_cell_params()

    rinplw = 1.d0/(dble(na*nb*nc))
    q0max = q0cut*1.01d0
!!!    q0min = 0.09d0
    nq0 = DINT(dLOG((q0max-q0min)*(lambda-1.d0)/dq_vdw+1.d0)/dLOG(lambda))+1
    maxk = dble(nr12)/r12max
    if(firstcall) call m_Parallel_init_mpi_nq(nfout,ipri,printable,nq0)

    call alloc_vdw()
    call m_CD_get_rspace_charge(nfout,na,nb,nc,rho,is)
    call print_vdw_parameters()
    call derivation(na,nb,nc,aa,rho,dv,grad)

    do i=1,nq0
       q = q0min + dq_vdw*(lambda**DBLE(i-1)-1.d0)/(lambda-1.d0)
       qar(i) = q
    enddo

    call spline0(nq0,qar,q2ar)

    if(eval_kernel_by_interpolation)then
       if(printable) write(nfout,'(a)') 'building the lookup table for the kernel function ...'
       call build_lookup_table(ndel,nphiD,phidD)
       if(printable) write(nfout,'(a)') '... done'
    endif
    if(printable) write(nfout,'(a)') '... done initialization'

  end subroutine initialize_vdwdf_oneshot

  subroutine do_cell_params()
    integer :: i
    na = fft_box_size_CD_3D(1,1)
    nb = fft_box_size_CD_3D(2,1)
    nc = fft_box_size_CD_3D(3,1)
    nabc = na * nb * nc
    do i=1,3
       aa(i,1:3) = altv(1:3,i)/dble(fft_box_size_CD(i,1))
    enddo
    dv = aa(1,1)*(aa(2,2)*aa(3,3)-aa(2,3)*aa(3,2)) &
         &      + aa(1,2)*(aa(2,3)*aa(3,1)-aa(2,1)*aa(3,3)) &
         &      + aa(1,3)*(aa(2,1)*aa(3,2)-aa(2,2)*aa(3,1))
  end subroutine do_cell_params

  function real_index(iq,pe) result(res)
    integer, intent(in) :: iq
    integer, intent(in),optional :: pe
    integer :: res
    integer :: ii,mpe
    mpe = myrank_ke
    if(present(pe))then
       mpe = pe
    endif
    if(mpe==0)then
       res = iq
       return
    endif
    res = 0
    do ii=0,mpe-1
       res = res+nel_nq(ii)
    enddo
    res = res+iq
    return
  end function real_index

  subroutine alloc_vdw()
    integer :: i,cix,ciy,ciz,itmp,c1,c2,c3,ca,cb,cc
    integer :: imin, imax, maxcount
    real(kind=DP) :: Ta,Tb,Tc, vec(3), bb(3,3)

    allocate(imp_fftcd(nabc));imp_fftcd=0
#ifdef FFT_3D_DIVISION_CD
    nfft_div_size = np_fftcd_x
    allocate(mp_fftcd(nfft_div_size));mp_fftcd=0
    mp_fftcd(:) = mp_fftcd_x(:)
#else
    if (sw_fft_xzy > 0) then
       if (kimg == 1) then
          nfft_div_size = np_fftcd_y / 2
       else
          nfft_div_size = np_fftcd_y
       end if
       allocate(mp_fftcd(nfft_div_size));mp_fftcd=0
       mp_fftcd(:) = mp_fftcd_y(:)
    else
       if (kimg == 1) then
          nfft_div_size = np_fftcd_z / 2
       else
          nfft_div_size = np_fftcd_z
       end if
       allocate(mp_fftcd(nfft_div_size));mp_fftcd=0
       mp_fftcd(:) = mp_fftcd_z(:)
    end if
!!$ T.Yamasaki 2019.06.12
!!$    write(nfout,'(" ** alloc_vdw ** nabc = ",i10)') nabc
!!$    write(nfout,'(" **   sw_fft_xzy = ",i8)') sw_fft_xzy
!!$    write(nfout,'(" **   nfft_div_size = ",i5)') nfft_div_size
    imin = 10000
    imax = -1
    maxcount = 0
    do i = 1, nfft_div_size
       if(imax < mp_fftcd(i)) imax = mp_fftcd(i)
       if(imin > mp_fftcd(i)) imin = mp_fftcd(i)
       if(mp_fftcd(i) > nabc) maxcount = maxcount+1
    end do
!!$    write(nfout,'(" **   mp_fftcd = [",i8," : ",i8,"], counter of exceeded = ",i8)') imin, imax,maxcount
!!$    call flush(nfout)

    if(maxcount > 0) then
       do i=1,nfft_div_size
          if(mp_fftcd(i) > nabc) cycle
          imp_fftcd(mp_fftcd(i)) = i
       enddo
    else
       do i=1,nfft_div_size
          imp_fftcd(mp_fftcd(i)) = i
       enddo
    end if
!!$ <==
#endif
    allocate(ica(nfft_div_size));ica=0
    allocate(icb(nfft_div_size));icb=0
    allocate(icc(nfft_div_size));icc=0
    call build_icaicbicc()
    Allocate(rho(nfft_div_size));rho=0.d0
    allocate(grad(nfft_div_size));grad=0.d0
    Allocate(theta_G(nfft_div_size))
    Allocate(theta_R(nfft_div_size))
    if(.not.oneshot)then
       Allocate(dtheta_R(np_nq,nfft_div_size))
       Allocate(ddtheta_R(np_nq,nfft_div_size))
    endif
    Allocate(theta_G_a(nfft_div_size))
    Allocate(theta_G_b(nfft_div_size))
    allocate(phi_ab(0:nr12))
    if(eval_kernel_by_interpolation.and.firstcall) allocate(phidD(0:ndel,-1:nphiD+1))
    allocate(qar(nq0));qar = 0.d0
    allocate(q2ar(nq0,nq0));q2ar = 0.d0
    !    allocate(theta_G_ab(nq0,-(na/2-1):na/2,-(nb/2-1):nb/2,-(nc/2-1):nc/2))
    if(sw_save_memory_vdw) then
       allocate(theta_G_ab(np_nq,nfft_div_size))
    else
       allocate(theta_G_ab(nq0,nfft_div_size))
    endif
    if(.not.oneshot)then
       !       allocate(ualpha_g(nq0,-(na/2-1):na/2,-(nb/2-1):nb/2,-(nc/2-1):nc/2))
       !       allocate(ualpha_r(nq0,na,nb,nc))
       allocate(dFdrho(nfft_div_size))
       allocate(dFddrho(nfft_div_size))
    endif

    allocate(rkar(nfft_div_size))
!    allocate(rkar(nabc))

#if 0
    Ta = na*DSQRT(aa(1,1)**2 + aa(1,2)**2 + aa(1,3)**2)
    Tb = nb*DSQRT(aa(2,1)**2 + aa(2,2)**2 + aa(2,3)**2)
    Tc = nc*DSQRT(aa(3,1)**2 + aa(3,2)**2 + aa(3,3)**2)
    Do ciz = -(nc/2-1),nc/2
       Do ciy = -(nb/2-1),nb/2
          Do cix = -(na/2-1),na/2
             i = get_index_1D(cix+na/2,ciy+nb/2,ciz+nc/2)
             rkar(i) = &
                  &    DSQRT((dble(cix)/Ta)**2 + (dble(ciy)/Tb)**2 + (dble(ciz)/Tc)**2)
          Enddo
       Enddo
    Enddo
#else
    Do i=1, 3
       bb(1,i) = rltv(i,1) /PAI /2.0d0
       bb(2,i) = rltv(i,2) /PAI /2.0d0
       bb(3,i) = rltv(i,3) /PAI /2.0d0
    End Do

    do i=1,nfft_div_size
       c1 = ica(i)
       c2 = icb(i)
       c3 = icc(i)
       ca = c1-na/2
       cb = c2-nb/2
       cc = c3-nc/2
       !cix = MOD(ca+na,na)+1
       !ciy = MOD(cb+nb,nb)+1
       !ciz = MOD(cc+nc,nc)+1
       !cir = get_index_1D(ca+na/2,cb+nb/2,cc+nc/2)
       vec(1:3) = dble(ca) *bb(1,1:3) &
            &    +dble(cb) *bb(2,1:3) &
            &    +dble(cc) *bb(3,1:3)
       rkar(i) = DSQRT( vec(1)**2 +vec(2)**2 +vec(3)**2 )
    enddo
!    Do ciz = -(nc/2-1),nc/2
!       Do ciy = -(nb/2-1),nb/2
!          Do cix = -(na/2-1),na/2
!             vec(1:3) = dble(cix) *bb(1,1:3) &
!                  &    +dble(ciy) *bb(2,1:3) &
!                  &    +dble(ciz) *bb(3,1:3)
!             i = get_index_1D(cix+na/2,ciy+nb/2,ciz+nc/2)
!             rkar(i) = DSQRT( vec(1)**2 +vec(2)**2 +vec(3)**2 )
!          End Do
!       End Do
!    End Do
#endif

  end subroutine alloc_vdw

  subroutine build_icaicbicc()
    integer :: i,j
    integer :: ind(3)
    do i=1,nfft_div_size
       j = mp_fftcd(i)
       ind = get_index_3D(j)
       ica(i) = ind(1)
       icb(i) = ind(2)
       icc(i) = ind(3)
    enddo
  end subroutine build_icaicbicc

  subroutine finalize_vdwdf()
    call dealloc()

  contains
    subroutine dealloc()
      deallocate(ica)
      deallocate(icb)
      deallocate(icc)
      deallocate(mp_fftcd)
      deallocate(imp_fftcd)
      deallocate(rho)
      deallocate(grad)
      deallocate(theta_G)
      deallocate(theta_R)
      if(.not.oneshot)then
         deallocate(dtheta_R)
         deallocate(ddtheta_R)
      endif
      deallocate(theta_G_a)
      deallocate(theta_G_b)
      deallocate(phi_ab)
      if(eval_kernel_by_interpolation.and.oneshot) deallocate(phidD)
      deallocate(qar)
      deallocate(q2ar)
      deallocate(theta_G_ab)
      if(.not.oneshot)then
         !         deallocate(ualpha_g)
         !         deallocate(ualpha_r)
         deallocate(dFdrho)
         deallocate(dFddrho)
      endif
      deallocate(rkar)
      !
      if ( allocated( cgrad ) ) deallocate( cgrad )
    end subroutine dealloc
  end subroutine finalize_vdwdf

  subroutine build_lookup_table(ndel,nphiD,phidD)
    use progress_bar
    implicit none
    integer, intent(in) :: ndel,nphiD
    real(kind=DP), intent(out), dimension(0:ndel,-1:nphiD+1) :: phidD
    integer cdel,cphiD
    real(kind=DP) :: del,phiD,di,dk,tmp
    real(kind=DP) :: ddel,dphiD
    integer :: ierr
    integer :: id_sname=-1

    call tstatc0_begin('build_lookup_table ',id_sname,1)

    phidD = 0.d0
    ddel = 1.d0 /dble(ndel)
    dphiD = (q0max*eta1-q0min*etai)/(dble(nphiD))

    if(lprog) then
       call reset_progress()
       call set_end(int(floor(dble(ndel+1)/dble(npes))))
    endif

    do cdel=0,ndel
       if(mod(cdel,npes)/=mype) cycle
       if(lprog) call progress()
       do cphiD = -1,nphiD+1
          del = dble(cdel)*ddel
          phiD = q0min*etai + dble(cphiD)*dphiD
          di = phiD*(1.d0+del)
          dk = phiD*(1.d0-del)
          Call kernel_phi(di,dk,tmp)
          phidD(cdel,cphiD) = tmp
       enddo
    enddo
    if(npes>1) &
         & call mpi_allreduce(MPI_IN_PLACE,phidD,(ndel+1)*(nphiD+2),mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
    call tstatc0_end(id_sname)
  end subroutine build_lookup_table

  subroutine build_theta()
    use progress_bar, only : reset_progress,progress,set_end
    integer :: cqa,ierr
!!$    integer :: imax, imin
    integer :: id_sname = -1
    real(kind=DP),allocatable,dimension(:) :: tmpdr,tmpddr

    call tstatc0_begin('build_theta ',id_sname,1)

    !+++++++++++++++ Execute FFT for theta_R_ab ++++++++++++++++++++
    if(printable.and.lprog) &
         & write(nfout,'(a)') 'building theta (spline coefficient x rho) and their Fourier transforms ...'
    if(lprog)then
       call reset_progress()
       call set_end(int(floor(dble(nq0)/dble(nrank_ke))))
    endif
    theta_G_ab(:,:) = (0.d0,0.d0)

    if(.not.oneshot)then
       allocate(tmpdr(nfft_div_size));tmpdr=0.d0
       allocate(tmpddr(nfft_div_size));tmpddr=0.d0
       dtheta_R=0.d0;ddtheta_R=0.d0
    else
       allocate(tmpdr(1))
       allocate(tmpddr(1))
    endif
!!$    write(nfout,'(" ** build_theta ** ista_nq,iend_nq = ",2i10)') ista_nq,iend_nq
!!$    imin = 1000
!!$    imax = -1
!!$    do cqa = ista_nq,iend_nq
!!$       if(imax<map_z_nq(cqa)) imax = map_z_nq(cqa)
!!$       if(imin>map_z_nq(cqa)) imin = map_z_nq(cqa)
!!$    end do
!!$    write(nfout,'(" **   map_z_nq = [",i8,", : ",i8,"]")') imin, imax
!!$    write(nfout,'(" **   nfft_div_size = ",i5)') nfft_div_size
!!$    if(sw_save_memory_vdw) then
!!$       write(nfout,'(" ** np_nq = ",i8)') np_nq
!!$    else
!!$       write(nfout,'(" ** nq0 = ",i8)') nq0
!!$    end if

    Do cqa = ista_nq,iend_nq
       if(lprog) call progress()
       Call theta_ab(nfft_div_size,cqa,nq0,q0min,q0max,rho,grad,rhomin,theta_R,tmpdr,tmpddr)
!!$       write(nfout,'(" ** build_theta na,nb,nc = ",3i8)') na,nb,nc
!!$       write(nfout,'("               size of theta_G, theta_R = ",2i8)') nfft_div_size,nfft_div_size
       Call RtoG(na,nb,nc,theta_R,theta_G)
       if(sw_save_memory_vdw)then
          theta_G_ab(map_z_nq(cqa),:) = theta_G(:)
       else
          theta_G_ab(cqa,:) = theta_G(:)
       endif
       if(.not.oneshot)then
          dtheta_R  (map_z_nq(cqa),:) = tmpdr  (:)
          ddtheta_R (map_z_nq(cqa),:) = tmpddr (:)
       endif
    Enddo

    if(.not.sw_save_memory_vdw)then
       call mpi_allreduce(MPI_IN_PLACE,theta_G_ab,nq0*nfft_div_size,mpi_double_complex, &
            &    mpi_sum,mpi_g_world,ierr)
    endif
    deallocate(tmpdr)
    deallocate(tmpddr)

    if(printable.and.lprog) write(nfout,'(a)') '... done'
    call tstatc0_end(id_sname)
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  end subroutine build_theta

  subroutine phiab_by_interpl(ca,cb,nr12,phi_ab)
    Implicit none
    integer, intent(in) :: ca,cb,nr12
    real(kind=DP), intent(out), dimension(0:nr12) :: phi_ab

    integer :: i,cr12,ierr
    real(kind=DP) :: r12,dr12,qa,qb,qab,phiD,del,ddel,phix,phiy
    integer :: cdel,cphiD

    real(kind=DP) :: rs,di,dk,phid_s,phid_s1,d_phid_s,phi2,phi4
    real(kind=DP) :: dr=0.001d0

    real(kind=DP) :: phi1u
    integer :: id_sname=-1

    call tstatc0_begin('phiab ',id_sname,1)

    phi_ab = 0.d0
    dr12 = r12max/DBLE(nr12)
    qa = qar(ca)
    qb = qar(cb)
    qab = DSQRT(qa**2+qb**2)

    ! Coefficients phi2 and phi4 in the local part is determined to
    ! match the non-local part in value and slope at d=d_s.

    rs = ds/qab

    di = qa*rs
    dk = qb*rs
    Call kernel_phi(di,dk,phid_s)

    di = qa*(rs+dr)
    dk = qb*(rs+dr)
    Call kernel_phi(di,dk,phid_s1)

    d_phid_s = (phid_s1 - phid_s)/dr

    phi2 = ( 2.d0/ds**2)*(phid_s-phi0) - (rs/(2.d0*ds**2))*d_phid_s
    phi4 = (-1.d0/ds**4)*(phid_s-phi0) + (rs/(2.d0*ds**4))*d_phid_s

    i = DINT((ds/qab)/dr12)
    ! Non-local part of phi_ab(r12)
    if(i.ge.nr12)return
    ddel = 1/dble(ndel)
    do cr12 = i+1,nr12
       if(mod(cr12,nrank_g) /= myrank_g ) cycle
       r12 = DBLE(cr12)*dr12

       di = qa*r12
       dk = qb*r12

       del  = dabs(di-dk)/(di+dk)
       phiD = 0.5d0*(di+dk)
       cdel = DINT(del/ddel)
       cphiD = DINT(DBLE(nphiD)*(phiD-q0min*etai)/(q0max*eta1-q0min*etai))
       !     write(6,*) 'cdel, cphiD: ',cdel,cphiD
       if(cdel.ge.ndel.or.cphiD.ge.nphiD.or.cdel.lt.0.or.cphiD.lt.-1)cycle
       phix = del/ddel - dble(cdel)
       phiy = dble(nphiD)*(phiD-q0min*etai)/(q0max*eta1-q0min*etai) - dble(cphiD)
       phi1u =  (1.d0-phix)*(1.d0-phiy) * phidD(cdel  ,cphiD  ) &
            &              +       phix *(1.d0-phiy) * phidD(cdel+1,cphiD  ) &
            &              +(1.d0-phix)*   phiy  * phidD(cdel  ,cphiD+1) &
            &              +      phix *   phiy  * phidD(cdel+1,cphiD+1)
       phi_ab(cr12) = phi1u
    enddo

    ! Local part of phi_ab(r12)
    Do cr12 = 0,i
       if(mod(cr12,nrank_g) /= myrank_g ) cycle
       r12 = DBLE(cr12)*dr12

       !d1 = qa*r12
       !d2 = qb*r12
       phiD = qab*r12
       phi_ab(cr12) = phi0 + phi2*phiD**2 + phi4*phiD**4
    Enddo

    if(nrank_g>1) &
    & call mpi_allreduce(MPI_IN_PLACE,phi_ab,nr12+1,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
    call tstatc0_end(id_sname)
  end subroutine phiab_by_interpl

  !** SUBROUTINE phiab **************************************************************************
  Subroutine phiab(ca,cb,nr12,phi_ab)
    Implicit none

    !++++++++++++++++++++++++++ VARIABLES ++++++++++++++++++++++++++
    Integer,intent(in) :: ca,cb,nr12
    Real(kind=DP),intent(out) :: phi_ab(0:nr12)

    ! Internal valuables
    Integer cr12,i,ierr
    Real(kind=DP) qab,qa,qb,r12,rs,phiD,dr12,phi2,phi4,phid_s,phid_s1,d_phid_s,d1,d2,abs_d
    Real(kind=DP) di,dk,tmp,dr
    Parameter(dr=0.001d0)
    !+++++++++++++++++++++ End VARIABLES +++++++++++++++++++++++++++
    integer :: id_sname=-1
    call tstatc0_begin('phiab ',id_sname,1)

    dr12 = r12max/DBLE(nr12)
    qa = qar(ca)
    qb = qar(cb)
    qab = DSQRT(qa**2+qb**2)

    ! Coefficients phi2 and phi4 in the local part is determined to
    ! match the non-local part in value and slope at d=d_s.

    rs = ds/qab

    di = qa*rs
    dk = qb*rs
    Call kernel_phi(di,dk,phid_s)

    di = qa*(rs+dr)
    dk = qb*(rs+dr)
    Call kernel_phi(di,dk,phid_s1)

    d_phid_s = (phid_s1 - phid_s)/dr

    phi2 = ( 2.d0/ds**2)*(phid_s-phi0) - (rs/(2.d0*ds**2))*d_phid_s
    phi4 = (-1.d0/ds**4)*(phid_s-phi0) + (rs/(2.d0*ds**4))*d_phid_s

    i = DINT((ds/qab)/dr12)
    ! Non-local part of phi_ab(r12)
    if(i.ge.nr12)return

    Do cr12 = i+1,nr12
       if(mod(cr12,nrank_g) /= myrank_g) cycle
       r12 = DBLE(cr12)*dr12

       di = qa*r12
       dk = qb*r12

       Call kernel_phi(di,dk,tmp)
       phi_ab(cr12) = tmp

    Enddo

    ! Local part of phi_ab(r12)
    Do cr12 = 0,i
       if(mod(cr12,nrank_g) /= myrank_g) cycle
       r12 = DBLE(cr12)*dr12

       d1 = qa*r12
       d2 = qb*r12
       phiD = qab*r12
       phi_ab(cr12) = phi0 + phi2*phiD**2 + phi4*phiD**4
    Enddo
    if(nrank_g>1) &
    & call mpi_allreduce(MPI_IN_PLACE,phi_ab,1+nr12,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
    call tstatc0_end(id_sname)
  End subroutine phiab
  !** End SUBROUTINE phiab **********************************************************************

  subroutine vdWdf()
    if(sw_save_memory_vdw) then
       call vdWdf_core()
       !call vdWdf_core2()
    else
       call vdWdf_core_org()
    endif
  end subroutine vdWdf

  subroutine vdWdf_core_org()
    use progress_bar
    real(kind=DP) :: fac
    integer :: cqaa,cqa,cqb,ierr
    integer :: i,ic,i2,i3
    complex(kind=DP), allocatable, dimension(:) :: tmpug,tmpur
    integer :: id_sname = -1
    integer :: id_sname3 = -1

    call tstatc0_begin('vdWdf_core ',id_sname,1)
    Ecnl_12 = 0.0d0
    if(printable.and.lprog) &
         & write(nfout,'(a)') 'performing the core operation : Ecnl = sum theta_a x phi_ab x theta_b ...'
    if(lprog)then
       call reset_progress()
       call set_end(int(floor(dble(nq0)/dble(nrank_ke))))
    endif
    if(.not.oneshot)then
       allocate(tmpug(nfft_div_size));tmpug=0.d0
       allocate(tmpur(nfft_div_size));tmpur=0.d0
       !       ualpha_g = (0.d0,0.d0)
    else
       allocate(tmpug(1))
    endif
    if(.not.oneshot)then
       dFdrho=0.d0;dFddrho=0.d0
    endif

    do cqaa=1,np_nq
       cqa = real_index(cqaa)
       theta_G_a(:) = theta_G_ab(cqa,:)
       if (.not.oneshot)  tmpug=(0.d0,0.d0)
       if (lprog) call progress()
       call tstatc0_begin('vdw_core_core ',id_sname3,1)
       do cqb=1,nq0
          theta_G_b(:) = theta_G_ab(cqb,:)
          fac = 1.d0
          if(eval_kernel_by_interpolation) then
             call phiab_by_interpl(cqa,cqb,nr12,phi_ab)
          else
             Call phiab(cqa,cqb,nr12,phi_ab)
          endif
          Call convolution_3d_by_fft(&
               &    na,nb,nc,cqa,cqb,nr12,phi_ab,theta_G_a,theta_G_b,fac,Ecnl_12_ab,tmpug)
          Ecnl_12 = Ecnl_12 + Ecnl_12_ab
       enddo
       call tstatc0_end(id_sname3)
       if(.not.oneshot)then
          call GtoR(na,nb,nc,tmpur,tmpug)
          dFdrho(:) = dFdrho(:) + dble(tmpur(:))*dtheta_R(cqaa,:)
          dFddrho(:) = dFddrho(:) + dble(tmpur(:))*ddtheta_R(cqaa,:)
       endif
    enddo

    if(.not.oneshot)then
       if(nrank_ke>1) then
          call mpi_allreduce(MPI_IN_PLACE,dFdrho, &
               &             nfft_div_size,mpi_double_precision,mpi_sum,mpi_g_world,ierr)
          call mpi_allreduce(MPI_IN_PLACE,dFddrho, &
               &             nfft_div_size,mpi_double_precision,mpi_sum,mpi_g_world,ierr)
          deallocate(tmpur)
       endif
       deallocate(tmpug)
    endif

    if(npes>1) then
       call mpi_allreduce(MPI_IN_PLACE,Ecnl_12,1,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
    endif
    if(printable.and.lprog) write(nfout,'(a)') '... done'
    call tstatc0_end(id_sname)

  end subroutine vdWdf_core_org

  subroutine vdWdf_core()
    use progress_bar
    real(kind=DP) :: fac, grho, rtmp, esum
    integer :: cqaa,cqa,cqb,ierr
    integer :: i,ic,i2,i3, cir, cix, ciy ,ciz
    complex(kind=DP), allocatable, dimension(:) :: tmpug,tmpur
    complex(kind=DP), allocatable, dimension(:,:) :: theta_buf_s,theta_buf_r
    integer :: ipos,i0,i1,isend,irecv,ireq,ireqr
    integer, allocatable, dimension(:) :: ista
    integer :: id_sname = -1, id_sname2 = -1, id_sname3 = -1

    Ecnl_12 = 0.0d0
    esum = 0.0d0

    call tstatc0_begin('vdWdf_core ',id_sname,1)
    if(printable.and.lprog) &
         & write(nfout,'(a)') 'performing the core operation : Ecnl = sum theta_a x phi_ab x theta_b ...'

    if(lprog)then
       call reset_progress()
       call set_end(int(floor(dble(nq0)/dble(nrank_ke))))
    endif

    if(nrank_ke>1)then
       allocate(ista(MPI_STATUS_SIZE))
       allocate(theta_buf_r(1:mp_nq,nfft_div_size));theta_buf_r=(0.d0,0.d0)
    endif

    if(.not.oneshot)then
       allocate(tmpug(nfft_div_size));tmpug=0.d0
       allocate(tmpur(nfft_div_size));tmpur=0.d0
       !       ualpha_g = (0.d0,0.d0)
    else
       allocate(tmpug(1))
       allocate(tmpur(1))
    endif
    if(.not.oneshot)then
       dFdrho=0.d0;dFddrho=0.d0
    endif

    irecv = myrank_ke+1
    if(irecv.ge.nrank_ke) irecv = irecv-nrank_ke
    isend = myrank_ke-1
    if(isend.lt.0)        isend = isend+nrank_ke

    do cqaa=1,mp_nq
       cqa = real_index(cqaa)
       if (nrank_ke>1)        theta_buf_r(1:np_nq,:) = theta_G_ab(1:np_nq,:)
       if (cqaa.le.np_nq) theta_G_a(:) = theta_G_ab(cqaa,:)

       if (.not.oneshot)  tmpug=(0.d0,0.d0)
       if (lprog)       call progress()

       call tstatc0_begin('vdw_core_core ',id_sname3,1)

       do i0=0,nrank_ke-1
          ipos = i0+myrank_ke
          if(ipos.ge.nrank_ke) ipos = ipos-nrank_ke
          if(ipos.lt.0)        ipos = ipos+nrank_ke

          if(cqaa.le.np_nq)then
             do i1=1,nel_nq(ipos)
                cqb = real_index(i1,ipos)
                if(nrank_ke>1)then
                   theta_G_b(:) = theta_buf_r(i1,:)
                else
                   theta_G_b(:) = theta_G_ab(i1,:)
                endif

                fac = 1.d0
                if(eval_kernel_by_interpolation) then
                   call phiab_by_interpl(cqa,cqb,nr12,phi_ab)
                else
                   Call phiab(cqa,cqb,nr12,phi_ab)
                endif
                Call convolution_3d_by_fft( na,nb,nc,cqa,cqb,nr12,phi_ab, &
                     &                      theta_G_a,theta_G_b,fac,Ecnl_12_ab,tmpug)
                Ecnl_12 = Ecnl_12 + Ecnl_12_ab
             enddo
          endif
          if(nrank_ke>1.and.i0.ne.nrank_ke-1)then
             call tstatc0_begin('vdWdf_core (comm) ',id_sname2,1)
             allocate(theta_buf_s(1:mp_nq,nfft_div_size))
             theta_buf_s = theta_buf_r
             call mpi_sendrecv(theta_buf_s,mp_nq*nfft_div_size,mpi_double_complex,isend,0, &
                  &            theta_buf_r,mp_nq*nfft_div_size,mpi_double_complex,irecv,0, &
                  &            mpi_g_world,ista,ierr)
             deallocate(theta_buf_s)
             call tstatc0_end(id_sname2)
          endif
       enddo
       call tstatc0_end(id_sname3)
       if(.not.oneshot) call GtoR(na,nb,nc,tmpur,tmpug)
       if(.not.oneshot.and.cqaa.le.np_nq)then
#if 1
          dFdrho(:) = dFdrho(:) + dble(tmpur(:))*dtheta_R(cqaa,:)
          if ( .not. grad_rho_eq_0 ) then
             do cir = 1,nfft_div_size
                grho = grad(cir)
                rtmp = rho(cir)
                if ( grho > 1.0d-6 ) then
                   dFddrho(cir) = dFddrho(cir) &
                        &                + dble(tmpur(cir)) &
                        &                 *ddtheta_R(cqaa,cir) /grho
                endif
             end do
          endif
#else
          dFdrho(:) = dFdrho(:) + dble(tmpur(:))*dtheta_R(cqaa,:)
          if ( .not. grad_rho_eq_0 ) then
             dFddrho(:) = dFddrho(:) + dble(tmpur(:))*ddtheta_R(cqaa,:)
          endif
#endif

          if ( iprixc >= 2 ) then
             do cir = 1,nabc
                esum = esum +dble(tmpur(cir)) *theta_R(cir)
             end do
          endif
       endif
    enddo

    if(.not.oneshot)then
       if(nrank_ke>1) then
          call mpi_allreduce(MPI_IN_PLACE,dFdrho, &
               &             nfft_div_size,mpi_double_precision,mpi_sum,mpi_g_world,ierr)
          call mpi_allreduce(MPI_IN_PLACE,dFddrho, &
               &             nfft_div_size,mpi_double_precision,mpi_sum,mpi_g_world,ierr)
          deallocate(tmpur)
       endif
       deallocate(tmpug)
    endif

    if(npes>1) then
       call mpi_allreduce(MPI_IN_PLACE,Ecnl_12,1,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
    endif

    if ( iprixc >=2 ) then
       esum = esum *univol *rinplw /2.0d0
       call mpi_allreduce(MPI_IN_PLACE,Esum,1,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
       write(nfout,*) "vdw nonolocal E12 "
       write(nfout,*) "evaluation in real       space : ", esum
       write(nfout,*) "           in reciprocal space : ", ecnl_12
    endif

    if(nrank_ke>1)then
       deallocate(ista)
       deallocate(theta_buf_r)
    endif
    if(printable.and.lprog) write(nfout,'(a)') '... done'
    call tstatc0_end(id_sname)
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  end subroutine vdWdf_core

  subroutine vdWdf_core2()
    use progress_bar
    real(kind=DP) :: fac
    integer :: cqbb,cqa,cqb,ierr
    integer :: i,ic,i2,i3
    complex(kind=DP), allocatable, dimension(:) :: tmpug,tmpur
    complex(kind=DP), allocatable, dimension(:,:) :: theta_buf_s,theta_buf_r
    real(kind=DP), allocatable, dimension(:,:) :: dtheta_buf
    integer :: ipos,i0,i1,isend,irecv,ireq,ireqr
    integer, allocatable, dimension(:) :: ista
    integer :: id_sname = -1
    integer :: id_sname2 = -1

    call tstatc0_begin('vdWdf_core2 ',id_sname,1)
    Ecnl_12 = 0.0d0
    if(printable.and.lprog) &
         & write(nfout,'(a)') 'performing the core operation : Ecnl = sum theta_a x phi_ab x theta_b ...'
    if(lprog)then
       call reset_progress()
       call set_end(int(floor(dble(nq0)/dble(nrank_ke))))
    endif
    if(nrank_ke>1)then
       allocate(ista(MPI_STATUS_SIZE))
       allocate(theta_buf_r(1:mp_nq,nfft_div_size));theta_buf_r=(0.d0,0.d0)
       if(.not.oneshot)then
          allocate(dtheta_buf(mp_nq,nfft_div_size));dtheta_buf=0.d0
       endif
    endif
    if(.not.oneshot)then
       allocate(tmpug(nfft_div_size));tmpug=0.d0
       allocate(tmpur(nfft_div_size));tmpur=0.d0
       !       ualpha_g = (0.d0,0.d0)
    else
       allocate(tmpug(1))
    endif
    if(.not.oneshot)then
       dFdrho=0.d0;dFddrho=0.d0
    endif

    irecv = myrank_ke+1
    if(irecv.ge.nrank_ke) irecv = irecv-nrank_ke
    isend = myrank_ke-1
    if(isend.lt.0)        isend = isend+nrank_ke

    if (nrank_ke>1) theta_buf_r(1:np_nq,:) = theta_G_ab(1:np_nq,:)
    do i0=0,nrank_ke-1
       ipos = i0+myrank_ke
       if(ipos.ge.nrank_ke) ipos = ipos-nrank_ke
       if(ipos.lt.0)        ipos = ipos+nrank_ke
       do i1=1,nel_nq(ipos)
          cqa = real_index(i1,ipos)
          if(nrank_ke>1)then
             theta_G_a(:) = theta_buf_r(i1,:)
          else
             theta_G_a(:) = theta_G_ab(i1,:)
          endif
          do cqbb=1,np_nq
             theta_G_b(:) = theta_G_ab(cqbb,:)
             fac = 1.d0
             cqb = real_index(cqbb)
             if(eval_kernel_by_interpolation) then
                call phiab_by_interpl(cqa,cqb,nr12,phi_ab)
             else
                Call phiab(cqa,cqb,nr12,phi_ab)
             endif
             Call convolution_3d_by_fft(&
                  &    na,nb,nc,cqa,cqb,nr12,phi_ab,theta_G_a,theta_G_b,fac,Ecnl_12_ab,tmpug)
             Ecnl_12 = Ecnl_12 + Ecnl_12_ab
          enddo
          if(.not.oneshot)then
             !             call mpi_allreduce(MPI_IN_PLACE,tmpug,na*nb*nc,mpi_double_complex,mpi_sum,MPI_CommGroup,ierr)
             call GtoR(na,nb,nc,tmpur,tmpug)
             dFdrho(:) = dFdrho(:) + dble(tmpur(:))*dtheta_R(i1,:)
             dFddrho(:) = dFddrho(:) + dble(tmpur(:))*ddtheta_R(i1,:)
          endif
       enddo
       if(nrank_ke>1)then
          call tstatc0_begin('vdWdf_core (comm) ',id_sname2,1)
          allocate(theta_buf_s(1:mp_nq,nfft_div_size))
          theta_buf_s = theta_buf_r
          call mpi_sendrecv(theta_buf_s,mp_nq*nfft_div_size,mpi_double_complex,isend,0, &
               &            theta_buf_r,mp_nq*nfft_div_size,mpi_double_complex,irecv,0, &
               &            mpi_g_world,ista,ierr)
          if(.not.oneshot)then
             dtheta_buf = dtheta_R
             call mpi_sendrecv(dtheta_buf,mp_nq*nfft_div_size,mpi_double_precision,isend,0, &
                  &            dtheta_R,  mp_nq*nfft_div_size,mpi_double_precision,irecv,0, &
                  &            mpi_g_world,ista,ierr)
             dtheta_buf = ddtheta_R
             call mpi_sendrecv(dtheta_buf,mp_nq*nfft_div_size,mpi_double_precision,isend,0, &
                  &            ddtheta_R, mp_nq*nfft_div_size,mpi_double_precision,irecv,0, &
                  &            mpi_g_world,ista,ierr)
          endif
          deallocate(theta_buf_s)
          call tstatc0_end(id_sname2)
       endif
    enddo

    if(.not.oneshot)then
       if(nrank_ke>1) then
          call mpi_allreduce(MPI_IN_PLACE,dFdrho, &
               &             nfft_div_size,mpi_double_precision,mpi_sum,mpi_g_world,ierr)
          call mpi_allreduce(MPI_IN_PLACE,dFddrho, &
               &             nfft_div_size,mpi_double_precision,mpi_sum,mpi_g_world,ierr)
          deallocate(tmpur)
       endif
       deallocate(tmpug)
    endif

    if(npes>1) then
       call mpi_allreduce(MPI_IN_PLACE,Ecnl_12,1,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
    endif
    if(nrank_ke>1)then
       deallocate(ista)
       deallocate(theta_buf_r)
       if(.not.oneshot) deallocate(dtheta_buf)
    endif
    if(printable.and.lprog) write(nfout,'(a)') '... done'
    call tstatc0_end(id_sname)
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  end subroutine vdWdf_core2

  subroutine corrections()
!!$if(printable) write(nfout,'(a)') 'calculating correction terms and contribution from the LDA ...'
    if(printable.and.lprog) write(nfout,'(a)') 'calculating correction terms ...'
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    Call piDphi(Ecnl_3,Ecnl_3s)
!!$    Call cLDA(na,nb,nc,rho,rhomin,dv,EcLDA)
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if(printable.and.lprog) write(nfout,'(a)') '... done'
  end subroutine corrections

  subroutine theta_ab(nxyz,ca,nq0,q0min,q0max,rho,grad,rhomin,theta_R,dtheta_R,ddtheta_R)
    Implicit none

    !++++++++++++++++++++++++++ VARIABLES ++++++++++++++++++++++++++
    Integer, intent(in) :: nxyz,ca,nq0
    Real(kind=DP), intent(in)  :: q0min,q0max,rhomin
    Real(kind=DP), intent(in)  :: rho(nxyz),grad(nxyz)
    real(kind=DP), intent(out) :: theta_R(nxyz),dtheta_R(nxyz),ddtheta_R(nxyz)

    ! Internal valuables
    Integer cir
    Real(kind=DP) ni,dni,q0,dqdn,dqddn

    real(kind=DP),allocatable,dimension(:) :: y,y2tmp
    real(kind=DP) :: yout,y1out
    integer :: id_sname = -1

    call tstatc0_begin('theta_ab ',id_sname,1)
    !+++++++++++++++++++++ End VARIABLES +++++++++++++++++++++++++++
    allocate(y(nq0));y=0.d0;y(ca)=1.d0
    allocate(y2tmp(nq0));y2tmp(:) = q2ar(ca,:)
    ! For the functions theta_R
    Do cir = 1,nxyz
       q0 = q0max
       ni = MAX(rho(cir),rhomin)
       dni = grad(cir)
       Call d_q0(ni,dni,q0min,q0max,q0,dqdn,dqddn)

       call cubic_spline(nq0,qar,y,y2tmp,q0,yout,y1out)
       theta_R  (cir) = yout*ni
       if(.not.oneshot)then
          dtheta_R (cir) = yout+ni*y1out*dqdn
          ddtheta_R(cir) = ni*y1out*dqddn
       endif
! ----
       if ( sw_use_WuGygi_method == ON ) then
          theta_R(cir) = theta_R(cir) *qar(ca) /q0
          if ( .not. oneshot ) then
             dtheta_R(cir) = dtheta_R(cir) *qar(ca) /q0 &
                  &                -theta_R(cir) *dqdn /q0
             ddtheta_R(cir) = ddtheta_R(cir) *qar(ca) /q0 &
                  &                -theta_R(cir) *dqddn /q0
          endif
       endif
! ----
    End do
    deallocate(y);   deallocate(y2tmp)
    call tstatc0_end(id_sname)
  end subroutine theta_ab

  !** SUBROUTINE piDphi **********************************************************************************
  Subroutine piDphi(Ecii,Ecii_s)
    implicit none


    !++++++++++++++++++++++++++ VARIABLES ++++++++++++++++++++++++++
    real(kind=DP), intent(out) :: Ecii,Ecii_s

    Real(kind=DP) da,db,a1,a2,dr
    Parameter (dr = 0.001d0)
    Integer  ci,cj,ck,ca,cb

    ! Gauss-Legendre integration
    Integer cD,nD
    Parameter (nD=10)
    Real(kind=DP) maxD,minD,dD,LD,PLD,LDxi(nD),LDwi(nD)
    Parameter (minD=0.d0)

    Integer  cix,ciy,ciz,cir,cjx,cjy,cjz,cjr
    Real(kind=DP)  di,dj,dk,d_di,d_dj,d_dk,eta,dx,dy,dz,n,ni,nj,nk
    Integer  i,j,k,zxp,zxm,zyp,zym,zzp,zzm
    Parameter (eta=0.00000001)

    Real(kind=DP) x,nnx,nny,nnz,nn2,r
    Real(kind=DP) nxp,nxm,nyp,nym,nzp,nzm
    Real(kind=DP) zx(-3:3),zy(-3:3),zz(-3:3),rn(3,-3:3)

    Real(kind=DP) phi

    Real(kind=DP) temp,rs,phid_s,phid_s1,d_phid_s,phi2,phi4

    ! The table of phi1D
    Integer c1D,n1D,ierr
    Parameter(n1D = 1000)
    Real(kind=DP) d1D,max1D,D,phix,phiy
    Real(kind=DP) phi1D(0:n1D+1)
    real(kind=DP) :: q0,dqdn,dqddn
    !++++++++++++++++++++++++ end VARIABLES ++++++++++++++++++++++++



    !---------------------- Calculation Start ----------------------
    ! Make the table of phi1D
    max1D = ds
    d1D = max1D/DBLE(n1D)
    Do c1D = 0,n1D+1
       D = DBLE(c1D)*d1D
       Call kernel_phi(D,D,phi)
       phi1D(c1D) = phi
    Enddo


    Ecii = 0
    Do cjr = 1,nfft_div_size
!       if(mod(cjr,npes)/=mype) cycle
!       cjx = 1+(cjr-1-MOD(cjr-1+nb*nc,nb*nc))/(nb*nc)
!       cjy = 1+(cjr-nb*nc*(cjx-1)-1-MOD(cjr-1+nc,nc))/nc
!       cjz = cjr-nc*(nb*(cjx-1)-1+cjy)

       n = MAX(rho(cjr),rhomin)

       Call d_q0(n,grad(cjr),q0min,q0max,q0,dqdn,dqddn)
       maxD=ds

       Call gauleg(minD,maxD,nD,LDxi,LDwi)

       temp = 0.d0
       Do cD=1,nD
          LD = LDxi(cD)

          c1D = DINT(LD/d1D)
          phix = (LD - d1D*DBLE(c1D))/d1D
          phiy = 1.d0 - phix
          temp = temp + LDwi(cD)*4.d0*pi*(LD**2) * (phiy*phi1D(c1D) + phix*phi1D(c1D+1))

       Enddo
       Ecii = Ecii + dv*0.5d0*(n**2)*temp/(q0**3)
    Enddo
    call mpi_allreduce(MPI_IN_PLACE,Ecii,1,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)

    Ecii_s = 0.d0
    Do cjr = 1,nfft_div_size
!       if(mod(cjr,npes)/=mype) cycle

       n = MAX(rho(cjr),rhomin)

       Call d_q0(n,grad(cjr),q0min,q0max,q0,dqdn,dqddn)
       rs = ds/DSQRT(q0**2+q0**2)
       di = q0*rs
       dk = q0*rs
       c1D = DINT(di/d1D)
       phix = (di - d1D*DBLE(c1D))/d1D
       phiy = 1.d0 - phix
       phid_s = phiy*phi1D(c1D) + phix*phi1D(c1D+1)

       di = q0*(rs+dr)
       dk = q0*(rs+dr)
       c1D = DINT(di/d1D)
       phix = (di - d1D*DBLE(c1D))/d1D
       phiy = 1.d0 - phix
       phid_s1 = phiy*phi1D(c1D) + phix*phi1D(c1D+1)

       d_phid_s = (phid_s1- phid_s)/dr

       phi2 = ( 2.d0/ds**2)*(phid_s-phi0) - (rs/(2.d0*ds**2))*d_phid_s
       phi4 = (-1.d0/ds**4)*(phid_s-phi0) + (rs/(2.d0*ds**4))*d_phid_s

       Ecii_s = Ecii_s + 0.5d0*4.d0*pi*dv*(n**2)*            &
            &                    (phi0*(rs**3)*(q0**0)/3.d0 +            &
            &                     phi2*(rs**5)*(q0**2)/5.d0 +            &
            &                     phi4*(rs**7)*(q0**4)/7.d0)


    Enddo
    call mpi_allreduce(MPI_IN_PLACE,Ecii_s,1,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
  End Subroutine piDphi
  !** End SUBROUTINE piDphi **********************

  !** SUBROUTINE cLDA ****************************
  Subroutine cLDA(na,nb,nc,rho,rhomin,dv,EcLDA)
    Implicit none

    !************************ Note *********************************
    ! This Algorism follows Dion's 1-shot method.
    !
    ! This program is a subroutine.
    ! This program calculates the correlation energy from LDA.
    ! The formula is given at Eq.(58) (p.93) of 'Theory of the
    !   Inhomogeneous Electron Gas' Lundqvist, March.
    !
    !
    ! Input
    !   rho(nrxyz,nsipn) : Total density
    !
    ! Output
    !   EcLDA : Correlation energy from LDA.
    !
    !
    !                            Written by Youky Ono in 2009/Jul.
    !***************************************************************

    !++++++++++++++++++++++++++ VARIABLES ++++++++++++++++++++++++++
    integer,intent(in)  :: na,nb,nc
    real(kind=DP),intent(in)  :: rho(nabc)
    real(kind=DP),intent(in)  :: rhomin
    real(kind=DP),intent(in)  :: dv
    real(kind=DP),intent(out) :: EcLDA
    integer  cjx,cjy,cjz,cjr,ierr
    integer  i,j,k

    real(kind=DP) rs,x,aB,ec,n
    parameter (aB=1.d0) ! aB=1[a.u.]

    real(kind=DP)  e,m
    parameter (e=1.d0,m=1.d0) ! Hatree atomic unit
    !++++++++++++++++++++++++ end VARIABLES ++++++++++++++++++++++++

    !---------------------- Calculation Start ----------------------
    EcLDA=0
    Do cjr = 1,nfft_div_size
!       if(mod(cjr,npes)/=mype) cycle
       n = MAX(rho(cjr),rhomin)

       rs = ((3.d0/(4.d0*pi*n))**(1.d0/3.d0))/aB
       x = rs/11.4d0
       ec = -0.0666d0*0.5d0*((1.d0+x**3)*DLOG(1.d0+1.d0/x)-x**2+x/2.d0-1.d0/3.d0)

       EcLDA = EcLDA + dv*n*ec
    Enddo
    if(npes>1) call mpi_allreduce(MPI_IN_PLACE,EcLDA,1,mpi_double_precision,mpi_sum,mpi_g_world,ierr)
  End Subroutine cLDA
  !** End SUBROUTINE cLDA **********************

  ! Execute FFT and transform theta_R to theta_G
  Subroutine RtoG(na,nb,nc,theta_R,theta_G)
    Implicit none

    !!include "fftw3.f"


    !++++++++++++++++++++++++++ VARIABLES ++++++++++++++++++++++++++
    Integer, intent(in) ::   na,nb,nc
    Integer :: cir
    Real(kind=DP) rx,ry,rz,kx,ky,kz,ra,rb,rc,ka,kb,kc,rk,r12,x,y,z,term,term1

    integer :: ca,cb,cc,cix,ciy,ciz,i,i2,c1,c2,c3,ierr
    integer :: m
    Real(kind=DP) theta_R(nfft_div_size)
    Complex(kind=DP) theta_G(nfft_div_size)

    ! FFTW3 !!!
    integer(kind=DP) :: plan
    Complex(kind=DP),allocatable :: temp_R(:),temp_G(:)

    real(kind=DP), allocatable, dimension(:) :: afft

    !+++++++++++++++++++++ End VARIABLES +++++++++++++++++++++++++++
    allocate(temp_R(nabc));temp_R=(0.d0,0.d0)
    allocate(temp_G(nabc));temp_G=(0.d0,0.d0)
    !***** FFT **************************************************
!!$    goto 1001
    do ca=1,nfft_div_size
       if(mp_fftcd(ca)>nabc) cycle
       temp_R(mp_fftcd(ca)) = DCMPLX(theta_R(ca))
    enddo
    call mpi_allreduce(MPI_IN_PLACE,temp_R,nabc,mpi_double_complex,mpi_sum,mpi_ke_world,ierr)

    allocate(afft(nfftp_nonpara));afft=0.d0


    call map_rho_to_afft(temp_R,afft)
    call m_FFT_CD0(nfout,afft,DIRECT)
    call map_afft_to_rho(afft,temp_G)

    ! FFTW3 !!!
!    call dfftw_plan_dft_3d(plan,na,nb,nc,temp_R(1),temp_G(1),FFTW_FORWARD,FFTW_ESTIMATE)
!    call dfftw_execute(plan)
!    call dfftw_destroy_plan(plan)
    do i=1,nfft_div_size
       c1 = ica(i)
       c2 = icb(i)
       c3 = icc(i)
       ca = c1-na/2
       cb = c2-nb/2
       cc = c3-nc/2
       cix = MOD(ca+na,na)+1
       ciy = MOD(cb+nb,nb)+1
       ciz = MOD(cc+nc,nc)+1
       theta_G(i) = temp_G(get_index_1D(cix,ciy,ciz))/dble(nabc)
    enddo

!    do ca = -(na/2-1),na/2
!       do cb = -(nb/2-1),nb/2
!          do cc = -(nc/2-1),nc/2
!             cix = MOD(ca+na,na)+1
!             ciy = MOD(cb+nb,nb)+1
!             ciz = MOD(cc+nc,nc)+1
!             i = get_index_1D(cix,ciy,ciz)
!             i2 = get_index_1D(ca+na/2,cb+nb/2,cc+nc/2)
!             if(imp_fftcd(i2)>0) theta_G(imp_fftcd(i2)) = temp_G(i) / DBLE(nabc)
!          enddo
!       enddo
!    enddo
    !***** END of FFT ******************************************

1001 continue
    if(allocated(temp_R)) deallocate(temp_R)
    if(allocated(temp_G)) deallocate(temp_G)
    if(allocated(afft))   deallocate(afft)
  End subroutine RtoG
  !** End SUBROUTINE RtoG ***********************************************************************

  subroutine GtoR(na,nb,nc,theta_R,theta_G)
    !    include "fftw3.f"
    integer,intent(in) :: na,nb,nc
    complex(kind=DP), dimension(nfft_div_size), intent(in) :: theta_G
    complex(kind=DP), dimension(nfft_div_size), intent(out) :: theta_R

    integer :: cir,ca,cb,cc,cix,ciy,ciz,i,i2,c1,c2,c3,ierr
    integer(kind=DP) :: plan
    complex(kind=DP),allocatable :: temp_R(:),temp_G(:)
    real(kind=DP),allocatable,dimension(:) :: afft
    integer :: id_sname=-1
    call tstatc0_begin('GtoR ',id_sname,1)

    allocate(temp_R(nabc))
    allocate(temp_G(nabc));temp_G = (0.d0,0.d0)
    do i=1,nfft_div_size
       c1 = ica(i)
       c2 = icb(i)
       c3 = icc(i)
       ca = c1-na/2
       cb = c2-nb/2
       cc = c3-nc/2
       cix = MOD(ca+na,na)+1
       ciy = MOD(cb+nb,nb)+1
       ciz = MOD(cc+nc,nc)+1
       temp_G(get_index_1D(cix,ciy,ciz)) = theta_G(i)
    enddo
    call mpi_allreduce(MPI_IN_PLACE,temp_G,nabc,mpi_double_complex,mpi_sum,mpi_ke_world,ierr)

!    Do ca = -(na/2-1),na/2
!       Do cb = -(nb/2-1),nb/2
!          Do cc = -(nc/2-1),nc/2
!             cix = MOD(ca+na,na)+1
!             ciy = MOD(cb+nb,nb)+1
!             ciz = MOD(cc+nc,nc)+1
!             i = get_index_1D(cix,ciy,ciz)
!             i2 = get_index_1D(ca+na/2,cb+nb/2,cc+nc/2)
!             temp_G(get_index_1D(cix,ciy,ciz)) = theta_G(get_index_1D(ca+na/2,cb+nb/2,cc+nc/2))
!          Enddo
!       Enddo
!    Enddo

    allocate(afft(nfftp_nonpara));afft = 0.d0
    call map_rhog_to_afft(temp_G,afft)
    call m_FFT_CD0(nfout,afft,INVERSE)
    call map_afft_to_rhog(afft,temp_R)
    do ca=1,nfft_div_size
       theta_R(ca) = temp_R(mp_fftcd(ca))
    enddo

    ! FFTW3 !!!
!    call dfftw_plan_dft_3d(plan,na,nb,nc,temp_G(1),temp_R(1),FFTW_BACKWARD,FFTW_ESTIMATE)
!    call dfftw_execute(plan)
!    call dfftw_destroy_plan(plan)
    !***** END of FFT ******************************************

    !    theta_R = temp_R/dble(na*nb*nc)
!    theta_R = temp_R

    deallocate(temp_R)
    deallocate(temp_G)
    deallocate(afft)
    call tstatc0_end(id_sname)
  end subroutine GtoR

  subroutine get_phi_ab_g(nr12,phi_ab,phi_ab_g)
    !    include "fftw3.f"
    integer, intent(in) :: nr12
    real(kind=DP), dimension(0:nr12), intent(in)  :: phi_ab
    real(kind=DP), dimension(0:nr12), intent(out) :: phi_ab_g
    complex(kind=DP), allocatable, dimension(:) :: cphiab_r
    complex(kind=DP), allocatable, dimension(:) :: cphiab_g
    integer :: cr,ck
    integer(kind=DP) :: plan
    real(kind=DP) :: dd,r12,term,rk
    real(kind=DP) :: rr
    complex(kind=DP) :: zsum
    complex(kind=DP), parameter :: zi = ( 0.0d0, 1.0d0 )

    integer :: id_sname = -1
    call tstatc0_begin('get_phi_ab_g ',id_sname,1)
    allocate(cphiab_r(0:nr12));cphiab_r=0.d0
    allocate(cphiab_g(0:nr12));cphiab_g=0.d0

    dd = r12max/dble(nr12)

    do cr=0,nr12
       rr = dble(cr) *dd
       cphiab_r(cr) = dcmplx( phi_ab(cr)*rr, 0.d0 )
    enddo

#if 0
    call dfftw_plan_dft_1d(plan,nr12+1,cphiab_r,cphiab_g,FFTW_BACKWARD,FFTW_ESTIMATE)
#else
    call dfftw_plan_dft_1d( plan, nr12, cphiab_r(0:nr12-1), cphiab_g(0:nr12-1), &
         &                  FFTW_BACKWARD, FFTW_ESTIMATE )
#endif
    call dfftw_execute(plan)
    call dfftw_destroy_plan(plan)

    do ck=1,nr12
       rk = dble(ck)/r12max
       phi_ab_g(ck) = 2.d0 *dd *dimag(cphiab_g(ck))/(dble(rk))
    enddo

    Do ck = 0,0
       rk = dd * dble(ck)
       term = 0.0d0
       Do cr = 0,nr12
          r12 = dd*dble(cr)
          term = term + phi_ab(cr) * (r12**2)
       Enddo
       phi_ab_g(ck) = 4.0d0*pi*dd * term
    End do

    deallocate(cphiab_r);   deallocate(cphiab_g)

    phi_ab_g = phi_ab_g / univol

    call tstatc0_end(id_sname)

  end Subroutine get_phi_ab_g

! ===============- For stress tensor =================
  subroutine initialize_vdwdf_stress( nspin, ispin, nfft_div_size, na, nb, nc, chgr, &
       &                              grad_rho, cgrad_rho, version_no )
    use progress_bar, only : set_printable

    integer, intent(in) :: nspin,ispin,nfft_div_size, na,nb,nc, version_no
    real(kind=DP), intent(in) :: chgr(nfft_div_size), grad_rho(nfft_div_size)
    real(kind=DP), intent(in) :: cgrad_rho(nfft_div_size,3)

    integer :: i,cir,cix,ciy,ciz,nrxyz,cir2
    real(kind=DP) :: q

    if ( version_no == 1 ) then
       Zab = Zab_0
    else if ( version_no == 2 ) then
       Zab = Zab_0 *2.222d0
    endif

    call set_printable(printable)
    call do_cell_params()

    rinplw = 1.d0/(dble(nabc))
    q0max = q0cut*1.01d0
!!!    q0min = 0.09d0
    nq0 = DINT(dLOG((q0max-q0min)*(lambda-1.d0)/dq_vdw+1.d0)/dLOG(lambda))+1
    maxk = dble(nr12)/r12max

    if(firstcall) call m_Parallel_init_mpi_nq(nfout,ipri,printable,nq0)
    if(firstcall) call print_vdw_parameters()

    call alloc_vdw()
    allocate( cgrad(nfft_div_size,3)); cgrad=0.d0
    rho(:) = chgr(:)
    grad(:) = grad_rho(:)
    cgrad(:,1:3) = cgrad_rho(:,1:3)

    if ( grad_rho_eq_0 )  grad = 0.0d0

    do i=1,nq0
       q = q0min + dq_vdw*(lambda**DBLE(i-1)-1.d0)/(lambda-1.d0)
       qar(i) = q
    enddo
    call spline0(nq0,qar,q2ar)

    if(eval_kernel_by_interpolation.and.firstcall) &
         &             call build_lookup_table(ndel,nphiD,phidD)
    firstcall = .false.
  end subroutine initialize_vdwdf_stress

  subroutine calc_stress_contrib2( na, nb, nc, cqa, cqb, nr12, phi_ab, &
       &                           theta_G_a, theta_G_b, fac, s_tensor, tmpug )
    Implicit none

    integer, intent(in) :: na,nb,nc,cqa,cqb,nr12
    complex(kind=DP), intent(in)  :: theta_G_a(nfft_div_size)
    complex(kind=DP), intent(in)  :: theta_G_b(nfft_div_size)
    real(kind=DP),    intent(in)  :: phi_ab(0:nr12)
    real(kind=DP),    intent(in)  :: fac
    real(kind=DP), intent(inout) :: s_tensor(3,3)
    complex(kind=DP), intent(inout) :: tmpug(nfft_div_size)

    Integer  i, u, v
    Integer cir,cix,ciy,ciz,cr,ck,c1,c2,c3

    Real(kind=DP) phix, phiy, Ta, Tb, Tc, bb(3,3), vec(3)
    Real(kind=DP) dk,rk,r12,term1, term2
    real(kind=DP) :: pi2rk, g1

    Real(kind=DP),allocatable :: core_G1(:), core_G2(:)
    Real(kind=DP),allocatable :: wk_x(:), wk_y(:)
    real(kind=DP) :: ctmp, z1

    integer :: id_sname = -1

    call tstatc0_begin('calc_stress_contrib2 ',id_sname,1)

    allocate(core_G1(0:nr12)); core_G1 = 0.d0
    allocate(core_G2(0:nr12)); core_G2 = 0.d0

    call get_phi_ab_g( nr12, phi_ab, core_G1 )

! ----
    dk = 1.0d0 /r12max

!    Do i=1, nr12 -1
!       z1 = ( core_G1(i) -core_G1(i-1) ) /dk
!       z2 = ( core_G1(i+1) -core_G1(i) ) /dk
!       core_G2(i) = ( z1 +z2 )/2.0d0
!    End Do

    allocate( wk_x(0:nr12) );     allocate( wk_y(0:nr12) )
    Do i=0, nr12
       wk_x(i) = dk *i
    End Do
    call init_cubic_spline( nr12+1, wk_x(0:nr12), core_G1(0:nr12), wk_y(0:nr12) )

    Do i=0, nr12
       call cubic_spline( nr12+1, wk_x(0:nr12), core_G1(0:nr12), wk_y(0:nr12), &
            &             wk_x(i), ctmp, core_G2(i) )
    End Do
    deallocate( wk_x ); deallocate( wk_y )

! ----
    s_tensor = 0.0d0

    Do i=1, 3
       bb(1,i) = rltv(i,1) /PAI /2.0d0
       bb(2,i) = rltv(i,2) /PAI /2.0d0
       bb(3,i) = rltv(i,3) /PAI /2.0d0
    End Do

    do i=1,nfft_div_size
       c1 = ica(i)
       c2 = icb(i)
       c3 = icc(i)
       cix = c1-na/2
       ciy = c2-nb/2
       ciz = c3-nc/2
       rk   = rkar(i)
       ck   = dint(rk/dk)
       phix = (rk - dk*DBLE(ck))/dk
       phiy = 1.0d0 - phix

       term1 = phiy*core_G1(ck) + phix*core_G1(ck+1)
       term2 = phiy*core_G2(ck) + phix*core_G2(ck+1)

       vec(1:3) = dble(cix) *bb(1,1:3) &
            &    +dble(ciy) *bb(2,1:3) &
            &    +dble(ciz) *bb(3,1:3)

       if ( ck > 0 ) then
          z1 = DCONJG( theta_G_a(i) ) *theta_G_b(i)  &
               &                                *term2 /rkar(i)
          Do v=1, 3
             Do u=1, 3
                s_tensor(u,v) = s_tensor(u,v) +z1 *vec(u) *vec(v)
             End Do
          End Do
       endif
       tmpug(i) = tmpug(i) &
            &              +theta_G_b(i) *term1 *univol
    enddo

    s_tensor = s_tensor *univol

    deallocate(core_G1);    deallocate(core_G2)
    call tstatc0_end(id_sname)

  end subroutine calc_stress_contrib2

  subroutine vdWdf_stress_tensor_core()
    use progress_bar

    real(kind=DP) :: fac
    integer :: cqaa,cqa,cqb,ierr
    integer :: i,ic,i2,i3
    complex(kind=DP), allocatable, dimension(:) :: tmpug,tmpur
    complex(kind=DP), allocatable, dimension(:,:) :: theta_buf_s,theta_buf_r

    real(kind=DP) :: s_tmp(3,3), grho, rtmp

    integer :: u, v, cix, ciy, ciz,cir
    integer :: ipos,i0,i1,isend,irecv,ireq,ireqr
    integer, allocatable, dimension(:) :: ista
    integer :: id_sname = -1
    integer :: id_sname2 = -1
    integer :: id_sname3 = -1

    call tstatc0_begin('vdWdf_stres_core ',id_sname,1)

    s_cnl1 = 0.0d0;  s_cnl2 = 0.0d0

    if(nrank_ke>1)then
       allocate(ista(MPI_STATUS_SIZE))
       allocate(theta_buf_r(1:mp_nq,nfft_div_size));
       theta_buf_r=(0.d0,0.d0)
    endif

    allocate(tmpug(nfft_div_size));tmpug=0.d0
    allocate(tmpur(nfft_div_size));tmpur=0.d0
    dFdrho=0.d0;dFddrho=0.d0

    irecv = myrank_ke+1
    if(irecv.ge.nrank_ke) irecv = irecv-nrank_ke
    isend = myrank_ke-1
    if(isend.lt.0)        isend = isend+nrank_ke

    do cqaa=1,mp_nq
       cqa = real_index(cqaa)
       if (nrank_ke>1)    theta_buf_r(1:np_nq,:) = theta_G_ab(1:np_nq,:)
       if (cqaa.le.np_nq) theta_G_a(:) = theta_G_ab(cqaa,:)

       if (.not.oneshot)  tmpug=(0.d0,0.d0)

       call tstatc0_begin('vdw_stress_core_core ',id_sname3,1)

       do i0=0,nrank_ke-1
          ipos = i0+myrank_ke
          if(ipos.ge.nrank_ke) ipos = ipos-nrank_ke
          if(ipos.lt.0)        ipos = ipos+nrank_ke

          if(cqaa.le.np_nq)then
             do i1=1,nel_nq(ipos)
                cqb = real_index(i1,ipos)
                if(nrank_ke>1)then
                   theta_G_b(:) = theta_buf_r(i1,:)
                else
                   theta_G_b(:) = theta_G_ab(i1,:)
                endif
                fac = 1.d0
                if(eval_kernel_by_interpolation) then
                  call phiab_by_interpl(cqa,cqb,nr12,phi_ab)
                else
                  Call phiab(cqa,cqb,nr12,phi_ab)
                endif
                Call calc_stress_contrib2( na, nb, nc, cqa, cqb, nr12, phi_ab, &
                     &                     theta_G_a, theta_G_b, fac, s_tmp, &
                     &                     tmpug )
                s_cnl2 = s_cnl2 +s_tmp
             enddo
          endif
          if(nrank_ke>1.and.i0.ne.nrank_ke-1)then
             call tstatc0_begin('vdWdf_stress_core (comm) ',id_sname2,1)
             allocate(theta_buf_s(1:mp_nq,nfft_div_size))
             theta_buf_s = theta_buf_r
             call mpi_sendrecv(theta_buf_s,mp_nq*nfft_div_size,mpi_double_complex,isend,0, &
                 &             theta_buf_r,mp_nq*nfft_div_size,mpi_double_complex,irecv,0, &
                 &             mpi_g_world,ista,ierr)
             deallocate(theta_buf_s)
             call tstatc0_end(id_sname2)
          endif
       enddo
       call tstatc0_end(id_sname3)
       if(.not.oneshot) call GtoR(na,nb,nc,tmpur,tmpug)
       if(.not.oneshot.and.cqaa.le.np_nq)then
          s_tmp = 0.0d0
          Do v=1, 3
             Do u=1, 3
                Do cir = 1,nfft_div_size
                   grho = grad(cir)
                   rtmp = rho(cir)
!                         if ( rtmp > 1.0d-6 .and. grho > 1.0d-6 ) then
                   if ( grho > 1.0d-6 ) then
                      s_tmp(u,v) = s_tmp(u,v) &
                           &      +dble(tmpur(cir)) &
                           &       *ddtheta_R(cqaa,cir) &
                           &       *cgrad(cir,u) &
                           &       *cgrad(cir,v) /grad(cir)
                   endif
                End Do
             End Do
          End Do
          s_cnl1(:,:) = s_cnl1(:,:) +s_tmp(:,:)
       endif
    enddo

    if(npes>1) then
       call mpi_allreduce( MPI_IN_PLACE, s_cnl1, 3*3, &
            &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
       call mpi_allreduce( MPI_IN_PLACE, s_cnl2, 3*3, &
            &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
       deallocate(tmpur)
    endif
    deallocate(tmpug)

    if(nrank_ke>1)then
       deallocate(ista);    deallocate(theta_buf_r)
    endif

    s_cnl1 = s_cnl1 *( univol *rinplw ) /univol
    s_cnl2 = s_cnl2 /2.0d0

    call tstatc0_end(id_sname)

  end subroutine vdWdf_stress_tensor_core
! ==== end for stress

  Subroutine convolution_3d_by_fft(na,nb,nc,cqa,cqb,nr12,phi_ab,theta_G_a,theta_G_b,fac,Ecnl_12_ab,tmpug)

    Implicit none

    !++++++++++++++++++++++++++ VARIABLES ++++++++++++++++++++++++++
    integer, intent(in) :: na,nb,nc,cqa,cqb,nr12
    complex(kind=DP), intent(in)  :: theta_G_a(nfft_div_size)
    complex(kind=DP), intent(in)  :: theta_G_b(nfft_div_size)
    real(kind=DP),    intent(in)  :: phi_ab(0:nr12)
    real(kind=DP),    intent(in)  :: fac
    real(kind=DP),    intent(out) :: Ecnl_12_ab
    complex(kind=DP), intent(inout) :: tmpug(nfft_div_size)
    complex(kind=DP) temp_c

    Real(kind=DP)  phix,phiy,Ta,Tb,Tc

    Integer cir,ck,i
    Integer :: c1,c2,c3,ca,cb,cc,cix,ciy,ciz
    Real(kind=DP) dk,rk,r12,term

    Real(kind=DP),allocatable :: core_G(:)

    real*8 :: dd

    real(kind=DP) :: pi2rk
    integer :: id_sname = -1

    !+++++++++++++++++++++ End VARIABLES +++++++++++++++++++++++++++
    call tstatc0_begin('convolution_3d_by_fft ',id_sname,1)
    allocate(core_G(0:nr12));core_G=0.d0

    !      dd = r12max/dble(nr12)

    call get_phi_ab_g(nr12,phi_ab,core_G)

    dk = 1.d0/r12max
    temp_c = (0.0d0,0.0d0)
    !tmpug = 0.d0

    if(oneshot)then
       do i = 1,nfft_div_size
          rk   = rkar(i)
          ck   = dint(rk/dk)

          phix = (rk - dk*DBLE(ck))/dk
          phiy = 1.0d0 -phix
          term = phiy*core_G(ck) + phix*core_G(ck+1)

          temp_c = temp_c + term *DCONJG(theta_G_a(i))   &
               &          *              theta_G_b(i)
       enddo
    else
       do i = 1,nfft_div_size
          rk   = rkar(i)
          ck   = dint(rk/dk)

          phix = (rk - dk*DBLE(ck))/dk
          phiy = 1.0d0 - phix
          term = phiy*core_G(ck) + phix*core_G(ck+1)

          temp_c = temp_c + term *DCONJG(theta_G_a(i))   &
               &           *             theta_G_b(i)

          tmpug(i) = tmpug(i) &
               &              +theta_G_b(i)*term *univol
       enddo
    endif

    temp_c = temp_c * univol**2
    Ecnl_12_ab = 0.5d0*DBLE(temp_c)*fac

    !***** End Calculate 'theta_G_a*core_G*theta_G_b'  ***********

    deallocate(core_G)

    call tstatc0_end(id_sname)

  End subroutine convolution_3d_by_fft
  !** End SUBROUTINE convolution_3d *************************************************************

!  Subroutine convolution_3d(na,nb,nc,cqa,cqb,nr12,phi_ab,theta_G_a,theta_G_b,Ecnl_12_ab)
!
!    Implicit none
!
!    !    include "fftw3.f"
!
!    !++++++++++++++++++++++++++ VARIABLES ++++++++++++++++++++++++++
!    integer, intent(in) :: na,nb,nc,cqa,cqb,nr12
!    complex(kind=DP), intent(in)  :: theta_G_a(-(na/2-1):na/2,-(nb/2-1):nb/2,-(nc/2-1):nc/2)
!    complex(kind=DP), intent(in)  :: theta_G_b(-(na/2-1):na/2,-(nb/2-1):nb/2,-(nc/2-1):nc/2)
!    real(kind=DP),    intent(in)  :: phi_ab(0:nr12)
!    real(kind=DP),    intent(out) :: Ecnl_12_ab
!    real(kind=DP) :: dr12
!    complex(kind=DP) temp_c
!
!    Integer  nabc
!    Real(kind=DP)  phix,phiy,Ta,Tb,Tc
!
!    Integer cix,ciy,ciz,cr,ck
!    Real(kind=DP) dk,rk,r12,term
!
!    Real(kind=DP),allocatable :: core_G(:)
!    Real(kind=DP),allocatable :: phiab_r(:)
!
!    real(kind=DP) :: pi2rk
!    !+++++++++++++++++++++ End VARIABLES +++++++++++++++++++++++++++
!
!    allocate(core_G(0:nk-1))
!    allocate(phiab_r(0:nr12));phiab_r=0.d0
!    nabc = na*nb*nc
!
!    !    dv = aa(1,1)*(aa(2,2)*aa(3,3)-aa(2,3)*aa(3,2)) &
!    !&      + aa(1,2)*(aa(2,3)*aa(3,1)-aa(2,1)*aa(3,3)) &
!    !&      + aa(1,3)*(aa(2,1)*aa(3,2)-aa(2,2)*aa(3,1))
!!
!    !    Ta = na*DSQRT(aa(1,1)**2 + aa(1,2)**2 + aa(1,3)**2)
!    !    Tb = nb*DSQRT(aa(2,1)**2 + aa(2,2)**2 + aa(2,3)**2)
!    !    Tc = nc*DSQRT(aa(3,1)**2 + aa(3,2)**2 + aa(3,3)**2)
!
!    !***** Make the core function and execute 3d-FFT by hand *****
!    dr12 = r12max/dble(nr12)
!    dk = maxk/dble(nk-1)
!
!    do cr=0,nr12
!       r12 = dr12*dble(cr)
!       phiab_r(cr) = phi_ab(cr)*r12
!    enddo
!
!    Do ck = 0,0
!       rk = dk * dble(ck)
!       term = 0
!       Do cr = 0,nr12
!          r12 = dr12*dble(cr)
!          term = term + phi_ab(cr) * (r12**2)
!       Enddo
!       core_G(ck) = 4.0d0*pi*dr12 * term
!    End do
!
!    Do ck = 1,nk-1
!       rk = dk * dble(ck)
!       pi2rk = 2.0*pi*rk
!       term = 0
!       Do cr = 0,nr12
!          r12 = dr12*dble(cr)
!          term = term + phiab_r(cr) * DSIN(pi2rk*r12)
!       Enddo
!       core_G(ck) = 2.d0 * dr12 * term / rk
!    End do
!    !***** END of Make the core function *************************
!
!    !    if(cqa==1)then
!    !       do ck=0,nr12
!    !          rk = dk * dble(ck)
!    !          write((cqa+1)*35+cqb,*) rk,core_G(ck)
!    !       enddo
!    !    endif
!    !***** Calculate 'theta_G_a*core_G*theta_G_b'  ***************
!    dk = maxk/dble(nk-1)
!    temp_c = (0.0,0.0)
!    Do cix = -(na/2-1),na/2
!       Do ciy = -(nb/2-1),nb/2
!          Do ciz = -(nc/2-1),nc/2
!
!             !rk = DSQRT((dble(cix)/Ta)**2 + (dble(ciy)/Tb)**2 + (dble(ciz)/Tc)**2)
!             rk = rkar(cix,ciy,ciz)
!             If(rk.LT.maxk) Then
!                ck = DINT(rk/dk)
!                phix = (rk - dk*DBLE(ck))/dk
!                phiy = 1.0 - phix
!                term = phiy*core_G(ck) + phix*core_G(ck+1)
!             Else
!                term = 0.0d0
!             Endif
!
!             temp_c = temp_c +                        &
!                  &             DCONJG(theta_G_a(cix,ciy,ciz)) *  &
!                  &                    theta_G_b(cix,ciy,ciz)  *  &
!                  &                  term
!          Enddo
!       Enddo
!    Enddo
!
!    temp_c = temp_c * dv*nabc
!    Ecnl_12_ab = 0.5d0*DBLE(temp_c)
!    !***** End Calculate 'theta_G_a*core_G*theta_G_b'  ***********
!
!    deallocate(core_G)
!    deallocate(phiab_r)
!
!
!  End subroutine convolution_3d
!  !** End SUBROUTINE convolution_3d *************************************************************

  Subroutine d_q0(n,dn,q0min,q0max,q0,dqdn,dqddn)
    Implicit None


    !************************ Note *********************************
    ! This program calculates q0.
    !
    ! Input
    !   rho(nrxyz,nsipn) : Total density
    !
    ! Output
    !   q0               :
    !
    !
    !
    !                            Written by Youky Ono
    !***************************************************************



    !++++++++++++++++++++++++++ VARIABLES ++++++++++++++++++++++++++
    real(kind=DP), intent(in)  :: n,dn,q0min,q0max
    real(kind=DP), intent(out) :: q0,dqdn,dqddn
    Double Precision rs,x,nn2,r,term1,term2,term3,rnn2

!    Double precision  kF,exc0,exLDA,excLDA,Zab,m,h,e,GxcLDA,eta,dq
!    Parameter (Zab=-0.8491d0)

    Double precision  kF,exc0,exLDA,excLDA,m,h,e,GxcLDA,eta,dq
    Parameter (e=1.0d0,m=1.0d0)

    integer :: i
    real(kind=DP) :: s2,s,qq
    real(kind=DP) :: drsdn,dxdn,dexcdn,dsdn,dkFdn
    !++++++++++++++++++++++++ end VARIABLES ++++++++++++++++++++++++



    !---------------------- Calculation Start ----------------------
    nn2 = dn*dn
    rnn2 = dsqrt(nn2)
    rs = (3.d0/(4.d0*pi*n))**(1.d0/3.d0)

    ! Eq.(58), (59) (p.93-94) of Theory of the Inhomogeneous Electron gas.
    ! S.Lundqvist and N.H.March 1983 Plenum Press, NY
    x = rs/11.4d0
    GxcLDA = 0.5d0*((1.d0+x**3)*DLog(1.d0+1.d0/x)-x**2+x/2.d0-1.d0/3.d0)
    excLDA = -0.458d0/rs-0.0666d0*GxcLDA

    kF = (3.d0*pi*pi*n)**(1.d0/3.d0)

    s2 = nn2/((2.d0*kF*n)**2);s=dsqrt(s2)
    qq = -(4.d0*pi/3.d0)*excLDA-(Zab/9.d0)*s2*kF

    call hxxc(qq,q0max,q0,dq)
    q0 = max(q0,q0min)

    if(q0.eq.q0min)then
       dqdn = 0.d0;dqddn=0.d0
    else
#if 0
       drsdn  = -(3.d0/(4.d0*PAI))**(1.d0/3.d0)/3.d0*n**(-4.d0/3.d0)
       dxdn   = drsdn/11.4d0
       dexcdn = (0.458d0/(rs*rs))*drsdn &
            &   -0.0333d0*(3*x*x*dlog(1+1.d0/x)-(1+x**3)/(x*(x+1))-2*x+0.5d0)*dxdn
       dkFdn  = (pi/kF)**2
!!!       dsdn   = -(rnn2/(4*kF*n*kF*n))*(dkFdn*n+kF)
       dsdn   = -(rnn2/(2*kF*n*kF*n))*(dkFdn*n+kF)
       dqdn   = -(4.d0*pi/3.d0)*dexcdn-(Zab/9.d0)*(2*s*dsdn*kF+s2*dkFdn)
       dqddn  = -(Zab/9.d0)*s/n

#else
       drsdn  = -rs /3.0d0 /n
       dxdn   = drsdn/11.4d0
       dexcdn = (0.458d0/(rs*rs))*drsdn &
            &   -0.0333d0*(3*x*x*dlog(1+1.d0/x)-(1+x**3)/(x*(x+1))-2*x+0.5d0)*dxdn
       dkFdn  = kF /3.0d0 /n
       dsdn   = -s /(kF*n) *(dkFdn*n+kF)
       dqdn   = -(4.d0*pi/3.d0)*dexcdn-(Zab/9.d0)*(2*s*dsdn*kF+s2*dkFdn)
!       dqddn  = -(Zab/9.d0)*s2 *2.0 /dn *KF
       dqddn  = -(Zab/9.d0)*s/n
#endif
       dqdn   = dqdn*dq
       dqddn  = dqddn*dq
    endif
  End Subroutine d_q0
  !** End SUBROUTINE d_q0 ***********************************************************************

  subroutine hxxc(x,xc,hx,dhx)
    real(kind=DP), intent(in)  :: x,xc
    real(kind=DP), intent(out) :: hx,dhx
    integer :: i
    integer :: mc = 12
    real(kind=DP) :: summ,dsumm
    hx = 0.d0;dhx = 0.d0
    summ=0.d0;dsumm=0.d0
    do i=1,mc
       summ = summ+(x/xc)**i/dble(i)
       dsumm = dsumm+x**(i-1)/xc**i
    enddo
    hx = xc*(1.d0-dexp(-summ))
    dhx = xc*dexp(-summ)*dsumm
  end subroutine hxxc

  !** SUBROUTINE derivation *********************************************************************
  Subroutine derivation(na,nb,nc,aa,rho,dv,grad)
    Implicit none

    !++++++++++++++++++++++++++ VARIABLES ++++++++++++++++++++++++++

    ! The unit cell and the electron density information
    integer,intent(in) :: na,nb,nc
    real(kind=DP), intent(in) :: rho(nabc)
    Real(kind=DP), intent(in) :: aa(3,3)
    real(kind=DP), intent(in) :: dv
    real(kind=DP), intent(out) :: grad(nabc)

    ! Integers
    Integer i,j,k,cx,cy,cz,cr

    ! Internal valuables
    Integer zx(-3:3),zy(-3:3),zz(-3:3)
    Real(kind=DP) rn(3,-3:3),detr,bb(3,3)
    Real(kind=DP),allocatable ::  darho(:),dbrho(:),dcrho(:)
    real(kind=DP) :: dx,dy,dz
    !+++++++++++++++++++++ End VARIABLES +++++++++++++++++++++++++++

    allocate(darho(nabc));darho=0.d0
    allocate(dbrho(nabc));dbrho=0.d0
    allocate(dcrho(nabc));dcrho=0.d0
    !    dv = aa(1,1)*(aa(2,2)*aa(3,3)-aa(2,3)*aa(3,2)) &
    !&      + aa(1,2)*(aa(2,3)*aa(3,1)-aa(2,1)*aa(3,3)) &
    !&      + aa(1,3)*(aa(2,1)*aa(3,2)-aa(2,2)*aa(3,1))


    Do cx = 1,na
       Do cy = 1,nb
          Do cz = 1,nc
             cr = get_index_1D(cx,cy,cz)
             Do j = -3,3
                zx(j) = MOD(2*na+(cx+j)-1,na)+1
                zy(j) = MOD(2*nb+(cy+j)-1,nb)+1
                zz(j) = MOD(2*nc+(cz+j)-1,nc)+1
                rn(1,j) = rho(get_index_1D(zx(j),cy,cz))
                rn(2,j) = rho(get_index_1D(cx,zy(j),cz))
                rn(3,j) = rho(get_index_1D(cx,cy,zz(j)))
             Enddo
             darho(cr) = (rn(1,3)-9.d0*rn(1,2)+45.d0*rn(1,1)-45.d0*rn(1,-1)+9.d0*rn(1,-2)-rn(1,-3))/(60.d0)
             dbrho(cr) = (rn(2,3)-9.d0*rn(2,2)+45.d0*rn(2,1)-45.d0*rn(2,-1)+9.d0*rn(2,-2)-rn(2,-3))/(60.d0)
             dcrho(cr) = (rn(3,3)-9.d0*rn(3,2)+45.d0*rn(3,1)-45.d0*rn(3,-1)+9.d0*rn(3,-2)-rn(3,-3))/(60.d0)
          Enddo
       Enddo
    Enddo

    detr = (aa(1,1)*aa(2,2)*aa(3,3)+aa(1,2)*aa(2,3)*aa(3,1)+aa(1,3)*aa(2,1)*aa(3,2)) &
         &- (aa(1,1)*aa(2,3)*aa(3,2)+aa(1,2)*aa(2,1)*aa(3,3)+aa(1,3)*aa(2,2)*aa(3,1))

    bb(1,1) =  (aa(2,2)*aa(3,3)-aa(2,3)*aa(3,2))/detr
    bb(2,1) = -(aa(2,1)*aa(3,3)-aa(2,3)*aa(3,1))/detr
    bb(3,1) =  (aa(2,1)*aa(3,2)-aa(2,2)*aa(3,1))/detr
    bb(1,2) = -(aa(1,2)*aa(3,3)-aa(1,3)*aa(3,2))/detr
    bb(2,2) =  (aa(1,1)*aa(3,3)-aa(1,3)*aa(3,1))/detr
    bb(3,2) = -(aa(1,1)*aa(3,2)-aa(1,2)*aa(3,1))/detr
    bb(1,3) =  (aa(1,2)*aa(2,3)-aa(1,3)*aa(2,2))/detr
    bb(2,3) = -(aa(1,1)*aa(2,3)-aa(1,3)*aa(2,1))/detr
    bb(3,3) =  (aa(1,1)*aa(2,2)-aa(1,2)*aa(2,1))/detr

    Do cr = 1,nabc
       dx = (bb(1,1)*darho(cr) + bb(1,2)*dbrho(cr) + bb(1,3)*dcrho(cr))
       dy = (bb(2,1)*darho(cr) + bb(2,2)*dbrho(cr) + bb(2,3)*dcrho(cr))
       dz = (bb(3,1)*darho(cr) + bb(3,2)*dbrho(cr) + bb(3,3)*dcrho(cr))
       grad(cr) = dsqrt(dx**2+dy**2+dz**2)
    Enddo

!    deallocate(darho)
!    deallocate(dbrho)
!    deallocate(dcrho)

  End subroutine derivation
  !** End SUBROUTINE derivation *****************************************************************

  subroutine spline0(nq0,x,y2)
    integer, intent(in) :: nq0
    real(kind=DP), intent(in),  dimension(nq0) :: x
    real(kind=DP), intent(out), dimension(nq0,nq0) :: y2
    real(kind=DP), allocatable,dimension(:) :: y,y2tmp
    integer :: i
    allocate(y(nq0));y=0.d0
    allocate(y2tmp(nq0));y2tmp=0.d0
    y2 = 0.d0
    do i=1,nq0
       y(:) = 0.d0;y(i) = 1.d0
       call init_cubic_spline(nq0,x,y,y2tmp)
       y2(i,:) = y2tmp(:)
    enddo
    deallocate(y)
    deallocate(y2tmp)
  end subroutine spline0

  function hofy(num,denom) result(res)
    Implicit none
    real(kind=DP), intent(in) :: num,denom
    real(kind=DP) :: res
    real(kind=DP) :: gam=1.396263402d0 ! 4pi/9
    if(abs(denom)<epsilon(denom))then
       res = 1.d0
    else
       res = 1.d0-dexp(-gam*((num/denom)**2))
    endif
  end function hofy

  function Wab(a,b) result (res)
    Implicit none
    real(kind=DP), intent(in) :: a,b
    real(kind=DP) :: res
    res = 2.d0*((3.d0-a*a)*b*DCOS(b)*DSIN(a) + (3.d0-b*b)*a*DCOS(a)*DSIN(b) &
         &       + (a*a+b*b-3.d0)*DSIN(a)*DSIN(b) &
         &       - 3.d0*a*b*DCOS(a)*DCOS(b))/((a*b)**3)
  end function Wab

  function Twxyz(w,x,y,z) result (res)
    Implicit none
    real(kind=DP),intent(in) :: w,x,y,z
    real(kind=DP) :: res
    res = 0.5d0*(1.d0/(w+x)+1.d0/(y+z))*(1.d0/((w+y)*(x+z))+1.d0/((w+z)*(y+x)))
  end function Twxyz

  !** SUBROUTINE kernel_phi *********************************************************************
  Subroutine kernel_phi(di,dk,phi)
    implicit none

    !************************ Note *********************************
    ! This Algorism follows Dion's 1-shot method.
    !
    ! This program is a subroutine.
    ! This program calculates the kernel function phi.
    !
    ! Input
    !
    !
    ! Output
    !   phi              :
    !
    !
    !
    !                            Written by Youky Ono in 2013/Jan.
    !***************************************************************



    !++++++++++++++++++++++++++ VARIABLES ++++++++++++++++++++++++++
    real(kind=DP), intent(in)  :: di,dk
    real(kind=DP), intent(out) :: phi
    Integer ca,cb,nb
    Real(kind=DP) a,b,h,v1,v2,v3,v4,W,T
    Real(kind=DP) e,m
    parameter (e=1.d0,m=1.d0)
    real(kind=DP) :: dr


    ! Gauss-Legendre integration
    Real(kind=DP),allocatable, dimension(:) :: xi,wi
    Real(kind=DP) :: fac
    !++++++++++++++++++++++++ end VARIABLES ++++++++++++++++++++++++

    !---------------------- Calculation Start ----------------------
    ! Call gauleg for Gauss-Legendre integral
    allocate(xi(na_gl));xi=0.d0
    allocate(wi(na_gl));wi=0.d0
    Call gauleg(a1,a2,na_gl,xi,wi)

    phi = 0.d0
    Do ca=1,na_gl
       Do cb=1,ca
          fac=2.d0
          if(ca.eq.cb) fac=1.d0
          a = xi(ca)
          b = xi(cb)

          v1 = (a**2)/(2.d0*hofy(a,di))
          v2 = (b**2)/(2.d0*hofy(b,di))
          v3 = (a**2)/(2.d0*hofy(a,dk))
          v4 = (b**2)/(2.d0*hofy(b,dk))
          phi = phi + fac*wi(ca)*wi(cb)*(a*a*b*b)*Wab(a,b)*Twxyz(v1,v2,v3,v4)
       End Do
    End Do

    phi = phi * 2.d0 *m *(e**4)/(pi**2)
    deallocate(xi)
    deallocate(wi)
  end Subroutine kernel_phi
  !** End SUBROUTINE kernel_phi *****************************************************************

  !** SUBROUTINE gauleg **********************************************************************************
  Subroutine gauleg(x1,x2,n,xi,wi)
    Implicit none
    real(kind=DP), intent(in) :: x1,x2
    integer, intent(in) :: n
    real(kind=DP), intent(out)  ::  xi(n),wi(n)

    !++++++++++++++++++++++++++ VARIABLES ++++++++++++++++++++++++++
    Integer  m,j,i
    REAL(kind=DP)  z1,z,xm,xl,pp,p3,p2,p1,eta
    Parameter (eta=0.0000000001d0)

    !++++++++++++++++++++++++ end VARIABLES ++++++++++++++++++++++++



    !---------------------- Calculation Start ----------------------
    m=(n+1)/2
    xm=0.5d0*(x2+x1)
    xl=0.5d0*(x2-x1)

    ! === DEBUG by tkato 2014/04/22 ================================================
    z1 = 0.0d0
    ! ==============================================================================
    Do i=1,m
       z=DCOS(pi*(i-0.25d0)/(n+0.5d0))
       Do While (ABS(z-z1).GT.eta)
          p1=1.0d0
          p2=0.0d0

          Do j=1,n
             p3=p2
             P2=p1
             p1=((2.0d0*j-1.d0)*z*p2-(j-1.d0)*p3)/dble(j)
          Enddo

          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp

       Enddo
       xi(i) = xm-xl*z
       xi(n+1-i) = xm+xl*z
       wi(i) = 2.d0*xl/((1.d0-z*z)*pp*pp)
       wi(n+1-i) = wi(i)
    Enddo
  End Subroutine gauleg
  !** End SUBROUTINE gauleg ******************************************************************************

  !** SUBROUTINE outputs ************************************************************************
  Subroutine outputs()
    Implicit none

    Ecnl = Ecnl_12 + Ecnl_3 - Ecnl_3s

    if(printable)then
       write(nfout,'(a)')  'Here are the results : '
       Write(nfout,*)  ' '
       Write(nfout,11) ExGGA
       Write(nfout,*)  ' '
       Write(nfout,12) EcLDA
       Write(nfout,13) Ecnl
       Write(nfout,14) EcLDA + Ecnl
       Write(nfout,*)  ' '
       Write(nfout,15) EcLDA + Ecnl + ExGGA
       Write(nfout,*)  ' '
11     Format('E_total(GGA exchange)      = ',F19.13)

12     Format('Ec(LDA)                    = ',F19.13)
13     Format('Ec(nl)                     = ',F19.13)
14     Format('Ec (= Ec(LDA) + Ec(nl) )   = ',F19.13)

15     Format('E_total(vdW-DF)            = ',F19.13)

       Write(nfout,*)  '                  Given in Hartree atomic units'
       Write(nfout,*)  ' '
    endif

  End Subroutine outputs
  !** End SUBROUTINE outputs ********************************************************************

  subroutine get_dFdrho_dFddrho(dFdrho_,dFddrho_)
    real(kind=DP), dimension(nfft_div_size), intent(out) :: dFdrho_,dFddrho_
    integer :: i,i1,i2,i3
    dFdrho_(:) = dFdrho(:)
    dFddrho_(:) = dFddrho(:)
    !    dFdrho = 0.0d0;dFddrho = 0.0d0 ! this can be calculated on the fly!!
    !    do i=1,nq0
    !       dFdrho (:,:,:) = dFdrho (:,:,:) + dble(ualpha_r(i,:,:,:))*dtheta_R (i,:,:,:)/(univol*rinplw)
    !       dFddrho(:,:,:) = dFddrho(:,:,:) + dble(ualpha_r(i,:,:,:))*ddtheta_R(i,:,:,:)/(univol*rinplw)
    !    enddo

  end subroutine get_dFdrho_dFddrho

  subroutine map_rho_to_afft(rho,afft)
    complex(kind=DP), intent(in),  dimension(nabc) :: rho
    real(kind=DP), intent(out), dimension(nfftp_nonpara) :: afft
    integer :: i,cix,ciy,ciz,cir,cir2
    afft = 0.d0
    do cix = 1,na
       do ciy = 1,nb
          do ciz = 1,nc
             cir = get_index_1D_kimg(cix,ciy,ciz)
             cir2 = get_index_1D(cix,ciy,ciz)
             if(kimg == 1)then
                afft(2*cir-1) = real(rho(cir2))
             else
                afft(2*cir-1) = real(rho(cir2))
                afft(2*cir)   = aimag(rho(cir2))
             endif
          enddo
       enddo
    enddo
  end subroutine map_rho_to_afft

  subroutine map_rhog_to_afft(rho,afft)
    complex(kind=DP), intent(in),  dimension(nabc) :: rho
    real(kind=DP), intent(out), dimension(nfftp_nonpara) :: afft
    integer :: i,cix,ciy,ciz,cir,cir2
    afft = 0.d0
    do cix = 1,na
       do ciy = 1,nb
          do ciz = 1,nc
             cir = get_index_1D_kimg(cix,ciy,ciz,.true.)
             cir2 = get_index_1D(cix,ciy,ciz)
             if(kimg == 1)then
                afft(cir) = real(rho(cir2))
             else
                afft(2*cir-1) = real(rho(cir2))
                afft(2*cir)   = aimag(rho(cir2))
             endif
          enddo
       enddo
    enddo
  end subroutine map_rhog_to_afft

  subroutine map_afft_to_rho(afft,rho)
    real(kind=DP), intent(in), dimension(nfftp_nonpara) :: afft
    complex(kind=DP), intent(out),  dimension(nabc) :: rho
    integer :: i,cix,ciy,ciz,cir,cir2
    integer :: m
    rho = (0.d0,0.d0)

    do cix = 1,na
       do ciy = 1,nb
          do ciz = 1,nc
             cir = get_index_1D_kimg(cix,ciy,ciz,.true.)
             cir2 = get_index_1D(cix,ciy,ciz)
             if(kimg == 1)then
                rho(cir2) = DCMPLX(afft(cir))
             else
                rho(cir2) = DCMPLX(afft(2*cir-1),afft(2*cir))
             endif
          enddo
       enddo
    enddo
  end subroutine map_afft_to_rho

  subroutine map_afft_to_rhog(afft,rho)
    real(kind=DP), intent(in), dimension(nfftp_nonpara) :: afft
    complex(kind=DP), intent(out),  dimension(nabc) :: rho
    integer :: i,cix,ciy,ciz,cir,cir2
    rho = (0.d0,0.d0)
    do cix = 1,na
       do ciy = 1,nb
          do ciz = 1,nc
             cir = get_index_1D_kimg(cix,ciy,ciz)
             cir2 = get_index_1D(cix,ciy,ciz)
             if(kimg == 1)then
                rho(cir2) = DCMPLX(afft(2*cir-1))
             else
                rho(cir2) = DCMPLX(afft(2*cir-1),afft(2*cir))
             endif
          enddo
       enddo
    enddo
  end subroutine map_afft_to_rhog

end module m_vdWDF
#endif
