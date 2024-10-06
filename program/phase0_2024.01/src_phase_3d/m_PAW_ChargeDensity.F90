!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 633 $)
!
!  MODULE: m_PAW_ChargeDensity
!
!  AUTHOR(S): T. Yamasaki, T. Yamamoto and T. Ohno November/2009
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
module m_PAW_ChargeDensity
  use m_db,                   only : getIntDB,getStringDB_TB,getIntDB_TB
  use m_Const_Parameters,     only : DP,PAI2,PAI4,BUCS,SphericalHarmonicsExpansion, GaussLegendre, LOWER, Bohr &
       &                           , FMAXVALLEN,NOCONV, ON, UP, OFF
  use m_Control_Parameters,   only : kimg,nspin,af,ipripaw,ipriinputfile,printable
  use m_Files,                only : nfout
  use m_Ionic_System,         only : ityp,natm,ntyp,pos,iwei, speciesname, amion, surface_integral_method_paw
  use m_Charge_Density,       only : hsr,chgq_l,chgsoft
  use m_PlaneWaveBasisSet,    only : igfp_l
  use m_PseudoPotential,      only : psirpw,phirpw,qrspspw &
       &                            ,rhcorpw,rhpcrpw &
       &                            ,ilmt,ltp,mtp,taup,radr_paw &
       &                            ,il2p,isph,dl2p,iqitg &
       &                            ,m_PP_find_maximum_l &
       &                            ,ipaw,wf_mnrc,flg_symmtry &
       &                            ,mmesh, nmesh
  use m_FFT,                  only : fft_box_size_CD,fft_box_size_CD_c, nfftp &
       , m_FFT_CD_inverse0 &
       , m_FFT_check_of_negative_CD &
       , m_FFT_cp_afft_CD &
       , m_FFT_alloc_CD_box &
       , m_FFT_dealloc_CD_box
  use m_Parallelization,      only : ista_kngp,iend_kngp,npes,mype,ierr &
       , ista_fftp,iend_fftp &
       , nis_fftp,nie_fftp,mp_fftp &
       , npes_cdfft, ista_fftph, iend_fftph &
       , MPI_CommGroup
  use m_Crystal_Structure,    only : altv,rltv,nopr,m_CS_op_in_PUCD,tau,op
  use m_PAW_Tecplot,          only : unusedUnit,output_box_tecplot
  use m_PlaneWaveBasisSet,    only : ngabc,kgp, kgp_reduced
  use m_Timing,               only : tstatc0_begin, tstatc0_end

! ============================ added by K. Tagami =============== 11.0
  use m_Control_Parameters,    only : ndim_magmom, noncol
  use m_Parallelization,      only : ista_atm, iend_atm
  use m_ES_NonCollinear,     only :  Global_Quantz_Axis_now
! =============================================================== 11.0

  use m_Control_Parameters,  only : sw_opencore, sw_xc_opencore_ae_only
  use m_PS_opencore,  only : rmag_opencore, has_opencore, mag_opencore_pol
  use mpi

  implicit none
!  include 'mpif.h'

  private

  integer:: ntheta,nphi
  real(DP),pointer,dimension(:)                 :: cos_theta
  real(DP),pointer,dimension(:)                 :: omg_wght
  real(DP),pointer,dimension(:,:,:)             :: ylm,dylm_dth,dylm_dph
  integer,pointer,dimension(:,:)               :: ia2ia_symmtry_op
  real(kind=DP), allocatable,dimension(:,:,:)  :: op_pr
!!$    character(len=64),pointer,dimension(:)       :: surface_integral_method_char
  integer,allocatable,dimension(:)             :: surface_integral_method

  real(kind=DP), pointer, dimension(:,:,:)      :: paw_cr2
  integer, pointer, dimension(:,:,:)           :: paw_isph2
  integer, pointer, dimension(:,:)             :: paw_mmt2
  integer, pointer, dimension(:)               :: paw_dnr

  character(len("paw_one_center_integral")),private,parameter :: tag_paw_integral="paw_one_center_integral"
  character(len("element_list")),private,parameter :: tag_element_list="element_list"
  character(len("element")),private,parameter :: tag_element = "element"
  character(len("surface_integral_method")),private,parameter :: tag_surface_integral_method="surface_integral_method"
  character(len("dnr")),private,parameter ::  tag_dnr = "dnr"
  character(len("SphericalHarmonicsExpansion")),private,parameter :: tag_Sphex ="sphericalharmonicsexpansion"
  character(len("SphericalHarmonics"))         ,private,parameter :: tag_Sphex2="sphericalharmonics"
  character(len("Sphex"))                      ,private,parameter :: tag_Sphex3="sphex"
  character(len("GaussLegendre"))              ,private,parameter :: tag_GL ="gausslegendre"
  character(len("GL"))                         ,private,parameter :: tag_GL2="gl"

  logical:: calcGaussLegendreIntegration = .false.
  logical:: calcSphericalHarmonicsExpansion = .false.

  public:: ntheta,nphi,omg_wght
  public:: m_PAWCD_rd_ntheta_nphi
  public:: m_PAWCD_init_surf_integration
  public:: m_PAWCD_set_ae_cd
  public:: m_PAWCD_set_ae_cd_sym
  public:: m_PAWCD_set_ae_der_cd
  public:: m_PAWCD_set_ae_der_cd_sym
  public:: m_PAWCD_set_ps_cd
  public:: m_PAWCD_set_ps_cd_sym
  public:: m_PAWCD_set_ps_der_cd
  public:: m_PAWCD_set_ps_der_cd_sym
  public:: ylm,dylm_dth,dylm_dph
  public:: cos_theta
  public:: m_PAWCD_dealloc
  public:: m_PAWCD_wd_cd

  public:: m_PAWCD_wd_number_of_state
  public:: surface_integral_method


! ======================== added by K. Tagami ================= 11.0
  integer max_sph_expansion
  parameter( max_sph_expansion = 25 )

  integer, parameter :: project_on_local_quantz_axis = 10
  integer, parameter :: project_on_global_quantz_axis = 15
!
  integer :: mode_of_projection = project_on_local_quantz_axis
!
  public :: m_PAWCD_ae_cd_sphex2_nonclA
  public :: m_PAWCD_ps_cd_sphex2_nonclA
  public :: max_sph_expansion
! ============================================================== 11.0

!!$    public:: m_PAWCD_set_ae_cd_sphex
  public:: m_PAWCD_set_ae_cd_sphex2
  public:: m_PAWCD_set_ae_cd_sphex2_2D
!!$    public:: m_PAWCD_set_ae_cd_sphex3
!!$    public:: m_PAWCD_set_ps_cd_sphex
  public:: m_PAWCD_set_ps_cd_sphex2
  public:: m_PAWCD_set_ps_cd_sphex2_2D
!!$    public:: m_PAWCD_set_ps_cd_sphex3
!!$    public:: m_PAWCD_set_sq_der_cd_sdphex
  public:: m_PAWCD_set_sq_der_cd_sdphex2
!!$    public:: m_PAWCD_set_sq_der_cd_sphex3
  public:: m_PAWCD_set_cr2_isph2_mmt2

  public:: m_PAW_dealloc
  public:: paw_cr2,paw_isph2,paw_mmt2,paw_dnr
  public:: calcGaussLegendreIntegration
  public:: calcSphericalHarmonicsExpansion


contains

  subroutine m_PAW_dealloc()
    call dealloc_integral_method_etc()
  end subroutine m_PAW_dealloc

  subroutine alloc_integral_method_etc
    allocate(surface_integral_method(ntyp));  surface_integral_method = SphericalHarmonicsExpansion
    allocate(paw_dnr(ntyp)); paw_dnr = 1
  end subroutine alloc_integral_method_etc

  subroutine dealloc_integral_method_etc()
    if(allocated(surface_integral_method)) deallocate(surface_integral_method)
    deallocate(paw_dnr)
  end subroutine dealloc_integral_method_etc

  subroutine m_PAWCD_rd_ntheta_nphi
    call getIntDB("paw surface_integration ntheta",ntheta,10)
    call getIntDB("paw surface_integration nphi",nphi,10)
    call alloc_integral_method_etc()
    surface_integral_method(:) = surface_integral_method_paw(:)
    deallocate(surface_integral_method_paw)
    call rd_surface_integral_method
    return
  end subroutine m_PAWCD_rd_ntheta_nphi

  subroutine rd_surface_integral_method
    integer :: i,j, npt, itr
    integer :: f_selectBlock, f_selectParentBlock
    integer :: f_selectFirstTableLine, f_selectNextTableLine, f_getStringValue, f_getIntValue
    logical :: found, tf
    character(len=FMAXVALLEN) :: rstr

    if(f_selectBlock(tag_paw_integral) == 0) then

       if(ipriinputfile >=2)  write(nfout,'(" !** -- tag_paw_one_center_integral --")')
       if( f_selectBlock(tag_element_list) == 0) then
          i = 1
          do while(.true.)
             if(i==1) then
                if( f_selectFirstTableLine() /= 0) then
                   exit
                end if
             else
                if( f_selectNextTableLine() /= 0) then
                   exit
                end if
             end if
             if(i > ntyp) exit
             found = .false.
             if(f_getStringValue(tag_element,rstr,NOCONV) == 0) then
                do j = 1, ntyp
                   if(trim(rstr).eq.trim(speciesname(j))) then
                      found = .true.
                      npt = j
                      exit
                   end if
                end do
             end if
             if(found) then
                if(npt >=1 .and. npt <=ntyp) then
                   if(f_getStringValue(tag_surface_integral_method,rstr,LOWER) == 0) then
                      call strncmp0(tag_Sphex,trim(rstr),tf)
                      if(.not.tf) call strncmp0(tag_Sphex2,trim(rstr),tf)
                      if(.not.tf) call strncmp0(tag_Sphex3,trim(rstr),tf)
                      if(tf) surface_integral_method(npt) = SphericalHarmonicsExpansion
                      call strncmp0(tag_GL,trim(rstr),tf)
                      if(.not.tf) call strncmp0(tag_GL2,trim(rstr),tf)
                      if(tf) surface_integral_method(npt) = GaussLegendre
                   end if
                   if(f_getIntValue(tag_dnr,itr) == 0) then
                      paw_dnr(npt) = itr
                   end if
                end if
             end if
             i = i+1
          end do
          itr = f_selectParentBlock()
       end if
       itr = f_selectParentBlock()
    end if
    if(ipriinputfile>=1) then
       write(nfout,*)
       write(nfout,*) '<***   PAW one center integral   ***>'
       do i=1,ntyp
          if(surface_integral_method(i) == SphericalHarmonicsExpansion) then
             write(nfout,'(i2,a,a,a,a,a,i2)') i,'-th element : ',trim(speciesname(i)) &
                  & ,' : ',trim(tag_Sphex),' : dnr = ',paw_dnr(i)
          else if(surface_integral_method(i) == GaussLegendre) then
             write(nfout,'(i2,a,a,a,a,a,i2)') i,'-th element : ',trim(speciesname(i)) &
                  & ,' : ',trim(tag_GL),' : dnr = ',paw_dnr(i)
          end if
          if(surface_integral_method(i) .ne. SphericalHarmonicsExpansion .and. &
               surface_integral_method(i) .ne. GaussLegendre) then
             write(nfout,*) 'Error in rd_surface_integral_method !'
             call phase_error_with_msg(nfout,'Error in rd_surface_integral_method !',__LINE__,__FILE__)
          end if
       end do
    end if

  end subroutine rd_surface_integral_method

  subroutine m_PAWCD_init_surf_integration
    logical:: doGaussLegendre
    doGaussLegendre=.false.
    call setDoGaussLegendre(doGaussLegendre)
    calcGaussLegendreIntegration=doGaussLegendre
    call setDoSphHarmonicsExpansion(calcSphericalHarmonicsExpansion)
    !        if(.not.doGaussLegendre) return
    allocate(cos_theta(ntheta));cos_theta=0
    allocate(omg_wght(ntheta));omg_wght=0
    call GaussLeg(-1.d0,1.d0,ntheta,cos_theta,omg_wght)
    omg_wght=omg_wght/2.d0/dble(nphi)
    call m_PAWCD_set_ylms(25)
    allocate(ia2ia_symmtry_op(natm,nopr+af))
    ia2ia_symmtry_op=0
    call set_ia2ia_symmtry_op
    return
  end subroutine m_PAWCD_init_surf_integration

  subroutine setDoGaussLegendre(flg)
    logical,intent(out):: flg
    integer:: i

    flg=.false.
    do i=1,ntyp
       if(ipaw(i)/=1) cycle
!!$            if(surface_integral_method(i) .eq. "GaussLegendre") then
       if(surface_integral_method(i) .eq. GaussLegendre) then
          flg=.true.
          exit
       end if
    end do
    return
  end subroutine setDoGaussLegendre

  subroutine setDoSphHarmonicsExpansion(flg)
    logical,intent(out):: flg
    integer:: i

    flg=.false.
    do i=1,ntyp
       if(ipaw(i)/=1) cycle
!!$            if(surface_integral_method(i) .eq. "SphericalHarmonicsExpansion") then
       if(surface_integral_method(i) .eq. SphericalHarmonicsExpansion) then
          flg=.true.
          exit
       end if
    end do
    return
  end subroutine setDoSphHarmonicsExpansion

  subroutine m_PAWCD_set_ylms(ilmax)
    integer,intent(in):: ilmax

    real(DP):: cos_the,sin_the
    real(DP):: phi,rr
    integer:: it,ip,itp,is
    real(DP),pointer,dimension(:)     :: gx,gy,gz,tylm
    real(DP),pointer,dimension(:,:)   :: tdylm

    allocate(gx(ntheta*nphi))
    allocate(gy(ntheta*nphi))
    allocate(gz(ntheta*nphi))
    allocate(tylm(ntheta*nphi))
    allocate(tdylm(ntheta*nphi,3))

    allocate(ylm(ntheta,nphi,ilmax))
    allocate(dylm_dth(ntheta,nphi,ilmax))
    allocate(dylm_dph(ntheta,nphi,ilmax))

    itp=0
    do it=1,ntheta
       cos_the=cos_theta(it)
       sin_the=sqrt(1.d0-cos_the*cos_the)
       do ip=1,nphi
          itp=itp+1
          phi=PAI2/dble(nphi)*dble(ip-1)
          gx(itp)=sin_the*cos(phi)
          gy(itp)=sin_the*sin(phi)
          gz(itp)=cos_the
       end do
    end do

    do is=1,ilmax
       call sphr(ntheta*nphi,is,gx,gy,gz,tylm)
       itp=0
       do it=1,ntheta
          do ip=1,nphi
             itp=itp+1
             ylm(it,ip,is)=tylm(itp)
          end do
       end do
    end do

    do is=1,ilmax
       call sphr_diff(ntheta*nphi,ntheta*nphi,is,gx,gy,gz,tdylm)
       itp=0
       do it=1,ntheta
          do ip=1,nphi
             itp=itp+1
!                    ylm(it,ip,is)=tylm(itp)
             rr=sqrt(gx(itp)**2+gy(itp)**2)
             dylm_dth(it,ip,is)= gz(itp)*gx(itp)/rr*tdylm(itp,1)+ &
                  gy(itp)*gz(itp)/rr*tdylm(itp,2)- &
                  rr*tdylm(itp,3)
             dylm_dph(it,ip,is)=-gy(itp)*tdylm(itp,1)+ &
                  gx(itp)*tdylm(itp,2)
          end do
       end do
    end do

    deallocate(gx,gy,gz,tylm,tdylm)
    return
  end subroutine m_PAWCD_set_ylms

  subroutine m_PAWCD_dealloc
! ======================= ASMS_DEBUG ================= 2013/02/07
    if (.not. associated(cos_theta)) return
                               !   sw_phonon == on .and. sw_calc_force == off
! ======================= ASMS_DEBUG ================= 2013/02/07
    deallocate(cos_theta)
    deallocate(omg_wght)
    deallocate(ylm,dylm_dth,dylm_dph)
    deallocate(ia2ia_symmtry_op)
    return
  end subroutine m_PAWCD_dealloc

  subroutine m_PAWCD_set_ae_cd(ith,iph,ia,nspin,nrc,n1nc,ista_nrc,iend_nrc)
    integer,intent(in) :: ith,iph,ia,nspin,nrc
    real(DP),intent(out):: n1nc(ista_nrc:iend_nrc,nspin)
    integer, intent(in) :: ista_nrc,iend_nrc

    real(DP):: n1(nrc,nspin),sum,fac, weight
    integer:: ir,it,is,lmt1,lmt2
    integer:: il1,il2,it1,it2,im1,im2
    integer:: ii,jj

    it=ityp(ia)
    n1nc=0.d0
    do lmt1=1,ilmt(it)
       il1=ltp(lmt1,it)
       im1=mtp(lmt1,it)
       it1=taup(lmt1,it)
       ii=(il1-1)**2+im1
       do lmt2=lmt1,ilmt(it)
          il2=ltp(lmt2,it)
          im2=mtp(lmt2,it)
          it2=taup(lmt2,it)
          jj=(il2-1)**2+im2
          fac=2.d0;if(lmt1.eq.lmt2) fac=1.d0

          do is=1,nspin
             sum=fac*hsr(ia,lmt1,lmt2,is)*ylm(ith,iph,ii)*ylm(ith,iph,jj)
!!$             do ir=1,nrc
             do ir=ista_nrc,iend_nrc
                n1nc(ir,is)=n1nc(ir,is)+sum*psirpw(ir,il1,it1,it)*psirpw(ir,il2,it2,it)
             end do
          end do
       end do
    end do

    do is=1,nspin
!!$       do ir=1,nrc
       do ir=ista_nrc,iend_nrc
          n1nc(ir,is)=(n1nc(ir,is)+rhcorpw(ir,it)/PAI4/dble(nspin)) &
               /radr_paw(ir,it)/radr_paw(ir,it)
       end do
    end do
    if ( sw_opencore == ON .and. nspin == 2 )then
       if ( has_opencore(it) == 1 ) then
          weight = 0.5d0 *mag_opencore_pol(ia,1)
          do ir=ista_nrc, iend_nrc
             n1nc(ir,1) = n1nc(ir,1) +rmag_opencore(ir,it)/PAI4 *weight &
                  &                        /radr_paw(ir,it)/radr_paw(ir,it)
             n1nc(ir,2) = n1nc(ir,2) -rmag_opencore(ir,it)/PAI4 *weight &
                  &                        /radr_paw(ir,it)/radr_paw(ir,it)
          end do
       endif
    endif

    return
  end subroutine m_PAWCD_set_ae_cd

  subroutine m_PAWCD_set_ae_cd_sym(vec0,ia,nspin,nrc,n1nc,ista_nrc,iend_nrc)
    integer,intent(in) :: ia,nspin,nrc
    real(DP),intent(in) :: vec0(3)
    real(DP),intent(out):: n1nc(ista_nrc:iend_nrc,nspin)
    integer, intent(in) :: ista_nrc,iend_nrc

    real(DP):: sum,fac,sum2, weight
    integer:: ir,it,is,lmt1,lmt2
    integer:: il1,il2,it1,it2,im1,im2
    integer:: ii,jj
    integer:: ja,iopr
    real(DP),pointer,dimension(:,:):: ylm_sym
    integer :: id_sname = -1
    call tstatc0_begin('m_PAWCD_set_ae_cd_sym ',id_sname)

    allocate(ylm_sym(nopr+af,25))
    ylm_sym=0.d0
    call set_ylms_symmtry_op(ia,25,vec0,ylm_sym)

    it=ityp(ia)

    n1nc=0.d0
    do lmt1=1,ilmt(it)
       il1=ltp(lmt1,it)
       im1=mtp(lmt1,it)
       it1=taup(lmt1,it)
       ii=(il1-1)**2+im1
       do lmt2=lmt1,ilmt(it)
          if(nspin==1 .and.dabs(hsr(ia,lmt1,lmt2,1)) < 1.d-10) cycle
          if(nspin==2) then
             if ( dabs(hsr(ia,lmt1,lmt2,1))+dabs(hsr(ia,lmt1,lmt2,2))<1.d-10 ) cycle
          endif
          il2=ltp(lmt2,it)
          im2=mtp(lmt2,it)
          it2=taup(lmt2,it)
          jj=(il2-1)**2+im2
          fac=2.d0;if(lmt1.eq.lmt2) fac=1.d0

          do is=1,nspin
             sum2=0.d0
             do iopr=1,nopr
                ja=abs(ia2ia_symmtry_op(ia,iopr))
                sum2=sum2+hsr(ja,lmt1,lmt2,is)*ylm_sym(iopr,ii)*ylm_sym(iopr,jj)
             end do
             sum2=sum2/dble(nopr)

!!$             do ir=1,nrc
             do ir=ista_nrc,min(nrc,iend_nrc)
                n1nc(ir,is)=n1nc(ir,is)+fac*psirpw(ir,il1,it1,it)*psirpw(ir,il2,it2,it)*sum2
             end do
          end do

       end do
    end do

    do is=1,nspin
!!$       do ir=1,nrc
       do ir=ista_nrc,min(nrc,iend_nrc)
          n1nc(ir,is)=(n1nc(ir,is)+rhcorpw(ir,it)/PAI4/dble(nspin)) &
               /radr_paw(ir,it)/radr_paw(ir,it)
       end do
    end do
    if ( sw_opencore == ON .and. nspin == 2 )then
       if ( has_opencore(it) == 1 ) then
          weight = 0.5d0 *mag_opencore_pol(ia,1)
          do ir=ista_nrc, min(nrc,iend_nrc)
             n1nc(ir,1) = n1nc(ir,1) +rmag_opencore(ir,it)/PAI4 *weight &
                  &                        /radr_paw(ir,it)/radr_paw(ir,it)
             n1nc(ir,2) = n1nc(ir,2) -rmag_opencore(ir,it)/PAI4 *weight &
                  &                        /radr_paw(ir,it)/radr_paw(ir,it)
          end do
       endif
    endif
    deallocate(ylm_sym)
    call tstatc0_end(id_sname)
    return
  end subroutine m_PAWCD_set_ae_cd_sym

  subroutine m_PAWCD_set_ae_cd2(ith,iph,ia,nspin,nrc,n1nc)
    integer,intent(in) :: ith,iph,ia,nspin,nrc
    real(8),intent(out):: n1nc(nrc,nspin)

    real(8):: n1(nrc,nspin),sum,fac, weight
    integer:: ir,it,is,lmt1,lmt2
    integer:: il1,il2,it1,it2,im1,im2
    integer:: ii,jj

    it=ityp(ia)
    do is=1,nspin
       do ir=1,nrc
          sum=0.d0
          do lmt1=1,ilmt(it)
             il1=ltp(lmt1,it)
             im1=mtp(lmt1,it)
             it1=taup(lmt1,it)
             ii=(il1-1)**2+im1
             do lmt2=lmt1,ilmt(it)
                il2=ltp(lmt2,it)
                im2=mtp(lmt2,it)
                it2=taup(lmt2,it)
                jj=(il2-1)**2+im2
                fac=2.d0;if(lmt1.eq.lmt2) fac=1.d0
                sum=sum+fac*hsr(ia,lmt1,lmt2,is)* &
                     psirpw(ir,il1,it1,it)* &
                     psirpw(ir,il2,it2,it)* &
                     ylm(ith,iph,ii)* &
                     ylm(ith,iph,jj)
             end do
          end do
          n1nc(ir,is)=(sum+rhcorpw(ir,it)/PAI4/dble(nspin)) &
               /radr_paw(ir,it)/radr_paw(ir,it)
       end do
    end do
    if ( sw_opencore == ON .and. nspin == 2 )then
       if ( has_opencore(it) == 1 ) then
          weight = 0.5d0 *mag_opencore_pol(ia,1)
          do ir=1,nrc
             n1nc(ir,1) = n1nc(ir,1) +rmag_opencore(ir,it)/PAI4 *weight &
                  &                        /radr_paw(ir,it)/radr_paw(ir,it)
             n1nc(ir,2) = n1nc(ir,2) -rmag_opencore(ir,it)/PAI4 *weight &
                  &                        /radr_paw(ir,it)/radr_paw(ir,it)
          end do
       endif
    endif

    return
  end subroutine m_PAWCD_set_ae_cd2

  subroutine m_PAWCD_set_ae_cd_sym2(vec0,ia,nspin,nrc,n1nc)
    integer,intent(in) :: ia,nspin,nrc
    real(DP),intent(in) :: vec0(3)
    real(8),intent(out):: n1nc(nrc,nspin)

    real(8):: n1(nrc,nspin),sum,fac,sum2, weight
    integer:: ir,it,is,lmt1,lmt2
    integer:: il1,il2,it1,it2,im1,im2
    integer:: ii,jj
    integer:: ja,iopr
    real(DP),pointer,dimension(:,:):: ylm_sym

    allocate(ylm_sym(nopr+af,25))
    ylm_sym=0.d0
    call set_ylms_symmtry_op(ia,25,vec0,ylm_sym)

    it=ityp(ia)
    do is=1,nspin
       do ir=1,nrc
          sum=0.d0
          do lmt1=1,ilmt(it)
             il1=ltp(lmt1,it)
             im1=mtp(lmt1,it)
             it1=taup(lmt1,it)
             ii=(il1-1)**2+im1
             do lmt2=lmt1,ilmt(it)
                il2=ltp(lmt2,it)
                im2=mtp(lmt2,it)
                it2=taup(lmt2,it)
                jj=(il2-1)**2+im2
                fac=2.d0;if(lmt1.eq.lmt2) fac=1.d0

                sum2=0.d0
                do iopr=1,nopr
                   ja=abs(ia2ia_symmtry_op(ia,iopr))
                   sum2=sum2+hsr(ja,lmt1,lmt2,is)* &
                        ylm_sym(iopr,ii)* &
                        ylm_sym(iopr,jj)
                end do
                sum2=sum2/dble(nopr)

                sum=sum+fac* &
                     psirpw(ir,il1,it1,it)* &
                     psirpw(ir,il2,it2,it)* &
                     sum2
             end do
          end do
          n1nc(ir,is)=(sum+rhcorpw(ir,it)/PAI4/dble(nspin)) &
               /radr_paw(ir,it)/radr_paw(ir,it)
       end do
    end do
    if ( sw_opencore == ON .and. nspin == 2 )then
       if ( has_opencore(it) == 1 ) then
          weight = 0.5d0 *mag_opencore_pol(ia,1)
          do ir=1,nrc
             n1nc(ir,1) = n1nc(ir,1) +rmag_opencore(ir,it)/PAI4 *weight &
                  &                        /radr_paw(ir,it)/radr_paw(ir,it)
             n1nc(ir,2) = n1nc(ir,2) -rmag_opencore(ir,it)/PAI4 *weight &
                  &                        /radr_paw(ir,it)/radr_paw(ir,it)
          end do
       endif
    endif

    deallocate(ylm_sym)
    return
  end subroutine m_PAWCD_set_ae_cd_sym2

  subroutine m_PAWCD_set_ae_vd_sym(vec0,ia,nspin,nrc,n1nc)
    integer,intent(in) :: ia,nspin,nrc
    real(DP),intent(in) :: vec0(3)
    real(8),intent(out):: n1nc(nrc,nspin)

    real(8):: n1(nrc,nspin),sum,fac,sum2
    integer:: ir,it,is,lmt1,lmt2
    integer:: il1,il2,it1,it2,im1,im2
    integer:: ii,jj
    integer:: ja,iopr
    real(DP),pointer,dimension(:,:):: ylm_sym

    allocate(ylm_sym(nopr+af,25))
    ylm_sym=0.d0
    call set_ylms_symmtry_op(ia,25,vec0,ylm_sym)

    it=ityp(ia)

    n1nc=0.d0
    do lmt1=1,ilmt(it)
       il1=ltp(lmt1,it)
       im1=mtp(lmt1,it)
       it1=taup(lmt1,it)
       ii=(il1-1)**2+im1
       do lmt2=lmt1,ilmt(it)
          il2=ltp(lmt2,it)
          im2=mtp(lmt2,it)
          it2=taup(lmt2,it)
          jj=(il2-1)**2+im2
          fac=2.d0;if(lmt1.eq.lmt2) fac=1.d0

          do is=1,nspin
             sum2=0.d0
             do iopr=1,nopr
                ja=abs(ia2ia_symmtry_op(ia,iopr))
                sum2=sum2+hsr(ja,lmt1,lmt2,is)* &
                     ylm_sym(iopr,ii)* &
                     ylm_sym(iopr,jj)
             end do
             sum2=sum2/dble(nopr)

             do ir=1,nrc
                n1nc(ir,is)=n1nc(ir,is)+fac* &
                     psirpw(ir,il1,it1,it)* &
                     psirpw(ir,il2,it2,it)* &
                     sum2
             end do
          end do

       end do
    end do

    do is=1,nspin
       do ir=1,nrc
          n1nc(ir,is)=n1nc(ir,is)/radr_paw(ir,it)/radr_paw(ir,it)
       end do
    end do
    deallocate(ylm_sym)
    return
  end subroutine m_PAWCD_set_ae_vd_sym

  subroutine m_PAWCD_set_ae_der_cd(ith,iph,ia,nspin,nrc,dn1_dth,dn1_dph,ista_nrc,iend_nrc)
    integer,intent(in) :: ith,iph,ia,nspin,nrc
    real(DP),intent(out):: dn1_dth(ista_nrc:iend_nrc,nspin)
    real(DP),intent(out):: dn1_dph(ista_nrc:iend_nrc,nspin)
    integer, intent(in) :: ista_nrc,iend_nrc

    real(DP):: sum1,sum2,fac,bb
    integer:: ir,it,is,lmt1,lmt2
    integer:: il1,il2,it1,it2,im1,im2
    integer:: ii,jj
    real(DP):: sin_theta
    integer :: id_sname = -1
    call tstatc0_begin('m_PAWCD_set_ae_der_cd ',id_sname)

    sin_theta=sqrt(1.d0-cos_theta(ith)**2)

    it=ityp(ia)
    dn1_dth=0.d0
    dn1_dph=0.d0
    do lmt1=1,ilmt(it)
       il1=ltp(lmt1,it)
       im1=mtp(lmt1,it)
       it1=taup(lmt1,it)
       ii=(il1-1)**2+im1
       do lmt2=lmt1,ilmt(it)
          il2=ltp(lmt2,it)
          im2=mtp(lmt2,it)
          it2=taup(lmt2,it)
          jj=(il2-1)**2+im2
          fac=2.d0;if(lmt1.eq.lmt2) fac=1.d0
          sum1=   dylm_dth(ith,iph,ii)* &
               ylm(ith,iph,jj)+ &
               ylm(ith,iph,ii)* &
               dylm_dth(ith,iph,jj)
          sum2=   dylm_dph(ith,iph,ii)* &
               ylm(ith,iph,jj)+ &
               ylm(ith,iph,ii)* &
               dylm_dph(ith,iph,jj)
          do is=1,nspin
             bb=fac*hsr(ia,lmt1,lmt2,is)
!!$             do ir=1,nrc
             do ir=ista_nrc,iend_nrc
                dn1_dth(ir,is)=dn1_dth(ir,is)+bb*sum1* &
                     psirpw(ir,il1,it1,it)* &
                     psirpw(ir,il2,it2,it)
                dn1_dph(ir,is)=dn1_dph(ir,is)+bb*sum2* &
                     psirpw(ir,il1,it1,it)* &
                     psirpw(ir,il2,it2,it)
             end do
          end do
       end do
    end do
    do is=1,nspin
!!$       do ir=1,nrc
       do ir=ista_nrc,iend_nrc
          dn1_dth(ir,is)=dn1_dth(ir,is)/radr_paw(ir,it)**3
          dn1_dph(ir,is)=dn1_dph(ir,is)/radr_paw(ir,it)**3
       end do
    end do
    dn1_dph=dn1_dph/sin_theta

    call tstatc0_end(id_sname)
    return
  end subroutine m_PAWCD_set_ae_der_cd

  subroutine m_PAWCD_set_ae_der_cd_sym(vec0,ith,iph,ia,nspin,nrc,dn1_dth,dn1_dph,ista_nrc,iend_nrc)
    real(DP),intent(in) :: vec0(3)
    integer,intent(in) :: ith,iph,ia,nspin,nrc
    real(DP),intent(out):: dn1_dth(ista_nrc:iend_nrc,nspin)
    real(DP),intent(out):: dn1_dph(ista_nrc:iend_nrc,nspin)
    integer, intent(in) :: ista_nrc,iend_nrc

!!$    real(DP):: n1(nrc,nspin),fac,bb,bb1,bb2
    real(DP):: fac,bb,bb1,bb2
    integer:: ir,it,is,lmt1,lmt2
    integer:: il1,il2,it1,it2,im1,im2
    integer:: ii,jj
    real(DP):: sin_theta
    integer:: ja,iopr,n
    real(DP),pointer,dimension(:,:):: ylm_sym
    real(DP),pointer,dimension(:,:):: dylm_dth_sym
    real(DP),pointer,dimension(:,:):: dylm_dph_sym
    integer :: id_sname = -1
    call tstatc0_begin('m_PAWCD_set_ae_der_cd_sym ',id_sname)

    call m_PP_find_maximum_l(n)
    n=(n-1)+(n-1)+1

    allocate(ylm_sym(nopr+af,n**2))
    allocate(dylm_dth_sym(nopr+af,n**2))
    allocate(dylm_dph_sym(nopr+af,n**2))
    ylm_sym=0.d0
    dylm_dth_sym=0.d0
    dylm_dph_sym=0.d0
    call set_diff_ylms_symmtry_op(ia,n**2,vec0,ylm_sym &
         & ,dylm_dth_sym,dylm_dph_sym)

    sin_theta=sqrt(1.d0-cos_theta(ith)**2)

    it=ityp(ia)
    dn1_dth=0.d0
    dn1_dph=0.d0
    do lmt1=1,ilmt(it)
       il1=ltp(lmt1,it)
       im1=mtp(lmt1,it)
       it1=taup(lmt1,it)
       ii=(il1-1)**2+im1
       do lmt2=lmt1,ilmt(it)
          if(nspin==1 .and.dabs(hsr(ia,lmt1,lmt2,1)) < 1.d-10) cycle
          if(nspin==2) then
             if ( dabs(hsr(ia,lmt1,lmt2,1))+dabs(hsr(ia,lmt1,lmt2,2))<1.d-10 ) cycle
          endif

          il2=ltp(lmt2,it)
          im2=mtp(lmt2,it)
          it2=taup(lmt2,it)
          jj=(il2-1)**2+im2
          fac=2.d0;if(lmt1.eq.lmt2) fac=1.d0

          do is=1,nspin
!                    bb1= fac*hsr(ia,lmt1,lmt2,is)* &
!                               (dylm_dth(ith,iph,ii)* &
!                                ylm(ith,iph,jj)+ &
!                                ylm(ith,iph,ii)* &
!                                dylm_dth(ith,iph,jj))
!                    bb2= fac*hsr(ia,lmt1,lmt2,is)*&
!                               (dylm_dph(ith,iph,ii)* &
!                                ylm(ith,iph,jj)+ &
!                                ylm(ith,iph,ii)* &
!                                dylm_dph(ith,iph,jj))

             bb1=0.d0
             do iopr=1,nopr
                ja=abs(ia2ia_symmtry_op(ia,iopr))
                bb1=bb1+hsr(ja,lmt1,lmt2,is)* &
                     (dylm_dth_sym(iopr,ii)* &
                     ylm_sym(iopr,jj)+ &
                     ylm_sym(iopr,ii)* &
                     dylm_dth_sym(iopr,jj))
             end do
             bb1=bb1/dble(nopr)*fac

             bb2=0.d0
             do iopr=1,nopr
                ja=abs(ia2ia_symmtry_op(ia,iopr))
                bb2=bb2+hsr(ja,lmt1,lmt2,is)* &
                     (dylm_dph_sym(iopr,ii)* &
                     ylm_sym(iopr,jj)+ &
                     ylm_sym(iopr,ii)* &
                     dylm_dph_sym(iopr,jj))
             end do
             bb2=bb2/dble(nopr)*fac

!!$             do ir=1,nrc
             do ir=ista_nrc,iend_nrc
                bb= psirpw(ir,il1,it1,it) * psirpw(ir,il2,it2,it)
                dn1_dth(ir,is)=dn1_dth(ir,is)+bb*bb1
                dn1_dph(ir,is)=dn1_dph(ir,is)+bb*bb2
             end do
          end do

       end do
    end do

    do is=1,nspin
!!$       do ir=1,nrc
       do ir=ista_nrc,iend_nrc
          dn1_dth(ir,is)=dn1_dth(ir,is)/radr_paw(ir,it)**3
          dn1_dph(ir,is)=dn1_dph(ir,is)/radr_paw(ir,it)**3
       end do
    end do
    dn1_dph=dn1_dph/sin_theta

    deallocate(ylm_sym,dylm_dth_sym,dylm_dph_sym)
    call tstatc0_end(id_sname)
    return
  end subroutine m_PAWCD_set_ae_der_cd_sym

  subroutine m_PAWCD_set_ae_der_cd2(ith,iph,ia,nspin,nrc,dn1_dth,dn1_dph)
    integer,intent(in) :: ith,iph,ia,nspin,nrc
    real(8),intent(out):: dn1_dth(nrc,nspin)
    real(8),intent(out):: dn1_dph(nrc,nspin)

    real(8):: n1(nrc,nspin),sum1,sum2,fac,bb
    integer:: ir,it,is,lmt1,lmt2
    integer:: il1,il2,it1,it2,im1,im2
    integer:: ii,jj
    real(8):: sin_theta

    sin_theta=sqrt(1.d0-cos_theta(ith)**2)

    it=ityp(ia)
    do is=1,nspin
       do ir=1,nrc
          sum1=0.d0
          sum2=0.d0
          do lmt1=1,ilmt(it)
             il1=ltp(lmt1,it)
             im1=mtp(lmt1,it)
             it1=taup(lmt1,it)
             ii=(il1-1)**2+im1
             do lmt2=lmt1,ilmt(it)
                il2=ltp(lmt2,it)
                im2=mtp(lmt2,it)
                it2=taup(lmt2,it)
                jj=(il2-1)**2+im2
                fac=2.d0;if(lmt1.eq.lmt2) fac=1.d0

                bb=fac*hsr(ia,lmt1,lmt2,is)* &
                     psirpw(ir,il1,it1,it)* &
                     psirpw(ir,il2,it2,it)
                sum1=sum1+bb* &
                     (dylm_dth(ith,iph,ii)* &
                     ylm(ith,iph,jj)+ &
                     ylm(ith,iph,ii)* &
                     dylm_dth(ith,iph,jj))
                sum2=sum2+bb* &
                     (dylm_dph(ith,iph,ii)* &
                     ylm(ith,iph,jj)+ &
                     ylm(ith,iph,ii)* &
                     dylm_dph(ith,iph,jj))
             end do
          end do
          dn1_dth(ir,is)=sum1/radr_paw(ir,it)**3
          dn1_dph(ir,is)=sum2/radr_paw(ir,it)**3
       end do
    end do
    dn1_dph=dn1_dph/sin_theta

    return
  end subroutine m_PAWCD_set_ae_der_cd2

  subroutine m_PAWCD_set_ps_cd(ith,iph,ia,nspin,nrc,nps,ista_nrc,iend_nrc)
    integer,intent(in) :: ith,iph,ia,nspin,nrc
    real(DP),intent(out):: nps(ista_nrc:iend_nrc,nspin)
    integer,intent(in) :: ista_nrc,iend_nrc

    real(DP):: sum,fac, weight
    integer:: ir,it,is,lmt1,lmt2
    integer:: il1,il2,it1,it2,im1,im2,l3
    integer:: ii,jj,n,ilm3,iiqitg
    integer,pointer,dimension(:):: il3
    real(DP),pointer,dimension(:):: qsum

    call m_PP_find_maximum_l(n)
    n=(n-1)+(n-1)+1
    allocate(il3(n**2));call substitute_il3(n**2,il3)
    allocate(qsum(ista_nrc:iend_nrc))

    it=ityp(ia)
    nps=0.d0
    do lmt1=1,ilmt(it)
       il1=ltp(lmt1,it)
       im1=mtp(lmt1,it)
       it1=taup(lmt1,it)
       ii=(il1-1)**2+im1
       do lmt2=lmt1,ilmt(it)
          il2=ltp(lmt2,it)
          im2=mtp(lmt2,it)
          it2=taup(lmt2,it)
          jj=(il2-1)**2+im2
          fac=2.d0;if(lmt1.eq.lmt2) fac=1.d0
          qsum=0.d0
          do n=1,il2p(lmt1,lmt2,it)
             ilm3=isph(lmt1,lmt2,n,it);l3=il3(ilm3)
             iiqitg=iqitg(il1,it1,il2,it2,l3+1,it)
             if(iiqitg==0) cycle
!!$             do ir=1,nrc
             do ir=ista_nrc,iend_nrc
                qsum(ir)=qsum(ir)+qrspspw(ir,iiqitg)* &
                     & dl2p(lmt1,lmt2,n,it)*ylm(ith,iph,ilm3)
             end do
          end do
          do is=1,nspin
!!$             do ir=1,nrc
             do ir=ista_nrc,iend_nrc
                nps(ir,is)=nps(ir,is)+fac*hsr(ia,lmt1,lmt2,is) &
                     & *(phirpw(ir,il1,it1,it)*phirpw(ir,il2,it2,it) &
                     &    *ylm(ith,iph,ii)*ylm(ith,iph,jj)+qsum(ir))
             end do
          end do
       end do
    end do
    do is=1,nspin
!!$       do ir=1,nrc
       do ir=ista_nrc,iend_nrc
          nps(ir,is)=(nps(ir,is)+rhpcrpw(ir,it)/PAI4/dble(nspin)) &
               /radr_paw(ir,it)/radr_paw(ir,it)
          if(nps(ir,is).lt.0.d0) nps(ir,is)=0.d0
       end do
    end do
    if ( sw_opencore == ON .and. nspin == 2 .and. sw_xc_opencore_ae_only ==OFF )then
       if ( has_opencore(it) == 1 ) then
          weight = 0.5d0 *mag_opencore_pol(ia,1)
          do ir=ista_nrc, iend_nrc
             nps(ir,1) = nps(ir,1) +rmag_opencore(ir,it)/PAI4 *weight &
                  &                        /radr_paw(ir,it)/radr_paw(ir,it)
             nps(ir,2) = nps(ir,2) -rmag_opencore(ir,it)/PAI4 *weight &
                  &                        /radr_paw(ir,it)/radr_paw(ir,it)
          end do
       endif
    endif

    deallocate(il3,qsum)
    return
  end subroutine m_PAWCD_set_ps_cd

  subroutine m_PAWCD_set_ps_cd_sym(vec0,ia,nspin,nrc,nps,ista_nrc,iend_nrc)
    integer,intent(in) :: ia,nspin,nrc
    real(DP),intent(in):: vec0(3)
    real(DP),intent(out):: nps(ista_nrc:iend_nrc,nspin)
    integer, intent(in) :: ista_nrc,iend_nrc

    real(DP):: sum,fac,sum2,sum3, weight
    integer:: ir,it,is,lmt1,lmt2
    integer:: il1,il2,it1,it2,im1,im2,l3
    integer:: ii,jj,n,ilm3,iiqitg
    integer,allocatable,dimension(:):: il3
    real(DP),allocatable,dimension(:):: qsum
    integer:: ja,iopr
    real(DP),allocatable,dimension(:,:):: ylm_sym
    integer :: id_sname = -1
    call tstatc0_begin('m_PAWCD_set_ps_cd_sym ',id_sname)

    call m_PP_find_maximum_l(n)
    n=(n-1)+(n-1)+1
    allocate(il3(n**2));call substitute_il3(n**2,il3)

    allocate(ylm_sym(nopr+af,n**2));   ylm_sym=0.d0
    call set_ylms_symmtry_op(ia,n**2,vec0,ylm_sym)

    allocate(qsum(ista_nrc:iend_nrc))

    it=ityp(ia)
    sum=0.d0
    nps=0.d0
    do lmt1=1,ilmt(it)
       il1=ltp(lmt1,it)
       im1=mtp(lmt1,it)
       it1=taup(lmt1,it)
       ii=(il1-1)**2+im1
       do lmt2=lmt1,ilmt(it)
          if(nspin==1 .and.dabs(hsr(ia,lmt1,lmt2,1)) < 1.d-10) cycle
          if(nspin==2) then
             if ( dabs(hsr(ia,lmt1,lmt2,1))+dabs(hsr(ia,lmt1,lmt2,2))<1.d-10 ) cycle
          endif

          il2=ltp(lmt2,it)
          im2=mtp(lmt2,it)
          it2=taup(lmt2,it)
          jj=(il2-1)**2+im2
          fac=2.d0;if(lmt1.eq.lmt2) fac=1.d0

          do is=1,nspin

             sum2=0.d0
             do iopr=1,nopr
                ja=abs(ia2ia_symmtry_op(ia,iopr))
                sum2=sum2+hsr(ja,lmt1,lmt2,is)*ylm_sym(iopr,ii)*ylm_sym(iopr,jj)
             end do
             sum2=sum2/dble(nopr)

             qsum=0.d0
             do n=1,il2p(lmt1,lmt2,it)
                ilm3=isph(lmt1,lmt2,n,it);l3=il3(ilm3)
                iiqitg=iqitg(il1,it1,il2,it2,l3+1,it)
                if(iiqitg==0) cycle

                sum3=0.d0
                do iopr=1,nopr
                   ja=abs(ia2ia_symmtry_op(ia,iopr))
                   sum3=sum3+hsr(ja,lmt1,lmt2,is)*ylm_sym(iopr,ilm3)
                end do
                sum3=sum3/dble(nopr)

!!$                do ir=1,nrc
                do ir=ista_nrc,min(nrc,iend_nrc)
                   qsum(ir)=qsum(ir)+qrspspw(ir,iiqitg)*dl2p(lmt1,lmt2,n,it)*sum3
                end do
             end do

!!$             do ir=1,nrc
             do ir=ista_nrc,min(nrc,iend_nrc)
                nps(ir,is)=nps(ir,is) &
                     & +fac*phirpw(ir,il1,it1,it)*phirpw(ir,il2,it2,it)*sum2+fac*qsum(ir)

             end do
          end do

       end do
    end do

    do is=1,nspin
!!$       do ir=1,nrc
       do ir=ista_nrc,min(nrc,iend_nrc)
          nps(ir,is)=(nps(ir,is)+rhpcrpw(ir,it)/PAI4/dble(nspin)) &
               /radr_paw(ir,it)/radr_paw(ir,it)
          if(nps(ir,is).lt.0.d0) nps(ir,is)=0.d0
       end do
    end do
    if ( sw_opencore == ON .and. nspin == 2 .and. sw_xc_opencore_ae_only ==OFF )then
       if ( has_opencore(it) == 1 ) then
          weight = 0.5d0 *mag_opencore_pol(ia,1)
          do ir=ista_nrc,min(nrc,iend_nrc)
             nps(ir,1) = nps(ir,1) +rmag_opencore(ir,it)/PAI4 *weight &
                  &                        /radr_paw(ir,it)/radr_paw(ir,it)
             nps(ir,2) = nps(ir,2) -rmag_opencore(ir,it)/PAI4 *weight &
                  &                        /radr_paw(ir,it)/radr_paw(ir,it)
          end do
       endif
    endif

    deallocate(il3,ylm_sym,qsum)
    call tstatc0_end(id_sname)
    return
  end subroutine m_PAWCD_set_ps_cd_sym

  subroutine m_PAWCD_set_ps_cd2(ith,iph,ia,nspin,nrc,nps)
    integer,intent(in) :: ith,iph,ia,nspin,nrc
    real(8),intent(out):: nps(nrc,nspin)

    real(8):: sum,fac,qsum, weight
    integer:: ir,it,is,lmt1,lmt2
    integer:: il1,il2,it1,it2,im1,im2,l3
    integer:: ii,jj,n,ilm3,iiqitg
    integer,pointer,dimension(:):: il3

    call m_PP_find_maximum_l(n)
    n=(n-1)+(n-1)+1
    allocate(il3(n**2));call substitute_il3(n**2,il3)

    it=ityp(ia)
    do is=1,nspin
       do ir=1,nrc
          sum=0.d0
          do lmt1=1,ilmt(it)
             il1=ltp(lmt1,it)
             im1=mtp(lmt1,it)
             it1=taup(lmt1,it)
             ii=(il1-1)**2+im1
             do lmt2=lmt1,ilmt(it)
                il2=ltp(lmt2,it)
                im2=mtp(lmt2,it)
                it2=taup(lmt2,it)
                jj=(il2-1)**2+im2
                fac=2.d0;if(lmt1.eq.lmt2) fac=1.d0
                qsum=0.d0
                do n=1,il2p(lmt1,lmt2,it)
                   ilm3=isph(lmt1,lmt2,n,it);l3=il3(ilm3)
                   iiqitg=iqitg(il1,it1,il2,it2,l3+1,it)
                   if(iiqitg==0) cycle
                   qsum=qsum+qrspspw(ir,iiqitg)* &
                        dl2p(lmt1,lmt2,n,it)* &
                        ylm(ith,iph,ilm3)
                end do
                sum=sum+fac*hsr(ia,lmt1,lmt2,is)* &
                     (phirpw(ir,il1,it1,it)* &
                     phirpw(ir,il2,it2,it)* &
                     ylm(ith,iph,ii)* &
                     ylm(ith,iph,jj)+qsum)
             end do
          end do
          nps(ir,is)=(sum+rhpcrpw(ir,it)/PAI4/dble(nspin)) &
               /radr_paw(ir,it)/radr_paw(ir,it)
          if(nps(ir,is).lt.0.d0) nps(ir,is)=0.d0
       end do
    end do
    if ( sw_opencore == ON .and. nspin == 2 .and. sw_xc_opencore_ae_only ==OFF )then
       if ( has_opencore(it) == 1 ) then
          weight = 0.5d0 *mag_opencore_pol(ia,1)
          do ir=1,nrc
             nps(ir,1) = nps(ir,1) +rmag_opencore(ir,it)/PAI4 *weight &
                  &                        /radr_paw(ir,it)/radr_paw(ir,it)
             nps(ir,2) = nps(ir,2) -rmag_opencore(ir,it)/PAI4 *weight &
                  &                        /radr_paw(ir,it)/radr_paw(ir,it)
          end do
       endif
    endif

    deallocate(il3)
    return
  end subroutine m_PAWCD_set_ps_cd2

  subroutine m_PAWCD_set_ps_cd_sym2(vec0,ia,nspin,nrc,nps)
    integer,intent(in) :: ia,nspin,nrc
    real(DP),intent(in):: vec0(3)
    real(8),intent(out):: nps(nrc,nspin)

    real(8):: sum,fac,qsum,sum2,sum3, weight
    integer:: ir,it,is,lmt1,lmt2
    integer:: il1,il2,it1,it2,im1,im2,l3
    integer:: ii,jj,n,ilm3,iiqitg
    integer,pointer,dimension(:):: il3
    integer:: ja,iopr
    real(DP),pointer,dimension(:,:):: ylm_sym

    call m_PP_find_maximum_l(n)
    n=(n-1)+(n-1)+1
    allocate(il3(n**2));call substitute_il3(n**2,il3)

    allocate(ylm_sym(nopr+af,25))
    ylm_sym=0.d0
    call set_ylms_symmtry_op(ia,25,vec0,ylm_sym)

    it=ityp(ia)
    do is=1,nspin
       do ir=1,nrc
          sum=0.d0
          do lmt1=1,ilmt(it)
             il1=ltp(lmt1,it)
             im1=mtp(lmt1,it)
             it1=taup(lmt1,it)
             ii=(il1-1)**2+im1
             do lmt2=lmt1,ilmt(it)
                il2=ltp(lmt2,it)
                im2=mtp(lmt2,it)
                it2=taup(lmt2,it)
                jj=(il2-1)**2+im2
                fac=2.d0;if(lmt1.eq.lmt2) fac=1.d0

                sum2=0.d0
                do iopr=1,nopr
                   ja=abs(ia2ia_symmtry_op(ia,iopr))
                   sum2=sum2+hsr(ja,lmt1,lmt2,is)* &
                        ylm_sym(iopr,ii)* &
                        ylm_sym(iopr,jj)
                end do
                sum2=sum2/dble(nopr)

                qsum=0.d0
                do n=1,il2p(lmt1,lmt2,it)
                   ilm3=isph(lmt1,lmt2,n,it);l3=il3(ilm3)
                   iiqitg=iqitg(il1,it1,il2,it2,l3+1,it)
                   if(iiqitg==0) cycle

                   sum3=0.d0
                   do iopr=1,nopr
                      ja=abs(ia2ia_symmtry_op(ia,iopr))
                      sum3=sum3+hsr(ja,lmt1,lmt2,is)* &
                           ylm_sym(iopr,ilm3)
                   end do
                   sum3=sum3/dble(nopr)

                   qsum=qsum+qrspspw(ir,iiqitg)* &
                        dl2p(lmt1,lmt2,n,it)* &
                        sum3
                end do

                sum=sum+fac* &
                     phirpw(ir,il1,it1,it)* &
                     phirpw(ir,il2,it2,it)* &
                     sum2 + &
                     fac*qsum

             end do
          end do
          nps(ir,is)=(sum+rhpcrpw(ir,it)/PAI4/dble(nspin)) &
               /radr_paw(ir,it)/radr_paw(ir,it)
          if(nps(ir,is).lt.0.d0) nps(ir,is)=0.d0
       end do
    end do
    if ( sw_opencore == ON .and. nspin == 2 .and. sw_xc_opencore_ae_only ==OFF )then
       if ( has_opencore(it) == 1 ) then
          weight = 0.5d0 *mag_opencore_pol(ia,1)
          do ir=1,nrc
             nps(ir,1) = nps(ir,1) +rmag_opencore(ir,it)/PAI4 *weight &
                  &                        /radr_paw(ir,it)/radr_paw(ir,it)
             nps(ir,2) = nps(ir,2) -rmag_opencore(ir,it)/PAI4 *weight &
                  &                        /radr_paw(ir,it)/radr_paw(ir,it)
          end do
       endif
    endif

    deallocate(il3,ylm_sym)
    return
  end subroutine m_PAWCD_set_ps_cd_sym2

  subroutine m_PAWCD_set_ps_der_cd(ith,iph,ia,nspin,nrc,dnps_dth,dnps_dph,ista_nrc,iend_nrc)
    integer,intent(in) :: ith,iph,ia,nspin,nrc
    real(DP),intent(out):: dnps_dth(nrc,nspin)
    real(DP),intent(out):: dnps_dph(nrc,nspin)
    integer, intent(in) :: ista_nrc,iend_nrc

    real(DP):: sum1,sum2,fac,qc,sin_theta
    integer:: ir,it,is,lmt1,lmt2
    integer:: il1,il2,it1,it2,im1,im2,l3
    integer:: ii,jj,n,ilm3,iiqitg
    integer,pointer,dimension(:):: il3
    real(DP),pointer,dimension(:):: qsum1,qsum2

    allocate(qsum1(nrc))
    allocate(qsum2(nrc))

    call m_PP_find_maximum_l(n)
    n=(n-1)+(n-1)+1
    allocate(il3(n**2));call substitute_il3(n**2,il3)

    sin_theta=sqrt(1.d0-cos_theta(ith)**2)

    it=ityp(ia)
    dnps_dth=0.d0
    dnps_dph=0.d0
    do lmt1=1,ilmt(it)
       il1=ltp(lmt1,it)
       im1=mtp(lmt1,it)
       it1=taup(lmt1,it)
       ii=(il1-1)**2+im1
       do lmt2=lmt1,ilmt(it)
          if(nspin==1 .and.dabs(hsr(ia,lmt1,lmt2,1)) < 1.d-10) cycle
          if(nspin==2) then
             if ( dabs(hsr(ia,lmt1,lmt2,1))+dabs(hsr(ia,lmt1,lmt2,2))<1.d-10 ) cycle
          endif

          il2=ltp(lmt2,it)
          im2=mtp(lmt2,it)
          it2=taup(lmt2,it)
          jj=(il2-1)**2+im2
          fac=2.d0;if(lmt1.eq.lmt2) fac=1.d0
          qsum1=0.d0
          qsum2=0.d0
          do n=1,il2p(lmt1,lmt2,it)
             ilm3=isph(lmt1,lmt2,n,it);l3=il3(ilm3)
             iiqitg=iqitg(il1,it1,il2,it2,l3+1,it)
             if(iiqitg==0) cycle
!!$             do ir=1,nrc
             do ir=ista_nrc,iend_nrc
                qc=qrspspw(ir,iiqitg)*dl2p(lmt1,lmt2,n,it)
                qsum1(ir)=qsum1(ir)+qc*dylm_dth(ith,iph,ilm3)
                qsum2(ir)=qsum2(ir)+qc*dylm_dph(ith,iph,ilm3)
             end do
          end do
          sum1=dylm_dth(ith,iph,ii)*ylm(ith,iph,jj)+ ylm(ith,iph,ii)*dylm_dth(ith,iph,jj)
          sum2=dylm_dph(ith,iph,ii)*ylm(ith,iph,jj)+ ylm(ith,iph,ii)*dylm_dph(ith,iph,jj)
          do is=1,nspin
!!$             do ir=1,nrc
             do ir=ista_nrc,iend_nrc
                dnps_dth(ir,is)=dnps_dth(ir,is)+fac*hsr(ia,lmt1,lmt2,is) &
                     & *(phirpw(ir,il1,it1,it)*phirpw(ir,il2,it2,it)*sum1+qsum1(ir))
                dnps_dph(ir,is)=dnps_dph(ir,is)+fac*hsr(ia,lmt1,lmt2,is) &
                     & *(phirpw(ir,il1,it1,it)*phirpw(ir,il2,it2,it)*sum2+qsum2(ir))
             end do
          end do
       end do
    end do
    do is=1,nspin
!!$       do ir=1,nrc
       do ir=ista_nrc,iend_nrc
          dnps_dth(ir,is)=dnps_dth(ir,is)/radr_paw(ir,it)**3
          dnps_dph(ir,is)=dnps_dph(ir,is)/radr_paw(ir,it)**3
       end do
    end do

    dnps_dph=dnps_dph/sin_theta

    deallocate(il3,qsum1,qsum2)
    return
  end subroutine m_PAWCD_set_ps_der_cd

  subroutine m_PAWCD_set_ps_der_cd_sym(vec0,ith,iph,ia,nspin,nrc,dnps_dth,dnps_dph,ista_nrc,iend_nrc)
    real(DP),intent(in) :: vec0(3)
    integer,intent(in) :: ith,iph,ia,nspin,nrc
    real(DP),intent(out):: dnps_dth(ista_nrc:iend_nrc,nspin)
    real(DP),intent(out):: dnps_dph(ista_nrc:iend_nrc,nspin)
    integer,intent(in) :: ista_nrc,iend_nrc

    real(DP):: fac,qc,bb1,bb2,sin_theta,cc1,cc2
    integer:: ir,it,is,lmt1,lmt2
    integer:: il1,il2,it1,it2,im1,im2,l3
    integer:: ii,jj,n,ilm3,iiqitg
    integer,pointer,dimension(:):: il3
    real(DP),pointer,dimension(:):: qsum1,qsum2
    integer:: ja,iopr
    real(DP),pointer,dimension(:,:):: ylm_sym
    real(DP),pointer,dimension(:,:):: dylm_dth_sym
    real(DP),pointer,dimension(:,:):: dylm_dph_sym
    integer :: id_sname = -1
    call tstatc0_begin('m_PAWCD_set_ps_der_cd_sym ',id_sname)

    call m_PP_find_maximum_l(n)
    n=(n-1)+(n-1)+1
    allocate(il3(n**2));call substitute_il3(n**2,il3)

    allocate(ylm_sym(nopr+af,n**2))
    allocate(dylm_dth_sym(nopr+af,n**2))
    allocate(dylm_dph_sym(nopr+af,n**2))
    ylm_sym=0.d0
    dylm_dth_sym=0.d0
    dylm_dph_sym=0.d0
    call set_diff_ylms_symmtry_op(ia,n**2,vec0,ylm_sym,dylm_dth_sym,dylm_dph_sym)

    allocate(qsum1(ista_nrc:iend_nrc))
    allocate(qsum2(ista_nrc:iend_nrc))

    sin_theta=sqrt(1.d0-cos_theta(ith)**2)

    it=ityp(ia)
    dnps_dth=0.d0
    dnps_dph=0.d0
    do lmt1=1,ilmt(it)
       il1=ltp(lmt1,it)
       im1=mtp(lmt1,it)
       it1=taup(lmt1,it)
       ii=(il1-1)**2+im1
       do lmt2=lmt1,ilmt(it)
          if(nspin==1 .and.dabs(hsr(ia,lmt1,lmt2,1)) < 1.d-10) cycle
          if(nspin==2) then
             if ( dabs(hsr(ia,lmt1,lmt2,1))+dabs(hsr(ia,lmt1,lmt2,2))<1.d-10 ) cycle
          endif

          il2=ltp(lmt2,it)
          im2=mtp(lmt2,it)
          it2=taup(lmt2,it)
          jj=(il2-1)**2+im2
          fac=2.d0;if(lmt1.eq.lmt2) fac=1.d0

          if(ipripaw >= 2) then
             iiqitg = 0
             do n=1,il2p(lmt1,lmt2,it)
                ilm3=isph(lmt1,lmt2,n,it);l3=il3(ilm3)
                if(iqitg(il1,it1,il2,it2,l3+1,it)>0) iiqitg=iiqitg+1
             end do
             write(nfout,'(" it, lmt1, lmt2, il2p iiqitg/=0 = ",5i8)') it, lmt1, lmt2,il2p(lmt1,lmt2,it),iiqitg
             if(iiqitg > 0) then
                do n=1,il2p(lmt1,lmt2,it)
                   ilm3=isph(lmt1,lmt2,n,it);l3=il3(ilm3)
!                         write(nfout,'("        il1,im1,it1, il2,im2,it2, n, iiqitg = ",8i6)') il1,im1,it1,il2,im2,it2, n, iqitg(il1,it1,il2,it2,l3+1,it)
                end do
             end if
          end if

          do is=1,nspin
!                    bb1=fac*hsr(ia,lmt1,lmt2,is)*( &
!                                dylm_dth(ith,iph,ii)* &
!                                ylm(ith,iph,jj)+ &
!                                ylm(ith,iph,ii)* &
!                                dylm_dth(ith,iph,jj))

             bb1=0.d0
             do iopr=1,nopr
                ja=abs(ia2ia_symmtry_op(ia,iopr))
                bb1=bb1+hsr(ja,lmt1,lmt2,is)*( &
                     dylm_dth_sym(iopr,ii)* &
                     ylm_sym(iopr,jj)+ &
                     ylm_sym(iopr,ii)* &
                     dylm_dth_sym(iopr,jj))
             end do
             bb1=bb1/dble(nopr)*fac

!                    bb2=fac*hsr(ia,lmt1,lmt2,is)*( &
!                                dylm_dph(ith,iph,ii)* &
!                                ylm(ith,iph,jj)+ &
!                                ylm(ith,iph,ii)* &
!                                dylm_dph(ith,iph,jj))

             bb2=0.d0
             do iopr=1,nopr
                ja=abs(ia2ia_symmtry_op(ia,iopr))
                bb2=bb2+hsr(ja,lmt1,lmt2,is)*( &
                     dylm_dph_sym(iopr,ii)* &
                     ylm_sym(iopr,jj)+ &
                     ylm_sym(iopr,ii)* &
                     dylm_dph_sym(iopr,jj))
             end do
             bb2=bb2/dble(nopr)*fac

             qsum1=0.d0
             qsum2=0.d0
             do n=1,il2p(lmt1,lmt2,it)
                ilm3=isph(lmt1,lmt2,n,it);l3=il3(ilm3)
                iiqitg=iqitg(il1,it1,il2,it2,l3+1,it)
                if(iiqitg==0) cycle

                cc1=0.d0
                do iopr=1,nopr
                   ja=abs(ia2ia_symmtry_op(ia,iopr))
                   cc1=cc1+dylm_dth_sym(iopr,ilm3)*hsr(ja,lmt1,lmt2,is)
                end do
                cc1=cc1/dble(nopr)*dl2p(lmt1,lmt2,n,it)

                cc2=0.d0
                do iopr=1,nopr
                   ja=abs(ia2ia_symmtry_op(ia,iopr))
                   cc2=cc2+dylm_dph_sym(iopr,ilm3)*hsr(ja,lmt1,lmt2,is)
                end do
                cc2=cc2/dble(nopr)*dl2p(lmt1,lmt2,n,it)

!!$                do ir=1,nrc
                do ir=ista_nrc,iend_nrc
                   qc=qrspspw(ir,iiqitg)
                   qsum1(ir)=qsum1(ir)+qc*cc1
                   qsum2(ir)=qsum2(ir)+qc*cc2
                end do
             end do
             qsum1=qsum1*fac
             qsum2=qsum2*fac

!!$             do ir=1,nrc
             do ir=ista_nrc,iend_nrc

                dnps_dth(ir,is)=dnps_dth(ir,is) &
                     & + bb1*phirpw(ir,il1,it1,it)*phirpw(ir,il2,it2,it)+qsum1(ir)

                dnps_dph(ir,is)=dnps_dph(ir,is) &
                     & + bb2*phirpw(ir,il1,it1,it)*phirpw(ir,il2,it2,it)+qsum2(ir)

             end do
          end do

       end do
    end do

    do is=1,nspin
!!$       do ir=1,nrc
       do ir=ista_nrc,iend_nrc
          dnps_dth(ir,is)=dnps_dth(ir,is)/radr_paw(ir,it)**3
          dnps_dph(ir,is)=dnps_dph(ir,is)/radr_paw(ir,it)**3
       end do
    end do

    dnps_dph=dnps_dph/sin_theta

    deallocate(il3,qsum1,qsum2)
    deallocate(ylm_sym,dylm_dth_sym,dylm_dph_sym)
    call tstatc0_end(id_sname)
    return
  end subroutine m_PAWCD_set_ps_der_cd_sym

  subroutine m_PAWCD_set_ps_der_cd2(ith,iph,ia,nspin,nrc,dnps_dth,dnps_dph)
    integer,intent(in) :: ith,iph,ia,nspin,nrc
    real(8),intent(out):: dnps_dth(nrc,nspin)
    real(8),intent(out):: dnps_dph(nrc,nspin)

    real(8):: sum1,sum2,fac,qsum1,qsum2,qc,sin_theta
    integer:: ir,it,is,lmt1,lmt2
    integer:: il1,il2,it1,it2,im1,im2,l3
    integer:: ii,jj,n,ilm3,iiqitg
    integer,pointer,dimension(:):: il3

    call m_PP_find_maximum_l(n)
    n=(n-1)+(n-1)+1
    allocate(il3(n**2));call substitute_il3(n**2,il3)

    sin_theta=sqrt(1.d0-cos_theta(ith)**2)

    it=ityp(ia)
    do is=1,nspin
       do ir=1,nrc
          sum1=0.d0
          sum2=0.d0
          do lmt1=1,ilmt(it)
             il1=ltp(lmt1,it)
             im1=mtp(lmt1,it)
             it1=taup(lmt1,it)
             ii=(il1-1)**2+im1
             do lmt2=lmt1,ilmt(it)
                il2=ltp(lmt2,it)
                im2=mtp(lmt2,it)
                it2=taup(lmt2,it)
                jj=(il2-1)**2+im2
                fac=2.d0;if(lmt1.eq.lmt2) fac=1.d0
                qsum1=0.d0
                qsum2=0.d0
                do n=1,il2p(lmt1,lmt2,it)
                   ilm3=isph(lmt1,lmt2,n,it);l3=il3(ilm3)
                   iiqitg=iqitg(il1,it1,il2,it2,l3+1,it)
                   if(iiqitg==0) cycle
                   qc=qrspspw(ir,iiqitg)* &
                        dl2p(lmt1,lmt2,n,it)
                   qsum1=qsum1+qc*dylm_dth(ith,iph,ilm3)
                   qsum2=qsum2+qc*dylm_dph(ith,iph,ilm3)
                end do
                sum1=sum1+fac*hsr(ia,lmt1,lmt2,is)* &
                     (phirpw(ir,il1,it1,it)* &
                     phirpw(ir,il2,it2,it)* &
                     (dylm_dth(ith,iph,ii)* &
                     ylm(ith,iph,jj)+ &
                     ylm(ith,iph,ii)* &
                     dylm_dth(ith,iph,jj))+qsum1)
                sum2=sum2+fac*hsr(ia,lmt1,lmt2,is)* &
                     (phirpw(ir,il1,it1,it)* &
                     phirpw(ir,il2,it2,it)* &
                     (dylm_dph(ith,iph,ii)* &
                     ylm(ith,iph,jj)+ &
                     ylm(ith,iph,ii)* &
                            dylm_dph(ith,iph,jj))+qsum2)
             end do
          end do
          dnps_dth(ir,is)=sum1/radr_paw(ir,it)**3
          dnps_dph(ir,is)=sum2/radr_paw(ir,it)**3
       end do
    end do

    dnps_dph=dnps_dph/sin_theta

    deallocate(il3)
    return
  end subroutine m_PAWCD_set_ps_der_cd2

  subroutine m_PAWCD_wd_cd(nfout)
    integer,intent(in):: nfout
    integer:: unit
    integer:: iloop !,ie,je,ke
    integer:: idp,nlp,nmp,nnp
    real(DP),allocatable,dimension(:) :: afft
    real(DP),pointer,dimension(:)      :: afft_mpi1
    real(DP),pointer,dimension(:)      :: afft_mpi2
    real(DP),pointer,dimension(:)      :: afft_mpi3
    real(DP),pointer,dimension(:,:,:)  :: chg_rs
    logical:: flg_wrtn(natm)

    call m_FFT_alloc_CD_box()
    allocate(afft(ista_fftp:iend_fftp))
    allocate(afft_mpi1(nfftp))
    if(npes >= 2) then
       allocate(afft_mpi2(mp_fftp))
       allocate(afft_mpi3(mp_fftp))
    end if
    if(mype==0) then
       flg_wrtn=.false.
!            ie=fft_box_size_CD(1)
!            je=fft_box_size_CD(2)
!            ke=fft_box_size_CD(3)
!!$            idp = fft_box_size_CD(0)
!!$            nlp = fft_box_size_CD(1)
!!$            nmp = fft_box_size_CD(2)
!!$            nnp = fft_box_size_CD(3)
#ifdef _MPIFFT_
       idp  = fft_box_size_CD_c(1,0)
#else
       idp  = fft_box_size_CD(1,0)
#endif
       nlp  = fft_box_size_CD(1,1)
       nmp  = fft_box_size_CD(2,1)
       nnp  = fft_box_size_CD(3,1)
!            allocate(chg_rs(0:ie,0:je,0:ke))
       allocate(chg_rs(0:nlp,0:nmp,0:nnp))
       chg_rs=0.d0
       unit=unusedUnit()
       open(unit,file='nfchgt_paw.data',status='unknown',form='formatted')
    end if

    do iloop=1,nspin
       call map_charge_onto_a_fft_box_paw('soft')
       call m_FFT_CD_inverse0(nfout,afft)
!!$            call m_FFT_check_of_negative_CD(afft,nfout,nspin,iloop)
!!$            call m_FFT_check_of_negative_CD(npes,ista_fftp,iend_fftp&
!!$               & ,ista_fftph,iend_fftph,afft,nfout,nspin,iloop)
       if(mype == 0) then
          call m_FFT_cp_afft_CD(afft_mpi1)
          call get_cd_in_rspace
       end if
    end do
!call tmp
!stop
    if(mype==0) then
       call output_box_tecplot(chg_rs,altv,'chgsoft','cd',unit,0)
       deallocate(chg_rs)
    end if

    call m_FFT_dealloc_CD_box()
    deallocate(afft)
!!$deallocate(afft_mpi1)
    if(npes >= 2) then
       deallocate(afft_mpi2);deallocate(afft_mpi3)
    end if

    if(mype==0) then
!            call wd_ae_cd(unit)
!            call wd_averaged_ae_cd(unit)
!            call wd_ps_cd(unit)
!            call wd_averaged_ps_cd(unit)
!            call wd_valence_cd(unit)
!            call wd_averaged_valence_cd(unit)
       if(flg_symmtry) then
          call wd_averaged_ps_cd(unit)
!                call wd_averaged_valence_cd(unit)
       else
          call wd_ps_cd(unit)
!                call wd_valence_cd(unit)
       end if
       close(unit)
    end if

    return

  contains

    subroutine map_charge_onto_a_fft_box_paw(s_or_a)
      character(len=*),intent(in):: s_or_a
      integer:: i,j,ip, ilast

      ilast = min(iend_kngp,kgp_reduced)
      afft_mpi1=0.d0
      select case(s_or_a)
      case('all')
         do j=1,kimg
!!$                    do i=ista_kngp,iend_kngp
            do i=ista_kngp,ilast
               ip=(igfp_l(i)-1)*kimg+j
               afft_mpi1(ip)=afft_mpi1(ip)+chgq_l(i,j,iloop)
            end do
         end do
      case('soft')
         do j=1,kimg
!!$                    do i=ista_kngp,iend_kngp
            do i=ista_kngp,ilast
               ip=(igfp_l(i)-1)*kimg+j
               afft_mpi1(ip)=afft_mpi1(ip)+chgsoft(i,j,iloop)
            end do
         end do
      case('hard')
         do j=1,kimg
!!$                    do i=ista_kngp,iend_kngp
            do i=ista_kngp,ilast
               ip=(igfp_l(i)-1)*kimg+j
               afft_mpi1(ip)=afft_mpi1(ip)+chgq_l(i,j,iloop)-chgsoft(i,j,iloop)
            end do
         end do
      case default
         write(6,*) 'Error in map_charge_onto_a_fft_box_paw !'
         write(6,*) 'Unknown a_or_s : ',s_or_a
         call phase_error_with_msg(nfout,'Error in map_charge_onto_a_fft_box_paw !',__LINE__,__FILE__)
      end select

      if(npes >= 2) then
         call mpi_barrier(MPI_CommGroup,ierr)
         do j=0,npes-1
            do i=nis_fftp(j),nie_fftp(j)
               afft_mpi2(i-nis_fftp(j)+1) = afft_mpi1(i)
            end do
            call mpi_allreduce(afft_mpi2,afft_mpi3,mp_fftp &
                 ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
            if(j==mype) then
               do i=ista_fftp,iend_fftp
                  afft(i)=afft_mpi3(i-ista_fftp+1)
               end do
            end if
         end do
      else
         afft=afft_mpi1
      end if
    end subroutine map_charge_onto_a_fft_box_paw

    subroutine get_cd_in_rspace
      integer:: i,j,k,ip
      integer:: ii,jj,kk
!real(8):: vec(3)
!            do k=1,ke
!                do j=1,je
!	                do i=1,ie
      do k=1,nnp
         do j=1,nmp
            do i=1,nlp
!                        ip = (ie+1)*je*(k-1) + (ie+1)*(j-1) + i
               ii=i
               jj=j
               kk=k
               if(kimg==1) then
                  if(2*(i-1).le.nlp) then
                     ii=2*i-1
                  else
                     ii=2*(-i+nlp+2)-1
                     jj=-j+nmp+2;if(j.eq.1) jj=1
                     kk=-k+nnp+2;if(k.eq.1) kk=1
                  end if
               end if
               ip = idp*nmp*(kk-1) + idp*(jj-1) + ii
               chg_rs(i-1,j-1,k-1) = chg_rs(i-1,j-1,k-1)+ &
                    afft_mpi1(ip*kimg+1-kimg)
!                                                afft_mpi1(2*ip-1)
!vec=(/dble(i)/dble(ie),dble(j)/dble(je),dble(k)/dble(ke)/)
!chg_rs(i-1,j-1,k-1)=getFFTValue(vec)
            end do
         end do
      end do
!            chg_rs(ie,:,:)=chg_rs(0,:,:)
!            chg_rs(:,je,:)=chg_rs(:,0,:)
!            chg_rs(:,:,ke)=chg_rs(:,:,0)
      chg_rs(nlp,:,:)=chg_rs(0,:,:)
      chg_rs(:,nmp,:)=chg_rs(:,0,:)
      chg_rs(:,:,nnp)=chg_rs(:,:,0)
      return
    end subroutine get_cd_in_rspace

    subroutine tmp
      integer:: idp,nlp,nmp,nnp
      integer:: i,j,k,n,nn


!!$            idp = fft_box_size_CD(0)
!!$            nlp = fft_box_size_CD(1)
!!$            nmp = fft_box_size_CD(2)
!!$            nnp = fft_box_size_CD(3)
#ifdef _MPIFFT_
      idp  = fft_box_size_CD_c(1,0)
#else
      idp  = fft_box_size_CD(1,0)
#endif
      nlp  = fft_box_size_CD(1,1)
      nmp  = fft_box_size_CD(2,1)
      nnp  = fft_box_size_CD(3,1)

      do n = 1,nfftp
         nn = (n+kimg-1)/kimg
         i  = mod(nn,idp)
         j  = mod((nn-1)/idp,nmp) + 1
         k  = (nn - (j-1)*idp - i)/(idp*nmp) + 1
!!$                print *,i,j,k,n,afft_mpi1(n)
      end do
      return
    end subroutine tmp

    function getFFTValue(vec,mode)
      real(8):: getFFTValue
      real(8),dimension(3),intent(inout):: vec
      character(len=*),intent(in):: mode
      integer:: i,is
      real(8):: value,phi,kk1,kk2,kk3

      value=0.d0
      select case(mode)
      case('soft')
         if(kimg==1) then
            do i=1,kgp
               kk1=dble(ngabc(i,1))
               kk2=dble(ngabc(i,2))
               kk3=dble(ngabc(i,3))
               phi=PAI2*(kk1*vec(1)+ &
                    kk2*vec(2)+ &
                    kk3*vec(3))
               do is=1,nspin
                  value=value+chgsoft(i,1,is)*cos(phi)
               end do
            end do
         else if(kimg==2) then
            do i=1,kgp
               kk1=dble(ngabc(i,1))
               kk2=dble(ngabc(i,2))
               kk3=dble(ngabc(i,3))
               phi=PAI2*(kk1*vec(1)+ &
                    kk2*vec(2)+ &
                    kk3*vec(3))
               do is=1,nspin
                  value=value+chgsoft(i,1,is)*cos(phi)- &
                       chgsoft(i,2,is)*sin(phi)
               end do
            end do
         end if
      end select
      getFFTValue=value
      return
    end function getFFTValue

    subroutine wd_ae_cd(unit)
      integer,intent(in):: unit
      integer:: ia,it,ith,iph,nrc,ir,is,ip
      integer:: lmt1,lmt2,il1,il2,im1,im2,it1,it2,ii,jj
      real(DP):: sint,cost,phi
      real(DP):: x0,y0,z0,xorg,yorg,zorg,pos1,pos2,pos3
      real(DP):: sum,fac,rr
      character(len=32):: title
      do ia=1,natm
         it=ityp(ia)
         if(ipaw(it)/=1) cycle
         nrc=wf_mnrc(it)
         pos1=pos(ia,1);if(pos1<0) pos1=pos1+1.d0
         pos2=pos(ia,2);if(pos2<0) pos2=pos2+1.d0
         pos3=pos(ia,3);if(pos3<0) pos3=pos3+1.d0
         xorg=altv(1,1)*pos1+altv(1,2)*pos2+altv(1,3)*pos3
         yorg=altv(2,1)*pos1+altv(2,2)*pos2+altv(2,3)*pos3
         zorg=altv(3,1)*pos1+altv(3,2)*pos2+altv(3,3)*pos3
         write(title,'("nae",i4.4)') ia
         write(unit,'(3a,i6,a,i6,a,i6,a)') &
              'ZONE T=',trim(title),'  I=',nrc &
              ,', J=',nphi+1 &
              ,', K=',ntheta &
              ,', F=POINT'
         do ith=1,ntheta
            cost=cos_theta(ith)
            sint=sqrt(1.d0-cost**2)
            do ip=1,nphi+1
               iph=ip
               if(ip.eq.nphi+1) iph=1
               phi=PAI2/dble(nphi)*dble(iph-1)
               x0=sint*cos(phi)
               y0=sint*sin(phi)
               z0=cost

               do ir=1,nrc
                  sum=0.d0
                  do is=1,nspin
                     do lmt1=1,ilmt(it)
                        il1=ltp(lmt1,it)
                        im1=mtp(lmt1,it)
                        it1=taup(lmt1,it)
                        ii=(il1-1)**2+im1
                        do lmt2=lmt1,ilmt(it)
                           il2=ltp(lmt2,it)
                           im2=mtp(lmt2,it)
                           it2=taup(lmt2,it)
                           jj=(il2-1)**2+im2
                           fac=2.d0;if(lmt1.eq.lmt2) fac=1.d0
                           sum=sum+fac*hsr(ia,lmt1,lmt2,is)* &
                                psirpw(ir,il1,it1,it)* &
                                psirpw(ir,il2,it2,it)* &
                                ylm(ith,iph,ii)* &
                                ylm(ith,iph,jj)
                        end do
                     end do
                  end do
                  rr=radr_paw(ir,it)
                  write(unit,'(4e19.6)') xorg+x0*rr, &
                       yorg+y0*rr, &
                       zorg+z0*rr,sum/rr/rr
               end do
            end do
         end do
      end do

      return
    end subroutine wd_ae_cd

    subroutine wd_ps_cd(unit)
      integer,intent(in):: unit
      integer:: ia,it,ith,iph,nrc,ir,is,ip
      integer:: lmt1,lmt2,il1,il2,im1,im2,it1,it2,ii,jj
      real(DP):: sint,cost,phi
      real(DP):: x0,y0,z0,xorg,yorg,zorg,pos1,pos2,pos3
      real(DP):: sum,fac,rr
      character(len=32):: title
      do ia=1,natm
         it=ityp(ia)
         if(ipaw(it)/=1) cycle
         nrc=wf_mnrc(it)
         pos1=pos(ia,1);if(pos1<0) pos1=pos1+1.d0
         pos2=pos(ia,2);if(pos2<0) pos2=pos2+1.d0
         pos3=pos(ia,3);if(pos3<0) pos3=pos3+1.d0
         xorg=altv(1,1)*pos1+altv(1,2)*pos2+altv(1,3)*pos3
         yorg=altv(2,1)*pos1+altv(2,2)*pos2+altv(2,3)*pos3
         zorg=altv(3,1)*pos1+altv(3,2)*pos2+altv(3,3)*pos3
         write(title,'("nps",i4.4)') ia
         write(unit,'(3a,i6,a,i6,a,i6,a)') &
              'ZONE T=',trim(title),'  I=',nrc &
              ,', J=',nphi+1 &
              ,', K=',ntheta &
              ,', F=POINT'
         do ith=1,ntheta
            cost=cos_theta(ith)
            sint=sqrt(1.d0-cost**2)
            do ip=1,nphi+1
               iph=ip
               if(ip.eq.nphi+1) iph=1
               phi=PAI2/dble(nphi)*dble(iph-1)
               x0=sint*cos(phi)
               y0=sint*sin(phi)
               z0=cost

               do ir=1,nrc
                  sum=0.d0
                  do is=1,nspin
                     do lmt1=1,ilmt(it)
                        il1=ltp(lmt1,it)
                        im1=mtp(lmt1,it)
                        it1=taup(lmt1,it)
                        ii=(il1-1)**2+im1
                        do lmt2=lmt1,ilmt(it)
                           il2=ltp(lmt2,it)
                           im2=mtp(lmt2,it)
                           it2=taup(lmt2,it)
                           jj=(il2-1)**2+im2
                           fac=2.d0;if(lmt1.eq.lmt2) fac=1.d0
                           sum=sum+fac*hsr(ia,lmt1,lmt2,is)* &
                                phirpw(ir,il1,it1,it)* &
                                phirpw(ir,il2,it2,it)* &
                                ylm(ith,iph,ii)* &
                                ylm(ith,iph,jj)
                        end do
                     end do
                  end do
                  rr=radr_paw(ir,it)
                  write(unit,'(4e19.6)') xorg+x0*rr, &
                       yorg+y0*rr, &
                       zorg+z0*rr,sum/rr/rr
               end do
            end do
         end do
      end do

      return
    end subroutine wd_ps_cd

    subroutine wd_valence_cd(unit)
      integer,intent(in):: unit
      integer:: ia,it,ith,iph,nrc,ir,is,ip
      integer:: lmt1,lmt2,il1,il2,im1,im2,it1,it2,ii,jj
      real(DP):: sint,cost,phi
      real(DP):: x0,y0,z0,xorg,yorg,zorg,pos1,pos2,pos3
      real(DP):: sum1,sum2,fac,rr,vec(3)
      character(len=32):: title
      do ia=1,natm
         it=ityp(ia)
         if(ipaw(it)/=1) cycle
         nrc=wf_mnrc(it)
         pos1=pos(ia,1);if(pos1<0) pos1=pos1+1.d0
         pos2=pos(ia,2);if(pos2<0) pos2=pos2+1.d0
         pos3=pos(ia,3);if(pos3<0) pos3=pos3+1.d0
         xorg=altv(1,1)*pos1+altv(1,2)*pos2+altv(1,3)*pos3
         yorg=altv(2,1)*pos1+altv(2,2)*pos2+altv(2,3)*pos3
         zorg=altv(3,1)*pos1+altv(3,2)*pos2+altv(3,3)*pos3
         write(title,'("chgval",i4.4)') ia
         write(unit,'(3a,i6,a,i6,a,i6,a)') &
              'ZONE T=',trim(title),'  I=',nrc &
              ,', J=',nphi+1 &
              ,', K=',ntheta &
              ,', F=POINT'
         do ith=1,ntheta
            cost=cos_theta(ith)
            sint=sqrt(1.d0-cost**2)
            do ip=1,nphi+1
               iph=ip
               if(ip.eq.nphi+1) iph=1
               phi=PAI2/dble(nphi)*dble(iph-1)
               x0=sint*cos(phi)
               y0=sint*sin(phi)
               z0=cost

               do ir=1,nrc
                  sum1=0.d0
                  sum2=0.d0
                  do is=1,nspin
                     do lmt1=1,ilmt(it)
                        il1=ltp(lmt1,it)
                        im1=mtp(lmt1,it)
                        it1=taup(lmt1,it)
                        ii=(il1-1)**2+im1
                        do lmt2=lmt1,ilmt(it)
                           il2=ltp(lmt2,it)
                           im2=mtp(lmt2,it)
                           it2=taup(lmt2,it)
                           jj=(il2-1)**2+im2
                           fac=2.d0;if(lmt1.eq.lmt2) fac=1.d0
                           sum1=sum1+fac*hsr(ia,lmt1,lmt2,is)* &
                                psirpw(ir,il1,it1,it)* &
                                psirpw(ir,il2,it2,it)* &
                                ylm(ith,iph,ii)* &
                                ylm(ith,iph,jj)
                           sum2=sum2+fac*hsr(ia,lmt1,lmt2,is)* &
                                phirpw(ir,il1,it1,it)* &
                                phirpw(ir,il2,it2,it)* &
                                ylm(ith,iph,ii)* &
                                ylm(ith,iph,jj)
                        end do
                     end do
                  end do
                  rr=radr_paw(ir,it)
                  vec=(/xorg+x0*rr,yorg+y0*rr,zorg+z0*rr/)
                  vec=matmul(rltv,vec)/PAI2
                  write(unit,'(4e19.6)') xorg+x0*rr, &
                       yorg+y0*rr, &
                       zorg+z0*rr, &
                       (sum1-sum2)/rr/rr+ &
                       getFFTValue(vec,'soft')
               end do
            end do
         end do
      end do

      return
    end subroutine wd_valence_cd

    subroutine wd_averaged_ae_cd(unit)
      integer,intent(in):: unit
      integer:: ia,it,ith,iph,nrc,ir,is,ip,iopr,ja
      integer:: lmt1,lmt2,il1,il2,im1,im2,it1,it2,ii,jj
      real(DP):: sint,cost,phi
      real(DP):: x0,y0,z0,xorg,yorg,zorg,pos1,pos2,pos3
      real(DP):: sum,fac,rr,vec0(3),sum2
      real(DP),pointer,dimension(:,:):: ylm_sym
      character(len=32):: title

      allocate(ylm_sym(nopr+af,25))
      ylm_sym=0.d0

      do ia=1,natm
         it=ityp(ia)
         if(ipaw(it)/=1) cycle
         nrc=wf_mnrc(it)
         pos1=pos(ia,1);if(pos1<0) pos1=pos1+1.d0
         pos2=pos(ia,2);if(pos2<0) pos2=pos2+1.d0
         pos3=pos(ia,3);if(pos3<0) pos3=pos3+1.d0
         xorg=altv(1,1)*pos1+altv(1,2)*pos2+altv(1,3)*pos3
         yorg=altv(2,1)*pos1+altv(2,2)*pos2+altv(2,3)*pos3
         zorg=altv(3,1)*pos1+altv(3,2)*pos2+altv(3,3)*pos3
         write(title,'("nae_sym",i4.4)') ia
         write(unit,'(3a,i6,a,i6,a,i6,a)') &
              'ZONE T=',trim(title),'  I=',nrc &
              ,', J=',nphi+1 &
              ,', K=',ntheta &
              ,', F=POINT'
         do ith=1,ntheta
            cost=cos_theta(ith)
            sint=sqrt(1.d0-cost**2)
            do ip=1,nphi+1
               iph=ip
               if(ip.eq.nphi+1) iph=1
               phi=PAI2/dble(nphi)*dble(iph-1)
               x0=sint*cos(phi)
               y0=sint*sin(phi)
               z0=cost
               vec0(1)=x0
               vec0(2)=y0
               vec0(3)=z0
               call set_ylms_symmtry_op(ia,25,vec0,ylm_sym)

               do ir=1,nrc
                  sum=0.d0
                  do is=1,nspin
                     do lmt1=1,ilmt(it)
                        il1=ltp(lmt1,it)
                        im1=mtp(lmt1,it)
                        it1=taup(lmt1,it)
                        ii=(il1-1)**2+im1
                        do lmt2=lmt1,ilmt(it)
                           il2=ltp(lmt2,it)
                           im2=mtp(lmt2,it)
                           it2=taup(lmt2,it)
                           jj=(il2-1)**2+im2
                           fac=2.d0;if(lmt1.eq.lmt2) fac=1.d0

                           sum2=0.d0
                           do iopr=1,nopr
                              ja=abs(ia2ia_symmtry_op(ia,iopr))
                              sum2=sum2+hsr(ja,lmt1,lmt2,is)* &
                                   ylm_sym(iopr,ii)* &
                                   ylm_sym(iopr,jj)
                           end do
                           sum2=sum2/dble(nopr)

                           sum=sum+fac* &
                                psirpw(ir,il1,it1,it)* &
                                psirpw(ir,il2,it2,it)* &
                                sum2
                        end do
                     end do
                  end do
                  rr=radr_paw(ir,it)
                  write(unit,'(4e19.6)') xorg+x0*rr, &
                       yorg+y0*rr, &
                       zorg+z0*rr,sum/rr/rr
               end do
            end do
         end do
      end do

      deallocate(ylm_sym)

      return
    end subroutine wd_averaged_ae_cd

    subroutine wd_averaged_ps_cd(unit)
      integer,intent(in):: unit
      integer:: ia,it,ith,iph,nrc,ir,is,ip,iopr,ja
      integer:: lmt1,lmt2,il1,il2,im1,im2,it1,it2,ii,jj
      real(DP):: sint,cost,phi
      real(DP):: x0,y0,z0,xorg,yorg,zorg,pos1,pos2,pos3
      real(DP):: sum,fac,rr,vec0(3),sum2
      real(DP),pointer,dimension(:,:):: ylm_sym
      character(len=32):: title

      allocate(ylm_sym(nopr+af,25))
      ylm_sym=0.d0

      do ia=1,natm
         it=ityp(ia)
         if(ipaw(it)/=1) cycle
         nrc=wf_mnrc(it)
         pos1=pos(ia,1);if(pos1<0) pos1=pos1+1.d0
         pos2=pos(ia,2);if(pos2<0) pos2=pos2+1.d0
         pos3=pos(ia,3);if(pos3<0) pos3=pos3+1.d0
         xorg=altv(1,1)*pos1+altv(1,2)*pos2+altv(1,3)*pos3
         yorg=altv(2,1)*pos1+altv(2,2)*pos2+altv(2,3)*pos3
         zorg=altv(3,1)*pos1+altv(3,2)*pos2+altv(3,3)*pos3
         write(title,'("nps_sym",i4.4)') ia
         write(unit,'(3a,i6,a,i6,a,i6,a)') &
              'ZONE T=',trim(title),'  I=',nrc &
              ,', J=',nphi+1 &
              ,', K=',ntheta &
              ,', F=POINT'
         do ith=1,ntheta
            cost=cos_theta(ith)
            sint=sqrt(1.d0-cost**2)
            do ip=1,nphi+1
               iph=ip
               if(ip.eq.nphi+1) iph=1
               phi=PAI2/dble(nphi)*dble(iph-1)
               x0=sint*cos(phi)
               y0=sint*sin(phi)
               z0=cost
               vec0(1)=x0
               vec0(2)=y0
               vec0(3)=z0
               call set_ylms_symmtry_op(ia,25,vec0,ylm_sym)

               do ir=1,nrc
                  sum=0.d0
                  do is=1,nspin
                     do lmt1=1,ilmt(it)
                        il1=ltp(lmt1,it)
                        im1=mtp(lmt1,it)
                        it1=taup(lmt1,it)
                        ii=(il1-1)**2+im1
                        do lmt2=lmt1,ilmt(it)
                           il2=ltp(lmt2,it)
                           im2=mtp(lmt2,it)
                           it2=taup(lmt2,it)
                           jj=(il2-1)**2+im2
                           fac=2.d0;if(lmt1.eq.lmt2) fac=1.d0

                           sum2=0.d0
                           do iopr=1,nopr
                              ja=abs(ia2ia_symmtry_op(ia,iopr))
                              sum2=sum2+hsr(ja,lmt1,lmt2,is)* &
                                   ylm_sym(iopr,ii)* &
                                   ylm_sym(iopr,jj)
                           end do
                           sum2=sum2/dble(nopr)

                           sum=sum+fac* &
                                phirpw(ir,il1,it1,it)* &
                                phirpw(ir,il2,it2,it)* &
                                sum2
                        end do
                     end do
                  end do
                  rr=radr_paw(ir,it)
                  write(unit,'(4e19.6)') xorg+x0*rr, &
                       yorg+y0*rr, &
                       zorg+z0*rr,sum/rr/rr
               end do
            end do
         end do
      end do

      deallocate(ylm_sym)
      return
    end subroutine wd_averaged_ps_cd

    subroutine wd_averaged_valence_cd(unit)
      integer,intent(in):: unit
      integer:: ia,it,ith,iph,nrc,ir,is,ip,iopr,ja
      integer:: lmt1,lmt2,il1,il2,im1,im2,it1,it2,ii,jj
      real(DP):: sint,cost,phi
      real(DP):: x0,y0,z0,xorg,yorg,zorg,pos1,pos2,pos3
      real(DP):: sum1,sum2,fac,rr,vec(3),vec0(3),sum3
      real(DP),pointer,dimension(:,:):: ylm_sym
      character(len=32):: title

      allocate(ylm_sym(nopr+af,25))
      ylm_sym=0.d0

      do ia=1,natm
         it=ityp(ia)
         if(ipaw(it)/=1) cycle
         nrc=wf_mnrc(it)
         pos1=pos(ia,1);if(pos1<0) pos1=pos1+1.d0
         pos2=pos(ia,2);if(pos2<0) pos2=pos2+1.d0
         pos3=pos(ia,3);if(pos3<0) pos3=pos3+1.d0
         xorg=altv(1,1)*pos1+altv(1,2)*pos2+altv(1,3)*pos3
         yorg=altv(2,1)*pos1+altv(2,2)*pos2+altv(2,3)*pos3
         zorg=altv(3,1)*pos1+altv(3,2)*pos2+altv(3,3)*pos3
         write(title,'("chgval_sym",i4.4)') ia
         write(unit,'(3a,i6,a,i6,a,i6,a)') &
              'ZONE T=',trim(title),'  I=',nrc &
              ,', J=',nphi+1 &
              ,', K=',ntheta &
              ,', F=POINT'
         do ith=1,ntheta
            cost=cos_theta(ith)
            sint=sqrt(1.d0-cost**2)
            do ip=1,nphi+1
               iph=ip
               if(ip.eq.nphi+1) iph=1
               phi=PAI2/dble(nphi)*dble(iph-1)
               x0=sint*cos(phi)
               y0=sint*sin(phi)
               z0=cost
               vec0(1)=x0
               vec0(2)=y0
               vec0(3)=z0
               call set_ylms_symmtry_op(ia,25,vec0,ylm_sym)

               do ir=1,nrc
                  sum1=0.d0
                  sum2=0.d0
                  do is=1,nspin
                     do lmt1=1,ilmt(it)
                        il1=ltp(lmt1,it)
                        im1=mtp(lmt1,it)
                        it1=taup(lmt1,it)
                        ii=(il1-1)**2+im1
                        do lmt2=lmt1,ilmt(it)
                           il2=ltp(lmt2,it)
                           im2=mtp(lmt2,it)
                           it2=taup(lmt2,it)
                           jj=(il2-1)**2+im2
                           fac=2.d0;if(lmt1.eq.lmt2) fac=1.d0

                           sum3=0.d0
                           do iopr=1,nopr
                              ja=abs(ia2ia_symmtry_op(ia,iopr))
                              sum3=sum3+hsr(ja,lmt1,lmt2,is)* &
                                   ylm_sym(iopr,ii)* &
                                   ylm_sym(iopr,jj)
                           end do
                           sum3=sum3/dble(nopr)

                           sum1=sum1+fac* &
                                psirpw(ir,il1,it1,it)* &
                                psirpw(ir,il2,it2,it)* &
                                sum3
                           sum2=sum2+fac* &
                                phirpw(ir,il1,it1,it)* &
                                phirpw(ir,il2,it2,it)* &
                                sum3
                        end do
                     end do
                  end do
                  rr=radr_paw(ir,it)
                  vec=(/xorg+x0*rr,yorg+y0*rr,zorg+z0*rr/)
                  vec=matmul(rltv,vec)/PAI2
                  write(unit,'(4e19.6)') xorg+x0*rr, &
                       yorg+y0*rr, &
                       zorg+z0*rr, &
                       (sum1-sum2)/rr/rr+ &
                       getFFTValue(vec,'soft')
               end do
            end do
         end do
      end do

      return
    end subroutine wd_averaged_valence_cd

  end subroutine m_PAWCD_wd_cd

  subroutine set_ia2ia_symmtry_op
    integer:: ia,it,no,ja,jt
    integer:: i,j,k
    real(DP):: pos0(3),pos1(3),pos2(3),pos3(3)
    real(DP):: distance

    allocate(op_pr(3,3,nopr+af))
    call m_CS_op_in_PUCD(nfout,op_pr,nopr+af)

    do ia=1,natm
       it=ityp(ia)
       if(ipaw(it)/=1) cycle
       pos0(1:3)=pos(ia,1:3)
!print *,'pos0=',pos0
       pos0(:) = pos0(:) - floor(pos0(:))
       do no=1,nopr
          pos1(:)=matmul(op_pr(:,:,no),pos0(:))+tau(:,no,BUCS)

          pos1(:) = pos1(:) - floor(pos1(:))
          KLoop: do k=-1,1
             do j=-1,1
                do i=-1,1
                   do ja=1,natm
                      jt=ityp(ja)
                      if(it/=jt) cycle
                      pos3(1:3) = pos(ja,1:3)
                      pos3(:) = pos3(:) - floor(pos3(:))
                      pos2(1:3)=pos3(1:3)+(/dble(i),dble(j),dble(k)/)

                      distance=abs(pos1(1)-pos2(1))+abs(pos1(2)-pos2(2)) &
                           +abs(pos1(3)-pos2(3))
                      if(distance < 1.d-5) then
                         ia2ia_symmtry_op(ia,no)=ja

                         exit KLoop
                      end if

                      if(kimg==1 .and. iwei(ja)==2) then
                         pos2(1:3)=-pos3(1:3)+(/dble(i),dble(j),dble(k)/)
                         distance=abs(pos1(1)-pos2(1))+abs(pos1(2)-pos2(2)) &
                              +abs(pos1(3)-pos2(3))
                         if(distance < 1.d-5) then
                            ia2ia_symmtry_op(ia,no)=-ja

                            exit KLoop
                         end if
                      end if

                      if(i.eq.1 .and. j.eq.1 .and. k.eq.1 .and. ja==natm) then
                         write(nfout,*) 'Error in set_ia2ia_symmtry_op (in m_PAW_ChargeDensity) !'
                         write(nfout,*) 'ia,no =', ia,no
                         write(nfout,'(" pos(ia)                     = ",3f10.4)') pos(ia,1:3)
                         write(nfout,'(" pos1(ia) = op*pos(ia)+tau   = ",3f10.4)') pos1(1:3)
                         call phase_error_with_msg(nfout,'Error in set_ia2ia_symmtry_op (in m_PAW_ChargeDensity) !'&
                                                  ,__LINE__,__FILE__)
                      end if

                   end do
                end do
             end do
          end do KLoop
       end do
    end do

    deallocate(op_pr)
    return
  end subroutine set_ia2ia_symmtry_op

  subroutine set_ylms_symmtry_op(ia,ilmax,vec0,ylm_sym)
    integer,intent(in):: ia,ilmax
    real(DP),intent(in):: vec0(3)
    real(DP),intent(out):: ylm_sym(nopr+af,ilmax)

    integer:: is,iopr
    real(DP),pointer,dimension(:)     :: gx,gy,gz,tylm
    real(DP):: vec(3)

    allocate(gx(nopr+af))
    allocate(gy(nopr+af))
    allocate(gz(nopr+af))
    allocate(tylm(nopr+af))

    do iopr=1,nopr
       vec(:)=matmul(op(:,:,iopr),vec0(:))
       if(ia2ia_symmtry_op(ia,iopr) < 0) vec=-vec
       gx(iopr)=vec(1)
       gy(iopr)=vec(2)
       gz(iopr)=vec(3)
    end do

    do is=1,ilmax
       call sphr(nopr,is,gx(1:nopr),gy(1:nopr),gz(1:nopr),tylm(1:nopr))
       do iopr=1,nopr
          ylm_sym(iopr,is)=tylm(iopr)
       end do
    end do

    deallocate(gx,gy,gz,tylm)
    return
  end subroutine set_ylms_symmtry_op

  subroutine set_diff_ylms_symmtry_op(ia,ilmax,vec0,ylm_sym,dylm_dth_sym,dylm_dph_sym)
    integer,intent(in):: ia,ilmax
    real(DP),intent(in):: vec0(3)
    real(DP),intent(out):: ylm_sym(nopr+af,ilmax)
    real(DP),intent(out):: dylm_dth_sym(nopr+af,ilmax)
    real(DP),intent(out):: dylm_dph_sym(nopr+af,ilmax)

    integer:: is,iopr
    real(DP),pointer,dimension(:)     :: gx,gy,gz,tylm
    real(DP),pointer,dimension(:,:)   :: tdylm
    real(DP):: vec(3),rr,vec1(3),vec2(3)

    allocate(gx(nopr+af))
    allocate(gy(nopr+af))
    allocate(gz(nopr+af))
    allocate(tylm(nopr+af))
    allocate(tdylm(nopr+af,3))

    do iopr=1,nopr
       vec(:)=matmul(op(:,:,iopr),vec0(:))
       if(ia2ia_symmtry_op(ia,iopr) < 0) vec=-vec
       gx(iopr)=vec(1)
       gy(iopr)=vec(2)
       gz(iopr)=vec(3)
    end do

    do is=1,ilmax
       call sphr(nopr,is,gx(1:nopr),gy(1:nopr),gz(1:nopr),tylm(1:nopr))
       do iopr=1,nopr
          ylm_sym(iopr,is)=tylm(iopr)
       end do
    end do

    rr=sqrt(vec0(1)**2+vec0(2)**2)
    vec1(1)=vec0(3)*vec0(1)/rr
    vec1(2)=vec0(2)*vec0(3)/rr
    vec1(3)=-rr
    vec2(1)=-vec0(2)
    vec2(2)=vec0(1)
    vec2(3)=0.d0

    do is=1,ilmax
       call sphr_diff(nopr,nopr,is,gx(1:nopr),gy(1:nopr),gz(1:nopr),tdylm(1:nopr,1:3))
       do iopr=1,nopr
          vec=matmul(op(:,:,iopr),vec1(:))
          dylm_dth_sym(iopr,is)=  tdylm(iopr,1)*vec(1)+ &
               tdylm(iopr,2)*vec(2)+ &
               tdylm(iopr,3)*vec(3)
          vec=matmul(op(:,:,iopr),vec2(:))
          dylm_dph_sym(iopr,is)=  tdylm(iopr,1)*vec(1)+ &
               tdylm(iopr,2)*vec(2)+ &
               tdylm(iopr,3)*vec(3)

          if(ia2ia_symmtry_op(ia,iopr) < 0) then
             dylm_dth_sym(iopr,is)=-dylm_dth_sym(iopr,is)
             dylm_dph_sym(iopr,is)=-dylm_dph_sym(iopr,is)
          end if
          !                ylm_sym(iopr,is)=tylm(iopr)
       end do
    end do

    deallocate(gx,gy,gz,tylm,tdylm)
    return
  end subroutine set_diff_ylms_symmtry_op

  subroutine m_PAWCD_wd_number_of_state(nfout)
    integer,intent(in):: nfout

    integer:: ia,ja,it,ith,iph
    integer:: ir,is,ier,nrc
    real(DP):: costh,sinth,phi,dvec(3)
    real(DP),pointer,dimension(:,:):: nae
    real(DP),pointer,dimension(:):: wos
    real(DP),pointer,dimension(:,:):: nsae
    logical,pointer,dimension(:):: flg_done

    allocate(nae(mmesh,nspin))
    allocate(wos(mmesh));wos=0.d0
    allocate(nsae(natm,nspin))
    if(af /= 0) then
       allocate(flg_done(natm))
       flg_done=.false.
    end if
    nsae=0.d0

    do ia=1,natm
       if(af /= 0) then
          if(flg_done(ia)) cycle
       end if
       it=ityp(ia)
       if(ipaw(it)/=1) cycle
       nrc=wf_mnrc(it)
       call set_weight_exp(ier,1,nrc,radr_paw(:,it),wos)
       do ir=1,nrc
          wos(ir)=wos(ir)*radr_paw(ir,it)**2
       end do
!            wos=wos*iwei(ia)
       do ith=1,ntheta
          costh=cos_theta(ith)
          sinth=sqrt(1.d0-costh**2)
          do iph=1,nphi
             phi=PAI2/dble(nphi)*dble(iph-1)
             dvec(1)=sinth*cos(phi)
             dvec(2)=sinth*sin(phi)
             dvec(3)=costh
             call m_PAWCD_set_ae_vd_sym &
                  (dvec,ia,nspin,nrc,nae(1:nrc,1:nspin))

             do is=1,nspin
                do ir=1,nrc
                   nsae(ia,is)=nsae(ia,is)+ &
                        nae(ir,is)*omg_wght(ith)*wos(ir)
                end do
             end do
          end do
       end do

       if(af /= 0) then
          ja=ia2ia_symmtry_op(ia,nopr+af)
          flg_done(ja)=.true.
       end if
    end do
    nsae=nsae*PAI4

    if(printable)then
       write(nfout,*) 'Site-decomposed number of states below the Fermi energy.'
       write(nfout,*) '      Site        Spin-up            Spin-down'
    endif

    flg_done=.false.

    do ia=1,natm
       if(af /= 0) then
          if(flg_done(ia)) cycle
       end if
       it=ityp(ia)
       if(ipaw(it)/=1) cycle
       if(nspin.eq.2) then
          if(printable) write(nfout,'(i10,2f19.6)')  ia,nsae(ia,1),nsae(ia,2)
       else
          if(printable) write(nfout,'(i10,2f19.6)')  ia,nsae(ia,1)/2.d0,nsae(ia,1)/2.d0
       end if
       if(af /= 0) then
          ja=ia2ia_symmtry_op(ia,nopr+af)
          if(printable) write(nfout,'(i10,2f19.6)')  ja,nsae(ia,2),nsae(ia,1)
          flg_done(ja)=.true.
       end if
    end do

    if(printable) write(nfout,*)

    deallocate(nae)
    deallocate(wos)
    deallocate(nsae)
    if(af /= 0) deallocate(flg_done)

  end subroutine m_PAWCD_wd_number_of_state


  subroutine m_PAWCD_set_ae_cd_sphex2_2D( ia,nspin,nrc,dnr &
       &                              ,msph,n1nc,msphmx &
       &                              ,num_isph,isph_chg )
    integer,intent(in) :: ia,nspin,nrc,dnr,msph
    real(8),intent(out):: n1nc(nrc,nspin,25)
    integer,intent(out):: msphmx
    integer,intent(out):: num_isph,isph_chg(25)

    real(8):: n1(nrc,nspin),sum,fac,isqrt4pi,eps,abshsr, weight
    integer:: ir,it,is,lmt1,lmt2
    integer:: il1,il2,it1,it2,im1,im2
    integer:: ii,jj,n,isp
    logical:: flg_isph(25)

    it=ityp(ia)
    n1nc=0.d0
    msphmx=0
    eps=1.d-10
    flg_isph=.false.
    do lmt1=1,ilmt(it)
       il1=ltp(lmt1,it)
       im1=mtp(lmt1,it)
       it1=taup(lmt1,it)
       ii=(il1-1)**2+im1
       do lmt2=lmt1,ilmt(it)
          abshsr=abs(hsr(ia,lmt1,lmt2,1))
          if(nspin.eq.2) abshsr=abshsr+abs(hsr(ia,lmt1,lmt2,2))
          if(abshsr .lt. eps) cycle
          il2=ltp(lmt2,it)
          im2=mtp(lmt2,it)
          it2=taup(lmt2,it)
          jj=(il2-1)**2+im2
          fac=2.d0;if(lmt1.eq.lmt2) fac=1.d0
          do n=1,il2p(lmt1,lmt2,it)
             isp=isph(lmt1,lmt2,n,it)
             if(isp.gt.msphmx) msphmx=isp
             flg_isph(isp)=.true.
             do is=1,nspin
                sum=fac*hsr(ia,lmt1,lmt2,is)*dl2p(lmt1,lmt2,n,it)

                do ir=1,nrc,dnr
                   n1nc(ir,is,isp)=n1nc(ir,is,isp)+sum* &
                        psirpw(ir,il1,it1,it)* &
                        psirpw(ir,il2,it2,it)
                end do
             end do
          end do
       end do
    end do

    num_isph=0
    isph_chg=0
    do isp=1,25
       if(flg_isph(isp)) then
          num_isph=num_isph+1
          isph_chg(num_isph)=isp
       end if
    end do

    isqrt4pi=1.d0/sqrt(PAI4)
    do is=1,nspin
       do ir=1,nrc,dnr
          n1nc(ir,is,1)=(n1nc(ir,is,1)*isqrt4pi+rhcorpw(ir,it)/PAI4/dble(nspin)) &
               /radr_paw(ir,it)/radr_paw(ir,it)
       end do
    end do

    if ( sw_opencore == ON .and. nspin == 2 ) then
       if ( has_opencore(it) == 1 ) then
          weight = 0.5d0 *mag_opencore_pol(ia,1)
          do ir=1,nrc,dnr
             n1nc(ir,1,1) = n1nc(ir,1,1) +rmag_opencore(ir,it)/PAI4 *weight &
                  &                        /radr_paw(ir,it)/radr_paw(ir,it)
             n1nc(ir,2,1) = n1nc(ir,2,1) -rmag_opencore(ir,it)/PAI4 *weight &
                  &                        /radr_paw(ir,it)/radr_paw(ir,it)
          end do
       endif
    endif

!        do n=2,msph
    do isp=2,num_isph
       n=isph_chg(isp)
       do is=1,nspin
          do ir=1,nrc,dnr
             n1nc(ir,is,n)=n1nc(ir,is,n)/radr_paw(ir,it)/radr_paw(ir,it)
          end do
       end do
    end do
!
!do ir=1,nrc
!print '(11e19.6)', radr_paw(ir,it),n1nc(ir,1,11:20)
!end do
!stop
    return
  end subroutine m_PAWCD_set_ae_cd_sphex2_2D

  subroutine m_PAWCD_set_ae_cd_sphex2( ia,nspin,nrc,dnr, msph,n1nc,msphmx &
       &                              ,num_isph,isph_chg,ista_nrc,iend_nrc,ist )
    integer,intent(in) :: ia,nspin,nrc,dnr,msph,ista_nrc,iend_nrc, ist
    real(DP),intent(out):: n1nc(ista_nrc:iend_nrc,nspin,25)
    integer,intent(out):: msphmx
    integer,intent(out):: num_isph,isph_chg(25)

    real(DP):: sum,fac,isqrt4pi,eps,abshsr, weight
    integer:: ir,it,is,lmt1,lmt2, ien
    integer:: il1,il2,it1,it2,im1,im2
    integer:: ii,jj,n,isp
    logical:: flg_isph(25)

    ien = min(iend_nrc,nrc)

    it=ityp(ia)
    n1nc=0.d0
    msphmx=0
    eps=1.d-10
    flg_isph=.false.
    do lmt1=1,ilmt(it)
       il1=ltp(lmt1,it)
       im1=mtp(lmt1,it)
       it1=taup(lmt1,it)
       ii=(il1-1)**2+im1
       do lmt2=lmt1,ilmt(it)
          abshsr=abs(hsr(ia,lmt1,lmt2,1))
          if(nspin.eq.2) abshsr=abshsr+abs(hsr(ia,lmt1,lmt2,2))
          if(abshsr .lt. eps) cycle
          il2=ltp(lmt2,it)
          im2=mtp(lmt2,it)
          it2=taup(lmt2,it)
          jj=(il2-1)**2+im2
          fac=2.d0;if(lmt1.eq.lmt2) fac=1.d0
          do n=1,il2p(lmt1,lmt2,it)
             isp=isph(lmt1,lmt2,n,it)
             if(isp.gt.msphmx) msphmx=isp
             flg_isph(isp)=.true.
             do is=1,nspin
                sum=fac*hsr(ia,lmt1,lmt2,is)*dl2p(lmt1,lmt2,n,it)

                do ir = ist, ien, dnr !   do ir=1,nrc,dnr
                   n1nc(ir,is,isp)=n1nc(ir,is,isp)+sum* &
                        psirpw(ir,il1,it1,it)* &
                        psirpw(ir,il2,it2,it)
                end do
             end do
          end do
       end do
    end do

    num_isph=0
    isph_chg=0
    do isp=1,25
       if(flg_isph(isp)) then
          num_isph=num_isph+1
          isph_chg(num_isph)=isp
       end if
    end do

    isqrt4pi=1.d0/sqrt(PAI4)
    do is=1,nspin
       do ir = ist,ien,dnr ! do ir=1,nrc,dnr
          n1nc(ir,is,1)=(n1nc(ir,is,1)*isqrt4pi+rhcorpw(ir,it)/PAI4/dble(nspin)) &
               /radr_paw(ir,it)/radr_paw(ir,it)
       end do
    end do

    if ( sw_opencore == ON .and. nspin == 2 ) then
       if ( has_opencore(it) == 1 ) then
          weight = 0.5d0 *mag_opencore_pol(ia,1)
          do ir=ist, ien, dnr
             n1nc(ir,1,1) = n1nc(ir,1,1) +rmag_opencore(ir,it)/PAI4 *weight &
                  &                        /radr_paw(ir,it)/radr_paw(ir,it)
             n1nc(ir,2,1) = n1nc(ir,2,1) -rmag_opencore(ir,it)/PAI4 *weight &
                  &                        /radr_paw(ir,it)/radr_paw(ir,it)
          end do
       endif
    endif

!        do n=2,msph
    do isp=2,num_isph
       n=isph_chg(isp)
       do is=1,nspin
          do ir=ist,ien,dnr !  do ir=1,nrc,dnr
             n1nc(ir,is,n)=n1nc(ir,is,n)/radr_paw(ir,it)/radr_paw(ir,it)
          end do
       end do
    end do
    return
  end subroutine m_PAWCD_set_ae_cd_sphex2

! ============================ added by K. Tagami ========== 11.0
  subroutine m_PAWCD_ae_cd_sphex2_nonclA( ia, nspin, nrc, dnr,  &
       &                                  msph, n1nc, msphmx, &
       &                                  num_isph, isph_chg, wos, &
       &                                  level_of_projection, &
       &                                  magmom_avg_spherical, rho_rad )
    integer, intent(in) :: ia,nspin,nrc,dnr,msph
    integer, intent(out):: msphmx
    integer, intent(out):: num_isph,isph_chg(max_sph_expansion)
    integer, intent(in) :: level_of_projection

    real(kind=DP), intent(out):: n1nc( nrc, nspin, max_sph_expansion )
    real(kind=DP), intent(in) :: wos(mmesh)
    real(kind=DP), intent(out) :: magmom_avg_spherical(3,max_sph_expansion)
    real(kind=DP), intent(out) :: rho_rad( nrc, ndim_magmom, max_sph_expansion )

    real(kind=DP) :: sum,fac,isqrt4pi,eps,abshsr, ctmp, weight

    integer:: ir,it,is,lmt1,lmt2
    integer:: il1,il2,it1,it2,im1,im2
    integer:: ii,jj,n,isp
    logical:: flg_isph( max_sph_expansion )

    it=ityp(ia)
    n1nc=0.d0; rho_rad = 0.0d0

    msphmx=0
    eps=1.d-10
    flg_isph=.false.

    do lmt1=1,ilmt(it)
       il1=ltp(lmt1,it);  im1=mtp(lmt1,it);   it1=taup(lmt1,it)

       ii=(il1-1)**2+im1

       do lmt2=lmt1,ilmt(it)
          abshsr = abs(hsr(ia,lmt1,lmt2,1))
          if(abshsr .lt. eps) cycle

          il2=ltp(lmt2,it); im2=mtp(lmt2,it);  it2=taup(lmt2,it)
          jj=(il2-1)**2+im2

          fac=2.d0;if(lmt1.eq.lmt2) fac=1.d0

          do n=1,il2p(lmt1,lmt2,it)
             isp = isph(lmt1,lmt2,n,it)
             if (isp.gt.msphmx ) msphmx=isp
             flg_isph(isp)=.true.

             do is=1,ndim_magmom
                sum = fac*hsr(ia,lmt1,lmt2,is)*dl2p(lmt1,lmt2,n,it)
                do ir=1,nrc,dnr
                   rho_rad(ir,is,isp) = rho_rad(ir,is,isp) &
                        &         +sum * psirpw(ir,il1,it1,it) &
                        &              * psirpw(ir,il2,it2,it)
                end do
             end do
          end do
       end do
    end do

    num_isph=0;   isph_chg=0

    do isp=1,max_sph_expansion
       if (flg_isph(isp)) then
          num_isph=num_isph+1
          isph_chg(num_isph)=isp
       end if
    end do

    if ( mode_of_projection == project_on_local_quantz_axis ) then
!
       select case( level_of_projection )
       case (1)
          call get_magmom_spherical( ia, nrc, dnr, rho_rad, wos, &
               &                     magmom_avg_spherical(:,1) )
          call radial_chg_along_LocalQuanzAxis( nrc, dnr, rho_rad, &
               &                                magmom_avg_spherical(:,1), n1nc )

       case(2,3)
          call get_magmom_spherical2( ia, nrc, dnr, rho_rad, wos, magmom_avg_spherical )
          call radial_chg_along_LocalQuanzAxis( nrc, dnr, rho_rad, &
               &                                magmom_avg_spherical(:,1), n1nc )

       case(4)
          call get_magmom_spherical( ia, nrc, dnr, rho_rad, wos, &
               &                     magmom_avg_spherical(:,1) )
          call radial_chg_alng_LocalQuanzAxis2( nrc, dnr, rho_rad, n1nc )

       end select

    else
       call phase_error_with_msg(nfout,'kt : under conctruction ',__LINE__,__FILE__)
    end if
!
    isqrt4pi=1.d0/sqrt(PAI4)
    do is=1,nspin
       do ir=1,nrc,dnr
          n1nc(ir,is,1) = (n1nc(ir,is,1)*isqrt4pi+rhcorpw(ir,it)/PAI4/dble(nspin)) &
               /radr_paw(ir,it)/radr_paw(ir,it)
       end do
    end do

    if ( sw_opencore == ON .and. nspin == 2 ) then
       if ( has_opencore(it) == 1 ) then
          ctmp = magmom_avg_spherical(1,1)**2 +magmom_avg_spherical(2,1)**2 &
               & +magmom_avg_spherical(3,1)**2
          ctmp = sqrt(ctmp)
          if ( ctmp > 1.0D-11 ) then
             weight = magmom_avg_spherical(1,1)/ctmp *mag_opencore_pol(ia,1) &
                  &  +magmom_avg_spherical(2,1)/ctmp *mag_opencore_pol(ia,2) &
                  &  +magmom_avg_spherical(3,1)/ctmp *mag_opencore_pol(ia,3)
             weight = 0.5d0 *weight
          else
             weight = 0.50d0
          endif

          do ir=1,nrc,dnr
             n1nc(ir,1,1) = n1nc(ir,1,1) +rmag_opencore(ir,it)/PAI4 *weight &
                  &                        /radr_paw(ir,it)/radr_paw(ir,it)
             n1nc(ir,2,1) = n1nc(ir,2,1) -rmag_opencore(ir,it)/PAI4 *weight &
                  &                        /radr_paw(ir,it)/radr_paw(ir,it)
          end do
       endif
    endif

!        do n=2,msph
    do isp=2,num_isph
       n=isph_chg(isp)
       do is=1,nspin
          do ir=1,nrc,dnr
             n1nc(ir,is,n)=n1nc(ir,is,n)/radr_paw(ir,it)/radr_paw(ir,it)
          end do
       end do
    end do

    return

  end subroutine m_PAWCD_ae_cd_sphex2_nonclA
! =========================================================== 11.0

! ============================== added by K. Tagami ================== 11.0
  subroutine radial_chg_along_LocalQuanzAxis( nrc, dnr, rho_rad, &
       &                                      magmom_avg_spherical, &
       &                                      rho_rad_along_quantz_axis )
! --
    integer, intent(in) :: nrc, dnr
    real(kind=DP), intent(in) :: rho_rad( nrc, ndim_magmom, max_sph_expansion )
    real(kind=DP), intent(in) :: magmom_avg_spherical( 3 )
    real(kind=DP), intent(out) :: &
         &                rho_rad_along_quantz_axis( nrc, nspin, max_sph_expansion )

    real(kind=DP) :: c_mx, c_my, c_mz, m_norm
    real(kind=DP) :: mmx, mmy, mmz

    real(kind=DP) :: c1, c2, d1, d2, ctmp2
    integer :: ir, isp
!    real(kind=DP) :: m_norm_criteria = 1.0D-20
    real(kind=DP) :: m_norm_criteria = 1.0D-18

    c_mx = magmom_avg_spherical(1)
    c_my = magmom_avg_spherical(2)
    c_mz = magmom_avg_spherical(3)
    m_norm = sqrt( c_mx**2 + c_my**2 + c_mz**2 )
!
    if ( m_norm > m_norm_criteria ) then
       mmx = c_mx / m_norm
       mmy = c_my / m_norm
       mmz = c_mz / m_norm
    else
       mmx = 0.0d0;  mmy = 0.0d0;  mmz = 1.0d0
       mmx = 0.0d0;  mmy = 0.0d0;  mmz = 0.0d0

!       ctmp2 =  c_mx *Global_Quantz_Axis_now(1) &
!            &  +c_my *Global_Quantz_Axis_now(2) &
!            &  +c_mz *Global_Quantz_Axis_now(3)
!
!       if ( ctmp2 > 0.0 ) then
!          mmx = Global_Quantz_Axis_now(1)
!          mmy = Global_Quantz_Axis_now(2)
!          mmz = Global_Quantz_Axis_now(3)
!       else
!          mmx = -Global_Quantz_Axis_now(1)
!          mmy = -Global_Quantz_Axis_now(2)
!          mmz = -Global_Quantz_Axis_now(3)
!       endif
    endif

    Do ir=1, nrc, dnr
       c1 = rho_rad(ir,1,1);
       c2 = rho_rad(ir,2,1) *mmx + rho_rad(ir,3,1) *mmy + rho_rad(ir,4,1) *mmz
       d1 = ( c1 + c2 ) /2.0d0;
       d2 = ( c1 - c2 ) /2.0d0;

       rho_rad_along_quantz_axis(ir,1,1) = d1
       rho_rad_along_quantz_axis(ir,2,1) = d2
    End do

    Do isp=2, max_sph_expansion
       Do ir=1, nrc, dnr
          c1 = rho_rad(ir,1,isp);
          c2 = rho_rad(ir,2,isp) *mmx +rho_rad(ir,3,isp) *mmy +rho_rad(ir,4,isp) *mmz
          d1 = ( c1 + c2 ) /2.0d0;
          d2 = ( c1 - c2 ) /2.0d0;

          rho_rad_along_quantz_axis(ir,1,isp) = d1
          rho_rad_along_quantz_axis(ir,2,isp) = d2
       End do
    End do
!
  end subroutine radial_chg_along_LocalQuanzAxis

  subroutine radial_chg_alng_LocalQuanzAxis2( nrc, dnr, rho_rad, &
       &                                      rho_rad_along_quantz_axis )
! --
    integer, intent(in) :: nrc, dnr
    real(kind=DP), intent(in) :: rho_rad( nrc, ndim_magmom, max_sph_expansion )
    real(kind=DP), intent(out) :: &
         &                rho_rad_along_quantz_axis( nrc, nspin, max_sph_expansion )

    real(kind=DP) :: c_mx, c_my, c_mz, m_norm
    real(kind=DP) :: mmx, mmy, mmz

    real(kind=DP) :: c1, c2, d1, d2, ctmp2, c3
    integer :: ir, isp
!    real(kind=DP) :: m_norm_criteria = 1.0D-20
    real(kind=DP) :: m_norm_criteria = 1.0D-18

    Do ir=1, nrc, dnr
       c1 = rho_rad(ir,1,1);
       c2 = sqrt( rho_rad(ir,2,1)**2 + rho_rad(ir,3,1)**2 + rho_rad(ir,4,1)**2 )

       d1 = ( c1 + c2 ) /2.0d0;
       d2 = ( c1 - c2 ) /2.0d0;

       rho_rad_along_quantz_axis(ir,1,1) = d1
       rho_rad_along_quantz_axis(ir,2,1) = d2

       Do isp=2, max_sph_expansion
          c1 = rho_rad(ir,1,isp);

          if ( c2 > m_norm_criteria ) then
             mmx = rho_rad(ir,2,1) / c2
             mmy = rho_rad(ir,3,1) / c2
             mmz = rho_rad(ir,4,1) / c2
          else
             mmx = 0.0;  mmy = 0.0;  mmz = 0.0
          endif

          c3 = rho_rad(ir,2,isp) *mmx +rho_rad(ir,3,isp) *mmy +rho_rad(ir,4,isp) *mmz

          d1 = ( c1 + c3 ) /2.0d0;
          d2 = ( c1 - c3 ) /2.0d0;

          rho_rad_along_quantz_axis(ir,1,isp) = d1
          rho_rad_along_quantz_axis(ir,2,isp) = d2
       End do
    End do
!
  end subroutine radial_chg_alng_LocalQuanzAxis2

  subroutine get_magmom_spherical( ia, nrc, dnr, rho_rad, wos, magmom_avg_spherical )
    integer, intent(in) :: nrc, ia, dnr
    real(kind=DP), intent(in) :: rho_rad(nrc,ndim_magmom,max_sph_expansion)
    real(kind=DP), intent(in) :: wos(nrc)
    real(kind=DP), intent(out) :: magmom_avg_spherical(3)

    integer :: is, ir
    real(kind=DP) :: csum

    Do is=2, ndim_magmom
       csum = 0.0d0
       Do ir=1, nrc, dnr
          csum = csum + rho_rad(ir,is,1) *wos(ir)
       End do
       magmom_avg_spherical( is-1 ) = csum
    End do
  end subroutine get_magmom_spherical

  subroutine get_magmom_spherical2( ia, nrc, dnr, rho_rad, wos, magmom_avg_spherical )
    integer, intent(in) :: nrc, ia, dnr
    real(kind=DP), intent(in) :: rho_rad(nrc,ndim_magmom,max_sph_expansion)
    real(kind=DP), intent(in) :: wos(nrc)
    real(kind=DP), intent(out) :: magmom_avg_spherical(3,max_sph_expansion)

    integer :: is, ir, ii
    real(kind=DP) :: csum

    magmom_avg_spherical = 0.0d0
    Do ii=1, max_sph_expansion
       Do is=2, ndim_magmom
          csum = 0.0d0
          Do ir=1, nrc, dnr
             csum = csum + rho_rad(ir,is,ii) *wos(ir)
          End do
          magmom_avg_spherical( is-1, ii ) = csum
       End do
    End do
  end subroutine get_magmom_spherical2
! ============================================================= 11.0
  subroutine m_PAWCD_set_ps_cd_sphex2_2D( ia,nspin,nrc,dnr &
       &                               ,msph,nps,msphmx &
       &                               ,num_isph,isph_chg )
    integer,intent(in) :: ia,nspin,nrc,dnr,msph
    real(8),intent(out):: nps(nrc,nspin,25)
    integer,intent(out):: msphmx
    integer,intent(out):: num_isph,isph_chg(25)

    real(8):: sum,fac,isqrt4pi,eps,abshsr, weight
    integer:: ir,it,is,lmt1,lmt2
    integer:: il1,il2,it1,it2,im1,im2,l3
    integer:: ii,jj,n,ilm3,iiqitg
    integer,pointer,dimension(:):: il3
    logical:: flg_isph(25)
!        real(8),pointer,dimension(:):: qsum

    call m_PP_find_maximum_l(n)
    n=(n-1)+(n-1)+1
    allocate(il3(25));call substitute_il3(25,il3)
!        allocate(qsum(nrc))
    it=ityp(ia)
    nps=0.d0
    msphmx=0
    eps=1.d-10
    flg_isph=.false.
    do lmt1=1,ilmt(it)
       il1=ltp(lmt1,it)
       im1=mtp(lmt1,it)
       it1=taup(lmt1,it)
       ii=(il1-1)**2+im1
       do lmt2=lmt1,ilmt(it)
          abshsr=abs(hsr(ia,lmt1,lmt2,1))
          if(nspin.eq.2) abshsr=abshsr+abs(hsr(ia,lmt1,lmt2,2))
          if(abshsr .lt. eps) cycle
          il2=ltp(lmt2,it)
          im2=mtp(lmt2,it)
          it2=taup(lmt2,it)
          jj=(il2-1)**2+im2
          fac=2.d0;if(lmt1.eq.lmt2) fac=1.d0
!                qsum=0.d0
          do n=1,il2p(lmt1,lmt2,it)
             ilm3=isph(lmt1,lmt2,n,it);l3=il3(ilm3)
             iiqitg=iqitg(il1,it1,il2,it2,l3+1,it)
!!$                    if(iiqitg == 0) cycle
             if(ilm3.gt.msphmx) msphmx=ilm3
             flg_isph(ilm3)=.true.
             do is=1,nspin
                sum=fac*hsr(ia,lmt1,lmt2,is)*dl2p(lmt1,lmt2,n,it)
                do ir=1,nrc,dnr
                   nps(ir,is,ilm3)=nps(ir,is,ilm3) + sum* &
                        & (phirpw(ir,il1,it1,it)*phirpw(ir,il2,it2,it) &
                        &  +qrspspw(ir,iiqitg))
                end do
             end do
          end do
       end do
    end do

    num_isph=0
    isph_chg=0
    do ilm3=1,25
       if(flg_isph(ilm3)) then
          num_isph=num_isph+1
          isph_chg(num_isph)=ilm3
       end if
    end do

    isqrt4pi=1.d0/sqrt(PAI4)
    do is=1,nspin
       do ir=1,nrc,dnr
          nps(ir,is,1)=(nps(ir,is,1)*isqrt4pi+rhpcrpw(ir,it)/PAI4/dble(nspin)) &
               /radr_paw(ir,it)/radr_paw(ir,it)
!            if(nps(ir,is,1).lt.0.d0) nps(ir,is,1)=0.d0
       end do
    end do

    if ( sw_opencore == ON .and. nspin == 2 .and. sw_xc_opencore_ae_only ==OFF ) then
       if ( has_opencore(it) == 1 ) then
          weight = 0.5d0 *mag_opencore_pol(ia,1)
          do ir=1,nrc,dnr
             nps(ir,1,1) = nps(ir,1,1) +rmag_opencore(ir,it)/PAI4 *weight &
                  &                        /radr_paw(ir,it)/radr_paw(ir,it)
             nps(ir,2,1) = nps(ir,2,1) -rmag_opencore(ir,it)/PAI4 *weight &
                  &                        /radr_paw(ir,it)/radr_paw(ir,it)
          end do
       endif
    endif

!        do n=2,msph
    do ilm3=2,num_isph
       n=isph_chg(ilm3)
       do is=1,nspin
          do ir=1,nrc,dnr
             nps(ir,is,n)=nps(ir,is,n)/radr_paw(ir,it)/radr_paw(ir,it)
!            if(nps(ir,is,n).lt.0.d0) nps(ir,is,n)=0.d0
          end do
       end do
    end do

!do ir=1,nrc
!print '(17e19.6)', radr_paw(ir,it),nps(ir,1,10:25)
!end do
!stop

    deallocate(il3)
    return
  end subroutine m_PAWCD_set_ps_cd_sphex2_2D

  subroutine m_PAWCD_set_ps_cd_sphex2( ia,nspin,nrc,dnr &
       &                               ,msph,nps,msphmx &
       &                               ,num_isph,isph_chg,ista_nrc,iend_nrc,ist )
    integer,intent(in) :: ia,nspin,nrc,dnr,msph,ista_nrc,iend_nrc,ist
    real(DP),intent(out):: nps(ista_nrc:iend_nrc,nspin,25) ! nps(nrc,nspin,25)
    integer,intent(out):: msphmx
    integer,intent(out):: num_isph,isph_chg(25)

    real(DP):: sum,fac,isqrt4pi,eps,abshsr, weight
    integer:: ir,it,is,lmt1,lmt2, ien
    integer:: il1,il2,it1,it2,im1,im2,l3
    integer:: ii,jj,n,ilm3,iiqitg
    integer,pointer,dimension(:):: il3
    logical:: flg_isph(25)

    ien = min(iend_nrc,nrc)

    call m_PP_find_maximum_l(n)
    n=(n-1)+(n-1)+1
    allocate(il3(25));call substitute_il3(25,il3)
    it=ityp(ia)
    nps=0.d0
    msphmx=0
    eps=1.d-10
    flg_isph=.false.
    do lmt1=1,ilmt(it)
       il1=ltp(lmt1,it)
       im1=mtp(lmt1,it)
       it1=taup(lmt1,it)
       ii=(il1-1)**2+im1
       do lmt2=lmt1,ilmt(it)
          abshsr=abs(hsr(ia,lmt1,lmt2,1))
          if(nspin.eq.2) abshsr=abshsr+abs(hsr(ia,lmt1,lmt2,2))
          if(abshsr .lt. eps) cycle
          il2=ltp(lmt2,it)
          im2=mtp(lmt2,it)
          it2=taup(lmt2,it)
          jj=(il2-1)**2+im2
          fac=2.d0;if(lmt1.eq.lmt2) fac=1.d0
          do n=1,il2p(lmt1,lmt2,it)
             ilm3=isph(lmt1,lmt2,n,it);l3=il3(ilm3)
             iiqitg=iqitg(il1,it1,il2,it2,l3+1,it)
!!$                    if(iiqitg == 0) cycle
             if(ilm3.gt.msphmx) msphmx=ilm3
             flg_isph(ilm3)=.true.
             do is=1,nspin
                sum=fac*hsr(ia,lmt1,lmt2,is)*dl2p(lmt1,lmt2,n,it)
                do ir=ist,ien,dnr ! do ir=1,nrc,dnr
                   nps(ir,is,ilm3)=nps(ir,is,ilm3) + sum* &
                        & (phirpw(ir,il1,it1,it)*phirpw(ir,il2,it2,it) &
                        &  +qrspspw(ir,iiqitg))
                end do
             end do
          end do
       end do
    end do

    num_isph=0
    isph_chg=0
    do ilm3=1,25
       if(flg_isph(ilm3)) then
          num_isph=num_isph+1
          isph_chg(num_isph)=ilm3
       end if
    end do

    isqrt4pi=1.d0/sqrt(PAI4)
    do is=1,nspin
       do ir=ist,ien,dnr !     do ir=1,nrc,dnr
          nps(ir,is,1)=(nps(ir,is,1)*isqrt4pi+rhpcrpw(ir,it)/PAI4/dble(nspin)) &
               /radr_paw(ir,it)/radr_paw(ir,it)
!            if(nps(ir,is,1).lt.0.d0) nps(ir,is,1)=0.d0
       end do
    end do

    if ( sw_opencore == ON .and. nspin == 2 .and. sw_xc_opencore_ae_only ==OFF ) then
       if ( has_opencore(it) == 1 ) then
          weight = 0.5d0 *mag_opencore_pol(ia,1)
          do ir=ist, ien, dnr
             nps(ir,1,1) = nps(ir,1,1) +rmag_opencore(ir,it)/PAI4 *weight &
                  &                        /radr_paw(ir,it)/radr_paw(ir,it)
             nps(ir,2,1) = nps(ir,2,1) -rmag_opencore(ir,it)/PAI4 *weight &
                  &                        /radr_paw(ir,it)/radr_paw(ir,it)
          end do
       endif
    endif


!        do n=2,msph
    do ilm3=2,num_isph
       n=isph_chg(ilm3)
       do is=1,nspin
          do ir=ist,ien,dnr ! do ir=1,nrc,dnr
             nps(ir,is,n)=nps(ir,is,n)/radr_paw(ir,it)/radr_paw(ir,it)
!            if(nps(ir,is,n).lt.0.d0) nps(ir,is,n)=0.d0
          end do
       end do
    end do

    deallocate(il3)
    return
  end subroutine m_PAWCD_set_ps_cd_sphex2


! ===================== added by K. Tagami ==================== 11.0
  subroutine m_PAWCD_ps_cd_sphex2_nonclA( ia, nspin, nrc, dnr, &
       &                                  msph, nps, msphmx, &
       &                                  num_isph, isph_chg, wos, &
       &                                  level_of_projection, &
       &                                  magmom_avg_spherical, rho_rad )

    integer, intent(in) :: ia,nspin,nrc,dnr,msph
    integer, intent(in) :: level_of_projection
    integer, intent(out):: msphmx
    integer, intent(out):: num_isph,isph_chg(max_sph_expansion)

    real(kind=DP), intent(out):: nps( nrc, nspin, max_sph_expansion )
    real(kind=DP), intent(in) :: wos(mmesh)
    real(kind=DP), intent(out) :: magmom_avg_spherical(3,max_sph_expansion)
    real(kind=DP), intent(out) :: rho_rad( nrc, ndim_magmom, max_sph_expansion )

    real(kind=DP) :: sum,fac,isqrt4pi,eps,abshsr, ctmp, weight

    integer:: ir,it,is,lmt1,lmt2
    integer:: il1,il2,it1,it2,im1,im2,l3
    integer:: ii,jj,n,ilm3,iiqitg
    integer,pointer,dimension(:):: il3

    logical:: flg_isph( max_sph_expansion )

    call m_PP_find_maximum_l(n)
    n=(n-1)+(n-1)+1
    allocate(il3(max_sph_expansion))
    call substitute_il3(max_sph_expansion,il3)

    it=ityp(ia)
    nps=0.d0; rho_rad = 0.0d0

    msphmx=0
    eps=1.d-10
    flg_isph=.false.

    do lmt1=1,ilmt(it)
       il1=ltp(lmt1,it); im1=mtp(lmt1,it); it1=taup(lmt1,it)
       ii=(il1-1)**2+im1

       do lmt2=lmt1,ilmt(it)
          abshsr =abs(hsr(ia,lmt1,lmt2,1))
          if(abshsr .lt. eps) cycle

          il2=ltp(lmt2,it);  im2=mtp(lmt2,it); it2=taup(lmt2,it)
          jj=(il2-1)**2+im2

          fac=2.d0;if(lmt1.eq.lmt2) fac=1.d0

          do n=1,il2p(lmt1,lmt2,it)
             ilm3=isph(lmt1,lmt2,n,it);l3=il3(ilm3)
             iiqitg=iqitg(il1,it1,il2,it2,l3+1,it)

             if(ilm3.gt.msphmx) msphmx=ilm3
             flg_isph(ilm3)=.true.

             do is=1, ndim_magmom
                sum=fac*hsr(ia,lmt1,lmt2,is)*dl2p(lmt1,lmt2,n,it)
                do ir=1,nrc,dnr
                   rho_rad(ir,is,ilm3) = rho_rad(ir,is,ilm3) &
                        &               + sum* ( phirpw(ir,il1,it1,it) &
                        &                       *phirpw(ir,il2,it2,it) &
                        &                        +qrspspw(ir,iiqitg) )
                end do
             end do
          end do
       end do
    end do

    num_isph=0; isph_chg=0

    do ilm3=1,max_sph_expansion
       if(flg_isph(ilm3)) then
          num_isph=num_isph+1
          isph_chg(num_isph)=ilm3
       end if
    end do

    if ( mode_of_projection == project_on_local_quantz_axis ) then

       select case( level_of_projection )
       case (1)
          call get_magmom_spherical( ia, nrc, dnr, rho_rad, wos, &
               &                     magmom_avg_spherical(:,1) )
          call radial_chg_along_LocalQuanzAxis( nrc, dnr, rho_rad, &
               &                                magmom_avg_spherical(:,1), nps )

       case(2,3)
          call get_magmom_spherical2( ia, nrc, dnr, rho_rad, wos, magmom_avg_spherical )
          call radial_chg_along_LocalQuanzAxis( nrc, dnr, rho_rad, &
               &                                magmom_avg_spherical(:,1), nps )

       case(4)
          call get_magmom_spherical( ia, nrc, dnr, rho_rad, wos, &
               &                     magmom_avg_spherical(:,1) )
          call radial_chg_alng_LocalQuanzAxis2( nrc, dnr, rho_rad, nps )

       end select

    else
       call phase_error_with_msg(nfout,'kt : under conctruction ',__LINE__,__FILE__)
    end if

    isqrt4pi=1.d0/sqrt(PAI4)
    do is=1,nspin
       do ir=1,nrc,dnr
          nps(ir,is,1) = (nps(ir,is,1)*isqrt4pi+rhpcrpw(ir,it)/PAI4/dble(nspin)) &
               /radr_paw(ir,it)/radr_paw(ir,it)
       end do
    end do

    if ( sw_opencore == ON .and. nspin == 2 .and. sw_xc_opencore_ae_only ==OFF ) then
       if ( has_opencore(it) == 1 ) then
          ctmp = magmom_avg_spherical(1,1)**2 +magmom_avg_spherical(2,1)**2 &
               & +magmom_avg_spherical(3,1)**2
          ctmp = sqrt(ctmp)
          if ( ctmp > 1.0D-11 ) then
             weight = magmom_avg_spherical(1,1)/ctmp *mag_opencore_pol(ia,1) &
                  &  +magmom_avg_spherical(2,1)/ctmp *mag_opencore_pol(ia,2) &
                  &  +magmom_avg_spherical(3,1)/ctmp *mag_opencore_pol(ia,3)
             weight = 0.5d0 *weight
          else
             weight = 0.50d0
          endif

          do ir=1,nrc,dnr
             nps(ir,1,1) = nps(ir,1,1) +rmag_opencore(ir,it)/PAI4 *weight &
                  &                        /radr_paw(ir,it)/radr_paw(ir,it)
             nps(ir,2,1) = nps(ir,2,1) -rmag_opencore(ir,it)/PAI4 *weight &
                  &                        /radr_paw(ir,it)/radr_paw(ir,it)
          end do
       endif
    endif

    do ilm3=2,num_isph
       n=isph_chg(ilm3)
       do is=1,nspin
          do ir=1,nrc,dnr
             nps(ir,is,n)=nps(ir,is,n)/radr_paw(ir,it)/radr_paw(ir,it)
          end do
       end do
    end do

    deallocate(il3)
    return
  end subroutine m_PAWCD_ps_cd_sphex2_nonclA
! ================================================================= 11.0

  subroutine m_PAWCD_set_cr2_isph2_mmt2()
    integer:: lcmax
    logical,save:: initialized=.false.
    if(.not.initialized) then
       lcmax=4
       allocate(paw_cr2(25,25,6))
       allocate(paw_isph2(25,25,6))
       allocate(paw_mmt2(25,25))
       call sphset3(nfout,ipripaw,paw_cr2,paw_isph2,paw_mmt2)
       initialized=.true.
    end if
    return
  end subroutine m_PAWCD_set_cr2_isph2_mmt2

  subroutine m_PAWCD_set_sq_der_cd_sdphex2( ia,nspin,nrc,dnr &    ! == For nrc decomposion. by takto 2012/12/05 ==
       &                    , num_isph_chg, isph_chg, nae_sph &
       &                    , dnae_dr_sph, grad_nae2, grad_tnae2, msphmx_grd &
       &                    , num_isph_grd ,isph_grd,ista_nrc,iend_nrc,ist )

    integer,intent(in) :: ia,nspin,nrc,dnr,num_isph_chg,isph_chg(25), ista_nrc,iend_nrc,ist

    real(DP),intent(in) :: nae_sph(ista_nrc:iend_nrc,nspin,25)     !  nae_sph(nrc,nspin,25)
    real(DP),intent(in) :: dnae_dr_sph(ista_nrc:iend_nrc,nspin,25) !  dnae_dr_sph(nrc,nspin,25)
    real(DP),intent(out):: grad_nae2(ista_nrc:iend_nrc,nspin,25)   !  grad_nae2(nrc,nspin,25)
    real(DP),intent(out):: grad_tnae2(ista_nrc:iend_nrc,25)        !  grad_tnae2(nrc,25)

    integer,intent(out):: msphmx_grd,num_isph_grd
    integer,intent(out):: isph_grd(25)

    real(DP):: fac
    integer:: ir,it,is,ien

    real(DP):: sum1(ista_nrc:iend_nrc,nspin),sum2(ista_nrc:iend_nrc,nspin,25) & ! sum1(nrc,nspin),sum2(nrc,nspin)
         &   , sum3(ista_nrc:iend_nrc),sum4(ista_nrc:iend_nrc,25)   ! sum3(nrc),sum4(nrc,25)

    real(DP):: dl,dl1,dl2,dl3,cijk
    integer:: isp,isp2,isp3,lcmax,n
    integer:: nsp,nsp2,nsp3
    logical:: flg_isp(25)
    integer,pointer,dimension(:):: il3
    logical,save:: initialized=.false.

    ien = min(iend_nrc, nrc)

    if(.not.initialized) then
       lcmax=4
!            allocate(paw_cr2(16,16,6))
!            allocate(paw_isph2(16,16,6))
!            allocate(paw_mmt2(16,16))
!            call sphset2(nfout,ipri,lcmax,paw_cr2,paw_isph2,paw_mmt2)
       allocate(paw_cr2(25,25,6))
       allocate(paw_isph2(25,25,6))
       allocate(paw_mmt2(25,25))
       call sphset3(nfout,ipripaw,paw_cr2,paw_isph2,paw_mmt2)
       initialized=.true.
    end if

    allocate(il3(25));call substitute_il3(25,il3)
    it=ityp(ia)
    grad_nae2=0.d0
    grad_tnae2=0.d0
    sum1=0.d0
    sum3=0.d0
    flg_isp=.false.

    do nsp=2,num_isph_chg !      do isp=2,msphmx
       isp=isph_chg(nsp)
       dl=dble(il3(isp))
       dl=dl*(dl+1.d0)
       do is=1,nspin

          do ir=ist,ien,dnr !      do ir=1,nrc,dnr
             sum1(ir,is)=sum1(ir,is) + &
                  dnae_dr_sph(ir,is,isp)**2 + &
                  dl*nae_sph(ir,is,isp)**2/radr_paw(ir,it)/radr_paw(ir,it)
          end do
       end do
       if(nspin.eq.2) then
          do ir=ist,ien,dnr !       do ir=1,nrc,dnr
             sum3(ir)=sum3(ir) + &
                  (dnae_dr_sph(ir,1,isp)+dnae_dr_sph(ir,2,isp))**2 + &
                  dl*(nae_sph(ir,1,isp)+nae_sph(ir,2,isp))**2/radr_paw(ir,it)/radr_paw(ir,it)
          end do
       end if
    end do
    do is=1,nspin
       do ir=ist,ien,dnr !   do ir=1,nrc,dnr
          grad_nae2(ir,is,1)=dnae_dr_sph(ir,is,1)**2 + sum1(ir,is)/PAI4
       end do
    end do
    if(nspin.eq.2) then
       do ir=ist,ien,dnr !   do ir=1,nrc,dnr
          grad_tnae2(ir,1)=(dnae_dr_sph(ir,1,1)+dnae_dr_sph(ir,2,1))**2+sum3(ir)/PAI4
       end do
    end if

    flg_isp(1)=.true.

    sum2=0.d0
    sum4=0.d0
    msphmx_grd=0
!        do isp2=2,min(16,msphmx)                                ! sphset2
!        do isp2=2,msphmx                                        ! sphset3
    do nsp2=2,num_isph_chg
       isp2=isph_chg(nsp2)
       dl2=dble(il3(isp2))
       dl2=dl2*(dl2+1.d0)
!            do isp3=isp2,min(16,msphmx)                         ! sphset2
!            do isp3=isp2,msphmx                                 ! sphset3
       do nsp3=nsp2,num_isph_chg
          isp3=isph_chg(nsp3)
          dl3=dble(il3(isp3))
          dl3=dl3*(dl3+1.d0)
          fac=2.d0;if(isp2.eq.isp3) fac=1.d0
          do n=1,paw_mmt2(isp2,isp3)
             isp=paw_isph2(isp2,isp3,n)
             if(isp.eq.1) cycle
             if(isp.gt.msphmx_grd) msphmx_grd=isp
             flg_isp(isp)=.true.
             cijk=paw_cr2(isp2,isp3,n)
             dl=dble(il3(isp))
             dl=dl*(dl+1.d0)
             do is=1,nspin

                do ir=ist,ien,dnr !               do ir=1,nrc,dnr
                   sum2(ir,is,isp)=sum2(ir,is,isp) &
                        & + fac*cijk*( &
                        &    dnae_dr_sph(ir,is,isp2)*dnae_dr_sph(ir,is,isp3) &
                        &  + 0.5d0*(dl2+dl3-dl)*nae_sph(ir,is,isp2)*nae_sph(ir,is,isp3) &
                        &     /radr_paw(ir,it)/radr_paw(ir,it) )
                end do
             end do

             if(nspin.eq.2) then
                do ir=ist,ien,dnr    !            do ir=1,nrc,dnr
                   sum4(ir,isp)=sum4(ir,isp) &
                        & + fac*cijk*( &
                        &    (dnae_dr_sph(ir,1,isp2)+dnae_dr_sph(ir,2,isp2)) &
                        &    * (dnae_dr_sph(ir,1,isp3)+dnae_dr_sph(ir,2,isp3)) &
                        &   +  0.5d0*(dl2+dl3-dl) &
                        &     * (nae_sph(ir,1,isp2)+nae_sph(ir,2,isp2)) &
                        &     * (nae_sph(ir,1,isp3)+nae_sph(ir,2,isp3)) &
                        &      /radr_paw(ir,it)/radr_paw(ir,it) )
                end do
             end if

          end do
       end do
    end do

    num_isph_grd=0
    isph_grd=0
    do isp=1,25
       if(flg_isp(isp)) then
          num_isph_grd=num_isph_grd+1
          isph_grd(num_isph_grd)=isp
       end if
    end do

    do nsp=2,num_isph_grd   !     do isp=2,msphmx_grd
       isp=isph_grd(nsp)
       do is=1,nspin
          do ir=ist,ien,dnr !      do ir=1,nrc,dnr
             grad_nae2(ir,is,isp)=2.d0*dnae_dr_sph(ir,is,1)*dnae_dr_sph(ir,is,isp) &
                  & + sum2(ir,is,isp)
          end do
       end do

       if(nspin.eq.2) then
          do ir=ist,ien,dnr !      do ir=1,nrc,dnr
             grad_tnae2(ir,isp)=2.d0*(dnae_dr_sph(ir,1,1)+dnae_dr_sph(ir,2,1)) * &
                  &               (dnae_dr_sph(ir,1,isp)+dnae_dr_sph(ir,2,isp))  &
                  &                + sum4(ir,isp)
          end do
       end if

    end do

    if(nspin.eq.1) then
!            do isp=1,msphmx_grd
       do nsp=1,num_isph_grd
          isp=isph_grd(nsp)
          do ir=ist,ien,dnr   !     do ir=1,nrc,dnr
             grad_tnae2(ir,isp)=grad_nae2(ir,1,isp)
          end do
       end do
    end if

    deallocate(il3)

    return
  end subroutine m_PAWCD_set_sq_der_cd_sdphex2
! ==============================================================================


end module m_PAW_ChargeDensity
