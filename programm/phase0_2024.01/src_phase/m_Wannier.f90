!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  MODULE: m_Wannier
!
!  AUTHOR(S): T. Yamamoto   May/05/2007
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
module m_Wannier
!! Maximally localized Wannier functions
!! Ref.
!! [1] N. Marzari and D. Vanderbilt, Phys. Rev. B 56 (1997) 12847.
!! [2] G. Berghold, et.al., Phys. Rev. B 61 (2000) 10040.
  use m_Const_Parameters,   only : DP,PAI2,SD,CG,CGPRC,BINARY,CUBE,VTK,DENSITY_ONLY &
       &                         , SKIP,EXECUT,GAMMA,ON
  use m_Control_Parameters, only : printable, kimg, neg, eps_wan, dt_wan &
       &                         , wannier_opt_method, max_iter_wan &
       &                         , wannier_filetype, ipri, sw_random_wannier &
       &                         , sw_potential_wannier, sw_continue_wannier
  use m_Crystal_Structure,  only : altv, rltv, univol
  use m_FFT,                only : nfft, fft_box_size_WF &
       &                         , m_FFT_alloc_WF_work &
       &                         , m_FFT_dealloc_WF_work
  use m_Electronic_Structure,only : totch, m_ES_WF_in_Rspace, neordr &
       &                         , fsr_l, fsi_l
  use m_Timing,             only : tstatc0_begin, tstatc0_end
  use m_Ionic_System,       only : ntyp,natm,natm2,ityp,iatomn,m_IS_pack_all_ions_in_uc &
       &                         , iwei,pos
  use m_PseudoPotential,    only : ival,modnrm,dk_wan,qitg_wan &
       &                         , m_PP_include_vanderbilt_pot &
       &                         , n_non0_lmtxlmt,index_lmt1_lmt2,lmta,nlmt,ilmt &
       &                         , ltp,taup,il2p,isph,iqitg,dl2p &
       &                         , m_PP_find_maximum_l
  use m_Files,              only : m_Files_open_nfwannier,nfwannier &
       &                         , m_Files_open_nfcntn_wannier &
       &                         , m_Files_open_nfpot_wannier &
       &                         , nfcntn_wannier,nfpot_wannier
  use m_Kpoints,            only : k_symmetry
  implicit none
  
  integer, public                     :: nwght
  real(kind=DP), private, allocatable :: wght(:) ! d(nwght)
  integer, public, allocatable        :: ghat(:,:) ! d(3,nwght) 

  integer, private                       :: nocc
  complex(kind=DP), private, allocatable :: zmat(:,:,:) ! d(nocc,nocc,nwght)
  complex(kind=DP), private, allocatable :: zn(:,:) ! d(nocc,nwght)

  real(kind=DP), private, allocatable :: umat(:,:) ! d(nocc,nocc)
  complex(kind=DP), private, allocatable :: bmat(:,:) ! d(nocc,nocc)
  real(kind=DP), private, allocatable :: mmat(:,:) ! d(nocc,nocc)

  integer, private                    :: ntri ! (nocc-1)*nocc/2
  real(kind=DP), private, allocatable :: amat(:) ! d(ntri)
  real(kind=DP), private, allocatable :: grad(:) ! d(ntri)
  real(kind=DP), private :: omega
  real(kind=DP), private :: gradmax

  complex(kind=DP), private, allocatable :: rmat(:,:) ! d(nocc,nocc)
  real(kind=DP), private, allocatable    :: lam(:) ! d(nocc)

  integer, private                    :: iteration

  complex(kind=DP), private, allocatable :: ftqe(:,:,:,:)!d(nlmt,nlmt,natm,nwght)

contains
  subroutine m_Wan_gen_weight(nfout)
    integer, intent(in) :: nfout
    integer :: i
    real(kind=DP) :: g(3,3),w
    real(kind=DP), parameter :: eps = 1.d-10

    nwght = 6
    allocate(wght(nwght),ghat(3,nwght))

    g = matmul(transpose(altv),altv)

    nwght = 0
    ! g^=(1,0,0)
    w = g(1,1)-g(1,2)-g(1,3)
    nwght = nwght + 1
    wght(nwght) = w
    ghat(1:3,nwght) = (/1,0,0/)
    ! g^=(0,1,0)
    w = g(2,2)-g(1,2)-g(2,3)
    nwght = nwght + 1
    wght(nwght) = w
    ghat(1:3,nwght) = (/0,1,0/)
    ! g^=(0,0,1)
    w = g(3,3)-g(1,3)-g(2,3)
    nwght = nwght + 1
    wght(nwght) = w
    ghat(1:3,nwght) = (/0,0,1/)
    ! g^=(1,1,0)
    w = g(1,2)
    if(abs(w) > eps) then
       nwght = nwght + 1
       wght(nwght) = w
       ghat(1:3,nwght) = (/1,1,0/)
    end if
    ! g^=(1,0,1)
    w = g(1,3)
    if(abs(w) > eps) then
       nwght = nwght + 1
       wght(nwght) = w
       ghat(1:3,nwght) = (/1,0,1/)
    end if
    ! g^=(1,1,0)
    w = g(2,3)
    if(abs(w) > eps) then
       nwght = nwght + 1
       wght(nwght) = w
       ghat(1:3,nwght) = (/0,1,1/)
    end if

    if(printable) then
       write(nfout,*) '!Wannier: nwght=',nwght
       do i=1,nwght
          write(nfout,'(" !Wannier: i=",i1," wght=",f18.5," g^=",3i2)') &
          & i,wght(i),ghat(1:3,i)
       end do
    end if

  end subroutine m_Wan_gen_weight

  subroutine m_Wan_opt_Omega(nfout)
    integer, intent(in) :: nfout
    integer :: id_sname = -1
    integer :: i
    call tstatc0_begin('m_Wan_opt_Omega ',id_sname)

    call alloc_zmat(nfout)
    if(sw_continue_wannier == ON) then
       call m_Files_open_nfcntn_wannier
       call rd_zmat_and_umat(nfout,nfcntn_wannier) ! => Z_Iij, U
    else
       call cal_zmat(nfout) ! => Z_Iij
    end if
    if(sw_random_wannier==ON) call random_number(amat)
    iteration = 0
    do 
       iteration = iteration + 1
       call diag_amat() ! => R+,L
       call cal_umat_and_bmat() ! => U,B
       amat(1:ntri) = 0.d0
       call update_zn() ! => Zn
       call cal_mmat() ! => M
       call cal_grad() ! => g = dOmega/dA, gradmax
       call cal_Omega() ! => Omega
       if(Omega_is_max(nfout)) exit
       call update_amat() ! => A
    end do
    call wannier_center(nfout)
    call wannier_function(nfout)
    call m_Files_open_nfcntn_wannier
    call wd_zmat_and_umat(nfout,nfcntn_wannier)! <= Z_Iij, U
    call dealloc_zmat()

    call tstatc0_end(id_sname)
  end subroutine m_Wan_opt_Omega

  subroutine alloc_zmat(nfout)
    integer, intent(in) :: nfout
    integer :: i
    nocc = int(totch+1.d-4)/2
    ntri = nocc*(nocc-1)/2
    if(printable) then
       write(nfout,*) '!Wannier: nocc=',nocc
    end if
    allocate(zmat(nocc,nocc,nwght))
    zmat(1:nocc,1:nocc,1:nwght) = (0.d0,0.d0)
    allocate(zn(nocc,nwght))
    zn(1:nocc,1:nwght) = (0.d0,0.d0)
    allocate(umat(nocc,nocc))
    umat(1:nocc,1:nocc) = 0.d0
    do i=1,nocc
       umat(i,i) = 1.d0
    end do
    allocate(bmat(nocc,nocc))
    allocate(mmat(nocc,nocc))
    allocate(amat(ntri))
    amat(1:ntri) = 0.d0
    allocate(grad(ntri))
    allocate(rmat(nocc,nocc))
    allocate(lam(nocc))
  end subroutine alloc_zmat

  subroutine dealloc_zmat
    deallocate(zmat)
    deallocate(zn)
    deallocate(umat)
    deallocate(bmat)
    deallocate(mmat)
    deallocate(amat)
    deallocate(grad)
  end subroutine dealloc_zmat

  subroutine cal_zmat(nfout)
    integer, intent(in) :: nfout
    integer :: id_sname = -1
    integer :: i,ip,iw,ik
    integer :: ib1,ib2,jb1,jb2
    integer :: nffth,ngrid
    integer :: id1,id2,id12
    integer :: r1,r2,r3
    real(kind=DP), allocatable :: rhat(:,:) !d(nfft/2,3)
    real(kind=DP), allocatable :: wf1(:) !d(nfft)
    real(kind=DP), allocatable :: wf2(:) !d(nfft)
    real(kind=DP), allocatable :: chgr(:) !d(nfft/2)
    real(kind=DP), allocatable :: chgi(:) !d(nfft/2)
    real(kind=DP), allocatable :: zcos(:,:) !d(nfft/2)
    real(kind=DP), allocatable :: zsin(:,:) !d(nfft/2)
    real(kind=DP) :: zr,zi,ph
    real(kind=DP) :: da(3)
    integer :: u,v,lmt1,lmt2,ia,it,mdvdb
    complex(kind=DP) :: defchg,fac,fs1u,fs1v,fs2u,fs2v
    integer :: kimg_t

    call tstatc0_begin('cal_zmat in m_Wan ',id_sname)

    ik = 1
    if(k_symmetry(ik) == GAMMA) then
       kimg_t = 1
    else
       kimg_t = kimg
    end if

    if(modnrm == EXECUT) call constract_of_ftq()

    ngrid = product(fft_box_size_WF(1:3,1))

    nffth = nfft/2
    allocate(wf1(nfft))
    allocate(wf2(nfft))
    allocate(chgr(nffth))
    if(k_symmetry(ik) /= GAMMA) allocate(chgi(nffth))
    call m_FFT_alloc_WF_work()

    allocate(zcos(nffth,nwght))
    allocate(zsin(nffth,nwght))

    ! Phase factors, exp(iGI*r)
    allocate(rhat(nffth,3))
    id1 = fft_box_size_WF(1,0)
    id2 = fft_box_size_WF(2,0)
    id12 = id1*id2
    da(1:3) = 1.d0/fft_box_size_WF(1:3,1)
    do i=1,nffth
       ip = i-1
       r3 = ip/id12
       r2 = (ip-r3*id12)/id1
       r1 = ip-r2*id1-r3*id12
       rhat(i,1) = dble(r1)*da(1)
       rhat(i,2) = dble(r2)*da(2)
       rhat(i,3) = dble(r3)*da(3)
    end do
    ! debug
    !!  do i=1,nffth
    !!     write(nfout,*) 'i=',i,' rhat=',rhat(i,1:3)
    !!  end do
    ! end debug
    do iw=1,nwght
       do i=1,nffth
          ph = PAI2*dot_product(ghat(:,iw),rhat(i,:))
          zcos(i,iw) = cos(ph)
          zsin(i,iw) = sin(ph)
       end do
    end do
    deallocate(rhat)

    ! ZI,ij = <i|exp(iGI*r)|j>
    ik=1
    do ib2 = 1,nocc
       jb2 = neordr(ib2,ik)
       call m_ES_WF_in_Rspace(ik,jb2,wf2)
       do ib1 = 1,ib2
          jb1 = neordr(ib1,ik)
          call m_ES_WF_in_Rspace(ik,jb1,wf1)
          if(k_symmetry(ik) == GAMMA) then
             do i=1,nffth
                ip = 2*i-1
                chgr(i) = wf1(ip)*wf2(ip)
             end do
          else
             do i=1,nffth
                ip = 2*i-1
                chgr(i) = wf1(ip)*wf2(ip)+wf1(ip+1)*wf2(ip+1)
                chgi(i) = wf1(ip)*wf2(ip+1)-wf1(ip+1)*wf2(ip)
             end do
          end if
          ! debug
          write(nfout,*) ib1,ib2,sum(chgr)/ngrid
          ! end debug
          do iw=1,nwght
             zr=0.d0; zi=0.d0
             if(k_symmetry(ik) == GAMMA) then
                do i=1,nffth
                  zr = zr + zcos(i,iw)*chgr(i)
                  zi = zi + zsin(i,iw)*chgr(i)
                end do
             else
                do i=1,nffth
                  zr = zr + zcos(i,iw)*chgr(i) - zsin(i,iw)*chgi(i)
                  zi = zi + zsin(i,iw)*chgr(i) + zcos(i,iw)*chgi(i)
                end do
             end if
             zmat(ib1,ib2,iw) = dcmplx(zr,zi)/ngrid
             if(modnrm == EXECUT) then
                defchg = (0.d0,0.d0)
                do ia=1,natm
                   it=ityp(ia)
                   mdvdb = m_PP_include_vanderbilt_pot(it)
                   if(mdvdb == SKIP) cycle
                   ph = PAI2*dot_product(ghat(:,iw),pos(ia,:))
                   if(kimg_t==1) then
                      do ip = 1, n_non0_lmtxlmt(it)
                         lmt1  = index_lmt1_lmt2(ip,it,1)
                         lmt2  = index_lmt1_lmt2(ip,it,2)
                         u     = lmta(lmt1,ia)
                         v     = lmta(lmt2,ia)
                         fac   = ftqe(lmt1,lmt2,it,iw)*iwei(ia)*exp(dcmplx(0.d0,ph))
                         !!$write(6,*) 'ip,lmt1,lmt2,u,v,ftqe,fac,=', &
                         !!$  & ip,lmt1,lmt2,u,v,ftqe(lmt1,lmt2,it,iw),fac
                         if(lmt1/=lmt2) then
                            defchg = defchg + &
                         & fac*(fsr_l(jb1,u,ik)*fsr_l(jb2,v,ik) + &
                         &      fsr_l(jb1,v,ik)*fsr_l(jb2,u,ik) )
                         else
                            defchg = defchg + &
                         & fac*fsr_l(jb1,u,ik)*fsr_l(jb2,v,ik)
                         end if
                      end do
                   else
                      do ip = 1, n_non0_lmtxlmt(it)
                         lmt1  = index_lmt1_lmt2(ip,it,1)
                         lmt2  = index_lmt1_lmt2(ip,it,2)
                         u     = lmta(lmt1,ia)
                         v     = lmta(lmt2,ia)
                         fac   = ftqe(lmt1,lmt2,it,iw)*iwei(ia)*exp(dcmplx(0.d0,ph))
                         if(lmt1/=lmt2) then
                            fs1u  = dcmplx(fsr_l(jb1,u,ik),fsi_l(jb1,u,ik))
                            fs1v  = dcmplx(fsr_l(jb1,v,ik),fsi_l(jb1,v,ik))
                            fs2u  = dcmplx(fsr_l(jb2,u,ik),fsi_l(jb2,u,ik))
                            fs2v  = dcmplx(fsr_l(jb2,v,ik),fsi_l(jb2,v,ik))
                            defchg = defchg + fac*(fs1u*fs2v+fs1v*fs2v)
                         else
                            fs1u  = dcmplx(fsr_l(jb1,u,ik),fsi_l(jb1,u,ik))
                            fs2v  = dcmplx(fsr_l(jb2,v,ik),fsi_l(jb2,v,ik))
                            defchg = defchg + fac*fs1u*fs2v
                         end if
                      end do
                   end if
                end do
                zmat(ib1,ib2,iw) = zmat(ib1,ib2,iw) + defchg
             end if
             zmat(ib2,ib1,iw) = zmat(ib1,ib2,iw)
          end do
       end do
    end do

    deallocate(wf1,wf2,chgr)
    if(allocated(chgi)) deallocate(chgi)
    deallocate(zcos,zsin)
    call m_FFT_dealloc_WF_work()

    ! debug
      do iw=1,nwght
         write(nfout,*) 'iw=',iw
         write(nfout,*) 'Real'
         do ib1=1,min(nocc,5)
            write(nfout,'(5(1x,e12.4))') dble(zmat(ib1,1:min(nocc,5),iw))
         end do
         write(nfout,*) 'Imag.'
         do ib1=1,min(nocc,5)
            write(nfout,'(5(1x,e12.4))') dimag(zmat(ib1,1:min(nocc,5),iw))
         end do
      end do
    ! end debug

    call tstatc0_end(id_sname)
  end subroutine cal_zmat

  subroutine diag_amat
    integer :: i,j,ip
    character(len=1), parameter :: JOBZ = 'V', UPLO = 'U'
    integer :: lwork,lrwork,liwork
    complex(kind=DP),allocatable,dimension(:) :: work
    real(kind=DP),allocatable,dimension(:) :: rwork
    integer,allocatable,dimension(:) :: iwork
    integer :: info

    rmat = (0.d0,0.d0)
    do j=2,nocc
       do i=1,j-1
          ip = i + (j-2)*(j-1)/2
          rmat(i,j) = dcmplx(0.d0,-amat(ip))
       end do
    end do

    lwork=nocc*(2+nocc)
    lrwork=1+nocc*(5+2*nocc)
    liwork=3+5*nocc

    allocate(work(lwork))
    allocate(rwork(lrwork))
    allocate(iwork(liwork))

    call zheevd(JOBZ,UPLO,nocc,rmat,nocc,lam,work,lwork,rwork,lrwork,iwork,liwork,info)
    if(info > 0) then
       write(6,*) 'info=',info
       stop 'zheevd error'
    end if

    deallocate(work)
    deallocate(rwork)
    deallocate(iwork)
  end subroutine diag_amat

  subroutine cal_umat_and_bmat()
    ! U = exp(A) = R+ * e^L * R
    ! B_kl = (e^lam_k-e^lam_l)/(lam_k-lam_l)
    complex(kind=DP) :: explam(nocc)
    complex(kind=DP) :: tmat(nocc,nocc)
    complex(kind=DP) :: ctmp
    real(kind=DP) :: dl,dl2,dl3,cr,ci
    real(kind=DP), parameter :: r6 = 1.d0/6.d0
    real(kind=DP), parameter :: r24 = 1.d0/24.d0
    real(kind=DP), parameter :: r120 = 1.d0/120.d0
    real(kind=DP), parameter :: r720 = 1.d0/720.d0
    real(kind=DP) :: umat0(nocc,nocc)
    integer :: i,j

    explam = dcmplx(cos(lam),sin(lam))
    do j=1,nocc
       do i=1,nocc
          tmat(i,j) = rmat(i,j) * explam(j)
       end do
    end do
    umat0 = matmul(tmat,transpose(dconjg(rmat)))
    umat = matmul(umat,umat0)
    do i=1,nocc
       bmat(i,i) = explam(i)
    end do
    do j=1,nocc
       do i=1,j-1
          dl = lam(i)-lam(j)
          if(abs(dl) > 1.d-3) then
             ctmp = explam(i)-explam(j)
             bmat(i,j) = dcmplx(dimag(ctmp),-dble(ctmp))/dl
          else
             dl2 = dl*dl
             dl3 = dl*dl2
             cr=1.d0-r6*dl2+r120*dl2*dl2
             ci=0.5d0*dl-r24*dl3+r720*dl2*dl3
             bmat(i,j) = explam(j)*dcmplx(cr,ci)
          end if
       end do
       bmat(j,i) = bmat(i,j)
    end do
  end subroutine cal_umat_and_bmat

  subroutine update_zn
    integer :: iw,i,j,n
    ! update Z_In
    do iw=1,nwght
       do n=1,nocc
          zn(n,iw) = (0.d0,0.d0)
          do j=1,nocc
             do i=1,nocc
                zn(n,iw) = zn(n,iw) + umat(i,n)*umat(j,n)*zmat(i,j,iw)
             end do
          end do
       end do
    end do
  end subroutine update_zn

  subroutine cal_mmat()
    integer :: i,j,iw
    real(kind=DP) :: fac
    complex(kind=DP) :: zz
    mmat = 0.d0
    do iw=1,nwght
       fac = 4.d0*wght(iw)
       do j=1,nocc
          do i=1,nocc
             !!zz = dot_product(umat(:,i),zmat(:,j,iw))
             zz = dot_product(umat(:,i),matmul(zmat(:,:,iw),umat(:,j)))
             mmat(i,j) = mmat(i,j) + &
             & fac*(dble(zz)*dble(zn(j,iw))+dimag(zz)*dimag(zn(j,iw)))
          end do
       end do
    end do
  end subroutine cal_mmat

  subroutine cal_grad()
    integer :: i,j,ip
    complex(kind=DP) :: tmat(nocc,nocc)
    ! R * M^T * R+
    !!tmat = matmul(matmul(transpose(dconjg(rmat)),transpose(mmat)),rmat)
    ! {(R * M^T * R+),B}
    !!tmat = tmat*bmat
    ! R+ * {(R * M^T * R+),B} * R
    !!tmat = matmul(matmul(rmat,tmat),transpose(dconjg(rmat)))

    tmat = transpose(mmat)
    do j=2,nocc
       do i=1,j-1
          ip = i + (j-2)*(j-1)/2
          grad(ip) = dble(tmat(j,i)) - dble(tmat(i,j))
       end do
    end do

    gradmax = 0.d0
    do i=1,ntri
       gradmax = max(gradmax,abs(grad(i)))
    end do
  end subroutine cal_grad

  subroutine update_amat()
    if(wannier_opt_method == SD) then 
       call sd_method
    else if(wannier_opt_method == CG) then
       call cg_method
    else if(wannier_opt_method == CGPRC) then
       call cgprc_method
    end if
  contains

    subroutine sd_method
      real(kind=DP), save :: dt
      real(kind=DP), allocatable, save :: amat0(:)
      real(kind=DP), allocatable, save :: grad0(:)
      integer, save :: nold = 0
      real(kind=DP) :: yd
      if(.not.allocated(amat0)) allocate(amat0(ntri))
      if(.not.allocated(grad0)) allocate(grad0(ntri))
      if(nold == 0) then
         !!dt = 0.d0
         dt = dt_wan
         amat0 = amat
         grad0 = grad
         nold = 1
      else
         dt = dt + dt_wan
      end if
      yd = -dot_product(grad,grad0)/PAI2**2
      write(6,'(3(1x,f15.5))') dt, omega, yd
      !!call line_search(dt,nold,grad0)
      !!amat = amat0 + dt*grad0
      !!write(6,*) 'dt = ',dt,' nold=',nold
      !!amat = amat0 + dt_wan*grad
      amat = dt_wan*grad
    end subroutine sd_method

    subroutine cg_method
    end subroutine cg_method

    subroutine cgprc_method
    end subroutine cgprc_method

    subroutine line_search(dt,nold,grad0)
      real(kind=DP), intent(inout) :: dt
      integer, intent(inout) :: nold
      real(kind=DP), intent(in) :: grad0(ntri)
      real(kind=DP), allocatable, save :: old_grad(:) ! d(ntri)
      real(kind=DP), save :: old_omega
      real(kind=DP), save :: x0,x1
      real(kind=DP) :: x_opt,y0,y1,yd0,yd1
      real(kind=DP), parameter :: eps = 1.d-2

      if(.not.allocated(old_grad)) allocate(old_grad(ntri))

      if(nold>=1) then
         if(nold==1.or.nold==2) then
            y0 = old_omega
            y1 = omega
            yd0 = -dt_wan*dot_product(old_grad,grad0)/PAI2**2
            yd1 = -dt_wan*dot_product(grad,grad0)/PAI2**2
         else if(nold==3) then
            y0 = omega
            y1 = old_omega
            yd0 = -dt_wan*dot_product(grad,grad0)/PAI2**2
            yd1 = -dt_wan*dot_product(old_grad,grad0)/PAI2**2
         end if
         if(yd0 < 0.d0 .and. yd1 > 0.d0) then
            call parabolic_fit(x_opt,y0,y1,yd0)
            write(6,*) 'x0,x1=',x0,x1
            write(6,*) 'y0,y1=',y0,y1
            write(6,*) 'yd0,yd1=',yd0,yd1
            write(6,*) 'x_opt=',x_opt
            dt = x_opt*(x1-x0)+x0
            if(y0<y1) then
               x1=dt
            else
               x0=dt
            end if
            dt = dt*dt_wan
            if(abs(x0-x1) < eps*max(x0,x1)) then
               nold = 0
               return
            else
               nold = 2
            end if
         else
            dt = dt + dt_wan
         end if
      else
         dt = dt + dt_wan
      end if

      if(nold==0) then
         x0 = 0.d0
         x1 = 1.d0
      end if
      if(old_omega > omega .or. nold==0 .or. nold==1) then
         old_omega  = omega
         old_grad = grad
         if(nold==0)then
            nold = 1
         else if(nold==2) then
            nold = 2
         end if
      else
         nold = 3
      end if
    end subroutine line_search

    subroutine parabolic_fit(xmin,y0,y1,yd0)
      real(kind=DP), intent(out) :: xmin
      real(kind=DP), intent(in)  :: y0,y1,yd0
      real(kind=DP) :: a, b, c
      c = y0
      b = yd0
      a = y1-(b+c)
      xmin = -0.5d0*b/a
    end subroutine parabolic_fit

  end subroutine update_amat

  logical function Omega_is_max(nfout)
    integer, intent(in) :: nfout

    Omega_is_max = .true. 
    if(gradmax > eps_wan) Omega_is_max = .false. 
    if(iteration >= max_iter_wan) Omega_is_max = .true.

    if(printable) then
       write(nfout,'("!Wannier: ",i5," Omega= ",f18.5," Gradmax= ",f18.5)') iteration,omega,gradmax
    end if
  end function Omega_is_max

  subroutine cal_Omega()
    integer :: i,iw
    real(kind=DP) :: zr,zi
    omega = 0.d0
    do iw=1,nwght
       do i=1,nocc
          zr = dble(zn(i,iw))
          zi = dimag(zn(i,iw))
          omega = omega + wght(iw)*(1.d0-(zr*zr+zi*zi))
       end do
    end do
    omega = omega/PAI2**2
  end subroutine cal_Omega

  subroutine wannier_center(nfout)
    integer, intent(in) :: nfout
    integer :: n,i
    real(kind=DP) :: wfc(nocc,3)

    ! wannier center
    wfc = atan2(dimag(zn(:,1:3)),dble(zn(:,1:3)))/PAI2
    do i=1,3
       do n=1,nocc
          if(wfc(n,i)<-1.d-4) wfc(n,i) = wfc(n,i)+1.d0
       end do
    end do
    if(printable) then
      write(nfout,*) '== Wannier centers (in internal) =='
      do n=1,nocc
         write(nfout,'(i5,3(1x,f18.10))') n,wfc(n,1:3)
      end do
    end if
    wfc = matmul(wfc,transpose(altv))
    if(printable) then
      write(nfout,*) '== Wannier centers (in Cartesian) =='
      do n=1,nocc
         write(nfout,'(i5,3(1x,f18.10))') n,wfc(n,1:3)
      end do
    end if

  end subroutine wannier_center

  subroutine wannier_function(nfout)
    integer, intent(in) :: nfout
    integer :: nspin,ispin,ib,ib1,jb1,i,ik
    integer :: ip,nffth
    real(kind=DP) :: wf(nfft/2), wf1(nfft)

    call m_FFT_alloc_WF_work()

    nffth = nfft/2

    nspin = 1
    do ispin=1,nspin
       ik = ispin
       do ib=1,nocc
          call m_Files_open_nfwannier(nspin,ispin,ib)
          wf = 0.d0
          do ib1 = 1,nocc
             jb1 = neordr(ib1,ik)
             call m_ES_WF_in_Rspace(ik,jb1,wf1)
             do i=1,nffth
                ip = 2*i-1
                wf(i) = wf(i) + wf1(ip)*umat(ib1,ib)
             end do
          end do
          call wd_wf(nfout,nfwannier,nspin,ispin,ib)
       end do
    end do

    call m_FFT_dealloc_WF_work()

  contains

    subroutine wd_wf(nfout,nfwf,nspin,ispin,ib)
      integer, intent(in) :: nfout,nfwf,nspin,ispin,ib
      integer :: i,j,k, id, nl, nm, nn, nlhf,inew,jnew,knew,ip,mm
      real(kind=DP),allocatable,dimension(:,:,:) :: wkwf
      real(kind=DP) ::      s1, s2, sratio,x,y,z
      integer, parameter :: UP = 1 , DOWN = 2
      integer ::            up_down
      real(kind=DP),allocatable,dimension(:,:) :: cps_full
      integer, allocatable,dimension(:) :: ityp_full
      real(kind=DP), dimension(3) :: r_wk    !!! K.Mae 040315
      integer :: ucret, m   !!! K.Mae 040315
      real(kind=DP) :: norm
      integer :: n1,n2,n3
      real(kind=DP) :: dn1,dn2,dn3
      real(kind=DP) :: eig
      integer :: nk 

      nk=1
      eig=0.d0

      id = fft_box_size_WF(1,0)
      mm = fft_box_size_WF(2,0)
      nl = fft_box_size_WF(1,1)
      nm = fft_box_size_WF(2,1)
      nn = fft_box_size_WF(3,1)

      if(kimg == 1) then
         nlhf = id/2
      else
         nlhf = id
      end if

      if(wannier_filetype == DENSITY_ONLY &
      & .or. wannier_filetype == VTK &
      & .or. wannier_filetype == BINARY) then
         allocate(wkwf(nl,nm,nn)); wkwf = 0.d0
      else if(wannier_filetype == CUBE) then
         allocate(wkwf(nn,nm,nl)); wkwf = 0.d0
      end if

      if(nspin==2) then
         if(ispin == 1) then
            up_down = UP
         else
            up_down = DOWN
         end if
      end if

      if(ipri >= 2) write(nfout,9001) nl*nm*nn, nl, nm, nn
9001  format(' Wannier function ',i8,'(',3i5,')')

      if(ipri >= 2) write(nfout,*) ' !D FFT cube mapping start'
      do i = 1, nm
         do j = 1, nn
            do k = 1, nl
               if(kimg == 1 .and. k > nlhf) then
                  knew = id - k
                  jnew = nn+2 - j
                  inew = nm+2 - i
                  if(jnew > nn) then
                     jnew = jnew - nn
                  end if
                  if(inew > nm) then
                     inew = inew - nm
                  end if
               else
                  knew = k; jnew = j; inew = i
               end if
               ip = nlhf*mm*(jnew-1) + nlhf*(inew-1) + knew
               if(wannier_filetype == DENSITY_ONLY &
               & .or. wannier_filetype == VTK &
               & .or. wannier_filetype == BINARY) then
                  wkwf(k,i,j) = wf(ip)
               else if(wannier_filetype == CUBE) then
                  wkwf(j,i,k) = wf(ip)
               end if
            end do
         end do
      end do
      norm = 1.d0/dsqrt(univol)
      do i = 1, nm
         do j = 1, nn
            do k = 1, nl
               if(wannier_filetype == DENSITY_ONLY &
               & .or. wannier_filetype == VTK &
               & .or. wannier_filetype == BINARY) then
                  wkwf(k,i,j) = norm*wkwf(k,i,j)
               else if(wannier_filetype == CUBE) then
                  wkwf(j,i,k) = norm*wkwf(j,i,k)
               end if
            end do
         end do
      end do

      if(wannier_filetype == DENSITY_ONLY) then
         write(nfwf,9001) nl*nm*nn, nl, nm, nn
         write(nfwf,'(6e13.5)') wkwf
      else if(wannier_filetype == BINARY) then
         write(nfwf) nl*nm*nn, nl, nm, nn
         write(nfwf) altv, nspin, up_down, nk, ib, eig
         write(nfwf) wkwf(:,:,:)
      else if(wannier_filetype == VTK) then
         write(nfwf,'("# vtk DataFile Version 2.0")')
         if(nspin == 2) then
            if(up_down == 1) then
               write(nfwf,'(" Wannier function UP : n=",i7)') ib
            else
               write(nfwf,'(" Wannier function DOWN : n=",i7)') ib
            end if
         else
            write(nfwf,'(" Wannier function : n=",i7)') ib
         end if
         write(nfwf,'("ASCII")')
         write(nfwf,'("DATASET STRUCTURED_GRID")')
         write(nfwf,'("DIMENSIONS",3(1x,i5))') nl+1,nm+1,nn+1
         write(nfwf,'("POINTS",1x,i10,1x,"float")') (nl+1)*(nm+1)*(nn+1)
         do n1=0,nl
            do n2=0,nm
               do n3=0,nn
                  dn1 = n1/dble(nl)
                  dn2 = n2/dble(nm)
                  dn3 = n3/dble(nn)
                  x = altv(1,1)*dn1 + altv(1,2)*dn2 + altv(1,3)*dn3
                  y = altv(2,1)*dn1 + altv(2,2)*dn2 + altv(2,3)*dn3
                  z = altv(3,1)*dn1 + altv(3,2)*dn2 + altv(3,3)*dn3
                  write(nfwf,'(3(1x,e13.5))') x,y,z
               end do
            end do
         end do
         write(nfwf,'("")')
         write(nfwf,'("POINT_DATA",1x,i10)') (nl+1)*(nm+1)*(nn+1)
         write(nfwf,'("SCALARS scalars float")')
         write(nfwf,'("LOOKUP_TABLE default")')
         do n1=0,nl
            i=n1+1
             if(n1==nl) i=1
            do n2=0,nm
               j=n2+1
               if(n2==nm) j=1
               do n3=0,nn
                  k=n3+1
                  if(n3==nn) k=1
                  write(nfwf,'(e13.5)') wkwf(i,j,k)
               end do
            end do
         end do
      else if(wannier_filetype == CUBE) then
         !!if(len_trim(wf_title) >= 1) then
         !!   write(nfwf,*) trim(wf_title)
         !!else
            write(nfwf,'(" Calculated by phase")')
         !!end if
         if(nspin == 2) then
            if(up_down == 1) then
               write(nfwf,'(" Wannier function UP : n=",i7)') ib
            else
               write(nfwf,'(" Wannier function DOWN : n=",i7)') ib
            end if
         else
            write(nfwf,'(" Wannier function : n=",i7)') ib
         end if
         x = 0.d0; y = 0.d0; z = 0.d0
         write(nfwf,'(i6,3f10.4)') -natm2, x,y,z
         do i = 1, 3
            write(nfwf,'(i6,3f10.6)') fft_box_size_WF(i,1), altv(1:3,i)/dble(fft_box_size_WF(i,1))
         end do

         allocate(cps_full(natm2,3))
         allocate(ityp_full(natm2))
         call m_IS_pack_all_ions_in_uc(ityp_full,cps_full)
         do i = 1, natm2
            m = ityp_full(i)
            write(nfwf,'(i6,4f10.6)') nint(iatomn(m)), ival(m), cps_full(i,1:3)
         end do
         deallocate(ityp_full,cps_full)

         write(nfwf,'(10i5)') 1,1
         write(nfwf,'(6e13.5)') wkwf(:,:,:)
  
      end if
      if(allocated(wkwf)) deallocate(wkwf)
    end subroutine wd_wf

  end subroutine wannier_function

  subroutine constract_of_ftq()
    ! local varialbes
    integer :: it,ik,ip,lmt1,lmt2,u,v
    integer :: il1,il2,tau1,tau2,l3
    integer :: n,ilm3,iiqitg,iw
    integer :: mdvdb
    integer, allocatable :: il3(:)
    real(kind=DP) :: fac
    real(kind=DP) :: dk(3),ylm,dga
    real(kind=DP) :: ftqr(nwght),ftqi(nwght),ftqb
    complex(kind=DP) :: zi = (0.d0,1.d0)

    call m_PP_find_maximum_l(n)   !  n-1: maximum l
    n = (n-1) + (n-1) + 1
    allocate(il3(n**2)); call substitute_il3(n**2,il3) ! -(b_Elec..)
    allocate(ftqe(nlmt,nlmt,ntyp,nwght)); ftqe=(0.d0,0.d0)
    TYPE: do it = 1, ntyp
       mdvdb = m_PP_include_vanderbilt_pot(it)
       if(mdvdb == SKIP) cycle TYPE
       LMT_1: do lmt1=1,ilmt(it)
          il1 = ltp(lmt1,it)
          tau1 = taup(lmt1,it)
          LMT_2: do lmt2=lmt1,ilmt(it)
             il2 = ltp(lmt2,it)
             tau2 = taup(lmt2,it)
             ftqr(1:nwght) = 0.d0
             ftqi(1:nwght) = 0.d0
             LM3: do n=1,il2p(lmt1,lmt2,it)
                ilm3 = isph(lmt1,lmt2,n,it); l3=il3(ilm3)
                iiqitg = iqitg(il1,tau1,il2,tau2,l3+1,it)
                if(iiqitg == 0) cycle LM3
                do iw=1,nwght
                   dk(1:3)=matmul(rltv,ghat(1:3,iw))
                   call sphrp2_for_Berry(ilm3,dk,ylm)
                   ftqb = qitg_wan(iiqitg,iw)*dl2p(lmt1,lmt2,n,it)*ylm
                   !!$ write(6,*) 'it,lmt1,lmt2,n,iw,ftqb=',it,lmt1,lmt2,n,iw,ftqb
                   if(mod(l3,2)==0) then
                      !!ftqr(iw)=ftqr(iw)+real(zi**(-l3))*ftqb
                      ftqr(iw)=ftqr(iw)+real(zi**(l3))*ftqb
                   else
                      !!ftqi(iw)=ftqi(iw)+dimag(zi**(-l3))*ftqb
                      ftqi(iw)=ftqi(iw)+dimag(zi**(l3))*ftqb
                   end if
                end do
             end do LM3
             do iw=1,nwght
                ftqe(lmt1,lmt2,it,iw) = dcmplx(ftqr(iw),ftqi(iw))
                ftqe(lmt2,lmt1,it,iw) = ftqe(lmt1,lmt2,it,iw)
             end do
          end do LMT_2
       end do LMT_1
    end do TYPE

    deallocate(il3)
  end subroutine constract_of_ftq

  subroutine rd_zmat_and_umat(nfout,nfcnt)
    integer, intent(in) :: nfout,nfcnt
    write(nfout,*) 'Reading Z and U matrices'
    rewind(nfcnt)
    read(nfcnt) zmat,umat
  end subroutine rd_zmat_and_umat

  subroutine wd_zmat_and_umat(nfout,nfcnt)
    integer, intent(in) :: nfout,nfcnt
    write(nfout,*) 'Writing Z and U matrices'
    rewind(nfcnt)
    write(nfcnt) zmat,umat
  end subroutine wd_zmat_and_umat

end module m_Wannier
