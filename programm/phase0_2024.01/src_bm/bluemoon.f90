!=======================================================================
!
!  PROGRAM  PHASE/0 2014.01 (rev.375)
!
!  "First-principles Electronic Structure Calculation Program"
!
!  AUTHOR(S): J. Koga 
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
!***************************************************************
module m_parameter
  use m_Const_Parameters, only : DP
  implicit none
  real(DP), parameter :: BOHR2ANG   = 0.5291772480d0
  real(DP), parameter :: Hartree2eV = 27.21139615d0
  real(DP), parameter :: autime2fs  = 2.418884327d-2
  real(DP), parameter :: HartreeBohr2eVAng = Hartree2eV/BOHR2ANG
  real(DP), parameter :: eV2kcal_mol = 627.509391d0/Hartree2eV
  real(DP), parameter :: eV2kJ_mol   = 2625.5d0/Hartree2eV
  real(DP), parameter :: erg2eV = 6.24150974d+11

  character(len=256) :: inpfile_name = 'nfinp.data'

end module m_parameter

! module for spline interpolation;
! from "numerical recipies in Fortran"
module m_spline
  use m_Const_Parameters, only : DP
  implicit none

  interface init_bicubic_spline
    module procedure init_bicubic_spline_noder
    module procedure init_bicubic_spline_der
  end interface init_bicubic_spline

  interface init_tricubic_spline
    module procedure init_tricubic_spline_noder
    module procedure init_tricubic_spline_der
  end interface init_tricubic_spline

  contains

  ! initialization routine for a 1D spline interpolation
  subroutine init_cubic_spline(x,y,n,yp1,ypn,y2)
    implicit none
    integer,intent(in) :: n
    real, intent(in) :: yp1,ypn
    real(DP), dimension(n), intent(in)  :: x,y
    real(DP), dimension(n), intent(out) :: y2
    integer, parameter :: NMAX=500
    integer i,k
    real(DP) p,qn,sig,un
    real(DP), dimension(NMAX) :: u

    y2=0.0d0
    u=0.0d0

    if(yp1.gt..99e30) then
        y2(1) = 0.
        u(1) = 0.
    else
        y2(1)=-0.5
        u(1) = (3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    endif

    do i=2,n-1
      sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
      p=sig*y2(i-1)+2
      y2(i)=(sig-1)/p
      u(i)=(6*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) &
      & /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    enddo

    if ( ypn.gt..99e30) then
        qn=0.
        un=0.
    else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
    endif

    y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1)
    do k=n-1,1,-1
      y2(k) = y2(k)*y2(k+1)+u(k)
    enddo
  end subroutine init_cubic_spline

  ! get the spline-interpolated value and its derivative
  ! using the second-derivative matrix
  ! obtained from init_cubic_spline.
  ! the return value is the interpolated value, while the last
  ! dummy argument is the derivative.
  real(DP) function cubic_spline(xa,ya,y2a,n,x,dy)
    integer, intent(in) :: n
    real(DP), intent(in) :: x
    real(DP), intent(out) :: dy
    real(DP), dimension(n), intent(in) :: xa,y2a,ya

    integer :: k,khi,klo
    real(DP) :: a,b,h

    klo=1
    khi=n

    do while ( (khi-klo).gt.1 ) 
      k=(khi+klo)/2
      if (xa(k).gt.x) then
          khi=k
      else
          klo=k
      endif
    enddo

    h=xa(khi)-xa(klo)
    if ( h.eq.0. ) stop 'bad xa input in splint'
    a=(xa(khi)-x)/h
    b=(x-xa(klo))/h
    dy = (ya(khi)-ya(klo))/h - (3*a**2-1)/6.0d0 * h * y2a(klo) + &
    & (3*b**2-1)/6.0d0 * h * y2a(khi)
    cubic_spline=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6
  end function cubic_spline

  ! initialization routine for a 2D cubic-spline interpolation
  subroutine init_bicubic_spline_noder(x1a,x2a,ya,m,n,y2a)
    implicit none
    integer, intent(in) :: m,n
    real(DP), dimension(m), intent(in) :: x1a
    real(DP), dimension(n), intent(in) :: x2a
    real(DP), dimension(m,n), intent(in) :: ya
    real(DP), dimension(m,n), intent(out) :: y2a
    integer :: i,j,k
    real(DP), dimension(n) :: y2tmp,ytmp

    y2tmp=0.0d0;ytmp=0.0d0;y2a=0.0d0

    do j=1,m
      do k=1,n
        ytmp(k) = ya(j,k)
      enddo
      call init_cubic_spline(x2a,ytmp,n,1.e30,1.e30,y2tmp)
      do k=1,n
        y2a(j,k)=y2tmp(k)
      enddo
    enddo
  end subroutine init_bicubic_spline_noder

  subroutine init_bicubic_spline_der(x1a,x2a,ya,yp,m,n,y2a)
    integer, intent(in) :: m,n
    real(DP), dimension(m), intent(in) :: x1a
    real(DP), dimension(n), intent(in) :: x2a
    real(DP), dimension(m,n), intent(in) :: ya
    real(DP), dimension(m,n), intent(in) :: yp
    real(DP), dimension(m,n), intent(out) :: y2a
    real(DP), dimension(n) :: y2tmp,ytmp,yptmp,yp1,ypn
    integer :: i,j,k

    y2tmp=0.0d0;ytmp=0.0d0;y2a=0.0d0;yptmp=0.0d0

    do j=1,m
      do k=1,n
        ytmp(k) = ya(j,k)
      enddo
!      call init_cubic_spline(x2a,ytmp,n,,,y2tmp)
      do k=1,n
        y2a(j,k)=y2tmp(k)
      enddo
    enddo
  end subroutine init_bicubic_spline_der

  ! get the 2D-spline-interpolated value using the second-derivative matrix
  ! obtained from init_bicubic_spline
  real(DP) function bicubic_spline_fast(x1a,x2a,ya,y2a,m,n,x1,x2,dy1,dy2)
    implicit none
    integer, intent(in) :: m,n
    real(DP), dimension(m), intent(in) :: x1a
    real(DP), dimension(n), intent(in) :: x2a
    real(DP), dimension(m,n), intent(in) :: ya,y2a
    real(DP), intent(in) :: x1,x2
    real(DP), intent(out) :: dy1,dy2

    integer j,k
    real(DP) :: dytmp
    real(DP), dimension(n) :: y2tmp, ytmp,yytmp,dytmpmat,y2tmptmp

    y2tmp=0.0d0;ytmp=0.0d0;yytmp=0.0d0;dytmpmat=0.0d0;y2tmptmp=0.0d0

    do j=1,m
      ytmp(:) = ya(j,:)
      y2tmp(:) = y2a(j,:)
      yytmp(j) = cubic_spline(x2a,ytmp,y2tmp,n,x2,dytmp)
      dytmpmat(j) = dytmp
      y2tmptmp(j) = cubic_spline(x2a,y2tmp,y2tmp,n,x2,dytmp)
    enddo

    dy1=0.0d0;dy2=0.0d0

!    y2tmp=0.0d0
    y2tmp = y2tmptmp
!    call init_cubic_spline(x1a,dytmpmat,m,1.e30,1.e30,y2tmp)
    dy2 = cubic_spline(x1a,dytmpmat,y2tmp,m,x1,dytmp)

!    y2tmp=0.0d0
!    call init_cubic_spline(x1a,yytmp,m,1.e30,1.e30,y2tmp)
    bicubic_spline_fast = cubic_spline(x1a,yytmp,y2tmp,m,x1,dy1)
  end function bicubic_spline_fast

  ! get the 2D-spline-interpolated value using the second-derivative matrix
  ! obtained from init_bicubic_spline
  real(DP) function bicubic_spline(x1a,x2a,ya,y2a,m,n,x1,x2,dy1,dy2)
    implicit none
    integer, intent(in) :: m,n
    real(DP), dimension(m), intent(in) :: x1a
    real(DP), dimension(n), intent(in) :: x2a
    real(DP), dimension(m,n), intent(in) :: ya,y2a
    real(DP), intent(in) :: x1,x2
    real(DP), intent(out) :: dy1,dy2

    integer j,k
    real(DP) :: dytmp
    real(DP), dimension(n) :: y2tmp, ytmp,yytmp,dytmpmat

    y2tmp=0.0d0;ytmp=0.0d0;yytmp=0.0d0;dytmpmat=0.0d0

    do j=1,m
      ytmp(:) = ya(j,:)
      y2tmp(:) = y2a(j,:)
      yytmp(j) = cubic_spline(x2a,ytmp,y2tmp,n,x2,dytmp)
      dytmpmat(j) = dytmp
    enddo

    dy1=0.0d0;dy2=0.0d0

!    y2tmp=0.0d0
    call init_cubic_spline(x1a,dytmpmat,m,1.e30,1.e30,y2tmp)
    dy2 = cubic_spline(x1a,dytmpmat,y2tmp,m,x1,dytmp)

!    y2tmp=0.0d0
    call init_cubic_spline(x1a,yytmp,m,1.e30,1.e30,y2tmp)
    bicubic_spline = cubic_spline(x1a,yytmp,y2tmp,m,x1,dy1)
  end function bicubic_spline

  subroutine init_tricubic_spline_der
  end subroutine init_tricubic_spline_der

  ! initialization routine for a 3D cubic-spline interpolation 
  ! the routine is actually exactly the same for the 2D spline.
  subroutine init_tricubic_spline_noder(x1a,x2a,x3a,ya,m,n,o,y2a)
    implicit none
    integer, intent(in) :: m,n,o
    real(DP), dimension(m), intent(in) :: x1a
    real(DP), dimension(n), intent(in) :: x2a
    real(DP), dimension(o), intent(in) :: x3a
    real(DP), dimension(m,n,o), intent(in) :: ya
    real(DP), dimension(m,n,o), intent(out) :: y2a

    integer :: i
    real(DP), dimension(m,n) :: y2atmp,yatmp

    y2atmp=0.0d0;yatmp=0.0d0;y2a=0.0d0

    do i=1,o
      yatmp(:,:) = ya(:,:,i)
      call init_bicubic_spline(x1a,x2a,yatmp,m,n,y2atmp)
      y2a(:,:,i) = y2atmp(:,:)
    enddo
  end subroutine init_tricubic_spline_noder

  real(DP) function tricubic_spline_2d_precal(&
  & x3a,y2d,dy2d1,dy2d2,x3,dy1,dy2,dy3,o)
    implicit none
    integer, intent(in) :: o
    real(DP), intent(in), dimension(o) :: y2d,dy2d1,dy2d2
    real(DP), dimension(o), intent(in) :: x3a
    real(DP), intent(in) :: x3
    real(DP), intent(out) :: dy1,dy2,dy3

    real(DP), dimension(o) :: y2tmp
    real(DP) :: dy1tmp,dy2tmp

    y2tmp=0.0d0;dy3=0.0d0

    call init_cubic_spline(x3a,dy2d1,o,1.e30,1.e30,y2tmp)
    dy1 = cubic_spline(x3a,dy2d1,y2tmp,o,x3,dy1tmp)

    call init_cubic_spline(x3a,dy2d2,o,1.e30,1.e30,y2tmp)
    dy2 = cubic_spline(x3a,dy2d2,y2tmp,o,x3,dy2tmp)

    call init_cubic_spline(x3a,y2d,o,1.e30,1.e30,y2tmp)
    tricubic_spline_2d_precal = &
    & cubic_spline(x3a,y2d,y2tmp,o,x3,dy3)
  end function tricubic_spline_2d_precal

  ! get the 3D-spline-interpolated value using the second-derivative matrix
  ! obtained from init_tricubic_spline
  ! we first do bi-cubic spline interpolation for the first two index, and from
  ! the interpolated array thus obtained, perform a final interpolation by the
  ! third index.
  real(DP) function tricubic_spline(x1a,x2a,x3a,ya,y2a,m,n,o,x1,x2,x3,&
  & dy1,dy2,dy3)
    implicit none
    integer, intent(in) :: m,n,o
    real(DP), dimension(m), intent(in) :: x1a
    real(DP), dimension(n), intent(in) :: x2a
    real(DP), dimension(o), intent(in) :: x3a
    real(DP), dimension(m,n,o), intent (in) :: ya
    real(DP), dimension(m,n,o), intent(in) :: y2a
    real(DP), intent(in) :: x1,x2,x3
    real(DP), intent(out) :: dy1,dy2,dy3

    integer j,k
    real(DP), dimension(o) :: y2tmp
    real(DP), dimension(o) :: yytmp
    real(DP), dimension(o) :: dy1tmpmat
    real(DP), dimension(o) :: dy2tmpmat
    real(DP) :: dy1tmp,dy2tmp
    real(DP), dimension(m,n) :: ytmp,y2atmp

    y2tmp=0.0d0;yytmp=0.0d0;dy1tmpmat=0.0d0
    dy2tmpmat=0.0d0;ytmp=0.0d0;y2atmp=0.0d0

    dy1=0.0d0;dy2=0.0d0;dy3=0.0d0
    do j=1,o
      ytmp(:,:) = ya(:,:,j)
      y2atmp(:,:) = y2a(:,:,j)
      yytmp(j) = bicubic_spline(x1a,x2a,ytmp,y2atmp,m,n,x1,x2,&
      & dy1tmp,dy2tmp)
      dy1tmpmat(j) = dy1tmp
      dy2tmpmat(j) = dy2tmp
    enddo

    call init_cubic_spline(x3a,dy1tmpmat,o,1.e30,1.e30,y2tmp)
    dy1 = cubic_spline(x3a,dy1tmpmat,y2tmp,o,x3,dy1tmp)

    call init_cubic_spline(x3a,dy2tmpmat,o,1.e30,1.e30,y2tmp)
    dy2 = cubic_spline(x3a,dy2tmpmat,y2tmp,o,x3,dy2tmp)

    call init_cubic_spline(x3a,yytmp,o,1.e30,1.e30,y2tmp)
    tricubic_spline = cubic_spline(x3a,yytmp,y2tmp,o,x3,dy3)
  end function tricubic_spline

  ! get the 3D-spline-interpolated value using the second-derivative matrix
  ! obtained from init_tricubic_spline
  ! we first do bi-cubic spline interpolation for the first two index, and from
  ! the interpolated array thus obtained, perform a final interpolation by the
  ! third index.
  real(DP) function tricubic_spline_fast(x1a,x2a,x3a,ya,y2a,m,n,o,x1,x2,x3,&
  & dy1,dy2,dy3)
    implicit none
    integer, intent(in) :: m,n,o
    real(DP), dimension(m), intent(in) :: x1a
    real(DP), dimension(n), intent(in) :: x2a
    real(DP), dimension(o), intent(in) :: x3a
    real(DP), dimension(m,n,o), intent (in) :: ya
    real(DP), dimension(m,n,o), intent(in) :: y2a
    real(DP), intent(in) :: x1,x2,x3
    real(DP), intent(out) :: dy1,dy2,dy3

    integer j,k
    real(DP), dimension(o) :: y2tmp
    real(DP), dimension(o) :: y2tmptmp
    real(DP), dimension(o) :: yytmp
    real(DP), dimension(o) :: dy1tmpmat
    real(DP), dimension(o) :: dy2tmpmat
    real(DP) :: dy1tmp,dy2tmp
    real(DP), dimension(m,n) :: ytmp,y2atmp

    y2tmp=0.0d0;yytmp=0.0d0;dy1tmpmat=0.0d0
    dy2tmpmat=0.0d0;ytmp=0.0d0;y2atmp=0.0d0
    y2tmptmp=0.0d0

    dy1=0.0d0;dy2=0.0d0;dy3=0.0d0
    do j=1,o
      ytmp(:,:) = ya(:,:,j)
      y2atmp(:,:) = y2a(:,:,j)
      yytmp(j) = bicubic_spline_fast(x1a,x2a,ytmp,y2atmp,m,n,x1,x2,&
      & dy1tmp,dy2tmp)
      dy1tmpmat(j) = dy1tmp
      dy2tmpmat(j) = dy2tmp
      y2tmptmp(j) = bicubic_spline_fast(x1a,x2a,y2atmp,y2atmp,m,n,x1,x2,&
      & dy1tmp,dy2tmp)
    enddo

    y2tmp = y2tmptmp
!    call init_cubic_spline(x3a,dy1tmpmat,o,1.e30,1.e30,y2tmp)
    dy1 = cubic_spline(x3a,dy1tmpmat,y2tmp,o,x3,dy1tmp)

!    call init_cubic_spline(x3a,dy2tmpmat,o,1.e30,1.e30,y2tmp)
    dy2 = cubic_spline(x3a,dy2tmpmat,y2tmp,o,x3,dy2tmp)

!    call init_cubic_spline(x3a,yytmp,o,1.e30,1.e30,y2tmp)
    tricubic_spline_fast = cubic_spline(x3a,yytmp,y2tmp,o,x3,dy3)
  end function tricubic_spline_fast

end module m_spline

module m_bm
use m_Const_Parameters, only : NOCONV
use m_parameter, only : BOHR2ANG,Hartree2eV,autime2fs,HartreeBohr2eVAng, eV2kJ_mol, eV2kcal_mol, inpfile_name
use m_spline, only : init_cubic_spline,cubic_spline
implicit none

integer :: nsample = 1
integer :: ndim   = 1
integer :: nsteps = 2000
integer :: nsteps_per_sample=2000
integer :: nequib = 1000
integer :: nreac_coords = 1
integer :: max_reac_coords = 1
integer :: istart_reac_coords = 1
character(256) :: basedir='.'

real(8), allocatable, dimension(:,:) :: mean_force
real(8), allocatable, dimension(:) :: mean_force_av
real(8), allocatable, dimension(:) :: value

real(8) :: d_xsi
integer :: n_xsi = 1000

logical :: smooth = .false.

contains
  subroutine parse_input()
    integer :: iarg
    integer :: iargc
    integer :: iret
    integer :: dret
    character(len=256)    :: string
    logical :: ex
    integer :: f_openInputFile,f_getIntValue,f_selectBlock,f_selectParentBlock,f_closeInputFile,&
   & f_selectTop, f_getRealValue, f_getStringValue, f_getBoolValue
    logical :: bret
    iarg = iargc()
    if (iarg/=0) then
      call getarg(1,inpfile_name)
    endif

    inquire(file=trim(adjustl(inpfile_name)),exist=ex)

    if(.not.ex)then
      write(0,'(a)') 'input file : '//trim(adjustl(inpfile_name))//' does not exist.'
      stop
    endif
    iret = f_openInputFile(trim(adjustl(inpfile_name)))
    iret = f_selectTop()
    if(f_selectBlock("thermodynamic_integration")==0)then
      if(f_getIntValue('ndim',iret)==0)then
        ndim = iret
      endif
      if(f_getIntValue('nsteps',iret)==0)then
        nsteps = iret
      endif
      if(f_getIntValue('nequib',iret)==0)then
        nequib = iret
        if(nequib.ge.nsteps)then
          write(0,'(a)') 'nequib must be smaller than nsteps'
          write(0,'(a)') 'setting nequib = 0.5 * nsteps'
          nequib = nint(0.5d0*real(nsteps))
        endif
      endif
      if(f_getIntValue('nsample',iret)==0)then
        nsample = iret
        if(nsample<1)then
          nsample=1
          write(0,*) 'nsample must be larger then 1'
        endif
      endif
      if(f_getIntValue('nreac_coords',iret)==0)then
        nreac_coords = iret
        max_reac_coords = nreac_coords
      endif
      if(f_getIntValue('max_nreac_coords',iret)==0)then
        if(iret>nreac_coords) max_reac_coords = iret
      endif
      if(f_getIntValue('istart_reac_coords',iret)==0)then
        istart_reac_coords = iret
        if(istart_reac_coords.ge.nreac_coords)then
          write(0,*) 'istart_reac_coords must be smaller than nreac_coords'
          istart_reac_coords=1
        endif
      endif
      if(f_getIntValue('n_xsi',iret)==0)then
        n_xsi = iret
      endif
      if(f_getBoolValue('smooth',bret)==0)then
        smooth = bret
      endif
      allocate(value(nreac_coords-(istart_reac_coords-1)))
      allocate(mean_force_av(nreac_coords-(istart_reac_coords-1)))
      allocate(mean_force(nsample,nreac_coords-(istart_reac_coords-1)))
      if(f_getStringValue('basedir',string,NOCONV)==0)then
        basedir = trim(adjustl(string))
      endif
      iret = f_selectParentBlock()
    else
      write(0,*) 'could not find the thermodynamic_integration block'
      stop
    endif
    open(unit=11,file=trim(basedir)//'/mean_force_raw.data',status='replace')
    if(smooth) open(unit=12,file=trim(basedir)//'/mean_force_smoothed.data',status='replace')
    open(unit=13,file=trim(basedir)//'/potential_of_mean_force.data',status='replace')
    write(13,'(a)') '#value, potetial of mean force in Hartree, eV, kcal/mol, kJ/mol'
    nsteps_per_sample = (nsteps-nequib)/nsample
  end subroutine parse_input

  subroutine do_mean_force(ireac)
    integer, intent(in) :: ireac
    integer :: i,j
    character(len=256) :: fname
    character(len=10) :: tmp
    logical :: ex
    integer :: ifile=10
    integer :: istep
    integer :: iindex
    real(8) :: sumlambda, lamb, val,sigma,dmetric,sumdmetric
    real(8) :: ndiv

    tmp = get_id_char(max_reac_coords,ireac)
    fname = trim(basedir)//'/nfbluemoon.data.reac'//trim(adjustl(tmp))
    inquire(file=trim(fname),exist=ex)
    if(.not.ex)then
      write(0,'(a)') 'data file : '//trim(adjustl(fname))//' does not exist.'
      stop
    endif
    open(ifile,file=trim(fname),status='old')

    do i=1,nequib
      read(ifile,fmt=*,err=2) istep,val,lamb,dmetric
2   continue
    enddo
    ! kara yomi
    iindex = ireac-(istart_reac_coords-1)
    do i=1,nsample
      sumlambda = 0.0d0
      sumdmetric = 0.0d0
      do j=1,nsteps_per_sample
        read(ifile,fmt=*,err=3) istep,val,lamb,dmetric
        sumlambda = sumlambda+lamb/(dsqrt(dmetric))
        sumdmetric = sumdmetric+1.0d0/(dsqrt(dmetric))
3       continue
      enddo
      sumlambda = sumlambda/real(nsteps_per_sample)
      sumdmetric = sumdmetric/real(nsteps_per_sample)
      mean_force(i,iindex) = sumlambda/sumdmetric
    enddo
    value(iindex) = val

    sumlambda = 0.d0
    do i=1,nsample
      sumlambda = sumlambda+mean_force(i,iindex)
    enddo

    mean_force_av(iindex) = sumlambda/real(nsample)

    sigma = 0.d0
    do i=1,nsample
      sigma = sigma+ (mean_force(i,iindex)-sumlambda)**2
    enddo
    sigma = sigma/real(nsample)
    write(11,'(3f20.10)') val,-sumlambda/real(nsample),dsqrt(sigma/real(nsample))

    return
    !!$write(0,'(a)') 'caught error while reading : ',trim(fname)
  end subroutine do_mean_force

  character(len=10) function get_id_char(imax,i)
    integer, intent(in) :: imax,i
    integer :: iketa=2
    integer :: itmp
    character(len=256) :: ch
    character(len=256) :: string
    itmp=int(log10(real(imax)))+1
    if(itmp>2)iketa=itmp
    write(ch,*) iketa
    write(get_id_char,'(i'//trim(adjustl(ch))//'.'//trim(adjustl(ch))//')') i
  end function get_id_char


  subroutine do_thermodynamic_integration()
    integer :: isamp
    integer :: i
    integer :: j
    real(8) :: tmpvi,tmpvj,tmpfi,tmpfj
    real(8),allocatable,dimension(:) :: smoothed_value,smoothed_mean_force
    real(8) :: tmptmp,sumsum,res,tmpval
    real(8) :: inival,sigma
    real(8), allocatable, dimension(:) :: mean_force_2
    real(8), allocatable, dimension(:,:) :: potential
    real(8), allocatable, dimension(:) :: potential_av
    if(smooth)then
      !the 'value' array must be in canonical order in order to perform a cubic spline interpl
      do i=1,nreac_coords-(istart_reac_coords-1)-1
        do j=i+1,nreac_coords-(istart_reac_coords-1)
          if(value(i)>value(j))then
            tmpvi=value(j)
            tmpvj=value(i)
            tmpfi=mean_force_av(j)
            tmpfj=mean_force_av(i)
            value(i)=tmpvi
            value(j)=tmpvj
            mean_force_av(i)=tmpfi
            mean_force_av(j)=tmpfj
          endif
        enddo
      enddo
      allocate(mean_force_2(nreac_coords-(istart_reac_coords-1)))
      d_xsi = (value(nreac_coords-(istart_reac_coords-1))-value(1))/real(n_xsi)
      call init_cubic_spline(value,mean_force_av,nreac_coords-(istart_reac_coords-1),1.e30,1.e30,mean_force_2)
      allocate(smoothed_value(n_xsi));allocate(smoothed_mean_force(n_xsi))
      do i=1,n_xsi+1
        tmpval = value(1) + d_xsi * (i-1)
        res    = cubic_spline(value,mean_force_av,mean_force_2,nreac_coords-(istart_reac_coords-1),tmpval,tmptmp)
        smoothed_value(i) = tmpval
        smoothed_mean_force(i) = res
        write(12,'(2f20.10)') tmpval,-res
      enddo
      do i=1,n_xsi
        sumsum = 0.d0
        do j=i+1,n_xsi-1
          sumsum = sumsum+smoothed_mean_force(j)
        enddo
        sumsum = sumsum+(smoothed_mean_force(i)+smoothed_mean_force(n_xsi))*0.5d0
        sumsum = sumsum * d_xsi
        write(13,'(5f20.10)') smoothed_value(i),sumsum,&
       & sumsum*Hartree2eV,sumsum*Hartree2eV*eV2kcal_mol,sumsum*Hartree2eV*eV2kJ_mol
      enddo
      deallocate(smoothed_value); deallocate(smoothed_mean_force)
      deallocate(mean_force_2)
    else
      n_xsi = nreac_coords-(istart_reac_coords-1)
      d_xsi = (value(nreac_coords-(istart_reac_coords-1))-value(1))/real(n_xsi)
      allocate(potential(nsample,n_xsi));potential=0.d0
      allocate(potential_av(n_xsi));potential_av=0.d0
      do isamp=1,nsample
        do i=1,n_xsi
          sumsum = 0.d0
          do j=i+1,n_xsi-1
            sumsum = sumsum+mean_force(isamp,j)
          enddo
          sumsum = sumsum+(mean_force(isamp,i)+mean_force(isamp,n_xsi))*0.5d0
          sumsum = sumsum * d_xsi
          potential(isamp,i) = sumsum
        enddo
      enddo
      do isamp=1,nsample
        potential_av(:) = potential_av(:) + potential(isamp,:)
      enddo
      potential_av(:) = potential_av(:)/real(nsample)

      do i=1,n_xsi
        sigma = 0.d0
        do j=1,nsample
          sigma = sigma+(potential(j,i)-potential_av(i))**2
        enddo
        sigma = dsqrt(sigma)/real(nsample)
        write(13,'(9f20.10)') value(i),potential_av(i),sigma,&
       & potential_av(i)*Hartree2eV,sigma*Hartree2eV, &
       & potential_av(i)*Hartree2eV*eV2kcal_mol,sigma*Hartree2eV*eV2kcal_mol, &
       & potential_av(i)*Hartree2eV*eV2kJ_mol,sigma*Hartree2eV*eV2kJ_mol
      enddo
      deallocate(potential)
    endif
  end subroutine do_thermodynamic_integration

  subroutine finalize()
    close(11)
    if (smooth) close(12)
    close(13)
    deallocate(value)
    deallocate(mean_force)
  end subroutine finalize

end module m_bm

program bluemoon
  use m_bm, only : parse_input, istart_reac_coords,nreac_coords,do_mean_force &
  &              , do_thermodynamic_integration, finalize
  implicit none
  integer :: i

  call parse_input()
  do i=istart_reac_coords,nreac_coords
    call do_mean_force(i)
  enddo

  call do_thermodynamic_integration()
  call finalize()

  stop

end program bluemoon

integer function f_getBoolValue(tag,ret)
  implicit none
  character(*),intent(in) :: tag
  logical,intent(out)     :: ret
  integer :: integ,iret,f_getIntValue
  iret = f_getIntValue(tag,integ)
  if ( integ==0 ) then
    ret = .false.
  else
    ret = .true.
  endif
  f_getBoolValue = iret
end function f_getBoolValue

