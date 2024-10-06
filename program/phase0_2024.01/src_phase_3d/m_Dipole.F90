!=======================================================================
!
!  SOFTWARE NAME : PHASE (ver. 6.01)
!
!  MODULE: m_Dipole
!
!  AUTHOR(S): T. Yamamoto   December/21/2005
!
!  Contact address :  Phase System Consortium
!                     E-mail: phase_system@nims.go.jp URL https://azuma.nims.go.jp
!
!
!
!=======================================================================
!
!
module m_Dipole
  use m_Const_Parameters,    only : DP,PAI4,PAI2,PAI,UP,DOWN,ON, speed_of_light,CONST_EV,BOHR_RADIUS
  use m_Parallelization,     only : MPI_CommGroup,ista_kngp,iend_kngp,npes,mype,ierr  &
      &                           , nrank_g, nrank_chg
  use m_Control_Parameters,  only : rvac,idir_dip,ndiv_dip,width_dip,kimg &
      &                           , printable,nspin,sw_dipole_correction &
      &                           , amix_dip, elec_field
  use m_Crystal_Structure,   only : altv,univol
  use m_Ionic_System,        only : natm, natm2, ityp, iwei, pos
  use m_PseudoPotential,     only : ival
  use m_PlaneWaveBasisSet,   only : ngabc
  use m_Charge_Density,      only : chgq_l
  use m_Electronic_Structure,only : vlhxc_l

! ================================ added by K. Tagami ============== 11.0
  use m_Control_Parameters,  only : noncol, ndim_magmom
! ================================================================== 11.0

!$$#ifdef PARA3D
  use m_Parallelization,     only : mpi_chg_world
  use m_PlaneWaveBasisSet,   only : ngabc_kngp_l, ngabc_kngp_B_l
  use mpi
!$$#endif
  implicit none
!  include 'mpif.h'

  real(kind=DP), private, parameter :: E30CM = CONST_EV * BOHR_RADIUS * 1.0d30
  real(kind=DP), private, parameter :: DEBYE = speed_of_light * CONST_EV * BOHR_RADIUS * 1.0d21
  real(kind=DP), public :: dipole(3),dipole_old(3)
  real(kind=DP), public :: field(3),field_old(3)
  real(kind=DP), public :: edip_ion,eext_ion

  real(kind=DP), allocatable :: vdip_l(:,:) !d(ista_kngp:iend_kngp,kimg)
  real(kind=DP), allocatable :: vext_l(:,:) !d(ista_kngp:iend_kngp,kimg)

contains


  subroutine ion_part(a1,a2,iax1,mui)
    real(kind=DP), intent(in) :: a1,a2
    integer, intent(in) :: iax1
    real(kind=DP), intent(out) :: mui

    real(kind=DP) :: ps(natm2),t
    integer :: ia,ja,ityp2(natm2)
#ifdef __TIMER_SUB__
  call timer_sta(1053)
#endif

    ja = 0
#ifdef __TIMER_DO__
  call timer_sta(1070)
#endif
    do ia=1,natm
       ja = ja + 1
       if(iwei(ia)==1) then
          ityp2(ja) = ityp(ia)
          ps(ja) = pos(ia,iax1)
       else
          ityp2(ja) = ityp(ia)
          ps(ja) = pos(ia,iax1)
          ja = ja + 1
          ityp2(ja) = ityp(ia)
          ps(ja) = -pos(ia,iax1)
       end if
    end do
#ifdef __TIMER_DO__
  call timer_end(1070)
#endif
    mui = 0.d0
#ifdef __TIMER_DO__
  call timer_sta(1071)
#endif
    do ia=1,natm2
       t = ps(ia) - a1
       t = mod(t,1.d0)
       t = t + a1
       if(t <= a2) mui = mui + ival(ityp2(ia))*t
    end do
#ifdef __TIMER_DO__
  call timer_end(1071)
#endif
    mui = mui * sqrt(sum(altv(1:3,iax1)**2))
#ifdef __TIMER_SUB__
  call timer_end(1053)
#endif
  end subroutine ion_part

  subroutine ion_part_gauss(a1,a2,iax1,mui)
    real(kind=DP), intent(in) :: a1,a2
    integer, intent(in) :: iax1
    real(kind=DP), intent(out) :: mui

    real(kind=DP) :: ps(natm2),t, derf
    integer :: ia,ja,ityp2(natm2)
    real(kind=DP) :: a,f1,f2,w,fac1,fac2,fac3,x,length

    length = sqrt(sum(altv(1:3,iax1)**2))
    w = width_dip
    a = 4.d0*log(2.d0)
    f1 = -a/w**2
    f2 = sqrt(a)/w
    fac1 = -0.5d0*length*w/sqrt(PAI*a)
    fac2 = 0.5d0*length

    ja = 0
    do ia=1,natm
       ja = ja + 1
       if(iwei(ia)==1) then
          ityp2(ja) = ityp(ia)
          ps(ja) = pos(ia,iax1)
       else
          ityp2(ja) = ityp(ia)
          ps(ja) = pos(ia,iax1)
          ja = ja + 1
          ityp2(ja) = ityp(ia)
          ps(ja) = -pos(ia,iax1)
       end if
    end do
    mui = 0.d0
    do ia=1,natm2
       t = ps(ia) - a1
       t = mod(t,1.d0)
       t = t + a1
       x = a2-t
       mui = mui + ival(ityp2(ia))*(fac1*dexp(f1*x**2)+t*fac2*(1.d0+derf(f2*x)))
    end do
  end subroutine ion_part_gauss


  subroutine edip_ion_part(a1,a2,iax1,zx)
    real(kind=DP), intent(in) :: a1,a2
    integer, intent(in) :: iax1
    real(kind=DP), intent(out) :: zx

    real(kind=DP) :: ps(natm2),t
    integer :: ia,ja,ityp2(natm2)

#ifdef __TIMER_SUB__
  call timer_sta(1055)
#endif

    ja = 0
#ifdef __TIMER_DO__
  call timer_sta(1074)
#endif
    do ia=1,natm
       ja = ja + 1
       if(iwei(ia)==1) then
          ityp2(ja) = ityp(ia)
          ps(ja) = pos(ia,iax1)
       else
          ityp2(ja) = ityp(ia)
          ps(ja) = pos(ia,iax1)
          ja = ja + 1
          ityp2(ja) = ityp(ia)
          ps(ja) = -pos(ia,iax1)
       end if
    end do
#ifdef __TIMER_DO__
  call timer_end(1074)
#endif
    zx = 0.d0
#ifdef __TIMER_DO__
  call timer_sta(1075)
#endif
    do ia=1,natm2
       t = ps(ia) - a1
       t = mod(t,1.d0)
       t = t - 0.5d0
       zx = zx + ival(ityp2(ia))*t
    end do
#ifdef __TIMER_DO__
  call timer_end(1075)
#endif
#ifdef __TIMER_SUB__
  call timer_end(1055)
#endif
  end subroutine edip_ion_part

  subroutine m_Dipole_wd_vlhxc(nfout,vxc_l)
    integer, intent(in) :: nfout

! ============================== modified by K. Tagami ================== 11.0
!    real(kind=DP), intent(in) :: vxc_l(ista_kngp:iend_kngp,kimg,nspin)
    real(kind=DP), intent(in) :: vxc_l(ista_kngp:iend_kngp,kimg,ndim_magmom)
! ======================================================================= 11.0

    integer :: g,ia,ix,i
    integer :: iax1,iax2,iax3
    real(kind=DP) :: a1,x,ph,dx
    real(kind=DP) :: vloc(0:ndiv_dip)
    real(kind=DP), allocatable :: vloc_mpi(:)
    real(kind=DP) :: vm,vp,vj

    integer, parameter :: ispin = 1

    if(idir_dip == 0) return

    iax1 = idir_dip
    iax2 = iax1 + 1
    if(iax2 >= 4) iax2 = iax2 - 3
    iax3 = iax2 + 1
    if(iax3 >= 4) iax3 = iax3 - 3

    vloc(:) = 0.d0
    a1 = rvac(idir_dip)
    dx = 1.d0/ndiv_dip
    do ix=0,ndiv_dip
       x = a1+dx*ix
       if(kimg==1) then
          do i = ista_kngp, iend_kngp
             if(ngabc_kngp_l(i,iax2) == 0 .and. ngabc_kngp_l(i,iax3) == 0) then
                g = ngabc_kngp_l(i,idir_dip)
                vloc(ix) = vloc(ix) + (vlhxc_l(i,kimg,ispin)-vxc_l(i,kimg,ispin))*cos(PAI2*g*x)
             end if
          end do
       else
          do i = ista_kngp, iend_kngp
             if(ngabc_kngp_l(i,iax2) == 0 .and. ngabc_kngp_l(i,iax3) == 0) then
                g = ngabc_kngp_l(i,idir_dip)
                ph = PAI2*g*x
                vloc(ix) = vloc(ix) + (vlhxc_l(i,1,ispin)-vxc_l(i,1,ispin))*cos(ph) &
                         &          - (vlhxc_l(i,kimg,ispin)-vxc_l(i,kimg,ispin))*sin(ph)
             end if
          end do
       end if
    end do
    if(nrank_chg > 1) then
       allocate(vloc_mpi(0:ndiv_dip))
       vloc_mpi = 0.d0
       call mpi_allreduce(vloc,vloc_mpi,(ndiv_dip+1),mpi_double_precision,mpi_sum,mpi_chg_world,ierr)
       vloc = vloc_mpi
       deallocate(vloc_mpi)
    end if
    if(printable) then
       vm = sum(vloc(1:5))/5
       vp = sum(vloc(ndiv_dip-6:ndiv_dip-1))/5
       vj = vm-vp
       write(nfout,'("Dipole: Vh; dir_dip=",i1)') idir_dip
       write(nfout,'("Dipole: Vh(-)=",f10.5," Vh(+)=",f10.5," Vh(-)-Vh(+)",f10.5)') vm,vp,vj
       write(nfout,'(" x, Vh(x)")')
       do ix=0,ndiv_dip
          x = a1+dx*ix
          write(nfout,'(f10.5,1x,f10.5)') x,vloc(ix)
       end do
    end if
  end subroutine m_Dipole_wd_vlhxc

!$$#ifdef PARA3D

!===============================================================================

  subroutine m_Dipole_vdip_alloc_3D
    if(allocated(vdip_l)) deallocate(vdip_l)
    if(allocated(vext_l)) deallocate(vext_l)
    allocate(vdip_l(ista_kngp:iend_kngp,kimg))
    allocate(vext_l(ista_kngp:iend_kngp,kimg))
  end subroutine m_Dipole_vdip_alloc_3D

!===============================================================================
  subroutine m_Dipole_calc_3D(nfout)
    integer, intent(in) :: nfout

    real(kind=DP) :: mue,mui,mu
    real(kind=DP) :: a1,a2
    integer :: id
#ifdef __TIMER_SUB__
  call timer_sta(1051)
#endif

    do id=1,3
       if(id /= idir_dip .and. idir_dip /= 0) cycle
       a1 = rvac(id)
       a2 = a1+1.d0
       call electron_part_3D(a1,a2,id,mue)
       call ion_part(a1,a2,id,mui)
       mu = mui - mue
       if(printable) then
          write(nfout,'("Dipole: rmin,rmax,idir =",2f15.8,i3)') a1,a2,id
          write(nfout,'("Dipole: Total =",f18.8,1x,"=",f13.4,1x,"Debye =",f13.4,"x10^-30 Cm")') mu,mu*DEBYE,mu*E30CM
          write(nfout,'("Dipole: Ion   =",f20.8)') mui
          write(nfout,'("Dipole: Elec  =",f20.8)') mue
       end if
       if(sw_dipole_correction == ON) then
          dipole(id) = mu/univol*sqrt(sum(altv(1:3,id)**2))
          field(id)  = -sum(altv(1:3,id)*elec_field(1:3))/PAI4
          if(printable) then
             write(nfout,'("Dipole: mix dipole and field with amix_dip=",f15.8)') amix_dip
             write(nfout,'("Dipole: (NOW) dipole,field =",2f15.8)') dipole(id),field(id)
             write(nfout,'("Dipole: (OLD) dipole,field =",2f15.8)') dipole_old(id),field_old(id)
          end if
          dipole(id) = (1.d0-amix_dip)*dipole_old(id) + amix_dip*dipole(id)
          field(id) = (1.d0-amix_dip)*field_old(id) + amix_dip*field(id)
          if(printable) then
             write(nfout,'("Dipole: (NEW) dipole,field =",2f15.8)') dipole(id),field(id)
          end if
          dipole_old(id) = dipole(id)
          field_old(id) = field(id)
       end if
    end do
#ifdef __TIMER_SUB__
  call timer_end(1051)
#endif
  end subroutine m_Dipole_calc_3D

!===============================================================================
  subroutine electron_part_3D(a1,a2,iax1,mue)
    real(kind=DP), intent(in) :: a1,a2
    integer, intent(in) :: iax1
    real(kind=DP), intent(out) :: mue

    real(kind=DP) :: ss,chgq0,mue_mpi,mue0,length
    integer :: ista,iax2,iax3,g,i,ispin
    complex(kind=DP) :: zi,zchg,zga1,zga2

#ifdef __TIMER_SUB__
  call timer_sta(1052)
#endif
    length = sqrt(sum(altv(1:3,iax1)**2))
    zi = cmplx(0.d0,1.d0)*PAI2

    iax2 = iax1 + 1
    if(iax2 >= 4) iax2 = iax2 - 3
    iax3 = iax2 + 1
    if(iax3 >= 4) iax3 = iax3 - 3

    ss = univol*length/PAI2**2

    ista = ista_kngp
    if(ista == 1) then
       chgq0 = sum(chgq_l(1,1,1:nspin))
       mue0  = 0.5d0*chgq0*univol*(a2**2-a1**2)*length
       ista = 2
    end if
    mue = 0.d0
#ifdef __TIMER_DO__
  call timer_sta(1067)
#endif
    if(kimg == 1) then
       do ispin=1,nspin
          do i = ista, iend_kngp
             if(ngabc_kngp_l(i,iax2) == 0 .and. ngabc_kngp_l(i,iax3) == 0) then
                g = ngabc_kngp_l(i,iax1)
                zga1 = zi*g*a1
                zga2 = zi*g*a2
                mue = mue  &
                  & + chgq_l(i,1,ispin)*( ( (zga1-1.d0)*cdexp(zga1) - (zga2-1.d0)*cdexp(zga2) )/(g*g) )
             end if
          end do
       end do
    else
       do ispin=1,nspin
          do i = ista, iend_kngp
             if(ngabc_kngp_l(i,iax2) == 0 .and. ngabc_kngp_l(i,iax3) == 0) then
                g = ngabc_kngp_l(i,iax1)
                zga1 = zi*g*a1
                zga2 = zi*g*a2
                zchg = cmplx(chgq_l(i,1,ispin),chgq_l(i,2,ispin))
                mue = mue &
                  & + dble(zchg * ( ( (zga1-1.d0)*cdexp(zga1) - (zga2-1.d0)*cdexp(zga2) )/(g*g) ) )
             end if
          end do
       end do
    end if
#ifdef __TIMER_DO__
  call timer_end(1067)
#endif
    if(nrank_chg > 1) then
#ifdef __TIMER_COMM__
  call timer_barrier(mpi_chg_world)
  call timer_sta(1068)
#endif
       call mpi_bcast(mue0,1,mpi_double_precision,0,mpi_chg_world,ierr)
#ifdef __TIMER_COMM__
  call timer_end(1068)
#endif
       mue_mpi = 0.d0
#ifdef __TIMER_COMM__
  call timer_barrier(mpi_chg_world)
  call timer_sta(1069)
#endif
       call mpi_allreduce(mue, mue_mpi, 1,mpi_double_precision,mpi_sum,mpi_chg_world,ierr)
#ifdef __TIMER_COMM__
  call timer_end(1069)
#endif
       mue = mue_mpi
    end if
    mue = mue * ss + mue0
#ifdef __TIMER_SUB__
  call timer_end(1052)
#endif
  end subroutine electron_part_3D

!===============================================================================
  subroutine m_Dipole_potential_3D(nfout,chg)

    integer, intent(in) :: nfout
    real(kind=DP), intent(in) :: chg(ista_kngp:iend_kngp,kimg,nspin)

    integer :: ista,iax1,iax2,iax3,g,ik,i,ispin,id
    real(kind=DP) :: mud,muf,a1,zx,fac
#ifdef __TIMER_SUB__
  call timer_sta(1054)
#endif

    vdip_l = 0.d0
    vext_l = 0.d0
    ista = ista_kngp
    if(ista == 1) then
       ista = 2
    end if
#ifdef __TIMER_DO__
  call timer_sta(1072)
#endif

    do id=1,3
       if(id /= idir_dip) cycle
       mud = dipole(id)
       muf = field(id)
       iax1 = id
       iax2 = iax1 + 1
       if(iax2 >= 4) iax2 = iax2 - 3
       iax3 = iax2 + 1
       if(iax3 >= 4) iax3 = iax3 - 3
       a1 = rvac(id)
       do ik=1,kimg
          do i = ista, iend_kngp
             if(ngabc_kngp_l(i,iax2) == 0 .and. ngabc_kngp_l(i,iax3) == 0) then
                g = ngabc_kngp_l(i,iax1)
                if(ik==1) then
                   fac = 2.d0*sin(PAI2*g*a1)/g
                else
                   fac = -2.d0*cos(PAI2*g*a1)/g
                end if
                vdip_l(i,ik) = vdip_l(i,ik) + mud*fac
                vext_l(i,ik) = vext_l(i,ik) + muf*fac
             end if
          end do
       end do
    end do
#ifdef __TIMER_DO__
  call timer_end(1072)
#endif

    edip_ion = 0.d0
    eext_ion = 0.d0
#ifdef __TIMER_DO__
  call timer_sta(1073)
#endif
    do id=1,3
       if(id /= idir_dip) cycle
       a1 = rvac(id)
       call edip_ion_part(a1,a1+1.d0,id,zx)
       edip_ion = edip_ion + PAI4*dipole(id)*zx
       eext_ion = eext_ion + PAI4*field(id)*zx
    end do
#ifdef __TIMER_DO__
  call timer_end(1073)
#endif
    edip_ion = 0.5d0*edip_ion
    if(printable) then
       write(nfout,'("Dipole: Edip(ion),Eext(ion)=",2f15.5)') edip_ion,eext_ion
       write(nfout,'("Dipole: potential jump (dip and ext)=",2f15.5)') PAI4*dipole(idir_dip),PAI4*field(idir_dip)
    end if
#ifdef __TIMER_SUB__
  call timer_end(1054)
#endif
  end subroutine m_Dipole_potential_3D
!$$#endif

end module m_Dipole
