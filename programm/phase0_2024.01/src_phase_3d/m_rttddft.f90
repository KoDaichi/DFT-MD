!!
!! Last updated: 01 Feb 2011 on asahi.nims.go.jp
!! real-time electron-ion dynamics simulation
!! time-dependent density functional theory
!!
!! 12 Sep 2012, some unstable routines removed, isir, osaka univ.
!!

module m_rttddft

use m_Const_Parameters,only : &
  &  DP,PAI,PAI2,PAI4,CMPLDP,zi,ELECTRON,DIRECT,ON,OFF,BUCS,CARTS,MAPPED &
  & ,SmallestPositiveNumber
use m_Parallelization,only : &
  &  MPI_CommGroup,ista_kngp,iend_kngp,npes,mype,ierr &
  & ,myrank_e,map_e,ista_e,iend_e,istep_e,myrank_k,map_k,ista_k,iend_k,np_e,map_z
use m_Control_Parameters,only : &
  &  neg,kimg,nspin,af &
  & ,sw_rttddft &
  & ,ext_ie_elec,ext_ie_hole &
  & ,ext_pulse_epsilon,ext_pulse_kx,ext_pulse_ky,ext_pulse_kz &
  & ,time_step_max,time_step_delta &
  & ,propagator_method,propagator_order
use m_Crystal_Structure,only : &
  &  altv,rltv,univol
use m_Ionic_System,only : &
  &  ntyp, iatom, ityp, iwei, natm, natm2, ityp, pos, cps
use m_FFT,only : &
  &  m_FFT_WF,m_FFT_Vlocal_W_3D,fft_box_size_WF,nfft
use m_PseudoPotential,only : &
  &  ival,ilmt,nlmt,nlmtt,nlmta,lmta,lmtt,ltp,mtp,taup,q,dion &
  & ,nloc,ntau,ntyp,mmesh,nmesh,betar,wos,lpsmax
use m_PlaneWaveBasisSet,only : &
  &  ngabc,kg1,m_pwBS_kinetic_energies,nbase,iba,igf
use m_Kpoints,only : &
  &  kv3,vkxyz,qwgt
use m_Charge_Density,only : &
  &  m_CD_softpart_3D,m_CD_hardpart,m_CD_rspace_charge &
  & ,m_CD_alloc_rspace_charge,m_CD_dealloc_rspace_charge &
  & ,chgq_l,hsr
!!---- phase_v800
!use m_Electronic_Structure,only : &
!  &  zaj_l,vlhxcQ,occup_l,eko_l,fsr_l,fsi_l,vlhxc_l &
!  & ,m_ES_Vlocal_in_Rspace,m_ES_WF_in_Rspace &
!  & ,m_ES_alloc_fft_related,m_ES_dealloc_fft_related,afft,bfft &
!  & ,m_ES_alloc_scss_etc,m_ES_dealloc_scss_etc,sc,ss,qc,qs &
!  & ,m_ES_Vnonlocal_W,vnlph_l,m_ES_betar_dot_WFs &
!  & ,m_ES_betar_dot_WFs_4_each_k,m_ES_betar_dot_WFs_4_lmta_k
!!---- phaseUnif
use m_Electronic_Structure,only : &
  &  zaj_l,vlhxcQ,occup_l,eko_l,fsr_l,fsi_l,vlhxc_l &
  & ,m_ES_Vlocal_in_Rspace_3D,m_ES_WF_in_Rspace_3D &
  & ,m_ES_alloc_fft_related,m_ES_dealloc_fft_related,afft,bfft &
  & ,vnlph_l
use m_ES_nonlocal,only : &
  &  sc,ss,qc,qs &
  & ,m_ES_Vnonlocal_W_3D &
  & ,m_ES_alloc_scss_etc_3D &
  & ,m_ES_dealloc_scss_etc &
  & ,m_ES_betar_dot_WFs_3D &
  & ,m_ES_betar_dot_WFs_4_each_k_3D &
  & ,m_ES_betar_dot_WFs_4_lmta_k_3D
!!----
use m_NonLocal_Potential,only : &
  &  snl
use m_Files,only : &
  &  nfout,nfchr,m_Files_open_nfchr
use m_Parallelization,only : mpi_kg_world, mpi_ke_world, mpi_ge_world, &
                             nel_fft_z, nel_fft_y, nel_fft_x, &
                             mp_g1k, myrank_g, ista_g1k, iend_g1k, xyz_fft_y, np_g1k
use m_FFT,            only : m_FFT_Direct_3D
use m_ES_WF_by_SDorCG,only : map_fft_to_WF_3D
use mpi

implicit none
!include 'mpif.h'

contains

!!
!!---- external impulse electric field at time=0+
!!     vector potential: A_ext(t) = [A0*r_i]*theta(t)
!!     electric field: E_ext(t) = -(dA(t)/dt)/c = -[A0*r_i]*delta(t)/c
!!     phi(r,t=0+) = exp(i*k_ext*r_i)*phi(r,t=0-)
!!
subroutine m_rttddft_init_impulse_field

  integer :: is,ik,ib,i1,ri,i,j,k,id,mm,nl,nm,nn,n1,n2,n3,nlhf,inew,jnew,knew,ip
  real(kind=DP) :: rx,ry,rz,dn1,dn2,dn3,denom
  complex(kind=CMPLDP) :: ctmp
  integer :: iadd, nx, ny, nz
! === FFT 3D ===================================================================
  integer :: ibsize, lsize, isrsize, fft_l_size
  real(kind=DP), allocatable, dimension(:,:) :: wk_bfft_l
  real(kind=DP), allocatable, dimension(:,:) :: bfft_l
! ==============================================================================

  if(ext_pulse_epsilon==0.0.and. &
  &  ext_pulse_kx==0.0.and.ext_pulse_ky==0.0.and.ext_pulse_kz==0.0) return
! === FFT 3D ===================================================================
  ibsize = 1
  lsize = max(maxval(nel_fft_x(:)),maxval(nel_fft_y(:)),maxval(nel_fft_z(:)))
  allocate(wk_bfft_l(lsize*kimg,ibsize), stat=ierr)
  allocate(bfft_l(lsize*kimg,ibsize) ,stat=ierr)
! ==============================================================================

  write(nfout,'(''# impulse electric field at time=0+'')')
  write(nfout,'(''# eps='',f6.3,'' kx='',f6.3,'' ky='',f6.3,'' kz='',f6.3)') &
  &  ext_pulse_epsilon,ext_pulse_kx,ext_pulse_ky,ext_pulse_kz

  call m_ES_alloc_fft_related()
  denom=1.0d0/product(fft_box_size_WF(1:3,1))
  nl=fft_box_size_WF(1,1); nm=fft_box_size_WF(2,1); nn=fft_box_size_WF(3,1)
  nx = xyz_fft_y(2,1)-xyz_fft_y(1,1)+1
  ny = xyz_fft_y(2,2)-xyz_fft_y(1,2)+1
  nz = xyz_fft_y(2,3)-xyz_fft_y(1,3)+1
! Now, kimg should not be 1!!!

  do ik=1,kv3
  if(map_k(ik) /= myrank_k) cycle
! === FFT 3D ===================================================================
  isrsize = min(lsize,mp_g1k(ik))
  fft_l_size  = nel_fft_x(myrank_g)
! ==============================================================================
  do ib=1,np_e
    call m_ES_WF_in_Rspace_3D(ik,ib,ib,ibsize,lsize,wk_bfft_l)
    do k=0,nz-1
    do j=0,ny-1
    do i=0,nx-1
! Now, kimg should not be 1!!!
      knew=k; jnew=j; inew=i
      ip=jnew+inew*ny+knew*ny*nx+1
      dn1=dble(i+xyz_fft_y(1,1)-1)/dble(nl)
      dn2=dble(j+xyz_fft_y(1,2)-1)/dble(nm)
      dn3=dble(k+xyz_fft_y(1,3)-1)/dble(nn)
      rx=altv(1,1)*dn1+altv(1,2)*dn2+altv(1,3)*dn3
      ry=altv(2,1)*dn1+altv(2,2)*dn2+altv(2,3)*dn3
      rz=altv(3,1)*dn1+altv(3,2)*dn2+altv(3,3)*dn3
      ctmp=cdexp(zi*ext_pulse_epsilon* &
        & (ext_pulse_kx*rx+ext_pulse_ky*ry+ext_pulse_kz*rz)) &
        & *dcmplx(wk_bfft_l(ip*2-1,1),wk_bfft_l(ip*2,1))
      wk_bfft_l(ip*2-1,1)=dreal(ctmp)
      wk_bfft_l(ip*2,  1)=dimag(ctmp)
    enddo
    enddo
    enddo
    call m_FFT_Direct_3D(nfout,wk_bfft_l,lsize,ibsize)
    call map_fft_to_WF_3D(ik,lsize,ibsize,wk_bfft_l,bfft_l,isrsize,fft_l_size)
    do ri=1,kimg
    do i=ista_g1k(ik),iend_g1k(ik)
      iadd = i - ista_g1k(ik) + 1
      i1=kimg*iadd+(ri-kimg)
      zaj_l(iadd,ib,ik,ri)=bfft_l(i1,1)*denom
    enddo
    enddo
  enddo
  enddo
  call m_ES_dealloc_fft_related()
! === FFT 3D ===================================================================
  deallocate(wk_bfft_l)
  deallocate(bfft_l)
! ==============================================================================
  do ik=1,kv3
    if(map_k(ik) /= myrank_k) cycle
    call m_ES_betar_dot_WFs_3D(nfout,ik)
  end do

end subroutine m_rttddft_init_impulse_field


!!
!!---- electron excitation: manual occupation control at time=0+
!!     occupation(filled-band,k) <--> occupation(empty-band,k)
!!
subroutine m_rttddft_init_occup_control

  integer :: ik,ie
  real(kind=DP) :: occup_tmp1,occup_tmp2
  real(kind=DP),allocatable,dimension(:,:) :: eko_mpi,occup_mpi,tmp_mpi

  if(ext_ie_elec==0.or.ext_ie_hole==0) return

  allocate(eko_mpi(neg,kv3)) ; eko_mpi=0.0d0
  allocate(occup_mpi(neg,kv3)) ; occup_mpi=0.0d0

  if(npes >= 2) then
    allocate(tmp_mpi(neg,kv3))
!-- occup_mpi
    tmp_mpi=0.0d0
    do ik=1,kv3
      if(map_k(ik) /= myrank_k) cycle
      do ie=1,neg
        if(map_e(ie) /= myrank_e) cycle
        tmp_mpi(ie,ik)=occup_l(map_z(ie),ik)
      enddo
    enddo
    call mpi_allreduce(tmp_mpi,occup_mpi,neg*kv3 &
      & ,mpi_double_precision,mpi_sum,mpi_kg_world,ierr)
    tmp_mpi = occup_mpi
    call mpi_allreduce(tmp_mpi,occup_mpi,neg*kv3 &
      & ,mpi_double_precision,mpi_sum,mpi_ge_world,ierr)
!-- eko_mpi
    tmp_mpi=0.0d0
    do ik=1,kv3
      if(map_k(ik) /= myrank_k) cycle
      do ie=1,neg
        if(map_e(ie) /= myrank_e) cycle
        tmp_mpi(ie,ik)=eko_l(map_z(ie),ik)
      enddo
    enddo
    call mpi_allreduce(tmp_mpi,eko_mpi,neg*kv3 &
      & ,mpi_double_precision,mpi_sum,mpi_kg_world,ierr)
    tmp_mpi = eko_mpi
    call mpi_allreduce(tmp_mpi,eko_mpi,neg*kv3 &
      & ,mpi_double_precision,mpi_sum,mpi_ge_world,ierr)
    deallocate(tmp_mpi)
  else
    occup_mpi=occup_l
    eko_mpi=eko_l
  endif

  if(mype==0) then
    write(nfout,'(''# time=0-'')')
    do ik=1,kv3
      write(nfout,'(''# occ ik='',i4)') ik
      write(nfout,'(10f8.4)') (occup_mpi(ie,ik),ie=1,neg)
      write(nfout,'(''# eko ik='',i4)') ik
      write(nfout,'(10f8.4)') (eko_mpi(ie,ik),ie=1,neg)
    enddo
  endif

!-- exchange occupation
  do ik=1,kv3
    occup_tmp1=occup_mpi(ext_ie_elec,ik)
    occup_tmp2=occup_mpi(ext_ie_hole,ik)
    occup_mpi(ext_ie_elec,ik)=occup_tmp2
    occup_mpi(ext_ie_hole,ik)=occup_tmp1
  enddo

  do ik=ista_k,iend_k
    do ie=1,neg
      if(map_e(ie) /= myrank_e) cycle
      occup_l(map_z(ie),ik)=occup_mpi(ie,ik)
    enddo
  enddo

  if(mype==0) then
    write(nfout,'(''# time=0+'')')
    do ik=1,kv3
      write(nfout,'(''# occ ik='',i4)') ik
      write(nfout,'(10f8.4)') (occup_mpi(ie,ik),ie=1,neg)
      write(nfout,'(''# eko ik='',i4)') ik
      write(nfout,'(10f8.4)') (eko_mpi(ie,ik),ie=1,neg)
    enddo
  endif

  deallocate(occup_mpi,eko_mpi)

end subroutine m_rttddft_init_occup_control


!!
!!---- propagator: Taylor expansion method
!!
subroutine m_rttddft_propagate_wf_by_taylor

  integer :: is,ik,ib,ie,iadd
  integer :: i,i1,iksnl,n_taylor,n_taylor_max
  real(kind=DP) :: fac,denom
  real(kind=DP),allocatable,dimension(:) :: ekin_tmp
  complex(kind=CMPLDP) :: cfac,ctmp
  complex(kind=CMPLDP),allocatable,dimension(:,:,:) :: czaj_l
  complex(kind=CMPLDP),allocatable,dimension(:) :: czaj_tmp
! === FFT 3D ===================================================================
  integer :: ibsize, lsize, isrsize, fft_l_size
  real(kind=DP), allocatable, dimension(:,:) :: wk_bfft_l
  real(kind=DP), allocatable, dimension(:,:) :: bfft_l
  real(kind=DP), allocatable, dimension(:) :: afft_l
! ==============================================================================

  allocate(czaj_l(maxval(np_g1k),np_e,ista_k:iend_k)); czaj_l=(0.0d0,0.0d0)
  allocate(czaj_tmp(maxval(np_g1k))); czaj_tmp=(0.0d0,0.0d0)
  allocate(ekin_tmp(maxval(np_g1k))); ekin_tmp=0.0d0

  call m_ES_alloc_fft_related()
! === FFT 3D ===================================================================
  ibsize = 1
  lsize = max(maxval(nel_fft_x(:)),maxval(nel_fft_y(:)),maxval(nel_fft_z(:)))
  allocate(wk_bfft_l(lsize*kimg,ibsize), stat=ierr)
  allocate(bfft_l(lsize*kimg,ibsize) ,stat=ierr)
  allocate(afft_l(lsize*kimg), stat=ierr)
! ==============================================================================

  if(propagator_order<0) then
    stop '** ERROR(RTTDDFT): propagator_order should be more than zero'
  endif
  n_taylor_max=propagator_order  !! input parameters (positive integer or zero)
  cfac=-zi*time_step_delta
  fac=1.0
  denom=1.0d0/product(fft_box_size_WF(1:3,1))

!!-- n=0 term
  do ik=1,kv3
    if(map_k(ik)/=myrank_k) cycle
    do ib=1,np_e
      do i=ista_g1k(ik),iend_g1k(ik)
        iadd = i - ista_g1k(ik) + 1
        czaj_l(iadd,ib,ik)=dcmplx(zaj_l(iadd,ib,ik,1),(kimg-1)*zaj_l(iadd,ib,ik,kimg))
      enddo
    enddo
  enddo

!!-- n>1 terms
  do n_taylor=1,n_taylor_max  !! serial processing
    fac=fac/dble(n_taylor)  !! inverse of n factorial = 1/(n!)
    do is=1,nspin,af+1
      call m_ES_Vlocal_in_Rspace_3D(is,afft_l,lsize,1,OFF)
      do ik=is,kv3-nspin+is,nspin
        if(map_k(ik)/=myrank_k) cycle
! === FFT 3D ===================================================================
        isrsize = min(lsize,mp_g1k(ik))
        fft_l_size  = nel_fft_x(myrank_g)
! ==============================================================================
        call m_pwBS_kinetic_energies(ik,vkxyz,ekin_tmp)
        iksnl=(ik-1)/nspin+1
        call m_ES_Vnonlocal_W_3D(ik,iksnl,is,OFF)
        do ie=1,np_e
          ib = ie
          czaj_tmp=(0.0d0,0.0d0)
!! T|wf>
          do i=ista_g1k(ik),iend_g1k(ik)
            iadd = i - ista_g1k(ik) + 1
            czaj_tmp(iadd)=czaj_tmp(iadd) &
              +ekin_tmp(iadd)*dcmplx(zaj_l(iadd,ib,ik,1),(kimg-1)*zaj_l(iadd,ib,ik,kimg))
          enddo
!! V_lhxc|wf>
          call m_ES_WF_in_Rspace_3D(ik,ib,ib,ibsize,lsize,wk_bfft_l)
          call m_FFT_Vlocal_W_3D(afft_l,wk_bfft_l,lsize,ibsize,nel_fft_y(myrank_g))
          call m_FFT_Direct_3D(nfout,wk_bfft_l,lsize,ibsize)
          call map_fft_to_WF_3D(ik,lsize,ibsize,wk_bfft_l,bfft_l,isrsize,fft_l_size)
          do i=ista_g1k(ik),iend_g1k(ik)
            iadd = i - ista_g1k(ik) + 1
            czaj_tmp(iadd)=czaj_tmp(iadd)+dcmplx(bfft_l(2*iadd-1,1),bfft_l(2*iadd,1))*denom
          enddo
!! V_nonlocal|wf>
          do i=ista_g1k(ik),iend_g1k(ik)
            iadd = i - ista_g1k(ik) + 1
            czaj_tmp(iadd)=czaj_tmp(iadd) &
              +dcmplx(vnlph_l(iadd,ib,1),(kimg-1)*vnlph_l(iadd,ib,kimg))
          enddo
! czaj_l=sum_n{(1/n!)(-idtH)^n|wf>} and zaj_l=(-idtH)^n|wf>
          do i=ista_g1k(ik),iend_g1k(ik)
            iadd = i - ista_g1k(ik) + 1
            czaj_l(iadd,ib,ik)=czaj_l(iadd,ib,ik)+fac*cfac*czaj_tmp(iadd)
            zaj_l(iadd,ib,ik,1)=dreal(cfac*czaj_tmp(iadd))
            if(kimg==2) zaj_l(iadd,ib,ik,kimg)=dimag(cfac*czaj_tmp(iadd))
          enddo
        enddo
        call m_ES_betar_dot_WFs_3D(nfout,ik)
      enddo
    enddo
  enddo

!-- zaj_l(t+dt)=exp(-idtH)*zaj_l(t)
  do is=1,nspin,af+1
    do ik=is,kv3-nspin+is,nspin
      if(map_k(ik) /= myrank_k) cycle
      do ib=1,np_e
        do i=ista_g1k(ik),iend_g1k(ik)
          iadd = i - ista_g1k(ik) + 1
          zaj_l(iadd,ib,ik,1)=dreal(czaj_l(iadd,ib,ik))
          if(kimg==2) zaj_l(iadd,ib,ik,kimg)=dimag(czaj_l(iadd,ib,ik))
        enddo
      enddo
      call m_ES_betar_dot_WFs_3D(nfout,ik)
    enddo
  enddo

  call m_ES_dealloc_fft_related()
! === FFT 3D ===================================================================
  deallocate(wk_bfft_l)
  deallocate(bfft_l)
  deallocate(afft_l)
! ==============================================================================
  deallocate(ekin_tmp)
  deallocate(czaj_tmp)
  deallocate(czaj_l)
end subroutine m_rttddft_propagate_wf_by_taylor



!!
!!---- dipole moment
!!     NOTE: Here is only valid for isolated systems,
!!     and a molecule should be located at the center of unit cell.
!!
subroutine m_rttddft_dipole
  integer :: i,ia,id,im,ip
  real(kind=DP) :: g(3),r0(3),len(3)
  real(kind=DP) :: d_ele(3),d_ion(3)
  real(kind=DP) :: d_ele_mpi(3)
  complex(kind=CMPLDP) :: ccha

  d_ele=0.0d0; d_ion=0.0d0
  r0(1)=0.0d0; r0(2)=0.0d0; r0(3)=0.0d0
  len(1)=sqrt(sum(altv(1:3,1)**2)); len(2)=sqrt(sum(altv(1:3,2)**2)); len(3)=sqrt(sum(altv(1:3,3)**2))

! electron
  do i=ista_kngp,iend_kngp
    ccha=dcmplx(chgq_l(i,1,1    ),(kimg-1)*chgq_l(i,kimg,1    )) &
      & +dcmplx(chgq_l(i,1,nspin),(kimg-1)*chgq_l(i,kimg,nspin))*(nspin-1)
    g(1)=rltv(1,1)*ngabc(i,1)+rltv(1,2)*ngabc(i,2)+rltv(1,3)*ngabc(i,3)
    g(2)=rltv(2,1)*ngabc(i,1)+rltv(2,2)*ngabc(i,2)+rltv(2,3)*ngabc(i,3)
    g(3)=rltv(3,1)*ngabc(i,1)+rltv(3,2)*ngabc(i,2)+rltv(3,3)*ngabc(i,3)
    do id=1,3
      if(id==1) then; ip=2; im=3; endif
      if(id==2) then; ip=3; im=1; endif
      if(id==3) then; ip=1; im=2; endif
      if(abs(g(ip))<1.d-10.and.abs(g(im))<1.d-10) then
        if(abs(g(id))<1.d-10) then
          d_ele(id)=d_ele(id)+len(ip)*len(im)*ccha*0.5*((r0(id)+len(id))**2.0-r0(id)**2.0)
        else
          d_ele(id)=d_ele(id)+len(ip)*len(im)*ccha &
        & *((zi*g(id)*(r0(id)+len(id))-1.0)*cdexp(zi*g(id)*(r0(id)+len(id))) &
        &  -(zi*g(id)*r0(id)-1.0)*cdexp(zi*g(id)*r0(id)))/(zi*g(id))**2.0
        endif
      endif
    enddo
  enddo
  if(npes > 1) then
    d_ele_mpi=0.0d0
    call mpi_allreduce(d_ele,d_ele_mpi,3,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
    d_ele=d_ele_mpi
    call mpi_allreduce(d_ele,d_ele_mpi,3,mpi_double_precision,mpi_sum,mpi_ge_world,ierr)
    d_ele=d_ele_mpi
  endif
  d_ele=-d_ele

! ion
  do ia=1,natm
  do id=1,3
    d_ion(id)=d_ion(id)+ival(ityp(ia))*cps(ia,id)
    if(iwei(ia)/=1) &
  & d_ion(id)=d_ion(id)-ival(ityp(ia))*cps(ia,id)
  enddo
  enddo

  write(nfout,'(''P='',3f20.10)') (d_ele(id)+d_ion(id),id=1,3)
! for debug
! write(nfout,'(''# dipole moment (tot) ='',3f20.10)') (d_ele(id)+d_ion(id),id=1,3)
! write(nfout,'(''# dipole moment (ele) ='',3f20.10)') (d_ele(id),id=1,3)
! write(nfout,'(''# dipole moment (ion) ='',3f20.10)') (d_ion(id),id=1,3)

end subroutine m_rttddft_dipole


!!
!!---- current density, induced polarization, & induced vector potential
!!
subroutine m_rttddft_current
  integer :: is,ik,ib,i,nb,id,iadd
  real(kind=DP) :: g(3),fac,tmp
  complex(kind=CMPLDP) :: j_ele(3),ctmp
  complex(kind=CMPLDP) :: j_ele_mpi(3)

  real(kind=DP),save :: p_ind(3)=0.0d0
  real(kind=DP),save :: a_ind(3)=0.0d0

  j_ele=(0.0d0,0.0d0)
  do is=1,nspin,af+1
  do ik=is,kv3-nspin+is,nspin
  if(map_k(ik) /= myrank_k) cycle
  do ib=1,np_e
  fac=dble(2/nspin)*qwgt(ik)*occup_l(ib,ik)
  do i=ista_g1k(ik),iend_g1k(ik)
    iadd = i - ista_g1k(ik) + 1
    nb=nbase(i,ik)
    g(1)=vkxyz(ik,1,CARTS)+rltv(1,1)*ngabc(nb,1)+rltv(1,2)*ngabc(nb,2)+rltv(1,3)*ngabc(nb,3)
    g(2)=vkxyz(ik,2,CARTS)+rltv(2,1)*ngabc(nb,1)+rltv(2,2)*ngabc(nb,2)+rltv(2,3)*ngabc(nb,3)
    g(3)=vkxyz(ik,3,CARTS)+rltv(3,1)*ngabc(nb,1)+rltv(3,2)*ngabc(nb,2)+rltv(3,3)*ngabc(nb,3)
    ctmp=fac*dconjg(dcmplx(zaj_l(iadd,ib,ik,1),(kimg-1)*zaj_l(iadd,ib,ik,kimg))) &
                 & *dcmplx(zaj_l(iadd,ib,ik,1),(kimg-1)*zaj_l(iadd,ib,ik,kimg))
    j_ele(1)=j_ele(1)+g(1)*ctmp
    j_ele(2)=j_ele(2)+g(2)*ctmp
    j_ele(3)=j_ele(3)+g(3)*ctmp
  enddo
  enddo
  enddo
  enddo
  if(npes > 1) then
    j_ele_mpi=0.0d0
    call mpi_allreduce(j_ele,j_ele_mpi,3,mpi_double_complex,mpi_sum,MPI_CommGroup,ierr)
    j_ele=j_ele_mpi
  endif
  j_ele=-j_ele

  write(nfout,'(''J='',3f20.10)') (dreal(j_ele(id)),id=1,3)
! for debug
! write(nfout,'(''# total current (ele) ='',3f20.10)') (dreal(j_ele(id)),id=1,3)

  tmp=sqrt(dimag(j_ele(1))**2.+dimag(j_ele(2))**2.+dimag(j_ele(3))**2.)
  if(tmp>1.d5*SmallestPositiveNumber) &
  &  write(nfout,'(''# WARNING: Im[J(t)]='',3f20.10)') (dimag(j_ele(id)),id=1,3)


!!-- induced polarization: dP(t)/dt = j(t)/V
!  p_ind=p_ind+dreal(j_ele)*time_step_delta/univol
!  write(nfout,'(''# induced polarization (ele) ='',3f20.10)') (p_ind(id),id=1,3)
!
!!-- **** here is not reliable ****
!!-- induced vector potential: dA(t)/dt = P(t)
!  a_ind=a_ind+p_ind*time_step_delta
!  write(nfout,'(''# induced vector potential (ele) ='',3f20.10)') (a_ind(id),id=1,3)

end subroutine m_rttddft_current


!!
!!---- for debug write: eigen values & occupations
!!
subroutine m_rttddft_print_eko_occup

  integer :: ie,ik
  real(kind=DP),allocatable,dimension(:,:) :: eko_mpi,occup_mpi,tmp_mpi

  allocate(eko_mpi(neg,kv3)) ; eko_mpi=0.0d0
  allocate(occup_mpi(neg,kv3)) ; occup_mpi=0.0d0

  if(npes >= 2) then
    allocate(tmp_mpi(neg,kv3)) ; tmp_mpi=0.0d0
!-- occup mpi
    do ik=1,kv3
      if(map_k(ik) /= myrank_k) cycle
      do ie=1,neg
        if(map_e(ie) /= myrank_e) cycle
        tmp_mpi(ie,ik)=occup_l(map_z(ie),ik)
      enddo
    enddo
    call mpi_allreduce(tmp_mpi,occup_mpi,neg*kv3 &
      & ,mpi_double_precision,mpi_sum,mpi_kg_world,ierr)
    tmp_mpi = occup_mpi
    call mpi_allreduce(tmp_mpi,occup_mpi,neg*kv3 &
      & ,mpi_double_precision,mpi_sum,mpi_ge_world,ierr)
!-- eko mpi
    tmp_mpi = 0.0d0
    do ik=1,kv3
      if(map_k(ik) /= myrank_k) cycle
      do ie=1,neg
        if(map_e(ie) /= myrank_e) cycle
        tmp_mpi(ie,ik)=eko_l(map_z(ie),ik)
      enddo
    enddo
    call mpi_allreduce(tmp_mpi,eko_mpi,neg*kv3 &
      & ,mpi_double_precision,mpi_sum,mpi_kg_world,ierr)
    tmp_mpi = eko_mpi
    call mpi_allreduce(tmp_mpi,eko_mpi,neg*kv3 &
      & ,mpi_double_precision,mpi_sum,mpi_ge_world,ierr)
    deallocate(tmp_mpi)
  else
    occup_mpi=occup_l
    eko_mpi=eko_l
  endif

  if(mype==0) then
    do ik=1,kv3
      write(nfout,'(''# occ ik='',i4)') ik
      write(nfout,'(10f8.4)') (occup_mpi(ie,ik),ie=1,neg)
      write(nfout,'(''# eko ik='',i4)') ik
      write(nfout,'(10f8.4)') (eko_mpi(ie,ik),ie=1,neg)
    enddo
  endif

  deallocate(occup_mpi,eko_mpi)

end subroutine m_rttddft_print_eko_occup



end module m_rttddft
