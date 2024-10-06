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
  &  m_FFT_WF,m_FFT_Vlocal_W,fft_box_size_WF,nfft
use m_PseudoPotential,only : &
  &  ival,ilmt,nlmt,nlmtt,nlmta,lmta,lmtt,ltp,mtp,taup,q,dion &
  & ,nloc,ntau,ntyp,mmesh,nmesh,betar,wos,lpsmax
use m_PlaneWaveBasisSet,only : &
  &  ngabc,kg1,m_pwBS_kinetic_energies,nbase,iba,igf
use m_Kpoints,only : &
  &  kv3,vkxyz,qwgt
use m_Charge_Density,only : &
  &  m_CD_softpart,m_CD_hardpart,m_CD_rspace_charge &
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
  & ,m_ES_Vlocal_in_Rspace,m_ES_WF_in_Rspace &
  & ,m_ES_alloc_fft_related,m_ES_dealloc_fft_related,afft,bfft &
  & ,vnlph_l
use m_ES_nonlocal,only : &
  &  sc,ss,qc,qs &
  & ,m_ES_Vnonlocal_W &
  & ,m_ES_alloc_scss_etc &
  & ,m_ES_dealloc_scss_etc &
  & ,m_ES_betar_dot_WFs &
  & ,m_ES_betar_dot_WFs_4_each_k &
  & ,m_ES_betar_dot_WFs_4_lmta_k
!!----
use m_NonLocal_Potential,only : &
  &  snl
use m_Files,only : &
  &  nfout,nfchr,m_Files_open_nfchr
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

  if(ext_pulse_epsilon==0.0.and. &
  &  ext_pulse_kx==0.0.and.ext_pulse_ky==0.0.and.ext_pulse_kz==0.0) return

  write(nfout,'(''# impulse electric field at time=0+'')')
  write(nfout,'(''# eps='',f6.3,'' kx='',f6.3,'' ky='',f6.3,'' kz='',f6.3)') &
  &  ext_pulse_epsilon,ext_pulse_kx,ext_pulse_ky,ext_pulse_kz

  call m_ES_alloc_fft_related()
  denom=1.0d0/product(fft_box_size_WF(1:3,1))
  id=fft_box_size_WF(1,0); mm=fft_box_size_WF(2,0)
  nl=fft_box_size_WF(1,1); nm=fft_box_size_WF(2,1); nn=fft_box_size_WF(3,1)
  if(kimg==1) then
    nlhf=id/2
  else
    nlhf=id
  endif

  do ik=1,kv3
  if(map_k(ik) /= myrank_k) cycle
  do ib=ista_e,iend_e,istep_e
    call m_ES_WF_in_Rspace(ik,ib,bfft)
    do k=1,nn
    do j=1,nm
    do i=1,nl
      if(kimg==1.and.i>nlhf) then
        inew=id-i
        jnew=nm+2-j
        knew=nn+2-k
        if(knew>nn) knew=knew-nn
        if(jnew>nm) jnew=jnew-nm
      else
        knew=k; jnew=j; inew=i
      endif
      ip=nlhf*mm*(knew-1)+nlhf*(jnew-1)+inew
      dn1=dble(i-1)/dble(nl)
      dn2=dble(j-1)/dble(nm)
      dn3=dble(k-1)/dble(nn)
      rx=altv(1,1)*dn1+altv(1,2)*dn2+altv(1,3)*dn3
      ry=altv(2,1)*dn1+altv(2,2)*dn2+altv(2,3)*dn3
      rz=altv(3,1)*dn1+altv(3,2)*dn2+altv(3,3)*dn3
      ctmp=cdexp(zi*ext_pulse_epsilon* &
        & (ext_pulse_kx*rx+ext_pulse_ky*ry+ext_pulse_kz*rz)) &
        & *dcmplx(bfft(ip*2-1),bfft(ip*2))
      bfft(ip*2-1)=dreal(ctmp)
      bfft(ip*2) =dimag(ctmp)
    enddo
    enddo
    enddo
    call m_FFT_WF(ELECTRON,nfout,bfft,DIRECT,ON)
    do ri=1,kimg
    do i=1,iba(ik)
      i1=kimg*igf(nbase(i,ik))+(ri-kimg)
      zaj_l(i,map_z(ib),ik,ri)=bfft(i1)*denom
    enddo
    enddo
  enddo
  enddo
  call m_ES_dealloc_fft_related()
  call m_ES_betar_dot_WFs(nfout)

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
      & ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
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
      & ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
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

  integer :: is,ik,ib,ie
  integer :: i,i1,iksnl,n_taylor,n_taylor_max
  real(kind=DP) :: fac,denom
  real(kind=DP),allocatable,dimension(:) :: ekin_tmp
  complex(kind=CMPLDP) :: cfac,ctmp
  complex(kind=CMPLDP),allocatable,dimension(:,:,:) :: czaj_l
  complex(kind=CMPLDP),allocatable,dimension(:) :: czaj_tmp

  allocate(czaj_l(kg1,np_e,ista_k:iend_k)); czaj_l=(0.0d0,0.0d0)
  allocate(czaj_tmp(kg1)); czaj_tmp=(0.0d0,0.0d0)
  allocate(ekin_tmp(kg1)); ekin_tmp=0.0d0

  call m_ES_alloc_fft_related()

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
      do i=1,iba(ik)
        czaj_l(i,ib,ik)=dcmplx(zaj_l(i,ib,ik,1),(kimg-1)*zaj_l(i,ib,ik,kimg))
      enddo
    enddo
  enddo

!!-- n>1 terms
  do n_taylor=1,n_taylor_max  !! serial processing
    fac=fac/dble(n_taylor)  !! inverse of n factorial = 1/(n!)
    do is=1,nspin,af+1
      call m_ES_Vlocal_in_Rspace(is,afft)
      do ik=is,kv3-nspin+is,nspin
        if(map_k(ik)/=myrank_k) cycle
        call m_pwBS_kinetic_energies(ik,vkxyz,ekin_tmp)
        iksnl=(ik-1)/nspin+1
        call m_ES_Vnonlocal_W(ik,iksnl,is,OFF)
        do ie=ista_e,iend_e,istep_e
          ib=map_z(ie)
          czaj_tmp=(0.0d0,0.0d0)
!! T|wf>
          do i=1,iba(ik)
            czaj_tmp(i)=czaj_tmp(i) &
              +ekin_tmp(i)*dcmplx(zaj_l(i,ib,ik,1),(kimg-1)*zaj_l(i,ib,ik,kimg))
          enddo
!! V_lhxc|wf>
          call m_ES_WF_in_Rspace(ik,ie,bfft)
          call m_FFT_Vlocal_W(afft,bfft)
          call m_FFT_WF(ELECTRON,nfout,bfft,DIRECT,ON)
          do i=1,iba(ik)
            i1=igf(nbase(i,ik))
            czaj_tmp(i)=czaj_tmp(i)+dcmplx(bfft(2*i1-1),bfft(2*i1))*denom
          enddo
!! V_nonlocal|wf>
          do i=1,iba(ik)
            czaj_tmp(i)=czaj_tmp(i) &
              +dcmplx(vnlph_l(i,ib,1),(kimg-1)*vnlph_l(i,ib,kimg))
          enddo
! czaj_l=sum_n{(1/n!)(-idtH)^n|wf>} and zaj_l=(-idtH)^n|wf>
          do i=1,iba(ik)
            czaj_l(i,ib,ik)=czaj_l(i,ib,ik)+fac*cfac*czaj_tmp(i)
            zaj_l(i,ib,ik,1)=dreal(cfac*czaj_tmp(i))
            if(kimg==2) zaj_l(i,ib,ik,kimg)=dimag(cfac*czaj_tmp(i))
          enddo
        enddo
      enddo
    enddo
    call m_ES_betar_dot_WFs(nfout)
  enddo

!-- zaj_l(t+dt)=exp(-idtH)*zaj_l(t)
  do is=1,nspin,af+1
    do ik=is,kv3-nspin+is,nspin
      if(map_k(ik) /= myrank_k) cycle
      do ib=1,np_e
        do i=1,iba(ik)
          zaj_l(i,ib,ik,1)=dreal(czaj_l(i,ib,ik))
          if(kimg==2) zaj_l(i,ib,ik,kimg)=dimag(czaj_l(i,ib,ik))
        enddo
      enddo
    enddo
  enddo
  call m_ES_betar_dot_WFs(nfout)

  call m_ES_dealloc_fft_related()
  deallocate(ekin_tmp)
  deallocate(czaj_tmp)
  deallocate(czaj_l)
end subroutine m_rttddft_propagate_wf_by_taylor

subroutine m_rttddft_T_wf(ik,ib,czaj_tmp)
  integer,intent(in) :: ik,ib
  integer :: i
  complex(kind=CMPLDP) :: ctmp
  real(kind=DP),allocatable,dimension(:) :: ekin_tmp
  complex(kind=CMPLDP),intent(inout),dimension(kg1) :: czaj_tmp
  allocate(ekin_tmp(kg1))
  call m_pwBS_kinetic_energies(ik,vkxyz,ekin_tmp)
  do i=1,iba(ik)
    czaj_tmp(i)=czaj_tmp(i)+ekin_tmp(i)*dcmplx(zaj_l(i,ib,ik,1),zaj_l(i,ib,ik,2))
  enddo
  deallocate(ekin_tmp)
end subroutine m_rttddft_T_wf

subroutine m_rttddft_Vlhxc_wf(is,ik,ib,czaj_tmp)
  integer,intent(in) :: is,ik,ib
  real(kind=DP) :: denom
  complex(kind=CMPLDP) :: ctmp
  complex(kind=CMPLDP),intent(inout),dimension(kg1) :: czaj_tmp
  integer :: i,i1
  denom=1.0d0/product(fft_box_size_WF(1:3,1))
  call m_ES_alloc_fft_related()
  call m_ES_Vlocal_in_Rspace(is,afft)
  call m_ES_WF_in_Rspace(ik,ib,bfft)
  call m_FFT_Vlocal_W(afft,bfft)
  call m_FFT_WF(ELECTRON,nfout,bfft,DIRECT,ON)
  do i=1,iba(ik)
    i1=igf(nbase(i,ik))
    czaj_tmp(i)=czaj_tmp(i)+dcmplx(bfft(2*i1-1),bfft(2*i1))*denom
  enddo
  call m_ES_dealloc_fft_related()
end subroutine m_rttddft_Vlhxc_wf

subroutine m_rttddft_Vnl_wf(is,ik,ib,czaj_tmp)
  integer,intent(in) :: is,ik,ib
  integer :: iksnl,i
  complex(kind=CMPLDP) :: ctmp
  complex(kind=CMPLDP),intent(inout),dimension(kg1) :: czaj_tmp
  iksnl=(ik-1)/nspin+1
  call m_ES_Vnonlocal_W(ik,iksnl,is,OFF)
  do i=1,iba(ik)
    czaj_tmp(i)=czaj_tmp(i)+dcmplx(vnlph_l(i,ib,1),vnlph_l(i,ib,2))
  enddo
end subroutine m_rttddft_Vnl_wf


!!
!!---- propagator: split operator method
!!
subroutine m_rttddft_propagate_wf_by_split(i_time)
  integer,intent(in) :: i_time
  integer :: i,n_split,n_split_max,split_order
  integer,save :: sw_tmp
  real(kind=DP) :: fac
  real(kind=DP),allocatable,dimension(:) :: pcof
  complex(kind=CMPLDP) :: cfac

  if(propagator_order/=2.and.propagator_order/=4) then
    stop '** ERROR(RT-TDDFT): invalid input propagator_order'
  endif
  split_order=propagator_order  !! input parameter (2 or 4)

  if(split_order==2) then
    n_split_max=1
    allocate(pcof(n_split_max)); pcof=0.0d0
    pcof(1)=1.0d0
  elseif(split_order==4) then
    n_split_max=5
    allocate(pcof(n_split_max)); pcof=0.0d0
    pcof(1)=1.0d0/(4.0d0-4.0d0**(1.0d0/3.0d0))
    pcof(2)=pcof(1)
    pcof(3)=1.0d0-4.0d0*pcof(1)
    pcof(4)=pcof(2)
    pcof(5)=pcof(1)
  endif

  do n_split=1,n_split_max  !! serial processing
    cfac=-0.5d0*zi*time_step_delta*pcof(n_split)
    call m_rttddft_expT_wf(cfac)
    if(i_time==1.and.n_split==1) then
      sw_tmp=1
    else
      sw_tmp=2
    endif
    call m_rttddft_expVnl_wf(cfac,sw_tmp,2)
    sw_tmp=2
    call m_rttddft_expVlhxc_wf(2.0d0*cfac)
    if(i_time==time_step_max.and.n_split==n_split_max) sw_tmp=3
    call m_rttddft_expVnl_wf(cfac,sw_tmp,1)
    call m_rttddft_expT_wf(cfac)
  enddo
  call m_ES_betar_dot_WFs(nfout)

  deallocate(pcof)
end subroutine m_rttddft_propagate_wf_by_split

subroutine m_rttddft_expT_wf(cfac)
  integer :: is,ik,ib,i
  real(kind=DP),allocatable,dimension(:) :: ekin_tmp
  complex(kind=CMPLDP),intent(in) :: cfac
  complex(kind=CMPLDP) :: ctmp
  allocate(ekin_tmp(kg1)); ekin_tmp=0.0d0

  do is=1,nspin,af+1
    do ik=is,kv3-nspin+is,nspin
      if(map_k(ik) /= myrank_k) cycle
      call m_pwBS_kinetic_energies(ik,vkxyz,ekin_tmp)
      do ib=1,np_e
        do i=1,iba(ik)
          ctmp=cdexp(cfac*ekin_tmp(i)) &
            *dcmplx(zaj_l(i,ib,ik,1),(kimg-1)*zaj_l(i,ib,ik,kimg))
          zaj_l(i,ib,ik,1)=dreal(ctmp)
          if(kimg==2) zaj_l(i,ib,ik,kimg)=dimag(ctmp)
        enddo
      enddo
    enddo
  enddo
  deallocate(ekin_tmp)

  call m_ES_betar_dot_WFs(nfout)

end subroutine m_rttddft_expT_wf

subroutine m_rttddft_expVlhxc_wf(cfac)
  integer :: is,ik,ib,i,i1,ri
  real(kind=DP) :: denom
  complex(kind=CMPLDP),intent(in) :: cfac
  complex(kind=CMPLDP) :: ctmp
  denom=1.0d0/product(fft_box_size_WF(1:3,1))
  call m_ES_alloc_fft_related()

  do is=1,nspin,af+1
    call m_ES_Vlocal_in_Rspace(is,afft)
    do ik=is,kv3-nspin+is,nspin
      if(map_k(ik)/=myrank_k) cycle
      do ib=ista_e,iend_e,istep_e
        call m_ES_WF_in_Rspace(ik,ib,bfft)

        do i=1,nfft-1,2
       !! ctmp=cdexp(cfac*dcmplx(afft(i),afft(i+1)))*dcmplx(bfft(i),bfft(i+1))
          ctmp=cdexp(cfac*afft(i))*dcmplx(bfft(i),bfft(i+1))
          bfft(i)=dreal(ctmp)
          bfft(i+1)=dimag(ctmp)
        enddo
        if(mod(nfft,2)==1) bfft(nfft)=dreal(cdexp(cfac*afft(nfft))*bfft(nfft))

        call m_FFT_WF(ELECTRON,nfout,bfft,DIRECT,ON)
        do i=1,iba(ik)
          i1=igf(nbase(i,ik))
          zaj_l(i,map_z(ib),ik,1)=bfft(2*i1-1)*denom
          if(kimg==2) zaj_l(i,map_z(ib),ik,kimg)=bfft(2*i1)*denom
        enddo

      enddo
    enddo
  enddo

  call m_ES_dealloc_fft_related()

  call m_ES_betar_dot_WFs(nfout)

end subroutine m_rttddft_expVlhxc_wf

subroutine m_rttddft_expVnl_wf(cfac,sw_tmp,as_or_des)
  integer,intent(in) :: as_or_des,sw_tmp
  integer :: it_ss,it_ee,it_ii,ia_ss,ia_ee,ia_ii,lmt_ss,lmt_ee,lmt_ii
  integer :: is,ik,ib,it,ia,i,j,nb,nr,iksnl
  integer :: lmt1,lmtt1,lmta1,il1,im1
  integer :: lmt2,lmtt2,lmta2,il2,im2
  complex(kind=CMPLDP),intent(in) :: cfac
  complex(kind=CMPLDP) :: ctmp,ccof,ccof_tmp,ccoffs
  real(kind=DP) :: ph,f1,f2,f3,brbr_tmp
  real(kind=DP),allocatable,dimension(:,:,:,:),save :: brbr
  real(kind=DP),allocatable,dimension(:,:,:),save :: zfsin_tmp,zfcos_tmp

  if(sw_tmp==1) then
    allocate(brbr(nlmtt,nlmtt,ntyp,kv3)); brbr=0.0d0
    allocate(zfsin_tmp(kg1,kv3,natm),zfcos_tmp(kg1,kv3,natm))
    zfsin_tmp=0.0d0; zfcos_tmp=0.0d0

!-- brbr=<betar_lmt2|betar_lmt1>
    if(mype==0) write(nfout,'(''# brbr=<beta_lmt2|beta_lmt1> G-space integration'')')
    do it=1,ntyp
      do lmt2=1,ilmt(it)
        il2=ltp(lmt2,it); lmtt2=lmtt(lmt2,it)
        do lmt1=1,ilmt(it)
          il1=ltp(lmt1,it); lmtt1=lmtt(lmt1,it)
          ccof=((-zi)**(-il2+1))*(zi**(-il1+1))
          do is=1,nspin,af+1
            do ik=is,kv3-nspin+is,nspin
            ! if(map_k(ik)/=myrank_k) cycle  ! comment for debug write
              iksnl=(ik-1)/nspin+1
              ctmp=(0.0d0,0.0d0)
              do i=1,iba(ik)
                ctmp=ctmp+ccof*snl(i,lmtt2,iksnl)*snl(i,lmtt1,iksnl)
              enddo
              brbr(lmtt2,lmtt1,it,iksnl)=dreal(ctmp)
              if(cdabs(ctmp) > 1.d-7.and.mype==0) &
                write(nfout,'(3i5,2f15.7)') it,lmtt1,lmtt2,dreal(ctmp),dimag(ctmp)
            enddo
          enddo
        enddo
      enddo
    enddo

!-- phase
    do ia=1,natm
      f1=pos(ia,1)*PAI2; f2=pos(ia,2)*PAI2; f3=pos(ia,3)*PAI2
      do is=1,nspin,af+1
        do ik=is,kv3-nspin+is,nspin
          if(map_k(ik)/=myrank_k) cycle
          do i=1,iba(ik)
            if(.not.allocated(nbase)) then
              stop '** ERROR(RT-TDDFT): nbase is not allocated at split-nonlocal'
            endif
            nb=nbase(i,ik)
            ph=ngabc(nb,1)*f1+ngabc(nb,2)*f2+ngabc(nb,3)*f3
            zfsin_tmp(i,ik,ia)=dsin(ph); zfcos_tmp(i,ik,ia)=dcos(ph)
          enddo
        enddo
      enddo
    enddo
  endif

  if(as_or_des==1) then
    it_ss=1; it_ee=ntyp; it_ii=+1
    ia_ss=1; ia_ee=natm; ia_ii=+1  !! ascending
  elseif(as_or_des==2) then
    it_ss=ntyp; it_ee=1; it_ii=-1
    ia_ss=natm; ia_ee=1; ia_ii=-1  !! descending
  endif

  do it=it_ss,it_ee,it_ii
    if(as_or_des==1) then
      lmt_ss=1; lmt_ee=ilmt(it); lmt_ii=+1  !! ascending
    elseif(as_or_des==2) then
      lmt_ss=ilmt(it); lmt_ee=1; lmt_ii=-1  !! descending
    endif
    do ia=ia_ss,ia_ee,ia_ii
      if(ityp(ia)/=it) cycle

      do lmt2=lmt_ss,lmt_ee,lmt_ii
        lmta2=lmta(lmt2,it); lmtt2=lmtt(lmt2,it)
        il2=ltp(lmt2,it); im2=mtp(lmt2,it)
        do lmt1=lmt_ss,lmt_ee,lmt_ii
          lmta1=lmta(lmt1,it); lmtt1=lmtt(lmt1,it)
          il1=ltp(lmt1,it); im1=mtp(lmt1,it)

          do is=1,nspin,af+1
            if(il1==il2.and.im1==im2) then
              ccof_tmp=cfac*iwei(ia)*(dion(lmt1,lmt2,it)+vlhxcQ(lmt1,lmt2,ia,is))
            else
              ccof_tmp=cfac*iwei(ia)*vlhxcQ(lmt1,lmt2,ia,is)
            endif
            do ik=is,kv3-nspin+is,nspin
              if(map_k(ik)/=myrank_k) cycle
              iksnl=(ik-1)/nspin+1
              brbr_tmp=brbr(lmtt2,lmtt1,it,iksnl)
              if(dabs(brbr_tmp) < 1.d-8) cycle
              ccof=(cdexp(ccof_tmp*brbr_tmp)-1.0d0)/brbr_tmp
              do ib=1,np_e
                ccoffs=ccof*dcmplx(fsr_l(ib,lmta2,ik),fsi_l(ib,lmta2,ik))
                do i=1,iba(ik)
                  ctmp=dcmplx(zaj_l(i,ib,ik,1),(kimg-1)*zaj_l(i,ib,ik,kimg)) &
                    +ccoffs*dcmplx(zfcos_tmp(i,ik,ia),-zfsin_tmp(i,ik,ia)) &
                    *snl(i,lmtt1,iksnl)*(zi**(-il1+1))
                  zaj_l(i,ib,ik,1)=dreal(ctmp)
                  if(kimg==2) zaj_l(i,ib,ik,kimg)=dimag(ctmp)
                enddo
              enddo
            enddo
          enddo
          call m_ES_betar_dot_WFs(nfout)

        enddo
      enddo
    enddo
  enddo

  if(sw_tmp==3) then
    deallocate(brbr)
    deallocate(zfsin_tmp,zfcos_tmp)
  endif
end subroutine m_rttddft_expVnl_wf


!!
!!---- for debug: check normality of wavefunctions & number of electrons
!!
subroutine m_rttddft_check_wf
  integer :: is,ik,ib,i,it,ia,lmt1,lmt2,lmta1,lmta2
  real(kind=DP) :: n_ele,n_ele_mpi
  complex(kind=CMPLDP) :: ctmp
  complex(kind=CMPLDP),allocatable,dimension(:,:) :: ctmp_soft,ctmp_hard
  complex(kind=CMPLDP),allocatable,dimension(:,:) :: ctmp_soft_mpi,ctmp_hard_mpi,ctmp_mpi

  if(mype==0) write(nfout,'(''# normality of Kohn-Sham wavefunction'')')

  allocate(ctmp_soft(np_e,ista_k:iend_k)) ; ctmp_soft=(0.0d0,0.0d0)
  allocate(ctmp_hard(np_e,ista_k:iend_k)) ; ctmp_hard=(0.0d0,0.0d0)
  n_ele=0.0d0

  do ik=1,kv3
    if(map_k(ik)/=myrank_k) cycle
    do ib=1,np_e
! soft part
      do i=1,iba(ik)
        ctmp_soft(ib,ik)=ctmp_soft(ib,ik) &
          +dconjg(dcmplx(zaj_l(i,ib,ik,1),(kimg-1)*zaj_l(i,ib,ik,kimg))) &
          *dcmplx(zaj_l(i,ib,ik,1),(kimg-1)*zaj_l(i,ib,ik,kimg))
      enddo
! hard part
      do it=1,ntyp
        do ia=1,natm ; if(ityp(ia)/=it) cycle
          do lmt1=1,ilmt(it) ; lmta1=lmta(lmt1,it)
            do lmt2=1,ilmt(it) ; lmta2=lmta(lmt2,it)
              ctmp_hard(ib,ik)=ctmp_hard(ib,ik) &
              +q(lmt1,lmt2,it)*iwei(ia) &
              *dconjg(dcmplx(fsr_l(ib,lmta1,ik),fsi_l(ib,lmta1,ik))) &
              *dcmplx(fsr_l(ib,lmta2,ik),fsi_l(ib,lmta2,ik))
            enddo
          enddo
        enddo
      enddo
      n_ele=n_ele+dble(2/nspin)*qwgt(ik)*occup_l(ib,ik) &
        *dreal(ctmp_soft(ib,ik)+ctmp_hard(ib,ik))
    enddo
  enddo

  allocate(ctmp_soft_mpi(neg,kv3)) ; ctmp_soft_mpi=(0.0d0,0.0d0)
  allocate(ctmp_hard_mpi(neg,kv3)) ; ctmp_hard_mpi=(0.0d0,0.0d0)

  if(npes >= 2) then
    allocate(ctmp_mpi(neg,kv3))
! soft
    ctmp_mpi=(0.0d0,0.0d0)
    do ik=1,kv3
      if(map_k(ik)/=myrank_k) cycle
      do ib=1,neg
        if(map_e(ib)/=myrank_e) cycle
        ctmp_mpi(ib,ik)=ctmp_soft(map_z(ib),ik)
      enddo
    enddo
    call mpi_allreduce(ctmp_mpi,ctmp_soft_mpi,neg*kv3 &
      ,mpi_double_complex,mpi_sum,MPI_CommGroup,ierr)
! hard
    ctmp_mpi=(0.0d0,0.0d0)
    do ik=1,kv3
      if(map_k(ik)/=myrank_k) cycle
      do ib=1,neg
        if(map_e(ib)/=myrank_e) cycle
        ctmp_mpi(ib,ik)=ctmp_hard(map_z(ib),ik)
      enddo
    enddo
    call mpi_allreduce(ctmp_mpi,ctmp_hard_mpi,neg*kv3 &
      ,mpi_double_complex,mpi_sum,MPI_CommGroup,ierr)
    deallocate(ctmp_mpi)
! number of electron
    call mpi_allreduce(n_ele,n_ele_mpi,1,mpi_double_precision &
      ,mpi_sum,MPI_CommGroup,ierr)
  else
    ctmp_soft_mpi=ctmp_soft
    ctmp_hard_mpi=ctmp_hard
    n_ele_mpi=n_ele
  endif

  deallocate(ctmp_soft,ctmp_hard)

  if(mype == 0) then
    do is=1,nspin,af+1
      do ik=is,kv3-nspin+is,nspin
        write(nfout,'(''# is='',i5,3x,''ik='',i5)') is,ik
        do ib=1,neg
! orthonormality of wavefunction: ctmp=<wf_ib|wf_ib>
          ctmp=ctmp_soft_mpi(ib,ik)+ctmp_hard_mpi(ib,ik)
          write(nfout,'(''# ib='',i5,'' tot='',2f10.5'' soft='',2f10.5,'' hard='',2f10.5)') &
            ib,ctmp,ctmp_soft_mpi(ib,ik),ctmp_hard_mpi(ib,ik)
        enddo
      enddo
    enddo
! number of electron
    write(nfout,'(''# number of electron  ='',f20.10)') n_ele_mpi
  endif

  deallocate(ctmp_soft_mpi,ctmp_hard_mpi)

end subroutine m_rttddft_check_wf


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
    call mpi_allreduce(d_ele,d_ele_mpi,3,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
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
  integer :: is,ik,ib,i,nb,id
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
  do i=1,iba(ik)
    nb=nbase(i,ik)
    g(1)=vkxyz(ik,1,CARTS)+rltv(1,1)*ngabc(nb,1)+rltv(1,2)*ngabc(nb,2)+rltv(1,3)*ngabc(nb,3)
    g(2)=vkxyz(ik,2,CARTS)+rltv(2,1)*ngabc(nb,1)+rltv(2,2)*ngabc(nb,2)+rltv(2,3)*ngabc(nb,3)
    g(3)=vkxyz(ik,3,CARTS)+rltv(3,1)*ngabc(nb,1)+rltv(3,2)*ngabc(nb,2)+rltv(3,3)*ngabc(nb,3)
    ctmp=fac*dconjg(dcmplx(zaj_l(i,ib,ik,1),(kimg-1)*zaj_l(i,ib,ik,kimg))) &
                 & *dcmplx(zaj_l(i,ib,ik,1),(kimg-1)*zaj_l(i,ib,ik,kimg))
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
      & ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
!-- eko mpi
    do ik=1,kv3
      if(map_k(ik) /= myrank_k) cycle
      do ie=1,neg
        if(map_e(ie) /= myrank_e) cycle
        tmp_mpi(ie,ik)=eko_l(map_z(ie),ik)
      enddo
    enddo
    call mpi_allreduce(tmp_mpi,eko_mpi,neg*kv3 &
      & ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
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


!!
!!---- charge density plot
!!
subroutine m_rttddft_wd_charge
  integer :: iloop
  call m_CD_alloc_rspace_charge()
  do iloop=1,nspin
    call m_Files_open_nfchr(nspin,iloop)
    call m_CD_rspace_charge(nspin,iloop,nfchr,nfout)
  enddo
  call m_CD_dealloc_rspace_charge()
end subroutine m_rttddft_wd_charge

end module m_rttddft
