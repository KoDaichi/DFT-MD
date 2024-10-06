
module m_ES_RSB
  use m_Files,                only : nfout
  use m_Const_Parameters,     only : DP,DIRECT,OFF,ELECTRON,EXECUT,NORMCONSERVATION,DELTA,ON,INVERSE
  use m_Kpoints,              only : kv3
  use m_FFT,                  only : nfft,fft_box_size_WF,m_FFT_WF
  use m_Parallelization,      only : myrank_k,map_k,np_e,ista_e,iend_e &
       &                           , istep_e, mpi_k_world, ista_k, iend_k, map_z &
       &                           , ista_ffth,iend_ffth,np_ffth,MPI_CommGroup,map_e,myrank_e,npes &
       &                           , nrank_e
  use m_Control_Parameters,   only : nspin,af,neg,printable,kimg,bisect_by,eps_rsb,lmax_rsb,ipri,iprirsb &
       &                           , sw_valence_electrons_only
  use m_PlaneWaveBasisSet,    only : kg1, iba, igf, nbase
  use m_Electronic_Structure, only : zaj_l,m_ES_WF_in_Rspace,fsr_l,fsi_l,occup_l,neordr,totch &
       &                           , m_ES_energy_eigen_values,m_ES_sort_eigen_values,eko_l
  use m_Timing,               only : tstatc0_begin,tstatc0_end
  use m_Crystal_Structure,    only : univol
  use m_ES_ortho,             only : m_ES_W_transpose_fft_cmplx, m_ES_W_transpose_back_fft_cmplx &
       &                           , m_ESortho_set_np_fs_x,np_fs_x, m_ES_F_transpose_r &
       &                           , m_ES_modified_gram_schmidt,np_g1k_x, m_ESortho_set_np_g1k_x &
       &                           , m_ES_W_transpose_r, m_ES_W_transpose_back_r,np_g1k
  use m_Realspace,            only : m_RS_get_coords
  use m_PseudoPotential,      only : modnrm,nlmta1_p,nlmta2_p,fqwei_p,nac_p,nlmta
  use m_ES_occup,             only : m_ESoc_fermi_parabolic
  use m_ES_nonlocal,          only : m_ES_betar_dot_WFs_4_each_k
  use mpi

  implicit none

  complex(kind=DP), private :: zero = (0.d0,0.d0)
  complex(kind=DP), private :: one  = (1.d0,0.d0)
  complex(kind=DP), private :: onei = (0.d0,1.d0)

  real(kind=DP), allocatable, dimension(:) :: volume
  real(kind=DP), allocatable, dimension(:,:,:) :: range_of_a_state
  real(kind=DP), allocatable, dimension(:,:,:,:) :: localized_region
  logical, allocatable, dimension(:,:) :: overlap
  integer, parameter :: nmaxregion=4
  integer :: nval,currnval

  integer :: noverl_old=0
  logical :: noverl_changed = .true.

  complex(kind=DP), allocatable, dimension(:,:) :: qout ! the unitary matrix which performs the RSB transformation

!  include 'mpif.h'                                      ! MPI

  contains

  function get_localization_volume(ndim,locrange) result(ret)
     integer, intent(in) :: ndim
     real(kind=DP), intent(in), dimension(ndim,3,nmaxregion,2) :: locrange
     real(kind=DP) :: ret
     real(kind=DP) :: dx,dy,dz,vol,v
     integer :: idi,ireg,ixyz
     vol = 0.d0
     do idi=1,ndim
        dx=0.d0;dy=0.d0;dz=0.d0
        do ireg=1,nmaxregion
           if(locrange(idi,1,ireg,1)>-0.1)then
              dx = dx + (locrange(idi,1,ireg,2)-locrange(idi,1,ireg,1))
           endif
           if(locrange(idi,2,ireg,1)>-0.1)then
              dy = dy + (locrange(idi,2,ireg,2)-locrange(idi,2,ireg,1))
           endif
           if(locrange(idi,3,ireg,1)>-0.1)then
              dz = dz + (locrange(idi,3,ireg,2)-locrange(idi,3,ireg,1))
           endif
        enddo
        v = dx*dy*dz
        if(printable.and.iprirsb>=3) write(nfout,'(a,i8,a,f10.5)') '!** volume of state ',idi,' : ',v
        vol = vol+v
     enddo
     ret = vol/real(ndim)
  end function get_localization_volume

  subroutine resolve_localization_region(ndim,nrecur,eigval,locrange)
     integer, intent(in) :: ndim,nrecur
     real(kind=DP), intent(in), dimension(ndim,nrecur) :: eigval
     real(kind=DP), intent(out), dimension(ndim,3,nmaxregion,2) :: locrange
     integer :: ir,irp,idi,ixyz,nr,iregion,iregionp,nregion,nregionp,icount
     real(kind=DP) :: c2,dr,rmin,rmax,rminp,rmaxp
     real(kind=DP),allocatable,dimension(:,:,:) :: region
     logical, allocatable, dimension(:,:) :: localized
     character(len=1),dimension(3) :: axi
     logical :: loc
     axi(1) = 'a'
     axi(2) = 'b'
     axi(3) = 'c'
     do idi=1,ndim
        do ixyz=1,3
           if(iprirsb>=3)write(nfout,'(a,i5,a)') '!-- resolving localization of state ',idi,', '//axi(ixyz)//'-axis'
           allocate(region(2**(nrecur-1)+1,2,nrecur));region=0.d0
           allocate(localized(2**(nrecur-1)+1,nrecur));localized=.false.
           do ir=1,nrecur
              nr = 2**ir
              nregion=2**(ir-1)+1
              dr = 1.d0/dble(nr)
              c2 = eigval(idi,(ir-1)*3+ixyz)
              do iregion=1,nregion
                 if(iregion.eq.1)then
                    region(iregion,1,ir) = 0.d0
                    region(iregion,2,ir) = dr
                 else if (iregion.eq.nregion)then
                    region(iregion,1,ir) = 1.d0-dr
                    region(iregion,2,ir) = 1.d0
                 else
                    region(iregion,1,ir) = (iregion-1)*2*dr-dr
                    region(iregion,2,ir) = iregion*2*dr-dr
                 endif
              enddo
              if(c2.lt.eps_rsb)then
                 do iregion=1,nregion,2
                    localized(iregion,ir) = .true.
                 enddo
              else if (1-c2.lt.eps_rsb)then
                 do iregion=2,nregion,2
                    localized(iregion,ir) = .true.
                 enddo
              endif
           enddo
           do ir=2,nrecur
              do irp=1,ir-1
                 nregion=2**(ir-1)+1;nregionp=2**(irp-1)+1
                 do iregionp=1,nregionp
                    if(localized(iregionp,irp))then
                       rminp = region(iregionp,1,irp)
                       rmaxp = region(iregionp,2,irp)
                       do iregion=1,nregion
                          rmin = region(iregion,1,ir)
                          rmax = region(iregion,2,ir)
                          if (rmaxp.lt.rmin .or. rminp.gt.rmax) then
                             localized(iregion,ir) = .false.
                          endif
                          if (rminp.gt.rmin.and.localized(iregion,ir))then
                             region(iregion,1,ir) = rminp
                          endif
                          if (rmaxp.lt.rmax.and.localized(iregion,ir))then
                             region(iregion,2,ir) = rmaxp
                          endif
                       enddo
                    endif
                 enddo
              enddo
           enddo
           do ir=1,nrecur
              c2 = eigval(idi,(ir-1)*3+ixyz)
              if(iprirsb>=3)write(nfout,'(a,i5)') '!-- level of recursion : ',ir
              nregion=2**(ir-1)+1
              do iregion=1,nregion
                 if(iprirsb>=3) write(nfout,'(a,i5,a,2f10.5,l,f10.5)') '!-- region ',iregion,' min, max, loc, cos^2 ', &
                   &  region(iregion,1,ir),region(iregion,2,ir),localized(iregion,ir),c2
              enddo
           enddo
           loc = .false.
           recur:do ir=nrecur,1,-1
              c2 = eigval(idi,(ir-1)*3+ixyz)
              nregion=2**(ir-1)+1
              icount = 0
              do iregion=1,nregion
                 if(localized(iregion,ir))then
                    if(iprirsb>=3) &
                    & write(nfout,'(a,2f10.5)') '!-- localized in region : ',region(iregion,1,ir),region(iregion,2,ir)
                    loc = .true.
                    icount = icount+1
                    locrange(idi,ixyz,icount,1) = region(iregion,1,ir)
                    locrange(idi,ixyz,icount,2) = region(iregion,2,ir)
                 endif
              enddo
              if(loc) exit recur
           enddo recur
           if(.not.loc.and.iprirsb>=3) write(nfout,'(a)') '!-- the state is delocalized'
           deallocate(region)
           deallocate(localized)
        enddo
     enddo
  end subroutine resolve_localization_region

  integer function get_nval()
     integer :: ibm,ik,ib,ib1,nv,irev,ierr
     ibm = 0
     do ik=1,kv3,af+1
        if(map_k(ik) /= myrank_k) cycle
        do ib=1,neg
           ib1 = neordr(ib,ik)
           if(map_e(ib1) == myrank_e) then
              if(occup_l(map_z(ib1),ik) > DELTA) then
                 ibm = max(ib,ibm)
              end if
           end if
        end do
     end do
     nv=0;irev=0
     if(npes>1) then
        call mpi_allreduce(ibm,irev,1,MPI_INTEGER,MPI_MAX,mpi_kg_world,ierr)
        nv= irev
        call mpi_allreduce(ibm,irev,1,MPI_INTEGER,MPI_MAX,mpi_ge_world,ierr)
        nv= irev
     else
        nv= ibm
     end if
     !!$call m_Parallel_init_mpi_nn(nfout,iprirsb,printable,nval*(nval-1)/2)
     get_nval = nv
     if(get_nval.eq.0) get_nval=totch*0.5d0

  end function get_nval

  subroutine resolve_overlap(ndim1,ndim,locregion,overl)
     integer, intent(in) :: ndim1,ndim
     real(kind=DP), intent(in), dimension(ndim,3,nmaxregion,2) :: locregion
     logical, intent(out), dimension(ndim1,ndim) :: overl
     integer :: i,j,iri,irj,ixyz,nover,ierr
     logical :: overlapping
     real(kind=DP) :: mini,maxi,minj,maxj
     logical, allocatable, dimension(:,:) :: overlt
     overl(:,:)=.true.
     nover=0
     do i=1,ndim
        if(map_e(i)/=myrank_e)cycle
        do j=1,ndim
           if (i.eq.j) cycle
           overl(map_z(i),j) = .true.
           do ixyz=1,3
              overlapping=.false.
              do iri=1,nmaxregion
                 if(locregion(i,ixyz,iri,1)>=-0.1)then
                    mini=locregion(i,ixyz,iri,1)
                    maxi=locregion(i,ixyz,iri,2)
                 else
                    cycle
                 endif
                 do irj=1,nmaxregion
                    if(locregion(j,ixyz,irj,1)>=-0.1)then
                       minj=locregion(j,ixyz,irj,1)
                       maxj=locregion(j,ixyz,irj,2)
                       if(.not.(maxi.le.minj.or.maxj.le.mini))then
                          overlapping=.true.
                       endif
                    else
                       cycle
                    endif
                 enddo
              enddo
              if(.not.overlapping)then
                 overl(map_z(i),j) = .false.
                 exit
              endif
           enddo
           if(.not.overl(map_z(i),j)) nover=nover+1
        enddo
     enddo
     call mpi_allreduce(MPI_IN_PLACE,nover,1,mpi_integer, MPI_SUM,mpi_k_world(myrank_k),ierr)
     if(printable) write(nfout,'(a,f10.5)') '!-- ratio of non-overlapping wf pairs    : ',dble(nover)/dble(ndim*ndim)

     noverl_changed = .false.
     if(noverl_old/=nover) then
        if(printable) write(nfout,'(a)') '!-- ratio of non-overlapping wf pairs has changed'
        noverl_changed = .true.
     endif
     noverl_old = nover

     allocate(overlt(ndim,ndim));overlt=.false.
     do i=1,ndim
        if(map_e(i)/=myrank_e)cycle
        overlt(i,:) = overl(map_z(i),:)
     enddo
     call mpi_allreduce(MPI_IN_PLACE,overlt,ndim*ndim,mpi_logical,mpi_lor,mpi_k_world(myrank_k),ierr)
     if(printable.and.iprirsb>=2)then
        do i=1,ndim
           do j=i+1,ndim
              write(nfout,'(a,2i8,a,l)') '!states ',i,j,' are overlapping : ',overlt(i,j)
           enddo
        enddo
     endif
     deallocate(overlt)
  end subroutine resolve_overlap

  subroutine m_ES_RSB_doit(verbose,fermi)
     logical, intent(in) :: verbose
     logical, intent(in), optional :: fermi
     integer :: ik,ib,ig,i1,ib1,ib2,i,j,k

     complex(kind=DP), allocatable, dimension(:,:) :: q,q1,qtran
     complex(kind=DP), allocatable, dimension(:,:,:) :: qoutt
     complex(kind=DP), allocatable, dimension(:) :: qtmp

     complex(kind=DP) :: ctmp
     real(kind=DP), allocatable, dimension(:) :: bfft
     real(kind=DP), allocatable, dimension(:,:) :: eigval
     real(kind=DP),allocatable,dimension(:,:,:) :: wf_t
     complex(kind=DP),allocatable,dimension(:,:) :: cwf_t,cwf_t2

     real(kind=DP) :: normr,normi,norm
     integer :: nffth
     integer :: modnrmt
     integer :: ir
     logical :: ferm
     integer :: id_sname=-1

     ferm = .true.
     if(present(fermi)) ferm = fermi

     call tstatc0_begin('m_ES_RSB ',id_sname,1)
     nffth = nfft/2
!     if(sw_valence_electrons_only==ON)then
!        nval = get_nval()
!     else
!        nval = neg
!     endif

!     if(ferm) call m_ESoc_fermi_parabolic(nfout)

     if(.not.allocated(qout)) then
        allocate(qout(nval,nval));qout=zero
        do i=1,nval
           qout(i,i) = one
        enddo
     else
!        do ik=1,kv3,af+1
!           if(map_k(ik)/=myrank_k) cycle
!           call m_ES_RSB_unitary_transform_wf(DIRECT,ik,zaj_l,.false.)
!        enddo
     endif

     allocate(qoutt(nval,nval,lmax_rsb*3));qoutt=zero
     allocate(qtmp(nffth));qtmp=zero
     allocate(bfft(nfft));bfft=0.d0
     allocate(eigval(nval,lmax_rsb*3));eigval=0.d0

     if(modnrm==EXECUT)then
        modnrmt = modnrm
        modnrm = NORMCONSERVATION
        call m_ES_modified_gram_schmidt(nfout)
        modnrm = modnrmt
     endif
     currnval = get_nval()
     do ik=1,kv3,af+1
        if(map_k(ik)/=myrank_k) cycle

        ! WF in rspace
        allocate(q(nffth,np_e));q=zero
        do ib=ista_e,iend_e,istep_e
           call m_ES_WF_in_Rspace(ik,ib,bfft)
           normr=0.d0;normi=0.d0;qtmp=zero
           do ig=1,nffth
              i1 = (ig-1)*kimg+1
              qtmp(ig) = dcmplx(bfft(i1),bfft(i1+1))
              normr = normr+bfft(i1)*bfft(i1)
              if(kimg==2) normi = normi+bfft(i1+1)*bfft(i1+1)
           enddo
           norm = normr+normi
           norm = 1.d0/dsqrt(norm)
           qtmp(1:nffth) = norm*qtmp(1:nffth)
           q(1:nffth,map_z(ib)) = qtmp(1:nffth)
        enddo
        allocate(qtran(np_ffth,neg));qtran=zero
        call m_ES_W_transpose_fft_cmplx(nffth,np_ffth,q,qtran)
        deallocate(q)

        ! recursive subspace bisection
        do ir=1,lmax_rsb
           call bisection_by_walsh(np_ffth,nval,0,qtran,qout,ir)
           qoutt(:,:,(ir-1)*3+1) = qout(:,:)
           call bisection_by_walsh(np_ffth,nval,1,qtran,qout,ir)
           qoutt(:,:,(ir-1)*3+2) = qout(:,:)
           call bisection_by_walsh(np_ffth,nval,2,qtran,qout,ir)
           qoutt(:,:,(ir-1)*3+3) = qout(:,:)
        enddo

        ! diagonalize the matrix qout
        call simultaneous_diagonalization(lmax_rsb*3,nval,qoutt,eigval,qout,verbose)

        localized_region=-1
        localized_region(:,:,1,1) = 0.d0
        localized_region(:,:,1,2) = 1.d0
        call resolve_localization_region(nval,lmax_rsb,eigval,localized_region)
        if(printable.and.verbose.or.iprirsb>=2)then
           do ib=1,nval
              write(nfout,'(a,i8)') '!-- localization info of state ',ib
              do ir=1,nmaxregion
                 if(localized_region(ib,1,ir,1)>=-0.1)then
                    write(nfout,'(a,f10.5,a,f10.5)') 'a-axis : ',localized_region(ib,1,ir,1),' - ',localized_region(ib,1,ir,2)
                 endif
              enddo
              do ir=1,nmaxregion
                 if(localized_region(ib,2,ir,1)>=-0.1)then
                    write(nfout,'(a,f10.5,a,f10.5)') 'b-axis : ',localized_region(ib,2,ir,1),' - ',localized_region(ib,2,ir,2)
                 endif
              enddo
              do ir=1,nmaxregion
                 if(localized_region(ib,3,ir,1)>=-0.1)then
                    write(nfout,'(a,f10.5,a,f10.5)') 'c-axis : ',localized_region(ib,3,ir,1),' - ',localized_region(ib,3,ir,2)
                 endif
              enddo
           enddo
        endif
        if(printable) write(nfout,'(a,f10.5)') '!-- compression ratio after localization : ',&
       & get_localization_volume(nval,localized_region)
        call resolve_overlap(np_e,nval,localized_region,overlap)

        ! unitary transformation
        !allocate(q1(np_ffth,neg));q1=zero
        !do i=1,nval
        !   do k=1,nval
        !      do j=1,np_ffth
        !         q1(j,i) = q1(j,i)+qtran(j,k)*(qout(k,i))
        !      enddo
        !   enddo
        !enddo
        !call zgemm('N','N',np_ffth,nval,nval,one,qtran,np_ffth,qout,nval,zero,q1,np_ffth)
        !if(neg>nval) q1(1:np_ffth,nval+1:neg) = qtran(1:np_ffth,nval+1:neg)
        !deallocate(qtran)
        deallocate(qtran)

        call m_ES_RSB_unitary_transform_wf(DIRECT,ik,zaj_l,.true.)
!        call m_ESortho_set_np_g1k_x(ik) !-> np_g1k_x
!        allocate(wf_t(np_g1k_x,neg,kimg));wf_t=0.d0
!        allocate(cwf_t(np_g1k_x,neg));cwf_t=0.d0
!        allocate(cwf_t2(np_g1k_x,neg));cwf_t2=0.d0
!        call m_ES_W_transpose_r(.false.,ista_k,iend_k,ik,zaj_l,wf_t)
!        if(kimg==1)then
!          cwf_t(:,:) = dcmplx(wf_t(:,:,1),0.d0)
!        else
!          cwf_t(:,:) = dcmplx(wf_t(:,:,1),wf_t(:,:,2))
!        endif
!        call zgemm('N','N',np_g1k_x,nval,nval,one,cwf_t,np_g1k_x,qout,nval,zero,cwf_t2,np_g1k_x)
!        wf_t(:,:,1) = dble(cwf_t2(:,:))
!        if(kimg>1)wf_t(:,:,2) = dimag(cwf_t2(:,:))
!        call m_ES_W_transpose_back_r(.false.,ista_k,iend_k,ik,zaj_l,wf_t)
!        deallocate(wf_t)
!        deallocate(cwf_t,cwf_t2)
!        ! back to G-space
        !allocate(q(nffth,np_e));q=0.d0
        !call m_ES_W_transpose_back_fft_cmplx(nffth,np_ffth,q,q1)
        !deallocate(q1)
        !do ib=1,np_e
        !   do ig=1,nffth
        !      i1 = (ig-1)*kimg+1
        !      bfft(i1) = real(q(ig,ib))
        !      if(kimg==2) bfft(i1+1) = aimag(q(ig,ib))
        !   enddo
        !   call m_FFT_WF(ELECTRON,nfout,bfft,DIRECT,OFF)
        !   do ig=1,iba(ik)
        !      i1 = kimg*igf(nbase(ig,ik))-1
        !      zaj_l(ig,ib,ik,1) = bfft(i1)
        !      if(kimg==2) zaj_l(ig,ib,ik,kimg) = bfft(i1+1)
        !   enddo
        !enddo
        !call m_ES_betar_dot_WFs_4_each_k(nfout,ik)   ! -> fsr_l,fsi_l
     enddo

     if(ferm)then
        !call m_ES_modified_gram_schmidt(nfout)
        !call ChargeDensity_Construction(1)
        !call Renewal_of_Potential()
        !call m_ES_energy_eigen_values(nfout,.true.)
        !call m_ESoc_fermi_parabolic(nfout)
     endif

     deallocate(qtmp)
     !!$deallocate(qout)
     deallocate(qoutt)
     deallocate(bfft)
     deallocate(eigval)
     call tstatc0_end(id_sname)

  end subroutine m_ES_RSB_doit

  real(kind=8) function get_sum_ek()
     integer :: i,ierr
     real(kind=8) :: sumek
     sumek=0.d0
     do i=1,np_e
        !if(occup_l(map_z(i),1) > DELTA) sumek = sumek+eko_l(i,1)
        sumek = sumek+eko_l(i,1)
     enddo
     call mpi_allreduce(MPI_IN_PLACE,sumek,1,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
     get_sum_ek = sumek
  end function get_sum_ek

  subroutine m_ES_RSB_alloc()
     if(sw_valence_electrons_only==ON)then
        nval = get_nval()
     else
        nval = neg
     endif
     if(.not.allocated(localized_region)) allocate(localized_region(nval,3,nmaxregion,2))
     localized_region=-1
     localized_region(:,:,1,1) = 0.d0
     localized_region(:,:,1,2) = 1.d0
     if(.not.allocated(overlap)) allocate(overlap(np_e,nval));overlap=.true.
  end subroutine m_ES_RSB_alloc

  subroutine m_ES_RSB_dealloc()
    deallocate(localized_region)
    deallocate(overlap)
  end subroutine m_ES_RSB_dealloc

  subroutine m_ES_RSB_test()
     integer :: ik,ib,ig,i1,ib1,ib2

     complex(kind=DP), allocatable, dimension(:,:) :: q,q1,qtran,qout
     complex(kind=DP), allocatable, dimension(:,:,:) :: qoutt
     complex(kind=DP), allocatable, dimension(:) :: qtmp

     real(kind=DP), allocatable, dimension(:) :: bfft
     real(kind=DP), allocatable, dimension(:) :: eigval
     real(kind=DP) :: normr,normi,norm
     integer :: nffth
     integer :: modnrmt
     integer, dimension(2) :: regionl,regionm,regionn

     integer :: id_sname=-1

     call tstatc0_begin('m_ES_RSB_test ',id_sname,1)
     nffth = nfft/2

     allocate(qout(neg,neg));qout=zero
     allocate(qoutt(neg,neg,1));qoutt=zero
     allocate(qtmp(nffth));qtmp=zero
     allocate(bfft(nfft));bfft=0.d0
     allocate(eigval(neg));eigval=0.d0

     if(modnrm==EXECUT)then
        modnrmt = modnrm
        modnrm = NORMCONSERVATION
        call m_ES_modified_gram_schmidt(nfout)
        modnrm = modnrmt
     endif

     regionl(1) = 1
     regionl(2) = fft_box_size_WF(1,1)
     regionm(1) = 1
     regionm(2) = fft_box_size_WF(2,1)
     regionn(1) = 1
     regionn(2) = fft_box_size_WF(3,1)

     do ik=1,kv3,af+1
        if(map_k(ik)/=myrank_k) cycle

        ! WF in rspace
        allocate(q(nffth,np_e));q=zero
        do ib=ista_e,iend_e,istep_e
           call m_ES_WF_in_Rspace(ik,ib,bfft)
           normr=0.d0;normi=0.d0;qtmp=zero
           do ig=1,nffth
              i1 = (ig-1)*kimg+1
              qtmp(ig) = dcmplx(bfft(i1),bfft(i1+1))
              normr = normr+bfft(i1)*bfft(i1)
              if(kimg==2) normi = normi+bfft(i1+1)*bfft(i1+1)
           enddo
           norm = normr+normi
           norm = 1.d0/dsqrt(norm)
           qtmp(1:nffth) = norm*qtmp(1:nffth)
           q(1:nffth,map_z(ib)) = qtmp(1:nffth)
        enddo
        allocate(qtran(np_ffth,neg));qtran=zero
        call m_ES_W_transpose_fft_cmplx(nffth,np_ffth,q,qtran)
        deallocate(q)

        ! bisection
        call bisection(regionl,regionm,regionn,np_ffth,bisect_by,qtran,qout,0)

        ! diagonalize the matrix qout
        !!call zheev_driver(neg,eigval,qout)
        qoutt(:,:,1) = qout(:,:)
        call simultaneous_diagonalization(1,neg,qoutt,eigval,qout)

        if(printable)then
           do ib=1,neg
              write(nfout,'(a,i8,a,f20.15)') &
             & ' !** eigen value of state ',ib,' : ',eigval(ib)
           enddo
        endif

        ! unitary transformation
        allocate(q1(np_ffth,neg));q1=zero
        call zgemm('N','N',np_ffth,neg,neg,one,qtran,np_ffth,qout,neg,zero,q1,np_ffth)
        deallocate(qtran)

        ! back to G-space
        allocate(q(nffth,np_e));q=0.d0
        call m_ES_W_transpose_back_fft_cmplx(nffth,np_ffth,q,q1)
        deallocate(q1)
        do ib=1,np_e
           do ig=1,nffth
              i1 = (ig-1)*kimg+1
              bfft(i1) = real(q(ig,ib))
              if(kimg==2) bfft(i1+1) = aimag(q(ig,ib))
           enddo
           call m_FFT_WF(ELECTRON,nfout,bfft,DIRECT,OFF)
           do ig=1,iba(ik)
              i1 = kimg*igf(nbase(ig,ik))-1
              zaj_l(ig,ib,ik,1) = bfft(i1)
              if(kimg==2) zaj_l(ig,ib,ik,kimg) = bfft(i1+1)
           enddo
        enddo
        deallocate(q)
     enddo

     call m_ES_modified_gram_schmidt(nfout)

     deallocate(qtmp)
     deallocate(qout)
     deallocate(bfft)
     deallocate(eigval)
     call tstatc0_end(id_sname)

  end subroutine m_ES_RSB_test

  subroutine m_ES_RSB_test_xyz()
     integer :: ik,ib,ig,i1,ib1,ib2

     complex(kind=DP), allocatable, dimension(:,:) :: q,q1,qtran
     complex(kind=DP), allocatable, dimension(:,:) :: qvec
     complex(kind=DP), allocatable, dimension(:,:,:) :: qout
     complex(kind=DP), allocatable, dimension(:) :: qtmp

     real(kind=DP), allocatable, dimension(:) :: bfft
     real(kind=DP), allocatable, dimension(:,:) :: eigval
     real(kind=DP) :: normr,normi,norm
     integer :: nffth
     integer :: modnrmt
     integer, dimension(2) :: regionl,regionm,regionn

     integer :: id_sname=-1

     call tstatc0_begin('m_ES_RSB_test ',id_sname,1)
     nffth = nfft/2

     allocate(qout(neg,neg,3));qout=zero
     allocate(qvec(neg,neg));qvec=zero
     allocate(qtmp(nffth));qtmp=zero
     allocate(bfft(nfft));bfft=0.d0
     allocate(eigval(neg,3));eigval=0.d0

     if(modnrm==EXECUT)then
        modnrmt = modnrm
        modnrm = NORMCONSERVATION
        call m_ES_modified_gram_schmidt(nfout)
        modnrm = modnrmt
     endif

     regionl(1) = 1
     regionl(2) = fft_box_size_WF(1,1)
     regionm(1) = 1
     regionm(2) = fft_box_size_WF(2,1)
     regionn(1) = 1
     regionn(2) = fft_box_size_WF(3,1)

     do ik=1,kv3,af+1
        if(map_k(ik)/=myrank_k) cycle

        ! WF in rspace
        allocate(q(nffth,np_e));q=zero
        do ib=ista_e,iend_e,istep_e
           call m_ES_WF_in_Rspace(ik,ib,bfft)
           normr=0.d0;normi=0.d0;qtmp=zero
           do ig=1,nffth
              i1 = (ig-1)*kimg+1
              qtmp(ig) = dcmplx(bfft(i1),bfft(i1+1))
              normr = normr+bfft(i1)*bfft(i1)
              if(kimg==2) normi = normi+bfft(i1+1)*bfft(i1+1)
           enddo
           norm = normr+normi
           norm = 1.d0/dsqrt(norm)
           qtmp(1:nffth) = norm*qtmp(1:nffth)
           q(1:nffth,map_z(ib)) = qtmp(1:nffth)
        enddo
        allocate(qtran(np_ffth,neg));qtran=zero
        call m_ES_W_transpose_fft_cmplx(nffth,np_ffth,q,qtran)
        deallocate(q)

        ! bisection
        call bisection_xyz(regionl,regionm,regionn,np_ffth,qtran,qout,1)

        ! diagonalize the matrix qout
        call simultaneous_diagonalization(3,neg,qout,eigval,qvec)

        !if(printable)then
        !   allocate(localization_info(lmax_rsb,3,neg))
        !   call build_localization_info(neg,lmax_rsb,eigval,localization_info)
        !   call print_localization_info(neg,lmax_rsb,localization_info,nfout)
        !   deallocate(localization_info)
        !endif

        ! unitary transformation
        allocate(q1(np_ffth,neg));q1=zero
        call zgemm('N','N',np_ffth,neg,neg,one,qtran,np_ffth,qvec,neg,zero,q1,np_ffth)
        deallocate(qtran)

        ! back to G-space
        allocate(q(nffth,np_e));q=0.d0
        call m_ES_W_transpose_back_fft_cmplx(nffth,np_ffth,q,q1)
        deallocate(q1)
        do ib=1,np_e
           do ig=1,nffth
              i1 = (ig-1)*kimg+1
              bfft(i1) = real(q(ig,ib))
              if(kimg==2) bfft(i1+1) = aimag(q(ig,ib))
           enddo
           call m_FFT_WF(ELECTRON,nfout,bfft,DIRECT,OFF)
           do ig=1,iba(ik)
              i1 = kimg*igf(nbase(ig,ik))-1
              zaj_l(ig,ib,ik,1) = bfft(i1)
              if(kimg==2) zaj_l(ig,ib,ik,kimg) = bfft(i1+1)
           enddo
        enddo
        deallocate(q)
     enddo

     call m_ES_modified_gram_schmidt(nfout)

     deallocate(qtmp)
     deallocate(qout)
     deallocate(qvec)
     deallocate(bfft)
     deallocate(eigval)
     call tstatc0_end(id_sname)
  end subroutine m_ES_RSB_test_xyz

  subroutine bisection_xyz(regionl,regionm,regionn,np_ffth,qtran,qout,lmax_rsb)
     integer, dimension(2), intent(in) :: regionl,regionm,regionn
     integer, intent(in) :: np_ffth
     complex(kind=DP),dimension(np_ffth,neg), intent(in) :: qtran
     complex(kind=DP), dimension(neg,neg,3), intent(out) :: qout
     integer, intent(in) :: lmax_rsb
!     call bisection(regionl,regionm,regionn,np_ffth,0,qtran,qout(:,:,1),phase)
!     call bisection(regionl,regionm,regionn,np_ffth,1,qtran,qout(:,:,2),phase)
!     call bisection(regionl,regionm,regionn,np_ffth,2,qtran,qout(:,:,3),phase)
     call bisection_by_walsh(np_ffth,neg,0,qtran,qout(:,:,1),lmax_rsb)
     call bisection_by_walsh(np_ffth,neg,1,qtran,qout(:,:,2),lmax_rsb)
     call bisection_by_walsh(np_ffth,neg,2,qtran,qout(:,:,3),lmax_rsb)
  end subroutine bisection_xyz

  subroutine bisection_by_walsh(np_ffth,nval,bisection_plane,qtran,qout,lmax_rsb)
     integer, intent(in) :: np_ffth,nval
     integer, intent(in) :: bisection_plane
     complex(kind=DP),dimension(np_ffth,nval), intent(in) :: qtran
     complex(kind=DP),dimension(nval,nval), intent(out) :: qout
     integer, intent(in) :: lmax_rsb

     integer :: ig
     integer :: ierr
     integer :: n1
     integer :: il,im,in
     integer :: fl,fm,fn
     integer :: lh,mh,nh
     complex(kind=DP), allocatable, dimension(:,:) :: q1
     real(kind=DP), allocatable, dimension(:,:) :: rq1,iq1,qout1,qout2,qout3,qout4
     complex(kind=DP) :: qtmp
     integer :: lmin,lmax,mmin,mmax,nmin,nmax
     integer :: id,mm,nl,nm,nn,nlhf
     integer :: ifft
     integer :: id_sname=-1
     real(kind=DP) :: x,denom

     call tstatc0_begin('bisection_by_walsh ',id_sname,1)

     id = fft_box_size_WF(1,0)
     mm = fft_box_size_WF(2,0)
     nl = fft_box_size_WF(1,1)
     nm = fft_box_size_WF(2,1)
     nn = fft_box_size_WF(3,1)
     if(kimg==1)then
       nlhf = id/2
     else
       nlhf = id
     endif
     n1 = 0
     allocate(q1(np_ffth,nval));q1=zero
     if(bisection_plane==0)then
        denom = 1.d0/dble(nl)
     else if (bisection_plane==1)then
        denom = 1.d0/dble(nm)
     else if (bisection_plane==2)then
        denom = 1.d0/dble(nn)
     endif

     do il=1,nl
        do im=1,nm
           do in=1,nn
              ig = nlhf*mm*(in-1)+nlhf*(im-1)+il
              if(ig.lt.ista_ffth.or.ig.gt.iend_ffth) cycle
              if(bisection_plane==0) then
                 x = dble(il)*denom
              else if(bisection_plane==1)then
                 x = dble(im)*denom
              else if(bisection_plane==2)then
                 x = dble(in)*denom
              endif
              if(walsh(lmax_rsb,x))then
                 n1=n1+1
                 q1(n1,1:nval) = qtran(ig-ista_ffth+1,1:nval)
              endif
           enddo
        enddo
     enddo

     ! qout = Q1H dot Q1
     qout=zero
#ifdef RSB_DEBUG
     do il=1,nval
        do im=1,nval
           qtmp=zero
           do in=1,n1
              qtmp = qtmp + dconjg(q1(in,il))*q1(in,im)
           enddo
           qout(il,im) = qtmp
        enddo
     enddo
#else
     if(n1>0)then
        call flush(0)
        allocate(rq1(n1,nval));rq1=0.d0
        allocate(iq1(n1,nval));iq1=0.d0
        allocate(qout1(nval,nval));qout1=0.d0
        allocate(qout2(nval,nval));qout2=0.d0
        allocate(qout3(nval,nval));qout3=0.d0
        allocate(qout4(nval,nval));qout4=0.d0
        rq1(1:n1,1:nval) = real(q1(1:n1,1:nval))
        iq1(1:n1,1:nval) = aimag(q1(1:n1,1:nval))
        call dgemm('T','N',nval,nval,n1,1.d0,rq1,n1,rq1,n1,0.d0,qout1,nval)
        call dgemm('T','N',nval,nval,n1,1.d0,iq1,n1,iq1,n1,0.d0,qout2,nval)
        call dgemm('T','N',nval,nval,n1,1.d0,iq1,n1,rq1,n1,0.d0,qout3,nval)
        call dgemm('T','N',nval,nval,n1,1.d0,rq1,n1,iq1,n1,0.d0,qout4,nval)
        qout = cmplx(qout1+qout2,-qout3+qout4)
        deallocate(rq1,iq1)
        deallocate(qout1,qout2,qout3,qout4)
     endif
#endif
     call mpi_allreduce(MPI_IN_PLACE,qout,nval*nval,mpi_double_complex,&
     &    mpi_sum,mpi_k_world(myrank_k),ierr)
     deallocate(q1)

     !currnval = 0.5*totch
        do il=1,nval
           do im=1,nval
              if(il.gt.currnval.and.im.le.currnval .or.  il.le.currnval.and.im.gt.currnval) qout(il,im)=zero
           enddo
        enddo

     call tstatc0_end(id_sname)
  end subroutine bisection_by_walsh

  subroutine bisection(regionl,regionm,regionn,np_ffth,bisection_plane,qtran,qout,phase)
     integer,dimension(2),intent(in) :: regionl,regionm,regionn
     integer, intent(in) :: np_ffth
     integer, intent(in) :: bisection_plane
     complex(kind=DP),dimension(np_ffth,neg), intent(in) :: qtran
     complex(kind=DP),dimension(neg,neg), intent(out) :: qout
     integer, intent(in) :: phase

     integer :: ig
     integer :: ierr
     integer :: n1
     integer :: il,im,in
     integer :: fl,fm,fn
     integer :: lh,mh,nh
     complex(kind=DP), allocatable, dimension(:,:) :: q1
     integer :: lmin,lmax,mmin,mmax,nmin,nmax
     integer :: id,mm,nl,nm,nn,nlhf

     integer :: id_sname=-1

     call tstatc0_begin('bisection ',id_sname,1)

     id = fft_box_size_WF(1,0)
     mm = fft_box_size_WF(2,0)
     nl = fft_box_size_WF(1,1)
     nm = fft_box_size_WF(2,1)
     nn = fft_box_size_WF(3,1)
     if(kimg==1)then
       nlhf = id/2
     else
       nlhf = id
     endif
     lmin=regionl(1);lmax=regionl(2)
     mmin=regionm(1);mmax=regionm(2)
     nmin=regionn(1);nmax=regionn(2)
     if(phase==0)then
        if(bisection_plane==0)then
           lmax = (lmax-lmin+1)/2
        else if (bisection_plane==1)then
           mmax = (mmax-mmin+1)/2
        else if (bisection_plane==2)then
           nmax = (nmax-nmin+1)/2
        endif
     else if (phase==1)then
        if(bisection_plane==0)then
           lmin = (lmax-lmin)/2
        else if (bisection_plane==1)then
           mmin = (mmax-mmin)/2
        else if (bisection_plane==2)then
           nmin = (nmax-nmin)/2
        endif
     endif
     n1 = 0
     allocate(q1(np_ffth,neg));q1=zero

     do il=lmin,lmax
        do im=mmin,mmax
           do in=nmin,nmax
              ig = nlhf*mm*(in-1)+nlhf*(im-1)+il
              if(ig.lt.ista_ffth.or.ig.gt.iend_ffth) cycle
              n1=n1+1
              q1(n1,1:neg) = qtran(ig-ista_ffth+1,1:neg)
           enddo
        enddo
     enddo

     ! qout = Q1H dot Q1
     qout=zero
     call zgemm('C','N',neg,neg,n1,one,q1,np_ffth,q1,np_ffth,zero,qout,neg)
     call mpi_allreduce(MPI_IN_PLACE,qout,neg*neg,mpi_double_complex,&
     &    mpi_sum,mpi_k_world(myrank_k),ierr)
     deallocate(q1)

     call tstatc0_end(id_sname)
  end subroutine bisection

  subroutine zheev_driver(ndim,eigval,cmtrx)
   integer, intent(in) :: ndim
   real(kind=DP), dimension(ndim), intent(out) :: eigval
   complex(kind=DP), dimension(ndim,ndim), intent(inout) :: cmtrx
   character(len=1) :: JOBZ,UPLO
!   integer, parameter :: lda=ndim
   integer :: lda
! --> T. Yamasaki 2009/07/24 revised according to a Dr. Katagiri's report
!!$   integer :: lwork,lwork_min,rwork_size
   integer :: lwork_min,rwork_size
   integer, save :: lwork = 0
! <--
   integer, parameter :: lwork_huge=10**8
   real(kind=DP),allocatable,dimension(:) :: work_lapack
   real(kind=DP),allocatable,dimension(:) :: rwork_lapack
   integer :: info

   integer :: id_sname = -1
   call tstatc0_begin('zheev_driver ', id_sname)

      lda=ndim
!(LAPACK)  JOBZ = N : eigenvalue, V : eigenvalue + eigenvector
     JOBZ = 'V'
!(LAPACK)  UPLO = U : upper triangle matrix,  L : lower triangle matrix
     UPLO = 'U'
      lwork_min=max(1,2*ndim-1)
      rwork_size=max(3*ndim-2,1)
     if(lwork < lwork_min)  lwork=lwork_min
     if(lwork > lwork_huge) lwork=lwork_min

     allocate(work_lapack(lwork*2))
     allocate(rwork_lapack(rwork_size))

     call zheev(JOBZ,UPLO,ndim,cmtrx,lda,eigval,work_lapack,lwork,rwork_lapack,info)
     lwork= int(work_lapack(1))+1

     deallocate(work_lapack)
     deallocate(rwork_lapack)

     call tstatc0_end(id_sname)
  end subroutine zheev_driver

  subroutine dsyev_driver(ndim,eig,w1hw2)
    integer, intent(in):: ndim
    real(kind=DP), intent(out) ,dimension(ndim) :: eig
    real(kind=DP), intent(inout) ,dimension(ndim,ndim) :: w1hw2

    integer, parameter :: PRINTLEVEL = 2
    character(len=1) :: JOBZ,UPLO
!   integer, parameter :: lda=ndim
    integer :: lda
!!$    integer :: lwork,lwork_min
    integer, save :: lwork=0
    integer       :: lwork_min
    integer,parameter :: nb = 64
    integer, parameter :: lwork_huge=10**8
    real(kind=DP),allocatable,dimension(:) :: work_lapack
    integer :: info

    integer :: id_sname = -1
    call tstatc0_begin('dsyev_driver ', id_sname)

    lda=ndim
!(LAPACK)  JOBZ = N : eigenvalue, V : eigenvalue + eigenvector
    JOBZ = 'V'
!(LAPACK)  UPLO = U : upper triangle matrix,  L : lower triangle matrix
    UPLO = 'U'
! (define lwork)
!!$    lwork = 0
!!$    lwork_min=max(1,3*ndim-1)
    lwork_min=max(1,(nb+2)*ndim)
    if(lwork == 0) then
       lwork = -1
       allocate(work_lapack(1))
       call dsyev(JOBZ,UPLO,ndim,w1hw2,lda,eig,work_lapack,lwork,info)
       lwork = int(work_lapack(1))+1
       deallocate(work_lapack)
    end if
    if(lwork < lwork_min) lwork=lwork_min
    if(lwork > lwork_huge) lwork=lwork_min
    allocate(work_lapack(lwork))

    call dsyev(JOBZ,UPLO,ndim,w1hw2,lda,eig,work_lapack,lwork,info)
!!$    write(*,*) 'lwork& info & suggested lwork = ',lwork,info,work_lapack(1)
!!$    lwork= int(work_lapack(1))+1
    deallocate(work_lapack)

    call tstatc0_end(id_sname)
  end subroutine dsyev_driver

  subroutine get_rot_mat(ndim,nmatrices,matrices,i,j,rii,rjj,rij,rji)
     integer, intent(in) :: ndim,nmatrices
     complex(kind=8), dimension(ndim,ndim,nmatrices) :: matrices
     integer, intent(in) :: i,j
     complex(kind=8), intent(out) :: rii,rjj,rij,rji
     integer :: k,ih,jh
     real(kind=DP), dimension(3,3) :: Gmat
     complex(kind=DP), dimension(3) :: hvec
     real(kind=DP), dimension(3) :: Geigval
     real(kind=DP), dimension(3,3) :: hhh
     complex(kind=DP) :: aii,ajj,aij,aji
     complex(kind=DP) :: c,s
     real(kind=DP) :: x,y,z

     ! calculate the rotation angle
     Gmat(:,:) = 0.d0
     do k=1,nmatrices
        aii=matrices(i,i,k)
        ajj=matrices(j,j,k)
        aij=matrices(i,j,k)
        aji=matrices(j,i,k)
        hvec(1) = aii-ajj
        hvec(2) = aij+aji
        hvec(3) = onei*(aji-aij)
        do ih=1,3
           do jh=1,3
              hhh(ih,jh) = dconjg(hvec(ih))*hvec(jh)
           enddo
        enddo
        Gmat(:,:) = Gmat(:,:)+real(hhh(:,:))
     enddo
     Geigval=0.d0
     call dsyev_driver(3,Geigval,Gmat)
     if(Gmat(1,3)<0) Gmat(:,:) = -Gmat(:,:)
     x = Gmat(1,3)
     y = Gmat(2,3)
     z = Gmat(3,3)

     c = dcmplx(dsqrt(0.5d0+0.5d0*x),0.d0)
     s = dcmplx(0.5d0*y,-0.5d0*z)/c

     rii = c
     rji = s
     rij = -dconjg(s)
     rjj = c
  end subroutine get_rot_mat

  subroutine gather_curr_indpair(i,j,currindpair)
     integer, intent(in) :: i,j
     integer, dimension(2,0:nrank_e-1), intent(out) :: currindpair
     integer, allocatable, dimension(:,:) :: currindpairb
     integer :: ierr
     if(nrank_e>1)then
        allocate(currindpairb(2,0:nrank_e-1));currindpairb=0
        currindpairb(1,myrank_e) = i
        currindpairb(2,myrank_e) = j
        call mpi_allgather(currindpairb(1,myrank_e),2,mpi_integer, &
                         & currindpair(1,0),2,mpi_integer,mpi_k_world(myrank_k),ierr)
        deallocate(currindpairb)
     else
        currindpair(1,0) = i
        currindpair(2,0) = j
     endif
  end subroutine gather_curr_indpair

  subroutine gather_rot_mat(rii,rjj,rij,rji,rotmat)
     complex(kind=DP), intent(in) :: rii,rjj,rij,rji
     complex(kind=DP), dimension(4,0:nrank_e-1), intent(out) :: rotmat
     complex(kind=DP), allocatable, dimension(:,:) :: rotmatb
     integer :: ierr
     if(nrank_e>1)then
        allocate(rotmatb(4,0:nrank_e-1));rotmatb=zero
        rotmatb(1,myrank_e) = rii
        rotmatb(2,myrank_e) = rjj
        rotmatb(3,myrank_e) = rij
        rotmatb(4,myrank_e) = rji
        rotmat=zero
        call mpi_allgather(rotmatb(1,myrank_e),4,mpi_double_complex, &
                         & rotmat(1,0),4,mpi_double_complex,mpi_k_world(myrank_k),ierr)
        deallocate(rotmatb)
     else
        rotmat(1,0) = rii
        rotmat(2,0) = rjj
        rotmat(3,0) = rij
        rotmat(4,0) = rji
     endif
  end subroutine gather_rot_mat

  ! subroutine which approximately diagonalizes a set of matrices simultaneously
  subroutine simultaneous_diagonalization_dbg(nmatrices,ndim,matrices,eigval,eigvec,verbose)
     integer, intent(in) :: nmatrices,ndim
     complex(kind=DP), dimension(ndim,ndim,nmatrices),intent(inout) :: matrices
     real(kind=DP),    dimension(ndim,nmatrices),     intent(out)   :: eigval
     complex(kind=DP), dimension(ndim,ndim),          intent(out)   :: eigvec
     logical, intent(in), optional                                  :: verbose
     real(kind=DP), dimension(3,3) :: Gmat
     real(kind=DP), dimension(3,3) :: hhh
     complex(kind=DP) :: aii,ajj,aij,aji,rii,rij,rji,rjj
     complex(kind=DP), dimension(3) :: hvec
     real(kind=DP), dimension(3) :: Geigval
     complex(kind=DP) :: c,s
     real(kind=DP) :: x,y,z
     integer :: isweep
     integer :: i,j,k,ih,jh,kh,ii,jj
     integer :: nsweep = 1000
     real(kind=DP) :: sum_diag,sum_diago,sum_offdiag,sum_offdiago
     real(kind=DP) :: eps_diag=1.d-10
     complex(kind=DP),allocatable,dimension(:) :: ctmpval1,ctmpval2
     integer :: nblock
     integer, allocatable, dimension(:,:) :: indpair
     complex(kind=DP), allocatable, dimension(:,:,:) :: matbuf
     logical :: verb=.false.
     integer :: id_sname=-1
     call tstatc0_begin('simultaneous diagonalization ',id_sname,1)
     if(present(verbose))then
        verb = verbose
     endif
     allocate(ctmpval1(ndim))
     allocate(ctmpval2(ndim))
     eigvec(1:ndim,1:ndim) = zero
     do ih=1,ndim
        eigvec(ih,ih) = one
     enddo
     nblock = (ndim/2)/nrank_e
     allocate(indpair(nblock,2))
     sum_diago=0.d0
     sum_offdiago=0.d0
     do k=1,nmatrices
        do ih=1,ndim-1
           sum_diago = sum_diago + real(dconjg(matrices(ih,ih,k))*matrices(ih,ih,k))
           do jh=ih+1,ndim
              sum_offdiago = sum_offdiago+real(dconjg(matrices(ih,jh,k))*matrices(ih,jh,k))
           enddo
        enddo
        sum_diago = sum_diago + real(dconjg(matrices(ndim,ndim,k))*matrices(ndim,ndim,k))
     enddo
     if(printable.and.verb.or.iprirsb>=2)then
        write(nfout,'(a)') '!-- simultaneous diagonalization'
        write(nfout,'(a,i5,i8)') '!-- number of matrices and their dimension : ',nmatrices,ndim
        write(nfout,'(a,3f20.10)') &
        & '!-- initial sum of the square of the diagonal, offdiagonal, and total elements :        ' &
        & ,sum_diago,sum_offdiago*2,sum_diago+sum_offdiago*2
     endif
     allocate(matbuf(ndim,ndim,nmatrices));matbuf=zero
     call initial_indices(ndim,nblock,myrank_e,indpair)
     do isweep=1,nsweep
        do ii=1,ndim-1
!           do j=i+1,ndim
           do jj=1,nblock
              do k=1,nblock
                 matbuf(:,indpair(k,1),:) = matrices(:,indpair(k,1),:)
                 matbuf(:,indpair(k,2),:) = matrices(:,indpair(k,2),:)
              enddo
              i = indpair(jj,1)
              j = indpair(jj,2)
              Gmat(:,:) = 0.d0
              do k=1,nmatrices
                 aii=matbuf(i,i,k)
                 ajj=matbuf(j,j,k)
                 aij=matbuf(i,j,k)
                 aji=matbuf(j,i,k)
                 hvec(1) = aii-ajj
                 hvec(2) = aij+aji
                 hvec(3) = onei*(aji-aij)
                 do ih=1,3
                    do jh=1,3
                       hhh(ih,jh) = dconjg(hvec(ih))*hvec(jh)
                    enddo
                 enddo
                 Gmat(:,:) = Gmat(:,:)+real(hhh(:,:))
              enddo
              Geigval=0.d0
              call dsyev_driver(3,Geigval,Gmat)
              if(Gmat(1,3)<0) Gmat(:,:) = -Gmat(:,:)
              x = Gmat(1,3)
              y = Gmat(2,3)
              z = Gmat(3,3)

              c = dcmplx(dsqrt(0.5d0+0.5d0*x),0.d0)
              s = dcmplx(0.5d0*y,-0.5d0*z)/c

              rii = c
              rji = s
              rij = -dconjg(s)
              rjj = c
              do k=1,nmatrices
                 do ih=1,ndim
                    ctmpval1(ih) = dconjg(rii)*matbuf(i,ih,k)+dconjg(rji)*matbuf(j,ih,k)
                    ctmpval2(ih) = dconjg(rij)*matbuf(i,ih,k)+dconjg(rjj)*matbuf(j,ih,k)
                 enddo
                 matrices(i,1:ndim,k) = ctmpval1(1:ndim)
                 matrices(j,1:ndim,k) = ctmpval2(1:ndim)
                 do ih=1,ndim
                    ctmpval1(ih) = matrices(ih,i,k)*rii+matrices(ih,j,k)*rji
                    ctmpval2(ih) = matrices(ih,i,k)*rij+matrices(ih,j,k)*rjj
                 enddo
                 matrices(1:ndim,i,k) = ctmpval1(1:ndim)
                 matrices(1:ndim,j,k) = ctmpval2(1:ndim)
              enddo
              do ih=1,ndim
                 ctmpval1(ih) = eigvec(ih,i)*rii+eigvec(ih,j)*rji
                 ctmpval2(ih) = eigvec(ih,i)*rij+eigvec(ih,j)*rjj
              enddo
              eigvec(1:ndim,i) = ctmpval1(1:ndim)
              eigvec(1:ndim,j) = ctmpval2(1:ndim)
           enddo
           call send_recv_and_rotate(myrank_e,nrank_e,nblock,indpair,ndim,nmatrices,matrices,eigvec)
        enddo
        sum_diag=0.d0
        sum_offdiag=0.d0
        do k=1,nmatrices
           do ih=1,ndim-1
              sum_diag = sum_diag + real(dconjg(matrices(ih,ih,k))*matrices(ih,ih,k))
              do jh=ih+1,ndim
                 sum_offdiag = sum_offdiag+real(dconjg(matrices(ih,jh,k))*matrices(ih,jh,k))
              enddo
           enddo
           sum_diag = sum_diag + real(dconjg(matrices(ndim,ndim,k))*matrices(ndim,ndim,k))
        enddo
        if(printable.and.verb.or.iprirsb>=2) then
        write(nfout,'(a,i5,a,3f20.10)') &
        & '!-- sweep : ',isweep,', sum of the square of the diagonal, offdiagonal, and total elements : ' &
        & ,sum_diag,sum_offdiag*2,sum_diag+sum_offdiag*2
        endif
        if(abs(sum_diago-sum_diag)/(sum_diag+sum_offdiag*2)<eps_diag) then
           if(printable.and.verb.or.iprirsb>=2) write(nfout,'(a)') '!-- convergence reached'
           exit
        endif
        sum_diago = sum_diag
     enddo
     do k=1,nmatrices
        do i=1,ndim
           eigval(i,k) = real(matrices(i,i,k))
        enddo
     enddo
     deallocate(ctmpval1)
     deallocate(ctmpval2)
     call tstatc0_end(id_sname)
  end subroutine simultaneous_diagonalization_dbg

  ! subroutine which approximately diagonalizes a set of matrices simultaneously
  subroutine simultaneous_diagonalization(nmatrices,ndim,matrices,eigval,eigvec,verbose)
     integer, intent(in) :: nmatrices,ndim
     complex(kind=DP), dimension(ndim,ndim,nmatrices),intent(inout) :: matrices
     real(kind=DP),    dimension(ndim,nmatrices),     intent(out)   :: eigval
     complex(kind=DP), dimension(ndim,ndim),          intent(out)   :: eigvec
     logical, intent(in),optional                                   :: verbose
     complex(kind=DP) :: rii,rij,rji,rjj,crii,crij,crji,crjj,mi,mj,ccrii,ccrij,ccrji,ccrjj
     complex(kind=DP) :: criic,crijc,crjic,crjjc
     complex(kind=DP) :: mii,mij,mji,mjj
     integer :: isweep
     integer :: i,j,k,ih,jh,kh,ii,jj
     integer :: nsweep = 1000
     real(kind=DP) :: sum_diag,sum_diago,sum_offdiag,sum_offdiago
     real(kind=DP) :: eps_diag=1.d-8
     complex(kind=DP),allocatable,dimension(:) :: ctmpval1,ctmpval2
     integer, allocatable, dimension(:,:) :: indpair
     complex(kind=DP), allocatable, dimension(:,:) :: rotmat
     complex(kind=DP), allocatable, dimension(:,:,:) :: mbuf
     integer, allocatable, dimension(:,:) :: currindpair
     complex(kind=DP), allocatable, dimension(:,:) :: eigvec_buf
     integer :: nblock,ierr,irank,itmp,jtmp,i1,j1,kk
     logical :: verb=.false.
     integer :: id_sname=-1
     call tstatc0_begin('simultaneous diagonalization ',id_sname,1)
     if(present(verbose))then
        verb = verbose
     endif
     allocate(ctmpval1(ndim))
     allocate(ctmpval2(ndim))

     nblock = (ndim/2)/nrank_e
     allocate(indpair(nblock,2))

     eigvec(1:ndim,1:ndim) = zero
     do ih=1,ndim
        eigvec(ih,ih) = one
     enddo
     sum_diago=0.d0
     sum_offdiago=0.d0
     do k=1,nmatrices
        do ih=1,ndim-1
           sum_diago = sum_diago + real(dconjg(matrices(ih,ih,k))*matrices(ih,ih,k))
           do jh=ih+1,ndim
              sum_offdiago = sum_offdiago+real(dconjg(matrices(ih,jh,k))*matrices(ih,jh,k))
           enddo
        enddo
        sum_diago = sum_diago + real(dconjg(matrices(ndim,ndim,k))*matrices(ndim,ndim,k))
     enddo
     if(printable.and.verb.or.iprirsb>=2)then
        write(nfout,'(a)') '!-- simultaneous diagonalization'
        write(nfout,'(a,i5,i8)') '!-- number of matrices and their dimension : ',nmatrices,ndim
        write(nfout,'(a,3f20.10)') &
        & '!-- initial sum of the square of the diagonal, offdiagonal, and total elements :        ' &
        & ,sum_diago,sum_offdiago*2,sum_diago+sum_offdiago*2
     endif

     call initial_indices(ndim,nblock,myrank_e,indpair)
     allocate(rotmat(4,0:nrank_e-1));rotmat=zero
     allocate(currindpair(2,0:nrank_e-1));currindpair=0
     allocate(mbuf(ndim,2,nmatrices));mbuf=zero
     do isweep=1,nsweep
        do ii=1,ndim-1
           do jj=1,nblock

              i = indpair(jj,1)
              j = indpair(jj,2)
              mbuf(:,1,:)=matrices(:,i,:)
              mbuf(:,2,:)=matrices(:,j,:)
              call gather_curr_indpair(i,j,currindpair)
              call get_rot_mat(ndim,nmatrices,matrices,i,j,rii,rjj,rij,rji)
              call gather_rot_mat(rii,rjj,rij,rji,rotmat)
              ccrii=dconjg(rii);ccrij=dconjg(rij);ccrji=dconjg(rji);ccrjj=dconjg(rjj)
              do k=1,nmatrices
                 do ih=1,ndim
                    mi=matrices(ih,i,k);mj=matrices(ih,j,k)
                    matrices(ih,i,k) = mi*rii+mj*rji
                    matrices(ih,j,k) = mi*rij+mj*rjj
                 enddo
              enddo

              do k=1,nmatrices
                 do irank=0,nrank_e-1
                    itmp = currindpair(1,irank)
                    jtmp = currindpair(2,irank)
                    crii = rotmat(1,irank);criic=dconjg(crii)
                    crjj = rotmat(2,irank);crjjc=dconjg(crjj)
                    crij = rotmat(3,irank);crijc=dconjg(crij)
                    crji = rotmat(4,irank);crjic=dconjg(crji)
                    mii=dconjg(mbuf(itmp,1,k));mij=dconjg(mbuf(jtmp,1,k))
                    mjj=dconjg(mbuf(jtmp,2,k));mji=dconjg(mbuf(itmp,2,k))
                    matrices(itmp,i,k) = dconjg((mii*crii+mij*crji)*ccrii+(mji*crii+mjj*crji)*ccrji)
                    matrices(jtmp,i,k) = dconjg((mii*crij+mij*crjj)*ccrii+(mji*crij+mjj*crjj)*ccrji)
                    matrices(itmp,j,k) = dconjg((mii*crii+mij*crji)*ccrij+(mji*crii+mjj*crji)*ccrjj)
                    matrices(jtmp,j,k) = dconjg((mii*crij+mij*crjj)*ccrij+(mji*crij+mjj*crjj)*ccrjj)
                    do kk=1,nblock
                       if(kk.eq.jj) cycle
                       i1 = indpair(kk,1);j1=indpair(kk,2)
                       mii=matrices(itmp,i1,k);mji=matrices(jtmp,i1,k)
                       mjj=matrices(jtmp,j1,k);mij=matrices(itmp,j1,k)
                       matrices(itmp,i1,k) = criic*mii+crjic*mji
                       matrices(itmp,j1,k) = criic*mij+crjic*mjj
                       matrices(jtmp,i1,k) = crijc*mii+crjjc*mji
                       matrices(jtmp,j1,k) = crijc*mij+crjjc*mjj
                    enddo
                 enddo
              enddo

              do ih=1,ndim
                 ctmpval1(ih) = eigvec(ih,i)*rii+eigvec(ih,j)*rji
                 ctmpval2(ih) = eigvec(ih,i)*rij+eigvec(ih,j)*rjj
              enddo
              eigvec(1:ndim,i) = ctmpval1(1:ndim)
              eigvec(1:ndim,j) = ctmpval2(1:ndim)
           enddo
           call send_recv_and_rotate(myrank_e,nrank_e,nblock,indpair,ndim,nmatrices,matrices,eigvec)
        enddo
        sum_diag=0.d0
        sum_offdiag=0.d0
        do k=1,nmatrices
           do jh=1,nblock
              i = indpair(jh,1)
              j = indpair(jh,2)
              sum_diag = sum_diag + real(dconjg(matrices(i,i,k))*matrices(i,i,k))
              sum_diag = sum_diag + real(dconjg(matrices(j,j,k))*matrices(j,j,k))
              do kh=1,ndim
                 if(i/=kh) sum_offdiag = sum_offdiag+real(dconjg(matrices(kh,i,k))*matrices(kh,i,k))
              enddo
              do kh=1,ndim
                 if(j/=kh) sum_offdiag = sum_offdiag+real(dconjg(matrices(kh,j,k))*matrices(kh,j,k))
              enddo
           enddo
        enddo
        if(nrank_e>1)then
           call mpi_allreduce(MPI_IN_PLACE,sum_diag,1,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr)
           call mpi_allreduce(MPI_IN_PLACE,sum_offdiag,1,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr)
        endif
        if(printable.and.verb.or.iprirsb>=2) then
           write(nfout,'(a,i5,a,3f20.10)') &
        & '!-- sweep : ',isweep,', sum of the square of the diagonal, offdiagonal, and total elements : ' &
        & ,sum_diag,sum_offdiag,sum_diag+sum_offdiag
        endif
        if(abs(sum_diago-sum_diag)/(sum_diag+sum_offdiag)<eps_diag) then
           if(printable.and.verb.or.iprirsb>=2) write(nfout,'(a)') '!-- convergence reached'
           exit
        endif
        sum_diago = sum_diag
     enddo

     if(printable)then
        write(nfout,'(a,i5,a)') '!-- simultaneous diagonalization converged after ',isweep,' sweeps'
     endif

     do k=1,nmatrices
        do i=1,ndim
           eigval(i,k) = real(matrices(i,i,k))
        enddo
     enddo

     allocate(eigvec_buf(ndim,ndim));eigvec_buf=zero
     do i=1,nblock
        eigvec_buf(1:ndim,indpair(i,1)) = eigvec(1:ndim,indpair(i,1))
        eigvec_buf(1:ndim,indpair(i,2)) = eigvec(1:ndim,indpair(i,2))
     enddo
     if(nrank_e>1) call mpi_allreduce(MPI_IN_PLACE,eigvec_buf,ndim*ndim,mpi_double_complex,mpi_sum,mpi_k_world(myrank_k),ierr)
     eigvec = eigvec_buf
     deallocate(eigvec_buf)

     deallocate(rotmat)
     deallocate(mbuf)
     deallocate(ctmpval1)
     deallocate(ctmpval2)
     deallocate(indpair)
     deallocate(currindpair)
     call tstatc0_end(id_sname)
  end subroutine simultaneous_diagonalization

  subroutine initial_indices(n,nblock,p,indpair)
     implicit none
     integer, intent(in) :: n,nblock,p
     integer, dimension(nblock,2) :: indpair
     integer :: i,istart,iend
     istart = nblock*p
     do i=1,nblock
        indpair(i,1) = istart+i
     enddo
     iend = n-p*nblock
     do i=1,nblock
        indpair(i,2) = iend-(i-1)
     enddo
  end subroutine initial_indices

  subroutine send_recv_and_rotate(mype,npes,nblock,indpair,ndim,nmatrices,matrices,eigvec)
     implicit none
     integer, intent(in) :: mype,npes,nblock
     integer, dimension(nblock,2), intent(inout) :: indpair
     integer, intent(in) :: ndim,nmatrices
     complex(kind=DP), dimension(ndim,ndim,nmatrices), intent(inout) :: matrices
     complex(kind=DP), dimension(ndim,ndim), intent(inout) :: eigvec
     integer :: ii,i,j
     integer :: indbuf1,indbufn,sendbuf,recvbuf,ierr
     complex(kind=DP), allocatable, dimension(:,:) :: cbuf,cbufr
     integer, allocatable, dimension(:) :: ista
     allocate(ista(MPI_STATUS_SIZE))
     allocate(cbuf(1:ndim,1:nmatrices+1));cbuf=zero
     allocate(cbufr(1:ndim,1:nmatrices+1));cbufr=zero

     if(mype<npes-1)then
        sendbuf = indpair(nblock,1)
        cbuf(1:ndim,1:nmatrices) = matrices(1:ndim,sendbuf,1:nmatrices)
        cbuf(1:ndim,nmatrices+1) = eigvec(1:ndim,sendbuf)
        call mpi_send(sendbuf,1,mpi_integer,mype+1,0,mpi_k_world(myrank_k),ierr)
        call mpi_recv(recvbuf,1,mpi_integer,mype+1,0,mpi_k_world(myrank_k),ista,ierr)
        call mpi_send(cbuf,ndim*(nmatrices+1),mpi_double_complex,mype+1,0,mpi_k_world(myrank_k),ierr)
        call mpi_recv(cbufr,ndim*(nmatrices+1),mpi_double_complex,mype+1,0,mpi_k_world(myrank_k),ista,ierr)
        indpair(nblock,1) = recvbuf
        matrices(1:ndim,recvbuf,1:nmatrices) = cbufr(1:ndim,1:nmatrices)
        eigvec(1:ndim,recvbuf) = cbufr(1:ndim,nmatrices+1)
     endif
     if(mype>0)then
        sendbuf = indpair(1,2)
        cbuf(1:ndim,1:nmatrices) = matrices(1:ndim,sendbuf,1:nmatrices)
        cbuf(1:ndim,nmatrices+1) = eigvec(1:ndim,sendbuf)
        call mpi_recv(recvbuf,1,mpi_integer,mype-1,0,mpi_k_world(myrank_k),ista,ierr)
        call mpi_send(sendbuf,1,mpi_integer,mype-1,0,mpi_k_world(myrank_k),ierr)
        call mpi_recv(cbufr,ndim*(nmatrices+1),mpi_double_complex,mype-1,0,mpi_k_world(myrank_k),ista,ierr)
        call mpi_send(cbuf,ndim*(nmatrices+1),mpi_double_complex,mype-1,0,mpi_k_world(myrank_k),ierr)
        indpair(1,2) = recvbuf
        matrices(1:ndim,recvbuf,1:nmatrices) = cbufr(1:ndim,1:nmatrices)
        eigvec(1:ndim,recvbuf) = cbufr(1:ndim,nmatrices+1)
     endif

     deallocate(ista)
     deallocate(cbuf)
     deallocate(cbufr)

     indbufn = indpair(nblock,1)
     if(mype==0)then
        do i=nblock,3,-1
           indpair(i,1) = indpair(i-1,1)
        enddo
     else
        do i=nblock,2,-1
           indpair(i,1) = indpair(i-1,1)
        enddo
     endif
     indbuf1 = indpair(1,2)
     do i=1,nblock-1
        indpair(i,2) = indpair(i+1,2)
     enddo
     if(mype==0)then
        indpair(2,1) = indbuf1
     else
        indpair(1,1) = indbuf1
     endif
     indpair(nblock,2) = indbufn
  end subroutine send_recv_and_rotate

  logical function walsh(l,x)
    implicit none
    integer, intent(in) :: l
    real(kind=DP), intent(in) :: x
    real(kind=DP) :: xx,dx,dx2
    real(kind=DP) :: wal
    if(l.eq.1)then
      if(x.lt.0.5d0) then
        walsh=.false.
      else
        walsh=.true.
      endif
      return
    endif
    dx = 1.d0/dble(2**l)
    if(x.lt.dx)then
      walsh=.false.
      return
    endif
    dx2 = 2.d0*dx
    xx = x+dx
    if(xx.gt.1) xx=xx-1
    wal = 2*(mod(floor(xx/dx2),2)-0.5)
    walsh = wal>0
  end function walsh

  !transform a vector to rsb space/submat space
  !imode==DIRECT  -> to rsb space,
  !imode==INVERSE -> to submat space
  subroutine m_ES_RSB_unitary_transform_vec(imode,ik,vec)
    integer,intent(in) :: imode,ik
    real(kind=DP), intent(inout), dimension(kg1,np_e,kimg) :: vec
    complex(kind=DP),allocatable,dimension(:,:) :: bufvec1,bufvec2
    integer :: npg1k
    npg1k = maxval(np_g1k)
    allocate(bufvec1(kg1,np_e));bufvec1=zero
    if(kimg==1)then
       bufvec1(:,:) = dcmplx(vec(:,:,1),0.d0)
    else
       bufvec1(:,:) = dcmplx(vec(:,:,1),vec(:,:,2))
    endif

    allocate(bufvec2(npg1k,neg));bufvec2=zero
    call m_ES_W_transpose_fft_cmplx(kg1,npg1k,bufvec1,bufvec2)
    deallocate(bufvec1)

    allocate(bufvec1(npg1k,neg));bufvec1=zero
    if(imode==DIRECT)then
        call zgemm('N','N',npg1k,nval,nval,one,bufvec2,npg1k,qout,nval,zero,bufvec1,npg1k)
    else if (imode==INVERSE)then
        call zgemm('N','C',npg1k,nval,nval,one,bufvec2,npg1k,qout,nval,zero,bufvec1,npg1k)
    endif
    deallocate(bufvec2)

    allocate(bufvec2(kg1,np_e));bufvec2=zero
    call m_ES_W_transpose_back_fft_cmplx(kg1,npg1k,bufvec2,bufvec1)
    deallocate(bufvec1)

    vec(:,:,1) = dble(bufvec2(:,:))
    if(kimg>1) vec(:,:,2) = dimag(bufvec2(:,:))
    deallocate(bufvec2)

  end subroutine m_ES_RSB_unitary_transform_vec

  !transform the wfs to rsb space/submat space
  !imode==DIRECT  -> to rsb space,
  !imode==INVERSE -> to submat space
  subroutine m_ES_RSB_unitary_transform_wf(imode,ik,zaj_l,do_betar)
    integer,intent(in) :: imode,ik
    logical,intent(in) :: do_betar
    real(kind=DP), intent(inout), dimension(kg1,np_e,ista_k:iend_k,kimg) :: zaj_l
    real(kind=DP),allocatable,dimension(:,:,:)  :: bufwf
    complex(kind=DP),allocatable,dimension(:,:) :: bufwf1,bufwf2

    call m_ESortho_set_np_g1k_x(ik) !-> np_g1k_x
    allocate(bufwf(np_g1k_x,neg,kimg));bufwf=0.d0
    call m_ES_W_transpose_r(.false.,ista_k,iend_k,ik,zaj_l,bufwf)
    allocate(bufwf1(np_g1k_x,neg));bufwf1=zero
    if(kimg==1)then
       bufwf1(:,:) = dcmplx(bufwf(:,:,1),0.d0)
    else
       bufwf1(:,:) = dcmplx(bufwf(:,:,1),bufwf(:,:,2))
    endif
    deallocate(bufwf)


    allocate(bufwf2(np_g1k_x,neg));bufwf2=zero
    if(imode==DIRECT)then
        call zgemm('N','N',np_g1k_x,nval,nval,one,bufwf1,np_g1k_x,qout,nval,zero,bufwf2,np_g1k_x)
    else if (imode==INVERSE)then
        call zgemm('N','C',np_g1k_x,nval,nval,one,bufwf1,np_g1k_x,qout,nval,zero,bufwf2,np_g1k_x)
    endif

    allocate(bufwf(np_g1k_x,neg,kimg));bufwf=0.d0
    bufwf(:,:,1) = dble(bufwf2(:,:))
    if(kimg>1) bufwf(:,:,2) = dimag(bufwf2(:,:))
    deallocate(bufwf2)

    call m_ES_W_transpose_back_r(.false.,ista_k,iend_k,ik,zaj_l,bufwf)
    deallocate(bufwf)
    if(do_betar) call m_ES_betar_dot_WFs_4_each_k(nfout,ik)   ! -> fsr_l,fsi_l
  end subroutine m_ES_RSB_unitary_transform_wf


end module m_ES_RSB

