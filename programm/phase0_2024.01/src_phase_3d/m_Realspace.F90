module m_Realspace
  use m_Const_Parameters,     only : DP, SOFTPART, HARDPART
  use m_PlaneWaveBasisSet,    only : kgp,kg1,ngabc
  use m_Crystal_Structure,    only : altv,univol
  use m_Ionic_System,         only : ntyp,pos,natm,ityp,iwei
  use m_Timing,               only : tstatc0_begin, tstatc0_end
  use m_Control_Parameters,   only : kimg,iprirs,printable,r0_factor,r0_factor_q
  use m_FFT,                  only : fft_box_size_WF,fft_box_size_CD,fft_box_size_CD_nonpara
  use m_Files,                only : nfout
  use m_PseudoPotential,      only : mmesh,nmesh,rmax,radr,wos,xh,ilmt,nqitg_sp &
  &                                , m_PP_tell_lmtt_l_m_tau, betar,qrspspw,il2p &
  &                                , ltp,taup,lmta,iqitg,m_PP_find_maximum_l    &
  &                                , isph, dl2p
  use m_Parallelization,      only : m_Parallel_init_mpi_rspace_aug

  use m_Const_Parameters,     only : ON
  use m_Control_Parameters,   only : sw_rspace_ekin_density
  use m_PseudoPotential,      only : wf_mnrc

  implicit none

  integer                                      :: nmesh_rs_max
  integer                                      :: nmesh_rs_max_h
  integer, allocatable, target, dimension(:)   :: nmesh_rs
  integer, allocatable, target, dimension(:)   :: nmesh_rs_h
  integer, allocatable, target, dimension(:,:) :: meshx_rs
  integer, allocatable, target, dimension(:,:) :: meshy_rs
  integer, allocatable, target, dimension(:,:) :: meshz_rs
  integer, allocatable, target, dimension(:,:) :: meshxyz_rs
  integer, allocatable, target, dimension(:,:) :: meshxyz_rs_h
  integer, allocatable, target, dimension(:,:) :: map_h
  integer, allocatable, target, dimension(:,:) :: map_h_i
  integer, allocatable, target, dimension(:,:) :: meshxyz_rs_conjg

  integer                                      :: nmesh_rs_aug_max
  integer                                      :: nmesh_rs_aug_max_h
  integer, allocatable, target, dimension(:)   :: nmesh_rs_aug
  integer, allocatable, target, dimension(:)   :: nmesh_rs_aug_h
  integer, allocatable, target, dimension(:,:) :: meshx_rs_aug
  integer, allocatable, target, dimension(:,:) :: meshy_rs_aug
  integer, allocatable, target, dimension(:,:) :: meshz_rs_aug
  integer, allocatable, target, dimension(:,:) :: meshxyz_rs_aug
  integer, allocatable, target, dimension(:,:) :: meshxyz_rs_aug_h
  integer, allocatable, target, dimension(:,:) :: map_aug_h
  integer, allocatable, target, dimension(:,:) :: map_aug_h_i
  integer, allocatable, target, dimension(:,:) :: meshxyz_rs_aug_conjg

  integer :: neix=2
  integer :: neiy=2
  integer :: neiz=2

  real(kind=DP), allocatable, dimension(:) :: rcut_betar, rcut_betar_rs
  real(kind=DP), allocatable, dimension(:) :: rcut_qr, rcut_qr_rs

  real(kind=DP), allocatable, dimension(:,:,:) :: qr_clm_ylm
  real(kind=DP), allocatable, dimension(:,:,:,:) :: dqr_clm_ylm
  
  integer, allocatable, dimension(:)   :: nlmtpair
  integer, allocatable, dimension(:,:) :: plmt1,plmt2

  interface m_RS_resolve_mesh_hard
    module procedure resolve_mesh_hard
    module procedure resolve_mesh_hard0
    module procedure resolve_mesh_hard1
  end interface m_RS_resolve_mesh_hard

  interface m_RS_build_qr_clm_ylm
    module procedure build_qr_clm_ylm0
    module procedure build_qr_clm_ylm1
  end interface m_RS_build_qr_clm_ylm

  contains

  subroutine m_RS_get_coords(nffth,flg,fmg,fng,ilg,img,ing)
     integer, intent(in) :: nffth 
     real(kind=DP), dimension(nffth), intent(out) :: flg,fmg,fng
     integer, dimension(nffth), intent(out) :: ilg,img,ing
     integer :: id,mm,nl,nm,nn
     real(kind=DP) :: invnl,invnm,invnn
     integer :: inl,inm,inn
     integer :: nlhf,ig
     flg=0.d0;fmg=0.d0;fng=0.d0
     id = fft_box_size_WF(1,0)
     mm = fft_box_size_WF(2,0)
     nl = fft_box_size_WF(1,1)
     nm = fft_box_size_WF(2,1)
     nn = fft_box_size_WF(3,1)
     invnl = 1.d0/dble(nl)
     invnm = 1.d0/dble(nm)
     invnn = 1.d0/dble(nn)
     if(kimg==1)then
       nlhf = id/2
     else
       nlhf = id
     endif
     do inl=1,nl
        do inm=1,nm
           do inn=1,nn
              ig = nlhf*mm*(inn-1)+nlhf*(inm-1)+inl
              flg(ig) = dble(inl)*invnl
              fmg(ig) = dble(inm)*invnm
              fng(ig) = dble(inn)*invnn
              ilg(ig) = inl
              img(ig) = inm
              ing(ig) = inn
           enddo
        enddo
     enddo
  end subroutine m_RS_get_coords

  subroutine resolve_cutoff_betar(nfout)
     integer, intent(in) :: nfout
     integer :: it,lmt,ilmtt,il,im,tau,nspher,imesh,icut,icutmax
     real(kind=DP)      :: hn,babs,eps,eps0
     real(kind=DP)      :: r0max,r1max
     if(.not.allocated(rcut_betar)) allocate(rcut_betar(ntyp));rcut_betar = 0.d0
     if(.not.allocated(rcut_betar_rs)) allocate(rcut_betar_rs(ntyp));rcut_betar_rs = 0.d0
     do it=1,ntyp
        call rmeshs(nmesh(it),nmesh(it),xh(it),rmax(it),radr,hn) ! -(b_PP)
        icutmax = -1
        do lmt=1,ilmt(it) 
           call m_PP_tell_lmtt_l_m_tau(lmt,it,ilmtt,il,im,tau,nspher)
           eps0 = 1.d0
           do imesh=nmesh(it),1,-1
              babs = abs(betar(imesh,il,tau,it))
              if(babs.eq.0.d0)then
                 eps = 1.d0
              else
                 eps = babs
              endif
              if(eps.gt.5.d0*eps0)then
                 if(babs.gt.1.d-10)then
                    if(imesh>=nmesh(it)-1)then
                       icut = imesh+1
                    else
                       icut = imesh+2
                    endif
                    exit
                 endif
              endif
              eps0 = eps
           enddo
           if(icutmax<icut)then
              icutmax = icut
           endif
        enddo
        rcut_betar(it) = radr(icutmax)
        rcut_betar_rs(it) = rcut_betar(it)*r0_factor
     enddo
     r0max=0.d0
     r1max=0.d0
     do it=1,ntyp
        if(r0max<rcut_betar(it))then
          r0max = rcut_betar(it)
        endif
     enddo
     !rcut_betar = r0max
     !rcut_betar_rs = r0max*r0_factor
     if(iprirs>=2.and.printable) then
       do it=1,ntyp
         write(nfout,'(a,i4,a,f10.5,a)') ' !RS cutoff for the projector of element           ',it,' : ',&
         & rcut_betar(it),' bohr'
         write(nfout,'(a,i4,a,f10.5,a)') ' !RS cutoff for the optimized projector of element ',it,' : ',&
         & rcut_betar_rs(it),' bohr'
       enddo
     endif
  end subroutine resolve_cutoff_betar

  subroutine resolve_cutoff_q(nfout)
     integer, intent(in) :: nfout
     integer :: it,lmt1,lmt2,np,iiqitg
     integer :: n,il1,tau1,ilmta1,il2,tau2,ilmta2,ilm3,l3
     integer :: imesh
     real(kind=DP) :: hn
     real(kind=DP),allocatable,dimension(:) :: radr
     real(kind=DP) :: rcut_tmp,qab,eps,eps0
     integer :: icut
     integer, allocatable :: il3(:)


     call m_PP_find_maximum_l(n)    ! n-1: maximum l
     n = (n-1) + (n-1) + 1
     allocate(il3(n**2))
     call substitute_il3(n**2,il3) ! -(b_Elec..)

     if(.not.allocated(rcut_qr)) allocate(rcut_qr(ntyp));rcut_qr=0.d0
     if(.not.allocated(rcut_qr_rs)) allocate(rcut_qr_rs(ntyp));rcut_qr_rs=0.d0
     allocate(radr(mmesh)) 
     do it=1,ntyp
        call rmeshs(nmesh(it),nmesh(it),xh(it),rmax(it),radr,hn) ! -(b_PP)
        rcut_tmp = 0.d0
        do lmt1=1,ilmt(it)
           il1  = ltp( lmt1, it)
           tau1 = taup(lmt1, it)
           do lmt2=1,ilmt(it)
              if(il2p(lmt1,lmt2,it)==0) cycle
              il2  = ltp( lmt2, it)
              tau2 = taup(lmt2, it)
              do np = 1, il2p(lmt1,lmt2,it)
                 ilm3 = isph(lmt1,lmt2,np,it)
                 l3 = il3(ilm3)
                 if(lmt2>=lmt1) then
                   iiqitg = iqitg(il1,tau1,il2,tau2,l3+1,it)
                 else
                   iiqitg = iqitg(il2,tau2,il1,tau1,l3+1,it)
                 end if
                 if(iiqitg == 0) cycle
                 eps0 = 1.d0
                 do imesh=nmesh(it),1,-1
                    if(iprirs>=3) write(500+iiqitg,'(2f20.10)') radr(imesh),qrspspw(imesh,iiqitg)
                    qab = abs(qrspspw(imesh,iiqitg))
                    if(qab.eq.0.d0)then
                       eps = 1.d0
                    else
                       eps = qab
                    endif
                    if(eps.gt.5.d0*eps0)then
                       if(qab.gt.1.d-10)then
                          if(imesh>=nmesh(it)-1)then
                             icut = imesh+1
                          else
                             icut = imesh+2
                          endif
                          exit
                       endif
                    endif
                    eps0 = eps
                 enddo
                 if(radr(icut)>rcut_tmp) rcut_tmp = radr(icut)
              enddo
           enddo
        enddo
        rcut_qr(it) = rcut_tmp
        rcut_qr_rs(it) = rcut_tmp * r0_factor_q
        if(iprirs>=2) write(nfout,'(a,i5,a,f10.5,a)') &
        & ' !RS cutoff of the aug. function for element ',it,' : ',rcut_qr(it),' bohr'
     enddo
     deallocate(il3)
     deallocate(radr)
  end subroutine resolve_cutoff_q

  subroutine resolve_cutoff_ekin_density(nfout)
    integer, intent(in) :: nfout

    integer :: it
    real(kind=DP) :: hn, rcut_tmp
    real(kind=DP),allocatable,dimension(:) :: radr

    if(.not.allocated(rcut_qr)) allocate(rcut_qr(ntyp));rcut_qr=0.d0
    if(.not.allocated(rcut_qr_rs)) allocate(rcut_qr_rs(ntyp));rcut_qr_rs=0.d0
    allocate(radr(mmesh))

    do it=1,ntyp
       call rmeshs(nmesh(it),nmesh(it),xh(it),rmax(it),radr,hn) ! -(b_PP)

       rcut_tmp = radr( wf_mnrc(it) )
       rcut_qr(it) = rcut_tmp
       rcut_qr_rs(it) = rcut_tmp * r0_factor_q
    end do
    deallocate(radr)

  end subroutine resolve_cutoff_ekin_density

  subroutine resolve_mesh_hard(nfout)
     integer, intent(in) :: nfout
     call resolve_mesh_hard1(nfout,fft_box_size_CD,fft_box_size_CD_nonpara)
  end subroutine resolve_mesh_hard

  subroutine resolve_mesh_hard1(nfout,box_size,box_size_nonpara)
     integer, intent(in) :: nfout
     integer, dimension(3,0:1), intent(in) :: box_size
     integer, dimension(3,0:0), intent(in) :: box_size_nonpara
     integer :: nl,nm,nn,mm,id
     id = box_size_nonpara(1,0)
     mm = box_size_nonpara(2,0)
     nl = box_size(1,1)
     nm = box_size(2,1)
     nn = box_size(3,1)
     call resolve_mesh_hard0(nfout,id,mm,nl,nm,nn)
  end subroutine resolve_mesh_hard1

  subroutine resolve_mesh_hard0(nfout,id,mm,nl,nm,nn)
     integer, intent(in) :: nfout
     integer, intent(in) :: id,mm,nl,nm,nn
     real(kind=DP) :: ntot
     integer :: id_sname = -1
     logical,save :: firstcall=.true.
     call tstatc0_begin('m_RS_resolve_mesh_hard ',id_sname,1)
     if (firstcall) then
#if 1
        call resolve_cutoff_q(nfout)
#else
        if ( sw_rspace_ekin_density == ON ) then
           call resolve_cutoff_ekin_density(nfout)
        else
          call resolve_cutoff_q(nfout)
        endif
#endif
        if(.not.allocated(nmesh_rs_aug)) allocate(nmesh_rs_aug(natm))
        if(kimg==1) then
           if(.not.allocated(nmesh_rs_aug_h)) allocate(nmesh_rs_aug_h(natm))
        endif
     endif

     ntot = dble(nl*nm*nn)

     if(.not.firstcall)then
       call m_RS_dealloc_hard()
     endif
     nmesh_rs_aug=0 
     if(kimg==1) nmesh_rs_aug_h=0 
     call resolve_atom_centered_mesh(&
     & prealloc=.true.,mode=HARDPART,rcutr=rcut_qr_rs, &
     & id=id,mm=mm,nl=nl,nm=nm,nn=nn,maxcount=nmesh_rs_aug_max,maxcount_h=nmesh_rs_aug_max_h)
     allocate(meshx_rs_aug(nmesh_rs_aug_max,natm));meshx_rs_aug = -1
     allocate(meshy_rs_aug(nmesh_rs_aug_max,natm));meshy_rs_aug = -1
     allocate(meshz_rs_aug(nmesh_rs_aug_max,natm));meshz_rs_aug = -1
     allocate(meshxyz_rs_aug(nmesh_rs_aug_max,natm));meshxyz_rs_aug = -1

     if(kimg==1) then
        allocate(meshxyz_rs_aug_conjg(nmesh_rs_max,natm));meshxyz_rs_aug_conjg = 1
        allocate(meshxyz_rs_aug_h(nmesh_rs_max_h,natm));meshxyz_rs_aug_h = -1
        allocate(map_aug_h(nmesh_rs_aug_max_h,natm));map_aug_h = -1
        allocate(map_aug_h_i(nmesh_rs_aug_max,natm));map_aug_h_i = -1
     endif

     call resolve_atom_centered_mesh(prealloc=.false.,mode=HARDPART,rcutr=rcut_qr_rs,id=id,mm=mm,nl=nl,nm=nm,nn=nn)

     if(iprirs>=2.and.printable)then
         call print_mesh(ntot,HARDPART)
     endif
     firstcall = .false.

     call m_Parallel_init_mpi_rspace_aug(nfout,iprirs,printable,natm,nmesh_rs_aug)

     call tstatc0_end(id_sname)

  end subroutine resolve_mesh_hard0

  subroutine m_RS_resolve_mesh_soft(nfout)
     integer, intent(in) :: nfout
     integer :: nl,nm,nn,mm,id
     real(kind=DP) :: ntot
     integer :: id_sname = -1
     logical,save :: firstcall=.true.
     call tstatc0_begin('m_RS_resolve_mesh_soft ',id_sname,1)
     if (firstcall) then
        call resolve_cutoff_betar(nfout)
        if(.not.allocated(nmesh_rs)) allocate(nmesh_rs(natm))
        if(kimg==1) then
           if(.not.allocated(nmesh_rs_h)) allocate(nmesh_rs_h(natm))
        endif
     endif
     id = fft_box_size_WF(1,0)
     mm = fft_box_size_WF(2,0)
     nl = fft_box_size_WF(1,1)
     nm = fft_box_size_WF(2,1)
     nn = fft_box_size_WF(3,1)
     ntot = dble(nl*nm*nn)

     if(.not.firstcall)then
       call m_RS_dealloc_soft()
     endif
     nmesh_rs=0 
     if(kimg==1) nmesh_rs_h=0 
     call resolve_atom_centered_mesh(prealloc=.true.,mode=SOFTPART,rcutr=rcut_betar_rs,&
        & id=id,mm=mm,nl=nl,nm=nm,nn=nn,maxcount=nmesh_rs_max,maxcount_h=nmesh_rs_max_h)
     allocate(meshx_rs(nmesh_rs_max,natm));meshx_rs = -1
     allocate(meshy_rs(nmesh_rs_max,natm));meshy_rs = -1
     allocate(meshz_rs(nmesh_rs_max,natm));meshz_rs = -1
     allocate(meshxyz_rs(nmesh_rs_max,natm));meshxyz_rs = -1

     if(kimg==1) then
        allocate(meshxyz_rs_conjg(nmesh_rs_max,natm));meshxyz_rs_conjg = 1
        allocate(meshxyz_rs_h(nmesh_rs_max_h,natm));meshxyz_rs_h = -1
        allocate(map_h(nmesh_rs_max_h,natm));map_h = -1
        allocate(map_h_i(nmesh_rs_max,natm));map_h_i = -1
     endif

     call resolve_atom_centered_mesh(prealloc=.false.,mode=SOFTPART,rcutr=rcut_betar_rs,id=id,mm=mm,nl=nl,nm=nm,nn=nn)

     if(iprirs>=2.and.printable)then
         call print_mesh(ntot,SOFTPART)
     endif
     firstcall = .false.

     call tstatc0_end(id_sname)

  end subroutine m_RS_resolve_mesh_soft

  subroutine print_mesh(ntot,mode)
     real(kind=DP), intent(in) :: ntot
     integer :: mode
     integer :: ia,nm
     integer :: totnmesh
     integer,pointer,dimension(:) :: nmesh_rs_p
     integer,pointer,dimension(:,:) :: meshx_rs_p,meshy_rs_p,meshz_rs_p,meshxyz_rs_p
     if(mode==SOFTPART)then
        nmesh_rs_p => nmesh_rs
        meshx_rs_p => meshx_rs
        meshy_rs_p => meshy_rs
        meshz_rs_p => meshz_rs
        meshxyz_rs_p => meshxyz_rs
     else if(mode==HARDPART) then
        nmesh_rs_p => nmesh_rs_aug
        meshx_rs_p => meshx_rs_aug
        meshy_rs_p => meshy_rs_aug
        meshz_rs_p => meshz_rs_aug
        meshxyz_rs_p => meshxyz_rs_aug
     endif
     totnmesh = 0
     do ia=1,natm
        write(nfout,'(a,i0,a,i0)') ' !RS number of mesh points associated to atom ',ia,' : ',nmesh_rs_p(ia)
        if(iprirs>=3)then
          write(nfout,'(a)')              ' !RS associated mesh points ... '
          do nm = 1,nmesh_rs_p(ia)
              write(nfout,'(a,4i8)')      ' !RS ',meshx_rs_p(nm,ia),meshy_rs_p(nm,ia),meshz_rs_p(nm,ia),meshxyz_rs_p(nm,ia)
          enddo
        endif
        totnmesh = totnmesh+nmesh_rs_p(ia)
     enddo
     totnmesh = totnmesh/natm
     write(nfout,'(a,i0,a,i0,a,f10.5)') ' !RS average number of mesh points associated to an atom / total FFT mesh = ' &
     &     ,totnmesh,'/',int(ntot),' = ',dble(totnmesh)/ntot
  end subroutine print_mesh

  subroutine resolve_atom_centered_mesh(prealloc,mode,rcutr,id,mm,nl,nm,nn,maxcount,maxcount_h)
     logical, intent(in) :: prealloc
     integer, intent(in) :: mode
     real(kind=DP), dimension(ntyp), intent(in) :: rcutr
     integer, intent(in) :: id,mm,nl,nm,nn
     integer, optional, intent(out) :: maxcount,maxcount_h
     real(kind=DP) :: rcut_max,rcut,rcut2,rr,rcutmax
     real(kind=DP) :: ex,ey,ez
     integer :: nl_per_atm,nm_per_atm,nn_per_atm
     integer :: nnx,nny,nnz
     integer :: it,ia
     integer :: i1,j1,k1,ii1,jj1,kk1,iil,iim,iin,cl,cm,cn
     integer :: i1min,i1max,j1min,j1max,k1min,k1max
     integer :: icount,icounth
     real(kind=DP) :: fl,fm,fn
     real(kind=DP) :: dx,dy,dz,r2,r,cx,cy,cz
     integer :: maxc,maxch
     integer :: pm
     logical :: smallx,smally,smallz,tof
     real(kind=DP) :: inl,inm,inn,nlhf

     integer,pointer,dimension(:) :: nmesh_rs_p
     integer,pointer,dimension(:) :: nmesh_rs_h_p
     integer,pointer,dimension(:,:) :: meshx_rs_p
     integer,pointer,dimension(:,:) :: meshy_rs_p
     integer,pointer,dimension(:,:) :: meshz_rs_p
     integer,pointer,dimension(:,:) :: meshxyz_rs_p
     integer,pointer,dimension(:,:) :: meshxyz_rs_h_p
     integer,pointer,dimension(:,:) :: map_h_p
     integer,pointer,dimension(:,:) :: map_h_i_p
     integer,pointer,dimension(:,:) :: meshxyz_rs_conjg_p

     if(mode==SOFTPART)then
        nmesh_rs_p => nmesh_rs
        nmesh_rs_h_p => nmesh_rs_h
        meshx_rs_p => meshx_rs
        meshy_rs_p => meshy_rs
        meshz_rs_p => meshz_rs
        meshxyz_rs_p => meshxyz_rs
        meshxyz_rs_h_p => meshxyz_rs_h
        map_h_p => map_h
        map_h_i_p => map_h_i 
        meshxyz_rs_conjg_p => meshxyz_rs_conjg
     else if(mode==HARDPART)then
        nmesh_rs_p => nmesh_rs_aug
        nmesh_rs_h_p => nmesh_rs_aug_h
        meshx_rs_p => meshx_rs_aug
        meshy_rs_p => meshy_rs_aug
        meshz_rs_p => meshz_rs_aug
        meshxyz_rs_p => meshxyz_rs_aug
        meshxyz_rs_h_p => meshxyz_rs_aug_h
        map_h_p => map_aug_h
        map_h_i_p => map_aug_h_i 
        meshxyz_rs_conjg_p => meshxyz_rs_aug_conjg
     endif

     rcut_max = 0.d0
     do it=1,ntyp 
        if(rcut_max<rcutr(it)) rcut_max = rcutr(it)
     enddo

     inl = 1.d0/dble(nl)
     inm = 1.d0/dble(nm)
     inn = 1.d0/dble(nn)

     if(kimg==1)then
       nlhf = id/2
     else
       nlhf = id
     endif

     ex = neix * dsqrt(altv(1,1)**2+altv(1,2)**2+altv(1,3)**2)
     ey = neiy * dsqrt(altv(2,1)**2+altv(2,2)**2+altv(2,3)**2)
     ez = neiz * dsqrt(altv(3,1)**2+altv(3,2)**2+altv(3,3)**2)

     nnx = floor(ex/rcut_max)
     nny = floor(ey/rcut_max)
     nnz = floor(ez/rcut_max)
     smallx = nnx<2
     smally = nny<2
     smallz = nnz<2

     nl_per_atm = floor(dble((2*neix+1)*nl)/dble(nnx))+1
     nm_per_atm = floor(dble((2*neiy+1)*nm)/dble(nny))+1
     nn_per_atm = floor(dble((2*neiz+1)*nn)/dble(nnz))+1
     if(iprirs>=2)then
       write(nfout,'(a,3i8)') ' !RS number of fft elements surrounding an atom  : ', &
       & nl_per_atm,nm_per_atm,nn_per_atm
     endif
     maxc = 0
     if(kimg==1) maxch = 0
     do ia=1,natm
        if(smallx)then
           i1min = -nl*neix
           i1max = +nl*(neix+1)
        else
           cl = -neix*nl+floor((2*neix+1)*nl*dble(neix+pos(ia,1))/dble(2*neix+1))
           i1min = cl-nl_per_atm
           i1max = cl+nl_per_atm
        endif
        if(smally)then
           j1min = -nm*neiy
           j1max = +nm*(neiy+1)
        else
           cm = -neiy*nm+floor((2*neiy+1)*nm*dble(neiy+pos(ia,2))/dble(2*neiy+1))
           j1min = cm-nm_per_atm
           j1max = cm+nm_per_atm
        endif
        if(smallz)then
           k1min = -nn*neiz
           k1max = +nn*(neiz+1)
        else
           cn = -neiz*nn+floor((2*neiz+1)*nn*dble(neiz+pos(ia,3))/dble(2*neiz+1))
           k1min = cn-nn_per_atm
           k1max = cn+nn_per_atm
        endif
        if(iprirs>=2.and. .not. (smallx.or.smally.or.smallz)) &
        &  write(nfout,'(a,i8,a,3i8)') ' !RS center of mesh for atom ',ia,':',cl,cm,cn
        it = ityp(ia)
        rcut2=rcutr(it)*rcutr(it)
        icount = 0
        if(kimg==1) icounth = 0
        do i1=i1min,i1max
        do j1=j1min,j1max
        do k1=k1min,k1max
           iil = i1-1
           iim = j1-1
           iin = k1-1

           ii1 = mod(iil+(2*neix+1)*nl,nl)+1
           jj1 = mod(iim+(2*neiy+1)*nm,nm)+1
           kk1 = mod(iin+(2*neiz+1)*nn,nn)+1

           fl = dble(iil) * inl
           fm = dble(iim) * inm
           fn = dble(iin) * inn

           dx = fl-pos(ia,1)
           dy = fm-pos(ia,2)
           dz = fn-pos(ia,3)

           cx = altv(1,1)*dx+altv(1,2)*dy+altv(1,3)*dz
           cy = altv(2,1)*dx+altv(2,2)*dy+altv(2,3)*dz
           cz = altv(3,1)*dx+altv(3,2)*dy+altv(3,3)*dz
           r2 = cx*cx+cy*cy+cz*cz
           if(r2<rcut2)then
              tof = .false.
              icount = icount+1
              if(kimg==1.and.ii1<=nlhf) then
                 icounth = icounth+1
                 tof = .true.
              endif
              if(.not.prealloc)then
                 meshx_rs_p(icount,ia) = iil
                 meshy_rs_p(icount,ia) = iim
                 meshz_rs_p(icount,ia) = iin
                 if(kimg==1)then
                    if(ii1>nlhf)then
                       ii1 = id - ii1
                       jj1 = nm+2 - jj1
                       kk1 = nn+2 - kk1
                       if(jj1>nm) jj1 = jj1-nm
                       if(jj1.eq.0) jj1=1
                       if(kk1>nn) kk1 = kk1-nn
                       if(kk1.eq.0) kk1=1
                       pm = +1.0d0
                    else
                       pm = -1.d0
                    endif
                 endif
                 meshxyz_rs_p(icount,ia) = nlhf*mm*(kk1-1)+nlhf*(jj1-1)+ii1
                 if(kimg==1) then
                    meshxyz_rs_conjg_p(icount,ia) = pm
                    if(tof) then
                       meshxyz_rs_h_p(icounth,ia) = meshxyz_rs_p(icount,ia)
                       map_h_p(icounth,ia) = icount
                       map_h_i_p(icount,ia) = icounth
                    endif
                 endif
              endif
           endif
        enddo
        enddo
        enddo
        if(prealloc)then
          if(icount>maxc) maxc = icount
          nmesh_rs_p(ia) = icount
          if(kimg==1) then
             if(icounth>maxch) maxch = icounth
             nmesh_rs_h_p(ia) = icounth
          endif
        endif
     enddo
     if(prealloc.and.present(maxcount))then
        maxcount = maxc
     endif
     if(prealloc.and.kimg==1.and.present(maxcount_h))then
        maxcount_h = maxch
     endif
     
  end subroutine resolve_atom_centered_mesh

  subroutine m_RS_dealloc_soft()
     deallocate(meshx_rs)
     deallocate(meshy_rs)
     deallocate(meshz_rs)
     deallocate(meshxyz_rs)
     if(kimg==1) then
        deallocate(meshxyz_rs_h)
        deallocate(map_h)
        deallocate(map_h_i)
        deallocate(meshxyz_rs_conjg)
     endif
  end subroutine m_RS_dealloc_soft

  subroutine m_RS_dealloc_hard()
     deallocate(meshx_rs_aug)
     deallocate(meshy_rs_aug)
     deallocate(meshz_rs_aug)
     deallocate(meshxyz_rs_aug)
     if(kimg==1) then
        deallocate(meshxyz_rs_aug_h)
        deallocate(map_aug_h)
        deallocate(map_aug_h_i)
        deallocate(meshxyz_rs_aug_conjg)
     endif
  end subroutine m_RS_dealloc_hard

  subroutine m_RS_R_minus_pos(pos,ia,nma,inl,inm,inn,cx,cy,cz,rdiff,meshx,meshy,meshz)
     real(kind=DP), dimension(natm,3), intent(in) :: pos
     integer, intent(in) :: ia,nma
     real(kind=DP), intent(in) :: inl,inm,inn
     real(kind=DP), intent(out), dimension(nma) :: cx,cy,cz,rdiff
     integer, dimension(nma,natm), intent(in) :: meshx,meshy,meshz
     real(kind=DP) :: dx,dy,dz
     integer :: ii
     real(kind=DP) :: fl,fm,fn
     real(kind=DP) :: small_shift = 1.d-8
     cx=0.d0;cy=0.d0;cz=0.d0;rdiff=0.d0
     do ii=1,nma
        fl = dble(meshx(ii,ia))*inl
        fm = dble(meshy(ii,ia))*inm
        fn = dble(meshz(ii,ia))*inn
        dx = fl-pos(ia,1)
        dy = fm-pos(ia,2)
        dz = fn-pos(ia,3)
        cx(ii) = altv(1,1)*dx+altv(1,2)*dy+altv(1,3)*dz
        cy(ii) = altv(2,1)*dx+altv(2,2)*dy+altv(2,3)*dz
        cz(ii) = altv(3,1)*dx+altv(3,2)*dy+altv(3,3)*dz
        rdiff(ii) = dsqrt(cx(ii)*cx(ii)+cy(ii)*cy(ii)+cz(ii)*cz(ii))
        if(rdiff(ii).lt.1.d-15)then
           dx = fl-pos(ia,1)+small_shift
           dy = fm-pos(ia,2)+small_shift
           dz = fn-pos(ia,3)+small_shift
           cx(ii) = altv(1,1)*dx+altv(1,2)*dy+altv(1,3)*dz
           cy(ii) = altv(2,1)*dx+altv(2,2)*dy+altv(2,3)*dz
           cz(ii) = altv(3,1)*dx+altv(3,2)*dy+altv(3,3)*dz
           rdiff(ii) = dsqrt(cx(ii)*cx(ii)+cy(ii)*cy(ii)+cz(ii)*cz(ii))
        endif
     enddo
  end subroutine m_RS_R_minus_pos

  subroutine build_qr_clm_ylm1(deriv)
    logical, intent(in), optional :: deriv
    call build_qr_clm_ylm0(deriv,fft_box_size_CD)
  end subroutine build_qr_clm_ylm1

  subroutine build_qr_clm_ylm0(deriv,box_size)
    logical, intent(in), optional :: deriv
    integer, intent(in), dimension(3,0:1) :: box_size
    integer :: ia,it,lmt1,lmt2,np,imesh,iiqitg,lmtp
    integer :: il1,tau1,ilmta1,il2,tau2,ilmta2,l3,ilm3,n
    real(kind=DP) :: f0
    real(kind=DP), allocatable, dimension(:) :: qrsum
    real(kind=DP), allocatable, dimension(:,:) :: dqrsum
    real(kind=DP) :: inl,inm,inn
    real(kind=DP), allocatable, dimension(:) :: cx,cy,cz,rdiff
    real(kind=DP), allocatable, dimension(:) :: ylm
    real(kind=DP), allocatable, dimension(:,:) :: dylm
    real(kind=DP), allocatable, dimension(:) :: q2,qtmp
    real(kind=DP) :: hn,qq,dqq
    integer :: nma
    integer, allocatable :: il3(:)
    logical :: deri
    integer :: lmtcount
    integer :: id_sname=-1
    deri = .true.
    if(present(deriv))then
       deri = deriv
    endif

    call tstatc0_begin('m_RS_build_qr_clm_ylm ',id_sname,1)
    inl = 1.d0/dble(box_size(1,1))
    inm = 1.d0/dble(box_size(2,1))
    inn = 1.d0/dble(box_size(3,1))
    if(allocated(qrsum)) deallocate(qrsum)
    if(allocated(dqrsum)) deallocate(dqrsum)
    if(allocated(plmt1)) deallocate(plmt1)
    if(allocated(plmt2)) deallocate(plmt2)
    if(allocated(nlmtpair)) deallocate(nlmtpair)
    allocate(qrsum(nmesh_rs_aug_max));qrsum=0.d0
    if(deri) then
       allocate(dqrsum(nmesh_rs_aug_max,3));dqrsum=0.d0
    endif
    allocate(cx(nmesh_rs_aug_max));cx=0.d0
    allocate(cy(nmesh_rs_aug_max));cy=0.d0
    allocate(cz(nmesh_rs_aug_max));cz=0.d0
    allocate(rdiff(nmesh_rs_aug_max));rdiff=0.d0
    allocate(ylm(nmesh_rs_aug_max));ylm=0.d0
    if(deri) then
       allocate(dylm(nmesh_rs_aug_max,3));dylm=0.d0
    endif
    allocate(q2(mmesh));q2=0.d0
    allocate(qtmp(mmesh));qtmp=0.d0

    call m_PP_find_maximum_l(n)    ! n-1: maximum l
    n = (n-1) + (n-1) + 1
    allocate(il3(n**2))
    call substitute_il3(n**2,il3) ! -(b_Elec..)
    lmtcount=0
    allocate(nlmtpair(natm));nlmtpair=0
    do ia=1,natm
       it = ityp(ia)
       do lmt1=1,ilmt(it)
          do lmt2 = lmt1,ilmt(it)
             if(il2p(lmt1,lmt2,it)==0) cycle
             lmtcount = lmtcount+1
          enddo
       enddo
       nlmtpair(ia) = lmtcount
       lmtcount=0
    enddo
    allocate(plmt1(maxval(nlmtpair),natm));plmt1=0
    allocate(plmt2(maxval(nlmtpair),natm));plmt2=0
    lmtcount=0
    do ia=1,natm
       it = ityp(ia)
       do lmt1=1,ilmt(it)
          do lmt2 = lmt1,ilmt(it)
             if(il2p(lmt1,lmt2,it)==0) cycle
             lmtcount = lmtcount+1
             plmt1(lmtcount,ia) = lmt1
             plmt2(lmtcount,ia) = lmt2
          enddo
       enddo
       lmtcount=0
    enddo
    if(allocated(qr_clm_ylm)) deallocate(qr_clm_ylm)
    allocate(qr_clm_ylm(nmesh_rs_aug_max,natm,maxval(nlmtpair)));qr_clm_ylm=0.d0
    if(deri)then
       if(allocated(dqr_clm_ylm)) deallocate(dqr_clm_ylm)
       allocate(dqr_clm_ylm(nmesh_rs_aug_max,natm,maxval(nlmtpair),3));dqr_clm_ylm=0.d0
    endif
    do ia=1,natm
       it = ityp(ia)
       nma = nmesh_rs_aug(ia)
       call rmeshs(nmesh(it),nmesh(it),xh(it),rmax(it),radr,hn) ! -(b_PP)
       call m_RS_R_minus_pos(pos,ia,nmesh_rs_aug_max,inl,inm,inn,cx,cy,cz,rdiff,meshx_rs_aug,meshy_rs_aug,meshz_rs_aug)
       do lmtp=1,nlmtpair(ia)
          lmt1 = plmt1(lmtp,ia)
          lmt2 = plmt2(lmtp,ia)
          il1  = ltp( lmt1, it)
          tau1 = taup(lmt1, it)
          ilmta1 = lmta(lmt1,ia)
          il2  = ltp( lmt2, it)
          tau2 = taup(lmt2, it)
          ilmta2 = lmta(lmt2,ia)
          qrsum=0.d0
          if(deri) dqrsum=0.d0
          do np=1,il2p(lmt1,lmt2,it)
             ilm3 = isph(lmt1,lmt2,np,it)
             l3 = il3(ilm3)
             if(lmt2>=lmt1) then
                iiqitg = iqitg(il1,tau1,il2,tau2,l3+1,it)
             else
                iiqitg = iqitg(il2,tau2,il1,tau1,l3+1,it)
             end if
             if(iiqitg == 0) cycle
             do imesh=1,nmesh(it)
                qtmp(imesh) = qrspspw(imesh,iiqitg)/radr(imesh)/radr(imesh)
             enddo
             call sphr_general(nma,ilm3,cx(1:nma),cy(1:nma),cz(1:nma),ylm(1:nma))
             if(deri) then
                call sphr_diff_general(nma,nma,ilm3,&
                     &                 cx(1:nma),cy(1:nma),cz(1:nma),dylm(1:nma,1:3))
             endif
             call init_cubic_spline(nmesh(it),radr(1:nmesh(it)),&
                  & qtmp(1:nmesh(it)),q2(1:nmesh(it)))
             f0 = iwei(ia) * dl2p(lmt1,lmt2,np,it)
             do imesh=1,nma
                call cubic_spline(nmesh(it),radr,qtmp(1:nmesh(it)),&
                     & q2(1:nmesh(it)),rdiff(imesh),qq,dqq)
                qrsum(imesh) = qrsum(imesh)+univol*f0*ylm(imesh)*qq
                if(deri)then
                  dqrsum(imesh,1) = dqrsum(imesh,1) &
                & -univol*f0*(ylm(imesh)*(dqq*cx(imesh)/rdiff(imesh))+qq*dylm(imesh,1))
                  dqrsum(imesh,2) = dqrsum(imesh,2) &
                & -univol*f0*(ylm(imesh)*(dqq*cy(imesh)/rdiff(imesh))+qq*dylm(imesh,2))
                  dqrsum(imesh,3) = dqrsum(imesh,3) &
                & -univol*f0*(ylm(imesh)*(dqq*cz(imesh)/rdiff(imesh))+qq*dylm(imesh,3))
                endif
             enddo
          enddo
          qr_clm_ylm(1:nma,ia,lmtp) = qrsum(1:nma)
          if(deri) dqr_clm_ylm(1:nma,ia,lmtp,1:3) = dqrsum(1:nma,1:3)
       enddo
    enddo
    deallocate(il3)
    deallocate(qrsum)
    if(deri) deallocate(dqrsum)
    deallocate(cx,cy,cz,rdiff,ylm,q2,qtmp)
    if(deri) deallocate(dylm)

    call tstatc0_end(id_sname)

  end subroutine build_qr_clm_ylm0

end module m_Realspace

