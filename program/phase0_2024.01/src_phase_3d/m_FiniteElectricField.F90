#define FEF_ZGEMM
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 593 $)
!
!  MODULE: m_FiniteElectricField
!
!  AUTHOR(S): T. Yamamoto   Oct/01/2007
!
!  Contact address :  IIS,The University of Tokyo RSS21 project
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
module m_FiniteElectricField
! $ID: $
  use m_Const_Parameters,     only : DP, PAI, PAI2, PAI4, ON, OFF, CARTS, BUCS &
                                 & , MONKHORST_PACK, PARA, EXECUT, SKIP, DELTA
  use m_Control_Parameters,   only : iprifef, nspin, printable, elec_field &
                                 & , kimg, way_ksample, sw_check_polar, af, neg, ipriparallel
  use m_Files,                only : nfout
  use m_Crystal_Structure,    only : altv, rltv, univol, nopr, op, imag
  use m_PlaneWaveBasisSet,    only : kg1, ngabc, nbase, iba, kgp
  use m_Electronic_Structure, only : totch, zaj_l, m_ES_add_it_to_vnlph, nrvf_ordr, neordr &
                                 & , fsr_l,fsi_l, occup_l
  use m_Kpoints,              only : kv3, k_symmetry, vkxyz, mp_index, kshift
  use m_Ionic_System,         only : ntyp,natm,natm2,ityp,cps,iwei,pos
  use m_PseudoPotential,      only : ival,modnrm,dk_fef,qitg_fef &
                                 & , m_PP_include_vanderbilt_pot &
                                 & , n_non0_lmtxlmt,index_lmt1_lmt2,lmta,nlmt,ilmt, nlmta &
                                 & , ltp,taup,il2p,isph,iqitg,dl2p &
                                 & , m_PP_find_maximum_l, nlmtt
  use m_NonLocal_Potential,   only : snl
  use m_Parallelization,      only : npes,mype,map_ek,map_z, mpi_e_world, map_z, map_e, MPI_CommGroup &
                                 & , map_k, myrank_k, myrank_e, mpi_ge_world, nrank_k, nrank_g, myrank_g, np_nvale
  use m_Parallelization,      only : mp_e, myrank_nval, mpi_ke_world, mpi_ge_world, mpi_kg_world &
                                 & , ista_g1k, iend_g1k, ista_fs, iend_fs, np_g1k
  use m_Parallelization,      only : mp_e, myrank_nvale, mpi_k_world, ista_nvale, iend_nvale &
                                 & , m_Parallel_init_mpi_nvale, ista_snl, iend_snl
  use m_Timing,               only : tstatc0_begin, tstatc0_end
  use mpi
  implicit none

  integer :: ns ! the number of valence states
  integer :: msize ! maxval(mp_index)
  integer :: numef ! the number of conponents of the electric field
  integer, allocatable :: kgbpp(:,:,:,:) !d(0:n1-1,0:n2-1,0:n3-1,numef)
  integer, allocatable :: kgbpm(:,:,:,:) !d(0:n1-1,0:n2-1,0:n3-1,numef)
  integer, allocatable :: ista_kgbpp(:,:,:,:) !d(0:n1-1,0:n2-1,0:n3-1,numef)
  integer, allocatable :: iend_kgbpp(:,:,:,:) !d(0:n1-1,0:n2-1,0:n3-1,numef)
  integer, allocatable :: ista_kgbpm(:,:,:,:) !d(0:n1-1,0:n2-1,0:n3-1,numef)
  integer, allocatable :: iend_kgbpm(:,:,:,:) !d(0:n1-1,0:n2-1,0:n3-1,numef)
  integer, allocatable :: indgbp(:,:,:,:,:,:) !d(kg1,2,0:n1-1,0:n2-1,0:n3-1,numef)
  integer, allocatable :: indgbm(:,:,:,:,:,:) !d(kg1,2,0:n1-1,0:n2-1,0:n3-1,numef)
  integer :: elec_id(3)
  real(kind=DP) :: efa(3)
  complex(kind=DP), allocatable :: smat(:,:) !d(ns,ns)
  real(kind=DP), allocatable :: berry(:,:,:) !d(0:msize-1,0:msize-1,numef)
  complex(kind=DP), allocatable :: berry_det(:,:,:,:) !d(0:n1-1,0:n2-1,0:n3-1,numef)
  complex(kind=DP), allocatable, dimension(:,:,:) :: berry_dett

  integer :: nopr_trs
  integer, allocatable :: op_trs(:,:,:) !d(3,3,nopr_trs)

  integer, allocatable :: map_nnn2ik(:,:,:) !d(0:n1-1,0:n2-1,0:n3-1)
  integer, allocatable :: map_nnn2op(:,:,:) !d(0:n1-1,0:n2-1,0:n3-1)
  integer, allocatable :: map_nnn2gs(:,:,:,:) !d(3,0:n1-1,0:n2-1,0:n3-1)
  integer, allocatable :: map_ik2nnn(:,:,:) !d(3,nopr_trs,kv3) = {[n1,n2,n3]}
  logical, allocatable :: flag_op(:,:) !d(nopr_trs,kv3)

  integer :: ngmax
  integer, allocatable :: ngpt(:,:) !d(ngmax,nopr_trs)

  real(kind=DP), allocatable :: grad_fef(:,:,:,:) !d(kg1,ns,kv3,kimg)

  complex(kind=DP), private, allocatable :: ftqe(:,:,:,:) !d(nlmt,nlmt,natm,numef)

  real(kind=DP), dimension(3) :: pmac,pmac_old

  real(kind=DP), allocatable :: fsr_fef(:,:,:) ! d(nval,nlmta,kv3)
  real(kind=DP), allocatable :: fsi_fef(:,:,:) ! d(nval,nlmta,kv3)
  real(kind=DP), allocatable :: neordr_val(:,:) ! d(nval,kv3)
  real(kind=DP), allocatable :: wfv(:,:,:,:) ! d(kg1,nval,kv3,kimg)
  real(kind=DP), allocatable :: snlt(:,:,:)

  complex(kind=DP) :: cone=dcmplx(1.d0,0.d0), czero=dcmplx(0.d0,0.d0)

!  include 'mpif.h'                                      ! MPI

contains

  function modp(k,m) result(ks)
    integer :: ks
    integer, intent(in) :: k,m
    if(k<0) then
       ks = k + m*(1-int(dble(k)/m-1.d-8))
    else
       ks = k - m*int(dble(k)/m+1.d-8)
    end if
  end function modp

  subroutine m_FEF_init(nfout)
    integer, intent(in) :: nfout
    integer :: ig,ik
    integer :: jj1,jj2,i1,i2,i3,j,jp,ik1,iopr1,ik2,iopr2
    integer :: ipm,ii(3)
    integer :: iopr,jg,ng(3)
    integer :: irot(3,3)
    integer :: m1,m2,m3
    integer :: igs(3)
    integer :: iwork,nn
    integer :: ierr

    call check_elec_field

    call check_kpoints

    write(nfout,*) 'efa=',efa(1:3)
    write(nfout,*) 'elec_id=',elec_id

    call alloc_ngpt !-> ngmax

    do iopr=1,nopr_trs
       irot(:,:) = nint(matmul(transpose(altv),matmul(op_trs(:,:,iopr),rltv))/PAI2)
       write(nfout,'("irot=",9(1x,i3))') irot(1:3,1:3)
       IG_LOOP: do ig=1,ngmax
          ng(1:3) = matmul(irot,ngabc(ig,1:3))
          do jg=1,kgp
             if(ngabc(jg,1)==ng(1).and.ngabc(jg,2)==ng(2).and.ngabc(jg,3)==ng(3)) then
                ngpt(ig,iopr) = jg
                cycle IG_LOOP
             end if
          end do
          write(nfout,*) 'Error: NOT found the rotated G'
          call phase_error_with_msg(nfout,'Error: NOT found the rotated G',__LINE__,__FILE__)
       end do IG_LOOP
    end do

    call alloc_indgb

    m1 = mp_index(1)-1
    m2 = mp_index(2)-1
    m3 = mp_index(3)-1
    do ig=1,3
       if(elec_id(ig)==0) cycle
       do i1=0,mp_index(1)-1
          do i2=0,mp_index(2)-1
             do i3=0,mp_index(3)-1
                write(nfout,*) 'i1,i2,i3,ig=',i1,i2,i3,ig
                ik1 = map_nnn2ik(i1,i2,i3)
                iopr1 = map_nnn2op(i1,i2,i3)
                ii(1) = i1
                ii(2) = i2
                ii(3) = i3
                jp = ii(ig)+1
                ipm = 0
                if(jp>=mp_index(ig)) then
                   ipm = +1
                   jp = 0
                end if
                ii(ig) = jp
                ik2 = map_nnn2ik(ii(1),ii(2),ii(3))
                iopr2 = map_nnn2op(ii(1),ii(2),ii(3))
                igs(1:3) = map_nnn2gs(1:3,ii(1),ii(2),ii(3)) - map_nnn2gs(1:3,i1,i2,i3)
                igs(ig) = igs(ig) + ipm
                write(nfout,*) 'igs=',igs
                call get_ind_of_g_plus_b(igs,i1,i2,i3,ig,elec_id(ig),ik1,iopr1,ik2,iopr2 &
                                      & ,m1,m2,m3,indgbp,kgbpp)

                ii(1) = i1
                ii(2) = i2
                ii(3) = i3
                jp = ii(ig)-1
                ipm = 0
                if(jp<0) then
                   ipm = -1
                   jp = mp_index(ig)-1
                end if
                ii(ig) = jp
                ik2 = map_nnn2ik(ii(1),ii(2),ii(3))
                iopr2 = map_nnn2op(ii(1),ii(2),ii(3))
                igs(1:3) = map_nnn2gs(1:3,ii(1),ii(2),ii(3)) - map_nnn2gs(1:3,i1,i2,i3)
                igs(ig) = igs(ig) + ipm
                write(nfout,*) 'igs=',igs
                call get_ind_of_g_plus_b(igs,i1,i2,i3,ig,elec_id(ig),ik1,iopr1,ik2,iopr2 &
                                      & ,m1,m2,m3,indgbm,kgbpm)
             end do
          end do
       end do
    end do
    call dealloc_ngpt

    do ig=1,3
       if(elec_id(ig)==0) cycle
       do i1=0,mp_index(1)-1
          do i2=0,mp_index(2)-1
             do i3=0,mp_index(3)-1
                nn = kgbpp(i1,i2,i3,elec_id(ig))
                iwork = (nn-1)/nrank_g + 1
                ista_kgbpp(i1,i2,i3,elec_id(ig)) = min(myrank_g*iwork+1, nn+1)
                iend_kgbpp(i1,i2,i3,elec_id(ig)) = min(ista_kgbpp(i1,i2,i3,elec_id(ig))+iwork-1, nn)

                nn = kgbpm(i1,i2,i3,elec_id(ig))
                iwork = (nn-1)/nrank_g + 1
                ista_kgbpm(i1,i2,i3,elec_id(ig)) = min(myrank_g*iwork+1, nn+1)
                iend_kgbpm(i1,i2,i3,elec_id(ig)) = min(ista_kgbpm(i1,i2,i3,elec_id(ig))+iwork-1, nn)
             enddo
          enddo
       enddo
    enddo

    pmac = 0.d0;pmac_old = 0.d0
  end subroutine m_FEF_init

  subroutine check_elec_field
    integer :: ig
    real(kind=8), parameter :: eps = 1.d-12

    elec_id(:) = 0
    numef = 0
    do ig=1,3
       efa(ig) =  dot_product(elec_field(:),altv(:,ig))
       if(abs(efa(ig))>eps) then
          numef = numef + 1
          elec_id(ig) = numef
       end if
    end do
  end subroutine check_elec_field

  subroutine check_kpoints
    integer :: i,j,jj,ik,jk
    integer :: i1,i2,i3,jj1,jj2
    integer :: iopr
    integer :: n1,n2,n3
    integer :: kn(3),kns(3)
    real(kind=8) :: rotk(3),dn,tmpk(3)
    logical :: flag_trs

    if(way_ksample /= MONKHORST_PACK) then
       write(nfout,*) 'set way_ksample = MONKHORST_PACK.'
       call phase_error_with_msg(nfout,'way_ksample /= MONKHORST_PACK.',__LINE__,__FILE__)
    end if
    do i=1,3
       if(abs(kshift(i))>1.d-12) then
          write(nfout,*) 'set kshift = 0.'
          call phase_error_with_msg(nfout,'kshit > 0',__LINE__,__FILE__)
       end if
    end do

    if(imag==PARA) then
       flag_trs = .true.
       do i=1,nopr
          if(abs(op(1,1,i)+1)<1.d-12.and. &
           & abs(op(2,2,i)+1)<1.d-12.and. &
           & abs(op(3,3,i)+1)<1.d-12) then
             flag_trs = .false.
             exit
          end if
       end do
    else
       flag_trs = .false.
    end if
    if(flag_trs) then
       nopr_trs = 2*nopr
       allocate(op_trs(3,3,nopr_trs))
       op_trs(1:3,1:3,1:nopr) = nint(op(1:3,1:3,1:nopr))
       op_trs(1:3,1:3,nopr+1:nopr_trs) = -nint(op(1:3,1:3,1:nopr))
    else
       nopr_trs = nopr
       allocate(op_trs(3,3,nopr_trs))
       op_trs = nint(op)
    end if
    write(nfout,*) 'nopr=',nopr
    write(nfout,*) 'nopr_trs=',nopr_trs

    msize = maxval(mp_index)
    n1 = mp_index(1)
    n2 = mp_index(2)
    n3 = mp_index(3)
    allocate(map_nnn2ik(0:n1-1,0:n2-1,0:n3-1)); map_nnn2ik=0
    allocate(map_nnn2op(0:n1-1,0:n2-1,0:n3-1)); map_nnn2op=0
    allocate(map_nnn2gs(3,0:n1-1,0:n2-1,0:n3-1)); map_nnn2gs=0
    allocate(map_ik2nnn(3,nopr_trs,kv3)); map_ik2nnn=-1
    allocate(flag_op(nopr_trs,kv3)); flag_op=.false.

    write(nfout,'("mp_index=",3i6)') mp_index(1:3)
!!$    write(nfout,'("i,mp_index=",2i6)') i,mp_index(i)
    do ik=1,kv3
       LOOP1: do iopr=1,nopr_trs
          tmpk = matmul(op_trs(:,:,iopr),vkxyz(ik,:,CARTS))
          rotk = matmul(transpose(altv),matmul(op_trs(:,:,iopr),vkxyz(ik,:,CARTS)))/PAI2
          write(nfout,'("rotk=",3f20.10)') rotk
          do j=1,3
             dn = mp_index(j)*rotk(j)
             kn(j) = nint(dn)
             if(abs(dn-kn(j)) > 1.d-8) cycle LOOP1
             if(kn(j)<0) then
                kns(j) = kn(j)+mp_index(j)*(1-int(dble(kn(j))/mp_index(j)-1.d-8))
             else
                kns(j) = kn(j)-mp_index(j)*int(dble(kn(j))/mp_index(j)+1.d-8)
             end if
          end do
          write(nfout,'("kn=",3(1x,i4)," kns=",3(1x,i4)," iopr=",i4)') kn(1:3),kns(1:3),iopr
          if(map_nnn2ik(kns(1),kns(2),kns(3))>0) cycle LOOP1
          map_nnn2ik(kns(1),kns(2),kns(3)) = ik
          map_nnn2op(kns(1),kns(2),kns(3)) = iopr
          map_nnn2gs(1:3,kns(1),kns(2),kns(3)) = (kns(1:3) - kn(1:3))/mp_index(1:3)
          map_ik2nnn(1:3,iopr,ik) = kns(1:3)
          flag_op(iopr,ik) = .true.
       end do LOOP1
    end do

    ! Write map func.
    do ik=1,kv3
       do iopr=1,nopr_trs
          if(flag_op(iopr,ik)) then
             write(nfout,'("ik=",i3," iopr=",i3," n1=",i3," n2=",i3," n3=",i3)') ik,iopr,map_ik2nnn(1:3,iopr,ik)
          end if
       end do
    end do
    do i1=0,mp_index(1)-1
       do i2=0,mp_index(2)-1
          do i3=0,mp_index(3)-1
             ik=map_nnn2ik(i1,i2,i3)
             iopr=map_nnn2op(i1,i2,i3)
             write(nfout,'(" n1=",i3," n2=",i3," n3=",i3,"ik=",i3," iopr=",i3," gs=",3(1x,i3))') &
            &  i1,i2,i3,ik,iopr,map_nnn2gs(1:3,i1,i2,i3)
          end do
       end do
    end do

  end subroutine check_kpoints

  subroutine alloc_ngpt
    integer :: ik
    ngmax = nbase(iba(1),1)
    do ik=2,kv3
       ngmax = max(ngmax,nbase(iba(ik),ik))
    end do
    write(nfout,*) 'ngmax=',ngmax
    write(nfout,*) 'nopr_trs=',nopr_trs
    allocate(ngpt(ngmax,nopr_trs))
  end subroutine alloc_ngpt

  subroutine dealloc_ngpt
    deallocate(ngpt)
  end subroutine dealloc_ngpt

  subroutine alloc_indgb
    integer :: m1,m2,m3
    m1 = mp_index(1)-1
    m2 = mp_index(2)-1
    m3 = mp_index(3)-1
    allocate(indgbp(kg1,2,0:m1,0:m2,0:m3,numef)); indgbp = 0
    allocate(indgbm(kg1,2,0:m1,0:m2,0:m3,numef)); indgbm = 0
    allocate(kgbpp(0:m1,0:m2,0:m3,numef))
    allocate(kgbpm(0:m1,0:m2,0:m3,numef))
    allocate(berry_det(0:m1,0:m2,0:m3,numef))
    allocate(berry(0:msize-1,0:msize-1,numef))
    allocate(ista_kgbpp(0:m1,0:m2,0:m3,numef))
    allocate(iend_kgbpp(0:m1,0:m2,0:m3,numef))
    allocate(ista_kgbpm(0:m1,0:m2,0:m3,numef))
    allocate(iend_kgbpm(0:m1,0:m2,0:m3,numef))
  end subroutine alloc_indgb

  subroutine alloc_smat
    logical, save :: first = .true.
    if(first) then
       ns = int(totch + 1.d-13)/2
       !!$write(nfout,*) 'debug FEF: ns=',ns
       write(nfout,*) 'm_FEF detected the number of valence bands; ns=',ns
       first = .false.
       call alloc_grad
    end if
    allocate(smat(ns,ns));smat = 0.d0
    grad_fef = 0.d0
  end subroutine alloc_smat

  subroutine dealloc_smat
    deallocate(smat)
  end subroutine dealloc_smat

  subroutine get_ind_of_g_plus_b(igs,n1,n2,n3,ig,id,ik1,iopr1,ik2,iopr2,m1,m2,m3,indgb,kgbp)
    implicit none
    integer, intent(in) :: igs(3),n1,n2,n3,ig,id,ik1,iopr1,ik2,iopr2,m1,m2,m3
    integer, intent(inout) :: indgb(kg1,2,0:m1,0:m2,0:m3,numef)
    integer, intent(inout) :: kgbp(0:m1,0:m2,0:m3,numef)

    ! local variables
    integer :: jj,i,j,k,l
    integer :: iga,igb,igc,jga,jgb,jgc

    jj = 0
    loop_i: do i=1,iba(ik1)
       k=ngpt(nbase(i,ik1),iopr1)
       iga = ngabc(k,1)+igs(1)
       igb = ngabc(k,2)+igs(2)
       igc = ngabc(k,3)+igs(3)
       loop_j: do j=1,iba(ik2)
          l=ngpt(nbase(j,ik2),iopr2)
          jga = ngabc(l,1)
          jgb = ngabc(l,2)
          jgc = ngabc(l,3)
          if(iga == jga .and. igb == jgb .and. igc == jgc) then
             jj = jj +1
             indgb(jj,2,n1,n2,n3,id) = i
             indgb(jj,1,n1,n2,n3,id) = j
             exit loop_j
          end if
       end do loop_j
    end do loop_i
    kgbp(n1,n2,n3,id) = jj
    write(nfout,'(a,5i8)') 'ind_of_g_plus_b:ik1,ik2,id,kgbp,kg1=',ik1,ik2,id,jj,kg1
  end subroutine get_ind_of_g_plus_b

  subroutine calc_overlap(n1,n2,n3,ik0,iopr0,gs0,ik1,iopr1,gs1,ig,id)
    integer, intent(in) :: n1,n2,n3,ik0,iopr0,ik1,iopr1,ig,id
    integer, intent(in) :: gs0(3),gs1(3)
    integer :: n,m,i,ii,jj,mm,nn
    integer :: ipm0,ipm1
    integer :: ip,it,ia,lmt1,lmt2,u,v,mdvdb
    real(kind=8) :: fnr,fni,fmr,fmi,sr,si,ph,dk(3),vk0(3),vk1(3)
    complex(kind=8) :: defchg,fac,fs1u,fs1v,fs2u,fs2v
    complex(kind=8) :: expkt(natm),csum
    complex(kind=8),allocatable,dimension(:) :: smatt
    integer :: ierr,id_sname=-1,id_sname2=-1
    call tstatc0_begin('calc_overlap  ',id_sname,1)

    if(iopr0>nopr) then
       ipm0=-1
    else
       ipm0=1
    end if
    if(iopr1>nopr) then
       ipm1=-1
    else
       ipm1=1
    end if
    allocate(smatt(ns));smatt=0.d0
    do m=1,ns
       mm = neordr_val(m,ik1)
       smatt = 0.d0
       do n=ista_nvale,iend_nvale
          nn = neordr_val(n,ik0)
          sr = 0.d0
          si = 0.d0
          do i=ista_kgbpp(n1,n2,n3,id),iend_kgbpp(n1,n2,n3,id)
             ii=indgbp(i,2,n1,n2,n3,id)
             jj=indgbp(i,1,n1,n2,n3,id)
             fnr = wfv(ii,nn,ik0,1)
             fni = wfv(ii,nn,ik0,2) * ipm0
             fmr = wfv(jj,mm,ik1,1)
             fmi = wfv(jj,mm,ik1,2) * ipm1
             sr = sr + fnr*fmr + fni*fmi
             si = si + fnr*fmi - fni*fmr
          end do
          call mpi_allreduce(mpi_in_place,sr,1,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
          call mpi_allreduce(mpi_in_place,si,1,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
          !!d$write(nfout,*) 'n,m,sr,si=',n,m,sr,si
          !smat(n,m) = dcmplx(sr,si)
          smatt(n) = dcmplx(sr,si)
       end do
       call mpi_allreduce(mpi_in_place,smatt,ns,mpi_double_complex,mpi_sum,mpi_kg_world,ierr)
       smat(:,m) = smatt(:)
    end do

    if(modnrm /= EXECUT) return
    do ia=1,natm
       ph = PAI2*dot_product(gs0-gs1,pos(ia,1:3))
       if(iopr0>1.and.iopr0<=nopr) then
          vk0 = vkxyz(ik0,1:3,CARTS)
          dk = matmul(op_trs(:,:,iopr0),vk0) - vk0
          ph = ph + dot_product(dk,cps(ia,1:3))
       else if(iopr0>nopr+1) then
          vk0 = vkxyz(ik0,1:3,CARTS)
          dk = vk0 + matmul(op_trs(:,:,iopr0),vk0)
          ph = ph + dot_product(dk,cps(ia,1:3))
       end if
       if(iopr1>1.and.iopr1<=nopr) then
          vk1 = vkxyz(ik1,1:3,CARTS)
          dk = vk1 - matmul(op_trs(:,:,iopr1),vk1)
          ph = ph + dot_product(dk,cps(ia,1:3))
       else if(iopr1>nopr+1) then
          vk1 = vkxyz(ik1,1:3,CARTS)
          dk = -(vk1 + matmul(op_trs(:,:,iopr1),vk1))
          ph = ph + dot_product(dk,cps(ia,1:3))
       end if
       expkt(ia) = dcmplx(cos(ph),sin(ph))
    end do

    do m=1,ns
       mm = neordr_val(m,ik1)
       smatt(:) = 0.d0
       do n=ista_nvale,iend_nvale
          nn = neordr_val(n,ik0)
          defchg = dcmplx(0.d0,0.d0)
          do ia=1,natm
             it=ityp(ia)
             mdvdb = m_PP_include_vanderbilt_pot(it)
             if(mdvdb == SKIP) cycle
             csum = dcmplx(0.d0,0.d0)
             do lmt1=1,ilmt(it)
                do lmt2=1,ilmt(it)
                   u = lmta(lmt1,ia)
                   v = lmta(lmt2,ia)
                   fac   = ftqe(lmt1,lmt2,it,id)*iwei(ia)
                   fs1u  = dcmplx(fsr_fef(nn,u,ik0),-ipm0*fsi_fef(nn,u,ik0))
                   fs2v  = dcmplx(fsr_fef(mm,v,ik1),ipm1*fsi_fef(mm,v,ik1))
                   csum = csum + fac*fs1u*fs2v
                end do
             end do
             defchg = defchg +  csum * expkt(ia) * iwei(ia)
          end do
          !!smat(n,m) = smat(n,m) + defchg
          smatt(n) = defchg
       end do
       call mpi_allreduce(mpi_in_place,smatt,ns,mpi_double_complex,mpi_sum,mpi_kg_world,ierr)
       smat(:,m) = smat(:,m)+smatt(:)
    end do
    call tstatc0_end(id_sname)
  end subroutine calc_overlap

  subroutine calc_berry_phase(n1,n2,n3,id)
    integer, intent(in) :: n1,n2,n3,id

    integer :: i, info, ipiv(ns), lwork, ip, ierr
    real(kind=8) :: ph
    complex(kind=DP) :: p
    complex(kind=DP), allocatable :: work(:)
    integer :: j

    ! LU-facterize the overlap matrix
    call zgetf2(ns,ns,smat,ns,ipiv,info)
    if(info/=0) then
       write(nfout,*) 'Error(ZGETF2): INFO=',info
       call phase_error_with_msg(nfout,'Error(ZGETF2)',__LINE__,__FILE__)
    end if

    ! Get the determinant of the overlap matrix
    p = dcmplx(1.d0,0.d0)
!!    do i=1,ns
    do i=ista_nvale,iend_nvale
       p = p * smat(i,i)
    end do

    call mpi_allreduce(mpi_in_place,p,1,mpi_double_complex,mpi_prod,mpi_kg_world,ierr)
!    berry_det(n1,n2,n3,id) = p
    berry_dett(n1,n2,n3) = p
!debug
!!    write(nfout,*) 'n1,n2,n3,id,berry_det=',n1,n2,n3,id,p
!end debug

    ! Inverte the overlap matrix
    lwork = -1
    allocate(work(ns))
    call zgetri(ns,smat,ns,ipiv,work,lwork,info)
    lwork = work(1) ! optimal size of the work array
    deallocate(work)
    allocate(work(lwork))
    call zgetri(ns,smat,ns,ipiv,work,lwork,info)
    if(info/=0) then
       write(nfout,*) 'Error(ZGETRI): INFO=',info
       call phase_error_with_msg(nfout,'Error(ZGETRI)',__LINE__,__FILE__)
    end if
    deallocate(work)
  end subroutine calc_berry_phase

#ifdef FEF_ZGEMM
  subroutine  calc_grad(n1,n2,n3,ik0,iopr0,gs0,m1,m2,m3,ik1,iopr1,gs1,ig,id)
    integer, intent(in) :: n1,n2,n3,ik0,iopr0,gs0(3)
    integer, intent(in) :: m1,m2,m3,ik1,iopr1,gs1(3)
    integer, intent(in) :: ig,id

    integer :: m, n, i, ii, jj, mm, nn
    integer :: ipm0,ipm1
    integer :: ia,it,nb,mdvdb
    integer :: lmt1,lmt2,v,l1
    real(kind=DP) :: fac, sr, si, fmr, fmi
    real(kind=DP), allocatable, dimension(:,:) :: grad,gsmat
    real(kind=DP), allocatable, dimension(:,:,:) :: gsmatb
    !real(kind=DP) :: grad(kg1,kimg),gsmat(kg1,kimg)
    real(kind=DP) :: ph,ph0,ph1
    real(kind=DP) :: dk(3),vk0(3),vk1(3)
    complex(kind=DP) :: exp0,exp1,qf,ctmp
    complex(kind=DP), allocatable, dimension(:) :: expgr
    complex(kind=DP), parameter :: zi=(0.d0,1.d0)
    complex(kind=DP), allocatable, dimension(:,:) :: qfmat,snlmat,cmat
    integer :: ierr,icount,ibcount
    integer :: id_sname = -1
    call tstatc0_begin('calc_grad ',id_sname,1)

    fac = efa(ig)/PAI4*mp_index(ig)
    allocate(grad(kg1,kimg));grad = 0.d0
    if(modnrm == EXECUT) then
      allocate(gsmat(kg1,kimg));gsmat = 0.d0
      allocate(gsmatb(kg1,kimg,np_nvale));gsmatb = 0.d0
      allocate(expgr(maxval(np_g1k)));expgr = dcmplx(0.d0,0.d0)
    endif
    if(iopr0>nopr) then
       ipm0=-1
    else
       ipm0=1
    end if
    if(iopr1>nopr) then
       ipm1=-1
    else
       ipm1=1
    end if

    if(iopr0==1) then
       if(modnrm == EXECUT)then
          allocate(snlmat(np_g1k(ik0),nlmta));snlmat=dcmplx(0.d0,0.d0)
          allocate(qfmat(nlmta,np_nvale));qfmat=dcmplx(0.d0,0.d0)
          allocate(cmat(np_g1k(ik0),np_nvale));cmat=dcmplx(0.d0,0.d0)
          icount = 0
          do ia=1,natm
             it=ityp(ia)
             mdvdb = m_PP_include_vanderbilt_pot(it)
             if(mdvdb == SKIP) cycle
             do i=ista_g1k(ik0),iend_g1k(ik0)
                nb = nbase(i,ik0)
                ph = -PAI2*dot_product(ngabc(nb,1:3),pos(ia,1:3))
                expgr(i-ista_g1k(ik0)+1) = dcmplx(cos(ph),sin(ph))
             end do
             do lmt1=1,ilmt(it)
                icount = icount+1
                do i=1,np_g1k(ik0)
                   snlmat(i,icount) = expgr(i)*snlt(i,lmt1,(ik0-1)/nspin+1)
                end do
             end do
          enddo
          ibcount = 0
          do m=ista_nvale,iend_nvale
             mm = neordr_val(m,ik1)
             ibcount = ibcount+1
             icount = 0
             do ia=1,natm
                it=ityp(ia)
                mdvdb = m_PP_include_vanderbilt_pot(it)
                if(mdvdb == SKIP) cycle
                ph1 = -PAI2*dot_product(gs1,pos(ia,1:3))
                if(iopr1>1.and.iopr1<=nopr) then
                   vk1 = vkxyz(ik1,1:3,CARTS)
                   dk = vk1 - matmul(op_trs(:,:,iopr1),vk1)
                   ph1 = ph1 + dot_product(dk,cps(ia,1:3))
                else if(iopr1>nopr+1) then
                   vk1 = vkxyz(ik1,1:3,CARTS)
                   dk = -(vk1 + matmul(op_trs(:,:,iopr1),vk1))
                   ph1 = ph1 + dot_product(dk,cps(ia,1:3))
                end if
                exp1 = exp(dcmplx(0.d0,ph1))
                do lmt1=1,ilmt(it)
                   l1 = ltp(lmt1,it)

                   qf = dcmplx(0.d0,0.d0)
                   do lmt2=1,ilmt(it)
                      v = lmta(lmt2,ia)
                      qf = qf + ftqe(lmt1,lmt2,it,id)*iwei(ia) &
                          & * dcmplx(fsr_fef(mm,v,ik1),ipm1*fsi_fef(mm,v,ik1))*exp1
                   end do
                   qf = qf * zi**(-l1)
                   icount = icount+1
                   qfmat(icount,ibcount) = qf
                end do
             end do
          enddo
          call zgemm('N','N',np_g1k(ik0),ibcount,icount,cone,snlmat,np_g1k(ik0),qfmat,nlmta,czero,cmat,np_g1k(ik0))
          !do i=1,np_g1k(ik0)
          !   do m=1,ibcount
          !      do lmt1=1,icount
          !         cmat(i,m) = cmat(i,m) + snlmat(i,lmt1)*qfmat(lmt1,m)
          !      enddo
          !   enddo
          !enddo
       endif
       if(modnrm == EXECUT)then
          ibcount = 0
          gsmatb = 0.d0
          do m=ista_nvale,iend_nvale
             ibcount = ibcount+1
             gsmat = 0.d0
             do i=ista_g1k(ik0),iend_g1k(ik0)
                ctmp = cmat(i-ista_g1k(ik0)+1,ibcount)
                gsmat(i,1) = gsmat(i,1) + dble(ctmp)
                gsmat(i,2) = gsmat(i,2) + dimag(ctmp)
             end do
             call mpi_allreduce(mpi_in_place,gsmat,kg1*kimg,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
             gsmatb(:,1:2,ibcount) = gsmat(:,1:2)
          enddo
       endif
       do n=1,ns
          nn = neordr_val(n,ik0)
          grad = 0.d0
          !!do m=1,ns
          ibcount = 0
          do m=ista_nvale,iend_nvale
             sr = dble(smat(m,n))
             si = dimag(smat(m,n))
             mm = neordr_val(m,ik1)
             if(modnrm /= EXECUT) then
             !   do i=1,kgbpp(n1,n2,n3,id)
                do i=ista_kgbpp(n1,n2,n3,id),iend_kgbpp(n1,n2,n3,id)
                   ii=indgbp(i,2,n1,n2,n3,id)
                   jj=indgbp(i,1,n1,n2,n3,id)
                   fmr = wfv(jj,mm,ik1,1)
                   fmi = wfv(jj,mm,ik1,2) * ipm1
                   grad(ii,1) = grad(ii,1) - sr*fmi - si*fmr
                   grad(ii,2) = grad(ii,2) + sr*fmr - si*fmi
                end do
             else
                ibcount = ibcount+1
                do i=ista_kgbpp(n1,n2,n3,id),iend_kgbpp(n1,n2,n3,id)
                   ii=indgbp(i,2,n1,n2,n3,id)
                   jj=indgbp(i,1,n1,n2,n3,id)
                   fmr = wfv(jj,mm,ik1,1) + gsmatb(ii,1,ibcount)
                   fmi = wfv(jj,mm,ik1,2) * ipm1 + gsmatb(ii,2,ibcount)
                   grad(ii,1) = grad(ii,1) - sr*fmi - si*fmr
                   grad(ii,2) = grad(ii,2) + sr*fmr - si*fmi
                end do
             end if
          end do
          call mpi_allreduce(mpi_in_place,grad,kg1*kimg,mpi_double_precision,mpi_sum &
                &           ,mpi_ke_world,ierr)
          call mpi_allreduce(mpi_in_place,grad,kg1*kimg,mpi_double_precision,mpi_sum &
                &           ,mpi_kg_world,ierr)
          do i=1,iba(ik0)
             grad_fef(i,n,ik0,1) = grad_fef(i,n,ik0,1) + fac * grad(i,1)
             grad_fef(i,n,ik0,2) = grad_fef(i,n,ik0,2) + fac * grad(i,2)
          end do
       end do
       if(modnrm == EXECUT)then
         deallocate(snlmat)
         deallocate(qfmat)
         deallocate(cmat)
       endif
    end if

    if(iopr1==1) then
       if(modnrm == EXECUT)then
          allocate(snlmat(np_g1k(ik1),nlmta));snlmat=dcmplx(0.d0,0.d0)
          allocate(qfmat(nlmta,np_nvale));qfmat=dcmplx(0.d0,0.d0)
          allocate(cmat(np_g1k(ik1),np_nvale));cmat=dcmplx(0.d0,0.d0)
          icount = 0
          do ia=1,natm
             it=ityp(ia)
             mdvdb = m_PP_include_vanderbilt_pot(it)
             if(mdvdb == SKIP) cycle
             do i=ista_g1k(ik1),iend_g1k(ik1)
                nb = nbase(i,ik1)
                ph = -PAI2*dot_product(ngabc(nb,1:3),pos(ia,1:3))
                expgr(i-ista_g1k(ik1)+1) = dcmplx(cos(ph),sin(ph))
             end do
             do lmt1=1,ilmt(it)
                icount = icount+1
                do i=1,np_g1k(ik1)
                   snlmat(i,icount) = expgr(i)*snlt(i,lmt1,(ik1-1)/nspin+1)
                end do
             end do
          enddo
          ibcount = 0
          do m=ista_nvale,iend_nvale
             mm = neordr_val(m,ik0)
             ibcount = ibcount+1
             icount = 0
             do ia=1,natm
                it=ityp(ia)
                mdvdb = m_PP_include_vanderbilt_pot(it)
                if(mdvdb == SKIP) cycle
                ph1 = -PAI2*dot_product(gs0,pos(ia,1:3))
                if(iopr0>1.and.iopr0<=nopr) then
                   vk1 = vkxyz(ik0,1:3,CARTS)
                   dk = vk1 - matmul(op_trs(:,:,iopr0),vk1)
                   ph1 = ph1 + dot_product(dk,cps(ia,1:3))
                else if(iopr0>nopr+1) then
                   vk1 = vkxyz(ik0,1:3,CARTS)
                   dk = -(vk1 + matmul(op_trs(:,:,iopr0),vk1))
                   ph1 = ph1 + dot_product(dk,cps(ia,1:3))
                end if
                exp1 = exp(dcmplx(0.d0,ph1))
                do lmt1=1,ilmt(it)
                   l1 = ltp(lmt1,it)

                   qf = dcmplx(0.d0,0.d0)
                   do lmt2=1,ilmt(it)
                      v = lmta(lmt2,ia)
                      qf = qf + dconjg(ftqe(lmt1,lmt2,it,id))*iwei(ia) &
                          & * dcmplx(fsr_fef(mm,v,ik0),ipm0*fsi_fef(mm,v,ik0))*exp1
                   end do
                   qf = qf * zi**(-l1)
                   icount = icount+1
                   qfmat(icount,ibcount) = qf
                end do
             end do
          enddo
          call zgemm('N','N',np_g1k(ik1),ibcount,icount,cone,snlmat,np_g1k(ik1),qfmat,nlmta,czero,cmat,np_g1k(ik1))
          !do i=1,np_g1k(ik1)
          !   do m=1,ibcount
          !      do lmt1=1,icount
          !         cmat(i,m) = cmat(i,m) + snlmat(i,lmt1)*qfmat(lmt1,m)
          !      enddo
          !   enddo
          !enddo
       endif
       if(modnrm == EXECUT)then
          ibcount = 0
          gsmatb = 0.d0
          do m=ista_nvale,iend_nvale
             ibcount = ibcount+1
             gsmat = 0.d0
             do i=ista_g1k(ik1),iend_g1k(ik1)
                ctmp = cmat(i-ista_g1k(ik1)+1,ibcount)
                gsmat(i,1) = gsmat(i,1) + dble(ctmp)
                gsmat(i,2) = gsmat(i,2) + dimag(ctmp)
             end do
             call mpi_allreduce(mpi_in_place,gsmat,kg1*kimg,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
             gsmatb(:,1:2,ibcount) = gsmat(:,1:2)
          enddo
       endif
       do n=1,ns
          nn = neordr_val(n,ik1)
          grad = 0.d0
!!          do m=1,ns
          ibcount = 0
          do m=ista_nvale,iend_nvale
             sr =  dble(smat(n,m))
             si = -dimag(smat(n,m))
             mm = neordr_val(m,ik0)
             if(modnrm /= EXECUT) then
                !do i=1,kgbpm(m1,m2,m3,id)
                do i=ista_kgbpm(m1,m2,m3,id),iend_kgbpm(m1,m2,m3,id)
                   ii=indgbm(i,2,m1,m2,m3,id)
                   jj=indgbm(i,1,m1,m2,m3,id)
                   fmr = wfv(jj,mm,ik0,1)
                   fmi = wfv(jj,mm,ik0,2) * ipm0
                   grad(ii,1) = grad(ii,1) + sr*fmi + si*fmr
                   grad(ii,2) = grad(ii,2) - sr*fmr + si*fmi
                end do
             else
                ibcount = ibcount+1
                do i=ista_kgbpm(m1,m2,m3,id),iend_kgbpm(m1,m2,m3,id)
                   ii=indgbm(i,2,m1,m2,m3,id)
                   jj=indgbm(i,1,m1,m2,m3,id)
                   fmr = wfv(jj,mm,ik0,1) + gsmatb(ii,1,ibcount)
                   fmi = wfv(jj,mm,ik0,2) * ipm0 + gsmatb(ii,2,ibcount)
                   grad(ii,1) = grad(ii,1) + sr*fmi + si*fmr
                   grad(ii,2) = grad(ii,2) - sr*fmr + si*fmi
                end do
             end if
          end do
          call mpi_allreduce(mpi_in_place,grad,kg1*kimg,mpi_double_precision,mpi_sum &
                &           ,mpi_ke_world,ierr)
          call mpi_allreduce(mpi_in_place,grad,kg1*kimg,mpi_double_precision,mpi_sum &
                &           ,mpi_kg_world,ierr)
          do i=1,iba(ik1)
             grad_fef(i,n,ik1,1) = grad_fef(i,n,ik1,1) + fac * grad(i,1)
             grad_fef(i,n,ik1,2) = grad_fef(i,n,ik1,2) + fac * grad(i,2)
          end do
       end do
       if(modnrm == EXECUT)then
         deallocate(snlmat)
         deallocate(qfmat)
         deallocate(cmat)
       endif
    end if
    if(modnrm == EXECUT)then
      deallocate(gsmat)
      deallocate(gsmatb)
      deallocate(expgr)
    endif
    deallocate(grad)
    call tstatc0_end(id_sname)
  end subroutine  calc_grad
#else
  subroutine  calc_grad(n1,n2,n3,ik0,iopr0,gs0,m1,m2,m3,ik1,iopr1,gs1,ig,id)
    integer, intent(in) :: n1,n2,n3,ik0,iopr0,gs0(3)
    integer, intent(in) :: m1,m2,m3,ik1,iopr1,gs1(3)
    integer, intent(in) :: ig,id

    integer :: m, n, i, ii, jj, mm, nn
    integer :: ipm0,ipm1
    integer :: ia,it,nb,mdvdb
    integer :: lmt1,lmt2,v,l1
    real(kind=DP) :: fac, sr, si, fmr, fmi
    real(kind=DP) :: grad(kg1,kimg),gsmat(kg1,kimg)
    real(kind=DP) :: ph,ph0,ph1
    real(kind=DP) :: dk(3),vk0(3),vk1(3)
    complex(kind=DP) :: exp0,exp1,qf,ctmp,expgr(kg1)
    complex(kind=DP), parameter :: zi=(0.d0,1.d0)
    integer :: ierr
    integer :: id_sname = -1
    call tstatc0_begin('calc_grad ',id_sname,1)

    fac = efa(ig)/PAI4*mp_index(ig)

    if(iopr0>nopr) then
       ipm0=-1
    else
       ipm0=1
    end if
    if(iopr1>nopr) then
       ipm1=-1
    else
       ipm1=1
    end if

    if(iopr0==1) then
       do n=1,ns
          nn = neordr_val(n,ik0)
          grad = 0.d0
          !!do m=1,ns
          do m=ista_nvale,iend_nvale
             sr = dble(smat(m,n))
             si = dimag(smat(m,n))
             mm = neordr_val(m,ik1)
             if(modnrm /= EXECUT) then
                do i=ista_kgbpp(n1,n2,n3,id),iend_kgbpp(n1,n2,n3,id)
                   ii=indgbp(i,2,n1,n2,n3,id)
                   jj=indgbp(i,1,n1,n2,n3,id)
                   fmr = wfv(jj,mm,ik1,1)
                   fmi = wfv(jj,mm,ik1,2) * ipm1
                   grad(ii,1) = grad(ii,1) - sr*fmi - si*fmr
                   grad(ii,2) = grad(ii,2) + sr*fmr - si*fmi
                end do
             else
                gsmat = dcmplx(0.d0,0.d0)
                do ia=1,natm
                   it=ityp(ia)
                   mdvdb = m_PP_include_vanderbilt_pot(it)
                   if(mdvdb == SKIP) cycle
                   ph1 = -PAI2*dot_product(gs1,pos(ia,1:3))
                   if(iopr1>1.and.iopr1<=nopr) then
                      vk1 = vkxyz(ik1,1:3,CARTS)
                      dk = vk1 - matmul(op_trs(:,:,iopr1),vk1)
                      ph1 = ph1 + dot_product(dk,cps(ia,1:3))
                   else if(iopr1>nopr+1) then
                      vk1 = vkxyz(ik1,1:3,CARTS)
                      dk = -(vk1 + matmul(op_trs(:,:,iopr1),vk1))
                      ph1 = ph1 + dot_product(dk,cps(ia,1:3))
                   end if
                   exp1 = exp(dcmplx(0.d0,ph1))
                   do i=1,iba(ik0)
                      nb = nbase(i,ik0)
                      ph = -PAI2*dot_product(ngabc(nb,1:3),pos(ia,1:3))
                      expgr(i) = dcmplx(cos(ph),sin(ph))
                   end do
                   do lmt1=1,ilmt(it)
                      l1 = ltp(lmt1,it)
                      qf = dcmplx(0.d0,0.d0)
                      do lmt2=1,ilmt(it)
                         v = lmta(lmt2,ia)
                         qf = qf + ftqe(lmt1,lmt2,it,id)*iwei(ia) &
                             & * dcmplx(fsr_fef(mm,v,ik1),ipm1*fsi_fef(mm,v,ik1))*exp1
                      end do
                      qf = qf * zi**(-l1)
                      do i=ista_g1k(ik0),iend_g1k(ik0)
                         ctmp = expgr(i)*snlt(i-ista_g1k(ik0)+1,lmt1,(ik0-1)/nspin+1)*qf
                         gsmat(i,1) = gsmat(i,1) + dble(ctmp)
                         gsmat(i,2) = gsmat(i,2) + dimag(ctmp)
                      end do
                   end do
                end do
                call mpi_allreduce(mpi_in_place,gsmat,kg1*kimg,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
                do i=ista_kgbpp(n1,n2,n3,id),iend_kgbpp(n1,n2,n3,id)
                   ii=indgbp(i,2,n1,n2,n3,id)
                   jj=indgbp(i,1,n1,n2,n3,id)
                   fmr = wfv(jj,mm,ik1,1) + gsmat(ii,1)
                   fmi = wfv(jj,mm,ik1,2) * ipm1 + gsmat(ii,2)
                   grad(ii,1) = grad(ii,1) - sr*fmi - si*fmr
                   grad(ii,2) = grad(ii,2) + sr*fmr - si*fmi
                end do
             end if
          end do
          call mpi_allreduce(mpi_in_place,grad,kg1*kimg,mpi_double_precision,mpi_sum &
                &           ,mpi_ke_world,ierr)
          call mpi_allreduce(mpi_in_place,grad,kg1*kimg,mpi_double_precision,mpi_sum &
                &           ,mpi_kg_world,ierr)
          do i=1,iba(ik0)
             grad_fef(i,n,ik0,1) = grad_fef(i,n,ik0,1) + fac * grad(i,1)
             grad_fef(i,n,ik0,2) = grad_fef(i,n,ik0,2) + fac * grad(i,2)
          end do
       end do
    end if

    if(iopr1==1) then
       do n=1,ns
          nn = neordr_val(n,ik1)
          grad = 0.d0
!!          do m=1,ns
          do m=ista_nvale,iend_nvale
             sr =  dble(smat(n,m))
             si = -dimag(smat(n,m))
             mm = neordr_val(m,ik0)
             if(modnrm /= EXECUT) then
                do i=ista_kgbpm(m1,m2,m3,id),iend_kgbpm(m1,m2,m3,id)
                   ii=indgbm(i,2,m1,m2,m3,id)
                   jj=indgbm(i,1,m1,m2,m3,id)
                   fmr = wfv(jj,mm,ik0,1)
                   fmi = wfv(jj,mm,ik0,2) * ipm0
                   grad(ii,1) = grad(ii,1) + sr*fmi + si*fmr
                   grad(ii,2) = grad(ii,2) - sr*fmr + si*fmi
                end do
             else
                gsmat = dcmplx(0.d0,0.d0)
                do ia=1,natm
                   it=ityp(ia)
                   mdvdb = m_PP_include_vanderbilt_pot(it)
                   if(mdvdb == SKIP) cycle
                   ph0 = -PAI2*dot_product(gs0,pos(ia,1:3))
                   if(iopr0>1.and.iopr0<=nopr) then
                      vk0 = vkxyz(ik0,1:3,CARTS)
                      dk = vk0 - matmul(op_trs(:,:,iopr0),vk0)
                      ph0 = ph0 + dot_product(dk,cps(ia,1:3))
                   else if(iopr0>nopr+1) then
                      vk0 = vkxyz(ik0,1:3,CARTS)
                      dk = -(vk0 + matmul(op_trs(:,:,iopr0),vk0))
                      ph0 = ph0 + dot_product(dk,cps(ia,1:3))
                   end if
                   exp0 = exp(dcmplx(0.d0,ph0))
                   do i=1,iba(ik1)
                      nb = nbase(i,ik1)
                      ph = -PAI2*dot_product(ngabc(nb,1:3),pos(ia,1:3))
                      expgr(i) = dcmplx(cos(ph),sin(ph))
                   end do
                   do lmt1=1,ilmt(it)
                      l1 = ltp(lmt1,it)
                      qf = dcmplx(0.d0,0.d0)
                      do lmt2=1,ilmt(it)
                         v = lmta(lmt2,ia)
                         qf = qf + dconjg(ftqe(lmt1,lmt2,it,id))*iwei(ia) &
                             & * dcmplx(fsr_fef(mm,v,ik0),ipm0*fsi_fef(mm,v,ik0))*exp0
                      end do
                      qf = qf * zi**(-l1)
                      do i=ista_g1k(ik1),iend_g1k(ik1)
                         ctmp = expgr(i)*snlt(i-ista_g1k(ik1)+1,lmt1,(ik1-1)/nspin+1)*qf
                         gsmat(i,1) = gsmat(i,1) + dble(ctmp)
                         gsmat(i,2) = gsmat(i,2) + dimag(ctmp)
                      end do
                   end do
                end do
                call mpi_allreduce(mpi_in_place,gsmat,kg1*kimg,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
                do i=ista_kgbpm(m1,m2,m3,id),iend_kgbpm(m1,m2,m3,id)
                   ii=indgbm(i,2,m1,m2,m3,id)
                   jj=indgbm(i,1,m1,m2,m3,id)
                   fmr = wfv(jj,mm,ik0,1) + gsmat(ii,1)
                   fmi = wfv(jj,mm,ik0,2) * ipm0 + gsmat(ii,2)
                   grad(ii,1) = grad(ii,1) + sr*fmi + si*fmr
                   grad(ii,2) = grad(ii,2) - sr*fmr + si*fmi
                end do
             end if
          end do
          call mpi_allreduce(mpi_in_place,grad,kg1*kimg,mpi_double_precision,mpi_sum &
                &           ,mpi_ke_world,ierr)
          call mpi_allreduce(mpi_in_place,grad,kg1*kimg,mpi_double_precision,mpi_sum &
                &           ,mpi_kg_world,ierr)
          do i=1,iba(ik1)
             grad_fef(i,n,ik1,1) = grad_fef(i,n,ik1,1) + fac * grad(i,1)
             grad_fef(i,n,ik1,2) = grad_fef(i,n,ik1,2) + fac * grad(i,2)
          end do
       end do
    end if
    call tstatc0_end(id_sname)
  end subroutine  calc_grad
#endif

  subroutine m_FEF_polarization(nfout,eplr)
    integer, intent(in) :: nfout
    real(kind=DP), intent(out) :: eplr

    integer :: ig,i,ia,id,ii(3)
    integer :: jj1,jj2,m1,m2,i1,i2,j
    real(kind=DP) :: bp(3),pel(3),pion(3)
    real(kind=DP) :: epel,epion
    real(kind=DP) :: fpi,ph
    complex(kind=DP) :: p
    real(kind=8), parameter :: paid2 = PAI/2


    bp(1:3) = 0.d0
    epel = 0.d0
    do ig=1,3
       id = elec_id(ig)
       if(id==0) cycle
       jj1 = mod(ig+1-1,3)+1
       jj2 = mod(ig+2-1,3)+1
       m1 = mp_index(jj1)
       m2 = mp_index(jj2)
       do i1=0,m1-1
          do i2=0,m2-1
             p = dcmplx(1.d0,0.d0)
             do j=0,mp_index(ig)-1
                ii(jj1) = i1
                ii(jj2) = i2
                ii(ig) = j
                p = p * berry_det(ii(1),ii(2),ii(3),id)
!debug
!!    write(nfout,*) 'n1,n2,n3,ig,berry_det=',ii(1),ii(2),ii(3),ig,berry_det(ii(1),ii(2),ii(3),id)
!end debug
             end do
!debug
!!    write(nfout,*) 'p=',p
!end debug
             ph = dimag(log(p))
             if(ph>paid2) then
                ph = ph - PAI
             else if(ph<-paid2) then
                ph = ph + PAI
             end if
             berry(i1,i2,id) = ph
! debug
!!             write(nfout,*) 'i1,i2,ig,ph=',i1,i2,ig,ph
! end debug
          end do
       end do
       bp(ig) = sum(berry(0:m1-1,0:m2-1,id))/(m1*m2)
       epel = epel + efa(ig)*bp(ig)
    end do
    if(nspin==1) then
       fpi = PAI
    else
       fpi = PAI2
    end if
    epel = -epel/fpi
    pel = matmul(altv,bp)/(univol*fpi)

    do i=1,3
       pion(i) = 0.d0
       do ia=1,natm
          pion(i) = pion(i) + ival(ityp(ia)) * cps(ia,i)
       end do
       epion = -dot_product(elec_field,pion)
       pion(i) = pion(i)/univol
    end do

    if(sw_check_polar==ON) then
       eplr = 0.d0
    else
       eplr = epel + epion
    end if

    write(nfout,'("BP   =",3(1x,f20.9))') bp(1:3)
    write(nfout,'("Pel  =",3(1x,f20.9))') pel(1:3)
    write(nfout,'("Pion =",3(1x,f20.9))') pion(1:3)
    write(nfout,'("Pmac =",3(1x,f20.9))') pel(1:3)+pion(1:3)
    pmac_old = pmac;pmac = pel+pion
  end subroutine m_FEF_polarization

  function get_delta_pmac() result(res)
    real(kind=DP) :: res
    res = sum(pmac - pmac_old)/3.d0
  end function get_delta_pmac

  subroutine m_FEF_build_grad
    integer :: ig,id,ik,iopr
    integer :: n1,n2,n3,ip(3)
    integer :: ik0,ik1,iopr0,iopr1
    integer :: gs0(3),gs1(3)
    integer :: ierr
    integer :: m1,m2,m3
    m1 = mp_index(1)
    m2 = mp_index(2)
    m3 = mp_index(3)
    allocate(berry_dett(0:m1-1,0:m2-1,0:m3-1))
    call gather_valence_states(nfout)

    call alloc_smat
    do ig=1,3
       id = elec_id(ig)
       if(id==0) cycle
       berry_dett = (0.d0,0.d0)
       do ik=1,kv3
          if(map_k(ik) /= myrank_k) cycle
          do iopr=1,nopr_trs
             if(.not.flag_op(iopr,ik)) cycle
             n1=map_ik2nnn(1,iopr,ik)
             n2=map_ik2nnn(2,iopr,ik)
             n3=map_ik2nnn(3,iopr,ik)
             ip(1) = n1
             ip(2) = n2
             ip(3) = n3
             ip(ig) = ip(ig) + 1
             gs1 = 0
             if(ip(ig)>=mp_index(ig)) then
                ip(ig)=0
                gs1(ig) = 1
             end if
             ik0 = map_nnn2ik(n1,n2,n3)
             !!ik0 = ik
             ik1 = map_nnn2ik(ip(1),ip(2),ip(3))
             iopr0 = map_nnn2op(n1,n2,n3)
             iopr1 = map_nnn2op(ip(1),ip(2),ip(3))
             gs0 = map_nnn2gs(1:3,n1,n2,n3)
             gs1 = gs1 + map_nnn2gs(1:3,ip(1),ip(2),ip(3))
             call calc_overlap(n1,n2,n3,ik0,iopr0,gs0,ik1,iopr1,gs1,ig,id) !=> smat
             call calc_berry_phase(n1,n2,n3,id) !=> smat,berry_det
             if(sw_check_polar==OFF) &
              & call calc_grad(n1,n2,n3,ik0,iopr0,gs0,ip(1),ip(2),ip(3),ik1,iopr1,gs1,ig,id) !=> grad_fef
          end do
       end do
       if(nrank_k>1) &
       & call mpi_allreduce(MPI_IN_PLACE,berry_dett,m1*m2*m3,mpi_double_complex,mpi_sum,mpi_ge_world,ierr)
       berry_det(:,:,:,id) = berry_dett(:,:,:)
    end do
    if(nrank_k>1) &
    & call mpi_allreduce(MPI_IN_PLACE,grad_fef,kg1*ns*kv3*kimg,mpi_double_precision,mpi_sum,mpi_ge_world,ierr)

    call dealloc_smat
    deallocate(berry_dett)
  end subroutine m_FEF_build_grad

  subroutine m_FEF_add_grad_to_vnlph(ik)
    integer, intent(in) :: ik

    integer :: ib,jb,ig
    real(kind=DP),allocatable,dimension(:,:) :: gradt

    if(sw_check_polar==ON) return

    if(map_k(ik) /= myrank_k) return

    allocate(gradt(kg1,kimg));gradt=0.d0
    do ib=1,ns
       jb = neordr(ib,ik)
       if(map_e(jb) == myrank_e) then
         do ig=1,np_g1k(ik)
            gradt(ig,1:kimg) = grad_fef(ista_g1k(ik)-1+ig,ib,ik,1:kimg)
         enddo
         call m_ES_add_it_to_vnlph(ik,map_z(jb),gradt)
       endif
    end do
    deallocate(gradt)
  end subroutine m_FEF_add_grad_to_vnlph

  subroutine m_FEF_add_grad_to_vnlph_RMM(ik,jb)
    integer, intent(in) :: ik,jb

    integer :: ib

    if(sw_check_polar==ON) return

    if(map_ek(jb,ik) == mype) then
       ib = nrvf_ordr(jb,ik)
       if(ib>ns) return
       call m_ES_add_it_to_vnlph(ik,jb,grad_fef(:,ib,ik,:))
    end if
  end subroutine m_FEF_add_grad_to_vnlph_RMM

  subroutine alloc_grad
    allocate(grad_fef(kg1,ns,kv3,kimg))
  end subroutine alloc_grad

  subroutine dealloc_grad
    deallocate(grad_fef)
  end subroutine dealloc_grad

  subroutine m_FEF_Constract_of_ftq()
    ! local varialbes
    integer :: it,ik,ip,lmt1,lmt2,u,v
    integer :: il1,il2,tau1,tau2,l3
    integer :: n,ilm3,iiqitg,ig,id
    integer :: mdvdb
    integer, allocatable :: il3(:)
    real(kind=DP) :: fac
    real(kind=DP) :: dk(3),ylm,dga
    real(kind=DP) :: ftqr(numef),ftqi(numef),ftqb
    complex(kind=DP) :: zi = (0.d0,1.d0)

    if(modnrm /= EXECUT) return

    call m_PP_find_maximum_l(n)   !  n-1: maximum l
    n = (n-1) + (n-1) + 1
    allocate(il3(n**2)); call substitute_il3(n**2,il3) ! -(b_Elec..)
    allocate(ftqe(nlmt,nlmt,ntyp,numef)); ftqe=(0.d0,0.d0)
    TYPE: do it = 1, ntyp
       mdvdb = m_PP_include_vanderbilt_pot(it)
       if(mdvdb == SKIP) cycle TYPE
       LMT_1: do lmt1=1,ilmt(it)
          il1 = ltp(lmt1,it)
          tau1 = taup(lmt1,it)
          LMT_2: do lmt2=lmt1,ilmt(it)
             il2 = ltp(lmt2,it)
             tau2 = taup(lmt2,it)
             ftqr(1:numef) = 0.d0
             ftqi(1:numef) = 0.d0
             LM3: do n=1,il2p(lmt1,lmt2,it)
                ilm3 = isph(lmt1,lmt2,n,it); l3=il3(ilm3)
                iiqitg = iqitg(il1,tau1,il2,tau2,l3+1,it)
                if(iiqitg == 0) cycle LM3
                do ig=1,3
                   id = elec_id(ig)
                   if(id==0) cycle
                   dk(1:3)=-rltv(1:3,ig)/mp_index(ig)
                   call sphrp2_for_Berry(ilm3,dk,ylm)
                   ftqb = qitg_fef(iiqitg,id)*dl2p(lmt1,lmt2,n,it)*ylm
                   !!$ write(6,*) 'it,lmt1,lmt2,n,id,ftqb=',it,lmt1,lmt2,n,id,ftqb
                   if(mod(l3,2)==0) then
                      ftqr(id)=ftqr(id)+real(zi**(-l3))*ftqb
                   else
                      ftqi(id)=ftqi(id)+dimag(zi**(-l3))*ftqb
                   end if
                end do
             end do LM3
             do id=1,numef
                ftqe(lmt1,lmt2,it,id) = dcmplx(ftqr(id),ftqi(id))
                ftqe(lmt2,lmt1,it,id) = ftqe(lmt1,lmt2,it,id)
             end do
          end do LMT_2
       end do LMT_1
    end do TYPE

    deallocate(il3)
  end subroutine m_FEF_Constract_of_ftq

  subroutine gather_valence_states(nfout)
    implicit none
    integer, intent(in) :: nfout

    integer :: ik,ib,ib1,ibm,irev
    real(kind=DP), allocatable :: wfv_mpi(:,:,:,:) ! d(kg1,nval,kv3,kimg)
    real(kind=DP), allocatable :: occup_val_mpi(:,:) ! d(nval,kv3)
    real(kind=DP), allocatable :: neordr_mpi(:,:) ! d(nval,kv3)
    real(kind=DP), allocatable :: fsr_mpi(:,:,:) ! d(nval,nlmta,kv3)
    real(kind=DP), allocatable :: fsi_mpi(:,:,:) ! d(nval,nlmta,kv3)
    real(kind=DP), allocatable :: wfvt(:,:)

    real(kind=DP), allocatable :: zaj_buf(:,:),zaj_buf2(:,:)
    real(kind=DP), allocatable :: fsrt(:),fsit(:)
    integer :: ig,iadd
    integer :: ip0,ip1
    integer :: ilmta,lmt,ia,it

    integer, allocatable, dimension(:) :: ista
    integer :: ierr
    integer :: iksnl
    integer,save  :: id_sname = -1
    logical, save :: firstcall = .true.

    if(firstcall)then
      allocate(snlt(maxval(np_g1k),nlmtt,kv3));snlt=0.d0
      do iksnl=ista_snl,iend_snl
        do ig=1,np_g1k(ik)
          snlt(ig,1:nlmtt,ik) = snl(ig,1:nlmtt,ik)
        enddo
      enddo
      call mpi_allreduce(mpi_in_place,snlt,maxval(np_g1k)*nlmtt*kv3/nspin &
         & ,mpi_double_precision,mpi_sum,mpi_ge_world,ierr)
      firstcall = .false.
    endif

    allocate(ista(MPI_STATUS_SIZE))

    call tstatc0_begin('gather_valence_states ',id_sname,level=1)

    ns = int(totch + 1.d-13)/2

    call m_Parallel_init_mpi_nvale(nfout,ipriparallel,printable,ns)

    if(allocated(wfv)) deallocate(wfv)
    if(allocated(neordr_val)) deallocate(neordr_val)
    !allocate(wfv(maxval(np_g1k),ns,kv3,kimg))
    allocate(wfv(kg1,ns,kv3,kimg))
    allocate(wfvt(kg1,kimg))
    if(modnrm==EXECUT)then
      allocate(fsrt(nlmta));allocate(fsit(nlmta))
      fsrt=0.d0;fsit=0.d0
    endif
    allocate(neordr_val(ns,kv3))
    wfv = 0.d0
    neordr_val = 0.d0
    if(modnrm == EXECUT) then
       if(allocated(fsr_fef)) deallocate(fsr_fef)
       if(allocated(fsi_fef)) deallocate(fsi_fef)
       allocate(fsr_fef(ns,nlmta,kv3))
       allocate(fsi_fef(ns,nlmta,kv3))
       fsr_fef = 0.d0
       fsi_fef = 0.d0
    end if
    do ik=1,kv3,af+1
       if(map_k(ik) /= myrank_k) cycle
       do ib=1,ns
          ib1 = neordr(ib,ik)
          if(map_e(ib1) == myrank_e) then
             wfvt=0.d0
             do ig=ista_g1k(ik),iend_g1k(ik)
                wfv(ig,ib1,ik,1:kimg) = zaj_l(ig-ista_g1k(ik)+1,map_z(ib1),ik,1:kimg)
!                wfvt(ig,1:kimg) = zaj_l(ig-ista_g1k(ik)+1,map_z(ib1),ik,1:kimg)
             enddo
             !call mpi_allreduce(mpi_in_place,wfvt,kg1*kimg,mpi_double_precision, mpi_sum,mpi_ke_world,ierr)
             !wfv(:,ib1,ik,:) = wfvt(:,:)
             neordr_val(ib,ik) = neordr(ib,ik)
             if(modnrm == EXECUT) then
                do ilmta=ista_fs,iend_fs
                   iadd = ilmta - ista_fs + 1
                   fsr_fef(ib1,ilmta,ik) = fsr_l(map_z(ib1),iadd,ik)
                   fsi_fef(ib1,ilmta,ik) = fsi_l(map_z(ib1),iadd,ik)
                enddo
             end if
          end if
       end do
    end do
!    call mpi_allreduce(mpi_in_place,wfv,kg1*ns*kv3*kimg,mpi_double_precision, mpi_sum,mpi_kg_world,ierr)
!    call mpi_allreduce(mpi_in_place,wfv,kg1*ns*kv3*kimg,mpi_double_precision, mpi_sum,mpi_ge_world,ierr)
    call mpi_allreduce(mpi_in_place,wfv,kg1*ns*kv3*kimg,mpi_double_precision, mpi_sum,MPI_CommGroup,ierr)

    call mpi_allreduce(mpi_in_place,neordr_val,ns*kv3,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_kg_world,ierr)
    call mpi_allreduce(mpi_in_place,neordr_val,ns*kv3,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ge_world,ierr)
    deallocate(wfvt)
    if(modnrm == EXECUT) then
       allocate(fsr_mpi(ns,nlmta,kv3))
       allocate(fsi_mpi(ns,nlmta,kv3))
       call mpi_allreduce(fsr_fef,fsr_mpi,ns*nlmta*kv3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_CommGroup,ierr)
       call mpi_allreduce(fsi_fef,fsi_mpi,ns*nlmta*kv3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_CommGroup,ierr)
       fsr_fef = fsr_mpi
       fsi_fef = fsi_mpi
       deallocate(fsr_mpi)
       deallocate(fsi_mpi)
    end if

    if(modnrm==EXECUT)then
      deallocate(fsrt)
      deallocate(fsit)
    endif
    deallocate(ista)
    call tstatc0_end(id_sname)
  end subroutine gather_valence_states

end module m_FiniteElectricField
