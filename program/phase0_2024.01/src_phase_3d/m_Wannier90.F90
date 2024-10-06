#define WAN90_SAVE_MEMORY
#define WAN90_SKIP_FFT
#define WAN90_SPN_FORMATTED
!=======================================================================
!
!  PROGRAM  PHASE/0 2014.01 ($Rev: 110 $)
!
!  MODULE: m_Wannier90
!
!  AUTHOR(S): T. Yamamoto   May/05/2008
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
module m_Wannier90
  use m_Const_Parameters,   only : DP,PAI2,PAI4,SKIP,EXECUT,GAMMA,BUCS,CRDTYP, CMPLDP, Hartree
  use m_Control_Parameters, only : printable, kimg, nspin, neg, ipri &
    &                            , wan90_seedname, nb_wan90, LEN_TITLE
  use m_Files,              only : nfwannier
  use m_Timing,             only : tstatc0_begin, tstatc0_end
  use m_Kpoints,            only : kv3,vkxyz,k_symmetry
  use m_PlaneWaveBasisSet,  only : ngabc,kg1,kgp,nbase,iba,nbase_gamma
  use m_Electronic_Structure,only: neordr,zaj_l,eko_l,m_ES_WF_in_Rspace_3D
  use m_FFT,                only : fft_box_size_WF, nfft, m_FFT_alloc_WF_work &
       &                         , m_FFT_dealloc_WF_work
  use m_Crystal_Structure,  only : univol,rltv
  use m_Parallelization,    only : MPI_CommGroup,npes,mype &
       &                         , nrank_e,myrank_e,map_e,ista_e,iend_e,istep_e,idisp_e &
       &                         , map_z,np_e,mpi_k_world,myrank_k,map_k,ista_k,iend_k &
       &                         , ista_snl, iend_snl, ierr, map_ek, nrank_k
! === KT_add === 2015/02/23, 09/02
  use m_Const_Parameters,   only : Delta07, zi, ON, UP, DOWN, BOHR
  use m_Control_Parameters,  only : ndim_spinor, sw_use_hardpart_wan90, noncol, &
       &                            spin_component_wan90
  use m_Ionic_System,  only : natm, ityp, cps, pos
  use m_PseudoPotential,   only : dk_wan, nloc, nlmta, m_PP_include_vanderbilt_pot, &
       &                          ltp, mtp, taup, ilmt, lmta, il2p, isph, iqitg, dl2p, &
       &                          qitg_wan, nqitg, phirpw, psirpw, q
  use m_Electronic_Structure,only: fsr_l, fsi_l
! ============== 2015/02/23, 09/02

! === KT_add === 2015/09/14
  use m_Const_Parameters,  only : ELECTRON, DIRECT, INVERSE
  use m_FFT,  only : m_FFT_WF
  use m_PlaneWaveBasisSet,  only : kg1, kg
! ============== 2015/09/14
  use mpi

  implicit none

  character(len=LEN_TITLE) :: comment
  logical :: calc_only_A
  real(kind=DP) :: real_lattice(3,3), recip_lattice(3,3)
  integer :: num_kp, n_proj, n_exclude, nntot, num_bands
  integer, parameter :: lmax1_pj = 3

  real(kind=DP), allocatable :: kvec(:,:) ! d(3,num_kp)
  real(kind=DP), allocatable :: centre(:,:), zaxis(:,:), xaxis(:,:) ! d(3,n_proj)
  real(kind=DP), allocatable :: zona(:) ! d(n_proj)
  integer, allocatable :: lang(:), mr(:), irf(:) ! d(n_proj)
  integer, allocatable :: exclude_bands(:) ! d(n_exclude)
  integer, allocatable :: nnlist(:,:) ! d(num_kp,nntot)
  integer, allocatable :: nncell(:,:,:) ! d(3,num_kp,nntot)
  integer, allocatable :: ib_inc(:) ! d(num_bands)

  real(kind=DP), allocatable :: projfunc(:,:,:,:) ! d(kg1,n_proj,ista_snl:iend_snl,kimg)

! ==== KT_add === 2015/04/13 & 09/02
  real(kind=DP), allocatable :: dk_unit(:,:)
  logical, allocatable :: centre_on_atom(:)

  integer, allocatable :: spn_index(:)
  real(kind=DP), allocatable :: spn_quant_dir(:,:)
! =============== 2015/04/13 & 09/02

! ==== KT_add === 2015/09/14
  integer, allocatable :: igf_wan90(:,:,:,:)
! =============== 2015/09/14

!  include 'mpif.h'                                      ! MPI
  integer istatus(mpi_status_size)                      ! MPI

contains

  subroutine m_Wan90_init(nfout)
    implicit none
    integer, intent(in) :: nfout

    character(len=15) :: dummy

    integer :: i, j, n
    integer :: kp, nnkp, ng1,ng2,ng3

! == KT add == 2015/02/23
    if ( mype /= 0 ) goto 100
! ============== 2015/02/23

    open(nfwannier,file=trim(wan90_seedname)//".nnkp",status="old",form="formatted")
    read(nfwannier,'(a80)') comment

    read(nfwannier,*)

    read(nfwannier,'(a15,l)') dummy,calc_only_A

    read(nfwannier,*)

    read(nfwannier,*)
    do i=1,3
       read(nfwannier,*) real_lattice(1:3,i)
    end do
    read(nfwannier,*)

    read(nfwannier,*)

    read(nfwannier,*)
    do i=1,3
       read(nfwannier,*) recip_lattice(1:3,i)
    end do
    read(nfwannier,*)

    read(nfwannier,*)

    read(nfwannier,*)
    read(nfwannier,*) num_kp
    if(ipri>=2) write(*,'("num_kp=",i5)') num_kp
    allocate(kvec(3,num_kp))
    do i=1,num_kp
       read(nfwannier,*) kvec(1:3,i)
    end do
    read(nfwannier,*)

    read(nfwannier,*)

    read(nfwannier,*)
    read(nfwannier,*) n_proj
    if(ipri>=2) write(*,'("n_proj=",i5)') n_proj
    allocate(centre(3,n_proj))
    allocate(lang(n_proj))
    allocate(mr(n_proj))
    allocate(irf(n_proj))
    allocate(zaxis(3,n_proj))
    allocate(xaxis(3,n_proj))
    allocate(zona(n_proj))

    if ( noncol ) then
       allocate( spn_index(n_proj) )
       allocate( spn_quant_dir(3,n_proj) )
    endif

    if ( noncol ) then
       do i=1,n_proj
          read(nfwannier,*) centre(1:3,i), lang(i), mr(i), irf(i)
          read(nfwannier,*) zaxis(1:3,i), xaxis(1:3,i), zona(i)
!
          read(nfwannier,*) spn_index(i), spn_quant_dir(1:3,i)
       end do
    else
       do i=1,n_proj
          read(nfwannier,*) centre(1:3,i), lang(i), mr(i), irf(i)
          read(nfwannier,*) zaxis(1:3,i), xaxis(1:3,i), zona(i)
       end do
    endif

    read(nfwannier,*)

    read(nfwannier,*)

    read(nfwannier,*)
    read(nfwannier,*) nntot
    if(ipri>=2) write(*,'("nntot=",i5)') nntot
    allocate(nnlist(num_kp,nntot))
    allocate(nncell(3,num_kp,nntot))
    do j=1,num_kp
       do i=1,nntot
          read(nfwannier,*) kp, nnkp, ng1,ng2,ng3
          nnlist(kp,i) = nnkp
          nncell(1,kp,i) = ng1
          nncell(2,kp,i) = ng2
          nncell(3,kp,i) = ng3
      end do
    end do
    read(nfwannier,*)

    read(nfwannier,*)

    read(nfwannier,*)
    read(nfwannier,*) n_exclude
    if(ipri>=2) write(*,'("n_exclude=",i5)') n_exclude
    allocate(exclude_bands(n_exclude))
    do i=1,n_exclude
       read(nfwannier,*) exclude_bands(i)
    end do
    read(nfwannier,*)

    close(nfwannier)

! == KT add == 2015/02/23 & 09/02
100 continue
    if ( npes > 1 ) then
       call mpi_bcast( calc_only_A, 1, mpi_logical, 0, MPI_CommGroup, ierr )
       call mpi_bcast( real_lattice, 3*3, mpi_double_precision, 0, &
            &          MPI_CommGroup, ierr )
       call mpi_bcast( recip_lattice, 3*3, mpi_double_precision, 0, &
            &          MPI_CommGroup, ierr )

       call mpi_bcast( num_kp, 1, mpi_integer, 0, MPI_CommGroup, ierr )
       if ( mype /= 0 ) then
          allocate( kvec(3,num_kp) );  kvec = 0.0d0
       endif
       call mpi_bcast( kvec, 3*num_kp, mpi_double_precision, 0, &
            &          MPI_CommGroup, ierr )

       call mpi_bcast( n_proj, 1, mpi_integer, 0, MPI_CommGroup, ierr )
       if ( mype /= 0 ) then
          allocate(centre(3,n_proj)); allocate(zona(n_proj))
          allocate(zaxis(3,n_proj));  allocate(xaxis(3,n_proj))
          allocate(lang(n_proj));     allocate(mr(n_proj));    allocate(irf(n_proj))
       endif
       call mpi_bcast( centre, 3*n_proj, mpi_double_precision, 0, &
            &          MPI_CommGroup, ierr )
       call mpi_bcast( zaxis,  3*n_proj, mpi_double_precision, 0, &
            &          MPI_CommGroup, ierr )
       call mpi_bcast( xaxis,  3*n_proj, mpi_double_precision, 0, &
            &          MPI_CommGroup, ierr )
       call mpi_bcast( zona,     n_proj, mpi_double_precision, 0, &
            &          MPI_CommGroup, ierr )

       if ( noncol ) then
          if ( mype /= 0 ) then
             allocate( spn_index(n_proj) )
             allocate( spn_quant_dir(3,n_proj) )
          endif
          call mpi_bcast( spn_index, n_proj, mpi_integer, 0, &
               &          MPI_CommGroup, ierr )
          call mpi_bcast( spn_quant_dir, 3*n_proj, mpi_double_precision, 0, &
               &          MPI_CommGroup, ierr )
       endif

       call mpi_bcast( lang, n_proj, mpi_integer, 0, MPI_CommGroup, ierr )
       call mpi_bcast( mr,   n_proj, mpi_integer, 0, MPI_CommGroup, ierr )
       call mpi_bcast( irf,  n_proj, mpi_integer, 0, MPI_CommGroup, ierr )

       call mpi_bcast( nntot, 1, mpi_integer, 0, MPI_CommGroup, ierr )
       if ( mype /= 0 ) then
          allocate(nnlist(num_kp,nntot));    allocate(nncell(3,num_kp,nntot))
       endif
       call mpi_bcast( nnlist,   num_kp*nntot, mpi_integer, 0, MPI_CommGroup, ierr )
       call mpi_bcast( nncell, 3*num_kp*nntot, mpi_integer, 0, MPI_CommGroup, ierr )

       call mpi_bcast( n_exclude, 1, mpi_integer, 0, MPI_CommGroup, ierr )
       if ( mype /= 0 ) then
          allocate(exclude_bands(n_exclude));  exclude_bands = 0
       endif
       call mpi_bcast( exclude_bands, n_exclude, mpi_integer, 0, MPI_CommGroup, ierr )
    endif
! ============= 2015/02/23 & 09/02

    num_bands = nb_wan90 - n_exclude
    allocate(ib_inc(num_bands))
    n = 0
    do i=1,nb_wan90
       if(excluding_band(i)) cycle
       n = n + 1
       ib_inc(n) = i
    end do

    if ( mype == 0 ) call write_down_input(nfout)

  end subroutine m_Wan90_init

  subroutine write_down_input(nfout)
    implicit none
    integer, intent(in) :: nfout

    integer :: i,j

    write(nfout,'("Wannier90")')
    write(nfout,'("comment: ",a)') trim(comment)
    write(nfout,'("calc_only_A: ",l)') calc_only_A
    write(nfout,'("real_lattice")')
    do i=1,3
       write(nfout,'(3(1x,f20.8))') real_lattice(1:3,i)
    end do
    write(nfout,'("recip_lattice")')
    do i=1,3
       write(nfout,'(3(1x,f20.8))') recip_lattice(1:3,i)
    end do
    write(nfout,'("k vectors:")')
    write(nfout,'("num_kp: ",i5)') num_kp
    do i=1,num_kp
       write(nfout,'("kvec(",i5,") =",3(1x,f20.8))') i,kvec(1:3,i)
    end do
    write(nfout,'("n_proj: ",i5)') n_proj
    do i=1,n_proj
       write(nfout,'(3(1x,f20.8),3(1x,i5))') centre(1:3,i), lang(i), mr(i), irf(i)
       write(nfout,'(7(1x,f20.8))') zaxis(1:3,i), xaxis(1:3,i), zona(i)
    end do
    write(nfout,'("kp nnlist nncell(1:3)")')
    do i=1,num_kp
       do j=1,nntot
          write(nfout,'(i5,1x,i5,3(1x,i5))') i,nnlist(i,j),nncell(1:3,i,j)
       end do
    end do
    write(nfout,'("excluded bands ")')
    write(nfout,'("n_exclude: ",i5)') n_exclude
    do i=1,n_exclude
       write(nfout,'(i5)') exclude_bands(i)
    end do
    write(nfout,'("included bands ")')
    do i=1,num_bands
       write(nfout,'(i5)') ib_inc(i)
    end do

  end subroutine write_down_input

  subroutine m_Wan90_gen_amat(nfout)
!
!    Amn(k) = <Psi_mk|g_n>
!
    implicit none
    integer, intent(in) :: nfout

    real(kind=DP), allocatable :: a_mat(:,:,:,:) ! d(num_bands,n_proj,kv3,2)

    allocate(a_mat(num_bands,n_proj,kv3,2));  a_mat = 0.d0

    call contrib_softpart
    if ( sw_use_hardpart_wan90 == ON ) call contrib_hardpart

    call gather_matrix
    call print_mat

    deallocate(a_mat)

  contains

    subroutine contrib_softpart
      integer :: ik,ikpj,n,m,ib,i,mi
      integer :: ik_start, ik_skip
      real(kind=DP) :: pr,pi,zr,zi

      allocate(projfunc(kg1,n_proj,ista_snl:iend_snl,kimg))
      call m_Wan90_Projectors(nfout,kv3,vkxyz)

      ik_start = 1;  ik_skip = 1
      if ( .not. noncol ) then
         if ( nspin == 2 .and. spin_component_wan90 == DOWN ) then
            ik_start = 2
         endif
         if ( nspin == 2 ) ik_skip = 2
      endif

      do ik = ik_start, kv3, ik_skip
         if(map_k(ik) /= myrank_k) cycle ! MPI
         ikpj = (ik-1)/nspin + 1

         do n=1,n_proj
            do m=1,num_bands
               mi = ib_inc(m)
               ib = neordr(mi,ik)
               if(map_e(ib) /= myrank_e) cycle ! MPI
               if(ipri>=2) write(*,'("m=",i5," n=",i5," ik=",i5)') m,n,ik
               if(kimg==1) then
                  do i=1,iba(ik)
                     pr = projfunc(i,n,ikpj,1)
                     zr = zaj_l(i,map_z(ib),ik,1)
                     a_mat(m,n,ik,1) = a_mat(m,n,ik,1) + zr*pr
                  end do
               else
                  if(k_symmetry(ik) == GAMMA) then
                     do i=2,iba(ik)
                        pr = projfunc(i,n,ikpj,1)
                        pi = projfunc(i,n,ikpj,2)
                        zr = zaj_l(i,map_z(ib),ik,1)
                        zi = zaj_l(i,map_z(ib),ik,2)
                        a_mat(m,n,ik,1) = a_mat(m,n,ik,1) + zr*pr + zi*pi
                     end do
                     pr = projfunc(1,n,ikpj,1)
                     zr = zaj_l(1,map_z(ib),ik,1)
                     a_mat(m,n,ik,1) = a_mat(m,n,ik,1)*2.d0 + zr*pr
                  else
                     do i=1,iba(ik)
                        pr = projfunc(i,n,ikpj,1)
                        pi = projfunc(i,n,ikpj,2)
                        zr = zaj_l(i,map_z(ib),ik,1)
                        zi = zaj_l(i,map_z(ib),ik,2)
                        a_mat(m,n,ik,1) = a_mat(m,n,ik,1) + zr*pr + zi*pi
                        a_mat(m,n,ik,2) = a_mat(m,n,ik,2) + zr*pi - zi*pr
                     end do
                  end if
               end if
            end do
         end do
      end do
      deallocate(projfunc)

    end subroutine contrib_softpart

! === KT_add === 2015/02/23
    subroutine contrib_hardpart
      integer, parameter :: nmesh = 1501

      integer :: ik ,ikpj, n, m, mi, ib
      integer :: ia, it, lmt1, il1, im1, tau1, ir, nsphere, lmta1
      integer :: ik_start, ik_skip
      real(kind=DP) :: dx, dy, dz, dist
      real(kind=DP) :: coeff_w90(16), csum, ph
      complex(kind=CMPLDP) :: z1

      real(kind=DP), allocatable, dimension(:) :: radr,wos,radfunc
      complex(kind=CMPLDP), allocatable :: zph(:)

      allocate( zph(natm) ); zph = 0.0d0

      allocate(radr(nmesh),wos(nmesh),radfunc(nmesh))
      call radr_and_wos(nmesh,radr,wos) ! --> radr, wos

      ik_start = 1;  ik_skip = 1
      if ( .not. noncol ) then
         if ( nspin == 2 .and. spin_component_wan90 == DOWN ) then
            ik_start = 2
         endif
         if ( nspin == 2 ) ik_skip = 2
      endif

      do ik = ik_start, kv3, ik_skip
         if(map_k(ik) /= myrank_k) cycle ! MPI
         ikpj = (ik-1)/nspin + 1

         Do ia=1, natm
            ph = PAI2 * dot_product( pos(ia,:), vkxyz(ik,:,BUCS) )
            zph(ia) = dcmplx( cos(ph), -sin(ph) )
         End Do

         do n=1,n_proj
            call radial_function( irf(n), zona(n), nmesh, radr, radfunc )
            call deompose_sphr_into_lm( lang(n), mr(n), coeff_w90 )
!
            do m=1,num_bands
               mi = ib_inc(m)
               if ( map_e(mi) /= myrank_e ) cycle

               ib = neordr(mi,ik)

               Do ia=1, natm
                  dx = centre(1,n) -pos(ia,1)
                  dy = centre(2,n) -pos(ia,2)
                  dz = centre(3,n) -pos(ia,3)
                  dist = sqrt( dx**2 +dy**2 +dz**2 )
                  if ( dist > 1.0E-3 ) cycle
!
                  it = ityp(ia)
                  Do lmt1=1, ilmt(it)
                     il1 = ltp(lmt1,it);  im1 = mtp(lmt1,it);  tau1 = taup(lmt1,it)
                     lmta1 = lmta(lmt1,ia)
                     !
                     csum = 0.0d0
                     Do ir=1, nmesh
                        csum = csum + wos(ir) *radr(ir) *radfunc(ir) &
                             &   *( psirpw(ir,il1,tau1,it) -phirpw(ir,il1,tau1,it) )
                     End do

                     nsphere = ( il1 -1 )**2 +im1
                     csum = csum *coeff_w90( nsphere )

                     z1 = csum *zph(ia) &
                          &     *dcmplx( fsr_l(map_z(ib),lmta1,ik), &
                          &             -fsi_l(map_z(ib),lmta1,ik) )

                     a_mat(m,n,ik,1) = a_mat(m,n,ik,1) + real(z1)
                     a_mat(m,n,ik,2) = a_mat(m,n,ik,2) + aimag(z1)
                  End Do
               End Do
            End do
         End do
      End do
      deallocate( radr ); deallocate( wos ); deallocate( radfunc )
    end subroutine contrib_hardpart
! ============== 2015/02/23

    subroutine gather_matrix
      real(kind=DP), allocatable :: a_mat_mpi(:,:,:,:) ! d(num_bands,n_proj,kv3,2)

      if(npes>1) then
         allocate(a_mat_mpi(num_bands,n_proj,kv3,2))
         a_mat_mpi = a_mat
         a_mat = 0.d0
         call mpi_allreduce(a_mat_mpi,a_mat,num_bands*n_proj*kv3*2,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
         deallocate(a_mat_mpi)
      end if
    end subroutine gather_matrix

    subroutine print_mat
      integer :: ik ,m, n, ik_start, iktmp, is
      real(kind=DP) :: theta, phi
      complex(kind=CMPLDP) :: ztmp, z1
      complex(kind=CMPLDP), allocatable :: spn_weight(:,:)

      if ( mype /= 0 ) return

      open(nfwannier,file=trim(wan90_seedname)//".amn",form="formatted")
      if ( noncol ) then
         write(nfwannier,*) 'Generated by PHASE noncol'
      else
         if ( nspin == 2 ) then
            write(nfwannier,*) 'Generated by PHASE (nspin=2), comp: ', &
                 &              spin_component_wan90
         else
            write(nfwannier,*) 'Generated by PHASE'
         endif
      endif
      write(nfwannier,*) num_bands, kv3/nspin, n_proj

      if ( noncol ) then
         allocate( spn_weight(2,n_proj) ); spn_weight = 0.0d0

         Do n=1, n_proj
            theta = acos( spn_quant_dir(3,n) )
            phi = atan2( spn_quant_dir(2,n), spn_quant_dir(1,n) )
!
            if ( spn_index(n) == 1 ) then
               spn_weight(1,n) = exp( zi *phi /2.0d0 ) *cos( theta /2.0d0 )
               spn_weight(2,n) = -exp( zi *phi /2.0d0 ) *sin( theta /2.0d0 )
            else
               spn_weight(1,n) = exp( -zi *phi /2.0d0 ) *sin( theta /2.0d0 )
               spn_weight(2,n) = exp( -zi *phi /2.0d0 ) *cos( theta /2.0d0 )
            endif
         End do

         do ik = 1, kv3, ndim_spinor
            iktmp = ( ik -1 )/ndim_spinor +1
            do n=1,n_proj
               do m=1,num_bands
                  ztmp = 0.0d0
                  Do is=1, ndim_spinor
                     z1 = dcmplx( a_mat(m,n,ik+is-1,1), a_mat(m,n,ik+is-1,2 ) )
                     ztmp = ztmp +spn_weight(is,n) *z1
                  End do

                  write(nfwannier,'(3(1x,i5),2(1x,f18.12))') m, n, iktmp, &
                       &         real(ztmp), aimag(ztmp)
               end do
            end do
         end do
         deallocate( spn_weight )

      else
         ik_start = 1
         if ( nspin == 2 .and. spin_component_wan90 == DOWN ) ik_start = 2

         do ik = ik_start, kv3, nspin
            iktmp = ( ik -1 )/nspin +1
            do n=1,n_proj
               do m=1,num_bands
                  write(nfwannier,'(3(1x,i5),2(1x,f18.12))') m, n, iktmp, &
                       &                     a_mat(m,n,ik,1:2)
               end do
            end do
         end do
      endif

      close(nfwannier)

    end subroutine print_mat

  end subroutine m_Wan90_gen_amat

  subroutine m_Wan90_Projectors(nfout,kv3,vkxyz)
    implicit none
    integer, intent(in) :: nfout,kv3
    real(kind=DP), intent(in) :: vkxyz(kv3,3,CRDTYP)

    integer, parameter :: nmesh = 1501

    integer :: ik,ikpj,il1,ip,n,i, ia
    real(kind=DP) :: dx, dy, dz, dist
    real(kind=DP) :: fac, facr, dnorm
    real(kind=DP), allocatable, dimension(:) :: qx,qy,qz,vlength,wka,wkb
    real(kind=DP), allocatable, dimension(:) :: zcos,zsin
    real(kind=DP), allocatable, dimension(:,:) :: pf
    real(kind=DP), allocatable, dimension(:,:) :: angf
    real(kind=DP), allocatable, dimension(:) :: radr,wos,radfunc
    integer, allocatable :: lmin_pj(:) ! d(n_proj)
    integer, allocatable :: lmax_pj(:) ! d(n_proj)

    integer             :: id_sname = -1
    call tstatc0_begin('m_Wan90_Projectors ',id_sname,1)

    allocate(lmin_pj(n_proj))
    allocate(lmax_pj(n_proj))
    allocate(qx(kg1),qy(kg1),qz(kg1),vlength(kg1),wka(kg1),wkb(kg1))
    allocate(pf(kg1,lmax1_pj),angf(kg1,lmax1_pj))
    allocate(radr(nmesh),wos(nmesh),radfunc(nmesh))
    allocate(zcos(kg1),zsin(kg1))

    allocate( centre_on_atom(n_proj) ); centre_on_atom = .false.
    Do n=1, n_proj
       Loop_ia: Do ia=1, natm
          dx = centre(1,n) -pos(ia,1)
          dy = centre(2,n) -pos(ia,2)
          dz = centre(3,n) -pos(ia,3)
          dist = sqrt( dx**2 +dy**2 +dz**2 )
          if ( dist < 1.0E-3 ) then
             centre_on_atom(n) = .true.
             exit Loop_ia
          endif
       End do Loop_ia

       if ( sw_use_hardpart_wan90 == ON ) then
          if ( .not. centre_on_atom(n) ) then
             write(nfout,*) "** Warning : proj ",n
             write(nfout,*) "** In the case of sw_use_hardpart_wan90 = ON "
             write(nfout,*) "** We recommend that Localized orbitals are on atomic sites"
!             stop
          endif
       End if
    End Do

    do ip=1,n_proj
       call set_lmin_lmax(lang(ip),lmin_pj(ip),lmax_pj(ip))
       if(ipri>=2) write(*,'("ip=",i5," lmin=",i5," lmax=",i5)') ip, lmin_pj(ip), lmax_pj(ip)
    end do

    call radr_and_wos(nmesh,radr,wos) ! --> radr, wos
    fac = PAI4/dsqrt(univol)
    do ik = 1, kv3, nspin
       if(map_k(ik) /= myrank_k) cycle                     ! MPI
       call k_plus_G_vectors_3D(ik,kgp,kg1,kv3,iba,nbase,vkxyz,ngabc,rltv&
                           &,qx,qy,qz,vlength) ! ->(bottom_Subr.)
!!       call G_vectors(kgp,iba(ik),nbase(1,ik),ngabc,rltv,qx,qy,qz,vlength)
       ikpj = (ik-1)/nspin + 1
       do ip=1,n_proj
          if(ipri>=2) write(*,'("ikpj=",i5," ip=",i5)') ikpj,ip

          call radial_function(irf(ip),zona(ip),nmesh,radr,radfunc)
! Debug
         !! if(ik==1) then
         !! write(1000+ip,'("radr radfunc")')
         !! do n=1,nmesh
         !!    write(1000+ip,'(f20.8,1x,f20.8)') radr(n), radfunc(n)
         !! end do
         !! end if
! end Debug
          if(ipri>=2) write(*,'("debug 1")')
          angf = 0.d0
          do il1=lmin_pj(ip)+1,lmax_pj(ip)+1
             call angular_function(lang(ip),mr(ip),iba(ik),qx,qy,qz,angf)
          end do
          if(ipri>=2) write(*,'("debug 2")')
          pf = 0.d0
          do n=1,nmesh
             facr = fac*wos(n)*radr(n)**2*radfunc(n)
             do i=1,iba(ik)
                wka(i) = vlength(i)*radr(n)
             end do
             do il1=lmin_pj(ip)+1,lmax_pj(ip)+1
                call dsjnv(il1-1,iba(ik),wka,wkb)     ! -(bottom_Subr.)
                do i=1,iba(ik)
                   pf(i,il1) = pf(i,il1) + facr*wkb(i)*angf(i,il1)
                end do
             end do
          end do
          if(ipri>=2) write(*,'("debug 3")')

          !! Projector(k+G) = sum_l i^(-l) * exp(-i(k+G)*R) * PF_l(k+G)
          call exp_KpG_dot_R(centre(1,ip),ik,kgp,iba,nbase,ngabc,vkxyz,zcos,zsin)

          if(ipri>=2) write(*,'("debug 4")')
          projfunc(1:kg1,ip,ikpj,1:kimg) = 0.d0
          do il1=lmin_pj(ip)+1,lmax_pj(ip)+1
             if(il1==1) then
                if(kimg==1) then
                   do i=1,iba(ik)
                      projfunc(i,ip,ikpj,1) = projfunc(i,ip,ikpj,1) + pf(i,il1) * zcos(i)
                   end do
                else
                   do i=1,iba(ik)
                      projfunc(i,ip,ikpj,1) = projfunc(i,ip,ikpj,1) + pf(i,il1) * zcos(i)
                      projfunc(i,ip,ikpj,2) = projfunc(i,ip,ikpj,2) - pf(i,il1) * zsin(i)
                   end do
                end if
             else if(il1==2) then
                if(kimg==1) then
                   do i=1,iba(ik)
                      projfunc(i,ip,ikpj,1) = projfunc(i,ip,ikpj,1) - pf(i,il1) * zsin(i)
                   end do
                else
                   do i=1,iba(ik)
                      projfunc(i,ip,ikpj,1) = projfunc(i,ip,ikpj,1) - pf(i,il1) * zsin(i)
                      projfunc(i,ip,ikpj,2) = projfunc(i,ip,ikpj,2) - pf(i,il1) * zcos(i)
                   end do
                end if
             else if(il1==3) then
                if(kimg==1) then
                   do i=1,iba(ik)
                      projfunc(i,ip,ikpj,1) = projfunc(i,ip,ikpj,1) - pf(i,il1) * zcos(i)
                   end do
                else
                   do i=1,iba(ik)
                      projfunc(i,ip,ikpj,1) = projfunc(i,ip,ikpj,1) - pf(i,il1) * zcos(i)
                      projfunc(i,ip,ikpj,2) = projfunc(i,ip,ikpj,2) + pf(i,il1) * zsin(i)
                   end do
                end if
             else if(il1==4) then
                if(kimg==1) then
                   do i=1,iba(ik)
                      projfunc(i,ip,ikpj,1) = projfunc(i,ip,ikpj,1) + pf(i,il1) * zsin(i)
                   end do
                else
                   do i=1,iba(ik)
                      projfunc(i,ip,ikpj,1) = projfunc(i,ip,ikpj,1) + pf(i,il1) * zsin(i)
                      projfunc(i,ip,ikpj,2) = projfunc(i,ip,ikpj,2) + pf(i,il1) * zcos(i)
                   end do
                end if
             end if
          end do
          if(ipri>=2) write(*,'("debug 5")')
! === KT_add ==== 2015/04/13
!          if ( sw_use_hardpart_wan90 == ON ) goto 1000
          if ( sw_use_hardpart_wan90 == ON .and. centre_on_atom(ip) ) goto 1000
! =============== 2015/04/13
#if 1
! Normalization
          dnorm = 0.d0
          do i=1,iba(ik)
             dnorm = dnorm + projfunc(i,ip,ikpj,1)**2 + projfunc(i,ip,ikpj,2)**2
          end do
! Debug
!!          write(nfout,'("ip=",i5," ikpj=",i5," dnorm=",f25.10," centre=",3(1x,f10.5))') &
!!           & ip,ikpj,dnorm,centre(1:3,ip)
! End Debug
          dnorm = 1.d0/sqrt(dnorm)
          if(kimg==1) then
             do i=1,iba(ik)
                projfunc(i,ip,ikpj,1) = projfunc(i,ip,ikpj,1) * dnorm
             end do
          else
             do i=1,iba(ik)
                projfunc(i,ip,ikpj,1) = projfunc(i,ip,ikpj,1) * dnorm
                projfunc(i,ip,ikpj,2) = projfunc(i,ip,ikpj,2) * dnorm
             end do
          end if
#endif
! === KT_add ==== 2015/02/28
1000      continue
! =============== 2015/02/28
      end do
    end do

    deallocate(lmin_pj)
    deallocate(lmax_pj)
    deallocate(qx,qy,qz,vlength,wka,wkb)
    deallocate(pf,angf)
    deallocate(zcos,zsin)
    deallocate(radr,wos,radfunc)

    call tstatc0_end(id_sname)
  end subroutine m_Wan90_Projectors

  subroutine radr_and_wos(nm,radr,wos)
    integer, intent(in) :: nm
    real(kind=DP), intent(out) :: radr(nm), wos(nm)

    real(kind=DP), parameter :: xh = 96.d0
    real(kind=DP), parameter :: rmax = 60.d0
    real(kind=DP) :: hn

    call rmeshs(nm,nm,xh,rmax,radr,hn) ! -(b_PP)
    call coef_simpson_integration(nm,nm,xh,radr,wos) ! -(b_PP)
  end subroutine radr_and_wos

  subroutine spherical_function(l,m,nq,qx,qy,qz,sphfunc)
    integer, intent(in) :: l,m,nq
    real(kind=DP), intent(in) :: qx(nq), qy(nq), qz(nq)
    real(kind=DP), intent(out) :: sphfunc(nq)

    if(l==0.and.m==1) then ! s
       call sphr(nq,1,qx,qy,qz,sphfunc)        ! -(bottom_Subr.)
    else if(l==1.and.m==1) then ! pz
       call sphr(nq,4,qx,qy,qz,sphfunc)        ! -(bottom_Subr.)
    else if(l==1.and.m==2) then ! px
       call sphr(nq,2,qx,qy,qz,sphfunc)        ! -(bottom_Subr.)
    else if(l==1.and.m==3) then ! py
       call sphr(nq,3,qx,qy,qz,sphfunc)        ! -(bottom_Subr.)
    else if(l==2.and.m==1) then ! dz2
       call sphr(nq,5,qx,qy,qz,sphfunc)        ! -(bottom_Subr.)
    else if(l==2.and.m==2) then ! dxz
       call sphr(nq,9,qx,qy,qz,sphfunc)        ! -(bottom_Subr.)
    else if(l==2.and.m==3) then ! dyz
       call sphr(nq,8,qx,qy,qz,sphfunc)        ! -(bottom_Subr.)
    else if(l==2.and.m==4) then ! dx2-y2
       call sphr(nq,6,qx,qy,qz,sphfunc)        ! -(bottom_Subr.)
    else if(l==2.and.m==5) then ! dxy
       call sphr(nq,7,qx,qy,qz,sphfunc)        ! -(bottom_Subr.)
    else if(l==3.and.m==1) then ! fz3
       call sphr(nq,10,qx,qy,qz,sphfunc)        ! -(bottom_Subr.)
    else if(l==3.and.m==2) then ! fxz2
       call sphr(nq,11,qx,qy,qz,sphfunc)        ! -(bottom_Subr.)
    else if(l==3.and.m==3) then ! fyz2
       call sphr(nq,12,qx,qy,qz,sphfunc)        ! -(bottom_Subr.)
    else if(l==3.and.m==4) then ! fz(x2-y2)
       call sphr(nq,13,qx,qy,qz,sphfunc)        ! -(bottom_Subr.)
    else if(l==3.and.m==5) then ! fxyz
       call sphr(nq,14,qx,qy,qz,sphfunc)        ! -(bottom_Subr.)
    else if(l==3.and.m==6) then ! fx(x2-3y2)
       call sphr(nq,15,qx,qy,qz,sphfunc)        ! -(bottom_Subr.)
    else if(l==3.and.m==7) then ! fy(3x2-y2)
       call sphr(nq,16,qx,qy,qz,sphfunc)        ! -(bottom_Subr.)
    end if
  end subroutine spherical_function

  subroutine angular_function(l,m,nq,qx,qy,qz,angfunc)
    integer, intent(in) :: l,m,nq
    real(kind=DP), intent(in) :: qx(nq), qy(nq), qz(nq)
    real(kind=DP), intent(out) :: angfunc(kg1,lmax1_pj)

    real(kind=DP), allocatable :: sphfunc(:)

    allocate(sphfunc(nq))

    if(l>=0) then
      call spherical_function(l,m,nq,qx,qy,qz,sphfunc)
      angfunc(1:nq,1) = sphfunc(1:nq)
    else if(l==-1.and.m==1) then
      call spherical_function(0,1,nq,qx,qy,qz,sphfunc) ! s
      angfunc(1:nq,1) = sphfunc(1:nq)/sqrt(2.d0)
      call spherical_function(1,2,nq,qx,qy,qz,sphfunc) ! px
      angfunc(1:nq,2) = sphfunc(1:nq)/sqrt(2.d0)
    else if(l==-1.and.m==2) then
      call spherical_function(0,1,nq,qx,qy,qz,sphfunc) ! s
      angfunc(1:nq,1) = sphfunc(1:nq)/sqrt(2.d0)
      call spherical_function(1,2,nq,qx,qy,qz,sphfunc) ! px
      angfunc(1:nq,2) = -sphfunc(1:nq)/sqrt(2.d0)
    else if(l==-2.and.m==1) then
      call spherical_function(0,1,nq,qx,qy,qz,sphfunc) ! s
      angfunc(1:nq,1) = sphfunc(1:nq)/sqrt(3.d0)
      call spherical_function(1,2,nq,qx,qy,qz,sphfunc) ! px
      angfunc(1:nq,2) = -sphfunc(1:nq)/sqrt(6.d0)
      call spherical_function(1,3,nq,qx,qy,qz,sphfunc) ! py
      angfunc(1:nq,2) = angfunc(1:nq,2) + sphfunc(1:nq)/sqrt(2.d0)
    else if(l==-2.and.m==2) then
      call spherical_function(0,1,nq,qx,qy,qz,sphfunc) ! s
      angfunc(1:nq,1) = sphfunc(1:nq)/sqrt(3.d0)
      call spherical_function(1,2,nq,qx,qy,qz,sphfunc) ! px
      angfunc(1:nq,2) = -sphfunc(1:nq)/sqrt(6.d0)
      call spherical_function(1,3,nq,qx,qy,qz,sphfunc) ! py
      angfunc(1:nq,2) = angfunc(1:nq,2) - sphfunc(1:nq)/sqrt(2.d0)
    else if(l==-2.and.m==3) then
      call spherical_function(0,1,nq,qx,qy,qz,sphfunc) ! s
      angfunc(1:nq,1) = sphfunc(1:nq)/sqrt(3.d0)
      call spherical_function(1,2,nq,qx,qy,qz,sphfunc) ! px
      angfunc(1:nq,2) = sphfunc(1:nq)*(2.d0/sqrt(6.d0))
    else if(l==-3.and.m==1) then
      call spherical_function(0,1,nq,qx,qy,qz,sphfunc) ! s
      angfunc(1:nq,1) = sphfunc(1:nq)*0.5d0
      call spherical_function(1,2,nq,qx,qy,qz,sphfunc) ! px
      angfunc(1:nq,2) = sphfunc(1:nq)*0.5d0
      call spherical_function(1,3,nq,qx,qy,qz,sphfunc) ! py
      angfunc(1:nq,2) = angfunc(1:nq,2) + sphfunc(1:nq)*0.5d0
      call spherical_function(1,1,nq,qx,qy,qz,sphfunc) ! pz
      angfunc(1:nq,2) = angfunc(1:nq,2) + sphfunc(1:nq)*0.5d0
    else if(l==-3.and.m==2) then
      call spherical_function(0,1,nq,qx,qy,qz,sphfunc) ! s
      angfunc(1:nq,1) = sphfunc(1:nq)*0.5d0
      call spherical_function(1,2,nq,qx,qy,qz,sphfunc) ! px
      angfunc(1:nq,2) = sphfunc(1:nq)*0.5d0
      call spherical_function(1,3,nq,qx,qy,qz,sphfunc) ! py
      angfunc(1:nq,2) = angfunc(1:nq,2) - sphfunc(1:nq)*0.5d0
      call spherical_function(1,1,nq,qx,qy,qz,sphfunc) ! pz
      angfunc(1:nq,2) = angfunc(1:nq,2) - sphfunc(1:nq)*0.5d0
    else if(l==-3.and.m==3) then
      call spherical_function(0,1,nq,qx,qy,qz,sphfunc) ! s
      angfunc(1:nq,1) = sphfunc(1:nq)*0.5d0
      call spherical_function(1,2,nq,qx,qy,qz,sphfunc) ! px
      angfunc(1:nq,2) = -sphfunc(1:nq)*0.5d0
      call spherical_function(1,3,nq,qx,qy,qz,sphfunc) ! py
      angfunc(1:nq,2) = angfunc(1:nq,2) + sphfunc(1:nq)*0.5d0
      call spherical_function(1,1,nq,qx,qy,qz,sphfunc) ! pz
      angfunc(1:nq,2) = angfunc(1:nq,2) - sphfunc(1:nq)*0.5d0
    else if(l==-3.and.m==4) then
      call spherical_function(0,1,nq,qx,qy,qz,sphfunc) ! s
      angfunc(1:nq,1) = sphfunc(1:nq)*0.5d0
      call spherical_function(1,2,nq,qx,qy,qz,sphfunc) ! px
      angfunc(1:nq,2) = -sphfunc(1:nq)*0.5d0
      call spherical_function(1,3,nq,qx,qy,qz,sphfunc) ! py
      angfunc(1:nq,2) = angfunc(1:nq,2) - sphfunc(1:nq)*0.5d0
      call spherical_function(1,1,nq,qx,qy,qz,sphfunc) ! pz
      angfunc(1:nq,2) = angfunc(1:nq,2) + sphfunc(1:nq)*0.5d0
    else if(l<=-4) then
      call phase_error_with_msg(nfout,'Not supported for l<=-4',__LINE__,__FILE__)
    end if

    deallocate(sphfunc)
  end subroutine angular_function

  subroutine radial_function(irf,zona,nm,radr,radfunc)
    implicit none
    integer, intent(in) :: irf, nm
    real(kind=DP) :: zona
    real(kind=DP) :: radr(nm), radfunc(nm)

    integer :: i
    real(kind=DP) :: fac, ar, alph

    alph = zona * bohr

    if(irf==1) then
      fac = 2.d0 * alph**(3.d0/2.d0)
      do i=1,nm
         radfunc(i) = fac * exp(-alph*radr(i))
      end do
    else if(irf==2) then
      fac = alph**(3.d0/2.d0) / sqrt(8.d0)
      do i=1,nm
         ar = alph * radr(i)
         radfunc(i) = fac * ( 2.d0 - ar ) * exp(-0.5d0*ar)
      end do
    else if(irf==3) then
      fac = alph**(3.d0/2.d0) * sqrt(4.d0/27.d0)
      do i=1,nm
         ar = alph * radr(i)
         radfunc(i) = fac * ( 1.d0 - (2.d0/3.d0)*ar + (2.d0/27.d0)*ar*ar ) * exp(-ar/3.d0)
      end do
    end if
  end subroutine radial_function

  subroutine set_lmin_lmax(l,lmin_pj,lmax_pj)
    implicit none
    integer, intent(in) :: l
    integer, intent(out) :: lmin_pj, lmax_pj

    if(l>=0) then
       lmin_pj = 0
       lmax_pj = 0
    else if(l==-1.or.l==-2.or.l==-3) then
       lmin_pj = 0
       lmax_pj = 1
    else if(l<=-4) then
       call phase_error_with_msg(nfout,'Not supported for l<=-4',__LINE__,__FILE__)
    end if
  end subroutine set_lmin_lmax

  subroutine exp_KpG_dot_R(centre,ik,kgp,iba,nbase,ngabc,vkxyz,zcos,zsin)
    implicit none
    real(kind=DP), intent(in) :: centre(3)
    integer, intent(in) :: ik, kgp
    integer, intent(in) :: iba(kv3), nbase(kg1,kv3), ngabc(kgp,3)
    real(kind=DP), intent(in) :: vkxyz(kv3,3,CRDTYP)
    real(kind=DP), intent(out) :: zcos(kg1), zsin(kg1)

    integer :: i
    real(kind=DP) :: ph

    do i=1,iba(ik)
       ph = PAI2 * dot_product(centre,vkxyz(ik,1:3,BUCS)+ngabc(nbase(i,ik),1:3))
       zcos(i) = cos(ph)
       zsin(i) = sin(ph)
    end do
  end subroutine exp_KpG_dot_R

  subroutine m_Wan90_wd_eig(nfout)
    implicit none
    integer, intent(in) :: nfout

    integer :: ik,m,ib,mi
    real(kind=DP) :: ek

    !!$if(npes /= 1) return
    if(mype==0) then
      open(nfwannier,file=trim(wan90_seedname)//".eig",form="formatted")
      do ik = 1, kv3
         do m=1,num_bands
            mi = ib_inc(m)
            ib = neordr(mi,ik)
            ek = eko_l(map_z(ib),ik) * hartree
            write(nfwannier,'(2(1x,i5),1x,e25.16)') m, ik, ek
         end do
      end do
      close(nfwannier)
   end if
  end subroutine m_Wan90_wd_eig

  subroutine m_Wan90_gen_mmat(nfout)
    use m_Parallelization, only : nel_fft_x, nel_fft_y, nel_fft_z
    implicit none
    integer, intent(in) :: nfout

    integer :: i,ik1,ik2,nn,n,m,ni,mi,ib1,ib2,nnc(3),ir,im
    integer :: nffth,ip,r1,r2,r3,ngrid
    integer :: id1,id2,id12
    real(kind=DP) :: rgrid(3),da(3),wr,wi,dvol,ph,sumr,sumi
    real(kind=DP), allocatable :: wf1(:,:), wf2(:,:), zcos(:), zsin(:)
    real(kind=DP), allocatable :: m_mat(:,:,:,:,:) ! d(neg,neg,kv3,nntot,2)
    real(kind=DP), allocatable :: m_mat_mpi(:,:,:,:,:) ! d(neg,neg,kv3,nntot,2)
    logical :: shift_k
    integer :: ibsize = 1, lsize
    lsize = max(maxval(nel_fft_x(:)),maxval(nel_fft_y(:)),maxval(nel_fft_z(:)))

    nffth = nfft/2
    ngrid = product(fft_box_size_WF(1:3,1))
    dvol = 1.d0/dble(ngrid)

    allocate(m_mat(neg,neg,kv3,nntot,2)); m_mat = 0.d0
    allocate(wf1(lsize*kimg,ibsize),wf2(lsize*kimg,ibsize))
    allocate(zcos(nffth),zsin(nffth))

    call m_FFT_alloc_WF_work()

    do nn=1,nntot
       do ik1=1,kv3
          ik2 = nnlist(ik1,nn)
          nnc(1:3) = nncell(1:3,ik1,nn)
          if(all(nnc(1:3) == 0)) then
             shift_k = .false.
          else
             shift_k = .true.
             ! Phase factors, exp(iG*r)
             id1 = fft_box_size_WF(1,0)
             id2 = fft_box_size_WF(2,0)
             id12 = id1*id2
             da(1:3) = 1.d0/fft_box_size_WF(1:3,1)
             do i=1,nffth
                ip = i-1
                r3 = ip/id12
                r2 = (ip-r3*id12)/id1
                r1 = ip-r2*id1-r3*id12
                rgrid(1) = dble(r1)*da(1)
                rgrid(2) = dble(r2)*da(2)
                rgrid(3) = dble(r3)*da(3)
                ph = PAI2*dot_product(nnc,rgrid)
                zcos(i) = cos(ph)
                zsin(i) = sin(ph)
             end do
          end if
          do n=1,num_bands
             ni = ib_inc(n)
             ib2 = neordr(ni,ik2)
             call m_ES_WF_in_Rspace_3D(ik2,ib2,ib2,ibsize,lsize,wf2) !-> u_n(k+b) or u_n(k')
             if(shift_k) then
                ! u_n(k+b) = exp(-iG*r) * u_n(k')
                ! k+b = k' + G
                do i=1,nffth
                   ir = 2*i-1
                   im = ir+1
                   wr = wf2(ir,1)
                   wi = wf2(im,1)
                   wf2(ir,1) = zcos(i)*wr + zsin(i)*wi
                   wf2(im,1) = zcos(i)*wi - zsin(i)*wr
                end do
             end if
             do m=1,num_bands
                mi = ib_inc(m)
                ib1 = neordr(mi,ik1)
                call m_ES_WF_in_Rspace_3D(ik1,ib1,ib1,ibsize,lsize,wf1) !-> u_m(k)
                ! <u_m(k)|u_n(k+b)>
                sumr = 0.d0
                sumi = 0.d0
                do i=1,nffth
                   ir = 2*i-1
                   im = ir+1
                   sumr = sumr + wf1(ir,1)*wf2(ir,1) + wf1(im,1)*wf2(im,1)
                   sumi = sumi + wf1(ir,1)*wf2(im,1) - wf1(im,1)*wf2(ir,1)
                end do
                m_mat(m,n,ik1,nn,1) = sumr * dvol
                m_mat(m,n,ik1,nn,2) = sumi * dvol
             end do
          end do
       end do
    end do

    call m_FFT_dealloc_WF_work()

    deallocate(wf1,wf2)
    deallocate(zcos,zsin)

    if(npes>1) then
       allocate(m_mat_mpi(num_bands,num_bands,kv3,nntot,2))
       m_mat_mpi = m_mat
       m_mat = 0.d0
       call mpi_allreduce(m_mat_mpi,m_mat,num_bands*num_bands*kv3*nntot*2,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
       deallocate(m_mat_mpi)
    end if

    if(mype==0) then
       open(nfwannier,file=trim(wan90_seedname)//".mmn",form="formatted")
       write(nfwannier,*) 'Generated by PHASE'
       write(nfwannier,*) num_bands, kv3, nntot
       do ik1 = 1, kv3
          do nn=1,nntot
             ik2 = nnlist(ik1,nn)
             nnc(1:3) = nncell(1:3,ik1,nn)
             write(nfwannier,'(5(1x,i5))') ik1,ik2,nnc(1:3)
             do n=1,num_bands
                do m=1,num_bands
                   write(nfwannier,'(2(1x,f16.12))') m_mat(m,n,ik1,nn,1:2)
                end do
             end do
          end do
       end do
       close(nfwannier)
    end if

    deallocate(m_mat)
  end subroutine m_Wan90_gen_mmat

  subroutine m_Wan90_write_unk(nfout)
    use m_Parallelization, only : nel_fft_x, nel_fft_y, nel_fft_z
    implicit none
    integer, intent(in) :: nfout

    integer :: ispin, iksnl, ik, n, nffth, ngrid, ni, ib
    integer :: id1, id2, id3, id12, ip
    integer :: nd1, nd2, nd3
    integer :: i1, i2, i3
    character(len=10) :: wfnname
    real(kind=DP), allocatable :: wf(:,:)
    complex(kind=DP), allocatable :: zwf(:,:,:)
    integer :: ibsize = 1, lsize
    lsize = max(maxval(nel_fft_x(:)),maxval(nel_fft_y(:)),maxval(nel_fft_z(:)))

    id1 = fft_box_size_WF(1,0)
    id2 = fft_box_size_WF(2,0)
    id3 = fft_box_size_WF(3,0)
    id12 = id1 * id2

    nd1 = fft_box_size_WF(1,1)
    nd2 = fft_box_size_WF(2,1)
    nd3 = fft_box_size_WF(3,1)

    allocate(wf(lsize*kimg,ibsize))
    allocate(zwf(id1,id2,id3))

    if(mype==0) then
       do ispin = 1, nspin
          do iksnl = 1, kv3/nspin
             ik = nspin*(iksnl-1) + ispin
             write(wfnname,'("UNK",i5.5,".",i1)') iksnl,ispin
             open(nfwannier,file=trim(wfnname),form="unformatted")
             write(nfwannier) fft_box_size_WF(1:3,1), iksnl, num_bands
             do n=1,num_bands
                ni = ib_inc(n)
                ib = neordr(ni,ik)
                call m_ES_WF_in_Rspace_3D(ik,ib,ib,ibsize,lsize,wf)
                do i3=1,nd3
                   do i2=1,nd2
                      do i1=1,nd1
                         ip = i1 + (i2-1)*id1 + (i3-1) * id12
                         zwf(i1,i2,i3) = cmplx(wf(2*ip-1,1),wf(2*ip,1))
                      end do
                   end do
                end do
                write(nfwannier) zwf
             end do
             close(nfwannier)
          end do
       end do
    end if

    deallocate(wf)
    deallocate(zwf)

  end subroutine m_Wan90_write_unk


  logical function excluding_band(m)
    implicit none
    integer, intent(in) :: m

    integer :: i

    do i=1,n_exclude
       if(exclude_bands(i)==m) then
          excluding_band = .true.
          return
       end if
    end do
    excluding_band = .false.
  end function excluding_band

! ==== KT_add === 2015/02/23
  subroutine m_Wan90_set_dk_wan
    integer :: ik1, ik2, nn, nnc(3)
    real(kind=DP) :: vec(3)

    allocate( dk_unit(3,nntot) );  dk_unit = 0.0d0

    ik1 = 1
    Do nn=1, nntot
       ik2 = nspin*( nnlist(ik1,nn) -1 )+1
       nnc(1:3) = nncell(1:3,ik1,nn)

       dk_unit(1:3,nn) = vkxyz(ik2,1:3,BUCS) +nnc(1:3) -vkxyz(ik1,1:3,BUCS)
!
       vec(1) = rltv(1,1) *dk_unit(1,nn) &
            & + rltv(1,2) *dk_unit(2,nn) &
            & + rltv(1,3) *dk_unit(3,nn)
       vec(2) = rltv(2,1) *dk_unit(1,nn) &
            & + rltv(2,2) *dk_unit(2,nn) &
            & + rltv(2,3) *dk_unit(3,nn)
       vec(3) = rltv(3,1) *dk_unit(1,nn) &
            & + rltv(3,2) *dk_unit(2,nn) &
            & + rltv(3,3) *dk_unit(3,nn)
       dk_wan(nn) = sqrt( dot_product( vec(:), vec(:) ) )
    End Do
  end subroutine m_Wan90_set_dk_wan

  subroutine deompose_sphr_into_lm( l, m, vec )
    integer, intent(in) :: l, m
    real(kind=DP), intent(out) :: vec(16)
!
!    assuming xaxis and azxis  are default
!
    vec = 0.0d0

    if ( l == 0 ) then
       if ( m == 1 ) vec(1) = 1.0d0   ! s
    endif
    if ( l == 1 ) then
       if ( m == 1 ) vec(4) = 1.0d0   ! pz
       if ( m == 2 ) vec(2) = 1.0d0   ! px
       if ( m == 3 ) vec(3) = 1.0d0   ! py
    endif
    if ( l == 2 ) then
       if ( m == 1 ) vec(5) = 1.0d0   ! dz2
       if ( m == 2 ) vec(9) = 1.0d0   ! dzx
       if ( m == 3 ) vec(8) = 1.0d0   ! dyz
       if ( m == 4 ) vec(6) = 1.0d0   ! dx2-y2
       if ( m == 5 ) vec(7) = 1.0d0   ! dxy
    endif
    if ( l == 3 ) then
       if ( m == 1 ) vec(10) = 1.0d0   ! fz3
       if ( m == 2 ) vec(11) = 1.0d0   ! fxz2
       if ( m == 3 ) vec(12) = 1.0d0   ! fyz2
       if ( m == 4 ) vec(13) = 1.0d0   ! fz(x2-y2)
       if ( m == 5 ) vec(14) = 1.0d0   ! fxyz
       if ( m == 6 ) vec(15) = 1.0d0   ! fx(x2-3y2)
       if ( m == 7 ) vec(16) = 1.0d0   ! fy(3x2-y2)
    endif
    if ( l == -1 ) then
       if ( m == 1 ) then            ! sp-1
          vec(1) = 1.0d0 /sqrt(2.0d0)   ! s
          vec(2) = vec(1)               ! px
       endif
       if ( m == 2 ) then            ! sp-2
          vec(1) = 1.0d0 /sqrt(2.0d0)   ! s
          vec(2) = -vec(1)              ! px
       endif
    endif
    if ( l == -2 ) then
       if ( m == 1 ) then            ! sp2-1
          vec(1) =  1.0d0 /sqrt(3.0d0)   ! s
          vec(2) = -1.0d0 /sqrt(6.0d0)   ! px
          vec(3) =  1.0d0 /sqrt(2.0d0)   ! py
       endif
       if ( m == 2 ) then            ! sp2-2
          vec(1) =  1.0d0 /sqrt(3.0d0)   ! s
          vec(2) = -1.0d0 /sqrt(6.0d0)   ! px
          vec(3) = -1.0d0 /sqrt(2.0d0)   ! py
       endif
       if ( m == 3 ) then            ! sp2-3
          vec(1) =  1.0d0 /sqrt(3.0d0)   ! s
          vec(2) =  2.0d0 /sqrt(6.0d0)   ! px
       endif
    endif
    if ( l == -3 ) then
       if ( m == 1 ) then     ! sp3-1
          vec(1) =  1.0d0 /2.0d0    ! s
          vec(2) =  1.0d0 /2.0d0    ! px
          vec(3) =  1.0d0 /2.0d0    ! py
          vec(4) =  1.0d0 /2.0d0    ! pz
       endif
       if ( m == 2 ) then     ! sp3-2
          vec(1) =  1.0d0 /2.0d0    ! s
          vec(2) =  1.0d0 /2.0d0    ! px
          vec(3) = -1.0d0 /2.0d0    ! py
          vec(4) = -1.0d0 /2.0d0    ! pz
       endif
       if ( m == 3 ) then     ! sp3-3
          vec(1) =  1.0d0 /2.0d0    ! s
          vec(2) = -1.0d0 /2.0d0    ! px
          vec(3) =  1.0d0 /2.0d0    ! py
          vec(4) = -1.0d0 /2.0d0    ! pz
       endif
       if ( m == 4 ) then     ! sp3-4
          vec(1) =  1.0d0 /2.0d0    ! s
          vec(2) = -1.0d0 /2.0d0    ! px
          vec(3) = -1.0d0 /2.0d0    ! py
          vec(4) =  1.0d0 /2.0d0    ! pz
       endif
    endif
    if ( l == -4 ) then
       if ( m == 1 ) then      ! sp3d-1
          vec(1) =  1.0d0 /sqrt(3.0d0)
          vec(2) = -1.0d0 /sqrt(6.0d0)
          vec(3) =  1.0d0 /sqrt(2.0d0)
       endif
       if ( m == 2 ) then      ! sp3d-2
          vec(1) =  1.0d0 /sqrt(3.0d0)
          vec(2) = -1.0d0 /sqrt(6.0d0)
          vec(3) = -1.0d0 /sqrt(2.0d0)
       endif
       if ( m == 3 ) then      ! sp3d-3
          vec(1) =  1.0d0 /sqrt(3.0d0)
          vec(2) =  2.0d0 /sqrt(6.0d0)
       endif
       if ( m == 4 ) then      ! sp3d-4
          vec(4) =  1.0d0 /sqrt(2.0d0)
          vec(5) =  1.0d0 /sqrt(2.0d0)
       endif
       if ( m == 5 ) then      ! sp3d-5
          vec(4) = -1.0d0 /sqrt(2.0d0)
          vec(5) =  1.0d0 /sqrt(2.0d0)
       endif
    endif
    if ( l == -5 ) then
       if ( m == 1 ) then      ! sp3d2-1
          vec(1) =  1.0d0 /sqrt(6.0d0)
          vec(2) = -1.0d0 /sqrt(2.0d0)
          vec(5) = -1.0d0 /sqrt(12.0d0)
          vec(6) =  1.0d0 /2.0d0
       endif
       if ( m == 2 ) then      ! sp3d2-2
          vec(1) =  1.0d0 /sqrt(6.0d0)
          vec(2) =  1.0d0 /sqrt(2.0d0)
          vec(5) = -1.0d0 /sqrt(12.0d0)
          vec(6) =  1.0d0 /2.0d0
       endif
       if ( m == 3 ) then      ! sp3d2-3
          vec(1) =  1.0d0 /sqrt(6.0d0)
          vec(2) = -1.0d0 /sqrt(2.0d0)
          vec(5) = -1.0d0 /sqrt(12.0d0)
          vec(6) = -1.0d0 /2.0d0
       endif
       if ( m == 4 ) then      ! sp3d2-4
          vec(1) =  1.0d0 /sqrt(6.0d0)
          vec(2) =  1.0d0 /sqrt(2.0d0)
          vec(5) = -1.0d0 /sqrt(12.0d0)
          vec(6) = -1.0d0 /2.0d0
       endif
       if ( m == 5 ) then      ! sp3d2-5
          vec(1) =  1.0d0 /sqrt(6.0d0)
          vec(4) = -1.0d0 /sqrt(2.0d0)
          vec(5) =  1.0d0 /sqrt(3.0d0)
       endif
       if ( m == 6 ) then      ! sp3d2-6
          vec(1) =  1.0d0 /sqrt(6.0d0)
          vec(4) =  1.0d0 /sqrt(2.0d0)
          vec(5) =  1.0d0 /sqrt(3.0d0)
       endif
    endif
  end subroutine deompose_sphr_into_lm
! =========== 2015/02/23

! ==== KT_add === 2015/09/04
  subroutine m_Wan90_gen_mat_spn
    integer :: npauli = 3

    complex(kind=CMPLDP), allocatable :: spn_mat(:,:,:,:) ! d(neg,neg,kv3/nspin,3)

    allocate( spn_mat( npauli, num_bands, num_bands, kv3/ndim_spinor) ); spn_mat = 0.d0

    call contrib_softpart
    if ( sw_use_hardpart_wan90 == ON ) call contrib_hardpart

    call gather_matrix
    call print_mat

    deallocate( spn_mat )

  contains

    subroutine contrib_softpart
      integer :: ik0, ik1, ik2, iktmp, is1, is2
      integer :: ig
      integer :: m, mi, n, ni, ib1, ib2
      complex(kind=CMPLDP) :: zsum

      real(kind=DP), allocatable :: wk_zaj(:,:)
      real(kind=DP) :: c1, c2
      complex(kind=CMPLDP) :: z1, z2

      allocate( wk_zaj( kg1, kimg ) ); wk_zaj = 0.0d0

      do ik0=1, kv3, ndim_spinor
         if ( map_k(ik0) /= myrank_k ) cycle

         iktmp = ( ik0 -1 )/ndim_spinor +1

         Do is2=1, ndim_spinor
            ik2 = ik0 +is2 -1

            Do n=1, num_bands
               ni = ib_inc(n)
               ib2 = neordr(ni,ik2)

               if ( map_e(ib2) == myrank_e ) then
                  wk_zaj(1:iba(ik2),1:kimg) = zaj_l(1:iba(ik2),map_z(ib2),ik2,1:kimg)
               endif
               call mpi_bcast( wk_zaj, kg1*kimg, mpi_double_precision, map_e(ib2), &
                    &          mpi_k_world(myrank_k), ierr )

               Do is1=1, ndim_spinor
                  ik1 = ik0 +is1 -1

                  do m=1,num_bands
                     mi = ib_inc(m)
                     ib1 = neordr(mi,ik1)

                     if ( map_e(ib1) /= myrank_e ) cycle

                     zsum = 0.0d0

                     if ( kimg == 1 ) then
                        Do ig=1, iba(ik0)
                           c1 = zaj_l(ig,map_z(ib1),ik1,1)
                           c2 = c1 *wk_zaj(ig,1)
                           zsum = zsum +c2
                        End do
                     else
                        Do ig=1, iba(ik0)
                           z1 = dcmplx( zaj_l(ig,map_z(ib1),ik1,1), &
                                &       zaj_l(ig,map_z(ib1),ik1,2) )
                           z2 = conjg(z1) *dcmplx( wk_zaj(ig,1), wk_zaj(ig,2) )
                           zsum = zsum +z2
                        End Do
                     endif
!
                     if ( is1 == 1 .and. is2 == 1 ) then
                        spn_mat(3,m,n,iktmp) = spn_mat(3,m,n,iktmp) +zsum
                     else if ( is1 == 1 .and. is2 == 2 ) then
                        spn_mat(1,m,n,iktmp) = spn_mat(1,m,n,iktmp) +zsum
                        spn_mat(2,m,n,iktmp) = spn_mat(2,m,n,iktmp) -zsum *zi
                     else if ( is1 == 2 .and. is2 == 1 ) then
                        spn_mat(1,m,n,iktmp) = spn_mat(1,m,n,iktmp) +zsum
                        spn_mat(2,m,n,iktmp) = spn_mat(2,m,n,iktmp) +zsum *zi
                     else if ( is1 == 2 .and. is2 == 2 ) then
                        spn_mat(3,m,n,iktmp) = spn_mat(3,m,n,iktmp) -zsum
                     endif

                  end do
               end Do
            end do
         end do
      end do

    end subroutine contrib_softpart

    subroutine contrib_hardpart
      integer :: ik0, ik1, ik2, iktmp, is1, is2
      integer :: mi, m, ni, n, ib1, ib2
      integer :: ia, it, mdvdb, lmt1, lmt2, il1, il2, it1, it2
      integer :: lmta1, lmta2
      complex(kind=CMPLDP) :: zsum, wf1, wf2
      complex(kind=CMPLDP), allocatable :: wk_fsri(:)

      allocate( wk_fsri(nlmta) ); wk_fsri = 0.0d0

! -- start
      Do ik0=1, kv3, ndim_spinor
         if ( map_k(ik0) /= myrank_k ) cycle

         iktmp = ( ik0 -1 )/ndim_spinor +1

         Do is2=1, ndim_spinor
            ik2 = ik0 +is2 -1

            do n=1,num_bands
               ni = ib_inc(n)
               ib2 = neordr(ni,ik2)

               if ( map_e(ib2) == mype ) then
                  wk_fsri(:) = dcmplx( fsr_l( map_z(ib2),:,ik2 ), &
                       &               fsi_l( map_z(ib2),:,ik2 ) )
               endif
               call mpi_bcast( wk_fsri, 2*nlmta, mpi_double_precision, &
                    &          map_e(ib2), mpi_k_world(myrank_k), ierr )

               Do is1=1, ndim_spinor
                  ik1 = ik0 +is1 -1

                  do m=1,num_bands
                     mi = ib_inc(m)
                     ib1 = neordr(mi,ik1)

                     if ( map_e(ib1) /= myrank_e ) cycle

                     zsum = 0.d0

                     Do ia=1, natm
                        it = ityp(ia)
                        mdvdb = m_PP_include_vanderbilt_pot(it)
                        if ( mdvdb == SKIP ) cycle

                        Do lmt1=1, ilmt(it)
                           il1 = ltp(lmt1,it); it1 = taup(lmt1,it)
                           lmta1 = lmta( lmt1,ia )
                           wf1 = dcmplx( fsr_l( map_z(ib1), lmta1, ik1 ), &
                                &        fsi_l( map_z(ib1), lmta1, ik1 ) )

                           Do lmt2=1, ilmt(it)
                              il2 = ltp(lmt2,it); it2 = taup(lmt2,it)
                              lmta2 = lmta( lmt2,ia )
                              wf2 = wk_fsri(lmta2)

                              zsum = zsum +conjg(wf1) * wf2 *q(lmt1,lmt2,it)
                           End do
                        End Do
                     End Do

                     if ( is1 == 1 .and. is2 == 1 ) then
                        spn_mat(3,m,n,iktmp) = spn_mat(3,m,n,iktmp) +zsum
                     else if ( is1 == 1 .and. is2 == 2 ) then
                        spn_mat(1,m,n,iktmp) = spn_mat(1,m,n,iktmp) +zsum
                        spn_mat(2,m,n,iktmp) = spn_mat(2,m,n,iktmp) -zsum *zi
                     else if ( is1 == 2 .and. is2 == 1 ) then
                        spn_mat(1,m,n,iktmp) = spn_mat(1,m,n,iktmp) +zsum
                        spn_mat(2,m,n,iktmp) = spn_mat(2,m,n,iktmp) +zsum *zi
                     else if ( is1 == 2 .and. is2 == 2 ) then
                        spn_mat(3,m,n,iktmp) = spn_mat(3,m,n,iktmp) -zsum
                     endif

                  End do
               End Do
            End DO
         End Do
      End Do

      deallocate( wk_fsri )

    end subroutine contrib_hardpart

    subroutine gather_matrix
      complex(kind=CMPLDP), allocatable :: spn_mat_mpi(:,:,:,:) ! d(neg,neg,kv3/nspin,3)

      if (npes>1) then
         allocate(spn_mat_mpi(npauli,num_bands,num_bands,kv3/nspin) )
         spn_mat_mpi = spn_mat
         spn_mat = 0.d0
         call mpi_allreduce( spn_mat_mpi, spn_mat, &
              &              num_bands*num_bands*kv3/nspin*npauli*2, &
              &              mpi_double_precision, mpi_sum, MPI_CommGroup,ierr )
         deallocate(spn_mat_mpi)
      end if

      call mpi_barrier( mpi_comm_world, ierr )

    end subroutine gather_matrix

    subroutine print_mat
      integer :: ik, mi, ni, ispn, num, nsize
      complex(kind=CMPLDP), allocatable :: work(:,:)

      if ( mype /= 0 ) return

      open(nfwannier,file=trim(wan90_seedname)//".spn",form="formatted")

#ifdef WAN90_SPN_FORMATTED
      write(nfwannier,*) 'Generated by PHASE noncol'
      write(nfwannier,*) num_bands, kv3/nspin

      do ik=1, kv3/nspin
         Do mi=1, num_bands
            Do ni=1, mi
               Do ispn=1, 3
                  write(nfwannier,'(2f16.12)') spn_mat(ispn,mi,ni,ik)
               End Do
            End Do
         End do
      end do
#else
      write(nfwannier) 'Generated by PHASE noncol'
      write(nfwannier) num_bands, kv3/nspin

      nsize = num_bands *( num_bands +1 ) /2
      allocate( work(3, nsize) ); work = 0.0d0

      Do ik=1, kv3/nspin
         num = 0
         Do mi=1, num_bands
            Do ni=1, mi
               num = num +1
               work(:,num) = spn_mat(:,mi,ni,ik)
            End Do
         End Do
         write(nfwannier) ( (work(ispn,mi),ispn=1,3), mi=1, nsize )
      End do
      deallocate( work  )
#endif

      close(nfwannier)

    end subroutine print_mat

  end subroutine m_Wan90_gen_mat_spn
! ================== 2015/09/04

! ==== KT_add ==== 2015/09/14
  subroutine m_Wan90_set_igf_wan90
    integer :: id, i
    integer :: igf1, igf2, igf3
    integer :: nx, ny, nz
    integer :: nxmin, nymin, nzmin, nxmax, nymax, nzmax
!
    nxmin = minval( nncell(1,:,:) );   nxmax = maxval( nncell(1,:,:) );
    nymin = minval( nncell(2,:,:) );   nymax = maxval( nncell(2,:,:) );
    nzmin = minval( nncell(3,:,:) );   nzmax = maxval( nncell(3,:,:) );

    allocate( igf_wan90( kg, nxmin:nxmax, nymin:nymax, nzmin:nzmax ) )
    igf_wan90 = 0
!
    id = fft_box_size_WF(1,0)

    Do nx=nxmin, nxmax
       Do ny=nymin, nymax
          Do nz=nzmin, nzmax

             do i = 1, kg
                igf1 = ngabc(i,1) + 1 -nx
                igf2 = ngabc(i,2) + 1 -ny
                igf3 = ngabc(i,3) + 1 -nz
!
                if ( igf1 <= 0 ) igf1 = igf1 + fft_box_size_WF(1,1)
                if ( igf2 <= 0 ) igf2 = igf2 + fft_box_size_WF(2,1)
                if ( igf3 <= 0 ) igf3 = igf3 + fft_box_size_WF(3,1)

                igf_wan90(i,nx,ny,nz) = igf1 + (igf2-1)*id &
                     &                 + (igf3-1)*id*fft_box_size_WF(2,0)
             enddo
          End Do
       End Do
    End Do
  end subroutine m_Wan90_set_igf_wan90

  subroutine m_Wan90_map_WF_on_fftmesh(k1,k2,ik,psi_l,bfft,igf_in)
    integer, intent(in) :: k1,k2,ik
    integer, intent(in) :: igf_in(kg1)
    real(kind=DP), intent(in),dimension(kg1,1,k1:k2,kimg) :: psi_l
    real(kind=DP), intent(inout), dimension(nfft) :: bfft

    integer :: i,i1,ri, j, i2, ii

    bfft = 0.d0
    if(k_symmetry(ik) == GAMMA) then
       if(kimg == 1) then
          i1 = igf_in(1)
          bfft(i1) = psi_l(1,1,ik,1)
#ifdef NEC_TUNE_SMP
!CDIR NODEP
#endif
          do ii = 2, iba(ik)
!!$             i = nbase(ii,1)
             i = nbase(ii,ik)
             i1 = igf_in(i)
             bfft(i1) = psi_l(ii,1,ik,1)
             j = nbase_gamma(ii,2)
             i2 = igf_in(j)
             bfft(i2) =   psi_l(ii,1,ik,1)
          end do
       else if(kimg == 2) then
          i1 = 2*igf_in(1) - 1
          bfft(i1)   = psi_l(1,1,ik,1)
          bfft(i1+1) = psi_l(1,1,ik,2)
#ifdef NEC_TUNE_SMP
!CDIR NODEP
#endif
          do ii = 2, iba(ik)
!!$             i = nbase(ii,1)
             i = nbase(ii,ik)
             i1 = 2*igf_in(i)-1
             bfft(i1  ) = psi_l(ii,1,ik,1)
             bfft(i1+1) = psi_l(ii,1,ik,2)
             j = nbase_gamma(ii,2)
             i2 = 2*igf_in(j)-1
             bfft(i2  ) = psi_l(ii,1,ik,1)
             bfft(i2+1) = -psi_l(ii,1,ik,2)
          end do
       end if
    else
#ifdef NEC_TUNE_SMP
!CDIR NOLOOPCHG
#endif
       do ri = 1, kimg
#ifdef NEC_TUNE_SMP
!CDIR NODEP
#endif
          do i = 1, iba(ik)
             i1 = kimg*igf_in(nbase(i,ik)) + (ri - kimg)
             bfft(i1) = psi_l(i,1,ik,ri)   ! MPI
          end do
       end do
    end if
  end subroutine m_Wan90_map_WF_on_fftmesh
! ===== 2015/09/14

end module m_Wannier90
