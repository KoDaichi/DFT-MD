!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 606 $)
!
!  MODULE: m_CS_SpaceGroup
!
!  AUTHORS: T. Yamamoto,    October/16/2005
!
!  Contact address :  Phase System Consortium
!                     E-mail: phase_system@nims.go.jp URL https://azuma.nims.go.jp
!
!
!
!=======================================================================
!
!
module m_CS_SpaceGroup
!      (m_CS_SG)
! $Id: m_CS_SpaceGroup.f90 606 2020-04-15 06:45:49Z ktagami $
  use m_Const_Parameters,   only : DP, PAI, PAI2, BUCS, CARTS, ANTIFERRO &
       &                         , oh=>oh_symbol, d6h=>d6h_symbol, OFF, CONTINUATION &
       &                         , FIXED_CHARGE_CONTINUATION, COORDINATE_CONTINUATION
  use m_Timing,             only : tstatc0_begin, tstatc0_end, m_Timing_initialize_tmp_timer &
       &                         , m_Timing_get_elapsed_time_tmp, m_Timing_finalize_tmp_timer
  use m_Control_Parameters, only : ipri, af, m_CtrlP_set_af &
       &                         , ipriinputfile, ipri_spg, printable &
       &                         , icond
  use m_Crystal_Structure,  only : altv,rltv,nopr,op,tau,imag &
       &                         , ngen,igen,jgen,iaf,jaf &
       &                         , alloc_igen_jgen,dealloc_igen_jgen &
       &                         , b2pmat,ig01,il,pg_symbol_system &
       &                         , misalignment, p2bmat &
       &                         , m_CS_put_Bravais_lattice &
       &                         , sw_non_symmorphic, cpumax_symsearch &
       &                         , sw_read_symmetry_from_file
  use m_Ionic_System,       only : natm,natm2,ntyp,ityp,pos,iatomn,iwei &
       &                         , m_IS_put_lattice_system &
       &                         , m_IS_put_latconst_len_angle
!!$       &                         , latconst_len,latconst_angle

! ============================ added by K. Tagami ================== 11.0
  use m_Const_Parameters,   only : NONCOLLINEAR, DELTA10, ON
! ================================================================== 11.0

! ============================ added by K. Tagami ================== 11.0
  use m_Control_Parameters,   only : noncol
  use m_Ionic_System,       only : magmom_local_now
  use m_Crystal_Structure,  only : sw_use_magnetic_symmetry
! ================================================================== 11.0

  use m_Parallelization,    only : npes,mype,MPI_CommGroup
  use m_Files,              only : m_Files_file_exists, m_Files_open_nfcntn &
       &                         , m_Files_close_nfcntn, nfcntn, F_CNTN
  use mpi
  implicit none

  character(len=9), public :: lattice_system

! ========= Added by K. Tagami ========================
   real(kind=DP), dimension(3,3), public :: kt_b2pmat
   real(kind=DP), dimension(3,3), public :: kt_p2bmat
! =========

   real(kind=DP),dimension(3) :: Bravais_lattice_length, Bravais_lattice_angle

   logical, private :: cpumax_reached = .false.
   integer, private :: mnprim, mnsym
   real(kind=DP), private, allocatable, dimension(:,:) :: mtprim, mtau
   logical, private, allocatable, dimension(:) :: mfsym

   integer, private :: mnprim_af, mnsym_af
   real(kind=DP), private, allocatable, dimension(:,:) :: mtprim_af, mtau_af
   logical, private, allocatable, dimension(:) :: mfsym_af

   character(len("symmetry_info")), private, parameter ::    tag_symmetry_info    = "symmetry_info"
   character(len("symmetry_info_af")), private, parameter :: tag_symmetry_info_af = "symmetry_info_af"

!   include 'mpif.h'
contains
    subroutine mat3inv(mat,matinv)
    implicit none

    real(kind=DP), intent(in)  :: mat(3,3)
    real(kind=DP), intent(out) :: matinv(3,3)

    ! local variables
    real(kind=DP) :: det
    integer :: i,j,i1,i2,j1,j2

    det=0.d0
    do i=1,3
      do j=1,3
        i1=mod(i,3)+1
        i2=mod(i+1,3)+1
        j1=mod(j,3)+1
        j2=mod(j+1,3)+1
        matinv(j,i)=mat(i1,j1)*mat(i2,j2)-mat(i1,j2)*mat(i2,j1)
      end do
      det=det+mat(i,1)*matinv(1,i)
    end do
    do i=1,3
      do j=1,3
        matinv(i,j)=matinv(i,j)/det
      end do
    end do

    return
    end subroutine mat3inv

  subroutine set_lattice_system(system)
    character(len=9), intent(in) :: system
    lattice_system = system
    call m_IS_put_lattice_system(lattice_system)
    pg_symbol_system = system
  end subroutine set_lattice_system

  subroutine set_ig01(iopr,nopr)
    integer, intent(in) :: iopr
    integer, intent(in) :: nopr(49)

    ig01(1:iopr) = nopr(1:iopr)
    if(imag==ANTIFERRO) then
       iaf = nopr(iopr+af)
    end if
  end subroutine set_ig01

  subroutine m_CS_SG_auto_gnrt_sym_op(paramset,nfout)
    logical, intent(in) :: paramset
    integer, intent(in) :: nfout
    integer ::             i, j, iopr, ip
    integer ::             af_t
    integer ::             id_sname = -1
    call tstatc0_begin('m_CS_SG_gnrt_symmetry_operations ',id_sname)
! ----- iopr, op, tau

    call autogen_space_group(nfout,ipri_spg,iopr,af_t,paramset) !! defined in this file

    if(printable) write(nfout,'(" !! af_t = ",i10)') af_t
    if(af_t > 1) call phase_error_with_msg(nfout, ' ! illegal af_t value',__LINE__,__FILE__)
    call m_CtrlP_set_af(af_t)
    if(printable) write(nfout,'(" !! af_t (after m_CtrlP_set_af) = " ,i6)') af_t

    if(ipri >= 0 .and. printable) then
       write(nfout,*) ' << rltv, altv >>'
       write(nfout,'(3f8.4,11x,3f8.4)') &
            & ((rltv(i,j),j=1,3),(altv(i,j),j=1,3),i=1,3)
    endif
    if(paramset) then
       nopr = iopr
       goto 1001
    end if

    call tau_in_CARTS()

1001 continue

    call tstatc0_end(id_sname)
  contains
    subroutine tau_in_CARTS
      integer i, no
      do i = 1, 3
         do no = 1, iopr + af
            tau(i,no,CARTS) = altv(i,1)*tau(1,no,BUCS) &
                 &     +      altv(i,2)*tau(2,no,BUCS) &
                 &     +      altv(i,3)*tau(3,no,BUCS)
         enddo
      enddo
    end subroutine tau_in_CARTS
  end subroutine m_CS_SG_auto_gnrt_sym_op

  subroutine autogen_space_group(nfout,ipri,iopr,af,paramset)
    implicit none
    integer, intent(in) :: nfout,ipri
    integer, intent(out) :: iopr,af
    logical, intent(in) :: paramset

    real(kind=DP) :: r(3,3,49),t(3,49)
    real(kind=DP) :: r_af(3,3,49),t_af(3,49)
    real(kind=DP) :: pos2(natm2,3)
    real(kind=DP) :: tprim(3,natm2),tprim_af(3,natm2)
    integer :: nopr(49),nopr_af(49)
    integer :: i,ia
    integer :: ipri_t
    integer :: ityp2(natm2)
    integer :: ityp_af(natm2)
    integer :: iopr_af
    integer :: nprim,nprim_af
    integer, parameter :: DEBUGPRINTLEVEL = 2
    character(len=9) :: system
! =================================== Added by K. Tagami ==========
    integer j,k,l,m,n
    real(kind=DP) :: tmp
!!!!    real(kind=DP) :: bmat(3,3), binv(3,3), p2bmat(3,3)
    real(kind=DP) :: bmat(3,3), binv(3,3)

! ============================== added by K. Tagami ================ 11.0
    real(kind=DP) :: mag_loc(natm2,3)
! ==================================================================== 11.0

! ============================== Added by K. Tagami ==== patch 0.9 =========
       Do j=1, 3
          Do k=1, 3
             kt_b2pmat(j,k) = b2pmat(k,j)
             kt_p2bmat(j,k) = p2bmat(k,j)
          End Do
       End do
!       Call mat3inv( kt_b2pmat, kt_p2bmat )
! =========================================================

    if(ipri>2) then
       ipri_t = 2
    else if(paramset) then
       ipri_t = 0
    else
       ipri_t = ipri
    endif

    pos2(1:natm,1:3) = pos(1:natm,1:3)
    ityp2(1:natm) = ityp(1:natm)
    i=natm
    do ia=1,natm
       if(iwei(ia)/=1) then
          i = i+1
          pos2(i,1:3) = -pos(ia,1:3)
          ityp2(i) = ityp(ia)
       end if
    end do

    iopr = 48
    call space_group(nfout,ipri_t,altv(1,1),altv(1,2),altv(1,3), &
         &           ntyp,natm2,ityp2,pos2,system,iopr,nopr,r,t, &
         &           nprim,tprim,paramset)

!!!!!!    write(*,*) 'iopr a =', iopr

    if(ipri>=DEBUGPRINTLEVEL) then
       write(nfout,'(" number of symmetry = ",i4)') iopr
       write(nfout,'(" system = ",a)') system
       call wd_operations(nfout,r,t,nopr,iopr)
    end if

    if(imag==ANTIFERRO) then
      do ia=1,natm2
        ityp_af(ia) = nint(iatomn(ityp2(ia)))
      end do
      iopr_af = 48
      call space_group(nfout,0,altv(1,1),altv(1,2),altv(1,3), &
                     & ntyp,natm2,ityp_af,pos2,system,iopr_af,nopr_af,r_af,t_af, &
                     & nprim_af,tprim_af,paramset,af=.true.)
      if(ipri>=DEBUGPRINTLEVEL) then
         write(nfout,'(" number of symmetry when imag=ANTIFERRO = ",i8)') iopr_af
         write(nfout,'(" system = ",a)') system
         call wd_operations(nfout,r_af,t_af,nopr_af,iopr_af)
         write(nfout,'(" nprim_af, nprim = ",2i8)') nprim_af, nprim
      end if
      if(iopr_af.ne.2*iopr.and.nprim_af.le.nprim) then
         write(nfout,*) 'iopr_af=',iopr_af
         write(nfout,*) 'iopr   =',iopr
         write(nfout,*) 'nprim_af=',nprim_af
         write(nfout,*) 'nprim=',nprim
         call phase_error_with_msg(nfout, 'iopr_af .ne. 2*iopr and nprim_af .le. nprim',__LINE__,__FILE__)
      end if
!=========================== Modified by K. Tagami ======= 1.0 ===
!      call set_af_operator(nfout,ipri_t,system,iopr_af,nopr_af,r_af,t_af,nprim_af,tprim_af,iopr,nopr,r,t,nprim,tprim)
      call set_af_operator_kt(nfout,ipri_t,system,iopr_af,nopr_af,r_af,t_af,nprim_af,tprim_af,iopr,nopr,r,t,nprim,tprim)
! =================================================================
      af = 1
    else
      af = 0
    end if

    call print_prim_trans(nfout,ipri_t,nprim,tprim)

! ========================== added by K. Tagami =============================== 11.0
!    if ( noncol .and. sw_use_magnetic_symmetry == ON ) then
!       i = natm
!       mag_loc(1:natm,:) = magmom_local_now(1:natm,:)
!       do ia=1,natm
!          if (iwei(ia)/=1) then
!             i = i+1
!             mag_loc(i,1:3) = magmom_local_now(ia,1:3)
!          end if
!       end do
!       call set_magnetic_symm( natm2, pos2, iopr, nopr, r, t, mag_loc )
!    endif
! ============================================================================ 11.0

    if(.not.paramset) then
! =================================== Added by K. Tagami ==========
!!!!!!!!!!!!!!      op(:,:,:) = r(:,:,1:iopr+af)
!
!       Call mat3inv( b2pmat, p2bmat )
!       bmat = matmul( altv, p2bmat )
!       Call mat3inv( bmat, binv )

       bmat = matmul( altv, kt_p2bmat )
       Call mat3inv( bmat, binv )
! =============================================================
!
       do i=1, iopr+af
          do n=1,3
            do m=1,3
               tmp = 0.d0
                do k=1,3
                  do l=1,3
                     tmp = tmp  + bmat(m,k) * r(k,l,i) *binv(l,n)
                  end do
               end do
               op(m,n,i) = tmp
            end do
         end do
      end do
! ========================================================================
      do i=1,iopr+af
        tau(1:3,i,BUCS) = t(1:3,i)
      end do
    end if

    call set_tspace_generators(nfout,ipri_t,iopr,nopr,r,t,natm2,nprim,tprim)
    if(ipri_t>0) call write_tspace_generators(nfout,system,nprim)

    call set_lattice_system(system)

    call set_ig01(iopr,nopr)

  end subroutine autogen_space_group

! ================================ added by K. Tagami ================= 11.0
  subroutine set_magnetic_symm( natom, apos, nsym,  nopr, rot, tau, mag_loc )
    implicit none

    integer, intent(in) :: natom
    integer, intent(inout) :: nsym
    integer, intent(out) :: nopr(48)
    real(kind=DP), intent(inout) :: rot(3,3,48), tau(3,48)
    real(kind=DP), intent(in) :: apos(natom,3)
    real(kind=DP), intent(in) :: mag_loc(natom,3)

    real(kind=DP) :: ktrot_pr(3,3,48), ktrot_carts(3,3,48)
    real(kind=DP) :: bmat(3,3), binv(3,3)
    logical :: mag_sym_flag(48)

    integer :: i, ia, ja, ja_found
    integer :: isym, count
    real(kind=DP) :: coord_tmp(3), dcoord(3), dist
    real(kind=DP) :: mag_tmp(3), dmag_p(3), dmag_m(3)
    real(kind=DP) :: c1, c2
!
    real(kind=DP), parameter :: criterion = 1.0D-5
!
    real(kind=DP), allocatable :: apos_t(:,:)
!
    integer :: it
! -- init --
    call set_drot( kt_b2pmat, kt_p2bmat, nsym, rot, ktrot_pr )

    Call matrix_product( altv, kt_p2bmat, bmat )
    Call mat3inv( bmat, binv )
    call set_drot( binv, bmat, nsym, rot, ktrot_carts )

    mag_sym_flag = .false.

    allocate(apos_t(natom,3))
    apos_t(1:natm,1:3) = apos(1:natom,1:3)
    do ia = 1, natom
       do i = 1, 3
          apos_t(ia,i) = apos_t(ia,i) - floor(apos_t(ia,i))
       end do
    end do

! ---
!    Do ia=1, natom
!       write(*,*) 'ia ', ia
!       write(*,*) apos_t(ia,1), apos_t(ia,2), apos_t(ia,3)
!    End do
! ---
!    Do isym=1, nsym
!       write(*,*) 'isym ', isym
!       write(*,*) ktrot_carts(1,1,isym), ktrot_carts(1,2,isym), ktrot_carts(1,3,isym)
!       write(*,*) ktrot_carts(2,1,isym), ktrot_carts(2,2,isym), ktrot_carts(2,3,isym)
!       write(*,*) ktrot_carts(3,1,isym), ktrot_carts(3,2,isym), ktrot_carts(3,3,isym)
!    End do

! -- begin --
    Do isym=1, nsym

       Do it=1, ntyp
          Do ia=1, natom
             if ( ityp(ia) /=it ) cycle

             coord_tmp(1:3) = matmul( ktrot_pr(:,:,isym), apos_t(ia,1:3) ) &
                  &         + tau(1:3,isym)
             coord_tmp(1:3) = coord_tmp(1:3) - floor(coord_tmp(1:3))

             ja_found = 0
             Do ja=1, natom
                dcoord(:) = apos_t(ja,1:3) - coord_tmp(:)

                dist = dcoord(1)**2 + dcoord(2)**2 + dcoord(3)**2
                dist = sqrt(dist)
                dist = dist - floor( dist +DELTA10 )

!!             write(*,*) 'isym, ia ja, dist =', isym, ia, ja, dist
                if ( dist < criterion ) then
                   ja_found = ja
!                goto 100
                   exit
                endif
             End do

100          continue
!
             if ( ja_found ==0 ) then
                write(*,*) 'kt: Not found symmetry', isym, ia
                !stop
                call phase_error_with_msg(6,'ikt: Not found symmetry',__LINE__,__FILE__)
             endif

             mag_tmp(1:3) = matmul( ktrot_carts(:,:,isym), mag_loc(ia,1:3) )
          !
             dmag_m(1:3) = mag_tmp(1:3) - mag_loc(ja_found,1:3)
             dmag_p(1:3) = mag_tmp(1:3) + mag_loc(ja_found,1:3)
          !
             c1 = abs( dmag_m(1) ) +abs( dmag_m(2) ) +abs( dmag_m(3) )
             c2 = abs( dmag_p(1) ) +abs( dmag_p(2) ) +abs( dmag_p(3) )
          !
!!!          if ( c1 < criterion .or. c2 < criterion ) then
             if ( c1 < criterion ) then
                mag_sym_flag(isym) = .true.
                exit
             endif
!
             if ( mag_sym_flag(isym) ) then
                if ( c1 > criterion ) then
                   mag_sym_flag(isym) = .false.
                endif
             endif
          End do
       End do
    End do
! -----
    count = 0
    Do isym=1, nsym
       if ( mag_sym_flag(isym) ) then
          count = count + 1
          nopr(count) = nopr(isym)
          rot(:,:,count) = rot(:,:,isym)
          tau(:,count) = tau(:,isym)
       endif
    End do

    nsym = count

    deallocate( apos_t )

  contains

    subroutine set_drot( a,ainv,nsym,rot,drot )
      real(kind=DP), intent(in) :: a(3,3),ainv(3,3)
      integer, intent(in) :: nsym
      real(kind=DP), intent(in) :: rot(3,3,nsym)
      real(kind=DP), intent(out) :: drot(3,3,nsym)

      integer :: n
      do n=1,nsym
         drot(1:3,1:3,n) = matmul(ainv,matmul(rot(:,:,n),a))
      end do
    end subroutine set_drot

    subroutine matrix_product(a,b,c)
      real(kind=DP), intent(in) :: a(3,3),b(3,3)
      real(kind=DP), intent(out) :: c(3,3)

      ! local variables
      integer :: i,j

      do j=1,3
         do i=1,3
            c(i,j) = sum(a(i,1:3)*b(1:3,j))
         end do
      end do

    end subroutine matrix_product

  end subroutine set_magnetic_symm
!=========================================================================== 11.0

  subroutine wd_operations(nfout,r,t,nopr,iopr)
    integer, intent(in) :: nfout, iopr, nopr(iopr)
    real(kind=DP), intent(in) :: r(3,3,iopr),t(3,iopr)
    integer :: n,j, i

    write(nfout,'(" iopr = ",i8," at <<autogen_space_group>> in m_CS_spaceGroup")') iopr
    do n = 1, iopr
       write(nfout,'(" #op = ", i8)') nopr(n)
       do i = 1, 3
          write(nfout,'(4f8.4)') (r(i,j,n),j=1,3),t(i,n)
       end do
    end do
  end subroutine wd_operations

  subroutine space_group(nfout,ipri,a1,a2,a3,ntype,natom &
                      & ,itype,apos,system,nsym,nopr,rot,tau &
                      & ,nprim,tprim,paramset,af)
    implicit none

    integer, intent(in) :: nfout,ipri
    real(kind=DP), intent(in) :: a1(3),a2(3),a3(3)
    integer, intent(in) :: ntype,natom,itype(natom)
    real(kind=DP), intent(in) :: apos(natom,3)
    character(len=9), intent(out) :: system
    integer, intent(out) :: nsym
    integer, intent(out) :: nopr(48)
    real(kind=DP), intent(out) :: rot(3,3,48),tau(3,48)
    integer, intent(out) :: nprim
    real(kind=DP), intent(out) :: tprim(3,natom)
    logical, intent(in) :: paramset
    logical, intent(in), optional :: af

    logical :: flat(48)
    integer :: nlat
    character(len=3) :: pg_name
    character(len=5) :: pg_name_i
    character(len=2) :: bl_name
    character(len=10) :: sg_name_i
    logical :: from_af

    integer :: i,j
    flat = .true.
    from_af = .false.
    if(present(af)) from_af = af

    call determine_crystal_system(a1,a2,a3,system)
    if(system=='cubic')then
       kt_b2pmat=0.d0
       kt_p2bmat=0.d0
       do i=1,3
          kt_b2pmat(i,i) = 1.0d0
          kt_p2bmat(i,i) = 1.0d0
       enddo
    endif
10  continue ! changing system for a monoclinic crystal with a hexagonal cell
    call generate_point_group(system,nsym,rot)
    ! debug
    if(ipri>1) call print_rot(nsym,rot,flat)
    ! end debug
    call determine_lattice_point_group(system,a1,a2,a3,nsym,nlat,nopr,rot,flat)
    call get_point_group_name(nfout,system,nlat,nopr,pg_name,pg_name_i)
    call get_bravais_lattice(nfout,0,a1,a2,a3,system,nopr,pg_name_i,bl_name) ! ipri = 0
    ! debug
    if(ipri>1) then
       write(nfout,'(9x,"Lattice point group: ",a,1x,a)') pg_name,pg_name_i
       write(nfout,'(13x,"Bravais lattice: ",a,1x,a)') bl_name
       write(nfout,'("Elements of the lattice point group:")')
       call print_rot(nsym,rot,flat)
    end if
    ! end debug
    call determine_space_group(a1,a2,a3,ntype,natom,itype,apos,nsym,nopr,rot,tau,flat,nprim,tprim,paramset,from_af)
    if(system_is_wrong(system,nsym,nopr)) go to 10
    !!$call get_point_group_name(nfout,system,nsym,nopr,pg_name,pg_name_i)
    !!$call get_bravais_lattice(nfout,ipri,a1,a2,a3,system,nopr,pg_name_i,bl_name)
    !!$call get_space_group_name(nfout,system,nsym,nopr,tau,.false.,bl_name,pg_name_i,sg_name_i)
    if(ipri>1) then
   !!$    write(nfout,'(13x,"Bravais lattice: ",a,1x,a)') bl_name
   !!$    write(nfout,'("Crystallographic point group: ",a,1x,a)') pg_name,pg_name_i
   !!$    write(nfout,'(17x,"Space group: ",a)') sg_name_i
   !!$    write(nfout,'("Elements of the space group:")')
       call print_space_group(nfout,system,nsym,nopr,rot,tau)
    end if

  contains

    logical function system_is_wrong(system,nsym,nopr)
      character(len=9), intent(inout) :: system
      integer, intent(in) :: nsym
      integer, intent(in) :: nopr(48)

      logical :: fexist(48)
      integer :: i

      fexist=.false.
      do i=1,nsym
         fexist(nopr(i)) = .true.
      end do

      system_is_wrong = .false.
      if(system == 'hexagonal'.and.il>0) then
         if(.not.fexist(3)) then
            ! lattice system is monoclinic if C3+ doesn't exist.
            system_is_wrong = .true.
            system = 'cubic'
            if(ipri>0) then
               write(nfout,*) 'The lattice system was found monoclinic or orthorhombic, and'
               write(nfout,*) 'system was changed from hexagonal to cubic.'
            end if
         end if
      end if
    end function system_is_wrong

    subroutine determine_crystal_system(a1,a2,a3,system)
    implicit none

    real(kind=DP), intent(in) :: a1(3),a2(3),a3(3)
    character(len=9), intent(out) :: system

    ! local variables
    real(kind=DP) :: length(3),angle(3)
    real(kind=DP) :: eps = 1.d-6
    real(kind=DP) :: pi2,pi3,pi4
    real(kind=DP) :: axis(3),axis_angle
    integer :: i

    pi2 = atan(1.d0)*2.d0
    pi3 = pi2*4.d0/3.d0
    pi4 = pi3*0.5d0
    do i=1,3
       axis(i) = a1(i)+a2(i)+a3(i)
    end do

    call vectors_length_angle(a1,a2,a3,length,angle)
    call get_axis_angle(axis,axis_angle)

    if(ipri>0) then
      write(nfout,'("a1=",f10.5," a2=",f10.5," a3=",f10.5,1x," a.u.")') length(1:3)
      write(nfout,'("a =",f10.5," b =",f10.5," g =",f10.5,1x," deg.")') angle(1:3)/(2.d0*pi2)*180
      write(nfout,'("axis angle",f10.5,1x," deg.")') axis_angle/(2.d0*pi2)*180
    end if

    if(il==0) then
       system = 'hexagonal'
       if(ipri>0) write(nfout,'("System is hexagonal (manual)")')
       return
    else if(il==-1) then
       system = 'hexagonal'
       if(ipri>0) write(nfout,'("System is rhombohedral (manual)")')
       return
    end if

    if(abs(length(1)-length(2)) < eps .and. &
     & abs(angle(1)-angle(2)) < eps) then

      if(abs(length(2)-length(3)) > eps .and. &
       & abs(angle(2)-pi2) < eps .and. &
       & abs(angle(3)-pi3) < eps) then
        system = 'hexagonal'
        if(ipri>0) write(nfout,'("System is hexagonal")')
        return
      else if(abs(length(2)-length(3)) < eps .and. &
            & abs(angle(2)-angle(3)) < eps .and. &
            & abs(angle(3)-pi2) > eps .and. &
            & abs(angle(3)-pi4) > eps .and. &
            & abs(axis_angle) < eps ) then
        system = 'hexagonal'
        if(ipri>0) write(nfout,'("System is rhombohedral")')
        return
      end if

    end if

    system = 'cubic'
    if(ipri>0) write(nfout,'("System is not hexagonal nor rhombohedral")')

    end subroutine determine_crystal_system


    subroutine get_axis_angle(axis,axis_angle)
    implicit none
    real(kind=DP), intent(in) :: axis(3)
    real(kind=DP), intent(out) :: axis_angle

    axis_angle = sqrt(sum(axis(1:3)**2))
    axis_angle = axis(3)/axis_angle
    axis_angle = acos(axis_angle)

    end subroutine get_axis_angle

    subroutine determine_lattice_point_group(system,a1,a2,a3,nsym,nlat,nopr,rot,flat)
      implicit none
      character(len=9), intent(in) :: system
      integer, intent(in) :: nsym
      integer, intent(out) :: nlat
      integer, intent(out) :: nopr(48)
      real(kind=DP), intent(in) :: a1(3),a2(3),a3(3)
      real(kind=DP), intent(in) :: rot(3,3,48)
      logical, intent(out) :: flat(48)

      ! local variables
      real(kind=DP) :: grid(3,16)
      real(kind=DP) :: x(3)
      integer :: i,j,n,isum
      real(kind=DP), parameter :: eps = 1.d-6
!============================== Added by K. Tagami ==============
      real(kind=DP) :: ktrot(3,3,48)
      real(kind=DP) :: a(3,3),ainv(3,3)
!!!!!!!!!      real(kind=DP) :: bmat(3,3), binv(3,3), p2bmat(3,3)
      real(kind=DP) :: bmat(3,3), binv(3,3)
      integer k, l
!
      a(1:3,1) = a1(1:3)
      a(1:3,2) = a2(1:3)
      a(1:3,3) = a3(1:3)

!      Call mat3inv( b2pmat, p2bmat )
!      Call matrix_product( altv, p2bmat, bmat )
!      Call mat3inv( bmat, binv )
!      call set_drot( binv, bmat, nsym,rot,ktrot )

      Call matrix_product( altv, kt_p2bmat, bmat )
      Call mat3inv( bmat, binv )
      call set_drot( binv, bmat, nsym,rot,ktrot )
!=======================================================================

!      write(900,*) 'DEBUG a =='
!      write(900,*) ktrot(1,1,4), ktrot(1,2,4), ktrot(1,3,4)
!      write(900,*) ktrot(2,1,4), ktrot(2,2,4), ktrot(2,3,4)
!      write(900,*) ktrot(3,1,4), ktrot(3,2,4), ktrot(3,3,4)

! --

      flat = .true.

      grid(1:3,1) = a1(1:3)
      grid(1:3,2) = a2(1:3)
      grid(1:3,3) = a3(1:3)
      do i=1,3
         grid(1:3,i+3) = -grid(1:3,i)
      end do
      grid(1:3,7) = a1(1:3)-a2(1:3)
      grid(1:3,8) = a2(1:3)-a3(1:3)
      grid(1:3,9) = a3(1:3)-a1(1:3)
      do i=7,9
         grid(1:3,i+3) = -grid(1:3,i)
      end do
      grid(1:3,13) = a1(1:3)+a2(1:3)+a3(1:3)
      grid(1:3,14) = -grid(1:3,13)
      grid(1:3,15) = a1(1:3)+a2(1:3)
      grid(1:3,16) = -grid(1:3,15)

      do n=1,nsym
         isum = 0
         LOOP: do i=1,3
! ================================ added by K. Tagami ================
!!            x(1:3) = rot(1:3,1,n)*grid(1,i) &
!!                 & + rot(1:3,2,n)*grid(2,i) &
!!                 & + rot(1:3,3,n)*grid(3,i)
            if ( system == "cubic" ) then
!               x(1:3) = nint(ktrot(1:3,1,n)) *grid(1,i) &
!                    & + nint(ktrot(1:3,2,n)) *grid(2,i) &
!                    & + nint(ktrot(1:3,3,n)) *grid(3,i)
               x(1:3) = ktrot(1:3,1,n) *grid(1,i) &
                    & + ktrot(1:3,2,n) *grid(2,i) &
                    & + ktrot(1:3,3,n) *grid(3,i)
            else               ! hexagonal
               x(1:3) = ktrot(1:3,1,n) *grid(1,i) &
                    & + ktrot(1:3,2,n) *grid(2,i) &
                    & + ktrot(1:3,3,n) *grid(3,i)
            endif
! ====================================================================
            do j=1,16
               if(abs(x(1)-grid(1,j)) < eps .and. &
               &  abs(x(2)-grid(2,j)) < eps .and. &
               &  abs(x(3)-grid(3,j)) < eps ) then
                  isum = isum + i
                  cycle LOOP
! debug
!!               else
!!                 write(*,'("n=",i2,"i=",i2,"j=",i2," X=",3f10.5," G=",3f10.5)') n,i,j,x,grid(1:3,j)
! end debug
               end if
            end do
         end do LOOP
         if(isum/=6) flat(n) = .false.
      end do

      nlat=0
      do n=1,nsym
         if(flat(n)) then
            nlat = nlat+1
            nopr(nlat) = n
         end if
      end do

    end subroutine determine_lattice_point_group

    subroutine generate_point_group(system,nsym,rot)
    implicit none

    character(len=9), intent(in) :: system
    integer, intent(out) :: nsym
    real(kind=DP), intent(out) :: rot(3,3,48)

    if(system .eq. 'cubic') then
      call cubic_point_group(nsym,rot)
    else if(system .eq. 'hexagonal') then
! ================================= added by K. Tagami ===
!      call hexagonal_point_group(nsym,rot)
      call hexagonal_point_group_bravais(nsym,rot)
! =======================================================
    else
      if(ipri>0) write(nfout,*) 'generate_point_group: system must be cubic or hexagonal.'
    end if

    end subroutine generate_point_group

    subroutine cubic_point_group(nsym,rot)
    implicit none

    integer, intent(out) :: nsym
    real(kind=DP), intent(out) :: rot(3,3,48)
    real(kind=DP) :: tmprot(3,3,48),ainv(3,3)

    ! local variables
    real(kind=DP) :: one=1.d0, zero=0.d0
    integer :: i,j

    nsym=48

    rot(1:3,1:3,1:nsym)=zero

    ! C31+ = O5
    rot(1,3,5) = one
    rot(2,1,5) = one
    rot(3,2,5) = one

    ! C4x+ = O19
    rot(1,1,19) =  one
    rot(2,3,19) = -one
    rot(3,2,19) =  one

    do i=1,3
      rot(i,i,1) = one  ! E = O1
    end do
    do j=1,3
      do i=1,3
        rot(i,j,22) = rot(j,i,19)  ! C4x- = O22 = O19^t
        rot(i,j,9)  = rot(j,i,5)  ! C31- = O9 = O5^t
      end do
    end do

    call matrix_product(rot(1,1,19),rot(1,1,19),rot(1,1,2)) ! C2x = O2 = O19*O19
    call matrix_product(rot(1,1,5),rot(1,1,19),rot(1,1,13)) ! C2a = O13 = O5*O19
    call matrix_product(rot(1,1,19),rot(1,1,5),rot(1,1,15)) ! C2c = O15 = O19*O5

    call matrix_product(rot(1,1,2),rot(1,1,9),rot(1,1,10))  ! C32- = O10 = O2*O9
    call matrix_product(rot(1,1,13),rot(1,1,15),rot(1,1,11))! C33- = O11 = O13*O15
    call matrix_product(rot(1,1,9),rot(1,1,2),rot(1,1,12))  ! C34- = O12 = O9*O2
    call matrix_product(rot(1,1,2),rot(1,1,5),rot(1,1,8))   ! C34+ = O8 = O2*O5
    call matrix_product(rot(1,1,5),rot(1,1,2),rot(1,1,6))   ! C32+ = O6 = O5*O2
    call matrix_product(rot(1,1,15),rot(1,1,13),rot(1,1,7)) ! C31+ = O7 = O15*O13
    call matrix_product(rot(1,1,13),rot(1,1,2),rot(1,1,21)) ! C4z+ = O21 = O13*O2
    call matrix_product(rot(1,1,2),rot(1,1,13),rot(1,1,24)) ! C4z- = O24 = O2*O13
    call matrix_product(rot(1,1,15),rot(1,1,2),rot(1,1,23)) ! C4y- = O23 = O15*O2
    call matrix_product(rot(1,1,2),rot(1,1,15),rot(1,1,20)) ! C4y+ = O20 = O2*O15

    call matrix_product(rot(1,1,9),rot(1,1,7),rot(1,1,3))  ! C2y = O3 = O9*O7
    call matrix_product(rot(1,1,9),rot(1,1,8),rot(1,1,4))  ! C2z = O4 = O9*O8
    call matrix_product(rot(1,1,15),rot(1,1,6),rot(1,1,14)) ! C2b = O14 = O15*O6
    call matrix_product(rot(1,1,13),rot(1,1,7),rot(1,1,18)) ! C2f = O18 = O13*O7
    call matrix_product(rot(1,1,13),rot(1,1,8),rot(1,1,16)) ! C2d = O16 = O13*O8
    call matrix_product(rot(1,1,7),rot(1,1,24),rot(1,1,17)) ! C2e = O17 = O7*O24

    do i=1,24
      rot(1:3,1:3,i+24) = -rot(1:3,1:3,i) ! add inversion
    end do

!!  to fractional
    call mat3inv(altv,ainv)
    call set_drot(altv,ainv,nsym,rot,tmprot)
    rot = tmprot
    end subroutine cubic_point_group

    subroutine hexagonal_point_group(nsym,rot)
    implicit none

    integer, intent(out) :: nsym
    real(kind=DP), intent(out) :: rot(3,3,48)

    ! local variables
    real(kind=DP) :: one=1.d0, zero=0.d0, half=0.5d0, r32
    integer :: i,j

    r32=sqrt(3.d0)*0.5d0

    nsym=24

    rot(1:3,1:3,1:nsym)=zero

    ! C6+ = O2
    rot(1,1,2) =  half
    rot(1,2,2) = -r32
    rot(2,1,2) =  r32
    rot(2,2,2) =  half
    rot(3,3,2) =  one

    ! C222 = O11
    rot(1,1,11) = -half
    rot(1,2,11) = -r32
    rot(2,1,11) = -r32
    rot(2,2,11) =  half
    rot(3,3,11) = -one


    do i=1,3
      rot(i,i,1) = one  ! E = O1
    end do

    do j=1,3
      do i=1,3
        rot(i,j,6) = rot(j,i,2)  ! C6- = O6 = O2^t
      end do
    end do
    call matrix_product(rot(1,1,2),rot(1,1,2),rot(1,1,3))   ! C3+ = O3 = O2*O2
    call matrix_product(rot(1,1,2),rot(1,1,11),rot(1,1,9))  ! C231 = O9 = O2*O11
    call matrix_product(rot(1,1,11),rot(1,1,2),rot(1,1,7))  ! C211 = O7 = O11*O2

    do j=1,3
      do i=1,3
        rot(i,j,5) = rot(j,i,3)  ! C3- = O5 = O3^t
      end do
    end do
    call matrix_product(rot(1,1,2),rot(1,1,3),rot(1,1,4))  ! C2 = O4 = O2*O3
    call matrix_product(rot(1,1,2),rot(1,1,9),rot(1,1,10)) ! C212 = O10 = O2*O9
    call matrix_product(rot(1,1,7),rot(1,1,3),rot(1,1,8))  ! C221 = O8 = O7*O3
    call matrix_product(rot(1,1,7),rot(1,1,2),rot(1,1,12)) ! C232 = O12 = O7*O2

    do i=1,12
      rot(1:3,1:3,i+12) = -rot(1:3,1:3,i) ! add inversion
    end do

    end subroutine hexagonal_point_group

! ============================== Added by K. Tagami ============
    subroutine hexagonal_point_group_bravais(nsym,rot)  ! in Bravais Lattice
    implicit none

    integer, intent(out) :: nsym
    real(kind=DP), intent(out) :: rot(3,3,48)

    ! local variables
    real(kind=DP) :: one=1.d0, zero=0.d0
    integer :: i,j

    nsym=24

    rot(1:3,1:3,1:nsym)=zero

    do i=1,3
      rot(i,i,1) = one  ! E = O1
    end do

    ! C6+ = O2
    rot(1,1,2) =  one;    rot(1,2,2) = -one
    rot(2,1,2) =  one;    rot(3,3,2) = one

    !C3+ = O3
    rot(1,2,3) = -one;     rot(2,1,3) = one;
    rot(2,2,3) = -one;     rot(3,3,3) = one

    !C2 = O4
    rot(1,1,4) = -one;     rot(2,2,4) = -one;   rot(3,3,4) = one

    !C3- = O5
    rot(1,1,5) = -one;     rot(1,2,5) = one;
    rot(2,1,5) = -one;     rot(3,3,5) = one

    !C6- = O6
    rot(1,2,6) = one;   rot(2,1,6) = -one;
    rot(2,2,6) = one;    rot(3,3,6) = one

    !C211 = O7
    rot(1,1,7) = -one;   rot(1,2,7) = one
    rot(2,2,7) =  one;   rot(3,3,7) = -one

    !C221 = O8
    rot(1,1,8) = one;  rot(2,1,8) = one
    rot(2,2,8) = -one;  rot(3,3,8) = -one

    !C231 = O9
    rot(1,2,9) = -one;  rot(2,1,9) = -one;  rot(3,3,9) = -one

    !C212 = C10
    rot(1,1,10) = one;  rot(1,2,10) = -one
    rot(2,2,10) = -one;  rot(3,3,10) = -one

    !C222 = C11
    rot(1,1,11) = -one;  rot(2,1,11) = -one
    rot(2,2,11) = one;   rot(3,3,11) = -one

    !C232 = C12
    rot(1,2,12) = one;   rot(2,1,12) = one;     rot(3,3,12) = -one

    do i=1,12
      rot(1:3,1:3,i+12) = -rot(1:3,1:3,i) ! add inversion
    end do

    end subroutine hexagonal_point_group_bravais
! ===============================================================

    subroutine matrix_product(a,b,c)
    implicit none

    real(kind=DP), intent(in) :: a(3,3),b(3,3)
    real(kind=DP), intent(out) :: c(3,3)

    ! local variables
    integer :: i,j

    do j=1,3
      do i=1,3
        c(i,j) = sum(a(i,1:3)*b(1:3,j))
      end do
    end do

    return
    end subroutine matrix_product


    subroutine print_rot(nsym,rot,flat)
    implicit none

    integer, intent(in) :: nsym
    real(kind=DP), intent(in) :: rot(3,3,nsym)
    logical, intent(in) :: flat(48)

    !local variables
    integer :: i, n

    n = 0
    do i=1,nsym
      if(.not.flat(i)) cycle
      n = n+1
      write(nfout,'(5x,3(f10.5,1x))')             rot(1,1:3,i)
      write(nfout,'("n=",i2,1x,3(f10.5,1x))') n, rot(2,1:3,i)
      write(nfout,'(5x,3(f10.5,1x))')             rot(3,1:3,i)
      write(nfout,*)
    end do

    end subroutine print_rot

    subroutine determine_space_group(a1,a2,a3,ntype,natom,itype,apos,nsym,nopr,rot,tau,flat,nprim,tprim,paramset,af)
    implicit none

    real(kind=DP), intent(in) :: a1(3),a2(3),a3(3)
    integer, intent(in) :: ntype,natom,itype(natom)
    real(kind=DP), intent(in) :: apos(natom,3)
    integer, intent(inout) :: nsym
    integer, intent(out) :: nopr(48)
    real(kind=DP), intent(inout) :: rot(3,3,48)
    real(kind=DP), intent(out) :: tau(3,48)
    logical, intent(in) :: flat(48)
    integer, intent(out) :: nprim
    real(kind=DP), intent(out) :: tprim(3,natom)
    logical, intent(in) :: paramset
    logical, intent(in), optional :: af

    ! local variables
! ================================= Usami ===========
!    real(kind=DP) :: ft(3,natom,natom,48)  ! fractional trnsformations
!    integer       :: lp(3,natom,natom,48)  ! lattice point vectors
!    integer :: f(natom,48),fi(natom,48),fp0(natom),fp1(natom) ! atom transformations
    real(kind=DP), allocatable :: ft(:,:,:,:)
    integer, allocatable       :: ift(:,:,:,:)
    integer, allocatable       :: itprim(:,:)
    integer      , allocatable :: lp(:,:,:,:)
    integer, allocatable ::   f(:,:),fi(:,:), fp0(:), fp1(:)
! ===================================================
    real(kind=DP) :: sympos(natom,3)
    integer :: i,j,k,l,n,m,npg
    real(kind=DP) :: x(3),y(3),z(3)
    real(kind=DP) :: a(3,3),ainv(3,3)
    real(kind=DP) :: trans(3)
    real(kind=DP) :: r(3,3),t(3),tmp
    real(kind=DP) :: tt(3)
    integer :: ir(3,3)
    integer :: ntrans(3)
    integer :: isum,isum1,isum0
    logical :: fsym(48)
    integer :: icount
    integer :: iptab(48,48)
    integer :: mpos, isym
    logical :: fpos(48)
    integer :: irot(3,3,48)

    real(kind=DP) :: eps=1.d-8
    real(kind=DP) :: eps2=1.d-6
    real(kind=DP) :: mis2,d(3)
! ============================ added by K. Tagami ==========
!!!!!!    real(kind=DP) :: ktrot(3,3,48) , p2bmat(3,3)
    real(kind=DP) :: ktrot(3,3,48)
    real(kind=DP) :: d1
    real(kind=DP) :: c1, c2, c3
    integer :: ierr1, ierr2, ierr3
    real(kind=DP) :: epsi=1.d+6
    integer :: one = 1
    integer :: ierr
    integer :: isatm,ieatm,iwork
    real(kind=DP) :: elapsed_time
    logical :: read_from_file = .false.
    logical :: logica
    character(len=256) :: tmpstr
    integer :: id_sname = -1
    logical             :: tag_is_found, EOF_reach
    integer,      parameter   :: len_str = 132
    character(len=len_str)       :: str
    logical :: from_af
    from_af = .false.
    if(present(af)) from_af = af
! ============================ Usami ============
    allocate( ft(3,natom,natom,48), lp(3,natom,natom,48) )
    allocate( f(natom,48), fi(natom,48),fp0(natom),fp1(natom) )
! =============================================================
    allocate(ift(3,natom,natom,48))
    allocate(itprim(3,natom))

    mis2 = misalignment**2

    a(1:3,1) = a1(1:3)
    a(1:3,2) = a2(1:3)
    a(1:3,3) = a3(1:3)
    call mat3inv(a,ainv)

!================================ added by Katunori Tagami =========
!!    call set_irot(a,ainv,nsym,rot,irot)
!!!!!!!!!!!!!    Call mat3inv( b2pmat, p2bmat )
!!!!!!    call set_drot( b2pmat, p2bmat, nsym,rot,ktrot )
    call set_drot( kt_b2pmat, kt_p2bmat, nsym,rot,ktrot )
!    if ( system == "cubic" ) then
!       ktrot(:,:,:) = nint(ktrot(:,:,:))
!    endif
! ==================================================================

!    write(910,*) 'DEBUG a =='
!    write(910,*) ktrot(1,1,4), ktrot(1,2,4), ktrot(1,3,4)
!    write(910,*) ktrot(2,1,4), ktrot(2,2,4), ktrot(2,3,4)
!    write(910,*) ktrot(3,1,4), ktrot(3,2,4), ktrot(3,3,4)

    ! make fractional trnsformations, ft(), and lattice point vectors, lp()
    do j=1,natom
      !!$y(1:3) = a1(1:3)*apos(j,1) + a2(1:3)*apos(j,2) + a3(1:3)*apos(j,3)
      y = apos(j,1:3)
      do i=1,natom
        !!$x(1:3) = a1(1:3)*apos(i,1) + a2(1:3)*apos(i,2) + a3(1:3)*apos(i,3)
        x = apos(i,1:3)
        do n=1,nsym
          if(.not.flat(n)) cycle

          !!z(1:3) = rot(1:3,1,n)*x(1)+rot(1:3,2,n)*x(2)+rot(1:3,3,n)*x(3)

!================================ added by Katunori Tagami =========
!!          z(1:3) = matmul(irot(:,:,n),x)
          z(1:3) = matmul( ktrot(:,:,n),x)
! ==================================================================
          !!z(1:3) = y(1:3)-z(1:3)

          ! remove primitive transformation
          !!trans(1:3) = ainv(1:3,1)*z(1)+ainv(1:3,2)*z(2)+ainv(1:3,3)*z(3)

          trans(1:3) = y(1:3)-z(1:3)
          ft(1:3,i,j,n) = trans(1:3) ! in the system a1,a2,a3
          call mod1(ft(1,i,j,n))
          lp(1:3,i,j,n) = nint(trans(1:3)-ft(1:3,i,j,n))
        end do
      end do
    end do

    ! debug
    !!write(nfout,*) 'fractional trasformations and lattice point vectors:'
    !!do j=1,natom
    !!  do i=1,natom
    !!    do n=1,nsym
    !!      if(.not.flat(n)) cycle
    !!      write(nfout,'(3(1x,i2),3(1x,f10.5),3(1x,3i))') i,j,n,ft(1:3,i,j,n),lp(1:3,i,j,n)
    !!    end do
    !!  end do
    !!end do
    ! end debug

    !!$write(*,'("max of isum=",i4)') natom*(natom+1)/2

    if((icond==CONTINUATION .or. icond==COORDINATE_CONTINUATION .or. icond==FIXED_CHARGE_CONTINUATION) &
    &  .and. sw_read_symmetry_from_file == ON) then
      if(m_Files_file_exists(F_CNTN)) then
         call m_Files_open_nfcntn()
         if(mype==0) then
           if(.not. from_af) then
             call rewind_to_tag0(nfcntn,len(tag_symmetry_info),tag_symmetry_info, &
            &     EOF_reach, tag_is_found, str,len_str)
           else
             call rewind_to_tag0(nfcntn,len(tag_symmetry_info_af),tag_symmetry_info_af, &
            &     EOF_reach, tag_is_found, str,len_str)
           endif
           if(tag_is_found) then
             read(nfcntn,*,err=50) nprim,nsym
             do i=1,nprim
                read(nfcntn,*,err=50) tprim(1:3,i)
             enddo
             do i=1,nsym
                read(nfcntn,*,err=50) tau(1:3,i),fsym(i)
             enddo
           endif
         endif
         call mpi_bcast(tag_is_found,1,mpi_logical,0,MPI_CommGroup,ierr)
         if(tag_is_found) then
           call mpi_bcast(nprim,1,mpi_integer,0,MPI_CommGroup,ierr)
           call mpi_bcast(nsym,1,mpi_integer,0,MPI_CommGroup,ierr)
           call mpi_bcast(tprim,3*natom,mpi_double_precision,0,MPI_CommGroup,ierr)
           call mpi_bcast(tau,3*nsym,mpi_double_precision,0,MPI_CommGroup,ierr)
           call mpi_bcast(fsym,nsym,mpi_logical,0,MPI_CommGroup,ierr)
           if(printable) write(nfout,'(a)') ' !** successfully read symmetry-info from file'
           read_from_file = .true.
         endif
50       continue
         if(read_from_file .and. .not.paramset)then
           if(.not. from_af) then
             if(allocated(mtprim)) deallocate(mtprim)
             if(allocated(mtau))   deallocate(mtau)
             if(allocated(mfsym))  deallocate(mfsym)
             allocate(mtprim(3,nprim))
             allocate(mtau(3,48))
             allocate(mfsym(48))
             mnprim = nprim;mnsym = nsym
             mtprim(1:3,1:mnprim) = tprim(1:3,1:mnprim)
             mtau(:,:) = tau(:,:)
             mfsym(:) = fsym(:)
           else
             if(allocated(mtprim_af)) deallocate(mtprim_af)
             if(allocated(mtau_af))   deallocate(mtau_af)
             if(allocated(mfsym_af))  deallocate(mfsym_af)
             allocate(mtprim_af(3,nprim))
             allocate(mtau_af(3,48))
             allocate(mfsym_af(48))
             mnprim_af = nprim;mnsym_af = nsym
             mtprim_af(1:3,1:mnprim_af) = tprim(1:3,1:mnprim_af)
             mtau_af(:,:) = tau(:,:)
             mfsym_af(:) = fsym(:)
           endif
         endif
      endif
    endif

    if(sw_non_symmorphic == ON .and. .not. cpumax_reached .and. .not. read_from_file) then
      if(printable) write(nfout,'(a)') &
      ' !** sw_non_symmorphic == ON; all space group will be taken into account'
      npg=0
      fsym(1:nsym)=.false.
      f(1:natom,1:nsym) = 0
      tau(1:3,1:nsym) = 0.d0
      nprim = 0
      ift = nint(ft*epsi)
      itprim = 0
      n=1
      iwork = (natom-1)/npes + 1
      isatm = min(mype*iwork+1, natom+1)
      ieatm = min(isatm+iwork-1,natom)
      call m_Timing_initialize_tmp_timer()
      if(flat(n)) then
        loop_l: do l = 1, natom
          elapsed_time = m_Timing_get_elapsed_time_tmp()
          if(cpumax_symsearch>0 .and. elapsed_time>cpumax_symsearch) then
             cpumax_reached = .true.
             exit loop_l
          endif
          loop_k: do k = 1, natom
            if(itype(k).ne.itype(l)) cycle loop_k
            do m=1,nprim
              if(sum(abs(ift(1:3,k,l,n)-itprim(1:3,m)))<=one) then
                 cycle loop_k
              end if
            end do
            isum=0
            isum1=0
            loop_j: do j=isatm,ieatm
              loop_i: do i=1,natom
                if(itype(i).ne.itype(j)) cycle loop_i
                if(abs(ift(1,i,j,n)-ift(1,k,l,n)) > one .or. &
                 & abs(ift(2,i,j,n)-ift(2,k,l,n)) > one .or. &
                 & abs(ift(3,i,j,n)-ift(3,k,l,n)) > one ) then
                   cycle loop_i
                end if
                isum = isum + i
              end do loop_i
            end do loop_j
            call mpi_allreduce(isum,isum1,1,mpi_integer,mpi_sum,MPI_CommGroup,ierr)
            if(isum1 .eq. natom*(natom+1)/2) then
               nprim = nprim + 1
               tprim (1:3,nprim) = ft (1:3,k,l,n)
               itprim(1:3,nprim) = ift(1:3,k,l,n)
            end if
          enddo loop_k
        enddo loop_l
      endif

      if(.not. cpumax_reached) then
        loop_sym: do n=1,nsym
          if(mod(n-1,npes) /= mype) cycle
          if(.not.flat(n)) cycle
          loop_l2: do l=1,natom
            elapsed_time = m_Timing_get_elapsed_time_tmp()
            if(cpumax_symsearch>0 .and. elapsed_time>cpumax_symsearch) then
               cpumax_reached = .true.
               exit loop_sym
            endif
            loop_k2: do k=1,natom
              if(itype(k).ne.itype(l)) cycle loop_k2
              isum=0
              loop_j2: do j=1,natom
                loop_i2: do i=1,natom
                  if(itype(i).ne.itype(j)) cycle loop_i2
                  if(abs(ift(1,i,j,n)-ift(1,k,l,n)) > one .or. &
                   & abs(ift(2,i,j,n)-ift(2,k,l,n)) > one .or. &
                   & abs(ift(3,i,j,n)-ift(3,k,l,n)) > one ) then
                     cycle loop_i2
                  end if
                  isum = isum + i
                  if(isum .eq. natom*(natom+1)/2) then
                     fsym(n) = .true.
                     npg = npg + 1
                     tau(1:3,n) = ft(1:3,k,l,n)
                     exit loop_l2
                  end if
                end do loop_i2
              end do loop_j2
            end do loop_k2
          end do loop_l2
        end do loop_sym
        call mpi_allreduce(mpi_in_place,tau,3*nsym,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
        call mpi_allreduce(mpi_in_place,fsym,nsym,mpi_logical,mpi_lor,MPI_CommGroup,ierr)
      endif
      call m_Timing_finalize_tmp_timer

      if(.not.paramset) then
        if(.not. from_af) then
          if(allocated(mtprim)) deallocate(mtprim)
          if(allocated(mtau))   deallocate(mtau)
          if(allocated(mfsym))  deallocate(mfsym)
          allocate(mtprim(3,nprim))
          allocate(mtau(3,48))
          allocate(mfsym(48))
          mnprim = nprim;mnsym = nsym
          mtprim(1:3,1:mnprim) = tprim(1:3,1:mnprim)
          mtau(:,:) = tau(:,:)
          mfsym(:) = fsym(:)
        else
          if(allocated(mtprim_af)) deallocate(mtprim_af)
          if(allocated(mtau_af))   deallocate(mtau_af)
          if(allocated(mfsym_af))  deallocate(mfsym_af)
          allocate(mtprim_af(3,nprim))
          allocate(mtau_af(3,48))
          allocate(mfsym_af(48))
          mnprim_af = nprim;mnsym_af = nsym
          mtprim_af(1:3,1:mnprim_af) = tprim(1:3,1:mnprim)
          mtau_af(:,:) = tau(:,:)
          mfsym_af(:) = fsym(:)
        endif
      endif

      do i=1,nprim
         do j=1,3
            tprim(j,i) = rational(tprim(j,i))
         end do
      end do

      if(ipri>=1.and.nprim > 1) then
         write(nfout,'("=== Primitive translations ===")')
         do i=2,nprim
            write(nfout,'(i4,3(1x,f10.5),1x,a)') i-1,tprim(1:3,i),get_lattice_point_name(tprim(1,i))
         end do
      end if

    endif

    if(sw_non_symmorphic == OFF .or. cpumax_reached) then
      fsym  = .false.
      tau   = 0.d0
      tprim = 0.d0
      nprim = 1
      if(printable .and. sw_non_symmorphic==OFF) write(nfout,'(a)') &
      ' !*** sw_non_symmorphic == OFF; only point group will be taken into account'
      do n = 1, nsym
         do j = 1, natom
            if(sum(ft(1:3,1,j,n))<eps2) then
               fsym(n) = .true.
             cycle
            endif
         enddo
      enddo
    endif

    if(cpumax_reached) then
       write(nfout,'(a,f10.2,a)') ' !** the time spent in automatic search of space group symmetry exceeded ', &
                        & cpumax_symsearch,' (s), and will be aborted.'
    endif


    i=0
! =================================== Modified by K. Tagami =========
!    do n=1,nsym
!       if(.not.fsym(n)) cycle
!       i=i+1
!       nopr(i) = n
!       rot(1:3,1:3,i) = rot(1:3,1:3,n)
!       !!$tau(1:3,i) = tau(1:3,n)
!       do j=1,3
!          tau(j,i) = rational(tau(j,n))
!       end do
!    end do
! ================================
    do n=1,nsym
       if(.not.fsym(n)) cycle
!
       call ktrational( tau(1,n), c1, ierr1 )
       call ktrational( tau(2,n), c2, ierr2 )
       call ktrational( tau(3,n), c3, ierr3 )
!
       if ( ierr1+ierr2+ierr3 < 0 ) then
          fsym(n) = .false.
          cycle
       endif
!
       i=i+1
       nopr(i) = n
       rot(1:3,1:3,i) = rot(1:3,1:3,n)

       tau(1,i) = c1;        tau(2,i) = c2;        tau(3,i) = c3
    end do
! =========================================================================
    nsym = i
    if(nsym == 1) return

    call set_product_table(iptab,nprim,tprim)

    loop_0: do i=2,nsym
       fsym(1:nsym) = .false.
       do n=1,nsym
          if(iptab(n,i) > 0) fsym(n) = .true.
       end do
       do n=2,nsym
          if(.not.fsym(n)) cycle
          loop_2: do m=2,nsym
             if(.not.fsym(m)) cycle
             do l=1,nsym
                if(.not.fsym(l)) cycle
                if(l == iptab(n,m)) exit loop_2
             end do
             cycle loop_0
          end do loop_2
       end do
       fpos(i) = .true.
    end do loop_0

    mpos = 0
    do i=2,nsym
       if(.not.fpos(i)) cycle
       m = 0
       do n=1,nsym
          if(iptab(n,i)>0) m = m + 1
       end do
       if(m>mpos) isym = i
    end do

    i=0
    do n=1,nsym
       if(iptab(n,isym)<1) cycle
       i=i+1
       nopr(i) = nopr(n)
       rot(1:3,1:3,i) = rot(1:3,1:3,n)
       do m=1,nprim
          tt(1:3) = tau(1:3,n)-tprim(1:3,m)
          call mod1(tt)
          if(sum(abs(tt(1:3))) < eps2) then
          ! experimental
          !!$d = matmul(altv,tt)
          !!$if(sum(d(1:3)**2) < mis2) then
          ! end experimental
            tau(1:3,i) = 0.d0
          else
            tau(1:3,i) = tau(1:3,n)
          end if
       end do
    end do
    nsym = i

    !!$call set_product_table(iptab,nprim,tprim)
    !!$stop 'determine_space_group: debug'
! =============================================== Usami ==========
    deallocate( ft,lp )
    deallocate( f, fi, fp0, fp1 )
! ==============================================================
    end subroutine determine_space_group

    subroutine set_irot(a,ainv,nsym,rot,irot)
      real(kind=DP), intent(in) :: a(3,3),ainv(3,3)
      integer, intent(in) :: nsym
      real(kind=DP), intent(in) :: rot(3,3,nsym)
      integer, intent(out) :: irot(3,3,nsym)

      integer :: n
      do n=1,nsym
         irot(1:3,1:3,n) = nint(matmul(ainv,matmul(rot(:,:,n),a)))
      end do
    end subroutine set_irot

! ================================ Added by K. Tagami =========
    subroutine set_drot( a,ainv,nsym,rot,drot )
      real(kind=DP), intent(in) :: a(3,3),ainv(3,3)
      integer, intent(in) :: nsym
      real(kind=DP), intent(in) :: rot(3,3,nsym)
      real(kind=DP), intent(out) :: drot(3,3,nsym)

      integer :: n
      do n=1,nsym
         drot(1:3,1:3,n) = matmul(ainv,matmul(rot(:,:,n),a))
      end do
    end subroutine set_drot
! ==============================================================

    subroutine set_product_table(iptab,nprim,tprim)
      integer, intent(out) :: iptab(48,48)
      integer, intent(in)  :: nprim
      real(kind=DP), intent(in) :: tprim(3,natom)

      integer :: i,j,k,l,m,n
! ============================= Added by K. Tagami ==============
!      integer :: irot(3,3,48),ir(3,3)
!!!!!!!!      real(kind=DP) :: drot(3,3,48), dr(3,3),  p2bmat(3,3)
      real(kind=DP) :: drot(3,3,48), dr(3,3)
! ===============================================================
      real(kind=DP) :: t(3),tt(3),tmp

      real(kind=DP), parameter :: eps  = 1.d-8
      real(kind=DP), parameter :: eps2 = 1.d-6
      real(kind=DP) :: mis2,d(3)

      mis2 = misalignment**2

      if(ipri>1) then
! =================================== Added by K. Tagami ======
!         write(nfout,'("irot:")')
         write(nfout,'("drot:")')
! ====================================================================
      end if
! ================================ Added by K. Tagami =========
!!!!!!!!      Call mat3inv( b2pmat, p2bmat )
!!!      call set_drot( b2pmat, p2bmat, nsym,rot, drot )
      call set_drot( kt_b2pmat, kt_p2bmat, nsym,rot, drot )
! ============================================================
      do i=1,nsym
!         do n=1,3
!            do m=1,3
!               tmp = 0.d0
!               do k=1,3
!                  do l=1,3
!                     tmp = tmp + rltv(k,m)*rot(k,l,i)*altv(l,n)
!                  end do
!               end do
!!!               irot(m,n,i) = nint(tmp/PAI2)
!            end do
!         end do
         if(ipri>1) then
! =================================== Added by K. Tagami ======
!            write(nfout,'(10(1x,i2))') i,irot(1:3,1:3,i)
             write(nfout,'(i2,9F8.4)') i, drot(1:3,1:3,i)
! ====================================================================
         end if
      end do

      do j=1,nsym
         LOOP: do i=1,nsym
            do k=1,3
! =================================== Added by K. Tagami ======
!!!!!               t(k) = sum(irot(k,1:3,i)*tau(1:3,j)) + tau(k,i)
               t(k) = sum( drot(k,1:3,i)*tau(1:3,j) ) + tau(k,i)
! ====================================================================
               do l=1,3
! =================================== Added by K. Tagami ======
!!!                  ir(k,l) = sum(irot(k,1:3,i)*irot(1:3,l,j))
                  dr(k,l) = sum( drot(k,1:3,i)*drot(1:3,l,j) )
! ====================================================================
               end do
            end do
            call mod1(t)
            ! debug
            !write(nfout,'(11(1x,i2),3(1x,f10.5))') i,j,ir(1:3,1:3),t(1:3)
            ! end debug
            do k=1,nsym
! =================================== Added by K. Tagami ======
!!!               if(sum(abs(ir(1:3,1:3)-irot(1:3,1:3,k))) < eps2) then
               if ( sum(abs(dr(1:3,1:3)-drot(1:3,1:3,k))) < eps2) then
! ====================================================================
                do l=1,nprim
                   tt(1:3) = t(1:3)-tau(1:3,k)-tprim(1:3,l)
                   call mod1(tt)
                   if(sum(abs(tt(1:3))) < eps2) then
                   ! experimental
                   !!$d = matmul(altv,tt)
                   !!$if(sum(d(1:3)**2) < mis2) then
                   ! end experimental
                      iptab(i,j) = k
                      cycle LOOP
                   end if
                 end do
               end if
            end do
            !!$write(nfout,'("Missing operation?")')
            !!$write(nfout,'(11(1x,i2),3(1x,f10.5))') i,j,ir(1:3,1:3),t(1:3)
            !!$stop 'set_product_table: missing operation'
            iptab(i,j) = 0
         end do LOOP
      end do

      if(ipri>1) then
         write(nfout,'("Product table:")')
         do i=1,nsym
            write(nfout,'(48(1x,i2))') iptab(i,1:nsym)
         end do
      end if
            !!$stop 'set_product_table: debug'

    end subroutine set_product_table

    function get_lattice_point_name(tprim) result(name)
      real(kind=DP), intent(in) :: tprim(3)
      character(15) :: name

      real(kind=DP) :: eps = 1.d-6
      real(kind=DP) :: mis2,d(3)

      mis2 = misalignment**2

      name = ''
      if(abs(tprim(1)-0.5d0) < eps .and. &
       & abs(tprim(2)-0.5d0) < eps .and. &
       & abs(tprim(3)) < eps ) then
         name = 'c-face centered'
      else if(abs(tprim(1)-0.5d0) < eps .and. &
       & abs(tprim(2)) < eps .and. &
       & abs(tprim(3)-0.5d0) < eps ) then
         name = 'b-face centered'
      else if(abs(tprim(1)) < eps .and. &
       & abs(tprim(2)-0.5d0) < eps .and. &
       & abs(tprim(3)-0.5d0) < eps ) then
         name = 'a-face centered'
      else if(abs(tprim(1)-0.5d0) < eps .and. &
       & abs(tprim(2)-0.5d0) < eps .and. &
       & abs(tprim(3)-0.5d0) < eps ) then
         name = 'body centered'
      end if

      ! experimental
      !!$name = 'c-face centered'
      !!$d = (/ 0.5d0, 0.5d0, 0.d0 /)
      !!$d = matmul(altv,tprim-d)
      !!$if(sum(d(1:3)**2) < mis2) return
      !!$name = 'b-face centered'
      !!$d = (/ 0.5d0, 0.d0, 0.5d0 /)
      !!$d = matmul(altv,tprim-d)
      !!$if(sum(d(1:3)**2) < mis2) return
      !!$name = 'a-face centered'
      !!$d = (/ 0.d0, 0.5d0, 0.5d0 /)
      !!$d = matmul(altv,tprim-d)
      !!$if(sum(d(1:3)**2) < mis2) return
      !!$name = 'body centered'
      !!$d = (/ 0.5d0, 0.5d0, 0.5d0 /)
      !!$d = matmul(altv,tprim-d)
      !!$if(sum(d(1:3)**2) < mis2) return
      !!$name = ''
      ! end experimental

    end function get_lattice_point_name

  end subroutine space_group

  subroutine print_prim_trans(nfout,ipri,nprim,tprim)
    integer, intent(in) :: nfout,ipri,nprim
    real(kind=DP) :: tprim(3,nprim)

    integer :: i

    if(ipri>0 .and. nprim > 1) then
       write(nfout,'("=== Primitive translations ===")')
       do i=2,nprim
          write(nfout,'(i4,3(1x,f10.5),1x,a)') i-1,tprim(1:3,i),get_lattice_point_name(tprim(1,i))
       end do
    end if
  contains
    function get_lattice_point_name(tprim) result(name)
      real(kind=DP), intent(in) :: tprim(3)
      character(15) :: name

      real(kind=DP) :: eps = 1.d-6

      name = ''
      if(abs(tprim(1)-0.5d0) < eps .and. &
       & abs(tprim(2)-0.5d0) < eps .and. &
       & abs(tprim(3)) < eps ) then
         name = 'c-face centered'
      else if(abs(tprim(1)-0.5d0) < eps .and. &
       & abs(tprim(2)) < eps .and. &
       & abs(tprim(3)-0.5d0) < eps ) then
         name = 'b-face centered'
      else if(abs(tprim(1)) < eps .and. &
       & abs(tprim(2)-0.5d0) < eps .and. &
       & abs(tprim(3)-0.5d0) < eps ) then
         name = 'a-face centered'
      else if(abs(tprim(1)-0.5d0) < eps .and. &
       & abs(tprim(2)-0.5d0) < eps .and. &
       & abs(tprim(3)-0.5d0) < eps ) then
         name = 'body centered'
      end if

    end function get_lattice_point_name

  end subroutine print_prim_trans

  subroutine print_space_group(nfout,system,nsym,nopr,rot,tau)
    implicit none

    integer, intent(in) :: nfout
    character(len=9), intent(in) :: system
    integer, intent(in) :: nsym
    integer, intent(in) :: nopr(nsym)
    real(kind=DP), intent(in) :: rot(3,3,nsym)
    real(kind=DP), intent(in) :: tau(3,nsym)
! ============================ added by K. Tagami =======
    real(kind=DP) :: taub(3,nsym), p2bmat(3,3), c1
! ============================================================

    ! local variables
    integer :: i

    if(.not.printable) return
! ===================================== added by K. Tagami ======
!!!!!!!!!!!!!!    taub = matmul( b2pmat, tau )
    taub = matmul( kt_b2pmat, tau )
! ~==============================================================
    if(system .eq. 'cubic') then
      do i=1,nsym
! ===================================== added by K. Tagami ======
!        write(nfout,'(8x,3(1x,f10.5),2x,f10.5)') rot(1,1:3,i),tau(1,i)
!        write(nfout,'(i2,1x,a5,3(1x,f10.5),2x,f10.5)') i,oh(nopr(i)),rot(2,1:3,i),tau(2,i)
!        write(nfout,'(8x,3(1x,f10.5),2x,f10.5)') rot(3,1:3,i),tau(3,i)
! ==
        write(nfout,'(8x,3(1x,f10.5),2x,f10.5)') rot(1,1:3,i),taub(1,i)
        write(nfout,'(i2,1x,a5,3(1x,f10.5),2x,f10.5)') i,oh(nopr(i)),rot(2,1:3,i),taub(2,i)
        write(nfout,'(8x,3(1x,f10.5),2x,f10.5)') rot(3,1:3,i),taub(3,i)
! ========================================================
        write(nfout,*)
      end do
    else
      do i=1,nsym
! ===================================== added by K. Tagami ======
!        write(nfout,'(8x,3(1x,f10.5),2x,f10.5)') rot(1,1:3,i),tau(1,i)
!        write(nfout,'(i2,1x,a5,3(1x,f10.5),2x,f10.5)') i,d6h(nopr(i)),rot(2,1:3,i),tau(2,i)
!        write(nfout,'(8x,3(1x,f10.5),2x,f10.5)') rot(3,1:3,i),tau(3,i)
!
        write(nfout,'(8x,3(1x,f10.5),2x,f10.5)') rot(1,1:3,i),taub(1,i)
        write(nfout,'(i2,1x,a5,3(1x,f10.5),2x,f10.5)') i,d6h(nopr(i)),rot(2,1:3,i),taub(2,i)
        write(nfout,'(8x,3(1x,f10.5),2x,f10.5)') rot(3,1:3,i),taub(3,i)
! ===============================================================
        write(nfout,*)
      end do
    end if
  end subroutine print_space_group


  subroutine set_af_operator(nfout,ipri,system &
                            & ,nsym_af,nopr_af,rot_af,tau_af,nprim_af,tprim_af &
                            & ,nsym,nopr,rot,tau,nprim,tprim)
    implicit none
    integer, intent(in) :: nfout,ipri
    character(len=9), intent(in) :: system
    integer, intent(in) :: nsym_af,nsym,nprim_af,nprim
    integer, intent(in) :: nopr_af(48)
    real(kind=DP), intent(in) :: rot_af(3,3,49),tau_af(3,49)
    real(kind=DP), intent(in) :: tprim_af(3,nprim_af),tprim(3,nprim)
    real(kind=DP), intent(inout) :: rot(3,3,48),tau(3,48)
    integer, intent(out) :: nopr(48)

    integer :: i,j,k,iaf,n
    real(kind=DP) :: r(3,3,24),t(3,24)
    logical :: faf(nsym_af),flag(nsym_af)
    real(kind=DP), parameter :: eps=1.d-8

    faf(1:nsym_af) = .true.
    n=0
    do iaf=1,nsym_af
      do i=1,nsym
         if(sum(abs(rot_af(1:3,1:3,iaf)-rot(1:3,1:3,i))) < eps &
           & .and. sum(abs(tau_af(1:3,iaf)-tau(1:3,i))) < eps) then
           n = n + 1
           faf(iaf) = .false.
         end if
      end do
    end do
    if(nsym_af .ne. 2*n) then
       if(nprim_af .gt. nprim) then
          iaf = 1
          do i=1,nprim
             do j=1,nprim_af
                if(sum(abs(tprim(1:3,i)-tprim_af(1:3,j))) > eps) then
                   nopr(nsym+1) = nopr_af(iaf)
                   rot(:,:,nsym+1) = rot_af(:,:,iaf)
                   tau(:,nsym+1) = tprim_af(:,j)
                   go to 1000
                end if
             end do
          end do
          call phase_error_with_msg(nfout, 'AF translation was not found.',__LINE__,__FILE__)
       else
          call phase_error_with_msg(nfout, 'nsym_af .ne. 2*n and nprim_af .le. nprim',__LINE__,__FILE__)
       end if
1000   continue
    else
       LOOP: do iaf=1,nsym_af
          if(faf(iaf)) then
             do i=1,nsym
                do j=1,3
                   do k=1,3
                      r(j,k,i) = sum(rot_af(j,:,iaf)*rot(:,k,i))
                   end do
                   t(j,i) = sum(rot_af(j,:,iaf)*tau(:,i)) + tau_af(j,iaf)
                end do
             end do
             flag(1:nsym_af) = .true.
             do i=1,nsym
                do j=1,nsym_af
                   if(faf(j)) then
                      if(sum(abs(rot_af(1:3,1:3,iaf)-r(1:3,1:3,i))) < eps &
                   & .and. sum(abs(tau_af(1:3,iaf)-t(1:3,i))) < eps) then
                         flag(j) = .false.
                      end if
                   end if
                end do
             end do
             do j=1,nsym_af
                if(faf(j).and.flag(j)) cycle LOOP
             end do
             nopr(nsym+1) = nopr_af(iaf)
             rot(:,:,nsym+1) = rot_af(:,:,iaf)
             tau(:,nsym+1) = tau_af(:,iaf)
             if(ipri>0) then
                write(nfout,'(" af operator is found. iaf = ",i8)') iaf
             end if
             exit LOOP
          end if
       end do LOOP
    end if

    if(ipri>0) then
      write(nfout,'("AF operator:")')
      call print_space_group(nfout,system,1,nopr(nsym+1),rot(1,1,nsym+1),tau(1,nsym+1))
    end if

  end subroutine set_af_operator

! =============================== Added by K. Tagami ==== 1.0 ==
  subroutine set_af_operator_kt(nfout,ipri,system &
                       & ,nsym_af,nopr_af,rot_af,tau_af,nprim_af,tprim_af &
                       & ,nsym,nopr,rot,tau,nprim,tprim )
    implicit none
    integer, intent(in) :: nfout,ipri
    character(len=9), intent(in) :: system
    integer, intent(in) :: nsym_af,nsym,nprim_af,nprim
    integer, intent(in) :: nopr_af(49)
    real(kind=DP), intent(in) :: rot_af(3,3,49),tau_af(3,49)
    real(kind=DP), intent(in) :: tprim_af(3,nprim_af),tprim(3,nprim)
    real(kind=DP), intent(inout) :: rot(3,3,49),tau(3,49)
    integer, intent(out) :: nopr(49)
    integer, parameter :: DEBUGPRINTLEVEL = 2

    integer :: i,j,k,iaf,n
    real(kind=DP) :: r(3,3,24),t(3,24)
    logical :: faf(nsym_af),flag(nsym_af)
    real(kind=DP), parameter :: eps=1.d-8
! ------------------
    integer :: jfound
!!$    real(kind=DP) :: c1, kt_i, kt_j
    integer ::  kt_i, kt_j
! -------------------

    faf(1:nsym_af) = .true.
    n=0
    do iaf=1,nsym_af
      do i=1,nsym
         if(sum(abs(rot_af(1:3,1:3,iaf)-rot(1:3,1:3,i))) < eps &
           & .and. sum(abs(tau_af(1:3,iaf)-tau(1:3,i))) < eps) then
           n = n + 1
           faf(iaf) = .false.
         end if
      end do
    end do
    if(ipri>=DEBUGPRINTLEVEL) then
       do iaf = 1, nsym_af
          if(faf(iaf)) then
             write(nfout,'(" faf(",i3,") = .true.")') iaf
          else
             write(nfout,'(" faf(",i3,") = .false.")') iaf
          end if
       end do
    end if
    if(nsym_af .ne. 2*n) then
       if(nprim_af .gt. nprim) then
          iaf = 1
! ------------------------------------------
          do j=1,nprim_af
             jfound = 0
             do i=1,nprim
                if( abs(tprim(1,i)-tprim_af(1,j)) < eps &
                   .and.abs(tprim(2,i)-tprim_af(2,j)) < eps &
                   .and.abs(tprim(3,i)-tprim_af(3,j)) < eps ) then
                   jfound = 1
                endif
             end do
             if ( jfound == 0 ) then
                nopr(nsym+1) = nopr_af(iaf)
                rot(:,:,nsym+1) = rot_af(:,:,iaf)
                tau(:,nsym+1) = tprim_af(:,j)
                go to 1000
             end if
          end do
! ------------------------------------------
          call phase_error_with_msg(nfout, 'AF translation was not found.',__LINE__,__FILE__)
       else
          call phase_error_with_msg(nfout, 'nsym_af .ne. 2*n and nprim_af .le. nprim',__LINE__,__FILE__)
       end if
1000   continue
    else
! --------------------- first try : E-symmetry ----
      iaf = 1
      do kt_j=1,nprim_af
        jfound = 0
        do kt_i=1,nprim
          if( abs(tprim(1,kt_i)-tprim_af(1,kt_j)) < eps &
              .and.abs(tprim(2,kt_i)-tprim_af(2,kt_j)) < eps &
              .and.abs(tprim(3,kt_i)-tprim_af(3,kt_j)) < eps ) then
            jfound = 1
          endif
        end do
        if ( jfound == 0 ) then
          nopr(nsym+1) = nopr_af(iaf)
          rot(:,:,nsym+1) = rot_af(:,:,iaf)
          tau(:,nsym+1) = tprim_af(:,kt_j)
          goto 2000
        end if
      end do
! ------------------------- next try --------
       LOOP: do iaf=1, nsym_af
          if(faf(iaf)) then
             do i=1, nsym
                do j=1,3
                   do k=1,3
                      r(j,k,i) = sum(rot_af(j,:,iaf)*rot(:,k,i))
                   end do
                   t(j,i) = sum(rot_af(j,:,iaf)*tau(:,i)) + tau_af(j,iaf)
                end do
             end do
             flag(1:nsym_af) = .true.
             do i=1, nsym
                do j=1, nsym_af
                   if(faf(j)) then
                      if(sum(abs(rot_af(1:3,1:3,iaf)-r(1:3,1:3,i))) < eps &
                   & .and. sum(abs(tau_af(1:3,iaf)-t(1:3,i))) < eps) then
                         flag(j) = .false.
                      end if
                   end if
                end do
             end do
             do j=1, nsym_af
                if(faf(j).and.flag(j)) cycle LOOP
             end do
             nopr(nsym+1) = nopr_af(iaf)
             rot(:,:,nsym+1) = rot_af(:,:,iaf)
             tau(:,nsym+1) = tau_af(:,iaf)
             exit LOOP
          end if
       end do LOOP
! ----------------------------------------------
2000  continue
    end if
    if(ipri>0) then
      write(nfout,'("AF operator:")')

      call print_space_group(nfout,system,1,nopr(nsym+1),rot(1,1,nsym+1),tau(1,nsym+1))
    end if

  end subroutine set_af_operator_kt
! =================================================================


  subroutine set_tspace_generators(nfout,ipri,nsym,nopr,rot,tau,natom,nprim,tprim)
    implicit none
    integer, intent(in) :: nfout,ipri
    integer, intent(in) :: nsym,nopr(49)
    real(kind=DP), intent(in) :: rot(3,3,49),tau(3,49)
    integer, intent(in) :: natom,nprim
    real(kind=DP), intent(in) :: tprim(3,natom)

    integer :: iptab(nsym,nsym)
    integer :: ng,ig(3)

    call set_product_table
    call search_generators
    call set_generators
! ============================ K. Tagami ==== Unknown ============
!!!!!    if(nprim>1) call set_translations_zero
! ===========================================================

  contains
    subroutine set_product_table
      integer :: i,j,k,l,m,n
! =============================== Added by K. Tagami ========
!      integer :: irot(3,3,48),ir(3,3)
!!!!!!!!!      real(kind=DP) :: drot(3,3,48),dr(3,3), p2bmat(3,3)
      real(kind=DP) :: drot(3,3,48),dr(3,3)
! ===========================================================
      real(kind=DP) :: t(3),tt(3),tmp
      real(kind=DP), parameter :: eps = 1.d-6

      if(ipri>1) then
! =============================== Added by K. Tagami ========
!         write(nfout,'("irot:")')
         write(nfout,'("drot2:")')
! ===========================================================
      end if
! ============================== Added by K. Tagami ========
!!!!!!!!!!      Call mat3inv( b2pmat, p2bmat )
! =========================================================
      do i=1,nsym
         do n=1,3
            do m=1,3
               tmp = 0.d0
               do k=1,3
                  do l=1,3
! ============================== Added by K. Tagami ===============
!                     tmp = tmp + rltv(k,m)*rot(k,l,i)*altv(l,n)
!!!!!!!!                      tmp = tmp + p2bmat(m,k)*rot(k,l,i)*b2pmat(l,n)
                      tmp = tmp + kt_p2bmat(m,k)*rot(k,l,i)*kt_b2pmat(l,n)
! ================================================================
                  end do
               end do
! ============================== Added by K. Tagami ===============
!!               irot(m,n,i) = nint(tmp/PAI2)
               drot(m,n,i) = tmp
! ================================================================
            end do
         end do
         if(ipri>1) then
! ============================== Added by K. Tagami ===============
!            write(nfout,'(10(1x,i2))') i,irot(1:3,1:3,i)
            write(nfout,'(i2,9F8.4)') i, drot(1:3,1:3,i)
! ================================================================
         end if
      end do

      do j=1,nsym
         LOOP: do i=1,nsym
            do k=1,3
! =============================== Added by K. Tagami ======
!!!               t(k) = sum(irot(k,1:3,i)*tau(1:3,j)) + tau(k,i)
               t(k) = sum( drot(k,1:3,i)*tau(1:3,j) ) + tau(k,i)
! =========================================================
               do l=1,3
! =============================== Added by K. Tagami ======
!!!                  ir(k,l) = sum(irot(k,1:3,i)*irot(1:3,l,j))
                  dr(k,l) = sum( drot(k,1:3,i)*drot(1:3,l,j) )
! =========================================================
               end do
            end do
            call mod1(t)
            ! debug
            !write(nfout,'(11(1x,i2),3(1x,f10.5))') i,j,ir(1:3,1:3),t(1:3)
            ! end debug
            do k=1,nsym
! =============================== Added by K. Tagami ======
!!!               if(sum(abs(ir(1:3,1:3)-irot(1:3,1:3,k))) < eps) then
               if(sum(abs( dr(1:3,1:3)-drot(1:3,1:3,k))) < eps) then
! =========================================================
                do l=1,nprim
                   tt(1:3) = t(1:3)-tau(1:3,k)-tprim(1:3,l)
                   call mod1(tt)
                   if(sum(abs(tt(1:3))) < eps) then
                      iptab(i,j) = k
                      cycle LOOP
                   end if
                 end do
               end if
            end do
            write(nfout,'("Missing operation?")')
! =============================== Added by K. Tagami ======
!            write(nfout,'(11(1x,i2),3(1x,f10.5))') i,j,ir(1:3,1:3),t(1:3)
            write(nfout,'(2(1x,i2),12F8.4)') i,j, dr(1:3,1:3),t(1:3)
! =========================================================
            call phase_error_with_msg(nfout, 'set_product_table: missing operation',__LINE__,__FILE__)
            !!$iptab(i,j) = 0
         end do LOOP
      end do

      if(ipri>1) then
         write(nfout,'("Product table:")')
         do i=1,nsym
            write(nfout,'(48(1x,i2))') iptab(i,1:nsym)
         end do
      end if
            !!$stop 'set_product_table: debug'

    end subroutine set_product_table

    subroutine search_generators
      integer :: i,j,k
      logical :: fexist(1:nsym)
      integer :: ii
               ! debug for the function `group_is_the_same'
               !!fexist(1:nsym) = .false.
               !!fexist(5) = .true.
               !!fexist(19) = .true.
               !!fexist(25) = .true.
               !!if(group_is_the_same(fexist)) then
               !!end if
               !!  do ii=1,nsym
               !!    write(nfout,*) ii,fexist(ii)
               !!  end do
               ! end debug
      ng = 1
      do i=1,nsym
         fexist(1:nsym) = .false.
         fexist(i) = .true.
         if(group_is_the_same(fexist)) then
            ig(ng) = i
            return
         end if
         ! debug
         !!write(nfout,'(3(i2,1x),48l)') i,fexist(1:nsym)
         ! end  debug
      end do
      ng = 2
      do i=2,nsym
         do j=i+1,nsym
            fexist(1:nsym) = .false.
            fexist(i) = .true.
            fexist(j) = .true.
            if(group_is_the_same(fexist)) then
               ig(1)  = i
               ig(ng) = j
               return
            end if
            ! debug
            !!write(nfout,'(3(i2,1x),48l)') i,j,fexist(1:nsym)
            ! end  debug
         end do
      end do
      ng = 3
      do i=2,nsym
         do j=i+1,nsym
            do k=j+1,nsym
               fexist(1:nsym) = .false.
               fexist(i) = .true.
               fexist(j) = .true.
               fexist(k) = .true.
               if(group_is_the_same(fexist)) then
                  ig(1)  = i
                  ig(2)  = j
                  ig(ng) = k
                  return
               end if
               ! debug
               !!!$write(nfout,'(3(i2,1x),48l)') i,j,k,fexist(1:nsym)
               ! end  debug
            end do
         end do
      end do
      call phase_error_with_msg(nfout, 'Set of generators of the space group was not found.',__LINE__,__FILE__)
    end subroutine search_generators

    logical function group_is_the_same(fexist)
      logical, intent(inout) :: fexist(nsym)

      integer :: i,j
      integer :: n,no
      logical :: flag(nsym)

      n=0
      no=-1
      do while(n > no)
         flag(1:nsym) = .false.
         do j=1,nsym
            do i=1,nsym
               if(fexist(i).and.fexist(j)) flag(iptab(i,j)) = .true.
            end do
         end do
         do i=1,nsym
            if(flag(i)) fexist(i) = .true.
         end do
         no = n
         n = 0
         do i=1,nsym
            if(fexist(i)) n=n+1
         end do
      end do
      if(n==nsym) then
         group_is_the_same = .true.
      else
         group_is_the_same = .false.
      end if

    end function group_is_the_same

    subroutine set_generators
      integer :: i,j,k,n
      real(kind=DP) :: t,r,tb(3)
      real(kind=DP), parameter :: eps = 1.d-4
      ngen = ng
      call dealloc_igen_jgen
      call alloc_igen_jgen
      do i=1,ngen
         igen(i) = nopr(ig(i))
         do j=1,3
           tb(j) = sum(tau(1:3,ig(i))*b2pmat(1:3,j))
         end do
         do j=1,3
            t = tb(j)
            do k=1,20
               r=t*k
               if(abs(nint(r)-r)<eps) then
                  jgen(1,j,i) = nint(r)
                  jgen(2,j,i) = k
                  exit
               end if
            end do
         end do
      end do
      if(af>0) then
         iaf = nopr(nsym+af)
         do j=1,3
           tb(j) = sum(tau(1:3,nsym+af)*b2pmat(1:3,j))
         end do
         do j=1,3
            t = tb(j)
            do k=1,20
               r=t*k
               if(abs(nint(r)-r)<eps) then
                  jaf(1,j) = nint(r)
                  jaf(2,j) = k
                  exit
               end if
            end do
         end do
      end if

    end subroutine set_generators

    subroutine set_translations_zero
      integer :: i,j
      do i=1,ngen
         do j=1,3
            jgen(1,j,i) = 0
            jgen(2,j,i) = 1
         end do
      end do
    end subroutine set_translations_zero

  end subroutine set_tspace_generators

  subroutine write_tspace_generators(nfout,system,nprim)
    integer, intent(in) :: nfout
    character(len=9), intent(in) :: system
    integer, intent(in) :: nprim

    integer :: i

    write(nfout,'("TSPACE Generators:")')
    if(nprim>1) write(nfout,'("Translations are invalid.")')
!============================= Modified by K. Tagami =============
!    if(system == 'cubic') then
!      do i=1,ngen
!         write(nfout,'("igen,jgen(2,3)=",i2,"(",a5,")",3(1x,i1,"/",i1))') igen(i),oh(igen(i)),jgen(:,:,i)
!      end do
!      if(af>0) then
!         write(nfout,'("iaf ,jaf (2,3)=",i2,"(",a5,")",3(1x,i1,"/",i1))') iaf,oh(iaf),jaf(:,:)
!      end if
!    else
!      do i=1,ngen
!         write(nfout,'("igen,jgen(2,3)=",i2,"(",a5,")",3(1x,i1,"/",i1))') igen(i),d6h(igen(i)),jgen(:,:,i)
!      end do
!      if(af>0) then
!         write(nfout,'("iaf ,jaf (2,3)=",i2,"(",a5,")",3(1x,i1,"/",i1))') iaf,d6h(iaf),jaf(:,:)
!      end if
!    end if
    if(system == 'cubic') then
      do i=1,ngen
         write(nfout,'("igen,jgen(2,3)=",i2,"(",a5,")",3(3x,i2,"/",i2))') igen(i),oh(igen(i)),jgen(:,:,i)
      end do
      if(af>0) then
         write(nfout,'("iaf ,jaf (2,3)=",i2,"(",a5,")",3(3x,i2,"/",i2))') iaf,oh(iaf),jaf(:,:)
      end if
    else
      do i=1,ngen
         write(nfout,'("igen,jgen(2,3)=",i2,"(",a5,")",3(3x,i2,"/",i2))') igen(i),d6h(igen(i)),jgen(:,:,i)
      end do
      if(af>0) then
         write(nfout,'("iaf ,jaf (2,3)=",i2,"(",a5,")",3(3x,i2,"/",i2))') iaf,d6h(iaf),jaf(:,:)
      end if
    end if
! ==================================================
  end subroutine write_tspace_generators

  subroutine mod1(t)
    real(kind=DP), intent(inout) :: t(3)
    real(kind=DP), parameter :: eps = 1.d-6
    integer :: k
    t(1:3) = mod(t(1:3),1.d0)
    do k=1,3
       if(t(k) < -eps) t(k) = t(k) + 1.d0
       if(t(k) > 1.d0 - eps) t(k) = t(k) - 1.d0
    end do
  end subroutine mod1

  subroutine get_point_group_name(nfout,system,nsym,nopr,pg_name,pg_name_i)
    integer, intent(in) :: nfout
    character(len=9), intent(in) :: system
    integer, intent(in) :: nsym
    integer, intent(in) :: nopr(nsym)
    character(len=3), intent(out) :: pg_name
    character(len=5), intent(out) :: pg_name_i

    logical :: fexist(48)
    integer :: g,i

    fexist(1:48) = .false.
    do i=1,nsym
       fexist(nopr(i)) = .true.
    end do

    pg_name = 'C1 '
    pg_name_i = '1    '
    g = nsym
    if(system == 'cubic') then
       if(fexist(25)) then
       ! IE
          g = g/2
       end if
       if(g==1) then
          if(fexist(25)) then
          ! IE
             pg_name = 'Ci '
             pg_name_i = '-1   '
          else
             pg_name = 'C1 '
             pg_name_i = '1    '
          end if
       else if(g==2) then
          if(fexist(2).or.fexist(3).or.fexist(4).or. &
           & fexist(13).or.fexist(14).or.fexist(15).or. &
           & fexist(16).or.fexist(17).or.fexist(18)) then
          ! C2x or C2y or C2z or
          ! C2a or C2b or C2c or
          ! C2d or C2e or C2f
             if(fexist(25)) then
                ! IE
                pg_name = 'C2h'
                pg_name_i = '2/m  '
             else
                pg_name = 'C2 '
                pg_name_i = '2    '
             end if
          else if(fexist(26).or.fexist(27).or.fexist(28).or. &
                & fexist(37).or.fexist(38).or.fexist(39).or. &
                & fexist(40).or.fexist(41).or.fexist(42)) then
          ! IC2x or IC2y or IC2z or
          ! IC2a or IC2b or IC2c or
          ! IC2d or IC2e or IC2f
             pg_name = 'Cs '
             pg_name_i = 'm    '
          end if
       else if(g==3) then
          if(fexist(25)) then
             ! IE
             pg_name = 'S6 '
             pg_name_i = '-3   '
          else
             pg_name = 'C3 '
             pg_name_i = '3    '
          end if
       else if(g==4) then
          if((fexist(2).and.fexist(3).and.fexist(4)).or.&
           & (fexist(2).and.fexist(16).and.fexist(18)).or.&
           & (fexist(3).and.fexist(15).and.fexist(17)).or.&
           & (fexist(4).and.fexist(13).and.fexist(14))) then
          ! C2x and C2y and C2z or
          ! C2x and C2d and C2f or
          ! C2y and C2c and C2e or
          ! C2z and C2a and C2b
             if(fexist(25)) then
                ! IE
                pg_name = 'D2h'
                pg_name_i = 'mmm  '
             else
                pg_name = 'D2 '
                pg_name_i = '222  '
             end if
          else if(fexist(26).or.fexist(27).or.fexist(28).or. &
                & fexist(37).or.fexist(38).or.fexist(39).or. &
                & fexist(40).or.fexist(41).or.fexist(42).and.&
                & .not.(fexist(19).or.fexist(20).or.fexist(21))) then
                ! IC2x or IC2y or IC2z or
                ! IC2a or IC2b or IC2c or
                ! IC2d or IC2e or IC2f
             pg_name = 'C2v'
             pg_name_i = 'mm2  '
          else if((fexist(19).and.fexist(22)).or. &
                & (fexist(20).and.fexist(23)).or. &
                & (fexist(21).and.fexist(24))) then
                ! (C4x+ and C4x-) or (C4y+ and C4y-) or (C4z+ and C4z-)
             if(fexist(25)) then
                ! IE
                pg_name = 'C4h'
                pg_name_i = '4/m'
             else
                pg_name = 'C4 '
                pg_name_i = '4    '
             end if
          else if((fexist(43).and.fexist(46)).or. &
                & (fexist(44).and.fexist(47)).or. &
                & (fexist(45).and.fexist(48))) then
                ! (IC4x+ and IC4x-) or (IC4y+ and IC4y-) or (IC4z+ and IC4z-)
             pg_name = 'S4 '
             pg_name_i = '-4   '
          end if
       else if(g==6) then
          if((fexist(5).and.fexist(14).and.fexist(17).and.fexist(18)) .or. &
           & (fexist(6).and.fexist(14).and.fexist(15).and.fexist(16)) .or. &
           & (fexist(7).and.fexist(13).and.fexist(15).and.fexist(18)) .or. &
           & (fexist(8).and.fexist(13).and.fexist(16).and.fexist(17))) then
           ! C31+ and C2b and C2e and C2f
           ! C32+ and C2b and C2c and C2d
           ! C33+ and C2a and C2c and C2f
           ! C34+ and C2a and C2d and C2e
             if(fexist(25)) then
             ! IE
                pg_name = 'D3d'
                pg_name_i = '-3m  '
             else
                pg_name = 'D3 '
                pg_name_i = '32   '
             end if
          else
             pg_name = 'C3v'
             pg_name_i = '3m   '
          end if
       else if(g==8) then
          if((fexist(19).and.fexist(22)).or. &
           & (fexist(20).and.fexist(23)).or. &
           & (fexist(21).and.fexist(24))) then
           ! (C4x+ and C4x-) or (C4y+ and C4y-) or (C4z+ and C4z-)
             if(fexist(37).or.fexist(38).or.fexist(39).or. &
              & fexist(40).or.fexist(41).or.fexist(42)) then
              ! IC2a or IC2b or IC2c or IC2d or IC2e or IC2f
                if(fexist(25)) then
                ! IE
                   pg_name = 'D4h'
                   pg_name_i = '4/mmm'
                else
                   pg_name = 'C4v'
                   pg_name_i = '4mm  '
                end if
             else
                pg_name = 'D4 '
                pg_name_i = '422  '
             end if
          else
             pg_name = 'D2d'
             pg_name_i = '-42m '
          end if
       else if(g==12) then
          if(fexist(25)) then
          ! IE
             pg_name = 'Th '
             pg_name_i = 'm-3   '
          else
             pg_name = 'T  '
             pg_name_i = '23   '
          end if
       else if(g==24) then
          if(fexist(19).and.fexist(22).and. &
           & fexist(20).and.fexist(23).and. &
           & fexist(21).and.fexist(24)) then
           ! C4x+ and C4x- and C4y+ and C4y- and C4z+ and C4z-
             if(fexist(25)) then
             ! IE
                pg_name = 'Oh '
                pg_name_i = 'm-3m  '
             else
                pg_name = 'O  '
                pg_name_i = '432  '
             end if
          else
             pg_name = 'Td '
             pg_name_i = '-43m '
          end if
       end if
    else ! hexagonal
       if(fexist(13)) then
       ! IE
          g = g/2
       else if(fexist(16)) then
       ! IC2
          g = g/2
       end if
       if(g==2) then
           if(fexist(7).or.fexist(8).or.fexist(9).or. &
           & fexist(10).or.fexist(11).or.fexist(12).or.fexist(4)) then
           ! C211 or C221 or C231 or C212 or C222 or C232 or C2
             if(fexist(13)) then
                ! IE
                pg_name = 'C2h'
                pg_name_i = '2/m  '
             else
                pg_name = 'C2 '
                pg_name_i = '2    '
             end if
          else if(fexist(19).or.fexist(20).or.fexist(21).or. &
                & fexist(22).or.fexist(23).or.fexist(24).or.fexist(16)) then
                 ! IC211 or IC221 or IC231 or IC212 or IC222 or IC232 or IC2
             pg_name = 'Cs '
             pg_name_i = 'm    '
          end if
       else if(g==3) then
          if(fexist(13)) then
          ! IE
             pg_name = 'S6 '
             pg_name_i = '-3   '
          else if(fexist(16)) then
          ! IC2
             pg_name = 'C3h'
             pg_name_i = '-6   '
          else
             pg_name = 'C3 '
             pg_name_i = '3    '
          end if
       else if(g==4) then
          if(fexist(4).and.( &
            & (fexist(7).and.fexist(10)).or. &
            & (fexist(8).and.fexist(11)).or. &
            & (fexist(9).and.fexist(12)))) then
          ! C2 and ( C211 and C212 ) or ( C221 and C222 ) or ( C231 and C232 )
             if(fexist(13)) then
                ! IE
                pg_name = 'D2h'
                pg_name_i = 'mmm  '
             else
                pg_name = 'D2 '
                pg_name_i = '222  '
             end if
          else if( (fexist(4).and.( &
            & (fexist(19).and.fexist(22)).or. &
            & (fexist(20).and.fexist(23)).or. &
            & (fexist(21).and.fexist(24))) ) .or. &
            & (fexist(16).and.( &
            & (fexist(7).and.fexist(22)).or. &
            & (fexist(8).and.fexist(23)).or. &
            & (fexist(9).and.fexist(24)).or. &
            & (fexist(19).and.fexist(10)).or. &
            & (fexist(20).and.fexist(11)).or. &
            & (fexist(21).and.fexist(12))) ) ) then
          ! C2 and (IC211 and IC212) or (IC221 and IC222) or (IC231 and IC232)
          ! or IC2 and (C211 and IC212) or (C221 and IC222) or (C231 and IC232)
          !         or (IC211 and C212) or (IC221 and C222) or (IC231 and C232)
             pg_name = 'C2v'
             pg_name_i = 'mm2  '
          end if
       else if(g==6) then
          if(fexist(2).and.fexist(6)) then
          ! C6+ and C6-
             if(fexist(13)) then
             ! IE
                pg_name = 'C6h'
                pg_name_i = '6/m  '
             else
                pg_name = 'C6 '
                pg_name_i = '6    '
             end if
          else if((fexist(7).and.fexist(8).and.fexist(9)) .or. &
                & (fexist(10).and.fexist(11).and.fexist(12))) then
                 ! (C211 and C221 and C231) or (C212 and C222 and C232)
             if(fexist(13)) then
             ! IE
                pg_name = 'D3d'
                pg_name_i = '-3m  '
             else if(fexist(16)) then
             ! IC2
                pg_name = 'D3h'
                pg_name_i = '-6m2 '
             else
                pg_name = 'D3 '
                pg_name_i = '32   '
             end if
          else if((fexist(19).and.fexist(20).and.fexist(21)) .or. &
                & (fexist(22).and.fexist(23).and.fexist(24))) then
                 ! (IC211 and IC221 and IC231) or (IC212 and IC222 and IC232)
             pg_name = 'C3v'
             pg_name_i = '3m   '
          end if
       else if(g==12) then
          if(fexist(7).and.fexist(8).and.fexist(9).and. &
           & fexist(10).and.fexist(11).and.fexist(12)) then
          ! C211 and C221 and C231 and C212 and C222 and C232
             if(fexist(13)) then
             ! IE
                pg_name = 'D6h'
                pg_name_i = '6/mmm'
             else
                pg_name = 'D6 '
                pg_name_i = '622  '
             end if
          else if(fexist(19).and.fexist(20).and.fexist(21).and. &
                & fexist(22).and.fexist(23).and.fexist(24)) then
                 ! IC211 and IC221 and IC231 and IC212 and IC222 and IC232
             pg_name = 'C6v'
             pg_name_i = '6mm  '
          end if
       end if
    end if
    if(g>1.and.pg_name=='C1 ') then
      write(nfout,*) 'get_point_group_name: error'
      write(nfout,*) 'system=',system
      write(nfout,*) 'nsym=',nsym
      write(nfout,*) 'g=',g
      write(nfout,*) 'fexist=',fexist
      call phase_error_with_msg(nfout, 'get_point_group_name: error',__LINE__,__FILE__)
    end if
  end subroutine get_point_group_name

  subroutine get_bravais_lattice(nfout,ipri,a1,a2,a3,system,nopr,pg_name_i,bl_name)
    integer, intent(in) :: nfout,ipri
    real(kind=DP), intent(in) :: a1(3),a2(3),a3(3)
    character(len=9), intent(in) :: system
    integer, intent(in) :: nopr(48)
    character(len=5), intent(in) :: pg_name_i
    character(len=2), intent(out) :: bl_name

    real(kind=DP) :: length(3),angle(3)
    real(kind=DP), parameter :: eps = 1.d-3
    real(kind=DP) :: a(3,3),b(3,3)

    a(1:3,1) = a1(1:3)
    a(1:3,2) = a2(1:3)
    a(1:3,3) = a3(1:3)

    call vectors_length_angle(a1,a2,a3,length,angle)
    angle(1:3) = angle(1:3)/PAI*180.d0

    bl_name = 'No'
    select case(pg_name_i)
    case('-1','1') ! Trigonal
       bl_name = 'aP'
    case('2/m','m','2') ! Monoclinic
       if(cell_is_monoclinic(a,length,angle)) then
          bl_name = 'mP'
          go to 100
       end if

       call change_primitive_to_bravais('I',a,b)
       if(cell_is_monoclinic(b,length,angle)) then
          bl_name = 'mI'
          go to 100
       end if

! --- ad hoc ----
       goto 1000

       if(nopr(2)==2 .or. nopr(2)==26) then
          call change_primitive_to_bravais('A',a,b)
          if(cell_is_monoclinic(b,length,angle)) then
             bl_name = 'mA'
             go to 100
          end if

       else if(nopr(2)==3 .or. nopr(2)==27) then
          call change_primitive_to_bravais('B',a,b)
          if(cell_is_monoclinic(b,length,angle)) then
             bl_name = 'mB'
             go to 100
          end if

       else if(nopr(2)==4 .or. nopr(2)==28) then
          call change_primitive_to_bravais('C',a,b)

          if(cell_is_monoclinic(b,length,angle)) then
             bl_name = 'mC'
             go to 100
          end if
       end if
! -----------------------

1000   continue
       call change_primitive_to_bravais('A',a,b)
       if(cell_is_monoclinic(b,length,angle)) then
          bl_name = 'mA'
          go to 100
       end if

       call change_primitive_to_bravais('B',a,b)
       if(cell_is_monoclinic(b,length,angle)) then
          bl_name = 'mB'
          go to 100
       end if
       call change_primitive_to_bravais('C',a,b)
       if(cell_is_monoclinic(b,length,angle)) then
          bl_name = 'mC'
          go to 100
       end if
! ----

       ! Hexagonal or Rhombohedral
       bl_name = 'mP'
       call vectors_length_angle(a1,a2,a3,length,angle)
       angle(1:3) = angle(1:3)/PAI*180.d0
       if(cell_is_rhombohedral(length,angle)) then
          length(1) = sqrt(sum((a2-a1)**2))
          length(2) = sqrt(sum((a3-a2)**2))
          length(3) = sqrt(sum((a1+a2+a3)**2))
       end if
       angle(1:2) = 90.d0
       angle(3) = 120.d0

    case('mmm','mm2','222') ! Orthorhombic
       if(cell_is_orthorhombic(a,length,angle)) then
          bl_name = 'oP'
          go to 100
       end if
       call change_primitive_to_bravais('A',a,b)
       if(cell_is_orthorhombic(b,length,angle)) then
          bl_name = 'oA'
          go to 100
       end if
       call change_primitive_to_bravais('B',a,b)
       if(cell_is_orthorhombic(b,length,angle)) then
          bl_name = 'oB'
          go to 100
       end if
       call change_primitive_to_bravais('C',a,b)
       if(cell_is_orthorhombic(b,length,angle)) then
          bl_name = 'oC'
          go to 100
       end if
       call change_primitive_to_bravais('F',a,b)
       if(cell_is_orthorhombic(b,length,angle)) then
          bl_name = 'oF'
          go to 100
       end if
       call change_primitive_to_bravais('I',a,b)
       if(cell_is_orthorhombic(b,length,angle)) then
          bl_name = 'oI'
          go to 100
       end if
       ! Hexagonal or Rhombohedral
       bl_name = 'oP'
       call vectors_length_angle(a1,a2,a3,length,angle)
       angle(1:3) = angle(1:3)/PAI*180.d0
       if(cell_is_rhombohedral(length,angle)) then
          length(1) = sqrt(sum((a2-a1)**2))
          length(2) = sqrt(sum((a3-a2)**2))
          length(3) = sqrt(sum((a1+a2+a3)**2))
       end if
       length(2) = 0.5d0*sqrt(3.d0)*length(2)
       angle(1:3) = 90.d0
    case('4/mmm','-42m','4mm','422','4/m','-4','4') ! Tetragonal
       if(cell_is_tetragonal(a,length,angle)) then
          bl_name = 'tP'
          go to 100
       end if
       call change_primitive_to_bravais('I',a,b)
       if(cell_is_tetragonal(b,length,angle)) then
          bl_name = 'tI'
          go to 100
       end if
    case('-3m','3m','32','-3','3') ! Trigonal
       bl_name = 'hP'
       !!$if(abs(length(1)-length(2)) < eps .and. &
       !!$ & abs(length(2)-length(3)) < eps .and. &
       !!$ & abs(angle(1)-angle(2)) < eps .and. &
       !!$ & abs(angle(2)-angle(3)) < eps ) then ! cell is Rhombohedral
       if(cell_is_rhombohedral(length,angle)) then
          bl_name = 'hR'
       end if
    case('6/mmm','-6m2','6mm','622','6/m','-6','6') ! Hexagonal
       bl_name = 'hP'
    case('m-3m','-43m','432','m-3','23') ! Cubic
       if(cell_is_cubic(a,length,angle)) then
          bl_name = 'cP'
          go to 100
       end if
       call change_primitive_to_bravais('F',a,b)
       if(cell_is_cubic(b,length,angle)) then
          bl_name = 'cF'
          go to 100
       end if
       call change_primitive_to_bravais('I',a,b)
       if(cell_is_cubic(b,length,angle)) then
          bl_name = 'cI'
          go to 100
       end if
    case default
       bl_name = 'No'
    end select
100 continue
    if(ipri > 0) then
       write(nfout,'("=== Lattice parameters ===")')
       write(nfout,'("a    ,b   ,c     = ",3f15.8," Bohr")') length
       write(nfout,'("alpha,beta,gamma = ",3f15.8," degree")') angle
    end if
    call m_IS_put_latconst_len_angle(length,angle)
    call m_CS_put_Bravais_lattice(length,angle)
!!$    latconst_len = length
!!$    latconst_angle = angle
  contains
    subroutine change_primitive_to_bravais(l,a,b)
      character(len=1), intent(in) :: l
      real(kind=DP), intent(in)    :: a(3,3) ! primitive
      real(kind=DP), intent(out)   :: b(3,3) ! bravais

      integer :: i,j
      real(kind=DP) :: p2b(3,3),tt(3)

      p2b = 0.d0
      if(l=='R' .or. l=='P') then
         ! Hexagonal or primitive
         do i=1,3
            p2b(i,i)=1.d0
         end do
      else if(l=='A') then
         ! A-face centered
         p2b(1,1) =  1.d0
         p2b(2,2) =  1.d0
         p2b(3,3) =  1.d0
         p2b(3,2) = -1.d0
         p2b(2,3) =  1.d0
      else if(l=='B') then
         ! B-face centered
         p2b(2,2) =  1.d0
         p2b(1,1) =  1.d0
         p2b(3,3) =  1.d0
         p2b(3,1) = -1.d0
         p2b(1,3) =  1.d0
      else if(l=='C') then
         ! C-face centered
         p2b(3,3) =  1.d0
         p2b(1,1) =  1.d0
         p2b(2,2) =  1.d0
         p2b(2,1) = -1.d0
         p2b(1,2) =  1.d0
      else if(l=='F') then
         ! All-face centered
         p2b(1:3,1:3) =  1.d0
         do i=1,3
            p2b(i,i)=-1.d0
         end do
      else if(l=='I') then
         ! Body centered
         p2b(1:3,1:3) =  1.d0
         do i=1,3
            p2b(i,i)=0.d0
         end do
      end if

      do i=1,3
         b(1:3,i) = 0.d0
         do j=1,3
            b(1:3,i) = b(1:3,i) + p2b(i,j)*a(1:3,j)
         end do
      end do
    end subroutine change_primitive_to_bravais

    logical function cell_is_cubic(b,length,angle)
      real(kind=DP), intent(in) :: b(3,3)
      real(kind=DP), intent(out) :: length(3),angle(3)
      cell_is_cubic = .false.
      call vectors_length_angle(b(1,1),b(1,2),b(1,3),length,angle)
      angle(1:3) = angle(1:3)/PAI*180.d0
      if( abs(length(1)-length(2)) < eps .and. &
        & abs(length(2)-length(3)) < eps .and. &
        & abs(angle(1) - 90.d0) < eps .and. &
        & abs(angle(2) - 90.d0) < eps .and. &
        & abs(angle(3) - 90.d0) < eps) then
        cell_is_cubic = .true.
      end if
    end function cell_is_cubic

    logical function cell_is_tetragonal(b,length,angle)
      real(kind=DP), intent(in) :: b(3,3)
      real(kind=DP), intent(out) :: length(3),angle(3)
      cell_is_tetragonal = .false.
      call vectors_length_angle(b(1,1),b(1,2),b(1,3),length,angle)
      angle(1:3) = angle(1:3)/PAI*180.d0
      if((abs(length(1)-length(2)) < eps .or. &
        & abs(length(2)-length(3)) < eps .or. &
        & abs(length(3)-length(1)) < eps ) .and. &
        & abs(angle(1) - 90.d0) < eps .and. &
        & abs(angle(2) - 90.d0) < eps .and. &
        & abs(angle(3) - 90.d0) < eps) then
        cell_is_tetragonal = .true.
      end if
    end function cell_is_tetragonal

    logical function cell_is_orthorhombic(b,length,angle)
      real(kind=DP), intent(in) :: b(3,3)
      real(kind=DP), intent(out) :: length(3),angle(3)
      cell_is_orthorhombic = .false.
      call vectors_length_angle(b(1,1),b(1,2),b(1,3),length,angle)
      angle(1:3) = angle(1:3)/PAI*180.d0
      !!$if(abs(angle(1) - 90.d0) < eps .and. &
      !!$  & abs(angle(2) - 90.d0) < eps .and. &
      !!$  & abs(angle(3) - 90.d0) < eps ) then
      if(abs(b(1,1))>eps.and.abs(b(2,1))<eps.and.abs(b(3,1))<eps.and.&
       & abs(b(1,2))<eps.and.abs(b(2,2))>eps.and.abs(b(3,2))<eps.and.&
       & abs(b(1,3))<eps.and.abs(b(2,3))<eps.and.abs(b(3,3))>eps) then
        cell_is_orthorhombic = .true.
      end if
    end function cell_is_orthorhombic

    logical function cell_is_monoclinic(b,length,angle)
      real(kind=DP), intent(in) :: b(3,3)
      real(kind=DP), intent(out) :: length(3),angle(3)
      cell_is_monoclinic = .false.
      call vectors_length_angle(b(1,1),b(1,2),b(1,3),length,angle)
      angle(1:3) = angle(1:3)/PAI*180.d0
      if((abs(angle(1) - 90.d0) < eps .and. &
        & abs(angle(1) - angle(2)) < eps) .or. &
        &(abs(angle(2) - 90.d0) < eps .and. &
        & abs(angle(2) - angle(3)) < eps) .or. &
        &(abs(angle(3) - 90.d0) < eps .and. &
        & abs(angle(3) - angle(1)) < eps) ) then
        cell_is_monoclinic = .true.
      end if
    end function cell_is_monoclinic

    logical function cell_is_rhombohedral(length,angle)
      real(kind=DP), intent(in) :: length(3),angle(3)
      cell_is_rhombohedral = .false.
      if(abs(length(1)-length(2)) < eps .and. &
        & abs(length(2)-length(3)) < eps .and. &
        & abs(angle(1)-angle(2)) < eps .and. &
        & abs(angle(2)-angle(3)) < eps ) then ! cell is Rhombohedral
          cell_is_rhombohedral = .true.
       end if
    end function cell_is_rhombohedral

  end subroutine get_bravais_lattice


  subroutine get_bravais_lattice_old(nfout,a1,a2,a3,system,nopr,pg_name_i,bl_name)
    integer, intent(in) :: nfout
    real(kind=DP), intent(in) :: a1(3),a2(3),a3(3)
    character(len=9), intent(in) :: system
    integer, intent(in) :: nopr(48)
    character(len=5), intent(in) :: pg_name_i
    character(len=2), intent(out) :: bl_name

    real(kind=DP) :: length(3),angle(3)
    real(kind=DP), parameter :: eps = 1.d-3

    length(1) = sqrt(sum(a1(1:3)*a1(1:3)))
    length(2) = sqrt(sum(a2(1:3)*a2(1:3)))
    length(3) = sqrt(sum(a3(1:3)*a3(1:3)))
    angle(1) = acos(sum(a2(1:3)*a3(1:3))/(length(2)*length(3)))
    angle(2) = acos(sum(a3(1:3)*a1(1:3))/(length(3)*length(1)))
    angle(3) = acos(sum(a1(1:3)*a2(1:3))/(length(1)*length(2)))

    angle(1:3) = angle(1:3)*180/PAI

    select case(pg_name_i)
    case('-1','1') ! Trigonal
       bl_name = 'aP'
    case('2/m','m','2') ! Monoclinic
       bl_name = 'mP'
       if(nopr(2) == 2) then
       ! C2x -> a-unique setting
          if(abs(angle(3)-90.d0) > eps) then
             if(abs(angle(2)-90.d0) > eps) then
                bl_name = 'mI'
             else
                bl_name = 'mB'
             end if
          else if(abs(angle(1)-90.d0) > eps) then
             bl_name = 'mA'
          end if
       else if(nopr(2) == 3) then
       ! C2y -> b-unique setting
          if(abs(angle(1)-90.d0) > eps) then
             if(abs(angle(3)-90.d0) > eps) then
                bl_name = 'mI'
             else
                bl_name = 'mA'
             end if
          else if(abs(angle(3)-90.d0) > eps) then
             bl_name = 'mC'
          end if
       else if(nopr(2) == 4) then
       ! C2z -> c-unique setting
          if(abs(angle(2)-90.d0) > eps) then
             if(abs(angle(1)-90.d0) > eps) then
                bl_name = 'mI'
             else
                bl_name = 'mC'
             end if
          else if(abs(angle(2)-90.d0) > eps) then
             bl_name = 'mB'
          end if
       end if
    case('mmm','mm2','222') ! Orthorhombic
       bl_name = 'oP'
       if(abs(angle(1)-90.d0) > eps .and. &
        & abs(angle(2)-90.d0) < eps .and. &
        & abs(angle(3)-90.d0) < eps ) then
          if(abs(a1(1))>eps) then
             bl_name = 'oA'
          else if(abs(a1(2))>eps) then
             bl_name = 'oB'
          else if(abs(a1(3))>eps) then
             bl_name = 'oC'
          end if
       else if(abs(angle(2)-90.d0) > eps .and. &
        & abs(angle(3)-90.d0) < eps .and. &
        & abs(angle(1)-90.d0) < eps ) then
          if(abs(a2(1))>eps) then
             bl_name = 'oA'
          else if(abs(a2(2))>eps) then
             bl_name = 'oB'
          else if(abs(a2(3))>eps) then
             bl_name = 'oC'
          end if
       else if(abs(angle(3)-90.d0) > eps .and. &
          & abs(angle(1)-90.d0) < eps .and. &
          & abs(angle(2)-90.d0) < eps ) then
          if(abs(a3(1))>eps) then
             bl_name = 'oA'
          else if(abs(a3(2))>eps) then
             bl_name = 'oB'
          else if(abs(a3(3))>eps) then
             bl_name = 'oC'
          end if
       else if((abs(a1(1)) < eps .and. &
          &  abs(a1(2)) > eps .and. &
          &  abs(a1(3)) > eps) .or. &
          & (abs(a1(2)) < eps .and. &
          &  abs(a1(3)) > eps .and. &
          &  abs(a1(1)) > eps) .or. &
          & (abs(a1(3)) < eps .and. &
          &  abs(a1(1)) > eps .and. &
          &  abs(a1(2)) > eps) ) then
          bl_name = 'oF'
       else if(abs(length(1)-length(2)) < eps .and. &
          & abs(length(2)-length(3)) < eps .and. &
          & ( abs(angle(1)-90.d0) > eps .or. &
          &   abs(angle(2)-90.d0) > eps .or. &
          &   abs(angle(3)-90.d0) > eps ) ) then
          bl_name = 'oI'
       end if
    case('4/mmm','-42m','4mm','422','4/m','-4','4') ! Tetragonal
       bl_name = 'tP'
       if(abs(angle(1)-angle(2)) > eps .or. &
        & abs(angle(2)-angle(3)) > eps ) then
          bl_name = 'tI'
       end if
    case('-3m','3m','32','-3','3') ! Trigonal
       bl_name = 'hP'
       if(abs(length(1)-length(2)) < eps .and. &
        & abs(length(2)-length(3)) < eps .and. &
        & abs(angle(1)-angle(2)) < eps .and. &
        & abs(angle(2)-angle(3)) < eps ) then
          bl_name = 'hR'
       end if
    case('6/mmm','-6m2','6mm','622','6/m','-6','6') ! Hexagonal
       bl_name = 'hP'
    case('m-3m','-43m','432','m-3','23') ! Cubic
       bl_name = 'cP'
       if(abs(angle(1) - 60.d0) < eps) then
          bl_name = 'cF'
       else if(abs(angle(1) - 109.47122d0) < eps) then
          bl_name = 'cI'
       end if
    case default
       bl_name = 'No'
    end select

  end subroutine get_bravais_lattice_old

  subroutine get_space_group_name(nfout,system,nsym,nopr,tau,fbravais,bl_name,pg_name_i,sg_name_i)
    integer, intent(in) :: nfout
    character(len=9), intent(in) :: system
    integer, intent(in) :: nsym
    integer, intent(in) :: nopr(48)
    real(kind=DP), intent(in) :: tau(3,*)
    logical, intent(in) :: fbravais
    character(len=2), intent(in) :: bl_name
    character(len=5), intent(in) :: pg_name_i
    character(len=10), intent(out) :: sg_name_i

    character(len=1) :: s,l,sp
    character(len=2) :: sub
    character(len=9) :: sg_name0
    integer :: i,j
    logical :: fexist(48)
    real(kind=DP) :: t(3,48)
    real(kind=DP), parameter :: eps = 1.d-6

    fexist(1:48) = .false.
    t(1:3,1:48) = 0.d0
    do i=1,nsym
       fexist(nopr(i)) = .true.
       t(1:3,nopr(i)) = tau(1:3,i)
    end do

    sp=''
    s(1:1) = bl_name(1:1)
    l(1:1) = bl_name(2:2)

    if(.not.fbravais) call get_trans_in_bravais_system(l,t)

    select case(pg_name_i)
    case('1','-1')     ! Triclinic
       sg_name_i = l // trim(pg_name_i)
    case('2')     ! Monoclinic
       sg_name0= trim(pg_name_i)
       if(system == 'cubic') then
          if(fexist(2)) then
          ! C2x
             if(abs(t(1,2)) > eps) sg_name0= trim(pg_name_i) // '_1'
          else if(fexist(3)) then
          ! C2y
             if(abs(t(2,3)) > eps) sg_name0= trim(pg_name_i) // '_1'
          else if(fexist(4)) then
          ! C2z
             if(abs(t(3,4)) > eps) sg_name0= trim(pg_name_i) // '_1'
          end if
       else ! hexagonal
          if(fexist(4)) then
             ! C2
             if(sum(abs(t(1:3,i))) > eps) sg_name0= trim(pg_name_i) // '_1'
          else
             do i=7,12
                if(fexist(i)) then
             ! C211 or C221 or C231 or C212 or C222 or C232
                   if(sum(abs(t(1:3,i))) > eps) sg_name0= trim(pg_name_i) // '_1'
                end if
             end do
          end if
       end if
       sg_name_i = l // trim(sg_name0)
    case('m')     ! Monoclinic
       if(system == 'cubic') then
          do i=0,2
             if(fexist(26+i)) then
             ! IC2x or IC2y or IC2z
               call get_symbol_of_sym_plane(t(1,26+i),sp)
             end if
          end do
       else ! hexagonal
          if(fexist(16)) then
             ! IC2
             call get_symbol_of_sym_plane(t(1,i),sp)
          else
             do i=19,24
                if(fexist(i)) then
             ! IC211 or IC221 or IC231 or IC212 or IC222 or IC232
                   call get_symbol_of_sym_plane(t(1,i),sp)
                end if
             end do
          end if
       end if
       if(sp=='a' .or. sp=='b' .or. sp=='n') then
          sp = 'c'
       end if
       sg_name_i = l // trim(sp)
    case('2/m')   ! Monoclinic
       sg_name0 = '2'
       if(system == 'cubic') then
          do i=0,2
             if(fexist(2+i)) then
             ! C2x or C2y or C2z
                if(sum(abs(t(1:3,2+i))) > eps) sg_name0= trim(sg_name0) // '_1'
                ! IC2x or IC2y or IC2z
                call get_symbol_of_sym_plane(t(1,26+i),sp)
             end if
          end do
          do i=0,5
             if(fexist(13+i)) then
             ! C2a or C2b or C2c or C2d or C2e or C2f
                if(sum(abs(t(1:3,13+i))) > eps) sg_name0= trim(sg_name0) // '_1'
                ! IC2a or IC2b or IC2c or IC2d or IC2e or IC2f
                call get_symbol_of_sym_plane(t(1,37+i),sp)
             end if
          end do
       else ! hexagonal
          if(fexist(4)) then
             ! C2
             if(sum(abs(t(1:3,4))) > eps) sg_name0= trim(pg_name_i) // '_1'
             ! IC2
             call get_symbol_of_sym_plane(t(1,16),sp)
          else
             do i=7,12
                if(fexist(i)) then
                ! C211 or C221 or C231 or C212 or C222 or C232
                   if(sum(abs(t(1:3,i))) > eps) sg_name0= trim(pg_name_i) // '_1'
                ! IC211 or IC221 or IC231 or IC212 or IC222 or IC232
                   call get_symbol_of_sym_plane(t(1,i+12),sp)
                end if
             end do
          end if
       end if
       if(sp=='a' .or. sp=='b' .or. sp=='n') sp = 'c'
       sg_name0= trim(sg_name0) // '/' // sp
       sg_name_i = l // trim(sg_name0)
    case('222')   ! Orthorhombic
       if(system == 'cubic') then
          ! C2x
          sg_name0 = '2'
          if(abs(t(1,2)) > eps) sg_name0= trim(sg_name0) // '_1'
          ! C2y
          sg_name0 = trim(sg_name0) // '2'
          if(abs(t(2,3)) > eps) sg_name0= trim(sg_name0) // '_1'
          ! C2z
          sg_name0 = trim(sg_name0) // '2'
          if(abs(t(3,4)) > eps) sg_name0= trim(sg_name0) // '_1'
          sg_name_i = l // trim(sg_name0)
       else ! hexagonal
          do i=0,2
             if(fexist(10+i)) then
                ! C212 or C222 or C232 (C2x)
                sg_name0 = '2'
                if(sum(abs(t(1:3,10+i))) > eps) sg_name0= trim(sg_name0) // '_1'
                ! C211 or C221 or C231 (C2y)
                sg_name0 = trim(sg_name0) // '2'
                if(sum(abs(t(1:3,7+i))) > eps) sg_name0= trim(sg_name0) // '_1'
                exit
             end if
          end do
          ! C2 (C2z)
          sg_name0 = trim(sg_name0) // '2'
          if(sum(abs(t(1:3,4))) > eps) sg_name0= trim(sg_name0) // '_1'
          sg_name_i = l // trim(sg_name0)
       end if
    case('mm2')   ! Orthorhombic
       sg_name0 = ''
       if(system == 'cubic') then
          do i=0,2
             if(fexist(2+i)) then
             ! C2x or C2y or C2z
                j = mod(i+1,3)
                if(fexist(26+j)) then
                ! IC2y or IC2z or IC2x
                   call get_symbol_of_sym_plane(t(1,26+j),sp)
                   sg_name0= trim(sg_name0) // sp
                end if
                j = mod(i+2,3)
                if(fexist(26+j)) then
                ! IC2z or IC2x or IC2y
                   call get_symbol_of_sym_plane(t(1,26+j),sp)
                   sg_name0= trim(sg_name0) // sp
                end if
                if(fexist(40-i)) then
                ! IC2d or IC2c or IC2b
                   call get_symbol_of_sym_plane(t(1,40-i),sp)
                   sg_name0= trim(sg_name0) // sp
                   if(40-i==38.and.fexist(37)) then
                   ! IC2f
                      call get_symbol_of_sym_plane(t(1,37),sp)
                   else if(40-i==39.and.fexist(41)) then
                   ! IC2e
                      call get_symbol_of_sym_plane(t(1,41),sp)
                   else if(40-i==40.and.fexist(42)) then
                   ! IC2a
                      call get_symbol_of_sym_plane(t(1,42),sp)
                   end if
                   sg_name0= trim(sg_name0) // sp
                end if
                sg_name0 = trim(sg_name0) // '2'
                if(abs(t(1+i,2+i)) > eps) sg_name0= trim(sg_name0) // '_1'
             end if
          end do
       else ! hexagonal
          if(fexist(4)) then
          ! c2 (C2z)
             do i=0,2
                if(fexist(22+i)) then
                   ! IC212 or IC222 or IC232 (IC2x)
                   call get_symbol_of_sym_plane(t(1,22+i),sp)
                   sg_name0= trim(sg_name0) // sp
                   ! IC211 or IC221 or IC231 (IC2y)
                   call get_symbol_of_sym_plane(t(1,19+i),sp)
                   sg_name0= trim(sg_name0) // sp
                end if
             end do
             sg_name0 = trim(sg_name0) // '2'
             if(sum(abs(t(1:3,4))) > eps) sg_name0= trim(sg_name0) // '_1'
          else
             do i=0,5
                if(fexist(7+i)) then
                ! C211 or C221 or C231 (C2y)
                ! or C212 or C222 or C232 (C2x)
                   if(i<3) then
                      j=i+3+19
                   else
                      j=i-3+19
                   end if
                   ! IC212 or IC222 or IC232 (IC2x)
                   ! or IC211 or IC221 or IC231 (IC2y)
                   call get_symbol_of_sym_plane(t(1,j),sp)
                   sg_name0= trim(sg_name0) // sp
                   ! IC2 (IC2z)
                   call get_symbol_of_sym_plane(t(1,16),sp)
                   sg_name0= trim(sg_name0) // sp

                   sg_name0 = trim(sg_name0) // '2'
                   if(sum(abs(t(1:3,7+i))) > eps) &
                      & sg_name0= trim(sg_name0) // '_1'
                   exit
                end if
             end do
          end if
       end if
       sg_name_i = l // trim(sg_name0)
    case('mmm')   ! Orthorhombic
       sg_name0 = ''
       if(system == 'cubic') then
          do i=0,2
             ! IC2x or IC2y or IC2z
             call get_symbol_of_sym_plane(t(1,26+i),sp)
             sg_name0= trim(sg_name0) // sp
          end do
       else ! hexagonal
          do i=0,2
             if(fexist(22+i)) then
                ! IC212 or IC222 or IC232 (IC2x)
                call get_symbol_of_sym_plane(t(1,22+i),sp)
                sg_name0= trim(sg_name0) // sp
                ! IC211 or IC221 or IC231 (IC2y)
                call get_symbol_of_sym_plane(t(1,19+i),sp)
                sg_name0= trim(sg_name0) // sp
                exit
             end if
          end do
          ! IC2 (IC2z)
          call get_symbol_of_sym_plane(t(1,16),sp)
          sg_name0= trim(sg_name0) // sp
       end if
       sg_name_i = l // trim(sg_name0)
    case('4')     ! Tetragonal
       do i=0,2
         if(fexist(19+i)) then
          ! C4x+ or C4y+ or C4z+
            call get_screw_sub(t(1+i,19+i),4,sub)
            sg_name0= '4' // trim(sub)
            exit
         end if
       end do
       sg_name_i = l // trim(sg_name0)
    case('-4')    ! Tetragonal
       sg_name_i = l // trim(pg_name_i)
    case('4/m')   ! Tetragonal
       do i=0,2
         if(fexist(19+i)) then
          ! C4x+ or C4y+ or C4z+
            call get_screw_sub(t(1+i,19+i),4,sub)
            sg_name0= '4' // trim(sub)
          ! IC2x or IC2y or IC2z
            call get_symbol_of_sym_plane(t(1,26+i),sp)
            sg_name0= trim(sg_name0) // "/" // sp
            exit
         end if
       end do
       sg_name_i = l // trim(sg_name0)
    case('422')   ! Tetragonal
       do i=0,2
          if(fexist(19+i)) then
           ! C4x+ or C4y+ or C4z+
             call get_screw_sub(t(1+i,19+i),4,sub)
             sg_name0= '4' // trim(sub)
           ! C2x or C2y or C2z
             j = mod(i+1,3)
             call get_screw_sub(t(1+j,2+j),2,sub)
             sg_name0= trim(sg_name0) // '2' // trim(sub) // '2'
             exit
          end if
       end do
       sg_name_i = l // sg_name0
    case('4mm')   ! Tetragonal
       do i=0,2
          if(fexist(19+i)) then
           ! C4x+ or C4y+ or C4z+
             call get_screw_sub(t(1+i,19+i),4,sub)
             sg_name0= '4' // trim(sub)
           ! IC2x or IC2y or IC2z
             j = mod(i+1,3)
             call get_symbol_of_sym_plane(t(1,26+j),sp)
             sg_name0= trim(sg_name0) // sp
           ! IC2d(40) or IC2c(39) or IC2b(38)
             call get_symbol_of_sym_plane(t(1,40-i),sp)
             sg_name0= trim(sg_name0) // sp
             exit
          end if
       end do
       sg_name_i = l // sg_name0
    case('-42m') ! Tetragonal
       do i=0,2
          if(fexist(43+i)) then
           ! C4x+ or C4y+ or C4z+
             sg_name0= '-4'
             j = mod(i+1,3)
             if(fexist(2+i)) then
              ! C2x or C2y or C2z
                call get_screw_sub(t(1+j,2+j),2,sub)
                sg_name0= trim(sg_name0) // '2' // trim(sub)
              ! IC2d(40) or IC2c(39) or IC2b(38)
                call get_symbol_of_sym_plane(t(1,40-i),sp)
                sg_name0= trim(sg_name0) // sp
             else
              ! IC2x or IC2y or IC2z
                call get_symbol_of_sym_plane(t(1,26+i),sp)
                sg_name0= trim(sg_name0) // sp
              ! C2d(16) or C2c(15) or C2b(14)
                sg_name0= trim(sg_name0) // '2'
             end if
             exit
          end if
       end do
       sg_name_i = l // sg_name0
    case('4/mmm') ! Tetragonal
       do i=0,2
          if(fexist(19+i)) then
           ! C4x+ or C4y+ or C4z+
             call get_screw_sub(t(1+i,19+i),4,sub)
             sg_name0= '4' // trim(sub)
           ! IC2x+ or IC2y+ or IC2z+
             call get_symbol_of_sym_plane(t(1,26+i),sp)
             sg_name0= trim(sg_name0) // "/" // sp
           ! IC2x or IC2y or IC2z
             j = mod(i+1,3)
             call get_symbol_of_sym_plane(t(1,26+j),sp)
             sg_name0= trim(sg_name0) // sp
           ! IC2d(40) or IC2c(39) or IC2b(38)
             call get_symbol_of_sym_plane(t(1,40-i),sp)
             sg_name0= trim(sg_name0) // sp
             exit
          end if
       end do
       sg_name_i = l // sg_name0
    case('3')     ! Trigonal
       if(l=='R') then
          sg_name_i = l // pg_name_i
       else
          ! C3+
          call get_screw_sub(t(3,3),3,sub)
          sg_name_i = l // '3' // trim(sub)
       end if
    case('-3')    ! Trigonal
       sg_name_i = l // '-3'
    case('32')    ! Trigonal
       if(l=='R') then
          sg_name_i = l // pg_name_i
       else
          ! C3+
          call get_screw_sub(t(3,3),3,sub)
          sg_name0 = '3' // trim(sub)
          if(fexist(10)) then
          ! C212
             sg_name0 = trim(sg_name0) // '2'
          else
             sg_name0 = trim(sg_name0) // '1'
          end if
          if(fexist(7)) then
          ! C211
             sg_name0 = trim(sg_name0) // '2'
          else
             sg_name0 = trim(sg_name0) // '1'
          end if
          sg_name_i = l // trim(sg_name0)
       end if
    case('3m','-3m')    ! Trigonal
       if(pg_name_i == '3m') then
          sg_name0 = '3'
       else
          sg_name0 = '-3'
       end if
       if(system == 'cubic') then
          ! C31+ or C32+ or C33+ or C34+
          do i=0,4
             if(fexist(5+i)) then
                j=0
                if(i>1) j=1
                ! C2b for C31+ and C32+, or C2a for C33+ and C34+
                call get_symbol_of_sym_plane(t(1,38-j),sp)
                exit
             end if
          end do
          sg_name0 = trim(sg_name0) // sp
       else
          if(fexist(22)) then
          ! IC212
             call get_symbol_of_sym_plane(t(1,22),sp)
             sg_name0 = trim(sg_name0) // sp
          else if(l == 'P') then
             sg_name0 = trim(sg_name0) // '1'
          end if
          if(fexist(19)) then
          ! IC211
             call get_symbol_of_sym_plane(t(1,19),sp)
             sg_name0 = trim(sg_name0) // sp
          else if(l == 'P') then
             sg_name0 = trim(sg_name0) // '1'
          end if
       end if
       sg_name_i = l // trim(sg_name0)
    case('6')     ! Hexagonal
       call get_screw_sub(t(3,2),6,sub)
       sg_name_i = l // '6' // trim(sub)
    case('-6')    ! Hexagonal
       sg_name_i = l // pg_name_i
    case('6/m')   ! Hexagonal
       call get_screw_sub(t(3,2),6,sub)
       sg_name_i = l // '6' // trim(sub) // "/m"
    case('622')   ! Hexagonal
       call get_screw_sub(t(3,2),6,sub)
       sg_name_i = l // '6' // trim(sub) // "22"
    case('6mm','6/mmm')   ! Hexagonal
       call get_screw_sub(t(3,2),6,sub)
       sg_name0 = '6' // trim(sub)
       if(pg_name_i == '6/mmm') then
          sg_name0 = trim(sg_name0) // "/m"
       end if
       ! IC212
       call get_symbol_of_sym_plane(t(1,22),sp)
       sg_name0 = trim(sg_name0) // sp
       ! IC211
       call get_symbol_of_sym_plane(t(1,19),sp)
       sg_name0 = trim(sg_name0) // sp
       sg_name_i = l // trim(sg_name0)
    case('-6m2')  ! Hexagonal
       sg_name0 = '-6'
       if(fexist(10)) then
        ! C212
          sg_name0 = trim(sg_name0) // '2'
       else
        ! IC212
          call get_symbol_of_sym_plane(t(1,22),sp)
          sg_name0 = trim(sg_name0) // sp
       end if
       if(fexist(7)) then
         ! C211
          sg_name0 = trim(sg_name0) // '2'
       else
       ! IC211
          call get_symbol_of_sym_plane(t(1,19),sp)
          sg_name0 = trim(sg_name0) // sp
       end if
       sg_name_i = l // trim(sg_name0)
    case('23')    ! Cubic
       ! C2x
       call get_screw_sub(t(1,2),2,sub)
       if(sub=='') call get_screw_sub(t(2,2),2,sub)
       if(sub=='') call get_screw_sub(t(3,2),2,sub)
       sg_name0= '2' // trim(sub) // '3'
       sg_name_i = l // trim(sg_name0)
    case('m-3','m-3m')   ! Cubic
       ! IC2x
       call get_symbol_of_sym_plane(t(1,26),sp)
       sg_name0= sp // '-3'
       if(pg_name_i == 'm-3m') then
         ! IC2a
         call get_symbol_of_sym_plane(t(1,37),sp)
         sg_name0= trim(sg_name0) // sp
       end if
       sg_name_i = l // trim(sg_name0)
    case('432')   ! Cubic
       ! C4x+
       call get_screw_sub(t(1,19),4,sub)
       sg_name0= '4' // trim(sub) // '32'
       sg_name_i = l // trim(sg_name0)
    case('-43m')  ! Cubic
       sg_name0= '-43'
       ! IC2a
       if(sum(abs(t(1:3,37))) < eps) then
          sp = 'm'
       else if(l == 'P') then
          sp = 'n'
       else if(l == 'F') then
          sp = 'c'
       else if(l == 'I') then
          sp = 'd'
       end if
       sg_name0= trim(sg_name0) // sp
       sg_name_i = l // trim(sg_name0)
    case default
    end select

  contains
    subroutine get_trans_in_bravais_system(l,t)
      character(len=1), intent(in) :: l
      real(kind=DP), intent(inout) :: t(3,48)

      integer :: i,j
      real(kind=DP) :: p2b(3,3),tt(3)

      p2b = 0.d0
      if(l=='R' .or. l=='P') then
         ! Hexagonal or primitive
         do i=1,3
            p2b(i,i)=1.d0
         end do
      else if(l=='A') then
         ! A-face centered
         p2b(1,1) =  1.d0
         p2b(2,2) =  0.5d0
         p2b(3,3) =  0.5d0
         p2b(2,3) =  0.5d0
         p2b(3,2) = -0.5d0
      else if(l=='B') then
         ! B-face centered
         p2b(2,2) =  1.d0
         p2b(1,1) =  0.5d0
         p2b(3,3) =  0.5d0
         p2b(3,1) =  0.5d0
         p2b(1,3) = -0.5d0
      else if(l=='C') then
         ! C-face centered
         p2b(3,3) =  1.d0
         p2b(1,1) =  0.5d0
         p2b(2,2) =  0.5d0
         p2b(1,2) =  0.5d0
         p2b(2,1) = -0.5d0
      else if(l=='F') then
         ! All-face centered
         p2b(1:3,1:3) =  0.5d0
         do i=1,3
            p2b(i,i)=0.d0
         end do
      else if(l=='I') then
         ! Body centered
         p2b(1:3,1:3) =  0.5d0
         do i=1,3
            p2b(i,i)=-0.5d0
         end do
      end if

      do j=1,48
         do i=1,3
            tt(i) = sum(p2b(i,1:3)*t(1:3,j))
         end do
         t(1:3,j) = tt(1:3)
      end do
    end subroutine get_trans_in_bravais_system

    subroutine get_symbol_of_sym_plane(t,sp)
      real(kind=DP), intent(in) :: t(3)
      character(len=1), intent(out) :: sp

      real(kind=DP), parameter :: eps = 1.d-6
      integer :: i

      if(sum(abs(t(1:3))) < eps) then
         sp = 'm'
      else if(abs(t(1))>eps.and.abs(t(2))<eps.and.abs(t(3))<eps) then
         sp = 'a'
      else if(abs(t(1))<eps.and.abs(t(2))>eps.and.abs(t(3))<eps) then
         sp = 'b'
      else if(abs(t(1))<eps.and.abs(t(2))<eps.and.abs(t(3))>eps) then
         sp = 'c'
      else
         do i=1,3
            if(abs(t(i))>eps) then
               if(abs(abs(t(i))-0.5d0) < eps) then
                  sp = 'n'
               else ! if(abs(abs(t(i))-0.25d0) < eps) then
                  sp = 'd'
               end if
               exit
           end if
        end do
      end if

    end subroutine get_symbol_of_sym_plane

    subroutine get_screw_sub(t,nr,sub)
      real(kind=DP), intent(in) :: t
      integer, intent(in) :: nr
      character(len=2), intent(out) :: sub

      real(kind=DP), parameter :: eps = 1.d-6
      integer :: ns
      ns = int(t*dble(nr)+eps)
      sub = ''
      if(ns == 1) then
        sub = '_1'
      else if(ns == 2) then
        sub = '_2'
      else if(ns == 3) then
        sub = '_3'
      else if(ns == 4) then
        sub = '_4'
      else if(ns == 5) then
        sub = '_5'
      end if
    end subroutine get_screw_sub

  end subroutine get_space_group_name

  subroutine m_CS_SG_print_space_group_name(nfout)
    integer, intent(in) :: nfout

    character(len=3) :: pg_name
    character(len=5) :: pg_name_i
    character(len=2) :: bl_name
    character(len=10) :: sg_name_i
    character(len=9) :: system
    integer :: ipri_t

    if(il>0) then
       system = 'cubic'
    else
       system = 'hexagonal'
    end if
    if(printable) then
      ipri_t = 1
    else
      ipri_t = 0
    end if

    call get_point_group_name(nfout,system,nopr,ig01,pg_name,pg_name_i)
    call get_bravais_lattice(nfout,ipri_t,altv(1,1),altv(1,2),altv(1,3),system,ig01,pg_name_i,bl_name)
    !!$call get_space_group_name(nfout,system,nopr,ig01,tau(1,1,BUCS),.true.,bl_name,pg_name_i,sg_name_i)
!!$       call get_space_group_name(nfout,system,nsym,nopr,tau(1,1,BUCS),.false.,bl_name,pg_name_i,sg_name_i)
!!$    if(.not.allocated(tau)) then
!!$       write(nfout,'(" tau is not allocated")')
!!$    else
       call get_space_group_name(nfout,system,nopr,ig01,tau(1,1,BUCS),.false.,bl_name,pg_name_i,sg_name_i)
!!$    end if

    if(printable) then
       write(nfout,'(13x,"Bravais lattice: ",a,1x,a)') bl_name
       write(nfout,'("Crystallographic point group: ",a,1x,a)') pg_name,pg_name_i
       write(nfout,'(17x,"Space group: ",a)') sg_name_i
    end if
  end subroutine m_CS_SG_print_space_group_name

  subroutine get_kpoint_symmetry(nsym,nopr,rot,kvec,nsymo,nopro,nsymo2,nopro2)
    integer, intent(in) :: nsym
    integer, intent(in) :: nopr(nsym)
    real(kind=DP), intent(in) :: rot(3,3,nsym)
    real(kind=DP), intent(in) :: kvec(3)
    integer, intent(out) :: nsymo
    integer, intent(out) :: nopro(nsym)
    integer, intent(out) :: nsymo2
    integer, intent(out) :: nopro2(nsym)

    real(kind=DP) :: rkv(3),dk(3)
    integer :: i,j
    real(kind=DP) :: eps = 1.d-6

    nsymo = 0
    nsymo2 = 0
    do i=1,nsym
       do j=1,3
          rkv(j) = sum(rot(j,1:3,i)*kvec(1:3))
       end do
       do j=1,3
          dk(j) = sum(altv(1:3,j)*(rkv(1:3)-kvec(1:3)))/PAI2
       end do
       if(sum(abs(dk(1:3))) < eps) then
          nsymo2 = nsymo2 + 1
! ================================= ASMS DEBUG ================= 2013/10/03
!          nopro2(nsymo2) = nopr(i)
          nopro2(nsymo2) = i
! ================================= ASMS DEBUG ================= 2013/10/03
       end if
       call mod1(dk)
       if(sum(abs(dk(1:3))) < eps) then
          nsymo = nsymo + 1
          nopro(nsymo) = nopr(i)
       end if
    end do

  end subroutine get_kpoint_symmetry

  subroutine vectors_length_angle(a1,a2,a3,length,angle)
    implicit none

    real(kind=DP), intent(in) :: a1(3),a2(3),a3(3)
    real(kind=DP), intent(out) :: length(3),angle(3)

    length(1) = sqrt(sum(a1(1:3)**2))
    length(2) = sqrt(sum(a2(1:3)**2))
    length(3) = sqrt(sum(a3(1:3)**2))

    angle(1) = acos(sum(a2(1:3)*a3(1:3))/(length(2)*length(3)))
    angle(2) = acos(sum(a3(1:3)*a1(1:3))/(length(3)*length(1)))
    angle(3) = acos(sum(a1(1:3)*a2(1:3))/(length(1)*length(2)))

  end subroutine vectors_length_angle

  function rational(a) result(b)
    real(kind=DP), intent(in) :: a
    real(kind=DP) :: b

    integer :: i,r
    real(kind=DP) :: eps = 1.d-4
    do i=1,20
       r = nint(a*i)
       b = dble(r)/dble(i)
       if(abs(a-b) < eps) exit
    end do
  end function rational

! ================================== Added by K. Tagami =========
  subroutine ktrational( a,b,ifound )
    real(kind=DP), intent(in) :: a
    real(kind=DP), intent(out) ::b
    integer, intent(inout) :: ifound

    integer :: i,r
    real(kind=DP) :: eps = 1.d-4
    ifound  = -1
    do i=1,20
       r = nint(a*i)
       b = dble(r)/dble(i)
       if(abs(a-b) < eps) then
          ifound = 0
          exit
       endif
    end do

  end subroutine ktrational
! ==================================================================

  subroutine m_CS_gen_opr_rspace_full( system, nsym, rot )
    character(len=9), intent(in) :: system
    integer, intent(in) :: nsym
    real(kind=DP), intent(out) :: rot(3,3,nsym)

    rot = 0.0d0
    if(system .eq. 'cubic') then
       call cubic_point_group
    else if(system .eq. 'hexagonal') then
       call hexagonal_point_group
    end if

  contains

    subroutine cubic_point_group
      ! local variables
      real(kind=DP) :: one=1.d0, zero=0.d0
      integer :: i,j

      ! C31+ = O5
      rot(1,3,5) = one;      rot(2,1,5) = one;      rot(3,2,5) = one

      ! C4x+ = O19
      rot(1,1,19) =  one;    rot(2,3,19) = -one;    rot(3,2,19) =  one

      do i=1,3
         rot(i,i,1) = one  ! E = O1
      end do
      do j=1,3
         do i=1,3
            rot(i,j,22) = rot(j,i,19)  ! C4x- = O22 = O19^t
            rot(i,j,9)  = rot(j,i,5)  ! C31- = O9 = O5^t
         end do
      end do

      call matrix_product(rot(1,1,19),rot(1,1,19),rot(1,1,2)) ! C2x = O2 = O19*O19
      call matrix_product(rot(1,1,5),rot(1,1,19),rot(1,1,13)) ! C2a = O13 = O5*O19
      call matrix_product(rot(1,1,19),rot(1,1,5),rot(1,1,15)) ! C2c = O15 = O19*O5

      call matrix_product(rot(1,1,2),rot(1,1,9),rot(1,1,10))  ! C32- = O10 = O2*O9
      call matrix_product(rot(1,1,13),rot(1,1,15),rot(1,1,11))! C33- = O11 = O13*O15
      call matrix_product(rot(1,1,9),rot(1,1,2),rot(1,1,12))  ! C34- = O12 = O9*O2
      call matrix_product(rot(1,1,2),rot(1,1,5),rot(1,1,8))   ! C34+ = O8 = O2*O5
      call matrix_product(rot(1,1,5),rot(1,1,2),rot(1,1,6))   ! C32+ = O6 = O5*O2
      call matrix_product(rot(1,1,15),rot(1,1,13),rot(1,1,7)) ! C31+ = O7 = O15*O13
      call matrix_product(rot(1,1,13),rot(1,1,2),rot(1,1,21)) ! C4z+ = O21 = O13*O2
      call matrix_product(rot(1,1,2),rot(1,1,13),rot(1,1,24)) ! C4z- = O24 = O2*O13
      call matrix_product(rot(1,1,15),rot(1,1,2),rot(1,1,23)) ! C4y- = O23 = O15*O2
      call matrix_product(rot(1,1,2),rot(1,1,15),rot(1,1,20)) ! C4y+ = O20 = O2*O15

      call matrix_product(rot(1,1,9),rot(1,1,7),rot(1,1,3))  ! C2y = O3 = O9*O7
      call matrix_product(rot(1,1,9),rot(1,1,8),rot(1,1,4))  ! C2z = O4 = O9*O8
      call matrix_product(rot(1,1,15),rot(1,1,6),rot(1,1,14)) ! C2b = O14 = O15*O6
      call matrix_product(rot(1,1,13),rot(1,1,7),rot(1,1,18)) ! C2f = O18 = O13*O7
      call matrix_product(rot(1,1,13),rot(1,1,8),rot(1,1,16)) ! C2d = O16 = O13*O8
      call matrix_product(rot(1,1,7),rot(1,1,24),rot(1,1,17)) ! C2e = O17 = O7*O24

      do i=1,24
         rot(1:3,1:3,i+24) = -rot(1:3,1:3,i) ! add inversion
      end do
    end subroutine cubic_point_group

    subroutine hexagonal_point_group
      ! local variables
      real(kind=DP) :: one=1.d0, zero=0.d0, half=0.5d0, r32
      integer :: i,j

      r32=sqrt(3.d0)*0.5d0

      ! C6+ = O2
      rot(1,1,2) =  half;      rot(1,2,2) = -r32;      rot(2,1,2) =  r32
      rot(2,2,2) =  half;      rot(3,3,2) =  one

      ! C222 = O11
      rot(1,1,11) = -half;      rot(1,2,11) = -r32;      rot(2,1,11) = -r32
      rot(2,2,11) =  half;      rot(3,3,11) = -one

      do i=1,3
         rot(i,i,1) = one  ! E = O1
      end do
      do j=1,3
         do i=1,3
            rot(i,j,6) = rot(j,i,2)  ! C6- = O6 = O2^t
         end do
      end do
      call matrix_product(rot(1,1,2),rot(1,1,2),rot(1,1,3))   ! C3+ = O3 = O2*O2
      call matrix_product(rot(1,1,2),rot(1,1,11),rot(1,1,9))  ! C231 = O9 = O2*O11
      call matrix_product(rot(1,1,11),rot(1,1,2),rot(1,1,7))  ! C211 = O7 = O11*O2

      do j=1,3
         do i=1,3
            rot(i,j,5) = rot(j,i,3)  ! C3- = O5 = O3^t
         end do
      end do
      call matrix_product(rot(1,1,2),rot(1,1,3),rot(1,1,4))  ! C2 = O4 = O2*O3
      call matrix_product(rot(1,1,2),rot(1,1,9),rot(1,1,10)) ! C212 = O10 = O2*O9
      call matrix_product(rot(1,1,7),rot(1,1,3),rot(1,1,8))  ! C221 = O8 = O7*O3
      call matrix_product(rot(1,1,7),rot(1,1,2),rot(1,1,12)) ! C232 = O12 = O7*O2

      do i=1,12
         rot(1:3,1:3,i+12) = -rot(1:3,1:3,i) ! add inversion
      end do

    end subroutine hexagonal_point_group

    subroutine matrix_product(a,b,c)
      real(kind=DP), intent(in) :: a(3,3),b(3,3)
      real(kind=DP), intent(out) :: c(3,3)

      integer :: i,j
      do j=1,3
         do i=1,3
            c(i,j) = sum(a(i,1:3)*b(1:3,j))
         end do
      end do
    end subroutine matrix_product

  end subroutine m_CS_gen_opr_rspace_full

  subroutine m_CS_remove_frac_translation
    use m_Crystal_Structure,  only :  tau, op, ig01

    integer :: iopr, count
    real(kind=DP) :: c1
    integer, allocatable :: ig01_wk(:)
    real(kind=DP), allocatable :: op_wk(:,:,:), tau_wk(:,:,:)

    real(kind=DP), parameter :: eps = 1.0E-4

    allocate( op_wk(3,3,nopr) );   op_wk = 0.0d0
    allocate( tau_wk(3,nopr,2) ); tau_wk = 0.0d0
    allocate( ig01_wk(nopr) );   ig01_wk = 0
!
    count = 0
    Do iopr=1, nopr
       c1 = abs( tau(1,iopr,BUCS) ) + abs( tau(2,iopr,BUCS) ) &
            &                       + abs( tau(3,iopr,BUCS) )
       if ( c1 < eps ) then
          count = count + 1
          op_wk(:,:,count) = op(:,:,iopr)
          tau_wk(:,count,:)  = tau(:,iopr,:)
          ig01_wk(count)    = ig01(iopr)
       endif
    End do
    deallocate( op );   deallocate( tau )

! ---------------------------------------------
    nopr = count
    allocate( op(3,3,nopr) );   op = 0.0d0
    allocate( tau(3,nopr,2) );  tau = 0.0d0
    ig01 = 0

    op(:,:,1:nopr) =  op_wk(:,:,1:nopr)
    tau(:,1:nopr,:) = tau_wk(:,1:nopr,:)
    ig01(1:nopr)    = ig01_wk(1:nopr)
!
    deallocate( op_wk );    deallocate( tau_wk );   deallocate( ig01_wk )

  end subroutine m_CS_remove_frac_translation

  subroutine m_CS_SG_wd_cntn(nfcntn)
    integer, intent(in) :: nfcntn
    integer :: i
    if(.not.allocated(mtprim)) return
    if(mype==0) then
      write(nfcntn, '(a)') tag_symmetry_info
      write(nfcntn,'(2i8)') mnprim,mnsym
      do i=1,mnprim
         write(nfcntn,'(3f20.15)') mtprim(1:3,i)
      enddo
      do i=1,mnsym
         write(nfcntn,'(3f20.15,l2)') mtau(1:3,i),mfsym(i)
      enddo
      if (af/=0) then
        write(nfcntn, '(a)') tag_symmetry_info_af
        write(nfcntn,'(2i8)') mnprim_af,mnsym_af
        do i=1,mnprim_af
           write(nfcntn,'(3f20.15)') mtprim_af(1:3,i)
        enddo
        do i=1,mnsym_af
           write(nfcntn,'(3f20.15,l2)') mtau_af(1:3,i),mfsym_af(i)
        enddo
      endif
    endif
  end subroutine m_CS_SG_wd_cntn

end module m_CS_SpaceGroup
