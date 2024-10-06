module m_ValenceBand_Spectrum

  use m_Control_Parameters,  only : printable, ipriepsilon
  use m_Const_Parameters, only : DP
  use m_Files,  only : nfout, nfpot, m_Files_open_ps_files, m_Files_close_ps_files, &
       &               m_Files_open_ps_file, m_Files_close_ps_file
  use m_Ionic_System,  only : ityp, ivan, iatomn, ntyp
  use m_Parallelization,  only : mype, ierr, MPI_CommGroup

  use m_Orbital_QuantumNum, only : num_orb_index_data, &
       &                           qnum_n_orb_index, qnum_l_orb_index, &
       &                           qnum_t_orb_index, qnum_tau_orb_index
  use mpi

  implicit none
!  include 'mpif.h'

  integer,  allocatable, dimension(:)     :: nppc_data       ! nppc_data(ntyp)
  integer,  allocatable, dimension(:,:)   :: phase_ylm1,phase_ylm2
                                           ! phase_ylm1(ntyp,:),phase_ylm2(ntyp,:)
  integer,  allocatable, dimension(:,:)    :: dipole_tau1,dipole_tau2
                                            ! dipole_tau1(ntyp,:),dipole_tau2(ntyp,:)
  integer                                 :: mnppc
  real(DP), allocatable, dimension(:,:,:)   :: dipole_dxyz_us
                                            ! dipole_dxyz_us(ntyp,:,3)

contains

  subroutine m_VBS_dealloc_dipole_ppc
    deallocate(nppc_data);   deallocate(dipole_dxyz_us)
    deallocate(phase_ylm1);  deallocate(phase_ylm2)
    deallocate(dipole_tau1); deallocate(dipole_tau2)
  end subroutine m_VBS_dealloc_dipole_ppc

  subroutine m_VBS_set_data_ppc_from_pp
    call m_Files_open_ps_files(ivan,iatomn,ntyp,ierr)
    if (ierr/=0) call mpi_stop(nfout)

    if(.not.allocated(nppc_data)) allocate(nppc_data(ntyp)) ; nppc_data=0
    call find_ppc_data_number(nfout)

    call alloc_ptrans_data_array_ek
    call read_ptrans_data_ek

    call m_Files_close_ps_files

  end subroutine m_VBS_set_data_ppc_from_pp

  subroutine m_VBS_set_data_ppc_from_pp_it()
    integer :: it
    if(.not.allocated(nppc_data)) allocate(nppc_data(ntyp)) ; nppc_data=0
    do it=1,ntyp
       call m_Files_open_ps_file(ivan,iatomn,ntyp,it,ierr)
       if (ierr/=0) call mpi_stop(nfout)
       call find_ppc_data_number_it(nfout,it)
       call m_Files_close_ps_file(it)
    enddo
    call alloc_ptrans_data_array_ek()
    do it=1,ntyp
       call m_Files_open_ps_file(ivan,iatomn,ntyp,it,ierr)
       if (ierr/=0) call mpi_stop(nfout)
       call find_ppc_data_number_it(nfout,it)
       call read_ptrans_data_ek_it(it)
       call m_Files_close_ps_file(it)
    enddo
    if(printable) write(nfout,'(1x,"!* ntyp = ",i3)') ntyp
  end subroutine m_VBS_set_data_ppc_from_pp_it

  subroutine find_ppc_data_number(nfout)
!
!   find data number of core-repair term
!   nppc_data(it) : data number for it-th type pseudopotential
!
    integer,intent(in) :: nfout
    integer            :: it

    do it=1, ntyp
       call find_ppc_data_number_it(nfout,it)
    end do

 end subroutine find_ppc_data_number

 subroutine find_ppc_data_number_it(nfout,it)
    integer, intent(in) :: nfout,it
    integer, save :: nfpp=0
    integer :: nptrans
    if(mype==0)then
       nfpp=nfpp+1
       call find_dipole_section(nfpot(nfpp),nfout,it,nptrans)
       if(printable) &
       & write(nfout,'(1x,"! PP transition moment correction data : it = ",i3, &
       & 2x,"number of data read from PP file = ",i3)') it, nptrans
       nppc_data(it)=nptrans
    endif
    if(it==ntyp) then
       call mpi_bcast(nppc_data,ntyp,mpi_integer,0,MPI_CommGroup,ierr) ! MPI
       nfpp=0
    endif
 end subroutine find_ppc_data_number_it

 subroutine find_dipole_section(nfp,nfout,it,nptrans)
!
!   find dipole section in gncpp2 potential file
!
    integer,intent(in)   :: nfp, nfout,it
    integer, intent(out) :: nptrans
    integer              :: idipole
    character(len=10)    :: line

    if(mype /= 0) return

! == KT_add === 2014/07/23
    rewind(nfp)
! ============= 2014/07/23

10  read(nfp,'(a10)',end=20) line
    idipole=index(line,'DIPOLE')
    if(idipole/=0) then
       read(nfp,*) nptrans
       if(nptrans==0) then
          if(printable) &
          & write(nfout,'(1x,"!* number of dipole data is zero for atom type ",i3,"   m_Epsilon_ek STOP")') it
          call phase_error_with_msg(nfout,'number of dipole data is zero',__LINE__,__FILE__)
       end if
       return
    end if
    goto 10
20  if(printable) write(*,'(1x,"!* dipole section is not found in pseudopotential file")')
    call phase_error_with_msg(nfout,'dipole section is not found in pseudopotential file',__LINE__,__FILE__)
 end subroutine find_dipole_section

 subroutine alloc_ptrans_data_array_ek
    implicit none
!
!   allocate data array for core-repair term
!
    integer :: it
! find maximum of nppc_data(it): it=1-ntyp
    if(mype == 0) then
       if(ntyp==1) then
          mnppc=nppc_data(1)
       else
          mnppc=nppc_data(1)
          do it=1,ntyp
             if(mnppc<nppc_data(it)) mnppc=nppc_data(it)
          end do
       end if
    end if
    call mpi_bcast(mnppc,1,mpi_integer,0,MPI_CommGroup, ierr) ! MPI
    if(printable) write(nfout,'(1x,"!* mnppc = ",i3)') mnppc
! allocate data array
    if(.not.allocated(dipole_dxyz_us)) allocate(dipole_dxyz_us(ntyp,mnppc,3)); dipole_dxyz_us=0.0d0
    if(.not.allocated(dipole_tau1))    allocate(dipole_tau1(ntyp,mnppc)); dipole_tau1=0
    if(.not.allocated(dipole_tau2))    allocate(dipole_tau2(ntyp,mnppc)); dipole_tau2=0
    if(.not.allocated(phase_ylm1))     allocate(phase_ylm1(ntyp,mnppc)); phase_ylm1=0
    if(.not.allocated(phase_ylm2))     allocate(phase_ylm2(ntyp,mnppc)); phase_ylm2=0
 end subroutine alloc_ptrans_data_array_ek

 subroutine read_ptrans_data_ek_it(it)
    integer, intent(in) :: it
    integer, save :: nfpp=0
    if(mype==0) then
       nfpp=nfpp+1
       if(printable) &
       & write(nfout,'(1x,"! ppc data for it = ",i4," are read: nppc = ",i4)') it,nppc_data(it)
! ====================== KT_mod ================================ 13.0S +2014/07/23
!          call read_ptrans_data_ek_core(nfpot(nfpp),nfout,it)
!
       if ( num_orb_index_data(it) > 0 ) then
          call read_ptrans_data_ek_core( nfpot(nfpp), nfout, it )
       else
          call read_ptrans_data_ek_core_kt( nfpot(nfpp), nfout, it )
       endif
! ============================================================== 13.0S +2014/07/23
    endif

    if(it==ntyp) then
       call mpi_bcast(dipole_dxyz_us,ntyp*mnppc*3,mpi_double_precision,0,MPI_CommGroup,ierr) ! MPI
       call mpi_bcast(dipole_tau1,ntyp*mnppc,mpi_integer,0,MPI_CommGroup,ierr)      ! MPI
       call mpi_bcast(dipole_tau2,ntyp*mnppc,mpi_integer,0,MPI_CommGroup,ierr)      ! MPI
       call mpi_bcast(phase_ylm1,ntyp*mnppc,mpi_integer,0,MPI_CommGroup,ierr)       ! MPI
       call mpi_bcast(phase_ylm2,ntyp*mnppc,mpi_integer,0,MPI_CommGroup,ierr)       ! MPI
       nfpp=0
    endif
 end subroutine read_ptrans_data_ek_it

 subroutine read_ptrans_data_ek
    implicit none
!
!   read core-repair term
!
    integer :: it
!!!    if(mype == 0)  then      !!! KT_mod: 2014/12/13
       do it =1, ntyp
          call read_ptrans_data_ek_it(it)
       end do
!!!    end if                   !!! KT_mod : 2014/12/13

    if(printable) write(nfout,'(1x,"!* ntyp = ",i3)') ntyp
  end subroutine read_ptrans_data_ek

     subroutine read_ptrans_data_ek_core(nfp,nfout,it)
       implicit none
!
!      read core-repair term pij(I)
!      cf. H. Kageshima and K. Shiraishi,Phys. Rev. B vol.56, 14985 (1997)
!      pij(I)
!      i, j :: atomic orbital index (CIAO)
!      i ->  ltmltm(n1,l1,t1,m1)  j-> ltmltm(n2,l2,t1,m1)  ltm: compound index
!      n :: principal quantum number    l :: azimuthal quantum number
!      m :: magnetic quantum number     t :: energy reference index
!      dipole_dxyz_us(it,ltmltm,ixyz)     :: core repair term for it-th type pseudopotential in ixyz direction
!                                         :: ixyz=1 -> x; ixyz=2 -> y;  ixyz= 3 -> z
!      phase_ylm  :: Ylm index for phase
!      dipole_tau :: reference index
!
!
       integer, intent(in) :: nfp,nfout,it
       integer :: n1,l1,t1,m1,n2,l2,t2,m2
       integer :: ltmltm,lmax, j1

       do ltmltm = 1, nppc_data(it)
!         n1 = n1_dipole_lm_us(ltmltm) ; n2 = n2_dipole_lm_us(ltmltm)
!         l1 = l1_dipole_lm_us(ltmltm) ; l2 = l2_dipole_lm_us(ltmltm)
!         t1 = t1_dipole_lm_us(ltmltm) ; t2 = t2_dipole_lm_us(ltmltm)
!         m1 = m1_dipole_lm_us(ltmltm) ; m2 = m2_dipole_lm_us(ltmltm)
          read(nfp,53) n1,l1,t1,m1,n2,l2,t2,m2, &
               &    dipole_dxyz_us(it,ltmltm,1), &
               &    dipole_dxyz_us(it,ltmltm,2), &
               &    dipole_dxyz_us(it,ltmltm,3), &
               &    phase_ylm1(it,ltmltm), phase_ylm2(it,ltmltm)
          dipole_tau1(it,ltmltm)=t1
          dipole_tau2(it,ltmltm)=t2

! === KT_add === 2014/07/23
          Do j1=1, num_orb_index_data(it)
             if ( n1 == qnum_n_orb_index(it,j1) &
                  &    .and. l1 == qnum_l_orb_index(it,j1) &
                  &    .and. t1 == qnum_t_orb_index(it,j1) ) then
                dipole_tau1( it,ltmltm ) = qnum_tau_orb_index(it,j1)
                exit
             end if
          End Do
          Do j1=1, num_orb_index_data(it)
             if ( n2 == qnum_n_orb_index(it,j1) &
                  &    .and. l2 == qnum_l_orb_index(it,j1) &
                  &    .and. t2 == qnum_t_orb_index(it,j1) ) then
                dipole_tau2( it,ltmltm ) = qnum_tau_orb_index(it,j1)
                exit
             end if
          End Do
! ============ 2014/07/23

          if (ipriepsilon>=2.and.printable ) then
             write(nfout,53) n1,l1,dipole_tau1(it,ltmltm),m1,n2,l2,dipole_tau2(it,ltmltm),m2, &
                  & dipole_dxyz_us(it,ltmltm,1),dipole_dxyz_us(it,ltmltm,2),dipole_dxyz_us(it,ltmltm,3), &
                  & phase_ylm1(it,ltmltm),phase_ylm2(it,ltmltm)
          end if
       end do

53     format(1x,8i3,3e18.10,2i3)

     end subroutine read_ptrans_data_ek_core

! ========================== KT_add =========================== 13.0S
! ****************************************************
!  Note: In the DIPOLE section of pp files, the variable tau ("t1, "t2") is
!        defined for each principal quantum number. For example,
!               5s: tau=1, 2,   6s: tau=1, 2
!
!        On the other hand, in the PHASE, the principal quantum numbers are not
!        treated explicitly. In the above exmaple,
!               s-orbital:  tau=1, 2, 3, 4
!
!        The subroutine read_ptrans_data_ek_core is not correct in the cases
!        when an l-quantum number corresponds to multiple pricipal numbers.
!
!        In the following subroutine, we could overcome this limitation by some tricks..
! ****************************************************
!
     subroutine read_ptrans_data_ek_core_kt(nfp,nfout,it)
       implicit none
       integer, intent(in) :: nfp,nfout,it

       integer :: n1, l1, t1, m1, n2, l2, t2, m2
       integer :: ltmltm, lval, lmax, nmax
       integer :: ndata, count, i
       integer :: old_n1, old_t1, old_n2, old_t2
!
       integer, allocatable :: data_n1(:), data_n2(:), data_l1(:), data_l2(:)
       integer, allocatable :: data_t1(:), data_t2(:), data_m1(:), data_m2(:)
       integer, allocatable :: data_ylm1(:), data_ylm2(:), igeta(:,:)
       real(kind=DP), allocatable :: data_dipol(:,:)

       ndata = nppc_data(it)

       lmax = 3         ! s:0, p:1,  d:2,  f:3
       nmax = 8         ! 8s, 8p, 8d .. etc.

       allocate( data_n1(ndata) );  allocate( data_n2(ndata) )
       allocate( data_l1(ndata) );  allocate( data_l2(ndata) )
       allocate( data_m1(ndata) );  allocate( data_m2(ndata) )
       allocate( data_t1(ndata) );  allocate( data_t2(ndata) )
       allocate( data_dipol(ndata,3) );
       allocate( data_ylm1(ndata) ); allocate( data_ylm2(ndata) );

       allocate( igeta(nmax,0:lmax) );  igeta = -100

       do i=1, ndata
          read(nfp,53) data_n1(i), data_l1(i), data_t1(i), data_m1(i), &
               &       data_n2(i), data_l2(i), data_t2(i), data_m2(i), &
               &       data_dipol(i,1), data_dipol(i,2), data_dipol(i,3), &
               &       data_ylm1(i), data_ylm2(i)
       End do

       Do lval=0, lmax
          count = 0
          Do i=1, ndata
             if ( data_l1(i) == lval ) then
                n1 = data_n1(i);  t1 = data_t1(i)

                if ( count==0 ) then
                   count = count +1
                   if ( igeta(n1,lval) == -100 ) igeta(n1,lval) = 0
                else
                   if ( old_n1 /=n1 .and. igeta(n1,lval)== -100 ) igeta(n1,lval)=old_t1
                endif
                old_n1 = n1;  old_t1 = t1
             endif
          End do

          count = 0
          Do i=1, ndata
             if ( data_l2(i) == lval ) then
                n2 = data_n2(i);  t2 = data_t2(i)

                if ( count==0 ) then
                   count = count +1
                   if ( igeta(n2,lval) == -100 ) igeta(n2,lval) = 0
                else
                   if ( old_n2 /=n2 .and. igeta(n2,lval)== -100 ) igeta(n2,lval)=old_t2
                endif
                old_n2 = n2;  old_t2 = t2
             endif
          End do
       End do

       Do i=1, ndata
          n1 = data_n1(i);  n2 = data_n2(i)
          l1 = data_l1(i);  l2 = data_l2(i)
          t1 = data_t1(i);  t2 = data_t2(i)
          m1 = data_m1(i);  m2 = data_m2(i)

          dipole_tau1(it, i) = data_t1(i) +igeta(n1,l1)
          dipole_tau2(it, i) = data_t2(i) +igeta(n2,l2)

          dipole_dxyz_us(it, i, 1) = data_dipol(i,1)
          dipole_dxyz_us(it, i, 2) = data_dipol(i,2)
          dipole_dxyz_us(it, i, 3) = data_dipol(i,3)

          phase_ylm1(it, i) = data_ylm1(i)
          phase_ylm2(it, i) = data_ylm2(i)

          if (ipriepsilon>=2.and.printable) then
             write(nfout,53) n1, l1, dipole_tau1(it,i), m1, &
                  &          n2, l2, dipole_tau2(it,i), m2, &
                  &          dipole_dxyz_us(it,i,1), &
                  &          dipole_dxyz_us(it,i,2), &
                  &          dipole_dxyz_us(it,i,3), &
                  &          phase_ylm1(it,i), phase_ylm2(it,i)
          end if
       end do
53     format(1x,8i3,3e18.10,2i3)

       deallocate( data_n1 );  deallocate( data_n2 )
       deallocate( data_l1 );  deallocate( data_l2 )
       deallocate( data_m1 );  deallocate( data_m2 )
       deallocate( data_t1 );  deallocate( data_t2 )
       deallocate( data_dipol );
       deallocate( data_ylm1 ); deallocate( data_ylm2 );
       deallocate( igeta )

     end subroutine read_ptrans_data_ek_core_kt
! ================================================================== 13.0S

  subroutine m_VBS_find_ptrans_index_ek( it, nspher1, nspher2, tau1, tau2, index, ifact )
    integer, intent(in) :: it, nspher1, nspher2, tau1, tau2
    integer, intent(out) :: index, ifact

    integer :: i
    integer :: nspher10, nspher20, tau10, tau20

    index = 0; ifact = 0
    if ( nspher1 > nspher2 ) then
       nspher20 = nspher1;   tau20 = tau1
       nspher10 = nspher2;   tau10 = tau2;    ifact = -1
    else
       nspher10 = nspher1;   tau10 = tau1;
       nspher20 = nspher2;   tau20 = tau2;    ifact = 1
    endif

    do i=1, nppc_data(it)
       if ( phase_ylm1(it,i) == nspher10 .and. phase_ylm2(it,i) == nspher20 ) then
          if ( dipole_tau1(it,i) == tau10 .and. dipole_tau2(it,i) == tau20 ) then
             index = i
             exit
          end if
       end if
    end do

  end subroutine m_VBS_find_ptrans_index_ek

end module m_ValenceBand_Spectrum
