!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: gnrt_k0_n, get_trmat1, gnrt_k0, get_trmat1, read_coordsystem,
!            ch4ksample, read_nkpnt, hexm0, fccm0, bccm0, sccm0, kpmwbz0,
!            kpmsf0, set_kfromtemporaryarray, get_trmat, read_konly,
!            down2wayksample, readk0, get_trmat, gen_gamma_p0
!
!  AUTHOR(S): T. Yamasaki, H. Mizouchi,   August/20/2003
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
! $Id: b_Kpoints.F90 570 2017-04-21 20:34:50Z yamasaki $
subroutine gen_gamma_p0(paramset,knv3,kv3,vkxyz,qwgt)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: gen_gamma_p0
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

  use m_Const_Parameters, only : CARTS, CRDTYP, DP
  implicit none
  logical, intent(in)        :: paramset
  integer, intent(in)        :: knv3
  integer, intent(out)       :: kv3
  real(kind=DP), intent(out) :: vkxyz(knv3,3,CRDTYP)
  real(kind=DP), intent(out) :: qwgt(knv3)

  kv3 = 1
  if(paramset) return

  vkxyz(1,1:3,CARTS) = 0.d0
  qwgt (1) = 1.d0
end subroutine gen_gamma_p0

subroutine readk0(paramset,nfout,ipri_kp,knv3,rltv,nfkpoint,nfmatbp&
     & ,kv3,vkxyz,qwgt)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: readk0
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!  FURTHER MODIFICATION: T. Yamasaki, October/10/2012
!=======================================================================

  use m_Const_Parameters, only : CARTS, CRDTYP, DP
  use m_ErrorMessages,    only : EOF_REACHED, FORMAT_ERROR
  implicit none
  logical, intent(in)        :: paramset
  integer, intent(in)        :: nfout, ipri_kp, knv3
  real(kind=DP), intent(in)  :: rltv(3,3)
  integer, intent(in)        :: nfkpoint, nfmatbp
  integer, intent(out)       :: kv3
  real(kind=DP), intent(out) :: vkxyz(knv3,3,CRDTYP),qwgt(knv3)

  real(kind=DP),dimension(3) ::fv
  real(kind=DP)              :: weightk
  integer                    :: nn,nk,nw,nwei,nv(3)
  integer                    :: fn_number_of_words
  integer,     parameter     :: len_str = 80
  character(len=len_str)     :: str
  integer                    :: ierr

  real(kind=DP), pointer, dimension(:,:) :: trmat,trbp,trpb,mat1,mat2 

  interface
     subroutine phase_error_wo_filename(ierr,nfout,nf,line,modulefile)
       integer, intent(in) :: ierr,nfout
       integer, intent(in), optional :: nf
       integer, intent(in), optional :: line
       character(len=*),intent(in), optional :: modulefile
     end subroutine phase_error_wo_filename
  end interface

  rewind nfkpoint
  call skip_commentlines(nfkpoint,ipri_kp,nfout,ierr)
  if(ierr == EOF_REACHED) then
!!$     write(nfout,'(" __LINE__ = ",i8)')  __LINE__
!!$     write(nfout,'(" __FILE__ = ",a)') __FILE__
#ifdef DEBUG_ERRORS
     call phase_error_wo_filename(ierr,nfout,nf=nfkpoint,line=__LINE__,modulefile=__FILE__)
#else
     call phase_error_wo_filename(ierr,nfout,nf=nfkpoint)
#endif
     call phase_error_with_msg(nfout,' phase_error',__LINE__,__FILE__)
  end if

  read(nfkpoint,'(a80)') str

  select case (fn_number_of_words(str))
  case (0)
     call phase_error_with_msg(nfout,' error in reading of the nfkpoint file',__LINE__,__FILE__)
  case (1)
     read(str,*) kv3
     nwei = kv3
  case (2)
     read(str,*) kv3, nwei
  case (3)
     kv3 = -1
     nwei = -1
  case default
     kv3 = -1
     nwei = -1
!!$     read(str,*) kv3, nwei
  end select
!!$  read(nfkpoint,*) kv3, nwei
  if(paramset .and. kv3 == -1) then
     kv3 = 1
     do while (.true.) 
        if(kv3 > 1) read(nfkpoint,'(a80)',end=1002,err=1002) str
        select case (fn_number_of_words(str))
        case (0,1,2)
           goto 1002
        case (3)
           read(str,*) fv(1), fv(2), fv(3)
        case (4)
           read(str,*) nv(1),nv(2),nv(3),nk
        case (5)
           read(str,*) nv(1),nv(2),nv(3),nk,nw
        case default
           read(str,*) nv(1),nv(2),nv(3),nk,nw
        end select
        kv3 = kv3+1
     end do
1002 continue
     kv3 = kv3-1 
     if(kv3 < 1) then
        call phase_error_with_msg(nfout,' kv3 is smaller than 1 in reading of nfkpoint',__LINE__,__FILE__)
     end if
  end if

  if(paramset) return

  allocate(trmat(3,3),trbp(3,3),trpb(3,3),mat1(3,3),mat2(3,3))

  call get_trmat(rltv,trbp,trpb,mat1,mat2,trmat,nfmatbp)

  if(ipri_kp >= 1) then
     write(nfout,'(" <<readk0 >>")')
     write(nfout,'(" rltv")')
     do nk = 1, 3
        write(nfout,'(3f12.6)') (rltv(nn,nk),nn=1,3)
     end do
     write(nfout,'(" trmat")')
     do nk = 1, 3
        write(nfout,'(3f12.6)') (trmat(nn,nk),nn=1,3)
     end do
  end if

  if(kv3 == -1) then
     kv3 = 1
     do while (.true.) 
        if(kv3 > 1) read(nfkpoint,'(a80)',end=1001,err=1001) str
        write(nfout,'(a80)') str
        select case (fn_number_of_words(str))
        case (0,1,2)
           goto 1001
!!$           stop ' error in reading of k-coordinates from the nfkpoint file'
        case (3)
           fv = 0.d0
           read(str,*) fv(1), fv(2), fv(3)
           write(nfout,'(" fv(1:3) = ",3f8.4)') fv(1:3)
           vkxyz(kv3,1:3,CARTS) = matmul(trmat,fv)
           if(nwei >= 1) then
              qwgt(kv3) = 1.d0/dble(nwei)
           else
              qwgt(kv3) = 1.d0
           end if
        case (4)
           read(str,*) nv(1),nv(2),nv(3),nk
           vkxyz(kv3,1:3,CARTS) = matmul(trmat,nv)/dble(nk)
           if(nwei >= 1) then
              qwgt(kv3) = 1.d0/dble(nwei)
           else
              qwgt(kv3) = 1.d0
           end if
        case (5)
           read(str,*) nv(1),nv(2),nv(3),nk,nw
           vkxyz(kv3,1:3,CARTS) = matmul(trmat,nv)/dble(nk)
           if(nwei >= 1) then
              qwgt(kv3) = dble(nw)/dble(nwei)
           else
              qwgt(kv3) = dble(nw)
           end if
        case default
           read(str,*) nv(1),nv(2),nv(3),nk,nw
           vkxyz(kv3,1:3,CARTS) = matmul(trmat,nv)/dble(nk)
           if(nwei >= 1) then
              qwgt(kv3) = dble(nw)/dble(nwei)
           else
              qwgt(kv3) = dble(nw)
           end if
        end select
        kv3 = kv3+1
     end do
1001 continue
     kv3 = kv3-1 
     if(kv3 < 1) then
        call phase_error_with_msg(nfout,' kv3 is smaller than 1 in reading of nfkpoint',__LINE__,__FILE__)
     end if
     if(nwei < 1) then
        weightk = 0.d0
        do nn = 1, kv3
           weightk = weightk + qwgt(nn)
        end do
        do nn = 1, kv3
           qwgt(nn) = qwgt(nn)/weightk
        end do
     end if
  else
     do nn = 1, kv3
        read(nfkpoint,'(a80)',end=1004, err=1003) str
        select case (fn_number_of_words(str))
        case (0,1,2)
           call phase_error_with_msg(nfout,' error in reading of k-coordinates from the nfkpoint file',__LINE__,__FILE__)
        case (3)
           fv = 0.d0
           read(str,*) fv(1), fv(2), fv(3)
           vkxyz(nn,1:3,CARTS) = matmul(trmat,fv)
           qwgt(nn) = 1.d0/dble(nwei)
        case (4)
           read(str,*) nv(1),nv(2),nv(3),nk
           vkxyz(nn,1:3,CARTS) = matmul(trmat,nv)/dble(nk)
           qwgt(nn) = 1.d0/dble(nwei)
        case (5)
           read(str,*) nv(1),nv(2),nv(3),nk,nw
           vkxyz(nn,1:3,CARTS) = matmul(trmat,nv)/dble(nk)
           qwgt(nn) = dble(nw)/dble(nwei)
        case default
           read(str,*) nv(1),nv(2),nv(3),nk,nw
           vkxyz(nn,1:3,CARTS) = matmul(trmat,nv)/dble(nk)
           qwgt(nn) = dble(nw)/dble(nwei)
        end select
     enddo
     goto 1010
1004 ierr = EOF_REACHED
#ifdef DEBUG_ERRORS
     call phase_error_wo_filename(ierr,nfout,nf=nfkpoint,line=__LINE__,modulefile=__FILE__)
#else
     call phase_error_wo_filename(ierr,nfout,nf=nfkpoint)
#endif

1003 ierr = FORMAT_ERROR
#ifdef DEBUG_ERRORS
     call phase_error_wo_filename(ierr,nfout,nf=nfkpoint,line=__LINE__,modulefile=__FILE__)
#else
     call phase_error_wo_filename(ierr,nfout,nf=nfkpoint)
#endif
1010 continue
  end if

  deallocate(trmat,trbp,trpb,mat1,mat2 )

  if(ipri_kp >= 1) then
     write(nfout,'(/,"  << readk0 >>")')
     write(nfout,*) ' !Total Generated Kpoints = ',kv3
     do nk = 1, kv3
        write(nfout,'(i4," ",3f12.6)') nk,(vkxyz(nk,nn,CARTS),nn=1,3)
     end do
  end if

end subroutine readk0

subroutine readk0_for_band_unfolding( paramset, nfout, ipri_kp, knv3, &
     &                                rltv_refcell, altv_supercell, &
     &                                nfkpoint,nfmatbp, kv3, vkxyz, qwgt )
  use m_Const_Parameters, only : CARTS, CRDTYP, DP, PAI2, BUCS
  use m_ErrorMessages,    only : EOF_REACHED, FORMAT_ERROR

  implicit none

  logical, intent(in)        :: paramset
  integer, intent(in)        :: nfout, ipri_kp, knv3
  real(kind=DP), intent(in)  :: rltv_refcell(3,3), altv_supercell(3,3)
  integer, intent(in)        :: nfkpoint, nfmatbp
  integer, intent(out)       :: kv3
  real(kind=DP), intent(out) :: vkxyz(knv3,3,CRDTYP),qwgt(knv3)

  real(kind=DP),dimension(3) ::fv
  real(kind=DP)              :: weightk
  integer                    :: nn,nk,nw,nwei,nv(3)
  integer                    :: fn_number_of_words
  integer,     parameter     :: len_str = 80
  character(len=len_str)     :: str
  integer                    :: ierr

  integer :: i, j, k
  real*8 c1

  real(kind=DP), pointer, dimension(:,:) :: trmat,trbp,trpb,mat1,mat2

  interface
     subroutine phase_error_wo_filename(ierr,nfout,nf,line,modulefile)
       integer, intent(in) :: ierr,nfout
       integer, intent(in), optional :: nf
       integer, intent(in), optional :: line
       character(len=*),intent(in), optional :: modulefile
     end subroutine phase_error_wo_filename
  end interface

  rewind nfkpoint
  call skip_commentlines(nfkpoint,ipri_kp,nfout,ierr)
  if(ierr == EOF_REACHED) then
!!$     write(nfout,'(" __LINE__ = ",i8)')  __LINE__
!!$     write(nfout,'(" __FILE__ = ",a)') __FILE__
#ifdef DEBUG_ERRORS
     call phase_error_wo_filename(ierr,nfout,nf=nfkpoint,line=__LINE__,modulefile=__FILE__)
#else
     call phase_error_wo_filename(ierr,nfout,nf=nfkpoint)
#endif
     call phase_error_with_msg(nfout,' phase_error',__LINE__,__FILE__)
  end if

  read(nfkpoint,'(a80)') str

  select case (fn_number_of_words(str))
  case (0)
     call phase_error_with_msg(nfout,' error in reading of the nfkpoint file',__LINE__,__FILE__)
  case (1)
     read(str,*) kv3
     nwei = kv3
  case (2)
     read(str,*) kv3, nwei
  case (3)
     kv3 = -1
     nwei = -1
  case default
     kv3 = -1
     nwei = -1
!!$     read(str,*) kv3, nwei
  end select
!!$  read(nfkpoint,*) kv3, nwei
  if(paramset .and. kv3 == -1) then
     kv3 = 1
     do while (.true.)
        if(kv3 > 1) read(nfkpoint,'(a80)',end=1002,err=1002) str
        select case (fn_number_of_words(str))
        case (0,1,2)
           goto 1002
        case (3)
           read(str,*) fv(1), fv(2), fv(3)
        case (4)
           read(str,*) nv(1),nv(2),nv(3),nk
        case (5)
           read(str,*) nv(1),nv(2),nv(3),nk,nw
        case default
           read(str,*) nv(1),nv(2),nv(3),nk,nw
        end select
        kv3 = kv3+1
     end do
1002 continue
     kv3 = kv3-1
     if(kv3 < 1) then
        call phase_error_with_msg(nfout,' kv3 is smaller than 1 in reading of nfkpoint',__LINE__,__FILE__)
     end if
  end if

  if(paramset) return

 allocate(trmat(3,3),trbp(3,3),trpb(3,3),mat1(3,3),mat2(3,3))

  call get_trmat(rltv_refcell,trbp,trpb,mat1,mat2,trmat,nfmatbp)
!
  if(ipri_kp >= 2) then
     write(nfout,'(" <<readk0_band_unfolding >>")')
     write(nfout,'(" rltv")')
     do nk = 1, 3
        write(nfout,'(3f12.6)') (rltv_refcell(nn,nk),nn=1,3)
     end do
     write(nfout,'(" trmat")')
     do nk = 1, 3
        write(nfout,'(3f12.6)') (trmat(nn,nk),nn=1,3)
     end do
  end if

  if(kv3 == -1) then
     kv3 = 1
     do while (.true.)
        if(kv3 > 1) read(nfkpoint,'(a80)',end=1001,err=1001) str
        write(nfout,'(a80)') str
        select case (fn_number_of_words(str))
        case (0,1,2)
           goto 1001
!!$           stop ' error in reading of k-coordinates from the nfkpoint file'
        case (3)
           fv = 0.d0
           read(str,*) fv(1), fv(2), fv(3)
           write(nfout,'(" fv(1:3) = ",3f8.4)') fv(1:3)
           vkxyz(kv3,1:3,CARTS) = matmul(trmat,fv)
           if(nwei >= 1) then
              qwgt(kv3) = 1.d0/dble(nwei)
           else
              qwgt(kv3) = 1.d0
           end if
        case (4)
           read(str,*) nv(1),nv(2),nv(3),nk
           vkxyz(kv3,1:3,CARTS) = matmul(trmat,nv)/dble(nk)
           if(nwei >= 1) then
              qwgt(kv3) = 1.d0/dble(nwei)
           else
              qwgt(kv3) = 1.d0
           end if
        case (5)
           read(str,*) nv(1),nv(2),nv(3),nk,nw
           vkxyz(kv3,1:3,CARTS) = matmul(trmat,nv)/dble(nk)
           if(nwei >= 1) then
              qwgt(kv3) = dble(nw)/dble(nwei)
           else
              qwgt(kv3) = dble(nw)
           end if
        case default
           read(str,*) nv(1),nv(2),nv(3),nk,nw
           vkxyz(kv3,1:3,CARTS) = matmul(trmat,nv)/dble(nk)
           if(nwei >= 1) then
              qwgt(kv3) = dble(nw)/dble(nwei)
           else
              qwgt(kv3) = dble(nw)
           end if
        end select
        kv3 = kv3+1
     end do
1001 continue
     kv3 = kv3-1
     if(kv3 < 1) then
        call phase_error_with_msg(nfout,' kv3 is smaller than 1 in reading of nfkpoint',__LINE__,__FILE__)
     end if
     if(nwei < 1) then
        weightk = 0.d0
        do nn = 1, kv3
           weightk = weightk + qwgt(nn)
        end do
        do nn = 1, kv3
           qwgt(nn) = qwgt(nn)/weightk
        end do
     end if
  else
     do nn = 1, kv3
        read(nfkpoint,'(a80)',end=1004, err=1003) str
        select case (fn_number_of_words(str))
        case (0,1,2)
           call phase_error_with_msg(nfout,' error in reading of k-coordinates from the nfkpoint file',__LINE__,__FILE__)
        case (3)
           fv = 0.d0
           read(str,*) fv(1), fv(2), fv(3)
           vkxyz(nn,1:3,CARTS) = matmul(trmat,fv)
           qwgt(nn) = 1.d0/dble(nwei)
        case (4)
           read(str,*) nv(1),nv(2),nv(3),nk
           vkxyz(nn,1:3,CARTS) = matmul(trmat,nv)/dble(nk)
           qwgt(nn) = 1.d0/dble(nwei)
        case (5)
           read(str,*) nv(1),nv(2),nv(3),nk,nw
           vkxyz(nn,1:3,CARTS) = matmul(trmat,nv)/dble(nk)
           qwgt(nn) = dble(nw)/dble(nwei)
        case default
           read(str,*) nv(1),nv(2),nv(3),nk,nw
           vkxyz(nn,1:3,CARTS) = matmul(trmat,nv)/dble(nk)
           qwgt(nn) = dble(nw)/dble(nwei)
        end select
     enddo
     goto 1010
1004 ierr = EOF_REACHED
#ifdef DEBUG_ERRORS
     call phase_error_wo_filename(ierr,nfout,nf=nfkpoint,line=__LINE__,modulefile=__FILE__)
#else
     call phase_error_wo_filename(ierr,nfout,nf=nfkpoint)
#endif

1003 ierr = FORMAT_ERROR
#ifdef DEBUG_ERRORS
     call phase_error_wo_filename(ierr,nfout,nf=nfkpoint,line=__LINE__,modulefile=__FILE__)
#else
     call phase_error_wo_filename(ierr,nfout,nf=nfkpoint)
#endif
1010 continue
  end if

! -- back to supercell --
  Do i=1, kv3
     Do j=1, 3
        c1 = 0.0d0
        Do k=1, 3
           c1 = c1 + altv_supercell(k,j) *vkxyz(i,k,CARTS)
        End Do
        vkxyz(i,j,BUCS) = c1 /PAI2
     End do
  End Do

  deallocate(trmat,trbp,trpb,mat1,mat2 )

  if(ipri_kp >= 2) then
     write(nfout,'(/,"  << readk0_band_unfolding >>")')
     write(nfout,*) ' !Total Generated Kpoints = ',kv3
     do nk = 1, kv3
        write(nfout,'(i4," ",3f12.6)') nk,(vkxyz(nk,nn,CARTS),nn=1,3)
     end do
  end if

end subroutine readk0_for_band_unfolding


subroutine get_trmat(rltv,trbp,trpb,mat1,mat2,trmat,nfmatbp)
  use m_Const_Parameters, only : DP
  implicit none
  real(kind=DP), intent(in), dimension(3,3) :: rltv
  real(kind=DP), dimension(3,3) ::              trbp,trpb,mat1,mat2
  real(kind=DP), intent(out), dimension(3,3) :: trmat
  integer, intent(in) :: nfmatbp
  !    make translation matrix  trpb (P -> B)
  integer  :: i
  logical :: open_check

  trbp = 0.d0
  do i = 1,3
     trbp(i,i) = 1.d0
  enddo

  inquire(unit=nfmatbp, opened = open_check)
  if(open_check) then
     rewind nfmatbp
     do i = 1, 3
        read(nfmatbp,*,end=1,err=1) trbp(i,1),trbp(i,2),trbp(i,3)
     enddo
     goto 2
1    trbp = 0.d0
     do i = 1,3
        trbp(i,i) = 1.d0
     enddo
2    continue
  end if

  call inver3n(3,trbp,trpb)

  mat1 = transpose(trpb)
  call inver3n(3,mat1,mat2)
  call matpr3(rltv,mat2,trmat)
!!$  trmat = matmul(rltv,mat2)

end subroutine get_trmat

subroutine read_konly(paramset,nfinp,nfout,knv3 &
     & ,altv,rltv,kv3,vkxyz,qwgt,str,len_str)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: read_konly
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

  use m_Const_Parameters, only : CARTS, CRDTYP, DP, NONAME, PAI2 &
       &                       , NODATA, BUCS
  implicit none
  logical, intent(in)        :: paramset
  integer, intent(in)        :: nfinp, nfout, knv3
  real(kind=DP), intent(in)  :: altv(3,3),rltv(3,3)
  integer, intent(out)       :: kv3
  real(kind=DP), intent(out) :: vkxyz(knv3,3,CRDTYP),qwgt(knv3)
  integer, intent(in)        :: len_str
  character(len=len_str)        :: str

  integer                    :: sw_k_coord_system
  integer                    :: ik,i,incunt,to,from

  real(kind=DP), pointer, dimension(:)   :: wka
  real(kind=DP), pointer, dimension(:,:) :: trmat

  rewind nfinp
  call down2wayksample  !-(contained here)
  read(nfinp,'(a132)') str
  call read_nkpnt(str,len_str,kv3) !-(b_Kpoints) ->kv3
  if(paramset) return

  allocate(wka(4));  allocate(trmat(3,3))

  read(nfinp,'(a132)') str
  call read_coordsystem(str,len_str,sw_k_coord_system)
  if(sw_k_coord_system == NODATA) then
     sw_k_coord_system = CARTS
  else
     read(nfinp,'(a132)') str
  end if
  do ik = 1, kv3
     call chnnm(str,len_str,4,wka,incunt)
     write(nfout,*) ' coordinates = ', wka(1),wka(2),wka(3),wka(4)
     if(incunt < 4) call phase_error_with_msg(nfout,' ! data shortage (nfinp, kpoints)',__LINE__,__FILE__)
     vkxyz(ik,1:3,sw_k_coord_system) = wka(1:3)
     qwgt(ik) = wka(4)
     if(ik < kv3) read(nfinp,'(a132)') str
  end do


  if(sw_k_coord_system == CARTS) then
     to = BUCS; from = CARTS; trmat = transpose(altv)/PAI2
  else
     to = CARTS; from = BUCS; trmat = rltv
  end if
  do ik = 1, kv3
     vkxyz(ik,1:3,to) = matmul(trmat,vkxyz(ik,1:3,from))
  end do

  do ik = 1, kv3
     write(nfout,'(i3,3f8.4,3x,3f8.4,3x,f8.4)') ik &
          &, (vkxyz(ik,i,CARTS),i=1,3) &
          &, (vkxyz(ik,i,BUCS) ,i=1,3),qwgt(ik)
  end do
  deallocate(trmat); deallocate(wka)
contains
  subroutine down2wayksample
    integer :: way_ksample_tmp

1   read(nfinp,'(a132)') str
    call ch4ksample(str,len_str,way_ksample_tmp) !-(b_Kpoints)
    if(way_ksample_tmp == NONAME) goto 1
  end subroutine down2wayksample
end subroutine read_konly

subroutine set_kfromtemporaryarray(nfout,ipri,nfmatbp,kv3_adj,kv3,kx,ky,kz,w &
     & ,rltv,vkxyz,qwgt)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: set_kfromtemporaryarray
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

  use m_Const_Parameters, only : CARTS, CRDTYP, DP, NONAME, PAI2 &
       &                       , NODATA, BUCS
  implicit none
  integer, intent(in) ::        nfout, ipri, nfmatbp,kv3,kv3_adj
  real(kind=DP), intent(in), dimension(kv3) :: kx,ky,kz,w
  real(kind=DP), intent(in)  :: rltv(3,3)
  real(kind=DP), intent(out) :: vkxyz(kv3_adj,3,CRDTYP),qwgt(kv3_adj)

  integer                    :: nn,nk
  real(kind=DP) :: kv(3)

  real(kind=DP), pointer, dimension(:,:) :: trmat,trbp,trpb,mat1,mat2 

  allocate(trmat(3,3))
  allocate(trbp(3,3))
  allocate(trpb(3,3))
  allocate(mat1(3,3))
  allocate(mat2(3,3))

!!$  call get_trmat  !-(contained here) ->(trmat)
  call get_trmat(rltv,trbp,trpb,mat1,mat2,trmat,nfmatbp)

  do nn = 1, kv3
     kv(1) = kx(nn); kv(2) = ky(nn); kv(3) = kz(nn)
     vkxyz(nn,1:3,CARTS) = matmul(trmat,kv)
     qwgt(nn) = w(nn)
  enddo

  deallocate(trmat)
  deallocate(trbp )
  deallocate(trpb )
  deallocate(mat1 )
  deallocate(mat2 )

  if(ipri >= 2) then
     write(nfout,'(/,"  << readk0 >>")')
     write(nfout,*) ' !Total Generated Kpoints = ',kv3
     do nk = 1, kv3
        write(nfout,'(i4," ",3f12.6)') nk,(vkxyz(nk,nn,CARTS),nn=1,3)
     end do
  end if
!!$contains
!!$  subroutine get_trmat
!!$!    make translation matrix  trpb (P -> B)
!!$    integer  :: i
!!$
!!$    rewind nfmatbp
!!$    do i = 1, 3
!!$       read(nfmatbp,*,end=1,err=1) trbp(i,1),trbp(i,2),trbp(i,3)
!!$    enddo
!!$    goto 2
!!$
!!$1   trbp = 0.d0
!!$    do i = 1,3
!!$       trbp(i,i) = 1.d0
!!$    enddo
!!$
!!$2   continue
!!$    call inver3n(3,trbp,trpb)
!!$
!!$    mat1 = transpose(trpb)
!!$    call inver3n(3,mat1,mat2)
!!$    call matpr3(rltv,mat2,trmat)
!!$
!!$  end subroutine get_trmat
end subroutine set_kfromtemporaryarray

subroutine kpmsf0(paramset,nfout,altv,rltv,ksm,knv3,kv3,vkxyz,qwgt)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: kpmsf0
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

  use m_Const_Parameters, only : BUCS, CARTS, CRDTYP, DP
  implicit none
  logical, intent(in)       :: paramset
  integer, intent(in)       :: nfout
  real(kind=DP), intent(in) :: altv(3,3), rltv(3,3)
  integer, intent(in)       :: ksm(3,2)
  integer, intent(in)       :: knv3
  integer, intent(out)      :: kv3
  real(kind=DP), intent(out):: vkxyz(knv3,3,CRDTYP)
  real(kind=DP), intent(out):: qwgt(knv3)

  real(kind=DP) :: fvy,fvz,dgvy,dgvz, fv2,fv3
  integer i2, i3

  if(.not.paramset) then
     if(ksm(1,1) >= 1) &
          & write(nfout,*) ' ksm(1,1) (=',ksm(1,1),') has no meaning here(kpmsf0).'
     if(ksm(1,2) >= 2) &
          & write(nfout,*) ' ksm(1,2) (=',ksm(1,2),') has no meaning here(kpmsf0).'
  end if
  fvy = ksm(2,1)*2.d0
  fvz = ksm(3,1)*2.d0
  kv3 = 0
!!$  write(nfout,*) ' !D nky, nkz = ', ksm(2,1),ksm(3,1) !nky, nkz
  dgvy=0.5d0*(rltv(2,2)/fvy+rltv(2,3)/fvz)
  dgvz=0.5d0*(rltv(3,2)/fvy+rltv(3,3)/fvz)
  do i2 = 1,ksm(2,2) ! nky2
     fv2 = (i2-ksm(2,1)-1)/fvy
     do i3 = 1,ksm(3,2) ! nkz2
        fv3 = (i3-ksm(3,1)-1)/fvz
        kv3 = kv3 + 1
!----*----*   for parameter setting 5th June '94 by Y.M.
        if(.not.paramset .and. kv3 <= knv3) then
           vkxyz(kv3,1,CARTS) = 0.D0
           vkxyz(kv3,2,CARTS) = rltv(2,2)*fv2 + rltv(2,3)*fv3 + dgvy
           vkxyz(kv3,3,CARTS) = rltv(3,2)*fv2 + rltv(3,3)*fv3 + dgvz
        end if
     end do
  end do
  if(.not.paramset) then
     qwgt(1:kv3) = 1.d0/kv3
     write(nfout,*) 'number of generated k-points=',kv3
     if(kv3 /= knv3) write(nfout,*) '**warn mm should be knv3**'
  end if

end subroutine kpmsf0

subroutine kpmwbz0(paramset,nfout,altv,rltv,ksm,knv3,kv3,vkxyz,qwgt)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: kpmwbz0
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

  use m_Const_Parameters, only : BUCS, CARTS, CRDTYP, DP
  implicit none
  logical, intent(in)       :: paramset
  integer, intent(in)       :: nfout
  real(kind=DP), intent(in) :: altv(3,3), rltv(3,3)
  integer, intent(in)       :: ksm(3,2)
  integer, intent(in)       :: knv3
  integer, intent(out)      :: kv3
  real(kind=DP), intent(out):: vkxyz(knv3,3,CRDTYP)
  real(kind=DP), intent(out):: qwgt(knv3)

  real(kind=DP)             :: dgv(3),fv123(3),dfv(3)
  integer                   :: i1,i2,i3

  if(ksm(1,1) == 0 &
       & .or. (     ksm(1,2) == 0 &
       &      .and. ksm(2,2) == 0 &
       &      .and. ksm(3,2) == 0)&
       ) then
     kv3 = 1
     if(paramset) return
     vkxyz(1,1:3,CARTS) = 0.d0
     return
  endif

  if(ksm(1,2) == 0) call phase_error_with_msg(nfout," ksm(1,2) == 0",__LINE__,__FILE__)
  if(ksm(2,2) == 0) call phase_error_with_msg(nfout," ksm(2,2) == 0",__LINE__,__FILE__)
  if(ksm(3,2) == 0)  call phase_error_with_msg(nfout," ksm(3,2) == 0",__LINE__,__FILE__)

  dfv(1:3) = 1.d0/(ksm(1:3,1)*2)
  kv3 = 0
  dgv = 0.5d0 * matmul(rltv,dfv)

  do i1 = 1, ksm(1,2)
     fv123(1) = (i1-ksm(1,1)-1)*dfv(1)
     do i2 = 1, ksm(2,2)
        fv123(2) = (i2-ksm(2,1)-1)*dfv(2)
        do i3 = 1, ksm(3,2)
           fv123(3) = (i3-ksm(3,1)-1)*dfv(3)
           kv3 = kv3 + 1
           if(.not.paramset .and. kv3 <= knv3) then
              vkxyz(kv3,1:3,CARTS) = matmul(rltv,fv123) + dgv
           end if
        end do
     end do
  end do

  if(.not.paramset) then
     qwgt(1:kv3) = 1.d0/kv3
  end if

  write(nfout,*) 'number of generated k-points=',kv3
  if(kv3 /= knv3) write(nfout,*) '**warning mm should be knv3**'
end subroutine kpmwbz0

subroutine sccm0(paramset,nfout,altv,rltv,ksm,knv3,kv3,vkxyz,qwgt)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: sccm0
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

  use m_Const_Parameters, only : BUCS, CARTS, CRDTYP, DP, PAI2
  implicit none
  logical, intent(in)       :: paramset
  integer, intent(in)       :: nfout
  real(kind=DP), intent(in) :: altv(3,3), rltv(3,3)
  integer, intent(in)       :: ksm(3,2)
  integer, intent(in)       :: knv3
  integer, intent(out)      :: kv3
  real(kind=DP), intent(out):: vkxyz(knv3,3,CRDTYP)
  real(kind=DP), intent(out):: qwgt(knv3)

  real(kind=DP)             :: h0, swgt, wgt, q(3), b
  integer                   :: ikx, iky, ikz, iq

  if(.not.paramset) b = PAI2/(altv(1,1)*2)

  h0   = 1.0d0/ksm(1,1)
  swgt = 0.0d0
  kv3  = 0

  do ikx = 0, ksm(1,1) !nkx
     q(1) = ikx*h0
     do iky = 0, ikx
        q(2) = iky*h0
        do ikz = 0,iky
           q(3) = ikz*h0
           iq = 1
           if(ikz == 0) then
              if(ikx == iky) then
                 if(ikx == 0) then
                    iq=48
                 else if(ikx == ksm(1,1)) then
                    iq=16
                 else
                    iq=4
                 end if
              else if(iky.eq.0) then
                 if(ikx.eq.ksm(1,1)) then
                    iq=16
                 else
                    iq=8
                 end if
              else if(ikx == ksm(1,1)) then
                 iq=4
              else
                 iq=2
              end if
           else if(ikx == iky) then
              if(ikz == ikx) then
                 if(ikz == ksm(1,1)) then
                    iq=48
                 else
                    iq=6
                 end if
              else if(ikx.eq.ksm(1,1)) then
                 iq=8
              else
                 iq=2
              end if
           else if(iky == ikz) then
              if(ikx == ksm(1,1)) then
                 iq=4
              else
                 iq=2
              end if
           else if(ikx == ksm(1,1)) then
              iq=2
           else
              iq=1
           end if
           wgt=48.0d0/(8.d0*iq*ksm(1,1)**3)
           swgt=swgt+wgt
           kv3 = kv3 + 1
!----*----*   for parameter setting 5th June '94 by Y.M.
           if(.not.paramset .and. kv3 <= knv3) then
              vkxyz(kv3,1:3,CARTS) = q(1:3)*b
              qwgt (kv3)           = wgt
           end if
!----*----*
        end do
     end do
  end do
  write(nfout,201) kv3,swgt
  201 format('  nktot =',i5/'  swgt  =',f20.16)
end subroutine sccm0

subroutine bccm0(paramset,nfout,altv,rltv,ksm,knv3,kv3,vkxyz,qwgt)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: bccm0
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

  use m_Const_Parameters, only : BUCS, CARTS, CRDTYP, DP, PAI2
  implicit none
  logical, intent(in)       :: paramset
  integer, intent(in)       :: nfout
  real(kind=DP), intent(in) :: altv(3,3), rltv(3,3)
  integer, intent(in)       :: ksm(3,2)
  integer, intent(in)       :: knv3
  integer, intent(out)      :: kv3
  real(kind=DP), intent(out):: vkxyz(knv3,3,CRDTYP)
  real(kind=DP), intent(out):: qwgt(knv3)

  integer      :: nkxh, nkx, ikx,iky,ikz,iq
  real(kind=DP) :: h0, swgt, wgt, q(3), b

  if(.not.paramset) b = PAI2/(altv(1,1)*2)

  nkxh = ksm(1,1)/2
  nkx  = nkxh+nkxh
  h0=1.0d0/nkx
  swgt=0.0d0
  KV3=0
  Loop_ikx :do ikx = 0, nkx
     q(1) = ikx*h0
     do iky = 0,ikx
        if(ikx+iky > nkx) cycle Loop_ikx
        q(2) = iky*h0
        do ikz = 0,iky
           q(3) = ikz*h0
           iq = 1
           if(ikz == 0) then
              if(ikx == iky) then
                 if(ikx == 0) then
                    iq=48
                 else if(ikx == nkxh) then
                    iq=8
                 else
                    iq=4
                 end if
              else if(iky == 0) then
                 if(ikx == nkx) then
                    iq=48
                 else
                    iq=8
                 end if
              else if(ikx+iky == nkx) then
                 iq=4
              else
                 iq=2
              end if
           else if(ikx == iky) then
              if(ikz == ikx) then
                 if(ikz == nkxh) then
                    iq=24
                 else
                    iq=6
                 end if
              else if(ikx == nkxh) then
                 iq=4
              else
                 iq=2
              end if
           else if(iky == ikz) then
              if(ikx+iky == nkx) then
                 iq=6
              else
                 iq=2
              end if
           else if(ikx+iky == nkx) then
              iq=2
           else
              iq=1
           end if
           wgt=48.0d0/dble(2*iq*nkx*nkx*nkx)
           swgt=swgt+wgt
           kv3=kv3+1
!----*----*   for parameter setting 5th June '94 by Y.M.
           if(.not.paramset .and. kv3 <= knv3) then
              Vkxyz(1:3,KV3,CARTS)=Q(1:3)*b
              qwgt(kv3)=wgt
            end if
!----*----*
         end do
      end do
  end do Loop_ikx
  write(nfout,201) kv3,swgt
201 format('  kv3 =',i5/'  swgt  =',d20.13)
end subroutine bccm0

subroutine fccm0(paramset,nfout,altv,rltv,ksm,knv3,kv3,vkxyz,qwgt)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: fccm0
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

  use m_Const_Parameters, only : BUCS, CARTS, CRDTYP, DP, PAI2
  implicit none
  logical, intent(in)       :: paramset
  integer, intent(in)       :: nfout
  real(kind=DP), intent(in) :: altv(3,3), rltv(3,3)
  integer, intent(in)       :: ksm(3,2)
  integer, intent(in)       :: knv3
  integer, intent(out)      :: kv3
  real(kind=DP), intent(out):: vkxyz(knv3,3,CRDTYP)
  real(kind=DP), intent(out):: qwgt(knv3)

  real(kind=DP) :: q(3),qw, h0, swgt, wgt, b
  integer       :: nkxh, nkxq, nkxh3, nkxq3, nkx,ikx,iky,ikz

!!$  print *,' altv(1,1) = ', altv(1,1)
  h0 = (altv(1,1) + altv(2,1) + altv(3,1))*0.5
  if(.not.paramset) b = PAI2/(h0*2)
  NKXH=ksm(1,1)/2
  NKXQ=NKXH/2
  if(nkxq == 0) then
     print '(" ksm(1,1) = ",i4," nkxh = ",i4," nkxq = ",i4)', ksm(1,1),nkxh,nkxq
     call phase_error_with_msg(nfout,' !! ksm(1,1) should be larger or equal to 4',__LINE__,__FILE__)
  end if
  NKXH=NKXQ+NKXQ
  NKX=NKXH+NKXH
  NKXH3=NKXH+NKXH+NKXH
  NKXQ3=NKXQ+NKXQ+NKXQ
  print *,' nkxh = ', nkxh, ' nkx = ', nkx,' nkxh3 = ', nkxh3, ' nkxq3 = ', nkxq3
  H0=1.0D0/DBLE(NKX)
  SWGT=0.0D0


  KV3=0
  do ikx = 0,nkx
     q(1) = ikx*h0
     Loop_iky :do iky = 0,ikx
        q(2) = iky*h0
        do ikz = 0,iky
           if(ikx+iky+ikz.gt.nkxh3) cycle Loop_iky
           q(3) = ikz*h0
           qw = 1.d0
           if(ikz == 0) then
              if(iky == 0) then
                 if(ikx == 0) then
                    qw=1.0d0
                 else if(ikx == nkx) then
                    qw=3.d0
                 else
                    qw=6.0d0
                 end if
              else if(ikx == iky) then
                 if(ikx == nkxq3) then
                    qw=3.5d0
                 else
                    qw=12.0d0
                 end if
              else if(ikx == nkx) then
                 if (iky == nkxh) then
                    qw=6.0d0
                 else
                    qw=12.0d0
                 end if
              else if(ikx+iky == nkxh3) then
                 qw=7.0d0
              else
                 qw=24.0d0
              end if
           else if(ikx == iky) then
              if(ikz == ikx) then
                 if(ikz == nkxh) then
                    qw=4.0d0
                 else
                    qw=8.0d0
                 end if
              else if(ikx+ikx+ikz == nkxh3) then
                 qw=12.0d0
              else
                 qw=24.0d0
              end if
           else if(ikx == nkx) then
              if(iky == ikz) then
                 if(iky == nkxq) then
                    qw=8.5d0
                 else
                    qw=12.0d0
                 end if
              else if(iky+ikz == nkxh) then
                 qw=17.0d0
              else
                 qw=24.0d0
              end if
           else if(iky == ikz) then
              if(ikx+iky+iky == nkxh3) then
                 qw=12.0d0
              else
                 qw=24.0d0
              end if
           else if(ikx+iky+ikz == nkxh3) then
              qw=24.0d0
           else
              qw=48.0d0
           end if
           wgt=qw/dble(4*nkx*nkx*nkx)
           swgt=swgt+wgt
           kv3=kv3+1
!----*----*   for parameter setting 5th June '94 by Y.M.
           if(.not.paramset .and. kv3 <= knv3) then
              vkxyz(kv3,1:3,CARTS) = q(1:3)*b
              qwgt(kv3)=wgt
           end if
!----*----*
        end do
     end do Loop_iky
  end do
  write(nfout,201) kv3,swgt
201 FORMAT('  KV3 =',I5/'  SWGT  =',D20.13)
end subroutine fccm0

subroutine hexm0(paramset,nfout,altv,rltv,ksm,knv3,kv3,vkxyz,qwgt)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: hexm0
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

  use m_Const_Parameters, only : BUCS, CARTS, CRDTYP, DP
  implicit none
  logical, intent(in)       :: paramset
  integer, intent(in)       :: nfout
  real(kind=DP), intent(in) :: altv(3,3), rltv(3,3)
  integer, intent(in)       :: ksm(3,2)
  integer, intent(in)       :: knv3
  integer, intent(out)      :: kv3
  real(kind=DP), intent(out):: vkxyz(knv3,3,CRDTYP)
  real(kind=DP), intent(out):: qwgt(knv3)

  real(kind=DP) :: q(3), hx, hy, wgt, swgt
  integer       :: nkxh, nkytth, nkyh, ikx, iky, ikz, iq

  if((mod(ksm(1,1),2) /= 0) .or. (mod(ksm(2,1),6) /= 0)) then
     write(nfout,*) ' ksm(1,1), ksm(2,1) err , ksm(1:2,1) =' &
          &,ksm(1,1),ksm(2,1)
     call phase_error_with_msg(nfout,'ksm(1,1), ksm(2,1) err',__LINE__,__FILE__)
  end if
  nkxh   = ksm(1,1)/2
  nkytth = (ksm(2,1)/3)*2
  nkyh   =  ksm(2,1)/2
  hx=1.0d0/ksm(1,1)
  hy=1.0d0/ksm(2,1)
  swgt=0.0d0
  kv3=0
  do ikx = 0, nkxh
     q(1) = ikx*hx
     do iky = 0,nkytth
        q(2) = iky*hy
        do ikz = 0,iky
           if((iky > ikz*2).or.(iky > ksm(2,1)-ikz)) cycle
           q(3) = ikz*hy
           iq = 1
           if((ikx == 0).or.(ikx == nkxh)) then
              if(iky == ikz) then
                 if(iky == 0) then
                    iq=24
                 else if(iky == nkyh) then
                    iq=8
                 else
                    iq=4
                 end if
              else if(iky == 2*ikz) then
                 if(iky == nkytth) then
                    iq=12
                 else
                    iq=4
                 end if
              else if(iky == ksm(2,1)-ikz) then
                 iq=4
              else
                 iq=2
              end if
           else
              if(iky == ikz) then
                 if(iky == 0) then
                    iq=12
                 else if(iky == nkyh) then
                    iq=4
                 else
                    iq=2
                 end if
              else if(iky == 2*ikz) then
                 if(iky == nkytth) then
                    iq=6
                 else
                    iq=2
                 end if
              else if(iky == ksm(2,1)-ikz) then
                 iq=2
              else
                 iq=1
              end if
           end if
!!$           wgt=24.0d0/(iq*ksm(1,1)*ksm(2,1)*ksm(2,1))
           wgt = 24.d0/(iq*product(ksm(1:3,1)))
           swgt=swgt+wgt
           kv3=kv3+1
!----*----*   for parameter setting 5th June '94 by Y.M.
           if(.not.paramset  .and. kv3 <= knv3) then
              vkxyz(kv3,1:3,CARTS) = matmul(rltv,q)
              qwgt(kv3) = wgt
           end if
!----*----*
        end do
     end do
  end do
  write(nfout,201) kv3,swgt
201 format('  kv3 =',i5/'  swgt  =',f20.16)
end subroutine hexm0

subroutine read_nkpnt(str,len_str,ikpnt)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: read_nkpnt
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

  implicit none
  integer,             intent(in)  :: len_str
  character(len=len_str), intent(in)  :: str
  integer,             intent(out) :: ikpnt

  logical :: tf
  integer :: ic = 1

  call strncmp2(str,len_str,'nkpoints',8,tf)
  if(tf) then
     print *,' str(ic:len_str) = ',str(ic:len_str)
     call skip2next(str,len_str,'n',ic)
     print *,' str(ic:len_str) = ',str(ic:len_str)
     read(str(ic:len_str),*) ikpnt
  else
     call phase_error_with_msg(6,' invalid input of nkpoints',__LINE__,__FILE__)
  end if
end subroutine read_nkpnt

subroutine ch4ksample(str,len_str,way_ksample)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: ch4ksample
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

  use m_Const_Parameters, only : NONAME,MESH, SKPS_DIRECT_IN &
       &, GAMMA, FILE, MONKHORST_PACK
  implicit none
  character(len=*),  intent(in) :: str
  integer,        intent(in) :: len_str
  integer,        intent(out):: way_ksample

  logical :: tf
  integer :: ic

  call strncmp2(str,len_str,'way_ksample',11,tf)
  if(tf) then
     call skip2next(str,len_str,'W',ic)
     if(str(ic:ic) == 'm' .or. str(ic:ic) == 'M') then
         if(str(ic+1:ic+1) == 'o' .or. str(ic+1:ic+1) == 'O') then
            way_ksample = MONKHORST_PACK
         else if(str(ic+1:ic+1) == 'e' .or. str(ic+1:ic+1) == 'E') then
            way_ksample = MESH
         end if
     else if(str(ic:ic) == 's' .or. str(ic:ic) == 'S') then
        way_ksample = SKPS_DIRECT_IN
     else if(str(ic:ic) == 'g' .or. str(ic:ic) == 'G') then
        way_ksample = GAMMA
     else if(str(ic:ic) == 'f' .or. str(ic:ic) == 'F') then
        way_ksample = FILE
     endif
  else
     way_ksample = NONAME
  endif
end subroutine ch4ksample

subroutine read_coordsystem(str,len_str,sw_k_coord_system)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: read_coordsystem
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

  use m_Const_Parameters, only : NODATA, CARTS, BUCS
  implicit none
  integer,         intent(in)  :: len_str
  character(len=len_str), intent(in)  :: str
  integer,         intent(out) :: sw_k_coord_system

  logical      :: tf

  sw_k_coord_system = NODATA
  call strncmp2(str,len_str,'k_coord_system',14,tf)
  if(tf) then
     call strncmp2(str,len_str,'pucv',4,tf)
     if(tf) then
        sw_k_coord_system = BUCS
        goto 1001
     endif
     call strncmp2(str,len_str,'bucs',4,tf)
     if(tf) then
        sw_k_coord_system = BUCS
        goto 1001
     end if
     call strncmp2(str,len_str,'cart',4,tf)
     if(tf) then
        sw_k_coord_system = CARTS
     end if
  end if
1001 continue
end subroutine read_coordsystem


subroutine gnrt_k0(nbztyp_spg,altv,nx,ny,nz &
                   & ,paramset,nfout,ipri,knv3 &
                   & ,rltv,kv3,vkxyz,qwgt,nfkpgn,nfspg,ipri_kp,ipri_spg, itrs )
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: gnrt_k0
!
!  AUTHOR(S): H. Mizouchi   February/25/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

  use m_Const_Parameters, only : CARTS, CRDTYP, DP &
       &                       , WHOLE_BZ, SIMPLE_CUBIC &
       &                       , BCC, FCC, DIAMOND, HEXAGONAL &
       &                       , GENERAL, GENERAL_LARGER 
  implicit none
  logical, intent(in)        :: paramset
  integer, intent(in)        :: nfout, ipri, knv3 &
       &                      , nfkpgn,nfspg,ipri_kp,ipri_spg, itrs 
  real(kind=DP), intent(in)  :: rltv(3,3),altv(3,3)
  integer, intent(out)       :: kv3
  real(kind=DP), intent(out) :: vkxyz(knv3,3,CRDTYP),qwgt(knv3)

  integer                    :: nn,nk,nw,nwei,nv(3)

  real(kind=DP), pointer, dimension(:,:) :: trmat,trbp,trpb,mat1,mat2 
  integer   :: np2,np1,np0
  integer   :: nxx0,nyy0,nzz0,nxx,nyy,nzz
  integer   :: nbztyp_spg,nx,ny,nz
  integer :: lmnp0, lmnp1, lmnp2
!!$  integer :: idim, il, ngen, inv
  integer :: il
  integer :: nx1, ny1, nz1, nd
  real(kind=DP), pointer, dimension(:,:) :: pa0_wk,pb0_wk,pb_wk 
  integer,       pointer, dimension(:,:) :: ka0_wk,ka2_wk
  integer,       pointer, dimension(:)   :: ip10_wk,ip02_wk,ip12_wk &
   &                                       ,ip01_wk,ip21_wk,iu21_wk,iv21_wk,nstar2_wk &
   &                                       ,ip20_wk

  if(nbztyp_spg == GENERAL .or.nbztyp_spg == GENERAL_LARGER) then
      rewind nfkpgn
      read(nfkpgn,*) nxx0,nyy0,nzz0
      write(6,*) 'nxx0,nyy0,nzz0 desu', nxx0,nyy0,nzz0

      rewind nfspg
      read(nfspg,*) nn
      read(nfspg,*) 
!!$      read(nfspg,*) idim, il, ngen, inv
      read(nfspg,*) nw, il, nk, nn

      call nskma0(il,nxx0,nyy0,nzz0,nxx,nyy,nzz,nx1,ny1,nz1,nd)         
      lmnp0=(nxx+1)*(nyy+1)*(nzz+1)
      lmnp1=lmnp0
      lmnp2=lmnp0

      allocate(ip10_wk(lmnp0))  ; ip10_wk = 0        
      allocate(ip20_wk(lmnp0))  ; ip20_wk = 0        
      allocate(ip01_wk(lmnp1))    ; ip01_wk = 0
      allocate(ip02_wk(lmnp2))    ; ip02_wk = 0
      allocate(ip21_wk(lmnp1))    ; ip21_wk = 0
      allocate(ip12_wk(lmnp2))    ; ip12_wk = 0
      allocate(iu21_wk(lmnp1))    ; iu21_wk = 0
      allocate(iv21_wk(lmnp1))    ; iv21_wk = 0
      allocate(nstar2_wk(lmnp2))  ; nstar2_wk = 0
      allocate(pa0_wk(3,lmnp0))  ; pa0_wk = 0
      allocate(pb0_wk(3,lmnp0))  ; pb0_wk = 0
      allocate(pb_wk(3,lmnp2))  ; pb_wk = 0
      allocate(ka0_wk(4,lmnp0))  ; ka0_wk = 0
      allocate(ka2_wk(4,lmnp2))  ; ka2_wk = 0


      call setkp0(np2,np1,np0,lmnp0,lmnp1,lmnp2 &
                & ,nxx0,nyy0,nzz0,nxx,nyy,nzz &
                & ,ip10_wk,ip20_wk,ip01_wk,ip02_wk,ip21_wk,ip12_wk &
                & ,iu21_wk,iv21_wk &
                & ,nstar2_wk,pa0_wk,pb0_wk,pb_wk,ka0_wk,ka2_wk &
                & ,nfspg,ipri_kp,ipri_spg, itrs )

  else

      if(nbztyp_spg.eq.SIMPLE_CUBIC) il = 1
      if(nbztyp_spg.eq.BCC) il = 3
      if(nbztyp_spg.eq.FCC) il = 2
      if(nbztyp_spg.eq.DIAMOND) il = 2
      if(nbztyp_spg.eq.HEXAGONAL) il = 0

      call nskma0(il,nx,ny,nz,nxx,nyy,nzz,nx1,ny1,nz1,nd)

      lmnp0=(nxx+1)*(nyy+1)*(nzz+1)
      lmnp1=lmnp0
      lmnp2=lmnp0

      allocate(ip10_wk(lmnp0))  ; ip10_wk = 0        
      allocate(ip20_wk(lmnp0))  ; ip20_wk = 0        
      allocate(ip01_wk(lmnp1))    ; ip01_wk = 0
      allocate(ip02_wk(lmnp2))    ; ip02_wk = 0
      allocate(ip21_wk(lmnp1))    ; ip21_wk = 0
      allocate(ip12_wk(lmnp2))    ; ip12_wk = 0
      allocate(iu21_wk(lmnp1))    ; iu21_wk = 0
      allocate(iv21_wk(lmnp1))    ; iv21_wk = 0
      allocate(nstar2_wk(lmnp2))  ; nstar2_wk = 0
      allocate(pa0_wk(3,lmnp0))  ; pa0_wk = 0
      allocate(pb0_wk(3,lmnp0))  ; pb0_wk = 0
      allocate(pb_wk(3,lmnp2))  ; pb_wk = 0
      allocate(ka0_wk(4,lmnp0))  ; ka0_wk = 0
      allocate(ka2_wk(4,lmnp2))  ; ka2_wk = 0

      call setkp0_default(nbztyp_spg,altv,nx,ny,nz &
                  & ,np2,np1,np0,lmnp0,lmnp1,lmnp2 &
                  & ,nxx,nyy,nzz &
                  & ,ip10_wk,ip20_wk,ip01_wk,ip02_wk,ip21_wk,ip12_wk &
                  & ,iu21_wk,iv21_wk &
                  & ,nstar2_wk,pa0_wk,pb0_wk,pb_wk,ka0_wk,ka2_wk &
                  & ,ipri_kp,ipri_spg, itrs )

  end if

  kv3 = np2
  nwei = np1


  if(paramset) then
     deallocate(ip10_wk)
     deallocate(ip20_wk)
     deallocate(ip01_wk)
     deallocate(ip02_wk)
     deallocate(ip21_wk)
     deallocate(ip12_wk)
     deallocate(iu21_wk)
     deallocate(iv21_wk)
     deallocate(nstar2_wk)
     deallocate(pa0_wk)
     deallocate(pb0_wk)
     deallocate(pb_wk)
     deallocate(ka0_wk)
     deallocate(ka2_wk)
     return
  endif

  allocate(trmat(3,3))
  allocate(trbp(3,3))
  allocate(trpb(3,3))
  allocate(mat1(3,3))
  allocate(mat2(3,3))

  call get_trmat1  !-(contained here) ->(trmat)

  do nn = 1, kv3
     nv(1) = ka0_wk(1,ip02_wk(nn))
     nv(2) = ka0_wk(2,ip02_wk(nn))
     nv(3) = ka0_wk(3,ip02_wk(nn))
     nk = ka0_wk(4,ip02_wk(nn))
     nw = nstar2_wk(nn)
     vkxyz(nn,1:3,CARTS) = matmul(trmat,nv)/dble(nk)
     qwgt(nn) = dble(nw)/dble(nwei)
  enddo

  deallocate(trmat)
  deallocate(trbp )
  deallocate(trpb )
  deallocate(mat1 )
  deallocate(mat2 )


  deallocate(ip10_wk)
  deallocate(ip20_wk)
  deallocate(ip01_wk)
  deallocate(ip02_wk)
  deallocate(ip21_wk)
  deallocate(ip12_wk)
  deallocate(iu21_wk)
  deallocate(iv21_wk)
  deallocate(nstar2_wk)
  deallocate(pa0_wk)
  deallocate(pb0_wk)
  deallocate(pb_wk)
  deallocate(ka0_wk)
  deallocate(ka2_wk)



  if(ipri >= 2) then
     write(nfout,'(/,"  << gnrt_k0 >>")')
     write(nfout,*) ' !Total Generated Kpoints = ',kv3
     do nk = 1, kv3
        write(nfout,'(i4," ",3f12.6)') nk,(vkxyz(nk,nn,CARTS),nn=1,3)
     end do
  end if
contains
  subroutine get_trmat1
!    make translation matrix  trpb (P -> B)
    integer  :: i
    real(DP)           :: tab(3,3),ta1(3,48)
    integer            :: ng1,lra1(3,3,48)

!! setspg  !!

  if(nbztyp_spg == 100 .or.nbztyp_spg == 101) then
      call setspg(tab,ng1,ta1,lra1,nfspg,ipri_spg)

  else
     call setspg_default(nbztyp_spg,altv,tab,ng1,ta1,lra1,ipri_spg)
  end if

!     make translation matrix  trpb (P -> B)
  do i = 1,3
     trbp(i,1) = tab(i,1)
     trbp(i,2) = tab(i,2)
     trbp(i,3) = tab(i,3)
  end do
    goto 2

1   trbp = 0.d0
    do i = 1,3
       trbp(i,i) = 1.d0
    enddo

2   continue
    call inver3n(3,trbp,trpb)

    mat1 = transpose(trpb)
    call inver3n(3,mat1,mat2)
    call matpr3(rltv,mat2,trmat)

  end subroutine get_trmat1
end subroutine gnrt_k0


! -- subroutine gnrt_k0_n ---------------------

! ================================== modified by K. Tagami ================ 12.0A
!subroutine gnrt_k0_n(il,ngen,inv,igen,jgen,a,b,c,ca,cb,cc &
!      &              , nbztyp_spg,nx,ny,nz &
!      &              , paramset,nfout,ipri,knv3,rltv &
!      &              , kv3,vkxyz,qwgt,ipri_kp,trmat_out)
subroutine gnrt_k0_n(il,ngen,inv,igen,jgen,a,b,c,ca,cb,cc &
      &              , nbztyp_spg,nx,ny,nz &
      &              , paramset,nfout,ipri,knv3,rltv &
      &              , kv3,vkxyz,qwgt,ipri_kp,trmat_out, &
      &                gen_tetramesh_mode, use_altv_rltv, altv, itrs, &
      &                gen_name_in_carts, &
      &                knv3_fbz, kv3_fbz, vkxyz_fbz, to_ibz_from_fbz_for_kpoint )
! =========================================================================== 12.0A

!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: gnrt_k0_n
!
!  AUTHOR(S): T. Yamasaki   May/31/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

  use m_Const_Parameters, only : CARTS, CRDTYP, DP &
       &                       , WHOLE_BZ, SIMPLE_CUBIC &
       &                       , BCC, FCC, DIAMOND, HEXAGONAL &
       &                       , GENERAL, GENERAL_LARGER 
  implicit none
  integer, intent(in) ::        il,ngen,inv
  integer, intent(in) ::        igen(ngen),jgen(2,3,ngen)
  real(kind=DP), intent(in) ::  a,b,c,ca,cb,cc
  integer, intent(in) ::        nbztyp_spg
  integer, intent(in) ::        nx,ny,nz
  real(kind=DP), intent(in)  :: rltv(3,3)
  logical, intent(in) ::        paramset
  integer, intent(in) ::        nfout, ipri, knv3 &
       &                      , ipri_kp
  integer, intent(out)       :: kv3
  real(kind=DP), intent(out) :: vkxyz(knv3,3,CRDTYP),qwgt(knv3),trmat_out(3,3)

! ======================================== added by K. Tagami =========== 12.0A
  integer, intent(in) :: gen_tetramesh_mode
  integer, intent(in) :: itrs
  logical, intent(in) :: use_altv_rltv, gen_name_in_carts
  real(kind=DP), intent(in) :: altv(3,3)
! ======================================================================== 12.0A

! === KT_add === 2014/09/30
  integer, intent(in)        :: knv3_fbz
  integer, intent(out)       :: kv3_fbz
  real(kind=DP), intent(out) :: vkxyz_fbz(knv3_fbz,3,CRDTYP)
  integer, intent(out) :: to_ibz_from_fbz_for_kpoint(knv3_fbz)
! ============== 2014/09/30

  integer                    :: nn,nk,nw,nwei,nv(3),ix,iy,iz, iix,iiy,iiz,npx,npy
  integer :: icub, ni, ip0, kx,ky,kz
  integer, dimension(8) :: ip8

  real(kind=DP), pointer, dimension(:,:) :: trmat,trbp,trpb,mat1,mat2 
  integer :: np2,np1,np0
  integer :: nxx0,nyy0,nzz0,nxx,nyy,nzz
  integer :: lmnp0, lmnp1, lmnp2
  integer :: nx1, ny1, nz1, nd
  real(kind=DP), pointer, dimension(:,:) :: pa0_wk,pb0_wk,pb_wk 
  integer,       pointer, dimension(:,:) :: ka0_wk,ka2_wk
  integer,       pointer, dimension(:)   :: ip10_wk,ip02_wk,ip12_wk &
   &                                       ,ip01_wk,ip21_wk,iu21_wk,iv21_wk,nstar2_wk &
   &                                       ,ip20_wk

  if(ipri_kp >=2 ) write(nfout,'(" <<gnrt_k0_n>>")')

  if(nbztyp_spg == GENERAL .or.nbztyp_spg == GENERAL_LARGER) then
     nxx0 = nx; nyy0 = ny; nzz0 = nz 
     write(6,*) 'nxx0,nyy0,nzz0 desu', nxx0,nyy0,nzz0

     call nskma0(il,nxx0,nyy0,nzz0,nxx,nyy,nzz,nx1,ny1,nz1,nd)         
     lmnp0=(nxx+1)*(nyy+1)*(nzz+1)
     lmnp1=lmnp0
     lmnp2=lmnp0

     allocate(ip10_wk(lmnp0))  ; ip10_wk = 0        
     allocate(ip20_wk(lmnp0))  ; ip20_wk = 0        
     allocate(ip01_wk(lmnp1))    ; ip01_wk = 0
     allocate(ip02_wk(lmnp2))    ; ip02_wk = 0
     allocate(ip21_wk(lmnp1))    ; ip21_wk = 0
     allocate(ip12_wk(lmnp2))    ; ip12_wk = 0
     allocate(iu21_wk(lmnp1))    ; iu21_wk = 0
     allocate(iv21_wk(lmnp1))    ; iv21_wk = 0
     allocate(nstar2_wk(lmnp2))  ; nstar2_wk = 0
     allocate(pa0_wk(3,lmnp0))  ; pa0_wk = 0
     allocate(pb0_wk(3,lmnp0))  ; pb0_wk = 0
     allocate(pb_wk(3,lmnp2))  ; pb_wk = 0
     allocate(ka0_wk(4,lmnp0))  ; ka0_wk = 0
     allocate(ka2_wk(4,lmnp2))  ; ka2_wk = 0

! === KT_mod === 2015/02/16
!     call setkp0_n(il,ngen,inv,igen,jgen,a,b,c,ca,cb,cc  &
!          &          ,np2,np1,np0,lmnp0,lmnp1,lmnp2 &
!          &          ,nxx0,nyy0,nzz0,nxx,nyy,nzz &
!          &          ,ip10_wk,ip20_wk,ip01_wk,ip02_wk,ip21_wk,ip12_wk &
!          &          ,iu21_wk,iv21_wk &
!          &          ,nstar2_wk,pa0_wk,pb0_wk,pb_wk,ka0_wk,ka2_wk &
!          &          ,ipri_kp)
!
     if ( gen_tetramesh_mode == 0 ) then
        call setkp0_n(il,ngen,inv,igen,jgen,a,b,c,ca,cb,cc  &
             &          ,np2,np1,np0,lmnp0,lmnp1,lmnp2 &
             &          ,nxx0,nyy0,nzz0,nxx,nyy,nzz &
             &          ,ip10_wk,ip20_wk,ip01_wk,ip02_wk,ip21_wk,ip12_wk &
             &          ,iu21_wk,iv21_wk &
             &          ,nstar2_wk,pa0_wk,pb0_wk,pb_wk,ka0_wk,ka2_wk &
             &          ,ipri_kp, itrs )
     else if ( gen_tetramesh_mode == 1 ) then
        call setkp0_n_kt(il,ngen,inv,igen,jgen,a,b,c,ca,cb,cc  &
             &          ,np2,np1,np0,lmnp0,lmnp1,lmnp2 &
             &          ,nxx0,nyy0,nzz0,nxx,nyy,nzz &
             &          ,ip10_wk,ip20_wk,ip01_wk,ip02_wk,ip21_wk,ip12_wk &
             &          ,iu21_wk,iv21_wk &
             &          ,nstar2_wk,pa0_wk,pb0_wk,pb_wk,ka0_wk,ka2_wk &
             &          ,ipri_kp, &
             &           use_altv_rltv, altv, rltv, itrs, &
             &           gen_name_in_carts )
     endif

  else

     call nskma0(il,nx,ny,nz,nxx,nyy,nzz,nx1,ny1,nz1,nd)

     lmnp0=(nxx+1)*(nyy+1)*(nzz+1)
     lmnp1=lmnp0
     lmnp2=lmnp0

     allocate(ip10_wk(lmnp0))  ; ip10_wk = 0        
     allocate(ip20_wk(lmnp0))  ; ip20_wk = 0        
     allocate(ip01_wk(lmnp1))    ; ip01_wk = 0
     allocate(ip02_wk(lmnp2))    ; ip02_wk = 0
     allocate(ip21_wk(lmnp1))    ; ip21_wk = 0
     allocate(ip12_wk(lmnp2))    ; ip12_wk = 0
     allocate(iu21_wk(lmnp1))    ; iu21_wk = 0
     allocate(iv21_wk(lmnp1))    ; iv21_wk = 0
     allocate(nstar2_wk(lmnp2))  ; nstar2_wk = 0
     allocate(pa0_wk(3,lmnp0))  ; pa0_wk = 0
     allocate(pb0_wk(3,lmnp0))  ; pb0_wk = 0
     allocate(pb_wk(3,lmnp2))  ; pb_wk = 0
     allocate(ka0_wk(4,lmnp0))  ; ka0_wk = 0
!!$     allocate(ka2_wk(4,lmnp2))  ; ka2_wk = 0

! ===================================== modified by K. Tagami ========= 12.0A
!     call setkp0_default_n(il,ngen,inv,igen,jgen,a,b,c,ca,cb,cc &
!          & ,nx,ny,nz &
!          & ,np2,np1,np0,lmnp0,lmnp1,lmnp2 &
!          & ,nxx,nyy,nzz &
!          & ,ip10_wk,ip20_wk,ip01_wk,ip02_wk,ip21_wk,ip12_wk &
!          & ,iu21_wk,iv21_wk &
!          & ,nstar2_wk,pa0_wk,pb0_wk,pb_wk,ka0_wk &
!          & ,ipri_kp)
!
     if ( gen_tetramesh_mode == 0 ) then
        call setkp0_default_n(il,ngen,inv,igen,jgen,a,b,c,ca,cb,cc &
             & ,nx,ny,nz &
             & ,np2,np1,np0,lmnp0,lmnp1,lmnp2 &
             & ,nxx,nyy,nzz &
             & ,ip10_wk,ip20_wk,ip01_wk,ip02_wk,ip21_wk,ip12_wk &
             & ,iu21_wk,iv21_wk &
             & ,nstar2_wk,pa0_wk,pb0_wk,pb_wk,ka0_wk &
             & ,ipri_kp, itrs )
     else if ( gen_tetramesh_mode == 1 ) then
        call setkp0_default_n_kt(il,ngen,inv,igen,jgen,a,b,c,ca,cb,cc &
             & ,nx,ny,nz &
             & ,np2,np1,np0,lmnp0,lmnp1,lmnp2 &
             & ,nxx,nyy,nzz &
             & ,ip10_wk,ip20_wk,ip01_wk,ip02_wk,ip21_wk,ip12_wk &
             & ,iu21_wk,iv21_wk &
             & ,nstar2_wk,pa0_wk,pb0_wk,pb_wk,ka0_wk &
             & ,ipri_kp, &
             &  use_altv_rltv, altv, rltv, itrs, &
             &  gen_name_in_carts )
     else if ( gen_tetramesh_mode == 2 ) then
        call setkp0_default_n_kt2(il,ngen,inv,igen,jgen,a,b,c,ca,cb,cc &
             & ,nx,ny,nz &
             & ,np2,np1,np0,lmnp0,lmnp1,lmnp2 &
             & ,nxx,nyy,nzz &
             & ,ip10_wk,ip20_wk,ip01_wk,ip02_wk,ip21_wk,ip12_wk &
             & ,iu21_wk,iv21_wk &
             & ,nstar2_wk,pa0_wk,pb0_wk,pb_wk,ka0_wk &
             & ,ipri_kp, &
             &  use_altv_rltv, altv, rltv, itrs, &
             &  gen_name_in_carts )
     endif
! =========================================================================== 12.0A
  end if

  kv3 = np2
  nwei = np1

  kv3_fbz = np1

  if(paramset .and. ipri_kp >= 3) then
     allocate(trmat(3,3))
     allocate(trbp(3,3))
     allocate(trpb(3,3))
     allocate(mat1(3,3))
     allocate(mat2(3,3))

     call get_trmat1  !-(contained here) ->(trmat)

     write(nfout,'(" !! ALL vkxyz")')
     write(nfout,'(" !! nxx,nyy,nzz = ",3i8)') nxx,nyy,nzz
     do nn = 1, np0
        nv(1:3) = ka0_wk(1:3,nn)
        nk = ka0_wk(4,nn)
        ix = mod(nn-1,nxx+1)
        iy = mod((nn-ix-1)/(nxx+1),nyy+1)
        iz = (nn-ix-1-iy*(nxx+1))/((nxx+1)*(nyy+1))
        write(nfout,'(i6,6x,3f16.8,i6)') nn, matmul(trmat,nv)/dble(nk),ip20_wk(nn)
!!$        write(nfout,'(" vkxyz(1:3,",i5,") = (",3f10.6,") (",3i8,") (",3i8,") nk=",i8)') &
!!$             & nn,matmul(trmat,nv)/dble(nk), nv(1:3),ix*ny1*nz1,nx1*iy*nz1,nx1*ny1*iz,nk
     end do

     npx=nxx+1
     npy=nyy+1                                                                 
     icub = 0
     write(nfout,'(" tetrahedrons")')
     do iz = 0, nzz-1
        do iy = 0, nyy-1
           do ix = 0, nxx-1
              icub=icub+1
              ni=npx*(npy*iz+iy)+ix
              do kz=1,2
                 do ky=1,2
                    do kx=1,2
                       ip0 = ni+npx*(npy*(kz-1)+ky-1)+kx
                       ip8(kx+2*(ky-1)+4*(kz-1)) = ip0
                    end do
                 end do
              end do
!!$           write(nfout,'(" icub = ",i8)') icub
              write(nfout,'(i6," = ",i5,",",i5,",",i5,",",i5, "  (",i6,")" 8i5)') &
                   &                                            (icub-1)*6+1, ip8(1),ip8(2),ip8(4),ip8(8), icub, ip8(1:8)
              write(nfout,'(i6," = ",i5,",",i5,",",i5,",",i5)') (icub-1)*6+2, ip8(1),ip8(2),ip8(6),ip8(8)
              write(nfout,'(i6," = ",i5,",",i5,",",i5,",",i5)') (icub-1)*6+3, ip8(1),ip8(5),ip8(6),ip8(8)
              write(nfout,'(i6," = ",i5,",",i5,",",i5,",",i5)') (icub-1)*6+4, ip8(1),ip8(3),ip8(4),ip8(8)
              write(nfout,'(i6," = ",i5,",",i5,",",i5,",",i5)') (icub-1)*6+5, ip8(1),ip8(3),ip8(7),ip8(8)
              write(nfout,'(i6," = ",i5,",",i5,",",i5,",",i5)') (icub-1)*6+6, ip8(1),ip8(5),ip8(7),ip8(8)
!!$           do kz=1,2
!!$              do ky=1,2
!!$                 do kx=1,2
!!$                    ip0 = ni+npx*(npy*(kz-1)+ky-1)+kx
!!$                    iix = ix+kx-1
!!$                    iiy = iy+ky-1
!!$                    iiz = iz+kz-1
!!$                    write(nfout,'(" ip0 = ",i8,"   ->ip20 = ",i6,"   ka0_wk = ",4i6,"  (",3i6,")")') &
!!$                         & ip0,ip20_wk(ip0),ka0_wk(1:4,ip0),iix*ny1*nz1,nx1*iiy*nz1,nx1*ny1*iiz
!!$                    ip8(kx+2*(ky-1)+4*(kz-1)) = ip0
!!$                 end do
!!$              end do
!!$           end do
           end do
        end do
     end do
     write(nfout,'(" end of tetrahedrons")')
     deallocate(trmat)
     deallocate(trbp )
     deallocate(trpb )
     deallocate(mat1 )
     deallocate(mat2 )

  end if

  if(paramset) then
     deallocate(ip10_wk)
     deallocate(ip20_wk)
     deallocate(ip01_wk)
     deallocate(ip02_wk)
     deallocate(ip21_wk)
     deallocate(ip12_wk)
     deallocate(iu21_wk)
     deallocate(iv21_wk)
     deallocate(nstar2_wk)
     deallocate(pa0_wk)
     deallocate(pb0_wk)
     deallocate(pb_wk)
     deallocate(ka0_wk)
     return
  endif

  allocate(trmat(3,3))
  allocate(trbp(3,3))
  allocate(trpb(3,3))
  allocate(mat1(3,3))
  allocate(mat2(3,3))

  call get_trmat1  !-(contained here) ->(trmat)
  trmat_out = trmat

  do nn = 1, kv3
     nv(1) = ka0_wk(1,ip02_wk(nn))
     nv(2) = ka0_wk(2,ip02_wk(nn))
     nv(3) = ka0_wk(3,ip02_wk(nn))
     nk = ka0_wk(4,ip02_wk(nn))
     nw = nstar2_wk(nn)
     vkxyz(nn,1:3,CARTS) = matmul(trmat,nv)/dble(nk)
     qwgt(nn) = dble(nw)/dble(nwei)
  end do

  Do nn=1, np1
     nv(1) = ka0_wk(1,ip01_wk(nn))
     nv(2) = ka0_wk(2,ip01_wk(nn))
     nv(3) = ka0_wk(3,ip01_wk(nn))
     nk = ka0_wk(4,ip01_wk(nn))
     vkxyz_fbz(nn,1:3,CARTS) = matmul(trmat,nv)/dble(nk)
     to_ibz_from_fbz_for_kpoint(nn) = ip21_wk(nn)
  End do

  deallocate(trmat)
  deallocate(trbp )
  deallocate(trpb )
  deallocate(mat1 )
  deallocate(mat2 )


  deallocate(ip10_wk)
  deallocate(ip20_wk)
  deallocate(ip01_wk)
  deallocate(ip02_wk)
  deallocate(ip21_wk)
  deallocate(ip12_wk)
  deallocate(iu21_wk)
  deallocate(iv21_wk)
  deallocate(nstar2_wk)
  deallocate(pa0_wk)
  deallocate(pb0_wk)
  deallocate(pb_wk)
  deallocate(ka0_wk)
  if(nbztyp_spg == GENERAL .or.nbztyp_spg == GENERAL_LARGER) deallocate(ka2_wk)

  if(ipri >= 2) then
     write(nfout,'(/,"  << gnrt_k0_n >>")')
     write(nfout,*) ' !Total Generated Kpoints = ',kv3
     do nk = 1, kv3
        write(nfout,'(i4," ",3f12.6)') nk,(vkxyz(nk,nn,CARTS),nn=1,3)
     end do
  end if
contains
  subroutine get_trmat1
!    make translation matrix  trpb (P -> B)

    call getspgtab(trbp)  ! spg+tetra

    call inver3n(3,trbp,trpb)

    mat1 = transpose(trpb)
    call inver3n(3,mat1,mat2)
    call matpr3(rltv,mat2,trmat)

  end subroutine get_trmat1
end subroutine gnrt_k0_n
