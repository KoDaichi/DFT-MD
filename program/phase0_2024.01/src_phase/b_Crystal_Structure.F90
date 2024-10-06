!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE:  readop, altv_2_rltv, string2values, inver3n,
!              matpr3, inver2n, gnrt_op, gnrt_op_n,
!              bravais2primitive, primitive2bravais
!
!  AUTHOR(S): T. Yamasaki, H. Mizouchi, H.Sawada, K. Betsuyaku
!
!                                                August/20/2003
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
subroutine readop(ipri,kopr,altv,nfout,nfopgr,nfmatbp,Iparas,nopr,op,tau,af)
! $Id: b_Crystal_Structure.f90 570 2017-04-21 20:34:50Z yamasaki $
!
!      read operation matrix OP(3,3,kopr) and TAU(3,kopr)
!      from the file 'opgr.data'
!     you must prepare the file 'matrix.BP'
!       ( translation matrix from Bravais lattice coodinates to Primitive one)
!
!     #1) antiferromagnetic calculation is added on 9th Jul. 1996
!                                          by H.Sawada
!     
!     #x1) <rplcps> is moved to <b_Ionic_System.f90>, Jun. 2003
!                                      by T. Yamasaki(FUJITSU LABORATORIES Ltd.)
!     #x2) <gnrt_on_n> is newly added, Jun. 2003
!                                      by T. Yamasaki(FUJITSU LABORATORIES Ltd.)
!**********************************************************************
  use m_Const_Parameters, only : DP
  implicit none
  integer, intent(in)                       :: ipri,kopr
  real(DP), intent(in),dimension(3,3)       :: altv
  integer, intent(in)                       :: nfout,nfopgr,nfmatbp,iparas
  integer, intent(out)                      :: nopr,af
  real(DP), intent(out),dimension(3,3,kopr) :: op
  real(DP), intent(out),dimension(3,kopr)   :: tau

  integer, parameter ::  PARAMETER_SET = 0 , PARAMETER_NO_SET = 1
!     work
  integer            :: nop(3,3,48),ii,jj,iop
  real(DP)           :: taub(3,48),br(3,3),brinv(3,3),trbp(3,3),trpb(3,3)
  real(DP)           :: cc(3,3),rop(3,3)

  rewind nfopgr
  rewind nfmatbp

!     make translation matrix  trpb (P -> B)
  do ii = 1,3
     read(nfmatbp,*)  trbp(ii,1),trbp(ii,2),trbp(ii,3)
  end do

  call inver3n(3,trbp,trpb)
  write(nfout,*) 'TR B->P'
  do ii = 1,3
     write(nfout,100) trbp(ii,1),trbp(ii,2),trbp(ii,3)
  end do

  write(nfout,*) 'TR P->B'
  do ii = 1,3
     write(nfout,100) trpb(ii,1),trpb(ii,2),trpb(ii,3)
  end do

  call matpr3(altv,trpb,br)
  write(nfout,*) 'Bravais lattice'
  write(nfout,100) (br(1,ii),ii=1,3)
  write(nfout,100) (br(2,ii),ii=1,3)
  write(nfout,100) (br(3,ii),ii=1,3)
100 format(8x,3f15.6)

  call inver3n(3,br,brinv)

  read(nfopgr,*) nopr
  if(nopr < 0) then
     nopr=-nopr
     af=1
  else
     af=0
  end if
  if(Iparas == PARAMETER_NO_SET) then
     do iop = 1,nopr
        read(nfopgr,*) taub(1,iop),taub(2,iop),taub(3,iop)
     end do

     do iop = 1,nopr
        read(nfopgr,*) ((nop(ii,jj,iop),ii=1,3),jj=1,3)
     end do
  else
     do iop = 1, nopr
        read(nfopgr,*)
     end do
     do iop = 1, nopr
        read(nfopgr,*)
     enddo
  end if

  if(Iparas == PARAMETER_SET) then
     if(af /= 0) then
        read(nfopgr,*) af
     end if
     goto 1000
  else
     if(af /= 0) then
        read(nfopgr,*) af
        if(ipri >= 1) then 
           write(nfout,*) ' ANTI-FERROMAGNETIC CALCULATION'
           write(nfout,'(41h NUMBER OF ANTI-FERROMAGNETIC OPERATOR = ,i3)') af
        end if
        do iop=nopr+1,nopr+af
           read(nfopgr,*) taub(1,iop),taub(2,iop),taub(3,iop)
        end do
        do iop = nopr+1,nopr+af
           read(nfopgr,*) ((nop(ii,jj,iop),ii=1,3),jj=1,3)
        end do
     end if
  end if

  if(nopr+af < kopr) then 
     write(nfout,*) 'Warning ! kopr should be ',nopr+af
  else if(nopr+af.gt.kopr) then
     write(nfout,*) 'ERROR ! kopr must be ',nopr+af
  end if

  if(IPRI >= 1) then 
     write(nfout,*) 'OPERATION MATRIX ( Bravais lattice )'
     do iop=1,nopr
        write(nfout,*) iop
        write(nfout,300) (nop(1,jj,iop),jj=1,3),taub(1,iop)
        write(nfout,300) (nop(2,jj,iop),jj=1,3),taub(2,iop)
        write(nfout,300) (nop(3,jj,iop),jj=1,3),taub(3,iop)
     end do
     if(af /= 0) then
        write(nfout,*) 'ANTI-FERROMAGNETIC OPERATION'
        do iop = nopr+1,nopr+af
           write(nfout,300) (nop(1,jj,iop),jj=1,3),taub(1,iop)
           write(nfout,300) (nop(2,jj,iop),jj=1,3),taub(2,iop)
           write(nfout,300) (nop(3,jj,iop),jj=1,3),taub(3,iop)
        end do
     end if
  end if

  do iop=1,nopr+af
     do ii=1,3
        do jj=1,3
           rop(ii,jj)=dfloat(nop(ii,jj,iop))
        end do
     end do

     call matpr3(rop,brinv,cc)
     call matpr3(br,cc,op(1,1,iop))

     do ii=1,3
        tau(ii,iop)=br(ii,1)*taub(1,iop)+br(ii,2)*taub(2,iop) &
             & +br(ii,3)*taub(3,iop)
     end do
  end do

!     debug
  if(IPRI >= 1) then 
     write(nfout,*) 'OPERATION MATRIX ( Cartesian coodinates )'
     do iop = 1,nopr
        write(nfout,*) iop
        write(nfout,200) (op(1,jj,iop),jj=1,3),tau(1,iop)
        write(nfout,200) (op(2,jj,iop),jj=1,3),tau(2,iop)
        write(nfout,200) (op(3,jj,iop),jj=1,3),tau(3,iop)
     end do
     if(af /= 0) then
        write(nfout,*) 'ANTI-FERROMAGNETIC OPERATION'
        do iop=nopr+1,nopr+af
           write(nfout,200) (op(1,jj,iop),jj=1,3),tau(1,iop)
           write(nfout,200) (op(2,jj,iop),jj=1,3),tau(2,iop)
           write(nfout,200) (op(3,jj,iop),jj=1,3),tau(3,iop)
        end do
     end if
200  format(8x,3f12.5,5x,f12.5)
300  format(8x,3i4,f10.4)
  end if
      
1000 if(Iparas .eq. PARAMETER_SET) nopr=nopr+af

  return
end subroutine readop

subroutine altv_2_rltv(altv,rltv,univol,rvol)
  use m_Const_Parameters, only : PAI2, DP
  implicit none
  real(kind=DP), dimension(3,3), intent(in)  :: altv
  real(kind=DP), dimension(3,3), intent(out) :: rltv
  real(kind=DP), intent(out):: univol, rvol
  real(kind=DP), external :: deter3

  call inver3n(3, altv, rltv)
  rltv = PAI2*transpose(rltv)

  univol=dabs(deter3(altv))
  rvol = dabs(deter3(rltv))

end subroutine altv_2_rltv

function deter3(a)
  use m_Const_Parameters, only : DP
  implicit none
  real(kind=DP), intent(in), dimension(3,3) :: a
  real(kind=DP)                             :: deter3

  deter3 = a(1,1)*(a(2,2)*a(3,3)-a(3,2)*a(2,3)) &
       &     +   a(2,1)*(a(3,2)*a(1,3)-a(1,2)*a(3,3)) &
       &     +   a(3,1)*(a(1,2)*a(2,3)-a(2,2)*a(1,3))

end function deter3

function deter3n(k,a)
  use m_Const_Parameters, only : DP
  implicit none
  integer, intent(in)                       :: k
  real(kind=DP), intent(in), dimension(k,k) :: a
  real(kind=DP)                             :: deter3n

  deter3n = a(1,1)*(a(2,2)*a(3,3)-a(3,2)*a(2,3)) &
       &     +   a(2,1)*(a(3,2)*a(1,3)-a(1,2)*a(3,3)) &
       &     +   a(3,1)*(a(1,2)*a(2,3)-a(2,2)*a(1,3))

end function deter3n

subroutine string2values(len_str,str,n,values,nc)
  use m_Const_Parameters, only : DP
  implicit none
  integer, intent(in)             :: len_str
  character(len=len_str), intent(in) :: str
  integer, intent(in)             :: n
  real(kind=DP), intent(out),dimension(n) :: values
  integer, intent(out)            :: nc

  call chnnm(str,len_str,n,values,nc)

end subroutine string2values

subroutine inver3n(k,a,b)
  use m_Const_Parameters, only : DP
  implicit none
  integer, intent(in)                       :: k
  real(kind=DP), intent(in), dimension(k,k) :: a
  real(kind=DP), intent(out),dimension(k,k) :: b

  real(kind=DP) :: c,d,deter3n,cmax
  integer       :: i,j,l

  d=deter3n(k,a)

  if(dabs(d) < 1.d-20) then
     write(6,9001) abs(d)
9001 format(' ***** ERROR: (inver3) abs(det(a))'&
          &       ,d9.2,'<1.d-20')
     !stop ' ***** ERROR: (inver3) abs(det(a))<1.d-20'
     call phase_error_with_msg(6, '***** ERROR: (inver3) abs(det(a))<1.d-20',__LINE__,__FILE__)
  endif
  b(1,1)=(a(2,2)*a(3,3)-a(3,2)*a(2,3))/d
  b(2,1)=(a(2,3)*a(3,1)-a(3,3)*a(2,1))/d
  b(3,1)=(a(2,1)*a(3,2)-a(3,1)*a(2,2))/d
  b(1,2)=(a(3,2)*a(1,3)-a(1,2)*a(3,3))/d
  b(2,2)=(a(3,3)*a(1,1)-a(1,3)*a(3,1))/d
  b(3,2)=(a(3,1)*a(1,2)-a(1,1)*a(3,2))/d
  b(1,3)=(a(1,2)*a(2,3)-a(2,2)*a(1,3))/d
  b(2,3)=(a(1,3)*a(2,1)-a(2,3)*a(1,1))/d
  b(3,3)=(a(1,1)*a(2,2)-a(2,1)*a(1,2))/d

  cmax=0.d0
  do i=1,3
     do j=1,3
        c=0.d0
        do l=1,3
           c=c+b(i,l)*a(l,j)
        enddo
        if(i == j) then
           c=c-1.d0
        endif
        if(cmax < abs(c)) cmax = abs(c)
     enddo
  enddo
  if(cmax > 1.d-10) then
     write(6,9002) cmax
9002 format(' ***** ERROR: (inver3n) cmax=',d10.3,'>1.d-10')
     !stop ' ***** ERROR: (inver3n) cmax>1.d-10'
     call phase_error_with_msg(6, '***** ERROR: (inver3n) cmax>1.d-10',__LINE__,__FILE__)
  endif
end subroutine inver3n

subroutine matpr3(a,b,c)
  use m_Const_Parameters, only : DP
  implicit none
  real(kind=DP), intent(in),  dimension(3,3) :: a,b
  real(kind=DP), intent(out), dimension(3,3) :: c
!   C=A*B  A,B,C is 3x3 matrix

  c = matmul(a,b)

end subroutine matpr3

subroutine inver2n(k,a,b)
  use m_Const_Parameters, only : DP
  implicit none
  integer, intent(in)                       :: k
  real(kind=DP), intent(in), dimension(k,k) :: a
  real(kind=DP), intent(out),dimension(k,k) :: b

  real(kind=DP) :: c,d,deter2n,cmax
  integer       :: i,j,l

  d=deter2n(k,a)

  if(dabs(d) < 1.d-20) then
     write(6,9001) abs(d)
9001 format(' ***** ERROR: (inver2n) abs(det(a))'&
          &       ,d9.2,'<1.d-20')
     !stop ' ***** ERROR: (inver2n) abs(det(a))<1.d-20'
     call phase_error_with_msg(6, '***** ERROR: (inver2n) abs(det(a))<1.d-20',__LINE__,__FILE__)
  endif
  b(1,1) =  a(2,2)/d
  b(1,2) = -a(1,2)/d
  b(2,1) = -a(2,1)/d
  b(2,2) =  a(1,1)/d

  cmax=0.d0
  do i=1,2
     do j=1,2
        c=0.d0
        do l=1,2
           c=c+b(i,l)*a(l,j)
        enddo
        if(i == j) then
           c=c-1.d0
        endif
        if(cmax < abs(c)) cmax = abs(c)
     enddo
  enddo
  if(cmax > 1.d-10) then
     write(6,9002) cmax
9002 format(' ***** ERROR: (inver2n) cmax=',d10.3,'>1.d-10')
     !stop ' ***** ERROR: (inver2n) cmax>1.d-10'
     call phase_error_with_msg(6, '***** ERROR: (inver2n) cmax>1.d-10',__LINE__,__FILE__)
  endif
end subroutine inver2n

function deter2n(k,a)
  use m_Const_Parameters, only : DP
  implicit none
  integer, intent(in)                       :: k
  real(kind=DP), intent(in), dimension(k,k) :: a
  real(kind=DP)                             :: deter2n

  deter2n = a(1,1)*a(2,2) - a(1,2)*a(2,1)

end function deter2n

subroutine gnrt_op(kopr,altv,nfout,Iparas,nopr,op,tau,af,nbztyp_spg,nfspg,ipri_spg)
!
!      read bnprp4.i5 in setspg
!
!     #1) antiferromagnetic calculation is added on 9th Jul. 1996
!                                          by H.Sawada
!**********************************************************************
  use m_Const_Parameters, only : DP
  implicit none
  integer, intent(in)                       :: kopr, nbztyp_spg, ipri_spg
  real(DP), intent(in),dimension(3,3)       :: altv
  integer, intent(in)                       :: nfout,iparas,nfspg ! 
  integer, intent(out)                      :: nopr,af
  real(DP), intent(out),dimension(3,3,kopr) :: op
  real(DP), intent(out),dimension(3,kopr)   :: tau

  integer, parameter ::  PARAMETER_SET = 0 , PARAMETER_NO_SET = 1
!     work
  integer            :: nop(3,3,48),ii,jj,iop
  real(DP)           :: taub(3,48),br(3,3),brinv(3,3),trbp(3,3),trpb(3,3)
  real(DP)           :: cc(3,3),rop(3,3)
  real(DP)           :: tab(3,3),ta1(3,48)
  integer            :: ng1,lra1(3,3,48)
  integer            :: jpri_spg

  jpri_spg = ipri_spg
!! setspg  !!


  if(nbztyp_spg == 100 .or.nbztyp_spg == 101) then
      call setspg(tab,ng1,ta1,lra1,nfspg,jpri_spg)
  else
      call setspg_default(nbztyp_spg,altv,tab,ng1,ta1,lra1,jpri_spg)
  end if

!     make translation matrix  trpb (P -> B)
  do ii = 1,3
     trbp(ii,1) = tab(ii,1)
     trbp(ii,2) = tab(ii,2)
     trbp(ii,3) = tab(ii,3)
  end do

  call inver3n(3,trbp,trpb)

  if(jpri_spg >= 0) then
      write(nfout,*) 'TR B->P'
      do ii = 1,3
         write(nfout,100) trbp(ii,1),trbp(ii,2),trbp(ii,3)
      end do
      write(nfout,*) 'TR P->B'
      do ii = 1,3
         write(nfout,100) trpb(ii,1),trpb(ii,2),trpb(ii,3)
      end do
  end if
  call matpr3(altv,trpb,br)

  if(jpri_spg >= 0) then
      write(nfout,*) 'Bravais lattice'
      write(nfout,100) (br(1,ii),ii=1,3)
      write(nfout,100) (br(2,ii),ii=1,3)
      write(nfout,100) (br(3,ii),ii=1,3)
  end if

100 format(8x,3f15.6)

  call inver3n(3,br,brinv)

  nopr = ng1 
  if(nopr < 0) then
     nopr=-nopr
     af=1
  else
     af=0
  end if
  if(Iparas == PARAMETER_NO_SET) then
     do iop = 1,nopr
         taub(1,iop) = ta1(1,iop)
         taub(2,iop) = ta1(2,iop)
         taub(3,iop) = ta1(3,iop)
     end do

     do iop = 1,nopr
        do jj =1,3
            do ii = 1,3
                nop(ii,jj,iop) = lra1(ii,jj,iop) 
            end do
        end do
     end do
  end if

  if(Iparas == PARAMETER_SET) then
     if(af /= 0) then
        af = ng1
     end if
     goto 1000
  else
     if(af /= 0) then
        af = ng1
        if(jpri_spg >= -1) then 
           write(nfout,*) ' ANTI-FERROMAGNETIC CALCULATION'
           write(nfout,'(41h NUMBER OF ANTI-FERROMAGNETIC OPERATOR = ,i3)') af
        end if
        do iop=nopr+1,nopr+af
           taub(1,iop) = ta1(1,iop)
           taub(2,iop) = ta1(2,iop)
           taub(3,iop) = ta1(3,iop)
        end do
        do iop = nopr+1,nopr+af
           do jj =1,3
               do ii = 1,3
                   nop(ii,jj,iop) = lra1(ii,jj,iop) 
               end do
           end do
        end do
     end if
  end if

  if(nopr+af < kopr) then 
     write(nfout,*) 'Warning ! kopr should be ',nopr+af
  else if(nopr+af.gt.kopr) then
     write(nfout,*) 'ERROR ! kopr must be ',nopr+af
  end if

  if(jpri_spg >= -1) then 
     write(nfout,*) 'OPERATION MATRIX ( Bravais lattice )'
     do iop=1,nopr
        write(nfout,*) iop
        write(nfout,300) (nop(1,jj,iop),jj=1,3),taub(1,iop)
        write(nfout,300) (nop(2,jj,iop),jj=1,3),taub(2,iop)
        write(nfout,300) (nop(3,jj,iop),jj=1,3),taub(3,iop)
     end do
     if(af /= 0) then
        write(nfout,*) 'ANTI-FERROMAGNETIC OPERATION'
        do iop = nopr+1,nopr+af
           write(nfout,300) (nop(1,jj,iop),jj=1,3),taub(1,iop)
           write(nfout,300) (nop(2,jj,iop),jj=1,3),taub(2,iop)
           write(nfout,300) (nop(3,jj,iop),jj=1,3),taub(3,iop)
        end do
     end if
  end if

  do iop=1,nopr+af
     do ii=1,3
        do jj=1,3
           rop(ii,jj)=dfloat(nop(ii,jj,iop))
        end do
     end do

     call matpr3(rop,brinv,cc)
     call matpr3(br,cc,op(1,1,iop))

     do ii=1,3
        tau(ii,iop)=br(ii,1)*taub(1,iop)+br(ii,2)*taub(2,iop) &
             & +br(ii,3)*taub(3,iop)
     end do
  end do

  if(jpri_spg >= -1) then 
     write(nfout,*) 'OPERATION MATRIX ( Cartesian coodinates )'
     do iop = 1,nopr
        write(nfout,*) iop
        write(nfout,200) (op(1,jj,iop),jj=1,3),tau(1,iop)
        write(nfout,200) (op(2,jj,iop),jj=1,3),tau(2,iop)
        write(nfout,200) (op(3,jj,iop),jj=1,3),tau(3,iop)
     end do
     if(af /= 0) then
        write(nfout,*) 'ANTI-FERROMAGNETIC OPERATION'
        do iop=nopr+1,nopr+af
           write(nfout,200) (op(1,jj,iop),jj=1,3),tau(1,iop)
           write(nfout,200) (op(2,jj,iop),jj=1,3),tau(2,iop)
           write(nfout,200) (op(3,jj,iop),jj=1,3),tau(3,iop)
        end do
     end if
200  format(8x,3f12.5,5x,f12.5)
300  format(8x,3i4,f10.4)
  end if
      
1000 if(Iparas .eq. PARAMETER_SET) nopr=nopr+af

  return
end subroutine gnrt_op

! ================================ modified by K. Tagami ============== 12.0A
!subroutine gnrt_op_paramset(il,ngen,inv,igen,jgen,imag,iaf,jaf,a,b,c,ca,cb,cc &
!  &     , nfout,nopr,af,nbztyp_spg,ig01)
!
subroutine gnrt_op_paramset( il, ngen, inv, igen, jgen, imag, iaf, jaf, &
     &                       a, b, c, ca, cb, cc, &
     &                       nfout, nopr, af, nbztyp_spg, ig01, &
     &                       use_altv_rltv, altv, rltv, gen_name_in_carts  )
! ========================================================================= 12.0A

!
!     #1) antiferromagnetic calculation is added on 9th Jul. 1996
!                                          by H.Sawada
!
!!!  subroutine gnrt_op(ipri,kopr,altv,nfout,Iparas,nopr,op,tau,af,nbztyp_spg)
!!!       modified by mizouchi@adv 2003.2.21 !!!!!
!
!  Subroutine <gnrt_op_n> is revised on <gnrt_op> by T. Yamasaki (FUJITSU Lab.)
!                                                 31 May 2003
!**********************************************************************
  use m_Const_Parameters, only : DP
  implicit none
  integer, intent(in) ::                       il,ngen,inv
  integer, intent(in), dimension(ngen) ::      igen
  integer, intent(in), dimension(2,3,ngen) ::  jgen
  integer, intent(in) ::                       imag !if(imag==1) Antiferro
  integer, intent(in) ::                       iaf
  integer, intent(in), dimension(2,3) ::       jaf
  real(kind=DP),intent(in) ::                  a,b,c,ca,cb,cc
  integer, intent(in) ::                       nbztyp_spg
  integer, intent(in) ::                       nfout
  integer, intent(out) ::                      nopr,af
  integer, intent(out),dimension(48) ::   ig01

! ================================== added by K. Tagami ================ 12.0A
  logical, intent(in) :: use_altv_rltv
  logical, intent(in) :: gen_name_in_carts
  real(DP), intent(in) :: rltv(3,3), altv(3,3)
! ======================================================================= 12.0A

!!$  real(DP)           :: br(3,3),brinv(3,3),trbp(3,3),trpb(3,3)
!!$  real(DP)           :: tab(3,3),ta1(3,48)
!!$  integer            :: ng1,lra1(3,3,48)
  real(DP)           :: tab(3,3),ta1(3,49)
  integer            :: ng1,lra1(3,3,49)

  if(nbztyp_spg == 100 .or.nbztyp_spg == 101) then
      call setspg_n(il,ngen,inv,igen,jgen,a,b,c,ca,cb,cc,tab,ng1,ta1,lra1,ig01)
!           setspg_n is modified from <setspg> by T. Yamasaki
  else

! ================================= modified by K. Tagami ================= 12.0A
!      call setspg_default_n(il,ngen,inv,igen,jgen,imag,iaf,jaf,a,b,c,ca,cb,cc,tab,ng1,ta1,lra1,0,ig01)
!!          <setspg_default_n> is modified from <setspg_default> by T. Yamasaki
!
     if ( use_altv_rltv ) then
        ng1 = 0
        call setspg_default_n_kt( il, ngen, inv, igen, jgen, imag, iaf, jaf, &
             &                    a, b, c, ca, cb, cc, tab, ng1, ta1, lra1, &
             &                    0, ig01, use_altv_rltv, altv, rltv, &
             &                    gen_name_in_carts )
     else
        call setspg_default_n( il, ngen, inv, igen, jgen, imag, iaf, jaf, &
             &                 a, b, c, ca, cb, cc, tab, ng1, ta1, lra1, 0, ig01 )
     endif
! ============================================================================ 12.0A

  end if

!     make translation matrix  trpb (P -> B)
!!$  do ii = 1,3
!!$     trbp(ii,1) = tab(ii,1)
!!$     trbp(ii,2) = tab(ii,2)
!!$     trbp(ii,3) = tab(ii,3)
!!$  end do
!!$
!!$  call inver3n(3,trbp,trpb)
!!$
!!$  call matpr3(altv,trpb,br)
!!$
!!$  call inver3n(3,br,brinv)

  nopr = ng1 
  if(imag == 1) then
  !!if(nopr < 0) then
  !!   nopr=-nopr
     af=1
  else
     af=0
  end if

  return
end subroutine gnrt_op_paramset

! ====================================== modified by K. Tagami =========== 12.0A
!subroutine gnrt_op_n(il,ngen,inv,igen,jgen,imag,iaf,jaf,a,b,c,ca,cb,cc &
!  &     , kopr,altv,nfout,nopr,op,tau,af,nbztyp_spg,ipri_spg,ig01)
!
subroutine gnrt_op_n( il, ngen, inv, igen, jgen, imag, iaf, jaf, &
     &                a, b, c, ca, cb, cc, &
     &                kopr, altv, nfout, nopr, op, tau, af, &
     &                nbztyp_spg, ipri_spg, ig01, &
     &                use_altv_rltv, rltv, gen_name_in_carts )
! ========================================================================== 12.0A

!
!     #1) antiferromagnetic calculation is added on 9th Jul. 1996
!                                          by H.Sawada
!
!!!  subroutine gnrt_op(ipri,kopr,altv,nfout,Iparas,nopr,op,tau,af,nbztyp_spg)
!!!       modified by mizouchi@adv 2003.2.21 !!!!!
!
!  Subroutine <gnrt_op_n> is revised on <gnrt_op> by T. Yamasaki (FUJITSU Lab.)
!                                                 31 May 2003
!**********************************************************************
  use m_Const_Parameters, only : DP, oh_symbol, d6h_symbol
  implicit none
  integer, intent(in) ::                       il,ngen,inv
  integer, intent(in), dimension(ngen) ::      igen
  integer, intent(in), dimension(2,3,ngen) ::  jgen
  integer, intent(in) ::                       imag !if(imag==1) Antiferro
  integer, intent(in) ::                       iaf
  integer, intent(in), dimension(2,3) ::       jaf
  real(kind=DP),intent(in) ::                  a,b,c,ca,cb,cc
  integer, intent(in) ::                       kopr, nbztyp_spg, ipri_spg
  real(DP), intent(in),dimension(3,3) ::       altv
  integer, intent(in) ::                       nfout
  integer, intent(out) ::                      nopr,af
  real(DP), intent(out),dimension(3,3,kopr) :: op
  real(DP), intent(out),dimension(3,kopr) ::   tau
  integer, intent(out),dimension(48) ::   ig01

! ================================== added by K. Tagami ================ 12.0A
  logical, intent(in) :: use_altv_rltv
  logical, intent(in) :: gen_name_in_carts
  real(DP), intent(in) :: rltv(3,3)
! ======================================================================= 12.0A

  integer, parameter ::  PARAMETER_SET = 0 , PARAMETER_NO_SET = 1
!     work
!!$  integer            :: nop(3,3,48),ii,jj,iop
  integer            :: nop(3,3,49),ii,jj,iop
!!$  real(DP)           :: taub(3,48),br(3,3),brinv(3,3),trbp(3,3),trpb(3,3)
  real(DP)           :: taub(3,49),br(3,3),brinv(3,3),trbp(3,3),trpb(3,3)
  real(DP)           :: ccc(3,3),rop(3,3)
!!$  real(DP)           :: tab(3,3),ta1(3,48)
  real(DP)           :: tab(3,3),ta1(3,49)
!!$  integer            :: ng1,lra1(3,3,48)
  integer            :: ng1,lra1(3,3,49)


  if(nbztyp_spg == 100 .or.nbztyp_spg == 101) then
      call setspg_n(il,ngen,inv,igen,jgen,a,b,c,ca,cb,cc,tab,ng1,ta1,lra1,ig01)
!           setspg_n is modified from <setspg> by T. Yamasaki
  else

! ============================ modified by K. Tagami ===================== 12.0A
!      call setspg_default_n(il,ngen,inv,igen,jgen,imag,iaf,jaf,a,b,c,ca,cb,cc,tab,ng1,ta1,lra1,ipri_spg,ig01)
!!          <setspg_default_n> is modified from <setspg_default> by T. Yamasaki
!
     if ( use_altv_rltv ) then
        ng1 = 0
        call setspg_default_n_kt( il, ngen, inv, igen, jgen, imag, iaf, jaf, &
             &                    a, b, c, ca, cb, cc, tab, ng1, ta1, lra1, &
             &                    ipri_spg, ig01, &
             &                    use_altv_rltv, altv, rltv, gen_name_in_carts )
     else
        call setspg_default_n( il, ngen, inv, igen, jgen, imag, iaf, jaf, &
             &                 a, b, c, ca, cb, cc, tab, ng1, ta1, lra1, &
             &                 ipri_spg, ig01 )
     end if
! ========================================================================= 12.0A     

  endif

!---
!     make translation matrix  trpb (P -> B)
  do ii = 1,3
     trbp(ii,1) = tab(ii,1)
     trbp(ii,2) = tab(ii,2)
     trbp(ii,3) = tab(ii,3)
  end do

  call inver3n(3,trbp,trpb)

  if(ipri_spg >= 1) then
      write(nfout,*) 'TR B->P'
      do ii = 1,3
         write(nfout,100) trbp(ii,1),trbp(ii,2),trbp(ii,3)
      end do
      write(nfout,*) 'TR P->B'
      do ii = 1,3
         write(nfout,100) trpb(ii,1),trpb(ii,2),trpb(ii,3)
      end do
  end if

  call matpr3(altv,trpb,br)


  if(ipri_spg >= 1) then
      write(nfout,*) 'Bravais lattice'
      write(nfout,100) (br(1,ii),ii=1,3)
      write(nfout,100) (br(2,ii),ii=1,3)
      write(nfout,100) (br(3,ii),ii=1,3)
  end if
100 format(8x,3f15.6)

  call inver3n(3,br,brinv)

  nopr = ng1 
  if(imag == 1) then
  !!if(nopr < 0) then
  !!   nopr=-nopr
     af=1
  else
     af=0
  end if
!!$  if(Iparas == PARAMETER_NO_SET) then
  do iop = 1,nopr
     taub(1,iop) = ta1(1,iop)
     taub(2,iop) = ta1(2,iop)
     taub(3,iop) = ta1(3,iop)
  end do

  do iop = 1,nopr
     do jj =1,3
        do ii = 1,3
           nop(ii,jj,iop) = lra1(ii,jj,iop) 
        end do
     end do
  end do

!!$  if(Iparas == PARAMETER_SET) then
!!$     !!if(af /= 0) then
!!$     !!   af = ng1
!!$     !!end if
!!$     goto 1000
!!$  else
  if(af /= 0) then
     !!af = ng1
     if(ipri_spg >= 1) then 
        write(nfout,*) ' ANTI-FERROMAGNETIC CALCULATION'
        write(nfout,'(41h NUMBER OF ANTI-FERROMAGNETIC OPERATOR = ,i3)') af
     end if
     do iop=nopr+1,nopr+af
        taub(1,iop) = ta1(1,iop)
        taub(2,iop) = ta1(2,iop)
        taub(3,iop) = ta1(3,iop)
     end do
     do iop = nopr+1,nopr+af
        do jj =1,3
           do ii = 1,3
              nop(ii,jj,iop) = lra1(ii,jj,iop) 
           end do
        end do
     end do
  end if
!!$  end if

  if(nopr+af < kopr) then 
     write(nfout,*) 'Warning ! kopr should be ',nopr+af
  else if(nopr+af.gt.kopr) then
     write(nfout,*) 'ERROR ! kopr must be ',nopr+af
  end if

  if(ipri_spg >= 1) then 
     write(nfout,*) 'OPERATION MATRIX ( Bravais lattice )'
     do iop=1,nopr
        if(il<=0) then
           write(nfout,*) iop,d6h_symbol(ig01(iop))
        else
           write(nfout,*) iop,oh_symbol(ig01(iop))
        end if
        write(nfout,300) (nop(1,jj,iop),jj=1,3),taub(1,iop)
        write(nfout,300) (nop(2,jj,iop),jj=1,3),taub(2,iop)
        write(nfout,300) (nop(3,jj,iop),jj=1,3),taub(3,iop)
     end do
     if(af /= 0) then
        write(nfout,*) 'ANTI-FERROMAGNETIC OPERATION'
        do iop = nopr+1,nopr+af
           if(il<=0) then
              write(nfout,*) iop,d6h_symbol(iaf)
           else
              write(nfout,*) iop,oh_symbol(iaf)
           end if
           write(nfout,300) (nop(1,jj,iop),jj=1,3),taub(1,iop)
           write(nfout,300) (nop(2,jj,iop),jj=1,3),taub(2,iop)
           write(nfout,300) (nop(3,jj,iop),jj=1,3),taub(3,iop)
        end do
     end if
  end if

  do iop=1,nopr+af
     do ii=1,3
        do jj=1,3
           rop(ii,jj)=dfloat(nop(ii,jj,iop))
        end do
     end do

     call matpr3(rop,brinv,ccc)
     call matpr3(br,ccc,op(1,1,iop))

     do ii=1,3
        tau(ii,iop)=br(ii,1)*taub(1,iop)+br(ii,2)*taub(2,iop) &
             & +br(ii,3)*taub(3,iop)
     end do
  end do

  if(ipri_spg >= 1) then 
     write(nfout,*) 'OPERATION MATRIX ( Cartesian coodinates )'
     do iop = 1,nopr
        if(il>0) then
           write(nfout,*) iop,oh_symbol(ig01(iop))
        else
           write(nfout,*) iop,d6h_symbol(ig01(iop))
        end if
        write(nfout,200) (op(1,jj,iop),jj=1,3),tau(1,iop)
        write(nfout,200) (op(2,jj,iop),jj=1,3),tau(2,iop)
        write(nfout,200) (op(3,jj,iop),jj=1,3),tau(3,iop)
     end do
     if(af /= 0) then
        write(nfout,*) 'ANTI-FERROMAGNETIC OPERATION'
        do iop=nopr+1,nopr+af
           if(il>0) then
              write(nfout,*) iop,oh_symbol(iaf)
           else
              write(nfout,*) iop,d6h_symbol(iaf)
           end if
           write(nfout,200) (op(1,jj,iop),jj=1,3),tau(1,iop)
           write(nfout,200) (op(2,jj,iop),jj=1,3),tau(2,iop)
           write(nfout,200) (op(3,jj,iop),jj=1,3),tau(3,iop)
        end do
     end if
200  format(8x,3f12.5,5x,f12.5)
300  format(8x,3i4,f10.4)
  end if
      
  return
end subroutine gnrt_op_n


subroutine bravais2primitive(nfout,b2pmat,a,b,c,ca,cb,cc,avec,bvec,cvec,il)
! This subroutine was coded by BETSUYAKU, K. (Fuji Research Institute Co., Ltd.), July 2003
!
! [use m_Crystal_Structure] is commented out and b2pmat is included in the argument list,
! and write and stop statements have been revised.
!   by T. Yamasaki (FUJITSU LABORATORIES Ltd.), 10th Aug. 2003

  use m_Const_Parameters, only  : DP
!!$  use m_Crystal_Structure, only : b2pmat
  implicit none
  integer, intent(in) ::          nfout
  real(kind=DP), intent(in), dimension(3,3) :: b2pmat
  real(kind=DP), intent(in)    :: a, b, c, ca, cb, cc
  real(kind=DP), intent(inout) :: avec(3), bvec(3), cvec(3)
  integer,       intent(in)    :: il

  real(kind=DP)                :: abcmat(3,3) = 0.d0
  real(kind=DP), parameter     :: epsilon = 1.d-10
  integer                      :: ido, jdo

  abcmat(1,:) = avec(:)
  abcmat(2,:) = bvec(:)
  abcmat(3,:) = cvec(:)

  do ido = 1, 3
     do jdo = 1, 3
        if (abs(abcmat(ido,jdo)) < epsilon) abcmat(ido,jdo) = 0.d0
     end do
  end do

  if (il == -1) then                  ! Trigonal
     if (abs(ca)>epsilon .or. abs(cb)>epsilon .or. abs(abs(cc)-0.5d0)>epsilon) then
!!$        write(nfout,'(" You must have given wrong lattice.: ca, cb, cc = ",3f8.4," (il=-1)")') ca, cb, cc
        write(nfout,1001) ca,cb,cc,il
        call phase_error_with_msg(nfout, 'il = -1 <<bravais2primitive>>',__LINE__,__FILE__)
        !stop ' il = -1 <<bravais2primitive>>'
     end if
  else if (il == 0) then ! Primitive (Hexagonal)
     if (abs(ca)>epsilon .or. abs(cb)>epsilon .or. abs(abs(cc)-0.5d0)>epsilon) then
!!$        write(nfout,'(" You must have given wrong lattice.: ca, cb, cc = ",3f8.4,', ca, cb, cc
        write(nfout,1001) ca,cb,cc,il
        !stop ' il = 0  <<bravais2primitive>>'
        call phase_error_with_msg(nfout, 'il = 0 <<bravais2primitive>>',__LINE__,__FILE__)
     end if
  else if (il == 1) then ! Primitive

  else if (il == 2) then              ! Face-centered
     if (abs(ca)>epsilon .or. abs(cb)>epsilon .or. abs(cc)>epsilon) then
!!$        write(nfout,*) 'You must have given wrong lattice.', ca, cb, cc
        write(nfout,1001) ca,cb,cc,il
        !stop ' il = 2 <<bravais2primitive>>'
        call phase_error_with_msg(nfout, 'il = 2 <<bravais2primitive>>',__LINE__,__FILE__)
     end if
  else if (il == 3) then              ! Body-centered
     if (abs(ca)>epsilon .or. abs(cb)>epsilon .or. abs(cc)>epsilon) then
!!$        write(nfout,*) 'You must have given wrong lattice.', ca, cb, cc
        write(nfout,1001) ca,cb,cc,il
        !stop ' il = 3 <<bravais2primitive>>'
        call phase_error_with_msg(nfout, 'il = 3 <<bravais2primitive>>',__LINE__,__FILE__)
     end if
  else if (il == 4) then              ! Base-centered
     if ( (abs(cb)>epsilon .and. (abs(ca)>epsilon .or. abs(cc)>epsilon)) .or. &
          (abs(cc)>epsilon .and. (abs(cb)>epsilon .or. abs(ca)>epsilon)) .or. &
          (abs(ca)>epsilon .and. (abs(cc)>epsilon .or. abs(cb)>epsilon)) ) then
!!$        write(nfout,*) 'You must have given wrong lattice.', ca, cb, cc
        write(nfout,1001) ca,cb,cc,il
        !stop ' il = 4 <<bravais2primitive>>'
        call phase_error_with_msg(6, 'il = 4 <<bravais2primitive>>',__LINE__,__FILE__)
     end if
  else
     !write(nfout,*) 'The value of il is invalid. il =', il
     call phase_error_with_msg(nfout, 'The value of il is invalid.',__LINE__,__FILE__)
     !stop
  end if
1001 format(" You must have given wrong lattice. : ca, cb, cc = ",3f8.4," (il = ",i2,")")

  abcmat = matmul(b2pmat,abcmat)

  do ido = 1, 3
     do jdo = 1, 3
        if (abs(abcmat(ido,jdo)) < epsilon) abcmat(ido,jdo) = 0.d0
     end do
  end do

  avec(:) = abcmat(1,:)
  bvec(:) = abcmat(2,:)
  cvec(:) = abcmat(3,:)

end subroutine bravais2primitive


subroutine primitive2bravais(nfout,p2bmat,avec,bvec,cvec,a,b,c,ca,cb,cc,il)
! This subroutine was coded by BETSUYAKU, K. (Fuji Research Institute Co., Ltd.), July 2003
!
! [use m_Crystal_Structure] is commented out and b2pmat is included in the argument list
! and write and stop statements have been revised.
!   by T. Yamasaki (FUJITSU LABORATORIES Ltd.), 10th Aug. 2003
  use m_Const_Parameters, only  : DP
  implicit none
!!$  use m_Crystal_Structure, only : p2bmat
  integer, intent(in) ::          nfout
  real(kind=DP), intent(in), dimension(3,3) :: p2bmat
  real(kind=DP), intent(in)    :: avec(3), bvec(3), cvec(3)
  real(kind=DP), intent(inout) :: a, b, c, ca, cb, cc
  integer,       intent(in)    :: il

  real(kind=DP)                :: abcmat(3,3) = 0.d0
  real(kind=DP), parameter     :: epsilon = 1.d-10

  abcmat(1,:) = avec(:)
  abcmat(2,:) = bvec(:)
  abcmat(3,:) = cvec(:)

  abcmat = matmul(p2bmat,abcmat)

  a = sqrt( sum(abcmat(1,:)*abcmat(1,:)) )
  b = sqrt( sum(abcmat(2,:)*abcmat(2,:)) )
  c = sqrt( sum(abcmat(3,:)*abcmat(3,:)) )

  ca =  dot_product(abcmat(2,:),abcmat(3,:))/(b*c)
  cb =  dot_product(abcmat(3,:),abcmat(1,:))/(c*a)
  cc =  dot_product(abcmat(1,:),abcmat(2,:))/(a*b)

  if (il <= 0) then      ! Trigonal or Hexagonal
     if ( abs(b-a) < epsilon ) b = a
     if ( abs(ca) < epsilon ) ca = 0.d0
     if ( abs(cb) < epsilon ) cb = 0.d0
     if ( abs(cc+0.5d0) < epsilon ) cc = -0.5d0
  else if (il >= 1) then ! the others
     if ( abs(b-a) < epsilon ) b = a
     if ( abs(c-a) < epsilon ) c = a
     if ( abs(ca) < epsilon ) ca = 0.d0
     if ( abs(cb) < epsilon ) cb = 0.d0
     if ( abs(cc) < epsilon ) cc = 0.d0
  else
     write(nfout,'(" The value of il is invalid. il =",i3)') il
     !stop ' The value of il is invalid. <<primitive2bravais>>'
     call phase_error_with_msg(nfout,' The value of il is invalid. <<primitive2bravais>>',__LINE__,__FILE__)
  end if

end subroutine primitive2bravais

logical function is_hexagonal(a,b,ca,cb,cc)
  use m_Const_Parameters, only : DP
  implicit none
  real(kind=DP), intent(in) :: a,b
  real(kind=DP), intent(in) :: ca,cb,cc
  real(kind=DP), parameter     :: epsilon = 1.d-10
  is_hexagonal = .false.
  if (abs(ca)<epsilon .and. abs(cb)<epsilon .and. abs(abs(cc)-0.5d0)<epsilon .and. abs(a-b)<epsilon) &
  & is_hexagonal=.true.

end function is_hexagonal

subroutine get_latvec_from_brav_to_prim( b2pmat, altv_brav, altv_pr )
  use m_Const_Parameters, only  : DP
  implicit none

  real(kind=DP), intent(in), dimension(3,3) :: b2pmat
  real(kind=DP), intent(in)   :: altv_brav(3,3)
  real(kind=DP), intent(out)    :: altv_pr(3,3)

  real(kind=DP)                :: abcmat(3,3) = 0.d0

  abcmat(1,:) = altv_brav(:,1)
  abcmat(2,:) = altv_brav(:,2)
  abcmat(3,:) = altv_brav(:,3)

  abcmat = matmul(b2pmat,abcmat)
  altv_pr(:,1) = abcmat(1,:)
  altv_pr(:,2) = abcmat(2,:)
  altv_pr(:,3) = abcmat(3,:)

end subroutine get_latvec_from_brav_to_prim

subroutine get_latvec_from_prim_to_brav( p2bmat, altv_pr, altv_brav )
  use m_Const_Parameters, only  : DP
  implicit none

  real(kind=DP), intent(in), dimension(3,3) :: p2bmat
  real(kind=DP), intent(in)    :: altv_pr(3,3)
  real(kind=DP), intent(out)   :: altv_brav(3,3)

  real(kind=DP)                :: abcmat(3,3) = 0.d0

  abcmat(1,:) = altv_pr(:,1)
  abcmat(2,:) = altv_pr(:,2)
  abcmat(3,:) = altv_pr(:,3)

  abcmat = matmul(p2bmat,abcmat)
  altv_brav(:,1) = abcmat(1,:)
  altv_brav(:,2) = abcmat(2,:)
  altv_brav(:,3) = abcmat(3,:)

end subroutine get_latvec_from_prim_to_brav
