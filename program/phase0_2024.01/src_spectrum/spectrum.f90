!=======================================================================
!
!  PROGRAM  PHASE/0 2014.01 (rev.375)
!
!  "First-principles Electronic Structure Calculation Program"
!
!  PROGRAM: dipole_strength_tensor
!
!  AUTHOR(S):  H. Momida
!  
!  Contact address :  Phase System Consortium
!                     E-mail: phase_system@nims.go.jp URL https://azuma.nims.go.jp
!  
!
!
!=======================================================================
!
!   The original version of this set of the computer programs "PHASE" was developed by 
!  the members of the Theory Group of Joint Research Center for Atom Technology 
!  (JRCAT), based in Tsukuba, in the period 1993-2001.  
!   Since 2002, this program set had been intensively developed as a part of the following 
!  national projects supported by the Ministry of Education, Culture, Sports, Science and 
!  Technology (MEXT) of Japan; "Frontier Simulation Software for Industrial Science 
!  (FSIS)" from 2002 to 2005, "Revolutionary Simulation Software (RSS21)" from 2006 to 
!  2008. "Research and Development of Innovative Simulation Software (RISS)" from 2008 
!  to 2013. These projects is lead by the Center for Research on Innovative Simulation 
!  Software (CISS), the Institute of Industrial Science (IIS), the University of Tokyo.
!   Since 2013, this program set has been further developed centering on PHASE System 
!  Consortium. 
!   The activity of development of this program set has been supervised by Takahisa Ohno.
!
program dipole_strength_tensor
implicit none
integer :: nt,ndata,i,j,omega_max,id
real*8 :: epsilon,gamma,dt,time,dipl,tmp,omega_delta,time_delta,fac,pi
real*8,allocatable,dimension(:) :: t,omega
real*8,allocatable,dimension(:,:) :: c,d
complex*16,allocatable,dimension(:,:) :: d_w,alpha,S
complex*16,parameter :: zi=(0.0d0,1.0d0)
character(len=2) :: ctmp2

!==== parameters
pi=acos(-1.0d0)

!==== input parameters
!!time_delta=(0.002d-15)/(2.418884327d-17)    !! input parameter (time step)
!time_delta=0.1                             !! input parameter (time step)
!epsilon=0.01      !! input parameter (magnitude of impulse electric filed)

!==== input and output files
open(10,file='j.data',status='old')
open(11,file='j.out',status='unknown')
open(12,file='p.out',status='unknown')
open(13,file='pw.out',status='unknown')
open(14,file='abs.out',status='unknown')

!-- read input (current density)
!read(10,*)
!read(10,*) nt,ndata
ndata=3
read(10,*) time_delta,epsilon
read(10,*) nt
allocate(t(nt),c(nt,3),d(nt,3)); t=0.0d0; c=0.0d0; d=0.0d0
do i=1,nt
  time=(i-1)*time_delta
  t(i)=time
  read(10,*) ctmp2,(c(i,id),id=1,3)
enddo

!-- write input (current density)
write(11,'(''time(au) total current(x,y,z) (au)'')')
write(11,'(2i8)') nt,ndata+1
do i=1,nt
  write(11,'(4f15.8)') t(i),(c(i,id),id=1,3)
enddo

!-- calculate dipole difference d(t)-d(0)
!! trapezoidal
do id=1,3
  d(1,id)=c(1,id)*(t(2)-t(1))*0.5d0
  do i=2,nt-1
    d(i,id)=d(i-1,id)+c(i,id)*(t(i+1)-t(i))
  enddo
  d(nt,id)=d(nt-1,id)+c(nt,id)*(t(nt)-t(nt-1))*0.5d0
enddo
!! simpson
!do id=1,3
!  d(1,id)=c(1,id)*time_delta/3.0d0
!  do i=2,nt-1
!    if(mod(i,2)==1) fac=2.0d0
!    if(mod(i,2)==0) fac=4.0d0
!    d(i,id)=d(i-1,id)+fac*c(i,id)*time_delta/3.0d0
!  enddo
!  d(nt,id)=d(nt-1,id)+c(nt,id)*time_delta/3.0d0
!enddo

!-- shift dipole difference
write(*,'(''# induced average dipole'')')
do id=1,3
tmp=sum(d(1:nt,id))/dble(nt)
write(*,'(''id='',i2,'' P-shift='',f15.8)') id,tmp
do i=1,nt
  d(i,id)=d(i,id)-tmp
enddo
enddo

!-- write dipole difference
write(12,'(''time(au) dipole difference(x,y,z) (au)'')')
write(12,'(2i8)') nt,ndata+1
do i=1,nt
  write(12,'(4f15.8)') t(i),(d(i,id),id=1,3)
enddo

!-- set frequency window, dumping factor = gamma
gamma=0.00001d0
omega_delta=0.0005d0
omega_max=4000
allocate(omega(0:omega_max)); omega=0.0d0
do j=0,omega_max
  omega(j)=j*omega_delta
enddo

!-- Fourier transform from TIME(t) to FREQUENCY(omega)
allocate(d_w(0:omega_max,3)); d_w=(0.0d0,0.0d0)
do id=1,3
  do j=0,omega_max
  do i=1,nt
!! P(t) -> P(w) : P(w)=a(w)*E(w)
    d_w(j,id)=d_w(j,id) &
  & +cdexp(zi*omega(j)*t(i))*dexp(-gamma*t(i)*t(i))*d(i,id)*time_delta
!----
!! J(t) -> J(w) * J(w)=a(w)*E(w)
!   d_w(j,id)=d_w(j,id) &
! & +cdexp(zi*omega(j)*t(i))*dexp(-gamma*t(i)*t(i))*(-c(i,id))*time_delta
  enddo
  enddo
enddo

!-- write Fourier transform of dipole difference
write(13,'(''omega(au) FT of dipole difference(x,y,z) (au)'')')
write(13,'(2i8)') omega_max+1,2*ndata+1
do i=0,omega_max
  write(13,'(7f15.8)') omega(i),(d_w(i,id),id=1,3)
enddo

!-- calculate dynamic polarizability(omega)
allocate(alpha(0:omega_max,3)); alpha=0.0d0
do id=1,3
  do j=0,omega_max
    alpha(j,id)=-d_w(j,id)/epsilon
  enddo
enddo

!-- calculate dipole strength(omega)
allocate(S(0:omega_max,3)); S=0.0d0
do id=1,3
  do j=0,omega_max
    S(j,id)=2.0d0/pi*omega(j)*alpha(j,id)
  enddo
enddo

!-- write alpha(omega) or S(omega)
write(14,'(''omega(au) dynamic polarizability(x,y,z)(au)'')')
write(14,'(i8,''  7'')') omega_max+1
do j=0,omega_max
! write(14,'(7f15.8)') omega(j),(alpha(j,id),id=1,3)
! write(14,'(7f15.8)') omega(j),(S(j,id),id=1,3)
  write(14,'(4f15.8)') omega(j),(imag(S(j,id)),id=1,3)
enddo

deallocate(t,c,d,alpha,S,omega,d_w)
close(10); close(11); close(12); close(13)
end
