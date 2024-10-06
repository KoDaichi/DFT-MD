!================================================
!  Software name : STM
!  Subroutine(s) : jrcat_r2ft_of_3D_forwafd, jrcat_r1ft_of_3D_back,
!                  jrcat_r2ft_of_3D_back, jrcat_r2ft_of_3D,
!                  jrcat_c2ft_of_3D_back, jrcat_c2ft_of_3D_forward,
!                  jrcat_c2ft_of_3D
!  Author(s)     : Takahiro Yamasaki and Koichi Kato (June 7, 2004)
!
!  Contact address :  IIS,The University of Tokyo RSS21 project
!  
!  "Multiscale Simulation System for Function Analysis of Nanomaterials"  
!
!================================================
!
!     The original version of this program "STM" was developed in the
!  period 1998-2001 at JRCAT by Koichi Kato (Toshiba Corporate
!  Research and Development Center) and Takahiro Yamasaki (Fujitsu 
!  Laboratories Ltd.) through a joint research between JRCAT and
!  Toshiba Corporate Research and Development Center and another joint
!  research between JRCAT and FUJITSU Laboratories Ltd.
!     Since 2003, this program has been tuned and revised as it
!  matches to the first-principles molecular dynamics program "PHASE"
!  as a part of the national project "Frontier Simulation Software for
!  Industrial Science (FSIS)", which is supported by the IT program of
!  the Ministry of Education, Culture, Sports, Science and Technology
!  (MEXT) of Japan.
!     Since 2006, this program set has been developed as a part of the
!  national project "Revolutionary Simulation Software (RSS21)", which
!  is supported by the next-generation IT program of MEXT of Japan.

      subroutine jrcat_r2ft_of_3D_forward(ca,cb,id,n1,n2,n3,cw2,cw3)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $
      implicit none

      integer id,n1,n2,n3
      real*8 ca (2,0:id/2-1,0:n2-1,0:n3-1)
      real*8 cb (0:id/2-1,0:n2-1,0:n3-1,2)
      complex*16 cw2(0:n2-1)
      complex*16 cw3(0:n3-1)

      integer  n1h, idh
      integer lp2,lp3,lp4,lp5,lp8,key

      n1h = n1/2
      idh = id/2

      key = -1

c**  cb(i1,i2,i3,2) <- CONJG(ca(2,i1,i2,i3))
      call trans_f7(ca,cb,idh,n2,n3,key)

c** Length:n3 CFFT
c**     cb -> ca
c**      X(r1,r2,r3) -> X(r1,r2,k3)

      call setlp(n3,lp2,lp3,lp4,lp5,lp8)
      call fsb_lp235(cb,ca,cw3,n3,idh*n2,lp2,lp3,lp4,lp5,lp8)

c** cb(i1,i3,i2,2) <- ca(i1,i2,i3,2)
c**      X(r1,r2,k3) -> X(r1,k3,r2)

      call trans_b1(ca,cb,idh,n2,n3)

c** Length:n2 CFFT
c**     cb -> ca

      call setlp(n2,lp2,lp3,lp4,lp5,lp8)
      call fsb_lp235(cb,ca,cw2,n2,idh*n3,lp2,lp3,lp4,lp5,lp8) 

c**  cb(2,i1,i2,i3) <- CONJG(ca(i1,i3,i2,2))
      call trans_f8(ca,cb,idh,n1h,n2,n3,key)

c**  cb(2,i1,i2,i3) -> ca(2,i1,i2,i3)
      call cp_back(cb,ca,idh,n2,n3)

      return
      end

      subroutine trans_f7(ca,cb,idh,n2,n3,key)

      implicit real*8(a-h,o-z)

      real*8 ca(2,0:idh-1,0:n2-1,0:n3-1)
      real*8 cb(  0:idh-1,0:n2-1,0:n3-1,2)
      integer key
c***
      real*8 parity

      if(key.eq.-1) then
         parity = -1.0d0
      else if(key.eq.1) then
         parity = +1.0d0
      end if

c* Packing from ca to cb
      do 10 i3=0,n3-1
      do 10 i2=0,n2-1
      do 10 i1=0,idh-1

         cb(i1,i2,i3,1) =  ca(1,i1,i2,i3) 
         cb(i1,i2,i3,2) =  ca(2,i1,i2,i3)*parity 

 10   continue

      return
      end

      subroutine trans_b1(ca,cb,mm1,n2,n3)

      implicit real*8(a-h,o-z)
      real*8 ca(0:mm1-1,0:n2-1,0:n3-1,2)
      real*8 cb(0:mm1-1,0:n3-1,0:n2-1,2)

      do 10 i3=0,n3-1
      do 10 i2=0,n2-1
      do 10 i1=0,mm1-1

         cb(i1,i3,i2,1) = ca(i1,i2,i3,1)
         cb(i1,i3,i2,2) = ca(i1,i2,i3,2)

  10  continue

      return
      end

      subroutine trans_f8(cb,ca,idh,n1h,n2,n3,key)

      implicit real*8(a-h,o-z)
      real*8 ca(2,0:idh-1,0:n2-1,0:n3-1 )
      real*8 cb(  0:idh-1,0:n3-1,0:n2-1 ,2)
      integer n1h,n2,n3,key, nm
c*
      real*8 parity
      if(key.eq.+1) then
         parity = +1.0d0
         nm = idh - 1
      else if(key.eq.-1) then
         parity = -1.0d0
         nm = n1h
      end if

      do 10 i3=0,n3-1
      do 10 i2=0,n2-1
      do 10 i1=0,nm

         ca(1,i1,i2,i3) =   cb(i1,i3,i2,1)
         ca(2,i1,i2,i3) =   cb(i1,i3,i2,2)*parity

  10  continue

      if(key .eq. -1) then
         do i3=0,n3-1
            do i2=0, n2-1
               do i1=n1h+1,idh-1
                  ca(1,i1,i2,i3) = 0.d0
                  ca(2,i1,i2,i3) = 0.d0
               enddo
            enddo
         enddo
      end if

      return
      end

      subroutine jrcat_r1ft_of_3D_back(ca,cb,id,n1,n2,n3,cw1)
c                           @(#)jrcat_ft_stm.f 1.2 99/01/16 23:56:50
      implicit none

      integer id,n1,n2,n3
      real*8 ca (2, 0:id/2-1, 0:n2-1, 0:n3-1)
      real*8 cb (0:id/2-1,0:n2-1,0:n3-1,2)
      complex*16 cw1(0:n1/2-1,2)

      integer n1h,idh,mm1,kc1
      integer lp2,lp3,lp4,lp5,lp8
      
      n1h = n1/2
      idh = id/2
      mm1 = n1h +1
      kc1 = 0

c** ca(2,i1,i2,i3) -> cb(i3,i2,i1,2)
      call trans_in_b6(ca,cb,idh,n2,n3,n1h)

c**
      call trans_in_b3(cb,ca,cw1(0,2),n3*n2,n1h,kc1)

c** Length:n1h RCFFT  

      call setlp(n1h,lp2,lp3,lp4,lp5,lp8)
      call fsb_lp235(ca,cb,cw1,n1h,n3*n2,lp2,lp3,lp4,lp5,lp8)

c** ca(2,i1,i2,i3) <- cb(i3,i2,i1,2)

      call trans_in_b4(cb,ca,n3,n2,idh,n1h)

      return

      end  

      subroutine jrcat_r2ft_of_3D_back(ca,cb,id,n1,n2,n3,cw2,cw3)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $
      implicit none

      integer id,n1,n2,n3
      real*8 ca (2,0:id/2-1,0:n2-1,0:n3-1)
      real*8 cb (0:id/2-1,0:n2-1,0:n3-1,2)
      complex*16 cw2(0:n2-1)
      complex*16 cw3(0:n3-1)

      integer  n1h, idh
      integer lp2,lp3,lp4,lp5,lp8,key

      n1h = n1/2
      idh = id/2

      key = +1
c** Packin from ca(2,0:idh-1,0:n2-1,0:n3-1)
c**                         to cb(0:idh-1,0:n2-1,0:n3-1,2)
      call trans_f7(ca,cb,idh,n2,n3,key)

c** Length:n3 CFFT
c**     cb -> ca
c**      X(r1,r2,r3) -> X(r1,r2,k3)

      call setlp(n3,lp2,lp3,lp4,lp5,lp8)
      call fsb_lp235(cb,ca,cw3,n3,idh*n2,lp2,lp3,lp4,lp5,lp8)

c** cb(i1,i3,i2,2) <- ca(i1,i2,i3,2)
c**      X(r1,r2,k3) -> X(r1,k3,r2)

      call trans_b1(ca,cb,idh,n2,n3)

c** Length:n2 CFFT
c**     cb -> ca
c**      X(r1,k3,r2) -> X(r1,k3,k2)

      call setlp(n2,lp2,lp3,lp4,lp5,lp8)
      call fsb_lp235(cb,ca,cw2,n2,idh*n3,lp2,lp3,lp4,lp5,lp8) 

c**  cb(2,i1,i2,i3) <- ca(i1,i3,i2,2)
      call trans_f8(ca,cb,idh,n1h,n2,n3,key)

c**  cb(2,i1,i2,i3) -> ca(2,i1,i2,i3)
      call cp_back(cb,ca,idh,n2,n3)

      return
      end

      subroutine jrcat_r2ft_of_3D(ca,cb,id,n1,n2,n3,cw2,cw3,key)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $
c* * key == +1 (backward FFT)
c* * key == -1 (forward FFT)
      implicit none
      integer id,n1,n2,n3,key
      real*8 ca (2,0:id/2-1,0:n2-1,0:n3-1)
      real*8 cb (0:id/2-1,0:n2-1,0:n3-1,2)
      complex*16 cw2(0:n2-1)
      complex*16 cw3(0:n3-1)
      integer  n1h, idh
      integer lp2,lp3,lp4,lp5,lp8

      n1h = n1/2
      idh = id/2

c** (key == +1)
c**  Packin from ca(2,0:idh-1,0:n2-1,0:n3-1)
c**                         to cb(0:idh-1,0:n2-1,0:n3-1,2)
c** (key == -1)
c**  cb(i1,i2,i3,2) <- CONJG(ca(2,i1,i2,i3))
      call trans_f7(ca,cb,idh,n2,n3,key)
c** Length:n3 CFFT
c**     cb -> ca
c**      X(r1,r2,r3) -> X(r1,r2,k3)

      call setlp(n3,lp2,lp3,lp4,lp5,lp8)
      call fsb_lp235(cb,ca,cw3,n3,idh*n2,lp2,lp3,lp4,lp5,lp8)

c** cb(i1,i3,i2,2) <- ca(i1,i2,i3,2)
c**      X(r1,r2,k3) -> X(r1,k3,r2)

      call trans_b1(ca,cb,idh,n2,n3)

c** Length:n2 CFFT
c**     cb -> ca
c**      X(r1,k3,r2) -> X(r1,k3,k2)

      call setlp(n2,lp2,lp3,lp4,lp5,lp8)
      call fsb_lp235(cb,ca,cw2,n2,idh*n3,lp2,lp3,lp4,lp5,lp8) 

c** (key == +1)  cb(2,i1,i2,i3) <- ca(i1,i3,i2,2)
c** (key == -1)  cb(2,i1,i2,i3) <- CONJG(ca(i1,i3,i2,2))
      call trans_f8(ca,cb,idh,n1h,n2,n3,key)

c**  cb(2,i1,i2,i3) -> ca(2,i1,i2,i3)
      call cp_back(cb,ca,idh,n2,n3)

      return
      end

      subroutine cp_back(cb,ca,idh,n2,n3)
      implicit real*8 (a-h,o-z)
      
      real*8 ca(2,0:idh-1,0:n2-1, 0:n3-1)
      real*8 cb(2,0:idh-1,0:n2-1, 0:n3-1)

      do i3 = 0, n3-1
        do i2 = 0, n2-1
          do i1 = 0, idh-1
            ca(1,i1,i2,i3) = cb(1,i1,i2,i3)
            ca(2,i1,i2,i3) = cb(2,i1,i2,i3)
          end do
        end do
      end do

      return
      end

      subroutine trans_b5(cb,ca,idh,n2,n3)

      implicit real*8(a-h,o-z)

      real*8 ca(2,0:idh-1, 0:n2-1, 0:n3-1)
      real*8 cb(  0:idh-1, 0:n3-1, 0:n2-1 ,2)

      do i2 = 0, n2-1
      do i3 = 0, n3-1
      do i1 = 0, idh-1

        ca(1,i1,i2,i3) = cb(i1,i3,i2,1)
        ca(2,i1,i2,i3) = cb(i1,i3,i2,2)

      end do
      end do
      end do

      return
      end 

      subroutine trans_in_b6(ca,cb,idh,n2,n3,n1h)

      implicit real*8(a-h,o-z)

      real*8 ca(2,0:idh-1,0:n2-1, 0:n3-1)
      real*8 cb(0:n3-1, 0:n2-1, 0:n1h ,2)

      do i1 = 0, n1h-1
      do i2 = 0, n2-1
      do i3 = 0, n3-1

        cb(i3,i2,i1,1) = ca(1,i1,i2,i3)
        cb(i3,i2,i1,2) = ca(2,i1,i2,i3)

      end do
      end do
      end do

      return
      end
c --------------------------------------------------------
      subroutine jrcat_c2ft_of_3D_back(ca,cb,id,n1,n2,n3,cw2,cw3)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $
      implicit none
      integer id,n1,n2,n3
      real*8 ca (0:id-1,0:n2-1,0:n3-1,2)
      real*8 cb (0:id-1,0:n2-1,0:n3-1,2)
      complex*16 cw2(0:n2-1)
      complex*16 cw3(0:n3-1)

      integer key
      key = -1
      
      call jrcat_c2ft_of_3D(ca,cb,id,n1,n2,n3,cw2,cw3,key)
      return
      end
c --------------------------------------------------------
      subroutine jrcat_c2ft_of_3D_forward(ca,cb,id,n1,n2,n3,cw2,cw3)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $
      implicit none
      integer id,n1,n2,n3
      real*8 ca (0:id-1,0:n2-1,0:n3-1,2)
      real*8 cb (0:id-1,0:n2-1,0:n3-1,2)
      complex*16 cw2(0:n2-1)
      complex*16 cw3(0:n3-1)

      integer key
      key = +1
      
      call jrcat_c2ft_of_3D(ca,cb,id,n1,n2,n3,cw2,cw3,key)
      return
      end


      subroutine trans_bc1(ca,cb,n1,n2,n3)

      implicit real*8(a-h,o-z)
      real*8 ca(0:n1-1,0:n2-1,0:n3-1,2)
      real*8 cb(0:n1-1,0:n3-1,0:n2-1,2)

      do 10 i3=0,n3-1
      do 10 i2=0,n2-1
      do 10 i1=0,n1-1

         cb(i1,i3,i2,1) = ca(i1,i2,i3,1)
         cb(i1,i3,i2,2) = ca(i1,i2,i3,2)

  10  continue

      return
      end
c --------------------------------------------------------
      subroutine jrcat_c2ft_of_3D(ca,cb,id,n1,n2,n3,cw2,cw3,key)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $
      implicit none

      integer id,n1,n2,n3,key
      real*8 ca (0:id-1,0:n2-1,0:n3-1,2)
      real*8 cb (0:id-1,0:n2-1,0:n3-1,2)
      complex*16 cw2(0:n2-1)
      complex*16 cw3(0:n3-1)

      integer lp2,lp3,lp4,lp5,lp8

c** Packin from ca(2,0:id-1,0:n2-1,0:n3-1)
c**                         to cb(0:id-1,0:n2-1,0:n3-1,2)
c**  if(key == -1)
c**    cb(i1,i2,i3,2) <- CONJG(ca(2,i1,i2,i3)))
c**  else if(key == +1)
c**    cb(i1,i2,i3,2) <- ca(2,i1,i2,i3)
      call trans_fb0(ca,cb,id,n1,n2,n3,key)

c** Length:n3 CFFT
c**     cb -> ca
c**      X(r1,r2,r3) -> X(r1,r2,k3)

      call setlp(n3,lp2,lp3,lp4,lp5,lp8)
      call fsb_lp235(cb,ca,cw3,n3,n1*n2,lp2,lp3,lp4,lp5,lp8)

c** cb(i1,i3,i2,2) <- ca(i1,i2,i3,2)
c**      X(r1,r2,k3) -> X(r1,k3,r2)

      call trans_bc1(ca,cb,n1,n2,n3)

c** Length:n2 CFFT
c**     cb -> ca
c**      X(r1,k3,r2) -> X(r1,k3,k2)

      call setlp(n2,lp2,lp3,lp4,lp5,lp8)
      call fsb_lp235(cb,ca,cw2,n2,n1*n3,lp2,lp3,lp4,lp5,lp8) 
c**
c**  if(key == -1)
c**    cb(i1,i2,i3,2) <- CONJG(ca(i1,i3,i2,2))
c**  else if(key == +1)
c**    cb(i1,i2,i3,2) <- ca(i1,i3,i2,2)
c**  cb(2,i1,i2,i3) <- ca(i1,i3,i2,2)
      call trans_bc5(ca,cb,id,n1,n2,n3,key)

c**  cb(2,i1,i2,i3) -> ca(2,i1,i2,i3)
      call cp_back(cb,ca,id,n2,n3)

      return
      end

      subroutine trans_bc5(ca,cb,id,n1,n2,n3,key)
      implicit none

      integer id,n1,n2,n3,key
      real*8 ca(0:n1-1, 0:n3-1, 0:n2-1, 2)
      real*8 cb(2, 0:id-1, 0:n2-1, 0:n3-1)
c
      integer i3,i2,i1
      real*8 parity

      if(key.eq.+1) then
         parity = +1.0d0
      else if(key.eq.-1) then
         parity = -1.0d0
      end if

      do i3=0,n3-1
         do i2=0,n2-1
            do i1=0,id-1
               cb(1,i1,i2,i3) = 0.d0
               cb(2,i1,i2,i3) = 0.d0
            end do
         end do
      end do

      do i2 = 0, n2-1
      do i3 = 0, n3-1
      do i1 = 0, n1-1

        cb(1,i1,i2,i3) = ca(i1,i3,i2,1)
        cb(2,i1,i2,i3) = ca(i1,i3,i2,2)*parity

      end do
      end do
      end do

      return
      end 

      subroutine trans_fb0(ca,cb,id,n1,n2,n3,key)
      implicit none
      integer id,n1,n2,n3,key
      real*8 ca(2,id,n2*n3)
      real*8 cb(n1,n2*n3,2)
c***
      integer i23,i1
      real*8 parity

      if(key.eq.+1) then
         parity = +1.0d0
      else if(key.eq.-1) then
         parity = -1.0d0
      end if

      do i23 = 1,n2*n3
         do i1 = 1,n1

            cb(i1,i23,1) =   ca(1,i1,i23) 
            cb(i1,i23,2) =   ca(2,i1,i23)*parity

         end do
      end do

      return
      end
c --------------------------------------------------------
      subroutine jrcat_r3ft(ra,work,id,n1,n2,n3,cw1,cw2,cw3,
     &                      kc1,kc2,kc3,key)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $
c**
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)
c
c------------------>      for Estimation of CPU cost.
c                             by T.Yamasaki
c                               22th May. 1996
      real*8 CRTVL, PCPUDF
      parameter (CRTVL = 1.d-5, PCPUDF = 2.0)
c        CRTVL : CRItical VaLue for Division.
c        PCPUDF: Percent of CPU DiFference.
      real*8 cpu0ol,cpudif,rcpudf,ecpu,tcpu
      integer icount
      data icount/0/
      data ecpu/0.d0/
      data cpu0ol/0.d0/
      real*8 t_start, t_end, UMICRO
      parameter (UMICRO = 1.d-6)
c <-------------------------------------------
      real*8 ra   (id,n2,n3)
      real*8 work (id,n2,n3)

      complex*16 cw1(0:n1-1)
      complex*16 cw2(0:n2-1)
      complex*16 cw3(0:n3-1)

      integer kc1,kc2,kc3
      integer nbank

c**********************************************************************
c*   key = 0 : set up of FFT.
c*             ( Calculation of cw1,cw2 and cw3. )
c**********************************************************************
c*   n1,n2 and n3 must include radix 2, and they may include 3 or 5.
c*   id must be larger than or equal to n1 + 2.
c*
c**********************************************************************
c*  key = -1 & -2 : F.T. of real 3d array ra(0:n1-1,0:n2-1,0:n3-1) .
c**********************************************************************
c*
c*   The result(complex 'cb') is overwritten on 'ra'.
c*   ( i.e. equivalence(ra,cb) ).
c*   
c*     cb(k1,k2,k3) 
c*   = sum(ix=0:n1-1,iy=0:n2-1,iz=0:n3-1) 
c*     ra(i1,i2,i3) exp [- 2 * pai * ci *
c*                                 ( k1*i1/n1 + k2*i2/n2 + k3*i3/n3 ) ]
c*     
c*     where 0 <= k1 <= n1/2, 0 <= k2 <= n2-1, 0 <= k3 <= n3-1.
c*
c**********************************************************************
c*  key = +1 & +2: Inverse F.T. of copmplex 3d array cb(0:n1/2,0:n2-1,
c*                                                        0:n3-1) .
c**********************************************************************
c*   Transform from 'cb' to 'ra'.
c*
c**********************************************************************
c*  kc1,kc2 kc3 are integer parameters specifying the mask for
c*  input data or output data.
c*
c*     If kc1.ge.n1/2.or.kc1.lt.0, then kc1 is forced to be 0. 
c*     If kc2.ge.n2/2.or.kc2.lt.0, then kc2 is forced to be 0. 
c*     If kc3.ge.n3/2.or.kc3.lt.0, then kc3 is forced to be 0. 
c*
c**********************************************************************
c** REMARK on Mask.
c**********************************************************************
c*
c*  key = -1 : Mask is for Output data ( cb ),
c*  key = +1 : Mask is for Input  data ( cb ),
c*  key = -2 : Mask is for Input  data ( ca ),
c*  key = +2 : Mask is for Output data ( ca ),
c*
c*     where ca(0:idh/2,0:n2-1,0:n3-1) and equivalence(ra,ca).
c*
c**********************************************************************
c*     Mask for Output data:
c**********************************************************************
c*
c*     If kc1.eq.0.and.kc2.eq.0.and.kc3.eq.0, then the full wavenumber
c*     range is calculated.
c*
c*     If kc1.ne.0, then only "0 <= k1 <= kc1" range is calculated.
c*
c*     If kc2.ne.0, then only "0 <= k2 <= kc2" and "n2-kc2 <= k2<= n2-1" 
c*     range is calculated.
c*
c*     If kc3.ne.0, then only "0 <= k3 <= kc3" and "n3-kc3 <= k3<= n3-1" 
c*     range is calculated.
c*
c**********************************************************************
c*     Mask for Input data:
c**********************************************************************
c*
c*     If kc1.eq.0.and.kc2.eq.0.and.kc3.eq.0, then the full wavenumber
c*     range is used to obtain real array 'ra'.
c*
c*     If kc1.ne.0, then only "0 <= k1 <= kc1" range is used.
c*
c*     If kc2.ne.0, then only "0 <= k2 <= kc2" and "n2-kc2 <= k2<= n2-1" 
c*     range is used.
c*
c*     If kc3.ne.0, then only "0 <= k3 <= kc3" and "n3-kc3 <= k3<= n3-1" 
c*     range is used.
c*
c**********************************************************************
c*
c******************************************************************
c**   choose nbank = 0 or nbank = 1 
c******************************************************************
c*   nbank = 0:
c*      The transpose is normally implemented.
c*      Bank conflict may appear patriculary in Vector Computers.
c*   nbank = 1:
c*      The transpose is implemented so as to prevent
c*      the bank conflict in the case the conflict may appear.
c******************************************************************
c
      !call gettod(t_start)

      nbank = 0

c**

      if(key.eq.0) then

c** set up of FFT.

         call r3ft_0(cw1,cw2,cw3,n1,n2,n3)

         goto 1001
c$$$         return

      end if

c**
      if(kc1.ge.n1/2.or.kc1.lt.0) then 
         write(6,*) 'warning in FFT: kc1 is irrelevant.'
         write(6,*) 'kc1 has changed to be 0.'
         kc1 = 0
      end if

      if(kc2.ge.n2/2.or.kc2.lt.0) then 
         write(6,*) 'warning in FFT: kc2 is irrelevant.'
         write(6,*) 'kc2 has changed to be 0.'
         kc2 = 0
      end if

      if(kc3.ge.n3/2.or.kc3.lt.0) then 
         write(6,*) 'warning in FFT: kc3 is irrelevant.'
         write(6,*) 'kc3 has changed to be 0.'
         kc3 = 0
      end if
c**

      if(key.eq.-1) then

c** Forword FFT with mask for output.

         call r3ft_f_mask_for_out(ra,work,id,n1,n2,n3,cw1,cw2,cw3,
     &                     kc1,kc2,kc3,nbank)

      else if(key.eq.-2) then

c** Forword FFT with mask for input.

         call r3ft_f_mask_for_in (ra,work,id,n1,n2,n3,cw1,cw2,cw3,
     &                     kc1,kc2,kc3,nbank)

      else if(key.eq.+1) then

c** Backword FFT with mask for input.

         call r3ft_b_mask_for_in (ra,work,id,n1,n2,n3,cw1,cw2,cw3,
     &                     kc1,kc2,kc3,nbank)

      else if(key.eq.+2) then

c** Backword FFT with mask for output.

         call r3ft_b_mask_for_out(ra,work,id,n1,n2,n3,cw1,cw2,cw3,
     &                     kc1,kc2,kc3,nbank)

      end if

c**  
 1001 continue

c -------->  CPU cost estimation
      !call gettod(t_end)
      tcpu = (t_end - t_start)*UMICRO
      icount = icount + 1
      ecpu = ecpu + tcpu
c <-------
      
      return

      entry jrc_r3strt
      icount = 0
      ecpu = 0.d0
      return
c--------------------------
      entry jrc_r3end
      cpudif = dabs(cpu0ol - ecpu)
      if(ecpu.gt.CRTVL) then
         rcpudf = cpudif/ecpu * 100.0
      else if(cpu0ol.le.CRTVL) then
         rcpudf = 0.0
      else
         rcpudf = 100.0
      endif
c               <-- Ratio in percent of CPU time difference
c                 between previous one and present one.
      cpu0ol = ecpu

c$$$      call eqivvl(rcpudf)

      if(rcpudf.gt.PCPUDF) then
         if(icount.gt.0) then
            tcpu = ecpu/dfloat(icount)
         endif
         !!write(6,9001) ecpu, icount, tcpu
 9001    format(1H ,' <<CPU TIME for JRCAT_R3FT = ',f11.3
     &        ,'(Sec.) = (',i6,' *',f10.5,')>>')
      endif
      return
      end  

      subroutine r3ft_0(cw1,cw2,cw3,n1,n2,n3)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $
c**
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)

      complex*16 cw1(0:n1/2-1,2)
      complex*16 cw2(0:n2-1)
      complex*16 cw3(0:n3-1)

c**
      if(mod(n1,2).ne.0) then
         stop 'The first length is not an even number.'
      end if

         n1h = n1/2

         pai=3.14159265358979323846d0
         pc =2.0d0*pai/dble(n1)
         ci=(0.0d0,1.0d0)

         call setcw(cw1,n1h)
         call setcw(cw2,n2 )
         call setcw(cw3,n3 ) 

      do i=0,n1h-1
         tr=dcos(pc*dble(i))
         ti=dsin(pc*dble(i))
         cw1(i,2)=ci*dcmplx(tr,ti)
      end do   

      return
      end  

      subroutine r3ft_b_mask_for_in(ca,cb,id,n1,n2,n3,cw1,cw2,cw3,
     &                        kc1,kc2,kc3,nbank)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $
c**
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)

      real*8 ca (2,0:id/2-1,0:n2-1,0:n3-1)
      real*8 cb (0:id/2-1,0:n2-1,0:n3-1,2)
c     complex*16 ca (0:id/2-1,0:n2-1,0:n3-1)
c     complex*16 cb (0:id/2-1,0:n2-1,0:n3-1)

      complex*16 cw1(0:n1/2-1,2)
      complex*16 cw2(0:n2-1)
      complex*16 cw3(0:n3-1)

      integer kc1,kc2,kc3
      integer nbank
c**
         n1h = n1/2
         idh = id/2

      if(kc1.eq.0) then
           mm1 = n1h+1
c          mm1 = idh
      else
           mm1 = kc1+1
      end if

      if(kc2.eq.0) then
           mm2 = n2
      else
           mm2 = kc2*2+1
      end if

c** Packin from ca(2,0:idh-1,0:n2-1,0:n3-1)
c**        to   cb(0:mm1-1,0:mm2-1,0:n3-1,2)

      call trans_in_b0(ca,cb,idh,n2,n3,mm1,mm2,kc1,kc2,kc3)

c** Length:n3 CFFT  

      call setlp(n3,lp2,lp3,lp4,lp5,lp8)
      call fsb_lp235(cb,ca,cw3,n3,mm1*mm2,lp2,lp3,lp4,lp5,lp8)

c** cb(i1,i3,i2,2) <- ca(i1,i2,i3,2)

      call trans_in_b1(ca,cb,mm1,n2,mm2,n3,kc2)

c** Length:n2 CFFT    

      call setlp(n2,lp2,lp3,lp4,lp5,lp8)
      call fsb_lp235(cb,ca,cw2,n2,mm1*n3,lp2,lp3,lp4,lp5,lp8) 

c** cb(i3,i2,i1,2) <- ca(i1,i3,i2,2)

      call trans_in_b2(ca,cb,mm1,n2,n3,n1h)

c**
      call trans_in_b3(cb,ca,cw1(0,2),n3*n2,n1h,kc1)

c** Length:n1h RCFFT  

      call setlp(n1h,lp2,lp3,lp4,lp5,lp8)
      call fsb_lp235(ca,cb,cw1,n1h,n3*n2,lp2,lp3,lp4,lp5,lp8)

c** ca(2,i1,i2,i3) <- cb(i3,i2,i1,2)

      call trans_in_b4(cb,ca,n3,n2,idh,n1h)

      return

      end  

      subroutine r3ft_b_mask_for_out(ca,cb,id,n1,n2,n3,
     &                        cw1,cw2,cw3,
     &                        kc1,kc2,kc3,nbank)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $
c**
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)

      real*8 ca (2,0:id/2-1,0:n2-1,0:n3-1)
      real*8 cb (0:id/2-1,0:n2-1,0:n3-1,2)
c     complex*16 ca (0:id/2-1,0:n2-1,0:n3-1)
c     complex*16 cb (0:id/2-1,0:n2-1,0:n3-1)

      complex*16 cw1(0:n1/2-1,2)
      complex*16 cw2(0:n2-1)
      complex*16 cw3(0:n3-1)

      integer kc1,kc2,kc3
      integer nbank
c**
         n1h = n1/2
         idh = id/2

      if(kc1.eq.0) then
           mm1 = n1h+1 
      else
           mm1 = kc1+1
      end if

      if(kc2.eq.0) then
           mm2 = n2
      else
           mm2 = kc2*2+1
      end if

      if(kc3.eq.0) then
           mm3 = n3
      else
           mm3 = kc3*2+1
      end if

c** cb(i1,i2,i3,2) <- ca(2,i1,i2,i3)

      call trans_out_b0(ca,cb,idh,n2,n3,n1h)

c** Length:n3 CFFT  

      call setlp(n3,lp2,lp3,lp4,lp5,lp8)
      call fsb_lp235(cb,ca,cw3,n3,(n1h+1)*n2,lp2,lp3,lp4,lp5,lp8)

c** cb(i1,i3,i2,2) <- ca(i1,i2,i3,2)

      call trans_out_b1(ca,cb,n1h,n2,n3,mm3,kc3)

c** Length:n2 CFFT    

      call setlp(n2,lp2,lp3,lp4,lp5,lp8)
      call fsb_lp235(cb,ca,cw2,n2,(n1h+1)*mm3,lp2,lp3,lp4,lp5,lp8) 

c** cb(i3,i2,i1,2) <- ca(i1,i3,i2,2)

      call trans_out_b2(ca,cb,n1h,n2,mm2,mm3,kc2)

c**
      call trans_out_b3(cb,ca,cw1(0,2),mm3*mm2,n1h)

c** Length:n1h RCFFT  

      call setlp(n1h,lp2,lp3,lp4,lp5,lp8)
      call fsb_lp235(ca,cb,cw1,n1h,mm3*mm2,lp2,lp3,lp4,lp5,lp8)

c** ca(2,i1,i2,i3) <- cb(i3,i2,i1,2)

      call trans_out_b4(cb,ca,idh,n1h,n2,n3,mm2,mm3,kc1,kc2,kc3)

      return

      end  

      subroutine r3ft_f_mask_for_out(ca,cb,id,n1,n2,n3,
     &                        cw1,cw2,cw3,
     &                        kc1,kc2,kc3,nbank)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $
c**
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)

      real*8 ca (2,0:id/2-1,0:n2-1,0:n3-1)
      real*8 cb (0:id/2-1,0:n2-1,0:n3-1,2)
c     complex*16 ca (0:id/2-1,0:n2-1,0:n3-1)
c     complex*16 cb (0:id/2-1,0:n2-1,0:n3-1)

      complex*16 cw1(0:n1/2-1,2)
      complex*16 cw2(0:n2-1)
      complex*16 cw3(0:n3-1)

      integer kc1,kc2,kc3
      integer nbank
c**
         n1h = n1/2
         idh = id/2

      if(kc1.eq.0) then
           mm1 = n1h+1
      else
           mm1 = kc1+1
      end if

      if(kc2.eq.0) then
           mm2 = n2
      else
           mm2 = kc2*2+1
      end if
    
c**   cb(i3,i2,i1,2) <- CONJG(ca(2,i1,i2,i3))

         call trans_out_f1(ca,cb,idh,n1h,n2,n3)
 
c** Length:n1h RFFT

      call setlp(n1h,lp2,lp3,lp4,lp5,lp8)

      call fsb_lp235(cb,ca,cw1,n1h,n3*n2,lp2,lp3,lp4,lp5,lp8)
      
c** Transform from 'n1h' complex to 'n1' real array. 

         call trans_out_f2(ca,cb,cw1(0,2),n3,n2,n1h,kc1) 

c**   ca(i3,k1,i2,2) <- cb(i3,i2,k1,2)

         call trans_out_f3(cb,ca,n3,n2,mm1,n1h)

c** Length:n2 CFFT  

      call setlp(n2,lp2,lp3,lp4,lp5,lp8)

      call fsb_lp235(ca,cb,cw2,n2,n3*mm1,lp2,lp3,lp4,lp5,lp8)

c**   ca(k1,k2,i3,2) <- cb(i3,k1,k2,2)

            call trans_out_f4(cb,ca,n3,mm1,n2,mm2,kc2)

c** Length:n3 CFFT  

      call setlp(n3,lp2,lp3,lp4,lp5,lp8)

      call fsb_lp235(ca,cb,cw3,n3,mm1*mm2,lp2,lp3,lp4,lp5,lp8)

c**   ca(2,i1,i2,i3) <- CONJG(cb(i1,i2,i3,2))

         call trans_out_f5(cb,ca,idh,n1h,n2,n3,mm1,mm2,kc1,kc2,kc3)

      return

      end  

      subroutine r3ft_f_mask_for_in(ca,cb,id,n1,n2,n3,cw1,cw2,cw3,
     &                        kc1,kc2,kc3,nbank)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $
c**
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)

      real*8 ca (2,0:id/2-1,0:n2-1,0:n3-1)
      real*8 cb (0:id/2-1,0:n2-1,0:n3-1,2)

c     complex*16 ca (0:id/2-1,0:n2-1,0:n3-1)
c     complex*16 cb (0:id/2-1,0:n2-1,0:n3-1)

      complex*16 cw1(0:n1/2-1,2)
      complex*16 cw2(0:n2-1)
      complex*16 cw3(0:n3-1)

      integer kc1,kc2,kc3
      integer nbank
c**
         n1h = n1/2
         idh = id/2

      if(kc1.eq.0) then
           mm1 = n1h+1 
      else
           mm1 = kc1+1
      end if

      if(kc2.eq.0) then
           mm2 = n2
      else
           mm2 = kc2*2+1
      end if

      if(kc3.eq.0) then
           mm3 = n3
      else
           mm3 = kc3*2+1
      end if
    
c**   cb(i3,i2,i1,2) <- CONJG(ca(2,i1,i2,i3))

         call trans_in_f1(ca,cb,idh,n1h,n2,n3,mm2,mm3,kc1,kc2,kc3)

c** Length:n1h RFFT

      call setlp(n1h,lp2,lp3,lp4,lp5,lp8)

      call fsb_lp235(cb,ca,cw1,n1h,mm3*mm2,lp2,lp3,lp4,lp5,lp8)
      
c** Transform from 'n1h' complex to 'n1' real array. 

         call trans_in_f2(ca,cb,cw1(0,2),mm3,mm2,n1h) 

c**   ca(i3,k1,i2,2) <- cb(i3,i2,k1,2)

         call trans_in_f3(cb,ca,mm3,n2,mm2,n1h,kc2)

c** Length:n2 CFFT  

      call setlp(n2,lp2,lp3,lp4,lp5,lp8)

      call fsb_lp235(ca,cb,cw2,n2,mm3*(n1h+1),lp2,lp3,lp4,lp5,lp8)

c**   ca(k1,k2,i3,2) <- cb(i3,k1,k2,2)

            call trans_in_f4(cb,ca,n1h,n2,n3,mm3,kc3)

c** Length:n3 CFFT  

      call setlp(n3,lp2,lp3,lp4,lp5,lp8)

      call fsb_lp235(ca,cb,cw3,n3,(n1h+1)*n2,lp2,lp3,lp4,lp5,lp8)

c**   ca(2,i1,i2,i3) <- CONJG(cb(i1,i2,i3,2))

         call trans_in_f5(cb,ca,idh,n1h,n2,n3)

      return

      end  

      subroutine fsb_lp235(ca,cb,cw,n,ndim,lp2,lp3,lp4,lp5,lp8)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $
      implicit complex*16(a-h,o-z)

      real*8 ca(ndim*n*2),cb(ndim*n*2)
c     complex*16 ca(ndim,n),cb(ndim,n)

      complex*16 cw(0:n-1)

      k =1
      la=1

c*** RADIX 2 FFT.

      do 2 i=1,lp2
      if(k.eq.1) then
         call fsb_lp2(ca,cb,cw,la,n,ndim)
      else
         call fsb_lp2(cb,ca,cw,la,n,ndim)
      end if
         k=k*(-1)
         la=la*2
    2 continue

c*** RADIX 3 FFT.

      do 3 i=1,lp3
      if(k.eq.1) then
         call fsb_lp3(ca,cb,cw,la,n,ndim)
      else
         call fsb_lp3(cb,ca,cw,la,n,ndim)
      end if
         k=k*(-1)
         la=la*3
    3 continue

c*** RADIX 4 FFT.

      do 4 i=1,lp4
      if(k.eq.1) then
         call fsb_lp4(ca,cb,cw,la,n,ndim)
      else
         call fsb_lp4(cb,ca,cw,la,n,ndim)
      end if
         k=k*(-1)
         la=la*4
    4 continue

c*** RADIX 5 FFT.

      do 5 i=1,lp5
      if(k.eq.1) then
         call fsb_lp5(ca,cb,cw,la,n,ndim)
      else
         call fsb_lp5(cb,ca,cw,la,n,ndim)
      end if
         k=k*(-1)
         la=la*5
    5 continue

c*** RADIX 8 FFT.

      do 8 i=1,lp8
      if(k.eq.1) then
         call fsb_lp8(ca,cb,cw,la,n,ndim)
      else
         call fsb_lp8(cb,ca,cw,la,n,ndim)
      end if
         k=k*(-1)
         la=la*8
    8 continue
c***
      if(k.eq.+1) then

         do 100 ii=1,ndim*n*2
            cb(ii)=ca(ii)
  100    continue

      else

         do 101 ii=1,ndim*n*2
            ca(ii)=cb(ii)
  101    continue

      end if
c**
      return
      end
c---------------------------------------------
      subroutine fsb_lp2(ca,cb,cw,la,n,ndim)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $
c---------------------------------------------
      implicit real*8(a-h,o-z)

      real*8 ca(ndim,n,2),cb(ndim,n,2),cw(2,0:n-1)

      m=n/2
      np=(m-la)/la

      i1=m
      j1=la
      jump=la  

      do 20 l=1,la
      do 20 ii=1,ndim

         cb(ii,l   ,1)=(ca(ii,l,1)+ca(ii,l+i1,1))
         cb(ii,l   ,2)=(ca(ii,l,2)+ca(ii,l+i1,2))

         cb(ii,l+j1,1)=(ca(ii,l,1)-ca(ii,l+i1,1))
         cb(ii,l+j1,2)=(ca(ii,l,2)-ca(ii,l+i1,2))

   20 continue

      do 21 kp=1,np
      do 21 l=1,la
      do 21 ii=1,ndim

         i0=l+kp*la         
         j0=l+kp*(la+jump)         

         cb(ii,j0   ,1)=(ca(ii,i0,1)+ca(ii,i0+i1,1))
         cb(ii,j0   ,2)=(ca(ii,i0,2)+ca(ii,i0+i1,2))

         cc1r = ca(ii,i0,1)-ca(ii,i0+i1,1)
         cc1i = ca(ii,i0,2)-ca(ii,i0+i1,2)

         cb(ii,j0+j1,1)=cc1r*cw(1,kp*la) - cc1i*cw(2,kp*la)
         cb(ii,j0+j1,2)=cc1r*cw(2,kp*la) + cc1i*cw(1,kp*la)

   21 continue

      return
      end
c---------------------------------------------
      subroutine fsb_lp3(ca,cb,cw,la,n,ndim)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $
c---------------------------------------------
      implicit real*8(a-h,o-z)  

      real*8 ca(ndim,n,2),cb(ndim,n,2),cw(2,0:n-1)

      real*8 sin60

      m=n/3
      np=(m-la)/la

      i1=m
      i2=m*2 

      j1=la
      j2=la*2

      jump=la*2

      sin60=0.5d0*dsqrt(3.0d0)

c*    k=0 
      do 30 l=1,la
      do 30 ii=1,ndim

         ct1r=ca(ii,l+i1,1)+ca(ii,l+i2,1)
         ct1i=ca(ii,l+i1,2)+ca(ii,l+i2,2)

         ct2r=ca(ii,l   ,1)-0.5d0*ct1r       
         ct2i=ca(ii,l   ,2)-0.5d0*ct1i       

         ct3r= - sin60*(ca(ii,l+i1,2)-ca(ii,l+i2,2))
         ct3i=   sin60*(ca(ii,l+i1,1)-ca(ii,l+i2,1))

         cb(ii,l   ,1)=ca(ii,l,1)+ct1r
         cb(ii,l   ,2)=ca(ii,l,2)+ct1i

         cb(ii,l+j1,1)=(ct2r+ct3r)
         cb(ii,l+j1,2)=(ct2i+ct3i)

         cb(ii,l+j2,1)=(ct2r-ct3r)
         cb(ii,l+j2,2)=(ct2i-ct3i)

   30 continue

      do 31 kp=1,np
      do 31 l=1,la
      do 31 ii=1,ndim

         i0=l+kp*la         
         j0=l+kp*(la+jump)         

         ct1r=ca(ii,i0+i1,1)+ca(ii,i0+i2,1)
         ct1i=ca(ii,i0+i1,2)+ca(ii,i0+i2,2)

         ct2r=ca(ii,i0   ,1)-0.5d0*ct1r       
         ct2i=ca(ii,i0   ,2)-0.5d0*ct1i       

         ct3r= - sin60*(ca(ii,i0+i1,2)-ca(ii,i0+i2,2))
         ct3i=   sin60*(ca(ii,i0+i1,1)-ca(ii,i0+i2,1))

         cb(ii,j0   ,1)=ca(ii,i0,1)+ct1r
         cb(ii,j0   ,2)=ca(ii,i0,2)+ct1i

         ct4r = ct2r + ct3r
         ct4i = ct2i + ct3i

         ct5r = ct2r - ct3r
         ct5i = ct2i - ct3i

         cb(ii,j0+j1,1)=ct4r*cw(1,kp*j1) - ct4i*cw(2,kp*j1)
         cb(ii,j0+j1,2)=ct4r*cw(2,kp*j1) + ct4i*cw(1,kp*j1)

         cb(ii,j0+j2,1)=ct5r*cw(1,kp*j2) - ct5i*cw(2,kp*j2)
         cb(ii,j0+j2,2)=ct5r*cw(2,kp*j2) + ct5i*cw(1,kp*j2)

   31 continue

      return
      end

c---------------------------------------------
      subroutine fsb_lp4(ca,cb,cw,la,n,ndim)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $
c---------------------------------------------
      implicit real*8(a-h,o-z)

      real*8 ca(ndim,n,2),cb(ndim,n,2),cw(2,0:n-1)

         m   =n/4
         np  =(m-la)/la
         jump=la*3
   
      do 40 l=1,la
      do 40 ii=1,ndim

      cc0r = ca(ii,l  ,1)+ca(ii,l+2*M,1)
      cc0i = ca(ii,l  ,2)+ca(ii,l+2*M,2)

      cc1r = ca(ii,l+M,1)+ca(ii,l+3*M,1)
      cc1i = ca(ii,l+M,2)+ca(ii,l+3*M,2)

      cc2r = ca(ii,l  ,1)-ca(ii,l+2*M,1)
      cc2i = ca(ii,l  ,2)-ca(ii,l+2*M,2)

      cc3r = - ca(ii,l+M,2) + ca(ii,l+3*M,2)
      cc3i =   ca(ii,l+M,1) - ca(ii,l+3*M,1)

         cb(ii,l     ,1)= cc0r+cc1r
         cb(ii,l     ,2)= cc0i+cc1i

         cb(ii,l+  LA,1)= cc2r+cc3r
         cb(ii,l+  LA,2)= cc2i+cc3i

         cb(ii,l+2*LA,1)= cc0r-cc1r
         cb(ii,l+2*LA,2)= cc0i-cc1i

         cb(ii,l+3*LA,1)= cc2r-cc3r
         cb(ii,l+3*LA,2)= cc2i-cc3i

   40 continue
   
      DO 41 KP=1,NP
      DO 41 L=1,LA
      DO 41 ii=1,ndim

         I0=L+KP* LA
         J0=L+KP*(LA+JUMP)

      cc1r =  ca(ii,I0  ,1)+ca(ii,I0+2*M,1)
      cc1i =  ca(ii,I0  ,2)+ca(ii,I0+2*M,2)

      cc2r =  ca(ii,I0+M,1)+ca(ii,I0+3*M,1)
      cc2i =  ca(ii,I0+M,2)+ca(ii,I0+3*M,2)

      cc3r =  ca(ii,I0  ,1)-ca(ii,I0+2*M,1)
      cc3i =  ca(ii,I0  ,2)-ca(ii,I0+2*M,2)

      cc4r = -ca(ii,I0+M,2)+ca(ii,I0+3*M,2)
      cc4i =  ca(ii,I0+M,1)-ca(ii,I0+3*M,1)

         cc5r = cc3r + cc4r
         cc5i = cc3i + cc4i

         cc6r = cc1r - cc2r
         cc6i = cc1i - cc2i

         cc7r = cc3r - cc4r
         cc7i = cc3i - cc4i

         cb(ii,J0     ,1) = cc1r+cc2r
         cb(ii,J0     ,2) = cc1i+cc2i

         cb(ii,J0+  LA,1) = cc5r*cw(1,KP*LA  ) - cc5i*cw(2,KP*LA  )
         cb(ii,J0+  LA,2) = cc5r*cw(2,KP*LA  ) + cc5i*cw(1,KP*LA  )

         cb(ii,J0+2*LA,1) = cc6r*cw(1,KP*LA*2) - cc6i*cw(2,KP*LA*2)
         cb(ii,J0+2*LA,2) = cc6r*cw(2,KP*LA*2) + cc6i*cw(1,KP*LA*2)

         cb(ii,J0+3*LA,1) = cc7r*cw(1,KP*LA*3) - cc7i*cw(2,KP*LA*3)
         cb(ii,J0+3*LA,2) = cc7r*cw(2,KP*LA*3) + cc7i*cw(1,KP*LA*3)

   41 continue
         
      RETURN
      END

c---------------------------------------------
      subroutine fsb_lp5(ca,cb,cw,la,n,ndim)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $
c---------------------------------------------
      implicit real*8(a-h,o-z)

      real*8 ca(ndim,n,2),cb(ndim,n,2),cw(2,0:n-1)
      real*8 pai,sin36,sin72,r54

      m=n/5
      np  =(m-la)/la

      i1=m
      i2=m*2 
      i3=m*3 
      i4=m*4 

      j1=la
      j2=la*2
      j3=la*3
      j4=la*4

      jump=la*4

      pai  =3.14159265358979323846d0
      sin36=dsin(pai*0.2d0)
      sin72=dsin(pai*0.4d0)
      r54  =dsqrt(5.0d0)/4.0d0

c** k=0
      do 50 l=1,la
      do 50 ii=1,ndim

         ct1r=ca(ii,l+i1,1)+ca(ii,l+i4,1)
         ct1i=ca(ii,l+i1,2)+ca(ii,l+i4,2)

         ct2r=ca(ii,l+i2,1)+ca(ii,l+i3,1)
         ct2i=ca(ii,l+i2,2)+ca(ii,l+i3,2)

         ct3r=ca(ii,l+i1,1)-ca(ii,l+i4,1)
         ct3i=ca(ii,l+i1,2)-ca(ii,l+i4,2)

         ct4r=ca(ii,l+i2,1)-ca(ii,l+i3,1)
         ct4i=ca(ii,l+i2,2)-ca(ii,l+i3,2)

         ct5r=     ct1r+ct2r
         ct5i=     ct1i+ct2i

         ct6r=r54*(ct1r-ct2r)
         ct6i=r54*(ct1i-ct2i)

         ct7r=ca(ii,l   ,1)-0.25d0*ct5r
         ct7i=ca(ii,l   ,2)-0.25d0*ct5i

         ct8r=ct7r+ct6r
         ct8i=ct7i+ct6i

         ct9r=ct7r-ct6r
         ct9i=ct7i-ct6i

         ct10r=-(sin72*ct3i+sin36*ct4i)
         ct10i= (sin72*ct3r+sin36*ct4r)

         ct11r=-(sin36*ct3i-sin72*ct4i)
         ct11i= (sin36*ct3r-sin72*ct4r)

         cb(ii,l   ,1)=ca(ii,l,1)+ct5r
         cb(ii,l   ,2)=ca(ii,l,2)+ct5i

         cb(ii,l+j1,1)=(ct8r+ct10r)
         cb(ii,l+j1,2)=(ct8i+ct10i)

         cb(ii,l+j2,1)=(ct9r+ct11r)
         cb(ii,l+j2,2)=(ct9i+ct11i)

         cb(ii,l+j3,1)=(ct9r-ct11r)
         cb(ii,l+j3,2)=(ct9i-ct11i)

         cb(ii,l+j4,1)=(ct8r-ct10r)
         cb(ii,l+j4,2)=(ct8i-ct10i)

   50 continue

      do 51 kp=1,np
      do 51 l=1,la
      do 51 ii=1,ndim

         I0=L+KP* LA
         J0=L+KP*(LA+JUMP)

         ct1r=ca(ii,i0+i1,1)+ca(ii,i0+i4,1)
         ct1i=ca(ii,i0+i1,2)+ca(ii,i0+i4,2)

         ct2r=ca(ii,i0+i2,1)+ca(ii,i0+i3,1)
         ct2i=ca(ii,i0+i2,2)+ca(ii,i0+i3,2)

         ct3r=ca(ii,i0+i1,1)-ca(ii,i0+i4,1)
         ct3i=ca(ii,i0+i1,2)-ca(ii,i0+i4,2)

         ct4r=ca(ii,i0+i2,1)-ca(ii,i0+i3,1)
         ct4i=ca(ii,i0+i2,2)-ca(ii,i0+i3,2)

         ct5r=     ct1r+ct2r
         ct5i=     ct1i+ct2i

         ct6r=r54*(ct1r-ct2r)
         ct6i=r54*(ct1i-ct2i)

         ct7r=ca(ii,i0  ,1)-0.25d0*ct5r
         ct7i=ca(ii,i0  ,2)-0.25d0*ct5i

         ct8r=ct7r+ct6r
         ct8i=ct7i+ct6i

         ct9r=ct7r-ct6r
         ct9i=ct7i-ct6i

         ct10r=-(sin72*ct3i+sin36*ct4i)
         ct10i= (sin72*ct3r+sin36*ct4r)

         ct11r=-(sin36*ct3i-sin72*ct4i)
         ct11i= (sin36*ct3r-sin72*ct4r)

         ct12r = ct8r + ct10r
         ct12i = ct8i + ct10i

         ct13r = ct9r + ct11r
         ct13i = ct9i + ct11i

         ct14r = ct9r - ct11r
         ct14i = ct9i - ct11i

         ct15r = ct8r - ct10r
         ct15i = ct8i - ct10i

         cb(ii,j0   ,1)=ca(ii,i0,1)+ct5r
         cb(ii,j0   ,2)=ca(ii,i0,2)+ct5i

         cb(ii,j0+j1,1)=ct12r*cw(1,kp*j1)-ct12i*cw(2,kp*j1)
         cb(ii,j0+j1,2)=ct12r*cw(2,kp*j1)+ct12i*cw(1,kp*j1)

         cb(ii,j0+j2,1)=ct13r*cw(1,kp*j2)-ct13i*cw(2,kp*j2)
         cb(ii,j0+j2,2)=ct13r*cw(2,kp*j2)+ct13i*cw(1,kp*j2)

         cb(ii,j0+j3,1)=ct14r*cw(1,kp*j3)-ct14i*cw(2,kp*j3)
         cb(ii,j0+j3,2)=ct14r*cw(2,kp*j3)+ct14i*cw(1,kp*j3)

         cb(ii,j0+j4,1)=ct15r*cw(1,kp*j4)-ct15i*cw(2,kp*j4)
         cb(ii,j0+j4,2)=ct15r*cw(2,kp*j4)+ct15i*cw(1,kp*j4)

   51 continue

      return
      end
c---------------------------------------------
      subroutine fsb_lp8(ca,cb,cw,la,n,ndim)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $
c---------------------------------------------
      implicit real*8(a-h,o-z)

      real*8 CA(ndim,n,2),CB(ndim,n,2),CW(2,0:n-1)
   
         R2=0.5D0
         X=DSQRT(R2)

         M   =N/8
         NP  =(M-LA)/LA
         JUMP=LA*7

      DO 80 L=1,LA
      DO 80 ii=1,ndim

      CC0r = CA(ii,l    ,1)+CA(ii,l+4*M,1)
      CC0i = CA(ii,l    ,2)+CA(ii,l+4*M,2)

      CC1r = CA(ii,l+  M,1)+CA(ii,l+5*M,1)
      CC1i = CA(ii,l+  M,2)+CA(ii,l+5*M,2)

      CC2r = CA(ii,l+2*M,1)+CA(ii,l+6*M,1)
      CC2i = CA(ii,l+2*M,2)+CA(ii,l+6*M,2)

      CC3r = CA(ii,l+3*M,1)+CA(ii,l+7*M,1)
      CC3i = CA(ii,l+3*M,2)+CA(ii,l+7*M,2)

      CC4r = CA(ii,l    ,1)-CA(ii,l+4*M,1)
      CC4i = CA(ii,l    ,2)-CA(ii,l+4*M,2)

      CC5r = CA(ii,l+  M,1)-CA(ii,l+5*M,1)
      CC5i = CA(ii,l+  M,2)-CA(ii,l+5*M,2)

      CC6r = -CA(ii,l+2*M,2)+CA(ii,l+6*M,2)
      CC6i =  CA(ii,l+2*M,1)-CA(ii,l+6*M,1)

      CC7r = -CA(ii,l+3*M,2)+CA(ii,l+7*M,2)
      CC7i =  CA(ii,l+3*M,1)-CA(ii,l+7*M,1)

         CD0r = CC0r+CC2r
         CD0i = CC0i+CC2i

         CD1r = CC1r+CC3r
         CD1i = CC1i+CC3i

         CD2r = CC4r+CC6r
         CD2i = CC4i+CC6i

            CE1r = CC5r+CC7r
            CE1i = CC5i+CC7i

         CD3r = (CE1r-CE1i)*x
         CD3i = (CE1r+CE1i)*x

         CD4r = CC0r-CC2r
         CD4i = CC0i-CC2i

         CD5r =-(CC1i-CC3i)
         CD5i = (CC1r-CC3r)

         CD6r = CC4r-CC6r
         CD6i = CC4i-CC6i

            CE2r = CC5r-CC7r
            CE2i = CC5i-CC7i

         CD7r = (-CE2r-CE2i)*x
         CD7i = ( CE2r-CE2i)*x

      CB(ii,l     ,1) = CD0r+CD1r
      CB(ii,l     ,2) = CD0i+CD1i

      CB(ii,l+  LA,1) = CD2r+CD3r
      CB(ii,l+  LA,2) = CD2i+CD3i

      CB(ii,l+2*LA,1) = CD4r+CD5r
      CB(ii,l+2*LA,2) = CD4i+CD5i

      CB(ii,l+3*LA,1) = CD6r+CD7r
      CB(ii,l+3*LA,2) = CD6i+CD7i

      CB(ii,l+4*LA,1) = CD0r-CD1r
      CB(ii,l+4*LA,2) = CD0i-CD1i

      CB(ii,l+5*LA,1) = CD2r-CD3r
      CB(ii,l+5*LA,2) = CD2i-CD3i

      CB(ii,l+6*LA,1) = CD4r-CD5r
      CB(ii,l+6*LA,2) = CD4i-CD5i

      CB(ii,l+7*LA,1) = CD6r-CD7r
      CB(ii,l+7*LA,2) = CD6i-CD7i

   80 CONTINUE
   
      DO 81 KP=1,NP
      DO 81 L=1,LA
      DO 81 ii=1,ndim

         I0=L+KP*LA
         J0=L+KP*(LA+JUMP)

      CC0r = CA(ii,I0    ,1)+CA(ii,I0+4*M,1)
      CC0i = CA(ii,I0    ,2)+CA(ii,I0+4*M,2)

      CC1r = CA(ii,I0+  M,1)+CA(ii,I0+5*M,1)
      CC1i = CA(ii,I0+  M,2)+CA(ii,I0+5*M,2)

      CC2r = CA(ii,I0+2*M,1)+CA(ii,I0+6*M,1)
      CC2i = CA(ii,I0+2*M,2)+CA(ii,I0+6*M,2)

      CC3r = CA(ii,I0+3*M,1)+CA(ii,I0+7*M,1)
      CC3i = CA(ii,I0+3*M,2)+CA(ii,I0+7*M,2)

      CC4r = CA(ii,I0    ,1)-CA(ii,I0+4*M,1)
      CC4i = CA(ii,I0    ,2)-CA(ii,I0+4*M,2)

      CC5r = CA(ii,I0+  M,1)-CA(ii,I0+5*M,1)
      CC5i = CA(ii,I0+  M,2)-CA(ii,I0+5*M,2)

      CC6r = -CA(ii,I0+2*M,2)+CA(ii,I0+6*M,2)
      CC6i =  CA(ii,I0+2*M,1)-CA(ii,I0+6*M,1)

      CC7r = -CA(ii,I0+3*M,2)+CA(ii,I0+7*M,2)
      CC7i =  CA(ii,I0+3*M,1)-CA(ii,I0+7*M,1)

         CD0r = CC0r+CC2r
         CD0i = CC0i+CC2i

         CD1r = CC1r+CC3r
         CD1i = CC1i+CC3i

         CD2r = CC4r+CC6r
         CD2i = CC4i+CC6i

            CE1r = CC5r+CC7r
            CE1i = CC5i+CC7i

         CD3r = (CE1r-CE1i)*x
         CD3i = (CE1r+CE1i)*x

         CD4r = CC0r-CC2r
         CD4i = CC0i-CC2i

         CD5r =-(CC1i-CC3i)
         CD5i = (CC1r-CC3r)

         CD6r = CC4r-CC6r
         CD6i = CC4i-CC6i

            CE2r = CC5r-CC7r
            CE2i = CC5i-CC7i

         CD7r = (-CE2r-CE2i)*x
         CD7i = ( CE2r-CE2i)*x

      CF1r = CD2r+CD3r
      CF1i = CD2i+CD3i

      CF2r = CD4r+CD5r
      CF2i = CD4i+CD5i

      CF3r = CD6r+CD7r
      CF3i = CD6i+CD7i

      CF4r = CD0r-CD1r
      CF4i = CD0i-CD1i

      CF5r = CD2r-CD3r
      CF5i = CD2i-CD3i

      CF6r = CD4r-CD5r
      CF6i = CD4i-CD5i

      CF7r = CD6r-CD7r
      CF7i = CD6i-CD7i

      CB(ii,J0     ,1) = CD0r+CD1r
      CB(ii,J0     ,2) = CD0i+CD1i

      CB(ii,J0+  LA,1) = CF1r*CW(1,KP*LA  )-CF1i*CW(2,KP*LA  )
      CB(ii,J0+  LA,2) = CF1r*CW(2,KP*LA  )+CF1i*CW(1,KP*LA  )

      CB(ii,J0+2*LA,1) = CF2r*CW(1,KP*LA*2)-CF2i*CW(2,KP*LA*2)
      CB(ii,J0+2*LA,2) = CF2r*CW(2,KP*LA*2)+CF2i*CW(1,KP*LA*2)

      CB(ii,J0+3*LA,1) = CF3r*CW(1,KP*LA*3)-CF3i*CW(2,KP*LA*3)
      CB(ii,J0+3*LA,2) = CF3r*CW(2,KP*LA*3)+CF3i*CW(1,KP*LA*3)

      CB(ii,J0+4*LA,1) = CF4r*CW(1,KP*LA*4)-CF4i*CW(2,KP*LA*4)
      CB(ii,J0+4*LA,2) = CF4r*CW(2,KP*LA*4)+CF4i*CW(1,KP*LA*4)

      CB(ii,J0+5*LA,1) = CF5r*CW(1,KP*LA*5)-CF5i*CW(2,KP*LA*5)
      CB(ii,J0+5*LA,2) = CF5r*CW(2,KP*LA*5)+CF5i*CW(1,KP*LA*5)

      CB(ii,J0+6*LA,1) = CF6r*CW(1,KP*LA*6)-CF6i*CW(2,KP*LA*6)
      CB(ii,J0+6*LA,2) = CF6r*CW(2,KP*LA*6)+CF6i*CW(1,KP*LA*6)

      CB(ii,J0+7*LA,1) = CF7r*CW(1,KP*LA*7)-CF7i*CW(2,KP*LA*7)
      CB(ii,J0+7*LA,2) = CF7r*CW(2,KP*LA*7)+CF7i*CW(1,KP*LA*7)

   81 CONTINUE
C**
      RETURN
      END
      subroutine setcw(cw,n)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)
      complex*16 cw(0:n-1)

      pai=3.14159265358979323846d0
      pc =2.0d0*pai/dble(n)

      do 10 i=0,n-1
         tr=dcos(pc*dble(i))
         ti=dsin(pc*dble(i))
         cw(i)=dcmplx(tr,ti)
   10 continue

      return
      end
      subroutine setlp(n,lp2,lp3,lp4,lp5,lp8)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $
      implicit real*8(a-b,d-h,o-z)

         nn =n

         lp8=0
         lp5=0
         lp4=0
         lp3=0
         lp2=0

      do while(mod(nn,8).eq.0)
         lp8=lp8+1
         nn=nn/8
      end do

      do while(mod(nn,5).eq.0)
         lp5=lp5+1
         nn=nn/5
      end do

      do while(mod(nn,4).eq.0)
         lp4=lp4+1
         nn=nn/4
      end do

      do while(mod(nn,3).eq.0)
         lp3=lp3+1
         nn=nn/3
      end do

      do while(mod(nn,2).eq.0)
         lp2=lp2+1
         nn=nn/2
      end do

c     if(lp8+lp4+lp2.eq.0) stop 'n is not an even number.'
      if(nn.ne.1) stop 'n is an improper number.'

      return
      end
      subroutine trans_in_b0(ca,cb,idh,n2,n3,mm1,mm2,kc1,kc2,kc3)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $

      implicit real*8(a-h,o-z)

      real*8 ca(2,0:idh-1,0:n2-1 ,0:n3-1)
      real*8 cb(0:mm1-1,0:mm2-1,0:n3-1,2)

c     complex*16 ca(0:idh-1,0:n2-1 ,0:n3-1)
c     complex*16 cb(0:mm1-1,0:mm2-1,0:n3-1)
c*
      if(kc1.eq.0) then
         kk1 = mm1-1
      else
         kk1 = kc1
      end if 

      if(kc2.eq.0) then
         kk2 = n2-1
      else
         kk2 = kc2
      end if 

c*
c* Packing from ca to cb
  
      do 10 i3=0,n3-1
      do 10 i2=0,kk2
      do 10 i1=0,kk1

         cb(i1,i2,i3,1) = ca(1,i1,i2,i3) 
         cb(i1,i2,i3,2) = ca(2,i1,i2,i3) 

  10  continue

      if(kc2.ne.0) then

         nskip2 = kc2*2+1-n2

         do 20 i3=0,n3-1
         do 20 i2=n2-kc2,n2-1
         do 20 i1=0,kk1

            cb(i1,i2+nskip2,i3,1) = ca(1,i1,i2,i3) 
            cb(i1,i2+nskip2,i3,2) = ca(2,i1,i2,i3) 

  20     continue

      end if
c*
      if(kc3.ne.0) then

         do 30 i3=kc3+1,n3-kc3-1
         do 30 i2=0,mm2-1
         do 30 i1=0,mm1-1

            cb(i1,i2,i3,1) = 0.0d0
            cb(i1,i2,i3,2) = 0.0d0

  30     continue

      end if
       
c*
      return
      end
      subroutine trans_in_b1(ca,cb,mm1,n2,mm2,n3,kc2)

      implicit real*8(a-h,o-z)
      real*8 ca(0:mm1-1,0:mm2-1,0:n3-1,2)
      real*8 cb(0:mm1-1,0:n3-1 ,0:n2-1,2)

c     complex*16 ca(0:mm1-1,0:mm2-1,0:n3-1)
c     complex*16 cb(0:mm1-1,0:n3-1 ,0:n2-1)

c*
      if(kc2.eq.0) then
         kk2 = n2-1
      else
         kk2 = kc2
      end if

c*
      do 10 i3=0,n3-1
      do 10 i2=0,kk2
      do 10 i1=0,mm1-1

         cb(i1,i3,i2,1) = ca(i1,i2,i3,1)
         cb(i1,i3,i2,2) = ca(i1,i2,i3,2)

  10  continue

c*
      if(kc2.ne.0) then

         do 20 i2=kc2+1,n2-kc2-1
         do 20 i3=0,n3-1
         do 20 i1=0,mm1-1

            cb(i1,i3,i2,1) = 0.0d0 
            cb(i1,i3,i2,2) = 0.0d0 

  20     continue

         nstep2 = kc2*2-n2+1

         do 30 i3=0,n3-1
         do 30 i2=n2-kc2,n2-1
         do 30 i1=0,mm1-1

            cb(i1,i3,i2,1) = ca(i1,i2+nstep2,i3,1)
            cb(i1,i3,i2,2) = ca(i1,i2+nstep2,i3,2)

  30     continue

      end if
c*

      return
      end
      subroutine trans_in_b2(ca,cb,mm1,n2,n3,n1h)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $

      implicit real*8(a-h,o-z)

      real*8 ca(0:mm1-1  ,0:n3*n2-1 ,2)
      real*8 cb(0:n3*n2-1,0:n1h     ,2)

c     real*8 cb(0:n3*n2-1,0:mm1-1   ,2)
c     complex*16 ca(0:mm1-1  ,0:n3*n2-1 )
c     complex*16 cb(0:n3*n2-1,0:mm1-1   )

         do 10 i2 = 0,n3*n2-1
         do 10 i1 = 0,mm1-1

            cb(i2,i1,1) = ca(i1,i2,1)
            cb(i2,i1,2) = ca(i1,i2,2)

 10      continue

c*
         do 20 i1 = mm1,n1h
         do 20 i2 = 0,n3*n2-1

            cb(i2,i1,1) = 0.0d0
            cb(i2,i1,2) = 0.0d0

 20      continue

c*

      return
      end 
      subroutine trans_in_b3(cb,ca,cw1,nn,n1h,kc1)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $

      implicit real*8(a-h,o-z)

      real*8 cb(nn,0:n1h  ,2)
      real*8 ca(nn,0:n1h-1,2)

c     complex*16 cb(nn,0:idh-1)
c     complex*16 ca(nn,0:idh-1)

      real*8 cw1(2,0:n1h-1)
c     complex*16 cw1(0:n1h-1)

      real*8 r0,r1
      real*8 rr0,rr1,ri0,ri1

c*
      do 10 ii=1,nn

         r0 = cb(ii,  0,1)
         r1 = cb(ii,n1h,1)

c        r0 = dble(cb(ii,  0))
c        r1 = dble(cb(ii,n1h))

         ca(ii,0,1)= r0 + r1 
         ca(ii,0,2)= r0 - r1 

c        ca(ii,0)=dcmplx( r0 + r1 , r0 - r1 )

  10  continue

      do 20 k1=1,n1h-1
      do 20 ii=1,nn

         rr0 = cb(ii,k1,1)
         ri0 = cb(ii,k1,2)

         rr1 = cb(ii,n1h-k1,1)
         ri1 = cb(ii,n1h-k1,2)

c        rr2 = rr0 + rr1 
c        ri2 = ri0 - ri1 

         rr3 = rr0 - rr1 
         ri3 = ri0 + ri1 

         ca(ii,k1,1) = rr0 + rr1 + rr3*cw1(1,k1) - ri3*cw1(2,k1)
         ca(ii,k1,2) = ri0 - ri1 + rr3*cw1(2,k1) + ri3*cw1(1,k1)

c        ca(ii,k1)=   cb(ii,k1)+dconjg(cb(ii,n1h-k1))
c    &             + (cb(ii,k1)-dconjg(cb(ii,n1h-k1)))*cw1(k1) 

  20  continue 

      return
      end 
      subroutine trans_in_b4(cb,ca,n3,n2,idh,n1h) 
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $

c**
      implicit real*8(a-h,o-z)

      real*8 cb(0:n3-1 ,0:n2-1,0:n1h-1,2)
      real*8 ca(2,0:idh-1,0:n2-1,0:n3-1 )

c     complex*16 cb(0:n3-1 ,0:n2-1,0:idh-1)
c     complex*16 ca(0:idh-1,0:n2-1,0:n3-1 )
c**
         do 10 i1=0,n1h-1
         do 10 i2=0,n2-1
         do 10 i3=0,n3-1

            ca(1,i1,i2,i3) = cb(i3,i2,i1,1)
            ca(2,i1,i2,i3) = cb(i3,i2,i1,2)

  10     continue

         do 20 i3=0,n3-1
         do 20 i2=0,n2-1
         do 20 i1=n1h,idh-1

            ca(1,i1,i2,i3) = 0.0d0
            ca(2,i1,i2,i3) = 0.0d0

  20     continue

c**
      return
      end
      subroutine trans_in_f1(ca,cb,idh,n1h,n2,n3,
     &                       mm2,mm3,kc1,kc2,kc3)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $

      implicit real*8(a-h,o-z)

      real*8 ca(2,0:idh-1,0:n2-1,0:n3-1)
      real*8 cb(0:mm3-1,0:mm2-1,0:n1h-1,2)

c***
      if(kc1.eq.0) then
         kk1 = n1h-1
      else
         kk1 = kc1
      end if

      if(kc2.eq.0) then
         kk2 = n2-1
      else
         kk2 = kc2
      end if

      if(kc3.eq.0) then
         kk3 = n3-1
      else
         kk3 = kc3
      end if
c**

         do 10 i3=0,kk3
         do 10 i2=0,kk2
         do 10 i1=0,kk1

            cb(i3,i2,i1,1) =   ca(1,i1,i2,i3) 
            cb(i3,i2,i1,2) = - ca(2,i1,i2,i3) 

   10    continue

c**
      if(kc2.ne.0) then

         do 20 i3=0,kk3
         do 20 i2=n2-kc2,n2-1
         do 20 i1=0,kk1

            ii2 = i2 + mm2 - n2
            cb(i3,ii2,i1,1) =   ca(1,i1,i2,i3) 
            cb(i3,ii2,i1,2) = - ca(2,i1,i2,i3) 

   20    continue
  
      end if
c**
      if(kc3.ne.0) then

         do 30 i3=n3-kc3,n3-1
         do 30 i2=0,kk2
         do 30 i1=0,kk1

            ii3 = i3 + mm3 - n3
            cb(ii3,i2,i1,1) =   ca(1,i1,i2,i3) 
            cb(ii3,i2,i1,2) = - ca(2,i1,i2,i3) 

   30    continue

      end if
c**
      if(kc2.ne.0.and.kc3.ne.0) then

         do 40 i3=n3-kc3,n3-1
         do 40 i2=n2-kc2,n2-1
         do 40 i1=0,kk1

            ii2 = i2 + mm2 - n2
            ii3 = i3 + mm3 - n3
            cb(ii3,ii2,i1,1) =   ca(1,i1,i2,i3) 
            cb(ii3,ii2,i1,2) = - ca(2,i1,i2,i3) 

   40    continue

      end if
c**
         do 50 i1=kk1+1,n1h-1
         do 50 i2=0,mm2-1
         do 50 i3=0,mm3-1

            cb(i3,i2,i1,1) = 0.0d0
            cb(i3,i2,i1,2) = 0.0d0

   50    continue
c**

      return
      end
      subroutine trans_in_f2(ca,cb,cw1,mm3,mm2,n1h)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $

      implicit real*8 (a-h,o-z)

      real*8 ca(mm3*mm2,0:n1h-1,2)
      real*8 cb(mm3*mm2,0:n1h  ,2)
      real*8 cw1(2,0:n1h-1)

c     complex*16 ca(mm3*mm2,0:n1h)
c     complex*16 cb(mm3*mm2,0:n1h)
c     complex*16 cw1(0:n1h-1)

      real*8 rr,ri
c**
      do i32 = 1,mm3*mm2

         rr = ca(i32,0,1)
         ri = ca(i32,0,2)

c        rr = dble (ca(i32,0))
c        ri = dimag(ca(i32,0))
       
         cb(i32,0  ,1) = rr-ri
         cb(i32,n1h,1) = rr+ri

         cb(i32,0  ,2) = 0.0d0 
         cb(i32,n1h,2) = 0.0d0 

c        cb(i32,0  )=dcmplx(rr-ri,0.0d0)
c        cb(i32,n1h)=dcmplx(rr+ri,0.0d0)

      end do
c**   
      do k1  = 1,n1h-1
      do ii = 1,mm3*mm2

         rr0 = ca(ii,k1,1)
         ri0 = ca(ii,k1,2)

         rr1 = ca(ii,n1h-k1,1)
         ri1 = ca(ii,n1h-k1,2)

c        rr2 = rr0 + rr1
c        ri2 = ri0 - ri1

         rr3 = rr0 - rr1
         ri3 = ri0 + ri1

         cb(ii,k1,1) = (rr0 + rr1 + rr3*cw1(1,k1) - ri3*cw1(2,k1))*0.5d0
         cb(ii,k1,2) = (ri0 - ri1 + rr3*cw1(2,k1) + ri3*cw1(1,k1))*0.5d0


c     cb(i32,k1) = ( ca(i32,k1)+dconjg(ca(i32,n1h-k1))
c    &            +( ca(i32,k1)-dconjg(ca(i32,n1h-k1)) )*cw1(k1))
c    &            * 0.5d0

      end do
      end do

      return
      end
      subroutine trans_in_f3(cb,ca,mm3,n2,mm2,n1h,kc2)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $

      implicit real*8(a-h,o-z)

      real*8 ca(0:mm3-1,0:n1h,0:n2-1,2)
      real*8 cb(0:mm3-1,0:mm2-1,0:n1h,2)
c     complex*16 ca(0:mm3-1,0:idh-1,0:n2-1)
c     complex*16 cb(0:mm3-1,0:mm2-1,0:idh-1)

c***********************************************
c**   transpose between 2nd axis and 3rd axis.
c***********************************************

      if(kc2.eq.0) then
         kk2 = n2-1
      else
         kk2 = kc2 
      end if

      do 10 i1=0,n1h
      do 10 i2=0,kk2
      do 10 i3=0,mm3-1

         ca(i3,i1,i2,1) = cb(i3,i2,i1,1)
         ca(i3,i1,i2,2) = cb(i3,i2,i1,2)

  10  continue  

c**
      if(kc2.ne.0) then

         do 30 ii2=kk2+1,n2-kc2-1
         do 30 i1=0,n1h
         do 30 i3=0,mm3-1

            ca(i3,i1,ii2,1) = 0.0d0
            ca(i3,i1,ii2,2) = 0.0d0

  30     continue  

         do 20 i1=0,n1h
         do 20 i2=kk2+1,mm2-1
         do 20 i3=0,mm3-1

            ii2 = i2 + n2 - mm2
            ca(i3,i1,ii2,1) = cb(i3,i2,i1,1)
            ca(i3,i1,ii2,2) = cb(i3,i2,i1,2)

  20     continue  

      end if
c**
      return
      end
      subroutine trans_in_f4(cb,ca,n1h,n2,n3,mm3,kc3)

      implicit real*8(a-h,o-z)

      real*8 ca((n1h+1)*n2 ,0:n3-1,2)
      real*8 cb(0:mm3-1,(n1h+1)*n2,2)
c     complex*16 ca(0:idh-1,0:n2-1 ,0:n3-1)
c     complex*16 cb(0:mm3-1,0:idh-1,0:n2-1)

c*
      if(kc3.eq.0) then
         kk3 = n3 - 1
      else
         kk3 = kc3
      end if

c     do 10 i2 =0,n2-1  
      do 10 i12 =1,(n1h+1)*n2
      do 10 i3 =0,kk3
      
         ca(i12,i3,1) = cb(i3,i12,1)
         ca(i12,i3,2) = cb(i3,i12,2)

  10  continue

c**
      if(kc3.ne.0) then

         do 20 ii3 = kk3+1,n3-kc3-1
c        do 20 i2 = 0,n2-1
         do 20 i12 =1,(n1h+1)*n2
      
            ca(i12,ii3,1) = 0.0d0
            ca(i12,ii3,2) = 0.0d0

  20     continue

c        do 30 i2 = 0,n2-1
         do 30 i12 =1,(n1h+1)*n2
         do 30 i3 = kk3+1,mm3-1
      
            ii3 = i3 + n3 - mm3 
            ca(i12,ii3,1) = cb(i3,i12,1)
            ca(i12,ii3,2) = cb(i3,i12,2)

  30     continue

      end if
c**
      return
      end


      subroutine trans_in_f5(cb,ca,idh,n1h,n2,n3)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $

      implicit real*8(a-h,o-z)

      real*8 ca(2,0:idh-1,0:n2*n3-1)
      real*8 cb(  0:n1h  ,0:n2*n3-1,2)
c     complex*16 ca(0:idh-1,0:n2-1,0:n3-1)
c     complex*16 cb(0:idh-1,0:n2-1,0:n3-1)
c*
      do 10 i23=0,n2*n3-1 
      do 10 i1=0,n1h 

         ca(1,i1,i23) =   cb(i1,i23,1)
         ca(2,i1,i23) = - cb(i1,i23,2)
c        ca(i1,i2,i3) = dconjg(cb(i1,i2,i3))

  10  continue                        

      do 20 i1=n1h+1,idh-1 
      do 20 i23=0,n2*n3-1 

         ca(1,i1,i23) = 0.0d0
         ca(2,i1,i23) = 0.0d0 

  20  continue                        

c*

      return
      end
      subroutine trans_out_b0(ca,cb,idh,n2,n3,n1h)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $

      implicit real*8(a-h,o-z)

      real*8 ca(2,0:idh-1,n2*n3)
      real*8 cb(0:n1h,n2*n3,2)

c     complex*16 ca(0:idh*n2*n3-1)
c     complex*16 cb(0:idh*n2*n3-1)
c*
  
      do 10 ii=1,n2*n3
      do 10  i=0,n1h   

         cb(i,ii,1) = ca(1,i,ii) 
         cb(i,ii,2) = ca(2,i,ii) 

  10  continue

c*
      return
      end
      subroutine trans_out_b1(ca,cb,n1h,n2,n3,mm3,kc3)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $

      real*8 ca(0:n1h,0:n2-1  ,0:n3-1,2)
      real*8 cb(0:n1h,0:mm3-1 ,0:n2-1,2)

c     complex*16 ca(0:idh-1,0:n2-1  ,0:n3-1)
c     complex*16 cb(0:idh-1,0:mm3-1 ,0:n2-1)

c*
      if(kc3.eq.0) then
         kk3 = n3-1
      else
         kk3 = kc3
      end if

c*
      do 10 i3=0,kk3
      do 10 i2=0,n2-1
      do 10 i1=0,n1h  

         cb(i1,i3,i2,1) = ca(i1,i2,i3,1)
         cb(i1,i3,i2,2) = ca(i1,i2,i3,2)

  10  continue

c*
      if(kc3.ne.0) then

         do 30 i3=n3-kc3,n3-1
         do 30 i2=0,n2-1
         do 30 i1=0,n1h   

            ii3 = i3 + mm3 - n3
            cb(i1,ii3,i2,1) = ca(i1,i2,i3,1)
            cb(i1,ii3,i2,2) = ca(i1,i2,i3,2)

  30     continue

      end if
c*

      return
      end
      subroutine trans_out_b2(ca,cb,n1h,n2,mm2,mm3,kc2)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $

      implicit real*8(a-h,o-z)
      real*8 ca(0:n1h,0:mm3-1,0:n2-1 ,2)
      real*8 cb(0:mm3-1,0:mm2-1,0:n1h,2)

c     complex*16 ca(0:idh-1,0:mm3-1,0:n2-1 )
c     complex*16 cb(0:mm3-1,0:mm2-1,0:idh-1)

      if(kc2.eq.0) then
         kk2 = n2-1
      else
         kk2 = kc2 
      end if

         do 10 i2 = 0,kk2
         do 10 i3 = 0,mm3-1
         do 10 i1 = 0,n1h   

            cb(i3,i2,i1,1) = ca(i1,i3,i2,1)
            cb(i3,i2,i1,2) = ca(i1,i3,i2,2)

 10      continue

      if(kc2.ne.0) then

         do 20 i2 = n2-kc2,n2-1
         do 20 i3 = 0,mm3-1
         do 20 i1 = 0,n1h   

            ii2 = i2 + mm2 - n2
            cb(i3,ii2,i1,1) = ca(i1,i3,i2,1)
            cb(i3,ii2,i1,2) = ca(i1,i3,i2,2)

 20      continue

      end if

      return
      end
      subroutine trans_out_b3(cb,ca,cw1,nn,n1h)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $

      implicit real*8(a-h,o-z)

      real*8 cb(nn,0:n1h  ,2)
      real*8 ca(nn,0:n1h-1,2)
      real*8 cw1(2,0:n1h-1)

c     complex*16 cb(nn,0:n1h)
c     complex*16 ca(nn,0:n1h-1)
c     complex*16 cw1(0:n1h-1)

      real*8 r0,r1

c**
      do 10 ii=1,nn

         r0 = cb(ii,  0,1)
         r1 = cb(ii,n1h,1)

c        r0 = dble(cb(ii,  0))
c        r1 = dble(cb(ii,n1h))

         ca(ii,0,1)= r0 + r1 
         ca(ii,0,2)= r0 - r1 

c        ca(ii,0)=dcmplx( r0 + r1 , r0 - r1 )

  10  continue

      do 20 k1=1,n1h-1
      do 20 ii=1,nn

         rr0 = cb(ii,k1,1)
         ri0 = cb(ii,k1,2)

         rr1 = cb(ii,n1h-k1,1)
         ri1 = cb(ii,n1h-k1,2)

c        rr2 = rr0 + rr1
c        ri2 = ri0 - ri1

         rr3 = rr0 - rr1
         ri3 = ri0 + ri1

         ca(ii,k1,1) = rr0 + rr1 + rr3*cw1(1,k1) - ri3*cw1(2,k1)
         ca(ii,k1,2) = ri0 - ri1 + rr3*cw1(2,k1) + ri3*cw1(1,k1)

c        ca(ii,k1)=   cb(ii,k1)+dconjg(cb(ii,n1h-k1))
c    &             + (cb(ii,k1)-dconjg(cb(ii,n1h-k1)))*cw1(k1) 

  20  continue 

      return
      end 
      subroutine trans_out_b4(cb,ca,idh,n1h,n2,n3,mm2,mm3,
     &                        kc1,kc2,kc3)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $
c**
      implicit real*8(a-h,o-z)

      real*8 cb(0:mm3-1 ,0:mm2-1,0:n1h-1,2)
      real*8 ca(2,0:idh-1 ,0:n2-1 ,0:n3-1 )

c     complex*16 cb(0:mm3-1 ,0:mm2-1,0:idh-1)
c     complex*16 ca(0:idh-1 ,0:n2-1 ,0:n3-1 )
c**
      if(kc1.eq.0) then
         kk1 = n1h-1
      else
         kk1 = kc1
      end if

      if(kc2.eq.0) then
         kk2 = n2-1
      else
         kk2 = kc2
      end if

      if(kc3.eq.0) then
         kk3 = n3-1
      else
         kk3 = kc3
      end if
c**
      do 100 i3 = 0,n3-1
      do 100 i2 = 0,n2-1
      do 100 i1 = 0,idh-1

         ca(1,i1,i2,i3) = 0.0d0
         ca(2,i1,i2,i3) = 0.0d0

 100  continue

c**
      do 10 i1=0,kk1
      do 10 i2=0,kk2
      do 10 i3=0,kk3

         ca(1,i1,i2,i3) = cb(i3,i2,i1,1)
         ca(2,i1,i2,i3) = cb(i3,i2,i1,2)

  10  continue
c**
      if(kc2.ne.0) then

         do 20 i1=0,kk1
         do 20 i2=kc2+1,kc2+kc2
         do 20 i3=0,kk3

            ii2 = i2 + n2-1-kc2*2
            ca(1,i1,ii2,i3) = cb(i3,i2,i1,1)
            ca(2,i1,ii2,i3) = cb(i3,i2,i1,2)

  20     continue

      end if
c**
      if(kc3.ne.0) then

         do 30 i1=0,kk1
         do 30 i2=0,kk2
         do 30 i3=kc3+1,kc3+kc3

            ii3 = i3 + n3-1-kc3*2
            ca(1,i1,i2,ii3) = cb(i3,i2,i1,1)
            ca(2,i1,i2,ii3) = cb(i3,i2,i1,2)

  30     continue

      end if
c**
      if(kc2.ne.0.and.kc3.ne.0) then

         do 40 i1=0,kk1
         do 40 i2=kc2+1,kc2+kc2
         do 40 i3=kc3+1,kc3+kc3

            ii2 = i2 + n2-1-kc2*2
            ii3 = i3 + n3-1-kc3*2
            ca(1,i1,ii2,ii3) = cb(i3,i2,i1,1)
            ca(2,i1,ii2,ii3) = cb(i3,i2,i1,2)

  40     continue

      end if

c**
      return
      end  
      subroutine trans_out_f1(ca,cb,idh,n1h,n2,n3)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $

      implicit real*8(a-h,o-z)

      real*8 ca(2,0:idh-1,0:n2-1,0:n3-1)
      real*8 cb(0:n3-1 ,0:n2-1,0:n1h-1,2)

c     complex*16 ca(0:m1-1,0:m2-1,0:m3-1),cb(0:m3-1,0:m2-1,0:m1-1)

c***
         do 10 i3=0,n3-1
         do 10 i2=0,n2-1
         do 10 i1=0,n1h-1

            cb(i3,i2,i1,1) =   ca(1,i1,i2,i3) 
            cb(i3,i2,i1,2) = - ca(2,i1,i2,i3) 

   10    continue

      return
      end
      subroutine trans_out_f2(ca,cb,cw1,n3,n2,n1h,kc1)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $

      implicit real*8(a-h,o-z)

      real*8 ca(n3*n2,0:n1h-1,2)
      real*8 cb(n3*n2,0:n1h  ,2)
      real*8 cw1(2,0:n1h-1)
      real*8 rr,ri
c**
      do i32 = 1,n3*n2

         rr = ca(i32,0,1)
         ri = ca(i32,0,2)
c        rr = dble (ca(i32,0))
c        ri = dimag(ca(i32,0))

         cb(i32,0  ,1) = rr-ri
         cb(i32,n1h,1) = rr+ri

         cb(i32,0  ,2) = 0.0d0
         cb(i32,n1h,2) = 0.0d0

c        rr = dble (ca(i32,0))
c        ri = dimag(ca(i32,0))
       
c        cb(i32,0  )=dcmplx(rr-ri,0.0d0)
c        cb(i32,n1h)=dcmplx(rr+ri,0.0d0)

      end do
c**   
      if(kc1.eq.0) then
         mm = n1h-1
      else
         mm = kc1
      end if
c**
      do k1  = 1,mm
      do ii = 1,n3*n2

         rr0 = ca(ii,k1,1)
         ri0 = ca(ii,k1,2)

         rr1 = ca(ii,n1h-k1,1)
         ri1 = ca(ii,n1h-k1,2)

         rr2 = rr0 + rr1
         ri2 = ri0 - ri1

         rr3 = rr0 - rr1
         ri3 = ri0 + ri1

         cb(ii,k1,1) = ( rr2 + rr3*cw1(1,k1) - ri3*cw1(2,k1) )*0.5d0
         cb(ii,k1,2) = ( ri2 + rr3*cw1(2,k1) + ri3*cw1(1,k1) )*0.5d0

c     cb(i32,k1) = ( ca(i32,k1)+dconjg(ca(i32,n1h-k1))
c    &            +( ca(i32,k1)-dconjg(ca(i32,n1h-k1)) )*cw1(k1))
c    &            * 0.5d0

      end do
      end do

      return
      end
      subroutine trans_out_f3(cb,ca,n3,n2,mm1,n1h)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $

      implicit real*8(a-h,o-z)
      real*8 ca(n3,0:mm1-1,n2,2)
      real*8 cb(n3,n2,0:n1h,2)

c***********************************************
c**   transpose between 2nd axis and 3rd axis.
c***********************************************

      do 10 i1=0,mm1-1
      do 10 i2=1,n2
      do 10 i3=1,n3

         ca(i3,i1,i2,1) = cb(i3,i2,i1,1)
         ca(i3,i1,i2,2) = cb(i3,i2,i1,2)

  10  continue  

      return
      end
      subroutine trans_out_f4(cb,ca,n3,mm1,n2,mm2,kc2)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $

      implicit real*8(a-h,o-z)
      real*8 ca(0:mm1-1,0:mm2-1, 0:n3-1,2)
      real*8 cb(0:n3-1 ,0:mm1-1, 0:n2-1,2)

c*
      if(kc2.eq.0) then
         kk2 = n2 - 1
      else
         kk2 = kc2
      end if

      do 10 i2 =0,kk2  
      do 10 i1 =0,mm1-1
      do 10 i3 =0,n3-1
      
         ca(i1,i2,i3,1) = cb(i3,i1,i2,1)
         ca(i1,i2,i3,2) = cb(i3,i1,i2,2)

  10  continue

c*
      if(kc2.ne.0) then

         do 20 i2 = -kc2,-1
         do 20 i1 =    0,mm1-1
         do 20 i3 =    0,n3-1
      
            ii2 = i2 + kc2*2 + 1 
            ca(i1,ii2,i3,1) = cb(i3,i1,i2+n2,1)
            ca(i1,ii2,i3,2) = cb(i3,i1,i2+n2,2)

  20     continue

      end if

c*
      return
      end
      subroutine trans_out_f5(cb,ca,idh,n1h,n2,n3,mm1,mm2,kc1,kc2,kc3)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $

      implicit real*8(a-h,o-z)
      real*8 ca(2,0:idh-1,0:n2-1 ,0:n3-1 )
      real*8 cb(0:mm1-1,0:mm2-1,0:n3-1 ,2)
c*
      if(kc1.eq.0) then
           kk1 = n1h
      else
           kk1 = kc1
      end if

      if(kc2.eq.0) then
           kk2 = n2-1
      else
           kk2 = kc2
      end if

      if(kc3.eq.0) then
           kk3 = n3-1
      else
           kk3 = kc3
      end if
c*

      do 100 i3=0,n3-1 
      do 100 i2=0,n2-1 
      do 100 i1=0,idh-1 

         ca(1,i1,i2,i3) = 0.0d0
 100     ca(2,i1,i2,i3) = 0.0d0

c*

      do 10 i3=0,kk3 
      do 10 i2=0,kk2 
      do 10 i1=0,kk1 

         ca(1,i1,i2,i3) =    cb(i1,i2,i3,1)
         ca(2,i1,i2,i3) =  - cb(i1,i2,i3,2)

  10  continue 

      if(kc2.ne.0) then

         nskip2 = n2-1-kc2*2

         do 11 i3=0,kk3 
         do 11 i2=kc2+1,kc2*2 
         do 11 i1=0,kk1 

            ca(1,i1,i2+nskip2,i3) =   cb(i1,i2,i3,1)
            ca(2,i1,i2+nskip2,i3) = - cb(i1,i2,i3,2)

  11     continue 

      end if

c******************************
      if(kc3.ne.0) then
c******************************

      do 20 i3=n3-kc3,n3-1 
      do 20 i2=0,kk2 
      do 20 i1=0,kk1 

         ca(1,i1,i2,i3) =   cb(i1,i2,i3,1)
         ca(2,i1,i2,i3) = - cb(i1,i2,i3,2)

  20  continue 

      if(kc2.ne.0) then

         nskip2 = n2-1-kc2*2

         do 21 i3=n3-kc3,n3-1 
         do 21 i2=kc2+1,kc2*2 
         do 21 i1=0,kk1 

            ca(1,i1,i2+nskip2,i3) =   cb(i1,i2,i3,1)
            ca(2,i1,i2+nskip2,i3) = - cb(i1,i2,i3,2)

  21     continue 

      end if

c******************************
      end if
c******************************

      return
      end
      subroutine jrcat_c3ft(ca,cwork,id,n1,n2,n3,cw1,cw2,cw3,
     &                      kc1,kc2,kc3,key)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $

c**********************************************************************
c*    3 dimensional complex to complex FFT with mask
c*    either for input array or output array,
c*    developed by T.Sanada (sanada@think.com) on 1996.01.04.
c*    It is specified for Fujitsu VPP500.
c**********************************************************************

      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)

c------------------>      for Estimation of CPU cost.
c                             by T.Yamasaki
c                               22th May. 1996
      real*8 CRTVL, PCPUDF
      parameter (CRTVL = 1.d-5, PCPUDF = 2.0)
c        CRTVL : CRItical VaLue for Division.
c        PCPUDF: Percent of CPU DiFference.
      real*8 cpu0ol,cpudif,rcpudf,ecpu,tcpu
      integer icount
      data icount/0/
      data ecpu/0.d0/
      data cpu0ol/0.d0/
      real*8 t_start, t_end, UMICRO
      parameter (UMICRO = 1.d-6)
c <-------------------------------------------

      complex*16 ca    (id,n2,n3)
      complex*16 cwork (id,n2,n3)

      complex*16 cw1(0:n1-1)
      complex*16 cw2(0:n2-1)
      complex*16 cw3(0:n3-1)

      integer kc1,kc2,kc3
      integer nbank

c**********************************************************************
c*   key = 0 : set up of FFT.
c*             ( Calculation of cw1,cw2 and cw3. )
c**********************************************************************
c*   n1,n2 and n3 must include radix 2, and they may include 3 or 5.
c*   id must be n1 + 2.
c*
c**********************************************************************
c*  key = -1 & -2 : F.T. of complex 3d array ca(0:n1-1,0:n2-1,0:n3-1) .
c**********************************************************************
c*
c*   The result(complex 'cb') is overwritten on 'ca'.
c*   
c*     cb(k1,k2,k3) 
c*   = sum(ix=0:n1-1,iy=0:n2-1,iz=0:n3-1) 
c*     ca(i1,i2,i3) exp [- 2 * pai * ci
c*                             * ( k1*i1/n1 + k2*i2/n2 + k3*i3/n3 ) ]
c*     
c*     where 0 <= k1 <= n1/2, 0 <= k2 <= n2-1, 0 <= k3 <= n3-1.
c*
c**********************************************************************
c*  key = +1 & +2: Inverse F.T. of copmplex 3d array cb(0:n1/2,0:n2-1
c*                                                           ,0:n3-1) .
c**********************************************************************
c*   Transform from 'cb' to 'ra'.
c*
c**********************************************************************
c*  kc1,kc2 kc3 are integer parameters specifying the mask for
c*  input data or output data.
c*
c*     If kc1.ge.n1/2.or.kc1.lt.0, then kc1 is forced to be 0. 
c*     If kc2.ge.n2/2.or.kc2.lt.0, then kc2 is forced to be 0. 
c*     If kc3.ge.n3/2.or.kc3.lt.0, then kc3 is forced to be 0. 
c*
c**********************************************************************
c** REMARK on Mask.
c**********************************************************************
c*
c*  key = -1 : Mask is for Output data ( cb ),
c*  key = +1 : Mask is for Input  data ( cb ),
c*  key = -2 : Mask is for Input  data ( ca ),
c*  key = +2 : Mask is for Output data ( ca ),
c*
c*     where ca(0:idh/2,0:n2-1,0:n3-1) and equivalence(ra,ca).
c*
c**********************************************************************
c*     Mask for Output data:
c**********************************************************************
c*
c*     If kc1.eq.0.and.kc2.eq.0.and.kc3.eq.0, then the full wavenumber
c*     range is calculated.
c*
c*     If kc1.ne.0, then only "0 <= k1 <= kc1" range is calculated.
c*
c*     If kc2.ne.0, then only "0 <= k2 <= kc2" and "n2-kc2 <= k2<= n2-1" 
c*     range is calculated.
c*
c*     If kc3.ne.0, then only "0 <= k3 <= kc3" and "n3-kc3 <= k3<= n3-1" 
c*     range is calculated.
c*
c**********************************************************************
c*     Mask for Input data:
c**********************************************************************
c*
c*     If kc1.eq.0.and.kc2.eq.0.and.kc3.eq.0, then the full wavenumber
c*     range is used to obtain real array 'ra'.
c*
c*     If kc1.ne.0, then only "0 <= k1 <= kc1" range is used.
c*
c*     If kc2.ne.0, then only "0 <= k2 <= kc2" and "n2-kc2 <= k2<= n2-1" 
c*     range is used.
c*
c*     If kc3.ne.0, then only "0 <= k3 <= kc3" and "n3-kc3 <= k3<= n3-1" 
c*     range is used.
c*
c**********************************************************************
c*
c******************************************************************
c**   choose nbank = 0 or nbank = 1 
c******************************************************************
c*   nbank = 0:
c*      The transpose is normally implemented.
c*      Bank conflict may appear patriculary in Vector Computers.
c*   nbank = 1:
c*      The transpose is implemented so as to prevent
c*      the bank conflict in the case the conflict may appear.
c******************************************************************
c
      !call gettod(t_start)

      nbank = 0

c**

      if(key.eq.0) then

c** set up of FFT.

         call c3ft_0(cw1,cw2,cw3,n1,n2,n3)

         goto 1001
c$$$         return

      end if

c**
      if(kc1.ge.n1/2.or.kc1.lt.0) then 
         write(6,*) 'warning in FFT: kc1 is irrelevant.'
         write(6,*) 'kc1 has changed to be 0.'
         kc1 = 0
      end if

      if(kc2.ge.n2/2.or.kc2.lt.0) then 
         write(6,*) 'warning in FFT: kc2 is irrelevant.'
         write(6,*) 'kc2 has changed to be 0.'
         kc2 = 0
      end if

      if(kc3.ge.n3/2.or.kc3.lt.0) then 
         write(6,*) 'warning in FFT: kc3 is irrelevant.'
         write(6,*) 'kc3 has changed to be 0.'
         kc3 = 0
      end if

c******************************************
      if(key.eq.-1.or.key.eq.+2) then
c******************************************

c** Forword  FFT with mask for output.
c** Backword FFT with mask for output.

         call c3ft_mask_for_out(ca,cwork,id,n1,n2,n3,cw1,cw2,cw3,
     &                     kc1,kc2,kc3,nbank,key)

c******************************************
      else 
c******************************************

c** Forword  FFT with mask for input.
c** Backword FFT with mask for input.

         call c3ft_mask_for_in (ca,cwork,id,n1,n2,n3,cw1,cw2,cw3,
     &                     kc1,kc2,kc3,nbank,key)


c******************************************
      end if
c******************************************

 1001 continue

c -------->  CPU cost estimation
      !call gettod(t_end)
      tcpu = (t_end - t_start)*UMICRO
      icount = icount + 1
      ecpu = ecpu + tcpu
c <-------

      return
      entry jrc_c3strt
      icount = 0
      ecpu = 0.d0
      return
c--------------------------
      entry jrc_c3end
      cpudif = dabs(cpu0ol - ecpu)
      if(ecpu.gt.CRTVL) then
         rcpudf = cpudif/ecpu * 100.0
      else if(cpu0ol.le.CRTVL) then
         rcpudf = 0.0
      else
         rcpudf = 100.0
      endif
c               <-- Ratio in percent of CPU time difference
c                 between previous one and present one.
      cpu0ol = ecpu

c$$$      call eqivvl(rcpudf)

      if(rcpudf.gt.PCPUDF) then
         if(icount.gt.0) then
            tcpu = ecpu/dfloat(icount)
         endif
         !!write(6,9001) ecpu, icount, tcpu
 9001    format(1H ,' <<CPU TIME for JRCAT_C3FT = ',f11.3
     &        ,'(Sec.) = (',i6,' *',f10.5,')>>')
      endif
      return
      end  

      subroutine c3ft_0(cw1,cw2,cw3,n1,n2,n3)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $
c**
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)

      complex*16 cw1(0:n1-1)
      complex*16 cw2(0:n2-1)
      complex*16 cw3(0:n3-1)

c**
         call setcw(cw1,n1)
         call setcw(cw2,n2)
         call setcw(cw3,n3) 

      return
      end  

      subroutine c3ft_mask_for_in(ca,cb,id,n1,n2,n3,cw1,cw2,cw3,
     &                            kc1,kc2,kc3,nbank,key)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $
c**
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)

      real*8 ca (0:id-1,0:n2-1,0:n3-1,2)
      real*8 cb (0:id-1,0:n2-1,0:n3-1,2)

c     complex*16 ca (0:id-1,0:n2-1,0:n3-1)
c     complex*16 cb (0:id-1,0:n2-1,0:n3-1)

      complex*16 cw1(0:n1-1)
      complex*16 cw2(0:n2-1)
      complex*16 cw3(0:n3-1)

      integer kc1,kc2,kc3
      integer nbank
c**
      if(kc1.eq.0) then
           mm1 = n1
      else
           mm1 = kc1*2+1
      end if

      if(kc2.eq.0) then
           mm2 = n2
      else
           mm2 = kc2*2+1
      end if

      if(kc3.eq.0) then
           mm3 = n3
      else
           mm3 = kc3*2+1
      end if

c** Packin from ca(2,0:id-1,0:n2-1,0:n3-1) to cb(0:mm1-1,0:mm2-1,0:n3-1,2)

      call trans_in_c0(ca,cb,id,n1,n2,n3,mm1,mm2,kc1,kc2,kc3,key)

c** Length:n3 CFFT  

      call setlp(n3,lp2,lp3,lp4,lp5,lp8)

      call fsb_lp235(cb,ca,cw3,n3,mm1*mm2,lp2,lp3,lp4,lp5,lp8)

c** cb(i1,i3,i2,2) <- ca(i1,i2,i3,2)

      call trans_in_c1(ca,cb,mm1,n2,mm2,n3,kc2)

c** Length:n2 CFFT    

      call setlp(n2,lp2,lp3,lp4,lp5,lp8)

      call fsb_lp235(cb,ca,cw2,n2,mm1*n3,lp2,lp3,lp4,lp5,lp8) 

c** cb(i3,i2,i1,2) <- ca(i1,i3,i2,2)

      call trans_in_c2(ca,cb,n1,mm1,n2,n3,kc1)

c** Length:n1 CFFT  

      call setlp(n1,lp2,lp3,lp4,lp5,lp8)

      call fsb_lp235(cb,ca,cw1,n1,n3*n2,lp2,lp3,lp4,lp5,lp8)

c** cb(2,i1,i2,i3) <- ca(i3,i2,i1,2)

      call trans_in_c3(ca,cb,n3,n2,n1,id,key)

c** ca(2,i1,i2,i3) <- cb(2,i1,i2,i3)

      call trans_in_c4(cb,ca,id,n2,n3)   

      return

      end  

      subroutine c3ft_mask_for_out(ca,cb,id,n1,n2,n3,
     &                        cw1,cw2,cw3,
     &                        kc1,kc2,kc3,nbank,key)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $
c**
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)

      real*8 ca (0:id-1,0:n2-1,0:n3-1,2)
      real*8 cb (0:id-1,0:n2-1,0:n3-1,2)

c     complex*16 ca (0:id-1,0:n2-1,0:n3-1)
c     complex*16 cb (0:id-1,0:n2-1,0:n3-1)

      complex*16 cw1(0:n1-1)
      complex*16 cw2(0:n2-1)
      complex*16 cw3(0:n3-1)

      integer kc1,kc2,kc3
      integer nbank

      if(kc1.eq.0) then
           mm1 = n1
      else
           mm1 = kc1*2+1
      end if

      if(kc2.eq.0) then
           mm2 = n2
      else
           mm2 = kc2*2+1
      end if

      if(kc3.eq.0) then
           mm3 = n3
      else
           mm3 = kc3*2+1
      end if

c** cb(i1,i2,i3,2) <- ca(2,i1,i2,i3)

      call trans_out_c0(ca,cb,id,n1,n2,n3,key)

c** Length:n3 CFFT  

      call setlp(n3,lp2,lp3,lp4,lp5,lp8)
      call fsb_lp235(cb,ca,cw3,n3,n1*n2,lp2,lp3,lp4,lp5,lp8)

c** cb(i1,i3,i2,2) <- ca(i1,i2,i3,2)

      call trans_out_c1(ca,cb,n1,n2,n3,mm3,kc3)

c** Length:n2 CFFT    

      call setlp(n2,lp2,lp3,lp4,lp5,lp8)
      call fsb_lp235(cb,ca,cw2,n2,n1*mm3,lp2,lp3,lp4,lp5,lp8) 

c** cb(i3,i2,i1,2) <- ca(i1,i3,i2,2)

      call trans_out_c2(ca,cb,n1,n2,mm2,mm3,kc2)

c** Length:n1 CFFT  

      call setlp(n1,lp2,lp3,lp4,lp5,lp8)

      call fsb_lp235(cb,ca,cw1,n1,mm3*mm2,lp2,lp3,lp4,lp5,lp8)

c** cb(2,i1,i2,i3) <- ca(i3,i2,i1,2)

      call trans_out_c3(ca,cb,id,n1,n2,n3,mm2,mm3,kc1,kc2,kc3,key)

c** ca(2,i1,i2,i3) <- cb(2,i1,i2,i3)

      call trans_out_c4(cb,ca,id,n2,n3)

      return
      end
      subroutine trans_in_c0(ca,cb,id,n1,n2,n3,mm1,mm2,kc1,kc2,kc3,key)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $

      implicit real*8(a-h,o-z)

      real*8 ca(2,0:id-1 ,0:n2-1 ,0:n3-1)
      real*8 cb(0:mm1-1,0:mm2-1,0:n3-1,2)

      real*8 parity
c*
      if(kc1.eq.0) then
         kk1 = n1-1
      else
         kk1 = kc1
      end if

      if(kc2.eq.0) then
         kk2 = n2-1
      else
         kk2 = kc2
      end if

      if(kc3.eq.0) then
         kk3 = n3-1
      else
         kk3 = kc3
      end if
c*
      if(key.eq.+1.or.key.eq.+2) then
         parity = +1.0d0
      else
         parity = -1.0d0
      end if

c*
      do 10 i3=0,n3-1
      do 10 i2=0,mm2-1
      do 10 i1=0,mm1-1

         cb(i1,i2,i3,1)=0.0d0
         cb(i1,i2,i3,2)=0.0d0

 10   continue
c*
      do 20 i3=0,kk3
      do 20 i2=0,kk2
      do 20 i1=0,kk1

         cb(i1,i2,i3,1) = ca(1,i1,i2,i3)
         cb(i1,i2,i3,2) = ca(2,i1,i2,i3) * parity

 20   continue

      if(kc1.ne.0) then

         do 30 i3=0,kk3
         do 30 i2=0,kk2
         do 30 i1=n1-kc1,n1-1
     
            ii1 = i1 + mm1 - n1
            cb(ii1,i2,i3,1) = ca(1,i1,i2,i3)
            cb(ii1,i2,i3,2) = ca(2,i1,i2,i3) * parity

 30      continue

      end if

      if(kc2.ne.0) then

         do 40 i3=0,kk3
         do 40 i2=n2-kc2,n2-1
         do 40 i1=0,kk1
     
            ii2 = i2 + mm2 - n2
            cb(i1,ii2,i3,1) = ca(1,i1,i2,i3)
            cb(i1,ii2,i3,2) = ca(2,i1,i2,i3) * parity

 40      continue

      end if

      if(kc3.ne.0) then

         do 50 i3=n3-kc3,n3-1
         do 50 i2=0,kk2
         do 50 i1=0,kk1
     
            ii3 = i3 
            cb(i1,i2,ii3,1) = ca(1,i1,i2,i3)
            cb(i1,i2,ii3,2) = ca(2,i1,i2,i3) * parity

 50      continue

      end if

      if(kc1.ne.0.and.kc2.ne.0) then

         do 60 i3=0,kk3
         do 60 i2=n2-kc2,n2-1
         do 60 i1=n1-kc1,n1-1
     
            ii1 = i1 + mm1 - n1
            ii2 = i2 + mm2 - n2
            cb(ii1,ii2,i3,1) = ca(1,i1,i2,i3)
            cb(ii1,ii2,i3,2) = ca(2,i1,i2,i3) * parity

 60      continue

      end if

      if(kc1.ne.0.and.kc3.ne.0) then

         do 70 i3=n3-kc3,n3-1
         do 70 i2=0,kk2
         do 70 i1=n1-kc1,n1-1
     
            ii1 = i1 + mm1 - n1
            ii3 = i3 
            cb(ii1,i2,ii3,1) = ca(1,i1,i2,i3)
            cb(ii1,i2,ii3,2) = ca(2,i1,i2,i3) * parity

 70      continue

      end if

      if(kc2.ne.0.and.kc3.ne.0) then

         do 80 i3=n3-kc3,n3-1
         do 80 i2=n2-kc2,n2-1
         do 80 i1=0,kk1
     
            ii2 = i2 + mm2 - n2
            ii3 = i3 
            cb(i1,ii2,ii3,1) = ca(1,i1,i2,i3)
            cb(i1,ii2,ii3,2) = ca(2,i1,i2,i3) * parity

 80      continue

      end if

      if(kc1.ne.0.and.kc2.ne.0.and.kc3.ne.0) then

         do 90 i3=n3-kc3,n3-1
         do 90 i2=n2-kc2,n2-1
         do 90 i1=n1-kc1,n1-1
     
            ii1 = i1 + mm1 - n1
            ii2 = i2 + mm2 - n2
            ii3 = i3 
            cb(ii1,ii2,ii3,1) = ca(1,i1,i2,i3)
            cb(ii1,ii2,ii3,2) = ca(2,i1,i2,i3) * parity

 90      continue

      end if

      return
      end
      subroutine trans_in_c1(ca,cb,mm1,n2,mm2,n3,kc2)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $

      implicit real*8(a-h,o-z)

      real*8 ca(0:mm1-1,0:mm2-1,0:n3-1,2)
      real*8 cb(0:mm1-1,0:n3-1 ,0:n2-1,2)

      if(kc2.eq.0) then
         kk2 = n2-1
      else
         kk2 = kc2
      end if

      do 10 i3=0,n3-1
      do 10 i2=0,kk2
      do 10 i1=0,mm1-1

         cb(i1,i3,i2,1) = ca(i1,i2,i3,1)
         cb(i1,i3,i2,2) = ca(i1,i2,i3,2)

  10  continue

      if(kc2.ne.0) then

         do 20 i2=kk2+1,n2-kc2-1
         do 20 i3=0,n3-1
         do 20 i1=0,mm1-1

            cb(i1,i3,i2,1) = 0.0d0
            cb(i1,i3,i2,2) = 0.0d0

  20     continue

         do 30 i2=n2-kc2,n2-1
         do 30 i3=0,n3-1
         do 30 i1=0,mm1-1

            ii2 = i2 + mm2 - n2
            cb(i1,i3,i2,1) = ca(i1,ii2,i3,1)
            cb(i1,i3,i2,2) = ca(i1,ii2,i3,2)

  30     continue

      end if

      return
      end
      subroutine trans_in_c2(ca,cb,n1,mm1,n2,n3,kc1)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $

      implicit real*8(a-h,o-z)

      real*8 ca(0:mm1-1  ,0:n3*n2-1,2)
      real*8 cb(0:n3*n2-1,0:n1-1   ,2)

      if(kc1.eq.0) then
         kk1 = n1-1
      else
         kk1 = kc1
      end if

c**
      do 10 i1 = 0,kk1
      do 10 i32 = 0,n3*n2-1

         cb(i32,i1,1) = ca(i1,i32,1)
         cb(i32,i1,2) = ca(i1,i32,2)

  10  continue

      if(kc1.ne.0) then

         do 20 i1 = kk1+1,n1-kc1-1
         do 20 i32 = 0,n3*n2-1

            cb(i32,i1,1) = 0.0d0
            cb(i32,i1,2) = 0.0d0

  20     continue

         do 30 i1 = n1-kc1,n1-1
         do 30 i32 = 0,n3*n2-1

            ii1 = i1 + mm1 - n1
            cb(i32,i1,1) = ca(ii1,i32,1)
            cb(i32,i1,2) = ca(ii1,i32,2)

  30     continue

      end if

      return
      end
      subroutine trans_in_c3(ca,cb,n3,n2,n1,id,key)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $

      implicit real*8(a-h,o-z)

      real*8 ca(0:n3-1,0:n2-1,0:n1-1,2)
      real*8 cb(2,0:id-1,0:n2-1,0:n3-1)

      real*8 parity
c*
      if(key.eq.+1.or.key.eq.+2) then
         parity = +1.0d0
      else
         parity = -1.0d0
      end if
c*
      do 10 i3=0,n3-1
      do 10 i2=0,n2-1
      do 10 i1=n1,id-1

         cb(1,i1,i2,i3) = 0.0d0
         cb(2,i1,i2,i3) = 0.0d0

  10  continue

      do 20 i1=0,n1-1
      do 20 i2=0,n2-1
      do 20 i3=0,n3-1

         cb(1,i1,i2,i3) = ca(i3,i2,i1,1)
         cb(2,i1,i2,i3) = ca(i3,i2,i1,2) * parity

  20  continue

      return
      end

      subroutine trans_in_c4(cb,ca,id,n2,n3)   
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $

      implicit real*8(a-h,o-z)

      real*8 cb(id*n2*n3*2)
      real*8 ca(id*n2*n3*2)

      do 10 ii=1,id*n2*n3*2

         ca(ii) = cb(ii)

 10   continue

      return
      end

      subroutine trans_out_c0(ca,cb,id,n1,n2,n3,key)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $

      implicit real*8(a-h,o-z)

      real*8 ca(2,id,n2*n3)
      real*8 cb(n1,n2*n3,2)

      real*8 parity
c*
      if(key.eq.+1.or.key.eq.+2) then
         parity = +1.0d0
      else
         parity = -1.0d0
      end if
c*
      do 10 i23=1,n2*n3
      do 10 i1=1,n1

         cb(i1,i23,1) = ca(1,i1,i23)
         cb(i1,i23,2) = ca(2,i1,i23) * parity

 10   continue

      return
      end 

      subroutine trans_out_c1(ca,cb,n1,n2,n3,mm3,kc3)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $

      implicit real*8(a-h,o-z)

      real*8 ca(0:n1-1,0:n2-1 ,0:n3-1,2)
      real*8 cb(0:n1-1,0:mm3-1,0:n2-1,2)

      if(kc3.eq.0) then
         kk3 = n3-1
      else
         kk3 = kc3  
      end if

      do 10 i3=0,kk3
      do 10 i2=0,n2-1
      do 10 i1=0,n1-1

         cb(i1,i3,i2,1) = ca(i1,i2,i3,1)
         cb(i1,i3,i2,2) = ca(i1,i2,i3,2)

 10   continue

      if(kc3.ne.0) then

         do 20 i3=n3-kc3,n3-1
         do 20 i2=0,n2-1
         do 20 i1=0,n1-1

            ii3 = i3 + mm3 - n3
            cb(i1,ii3,i2,1) = ca(i1,i2,i3,1)
            cb(i1,ii3,i2,2) = ca(i1,i2,i3,2)

 20      continue

      end if

      return
      end
      subroutine trans_out_c2(ca,cb,n1,n2,mm2,mm3,kc2)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $

      implicit real*8(a-h,o-z)

      real*8 ca(0:n1-1 ,0:mm3-1,0:n2-1,2)
      real*8 cb(0:mm3-1,0:mm2-1,0:n1-1,2)

      if(kc2.eq.0) then
         kk2 = n2-1
      else
         kk2 = kc2
      end if

      do 10 i2 = 0,kk2
      do 10 i3 = 0,mm3-1
      do 10 i1 = 0,n1-1
   
         cb(i3,i2,i1,1) = ca(i1,i3,i2,1)
         cb(i3,i2,i1,2) = ca(i1,i3,i2,2)

  10  continue

      if(kc2.ne.0) then

      do 20 i2 = n2-kc2,n2-1
      do 20 i3 = 0,mm3-1
      do 20 i1 = 0,n1-1
   
         ii2 = i2 + mm2 - n2
         cb(i3,ii2,i1,1) = ca(i1,i3,i2,1)
         cb(i3,ii2,i1,2) = ca(i1,i3,i2,2)

  20  continue

      end if

      return
      end
      subroutine trans_out_c3(ca,cb,id,n1,n2,n3,mm2,mm3,kc1,kc2,kc3,key)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $

      implicit real*8(a-h,o-z)

      real*8 ca(0:mm3-1,0:mm2-1,0:n1-1,2)
      real*8 cb(2,0:id-1 ,0:n2-1 ,0:n3-1)

      real*8 parity
c*
      if(key.eq.+1.or.key.eq.+2) then
         parity = +1.0d0
      else
         parity = -1.0d0
      end if
c*

      do 1 i3=0,n3-1
      do 1 i2=0,n2-1
      do 1 i1=0,id-1

         cb(1,i1,i2,i3) = 0.0d0
         cb(2,i1,i2,i3) = 0.0d0

  1   continue

      if(kc1.eq.0) then
         kk1 = n1-1
      else
         kk1 = kc1
      end if

      if(kc2.eq.0) then
         kk2 = n2-1
      else
         kk2 = kc2
      end if

      if(kc3.eq.0) then
         kk3 = n3-1
      else
         kk3 = kc3
      end if

      do 10 i1 = 0,kk1
      do 10 i3 = 0,kk3
      do 10 i2 = 0,kk2

         cb(1,i1,i2,i3) = ca(i3,i2,i1,1)
         cb(2,i1,i2,i3) = ca(i3,i2,i1,2) * parity

  10  continue

      if(kc1.ne.0) then

      do 20 i1 = n1-kc1,n1-1
      do 20 i3 = 0,kk3
      do 20 i2 = 0,kk2

         ii1 = i1
         cb(1,i1,i2,i3) = ca(i3,i2,ii1,1)
         cb(2,i1,i2,i3) = ca(i3,i2,ii1,2) * parity

  20  continue

      end if

      if(kc2.ne.0) then

      do 30 i1 = 0,kk1
      do 30 i3 = 0,kk3
      do 30 i2 = n2-kc2,n2-1

         ii2 = i2 + mm2 - n2
         cb(1,i1,i2,i3) = ca(i3,ii2,i1,1)
         cb(2,i1,i2,i3) = ca(i3,ii2,i1,2) * parity

  30  continue

      end if

      if(kc3.ne.0) then

      do 40 i1 = 0,kk1
      do 40 i3 = n3-kc3,n3-1
      do 40 i2 = 0,kk2

         ii3 = i3 + mm3 - n3
         cb(1,i1,i2,i3) = ca(ii3,i2,i1,1)
         cb(2,i1,i2,i3) = ca(ii3,i2,i1,2) * parity

  40  continue

      end if

      if(kc1.ne.0.and.kc2.ne.0) then

      do 50 i1 = n1-kc1,n1-1
      do 50 i3 = 0,kk3
      do 50 i2 = n2-kc2,n2-1

         ii1 = i1
         ii2 = i2 + mm2 - n2
         cb(1,i1,i2,i3) = ca(i3,ii2,ii1,1)
         cb(2,i1,i2,i3) = ca(i3,ii2,ii1,2) * parity

  50  continue

      end if

      if(kc1.ne.0.and.kc3.ne.0) then

      do 60 i1 = n1-kc1,n1-1
      do 60 i3 = n3-kc3,n3-1
      do 60 i2 = 0,kk2

         ii1 = i1
         ii3 = i3 + mm3 - n3
         cb(1,i1,i2,i3) = ca(ii3,i2,ii1,1)
         cb(2,i1,i2,i3) = ca(ii3,i2,ii1,2) * parity

  60  continue

      end if

      if(kc2.ne.0.and.kc3.ne.0) then

      do 70 i1 = 0,kk1
      do 70 i3 = n3-kc3,n3-1
      do 70 i2 = n2-kc2,n2-1

         ii2 = i2 + mm2 - n2
         ii3 = i3 + mm3 - n3
         cb(1,i1,i2,i3) = ca(ii3,ii2,i1,1)
         cb(2,i1,i2,i3) = ca(ii3,ii2,i1,2) * parity

  70  continue

      end if

      if(kc1.ne.0.and.kc2.ne.0.and.kc3.ne.0) then

      do 80 i1 = n1-kc1,n1-1
      do 80 i3 = n3-kc3,n3-1
      do 80 i2 = n2-kc2,n2-1

         ii1 = i1
         ii2 = i2 + mm2 - n2
         ii3 = i3 + mm3 - n3
         cb(1,i1,i2,i3) = ca(ii3,ii2,ii1,1)
         cb(2,i1,i2,i3) = ca(ii3,ii2,ii1,2) * parity

  80  continue

      end if

      return
      end
      subroutine trans_out_c4(cb,ca,id,n2,n3)
c $Id: jrcat_ft_stm.f,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $

      implicit real*8(a-h,o-z)

      real*8 ca(id*n2*n3*2)
      real*8 cb(id*n2*n3*2)

c*
      do 10 ii=1,id*n2*n3*2

         ca(ii) = cb(ii)

 10   continue

      return
      end
