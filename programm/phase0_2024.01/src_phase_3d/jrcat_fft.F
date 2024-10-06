#ifdef HIUX
*OPTION MP(P(0))
#endif
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
#ifdef WF_JRCATFFT
      subroutine jrcat_r3ft(ra,work,id,n1,n2,n3,cw1,cw2,cw3,
     &                      kc1,kc2,kc3,key)
c $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $
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

c********************************************************************************
c*   key = 0 : set up of FFT.
c*             ( Calculation of cw1,cw2 and cw3. )
c********************************************************************************
c*   n1,n2 and n3 must include radix 2, and they may include 3 or 5.
c*   id must be larger than or equal to n1 + 2.
c*
c********************************************************************************
c*  key = -1 & -2 : F.T. of real 3d array ra(0:n1-1,0:n2-1,0:n3-1) .
c********************************************************************************
c*
c*   The result(complex 'cb') is overwritten on 'ra'.
c*   ( i.e. equivalence(ra,cb) ).
c*   
c*     cb(k1,k2,k3) 
c*   = sum(ix=0:n1-1,iy=0:n2-1,iz=0:n3-1) 
c*     ra(i1,i2,i3) exp [- 2 * pai * ci * ( k1*i1/n1 + k2*i2/n2 + k3*i3/n3 ) ]
c*     
c*     where 0 <= k1 <= n1/2, 0 <= k2 <= n2-1, 0 <= k3 <= n3-1.
c*
c********************************************************************************
c*  key = +1 & +2: Inverse F.T. of copmplex 3d array cb(0:n1/2,0:n2-1,0:n3-1) .
c********************************************************************************
c*   Transform from 'cb' to 'ra'.
c*
c********************************************************************************
c*  kc1,kc2 kc3 are integer parameters specifying the mask for
c*  input data or output data.
c*
c*     If kc1.ge.n1/2.or.kc1.lt.0, then kc1 is forced to be 0. 
c*     If kc2.ge.n2/2.or.kc2.lt.0, then kc2 is forced to be 0. 
c*     If kc3.ge.n3/2.or.kc3.lt.0, then kc3 is forced to be 0. 
c*
c********************************************************************************
c** REMARK on Mask.
c********************************************************************************
c*
c*  key = -1 : Mask is for Output data ( cb ),
c*  key = +1 : Mask is for Input  data ( cb ),
c*  key = -2 : Mask is for Input  data ( ca ),
c*  key = +2 : Mask is for Output data ( ca ),
c*
c*     where ca(0:idh/2,0:n2-1,0:n3-1) and equivalence(ra,ca).
c*
c********************************************************************************
c*     Mask for Output data:
c********************************************************************************
c*
c*     If kc1.eq.0.and.kc2.eq.0.and.kc3.eq.0, then the full wavenumber
c*     range is calculated.
c*
c*     If kc1.ne.0, then only "0 <= k1 <= kc1" range is calculated.
c*
c*     If kc2.ne.0, then only "0 <= k2 <= kc2" and "n2-kc2 <= k2 <= n2-1" 
c*     range is calculated.
c*
c*     If kc3.ne.0, then only "0 <= k3 <= kc3" and "n3-kc3 <= k3 <= n3-1" 
c*     range is calculated.
c*
c********************************************************************************
c*     Mask for Input data:
c********************************************************************************
c*
c*     If kc1.eq.0.and.kc2.eq.0.and.kc3.eq.0, then the full wavenumber
c*     range is used to obtain real array 'ra'.
c*
c*     If kc1.ne.0, then only "0 <= k1 <= kc1" range is used.
c*
c*     If kc2.ne.0, then only "0 <= k2 <= kc2" and "n2-kc2 <= k2 <= n2-1" 
c*     range is used.
c*
c*     If kc3.ne.0, then only "0 <= k3 <= kc3" and "n3-kc3 <= k3 <= n3-1" 
c*     range is used.
c*
c********************************************************************************
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
      call gettod(t_start)

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
      call gettod(t_end)
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
         write(6,9001) ecpu, icount, tcpu
 9001    format(1H ,' <<CPU TIME for JRCAT_R3FT = ',f11.3
     &        ,'(Sec.) = (',i6,' *',f10.5,')>>')
      endif
      return
      end  

      subroutine r3ft_0(cw1,cw2,cw3,n1,n2,n3)
c $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $
c**
      use m_Const_Parameters, only  : PAI
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

c**         pai=3.14159265358979323846d0
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
c $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $
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

c** Packin from ca(2,0:idh-1,0:n2-1,0:n3-1) to cb(0:mm1-1,0:mm2-1,0:n3-1,2)

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
c $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $
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

c$$$      if(kc1.eq.0) then
c$$$           mm1 = n1h+1 
c$$$      else
c$$$           mm1 = kc1+1
c$$$      end if

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
c $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $
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
c $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $
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

c$$$      if(kc1.eq.0) then
c$$$           mm1 = n1h+1 
c$$$      else
c$$$           mm1 = kc1+1
c$$$      end if

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
c $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $
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
c $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $
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
c $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $
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
c $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $
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
c $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $
c---------------------------------------------
      use m_Const_Parameters, only : PAI
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

c**      pai  =3.14159265358979323846d0
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
c $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $
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
c $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $
      use m_Const_Parameters, only : PAI
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)
      complex*16 cw(0:n-1)

c**      pai=3.14159265358979323846d0
      pc =2.0d0*pai/dble(n)

      do 10 i=0,n-1
         tr=dcos(pc*dble(i))
         ti=dsin(pc*dble(i))
         cw(i)=dcmplx(tr,ti)
   10 continue

      return
      end
      subroutine setlp(n,lp2,lp3,lp4,lp5,lp8)
c $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $
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
c $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $

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
c $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $

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
c $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $

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
c $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $

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
c $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $

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
c $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $

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
c $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $

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
c $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $

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
c $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $

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
c $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $

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
c $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $

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
c $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $

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
c $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $

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
c $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $

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
c $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $
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
c $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $

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
c $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $

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
! $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $

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
! $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $

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
! $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $

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
! $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $

c********************************************************************************
c*    3 dimensional complex to complex FFT with mask
c*    either for input array or output array,
c*    developed by T.Sanada (sanada@think.com) on 1996.01.04.
c*    It is specified for Fujitsu VPP500.
c********************************************************************************

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

c********************************************************************************
c*   key = 0 : set up of FFT.
c*             ( Calculation of cw1,cw2 and cw3. )
c********************************************************************************
c*   n1,n2 and n3 must include radix 2, and they may include 3 or 5.
c*   id must be n1 + 2.
c*
c********************************************************************************
c*  key = -1 & -2 : F.T. of complex 3d array ca(0:n1-1,0:n2-1,0:n3-1) .
c********************************************************************************
c*
c*   The result(complex 'cb') is overwritten on 'ca'.
c*   
c*     cb(k1,k2,k3) 
c*   = sum(ix=0:n1-1,iy=0:n2-1,iz=0:n3-1) 
c*     ca(i1,i2,i3) exp [- 2 * pai * ci * ( k1*i1/n1 + k2*i2/n2 + k3*i3/n3 ) ]
c*     
c*     where 0 <= k1 <= n1/2, 0 <= k2 <= n2-1, 0 <= k3 <= n3-1.
c*
c********************************************************************************
c*  key = +1 & +2: Inverse F.T. of copmplex 3d array cb(0:n1/2,0:n2-1,0:n3-1) .
c********************************************************************************
c*   Transform from 'cb' to 'ra'.
c*
c********************************************************************************
c*  kc1,kc2 kc3 are integer parameters specifying the mask for
c*  input data or output data.
c*
c*     If kc1.ge.n1/2.or.kc1.lt.0, then kc1 is forced to be 0. 
c*     If kc2.ge.n2/2.or.kc2.lt.0, then kc2 is forced to be 0. 
c*     If kc3.ge.n3/2.or.kc3.lt.0, then kc3 is forced to be 0. 
c*
c********************************************************************************
c** REMARK on Mask.
c********************************************************************************
c*
c*  key = -1 : Mask is for Output data ( cb ),
c*  key = +1 : Mask is for Input  data ( cb ),
c*  key = -2 : Mask is for Input  data ( ca ),
c*  key = +2 : Mask is for Output data ( ca ),
c*
c*     where ca(0:idh/2,0:n2-1,0:n3-1) and equivalence(ra,ca).
c*
c********************************************************************************
c*     Mask for Output data:
c********************************************************************************
c*
c*     If kc1.eq.0.and.kc2.eq.0.and.kc3.eq.0, then the full wavenumber
c*     range is calculated.
c*
c*     If kc1.ne.0, then only "0 <= k1 <= kc1" range is calculated.
c*
c*     If kc2.ne.0, then only "0 <= k2 <= kc2" and "n2-kc2 <= k2 <= n2-1" 
c*     range is calculated.
c*
c*     If kc3.ne.0, then only "0 <= k3 <= kc3" and "n3-kc3 <= k3 <= n3-1" 
c*     range is calculated.
c*
c********************************************************************************
c*     Mask for Input data:
c********************************************************************************
c*
c*     If kc1.eq.0.and.kc2.eq.0.and.kc3.eq.0, then the full wavenumber
c*     range is used to obtain real array 'ra'.
c*
c*     If kc1.ne.0, then only "0 <= k1 <= kc1" range is used.
c*
c*     If kc2.ne.0, then only "0 <= k2 <= kc2" and "n2-kc2 <= k2 <= n2-1" 
c*     range is used.
c*
c*     If kc3.ne.0, then only "0 <= k3 <= kc3" and "n3-kc3 <= k3 <= n3-1" 
c*     range is used.
c*
c********************************************************************************
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
      call gettod(t_start)

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
      call gettod(t_end)
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
         write(6,9001) ecpu, icount, tcpu
 9001    format(1H ,' <<CPU TIME for JRCAT_C3FT = ',f11.3
     &        ,'(Sec.) = (',i6,' *',f10.5,')>>')
      endif
      return
      end  

      subroutine c3ft_0(cw1,cw2,cw3,n1,n2,n3)
! $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $
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
! $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $
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

c$$$      if(kc3.eq.0) then
c$$$           mm3 = n3
c$$$      else
c$$$           mm3 = kc3*2+1
c$$$      end if

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
! $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $
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

c$$$      if(kc1.eq.0) then
c$$$           mm1 = n1
c$$$      else
c$$$           mm1 = kc1*2+1
c$$$      end if

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
! $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $

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
! $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $

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
! $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $

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
! $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $

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
! $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $

      implicit real*8(a-h,o-z)

      real*8 cb(id*n2*n3*2)
      real*8 ca(id*n2*n3*2)

      do 10 ii=1,id*n2*n3*2

         ca(ii) = cb(ii)

 10   continue

      return
      end
      subroutine trans_out_c0(ca,cb,id,n1,n2,n3,key)
! $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $

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
! $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $

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
! $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $

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
! $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $

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
! $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $

      implicit real*8(a-h,o-z)

      real*8 ca(id*n2*n3*2)
      real*8 cb(id*n2*n3*2)

c*
      do 10 ii=1,id*n2*n3*2

         ca(ii) = cb(ii)

 10   continue

      return
      end
#else
!
c $Id: jrcat_fft.F 376 2014-06-17 07:48:31Z jkoga $
c---------------------------------------------
      subroutine jrcat_c3ft(ca,cwork,id,n1,n2,n3,cw1,cw2,cw3,
     &                      kc1,kc2,kc3,key)
c---------------------------------------------
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
      dimension ca    (id,n2,n3)
      dimension cwork (id+n2+n3,2)

      dimension cw1(0:n1-1)
      dimension cw2(0:n2-1)
      dimension cw3(0:n3-1)

      integer kc1,kc2,kc3
      integer key

c*******************************************************************************
c*  3 DIMENSIONAL FAST FOURIER TRANSFORM OF COMPLEX ARRAY WITH MASK,
c*                developed on 1997.6.5.
c*                modified  on 1998.5.18.
c*******************************************************************************
c*  
c*  The 3-dimensional array ca(id,n2,n3) is fourier-transfomed and
c*  its result is overwritten on ca. The size of fourier-transformed is
c*  (n1,n2,n3) and hence id must be larger than or equal to n1.
c*  Also n1,n2 and n3 must be radix 2,3,5 integers.
c*
c*  cwork is complex work array and the size must be at least (id+n2+n3)*2
c*
c*******************************************************************************
c*  How to set the value of key.
c*
c*(1) key = 0       : set up of FFT 
c*    
c*     Calculations of cw1,cw2 and cw3.
c*
c*(2) key = -1 & -2 : 
c*
c*     Forward F.T. of complex 3d array ca(0:id-1,0:n2-1,0:n3-1),
c*
c*     here the definition of fourier-transform is
c*  
c*       cb(k1,k2,k3) 
c*     = sum(i1=0:n1-1,i2=0:n2-1,i3=0:n3-1) 
c*       ca(i1,i2,i3) exp [- 2 * pai * ci * ( k1*i1/n1 + k2*i2/n2 + k3*i3/n3 ) ]
c*     
c*        where 0 <= k1 <= n1-1, 0 <= k2 <= n2-1, 0 <= k3 <= n3-1,
c*              and ci is the imaginary unit.
c*
c*    The result(complex array cb(0:id-1,0:n2-1,0:n3-1)) is overwritten on 'ca'.
c*
c*(3) key = +1 & +2 : 
c*
c*     Inverse F.T. of complex 3d array cb(0:id-1,0:n2-1,0:n3-1) 
c*                   to ca(0:id-1,0:n2-1,0:n3-1).
c*
c*******************************************************************************
c*  How to set the mask.
c*
c*  kc1,kc2,kc3 and key are integer parameters specifying the mask for
c*  input data or output data.
c*
c*     key = -1 : Mask is for Output data ( cb ),
c*     key = +1 : Mask is for Input  data ( cb ),
c*     key = -2 : Mask is for Input  data ( ca ),
c*     key = +2 : Mask is for Output data ( ca ),
c*
c*     If kc1.ge.n1/2.or.kc1.lt.0, then kc1 is forced to be 0. 
c*     If kc2.ge.n2/2.or.kc2.lt.0, then kc2 is forced to be 0. 
c*     If kc3.ge.n3/2.or.kc3.lt.0, then kc3 is forced to be 0. 
c*
c*     where ca(0:id-1,0:n2-1,0:n3-1) is equivalent to cb.
c*
c*******************************************************************************
c*  Mask for Output data:
c*
c*     If kc1.eq.0.and.kc2.eq.0.and.kc3.eq.0, then the full 
c*     range is calculated.
c*
c*     If kc1.ne.0, then only "0 <= k1 <= kc1" and "n1-kc1 <= k1 <= n1-1" 
c*     range is calculated.
c*
c*     If kc2.ne.0, then only "0 <= k2 <= kc2" and "n2-kc2 <= k2 <= n2-1" 
c*     range is calculated.
c*
c*     If kc3.ne.0, then only "0 <= k3 <= kc3" and "n3-kc3 <= k3 <= n3-1" 
c*     range is calculated.
c*
c*******************************************************************************
c*  Mask for Input data:
c*
c*     If kc1.eq.0.and.kc2.eq.0.and.kc3.eq.0, then the full 
c*     range is used to obtain complex array 'ca'.
c*
c*     If kc1.ne.0, then only "0 <= k1 <= kc1" and "n1-kc1 <= k1 <= n1-1" 
c*     range is calculated.
c*
c*     If kc2.ne.0, then only "0 <= k2 <= kc2" and "n2-kc2 <= k2 <= n2-1" 
c*     range is used.
c*
c*     If kc3.ne.0, then only "0 <= k3 <= kc3" and "n3-kc3 <= k3 <= n3-1" 
c*     range is used.
c*
c*******************************************************************************
c
      call gettod(t_start)

c** Set up of FFT.
      if(key.eq.0) then
         call c3ft_0(cw1,cw2,cw3,n1,n2,n3)
         return
      end if
c**
      if(kc1.ge.n1/2.or.kc1.lt.0) then 
         write(6,*) 'warning in JRCAT_C3FT: kc1 is irrelevant.'
         write(6,*) 'kc1 has changed to be 0.'
         kc1 = 0
      end if
c*
      if(kc2.ge.n2/2.or.kc2.lt.0) then 
         write(6,*) 'warning in JRCAT_C3FT: kc2 is irrelevant.'
         write(6,*) 'kc2 has changed to be 0.'
         kc2 = 0
      end if
c*
      if(kc3.ge.n3/2.or.kc3.lt.0) then 
         write(6,*) 'warning in JRCAT_C3FT: kc3 is irrelevant.'
         write(6,*) 'kc3 has changed to be 0.'
         kc3 = 0
      end if
c**
      if(key.eq.-1.or.key.eq.-2) then
         call c3ft_conj(cw1,n1)
         call c3ft_conj(cw2,n2)
         call c3ft_conj(cw3,n3)
      end if
c********************************************************************
      if(key.eq.-1.or.key.eq.+2) then
c********************************************************************
c*  Forword or Backword FFT with mask for output.

         call c3ft_mask_for_out(ca,cwork,id,n1,n2,n3,cw1,cw2,cw3,
     &                     kc1,kc2,kc3,key)

c********************************************************************
      else 
c********************************************************************
c*  Forword or Backword FFT with mask for input.

         call c3ft_mask_for_in (ca,cwork,id,n1,n2,n3,cw1,cw2,cw3,
     &                     kc1,kc2,kc3,key)

c********************************************************************
      end if
c********************************************************************
      if(key.eq.-1.or.key.eq.-2) then
         call c3ft_conj(cw1,n1)
         call c3ft_conj(cw2,n2)
         call c3ft_conj(cw3,n3)
      end if
c**  
c -------->  CPU cost estimation
      call gettod(t_end)
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
         write(6,9001) ecpu, icount, tcpu
 9001    format(1H ,' <<CPU TIME for JRCAT_C3FT = ',f11.3
     &        ,'(Sec.) = (',i6,' *',f10.5,')>>')
      endif
      return
      end  

c---------------------------------------------
      subroutine jrcat_r3ft(ra,work,id,n1,n2,n3,cw1,cw2,cw3,
     &                      kc1,kc2,kc3,key)
c---------------------------------------------
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

      dimension ra   (id,n2,n3)
      dimension work (id+n2+n3,2)

      dimension cw1(0:n1-1)
      dimension cw2(0:n2-1)
      dimension cw3(0:n3-1)

      integer kc1,kc2,kc3
      integer key 

c*******************************************************************************
c*  3 DIMENSIONAL FAST FOURIER TRANSFORM OF REAL ARRAY WITH MASK,
c*                developed on 1997.6.5.
c*******************************************************************************
c*
c*  The 3-dimensional array ra(id,n2,n3) is fourier-transfomed and
c*  its result is overwritten on ra. The size of fourier-transformed is
c*  (n1,n2,n3) and n1 axis is the complex conjugate axis. id must be n1 + 2.
c*  Also n1 must be an even number and n1,n2 and n3 must be radix 2,3,5
c*  integers.
c*
c*  work is real work array and the size must be at least (id+n2+n3)*2.
c*
c*******************************************************************************
c*  How to set the value of key.
c*
c*(1) key = 0 : set up of FFT.
c*
c*    Calculation of cw1,cw2 and cw3. 
c*
c*(2) key = -1 & -2 :
c*
c*    F.T. of real 3d array ra(0:id-1,0:n2-1,0:n3-1),
c*
c*    here the definition of fourier-transform is
c*
c*     cb(k1,k2,k3) 
c*   = sum(i1=0:n1-1,i2=0:n2-1,i3=0:n3-1) 
c*     ra(i1,i2,i3) exp [- 2 * pai * ci * ( k1*i1/n1 + k2*i2/n2 + k3*i3/n3 ) ]
c*     
c*     where 0 <= k1 <= n1/2, 0 <= k2 <= n2-1, 0 <= k3 <= n3-1.
c*           and ci is the imaginary unit.
c*
c*  The result(complex array cb(0:id/2-1,0:n2-1,0:n3-1)) is overwritten on 'ra'.
c*
c*(3) key = +1 & +2 :
c*
c*     Inverse F.T. of complex 3d array cb(0:id/2-1,0:n2-1,0:n3-1)
c*                   to ra(0:id-1,0:n2-1,0:n3-1).
c*
c*******************************************************************************
c*  How to set the mask.
c*
c*  kc1,kc2 kc3 and key are integer parameters specifying the mask for
c*  input data or output data.
c*
c*     key = -1 : Mask is for Output data ( cb ),
c*     key = +1 : Mask is for Input  data ( cb ),
c*     key = -2 : Mask is for Input  data ( ca ),
c*     key = +2 : Mask is for Output data ( ca ),
c*
c*     If kc1.ge.n1/2.or.kc1.lt.0, then kc1 is forced to be 0. 
c*     If kc2.ge.n2/2.or.kc2.lt.0, then kc2 is forced to be 0. 
c*     If kc3.ge.n3/2.or.kc3.lt.0, then kc3 is forced to be 0. 
c*
c*  where ca(0:id/2-1,0:n2-1,0:n3-1) are equivalent to ra(0:id-1,0:n2-1,0:n3-1).
c*
c*******************************************************************************
c*  Mask for Output data:
c*
c*     If kc1.eq.0.and.kc2.eq.0.and.kc3.eq.0, then the full 
c*     range is calculated.
c*
c*     If kc1.ne.0, then only "0 <= k1 <= kc1" range is calculated.
c*
c*     If kc2.ne.0, then only "0 <= k2 <= kc2" and "n2-kc2 <= k2 <= n2-1" 
c*     range is calculated.
c*
c*     If kc3.ne.0, then only "0 <= k3 <= kc3" and "n3-kc3 <= k3 <= n3-1" 
c*     range is calculated.
c*
c*******************************************************************************
c*  Mask for Input data:
c*
c*     If kc1.eq.0.and.kc2.eq.0.and.kc3.eq.0, then the full wavenumber
c*     range is used to obtain real array 'ra'.
c*
c*     If kc1.ne.0, then only "0 <= k1 <= kc1" range is used.
c*
c*     If kc2.ne.0, then only "0 <= k2 <= kc2" and "n2-kc2 <= k2 <= n2-1" 
c*     range is used.
c*
c*     If kc3.ne.0, then only "0 <= k3 <= kc3" and "n3-kc3 <= k3 <= n3-1" 
c*     range is used.
c*
c*******************************************************************************
c
      call gettod(t_start)

      if(mod(n1,2).ne.0) then
         write(6,*) 'ERROR in JRCAT_R3FT: n1 is not an even number.'
         stop
      end if
c*******************************************************************************
      if(key.eq.0) then
c** set up of FFT.
         call r3ft_0(cw1,cw2,cw3,n1,n2,n3)
         return
      end if
c*******************************************************************************
      if(id.ne.n1+2) then
         write(6,*) 'ERROR in JRCAT_R3FT: id is not n1+2.'
         stop
      end if
c*******************************************************************************
      if(kc1.ge.n1/2.or.kc1.lt.0) then 
         write(6,*) 'warning in JRCAT_R3FT: kc1 is irrelevant.'
         write(6,*) 'kc1 has changed to be 0.'
         kc1 = 0
      end if

      if(kc2.ge.n2/2.or.kc2.lt.0) then 
         write(6,*) 'warning in JRCAT_R3FT: kc2 is irrelevant.'
         write(6,*) 'kc2 has changed to be 0.'
         kc2 = 0
      end if

      if(kc3.ge.n3/2.or.kc3.lt.0) then 
         write(6,*) 'warning in JRCAT_R3FT: kc3 is irrelevant.'
         write(6,*) 'kc3 has changed to be 0.'
         kc3 = 0
      end if
c*******************************************************************************
      if(key.eq.-1) then  
c** Forword FFT with mask for output.
         call r3ft_f_mask_for_out(ra,work,id,n1,n2,n3,cw1,cw2,cw3,
     &                            kc1,kc2,kc3)
      else if(key.eq.-2) then
c** Forword FFT with mask for input.
         call r3ft_f_mask_for_in (ra,work,id,n1,n2,n3,cw1,cw2,cw3,
     &                            kc1,kc2,kc3)
      else if(key.eq.+1) then
c** Backword FFT with mask for input.
         call r3ft_b_mask_for_in (ra,work,id,n1,n2,n3,cw1,cw2,cw3,
     &                            kc1,kc2,kc3)
      else if(key.eq.+2) then
c** Backword FFT with mask for output.
         call r3ft_b_mask_for_out(ra,work,id,n1,n2,n3,cw1,cw2,cw3,
     &                            kc1,kc2,kc3)
      end if
c*******************************************************************************
c -------->  CPU cost estimation
      call gettod(t_end)
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
         write(6,9001) ecpu, icount, tcpu
 9001    format(1H ,' <<CPU TIME for JRCAT_R3FT = ',f11.3
     &        ,'(Sec.) = (',i6,' *',f10.5,')>>')
      endif
      return
      end  
c---------------------------------------------
      subroutine c3ft_0(cw1,cw2,cw3,n1,n2,n3)
c---------------------------------------------
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)

      dimension cw1(0:n1-1)
      dimension cw2(0:n2-1)
      dimension cw3(0:n3-1)

         call setcw(cw1,n1)
         call setcw(cw2,n2)
         call setcw(cw3,n3) 

      return
      end  
c---------------------------------------------
      subroutine c3ft_conj(cw,nn)
c---------------------------------------------
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)
      dimension cw(nn) 

      do 10 ii=1,nn
         cw(ii) = dconjg(cw(ii))
  10  continue
c 10  end do

      return
      end 
c---------------------------------------------
      subroutine c3ft_mask_for_in(ca,cwork,id,n1,n2,n3,
     &                            cw1,cw2,cw3,kc1,kc2,kc3,key)
c---------------------------------------------
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)

      dimension ca   (0:id-1,0:n2-1,0:n3-1)
      dimension cwork(*)

      dimension cw1(0:n1-1)
      dimension cw2(0:n2-1)
      dimension cw3(0:n3-1)

      integer kc1,kc2,kc3
      integer key
      integer lp(8,3)
c**
      if(kc2.eq.0) then
           nk2 = n2
      else
           nk2 = kc2+1
      end if
      if(kc3.eq.0) then
           nk3 = n3
      else
           nk3 = kc3+1
      end if
c**
      call setlp(n1,lp(1,1))
      call setlp(n2,lp(1,2))
      call setlp(n3,lp(1,3))
c**
      call c3ft_zeroset(ca,id,n1,n2,n3,kc1,kc2,kc3)

c** Length:n1&n2 CFFT    
      do 10 i3=0,n3-1
         if(i3.le.nk3-1.or.i3.ge.n3-kc3) then
         do 11 i2=0,nk2-1
            call fsb_lp235(ca(0,i2,i3),cwork,cw1,n1,lp(1,1),key)
  11     continue
c 11     end do
         do 12 i2=n2-kc2,n2-1 
            call fsb_lp235(ca(0,i2,i3),cwork,cw1,n1,lp(1,1),key)
  12     continue
c 12     end do
         call fsb_lp235v(ca(0,0,i3),cwork,cw2,n2,id,1,n1,lp(1,2),key) 
         end if
  10  continue
c 10  end do

c** Length:n3 CFFT  
      call fsb_lp235v(ca,cwork,cw3,n3,id*n2,1,id*n2,lp(1,3),key)
c**
      return
      end
c---------------------------------------------
      subroutine c3ft_mask_for_out(ca,cwork,id,n1,n2,n3,
     &                             cw1,cw2,cw3,kc1,kc2,kc3,key)
c---------------------------------------------
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)

      dimension ca   (0:id-1,0:n2-1,0:n3-1)
      dimension cwork(*)

      dimension cw1(0:n1-1)
      dimension cw2(0:n2-1)
      dimension cw3(0:n3-1)

      integer kc1,kc2,kc3
      integer key         
      integer lp(8,3)
c**
      if(kc2.eq.0) then
           nk2 = n2
      else
           nk2 = kc2+1
      end if
      if(kc3.eq.0) then
           nk3 = n3
      else
           nk3 = kc3+1
      end if
c**
      call setlp(n1,lp(1,1))
      call setlp(n2,lp(1,2))
      call setlp(n3,lp(1,3))

c** Length:n3 CFFT  
      call fsb_lp235v(ca,cwork,cw3,n3,id*n2,1,id*n2,lp(1,3),key)

c** Length:n2 & n1 CFFT    
      do 10 i3=0,n3-1
         if(i3.le.nk3-1.or.i3.ge.n3-kc3) then
            call fsb_lp235v(ca(0,0,i3),cwork,cw2,n2,id,1,n1,lp(1,2),key) 
         do 11 i2=0,nk2-1
            call fsb_lp235(ca(0,i2,i3),cwork,cw1,n1,lp(1,1),key)
  11     continue
c 11     end do
         do 12 i2=n2-kc2,n2-1 
            call fsb_lp235(ca(0,i2,i3),cwork,cw1,n1,lp(1,1),key)
  12     continue
c 12     end do
         end if
  10  continue
c 10  end do
c**
      call c3ft_zeroset(ca,id,n1,n2,n3,kc1,kc2,kc3)
c**
      return
      end
c---------------------------------------------
      subroutine c3ft_zeroset(ca,id,n1,n2,n3,kc1,kc2,kc3)
c---------------------------------------------
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)
      dimension ca(0:id-1,0:n2-1,0:n3-1)

      do 10 i3=0,n3-1
      do 10 i2=0,n2-1
      do 10 i1=n1,id-1
         ca(i1,i2,i3) = (0.0d0,0.0d0)
 10   continue
c10   end do

      if(kc1.ne.0) then
      do 20 i3=0,n3-1
      do 20 i2=0,n2-1
      do 20 i1=kc1+1,n1-kc1-1
         ca(i1,i2,i3) = (0.0d0,0.0d0)
 20   continue
c20   end do
      end if

      if(kc2.ne.0) then
      do 30 i3=0,n3-1
      do 30 i2=kc2+1,n2-kc2-1
      do 30 i1=0,id-1
         ca(i1,i2,i3) = (0.0d0,0.0d0)
 30   continue
c30   end do
      end if

      if(kc3.ne.0) then
      do 40 i3=kc3+1,n3-kc3-1
      do 40 i2=0,n2-1
      do 40 i1=0,id-1
         ca(i1,i2,i3) = (0.0d0,0.0d0)
 40   continue
c40   end do
      end if

      return
      end 
c---------------------------------------------
      subroutine r3ft_0(cw1,cw2,cw3,n1,n2,n3)
c---------------------------------------------
      use m_Const_Parameters, only : PAI
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)

      dimension cw1(0:n1/2-1,2)
      dimension cw2(0:n2-1)
      dimension cw3(0:n3-1)
c**
      n1h = n1/2

      call setcw(cw1,n1h)
      call setcw(cw2,n2 )
      call setcw(cw3,n3 ) 

c**      pai=3.14159265358979323846d0
      pc =2.0d0*pai/dble(n1)
      ci=(0.0d0,1.0d0)

      do i=0,n1h-1
         tr=dcos(pc*dble(i))
         ti=dsin(pc*dble(i))
c        cw1(i,2)=ci*dcmplx(tr,ti)
         cw1(i,2)=dcmplx(-ti,tr)
      end do   
c**
      return
      end  

c---------------------------------------------
      subroutine r3ft_b_mask_for_in (ca,cwork,id,n1,n2,n3,
     &                               cw1,cw2,cw3,kc1,kc2,kc3)
c---------------------------------------------
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)

      dimension ca   (0:id/2-1,0:n2-1,0:n3-1)
      dimension cwork(*)

      dimension cw1(0:n1/2-1,2)
      dimension cw2(0:n2-1)
      dimension cw3(0:n3-1)

      integer kc1,kc2,kc3
      integer lp(8,3)
c**
      n1h = n1/2
      idh = id/2
c**
      if(kc1.eq.0) then
           nk1 = idh
      else
           nk1 = kc1+1
      end if

      if(kc3.eq.0) then
           nk3 = n3
      else
           nk3 = kc3+1
      end if
c**
      call setlp(n1h,lp(1,1))
      call setlp(n2 ,lp(1,2))
      call setlp(n3 ,lp(1,3))
c**
      call r3ft_zeroset(ca,id,n1,n2,n3,kc1,kc2,kc3)
    
c** Length:n2 CFFT
      do 10 i3=0,n3-1
         if(i3.le.nk3-1.or.i3.ge.n3-kc3) then
      call fsb_lp235v(ca(0,0,i3),cwork,cw2,n2,idh,1,nk1,lp(1,2),+1)
         end if
  10  end do

c** Length:n3 CFFT
      if(kc1.eq.0) then
         call fsb_lp235v(ca,cwork,cw3,n3,idh*n2,1,idh*n2,lp(1,3),+1)
      else
         do 20 nn1=1,idh*(n2-1)+1,idh
            nn2=nn1+nk1-1
         call fsb_lp235v(ca,cwork,cw3,n3,idh*n2,nn1,nn2,lp(1,3),+1)
  20     continue
c 20     end do
      end if

c** Length:n1h RFFT
      do 30 i3=0,n3-1
      do 30 i2=0,n2-1
         call trans_b(ca(0,i2,i3),cw1(0,2),n1h,0)
         call fsb_lp235(ca(0,i2,i3),cwork,cw1,n1h,lp(1,1),+1)
  30  continue
c 30  end do
c**  
      return
      end  

c---------------------------------------------
      subroutine r3ft_b_mask_for_out(ca,cwork,id,n1,n2,n3,
     &                               cw1,cw2,cw3,kc1,kc2,kc3)
c---------------------------------------------
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)

      dimension ca   (0:id/2-1,0:n2-1,0:n3-1)
      dimension cwork(*)

      dimension cw1(0:n1/2-1,2)
      dimension cw2(0:n2-1)
      dimension cw3(0:n3-1)

      integer kc1,kc2,kc3
      integer lp(8,3)
c**
      n1h = n1/2
      idh = id/2

      if(kc2.eq.0) then
           nk2 = n2
      else
           nk2 = kc2+1
      end if

      if(kc3.eq.0) then
           nk3 = n3
      else
           nk3 = kc3+1
      end if
c**
      call setlp(n1h,lp(1,1))
      call setlp(n2 ,lp(1,2))
      call setlp(n3 ,lp(1,3))

c** Length:n3 CFFT
      call fsb_lp235v(ca,cwork,cw3,n3,idh*n2,1,idh*n2,lp(1,3),+1)

      do 30 i3=0,n3-1
         if(i3.le.nk3-1.or.i3.ge.n3-kc3) then
c** Length:n2 CFFT
         call fsb_lp235v(ca(0,0,i3),cwork,cw2,n2,idh,1,idh,lp(1,2),+1)
      do 10 i2=0,nk2-1
c** Length:n1h RFFT
         call trans_b(ca(0,i2,i3),cw1(0,2),n1h,0)
         call fsb_lp235(ca(0,i2,i3),cwork,cw1,n1h,lp(1,1),+1)
  10  end do
      do 20 i2=n2-kc2,n2-1
c** Length:n1h RFFT
         call trans_b(ca(0,i2,i3),cw1(0,2),n1h,0)
         call fsb_lp235(ca(0,i2,i3),cwork,cw1,n1h,lp(1,1),+1)
  20  end do
         end if
  30  end do
c**
      call r3ft_zeroset(ca,id,n1,n2,n3,kc1,kc2,kc3)
c**
      return
      end  

c---------------------------------------------
      subroutine r3ft_f_mask_for_in (ca,cwork,id,n1,n2,n3,
     &                               cw1,cw2,cw3,kc1,kc2,kc3)
c---------------------------------------------
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)

      dimension ca   (0:id/2-1,0:n2-1,0:n3-1)
      dimension cwork(*)

      dimension cw1(0:n1/2-1,2)
      dimension cw2(0:n2-1)
      dimension cw3(0:n3-1)

      integer kc1,kc2,kc3
      integer lp(8,3)
c**
      n1h = n1/2
      idh = id/2

      if(kc2.eq.0) then
           nk2 = n2
      else
           nk2 = kc2+1
      end if

      if(kc3.eq.0) then
           nk3 = n3
      else
           nk3 = kc3+1
      end if
c**   
      call setlp(n1h,lp(1,1))
      call setlp(n2 ,lp(1,2))
      call setlp(n3 ,lp(1,3))
c**   
      call c3ft_conj(cw1,n1)
      call c3ft_conj(cw2,n2)
      call c3ft_conj(cw3,n3)
c**   
      call r3ft_zeroset(ca,id,n1,n2,n3,kc1,kc2,kc3)
c**   
      do 30 i3=0,n3-1
         if(i3.le.nk3-1.or.i3.ge.n3-kc3) then
      do 10 i2=0,nk2-1
c** Length:n1h RFFT
         call fsb_lp235(ca(0,i2,i3),cwork,cw1,n1h,lp(1,1),-1)
         call trans_f(ca(0,i2,i3),cw1(0,2),n1h,0)
  10  end do
      do 20 i2=n2-kc2,n2-1
c** Length:n1h RFFT
         call fsb_lp235(ca(0,i2,i3),cwork,cw1,n1h,lp(1,1),-1)
         call trans_f(ca(0,i2,i3),cw1(0,2),n1h,0)
  20  end do
c** Length:n2 CFFT  
         call fsb_lp235v(ca(0,0,i3),cwork,cw2,n2,idh,1,idh,lp(1,2),-1)   
         end if
  30  end do

c** Length:n3 CFFT  
      call fsb_lp235v(ca,cwork,cw3,n3,idh*n2,1,idh*n2,lp(1,3),-1)
c**   
      call c3ft_conj(cw1,n1)
      call c3ft_conj(cw2,n2)
      call c3ft_conj(cw3,n3)
c**   
      return
      end  

c---------------------------------------------
      subroutine r3ft_f_mask_for_out(ca,cwork,id,n1,n2,n3,
     &                               cw1,cw2,cw3,kc1,kc2,kc3)
c---------------------------------------------
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)

      dimension ca   (0:id/2-1,0:n2-1,0:n3-1)
      dimension cwork(*)

      dimension cw1(0:n1/2-1,2)
      dimension cw2(0:n2-1)
      dimension cw3(0:n3-1)

      integer kc1,kc2,kc3
      integer lp(8,3)
c**
      n1h = n1/2
      idh = id/2
c**
      if(kc1.eq.0) then
           nk1 = idh
      else
           nk1 = kc1+1
      end if

      if(kc3.eq.0) then
           nk3 = n3
      else
           nk3 = kc3+1
      end if
c**   
      call setlp(n1h,lp(1,1))
      call setlp(n2 ,lp(1,2))
      call setlp(n3 ,lp(1,3))
c**   
      call c3ft_conj(cw1,n1)
      call c3ft_conj(cw2,n2)
      call c3ft_conj(cw3,n3)
      
c** Length:n1h RFFT
      do 10 i3=0,n3-1
      do 10 i2=0,n2-1
         call fsb_lp235(ca(0,i2,i3),cwork,cw1,n1h,lp(1,1),-1)
         call trans_f(ca(0,i2,i3),cw1(0,2),n1h,kc1)
  10  continue
c 10  end do

c** Length:n3 CFFT  
      if(kc1.eq.0) then
         call fsb_lp235v(ca,cwork,cw3,n3,idh*n2,1,idh*n2,lp(1,3),-1)
      else
         do 20 nn1=1,idh*(n2-1)+1,idh
            nn2=nn1+nk1-1
         call fsb_lp235v(ca,cwork,cw3,n3,idh*n2,nn1,nn2,lp(1,3),-1)
  20     continue
c 20     end do
      end if

c** Length:n2 CFFT  
      do 30 i3=0,n3-1
         if(i3.le.nk3-1.or.i3.ge.n3-kc3) then
      call fsb_lp235v(ca(0,0,i3),cwork,cw2,n2,idh,1,nk1,lp(1,2),-1)   
         end if
  30  end do
c**   
      call r3ft_zeroset(ca,id,n1,n2,n3,kc1,kc2,kc3)
c**   
      call c3ft_conj(cw1,n1)
      call c3ft_conj(cw2,n2)
      call c3ft_conj(cw3,n3)
c**   
      return
      end  

c---------------------------------------------
      subroutine r3ft_zeroset(ca,id,n1,n2,n3,kc1,kc2,kc3)
c---------------------------------------------
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)
      dimension ca(0:id/2-1,0:n2-1,0:n3-1)

      idh = id/2

      if(kc1.ne.0) then
      do 20 i3=0,n3-1
      do 20 i2=0,n2-1
      do 20 i1=kc1+1,idh-1
         ca(i1,i2,i3) = (0.0d0,0.0d0)
 20   continue
c20   end do
      end if

      if(kc2.ne.0) then
      do 30 i3=0,n3-1
      do 30 i2=kc2+1,n2-kc2-1
      do 30 i1=0,idh-1
         ca(i1,i2,i3) = (0.0d0,0.0d0)
 30   continue
c30   end do
      end if

      if(kc3.ne.0) then
      do 40 i3=kc3+1,n3-kc3-1
      do 40 i2=0,n2-1
      do 40 i1=0,idh-1
         ca(i1,i2,i3) = (0.0d0,0.0d0)
 40   continue
c40   end do
      end if

      return
      end 
c---------------------------------------------
      subroutine fsb_lp235(ca,cb,cw,n,lp,key)
c---------------------------------------------
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)
      dimension ca(0:n-1),cb(0:n-1)
      dimension cw(0:n-1)
      dimension lp(8)

      k =1
      la=1

c*** radix 2 fft.

      do 2 i=1,lp(2)
      if(k.eq.1) then
         call fsb_lp2(ca,cb,cw,la,n,key)
      else
         call fsb_lp2(cb,ca,cw,la,n,key)
      end if
         k=k*(-1)
         la=la*2
    2 continue

c*** radix 3 fft.

      do 3 i=1,lp(3)
      if(k.eq.1) then
         call fsb_lp3(ca,cb,cw,la,n,key)
      else
         call fsb_lp3(cb,ca,cw,la,n,key)
      end if
         k=k*(-1)
         la=la*3
    3 continue

c*** radix 4 fft.

      do 4 i=1,lp(4)
      if(k.eq.1) then
         call fsb_lp4(ca,cb,cw,la,n,key)
      else
         call fsb_lp4(cb,ca,cw,la,n,key)
      end if
         k=k*(-1)
         la=la*4
    4 continue

c*** radix 5 fft.

      do 5 i=1,lp(5)
      if(k.eq.1) then
         call fsb_lp5(ca,cb,cw,la,n,key)
      else
         call fsb_lp5(cb,ca,cw,la,n,key)
      end if
         k=k*(-1)
         la=la*5
    5 continue

c*** radix 6 fft.

      do 6 i=1,lp(6)
      if(k.eq.1) then
         call fsb_lp6(ca,cb,cw,la,n,key)
      else
         call fsb_lp6(cb,ca,cw,la,n,key)
      end if
         k=k*(-1)
         la=la*6
    6 continue

c*** radix 8 fft.

      do 8 i=1,lp(8)
      if(k.eq.1) then
         call fsb_lp8(ca,cb,cw,la,n,key)
      else
         call fsb_lp8(cb,ca,cw,la,n,key)
      end if
         k=k*(-1)
         la=la*8
    8 continue
c***
      if(k.eq.-1) then
         do 10 i=0,n-1
            ca(i)=cb(i)
   10    continue
      end if
c***
      return
      end
c---------------------------------------------
      subroutine fsb_lp2(ca,cb,cw,la,n,key)
c---------------------------------------------
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)
      dimension ca(0:n-1),cb(0:n-1)
      dimension cw(0:n-1)

      m=n/2

      if(la.eq.1) go to 22

      do 20 l=0,la-1

         cb(l   )=(ca(l)+ca(l+m))
         cb(l+la)=(ca(l)-ca(l+m))

   20 continue

      do 21 kpla=la,m-la,la
      do 21 l=kpla,kpla+la-1

         j0=l+kpla         

         cb(j0   )=(ca(l)+ca(l+m))
         cb(j0+la)=(ca(l)-ca(l+m))*cw(kpla)

   21 continue

      return

   22 cb(0)=(ca(0)+ca(m))
      cb(1)=(ca(0)-ca(m))

      do 23 l=1,m-1

         j0=l*2         

         cb(j0  )=(ca(l)+ca(l+m))
         cb(j0+1)=(ca(l)-ca(l+m))*cw(l)

   23 continue

      return
      end
c---------------------------------------------
      subroutine fsb_lp3(ca,cb,cw,la,n,key)
c---------------------------------------------
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)
      dimension ca(0:n-1),cb(0:n-1)
      dimension cw(0:n-1)

      m=n/3

      sin60=0.5d0*dsqrt(3.0d0)
      if(key.gt.0) then
         cisin60=dcmplx(0.0d0,+sin60)
      else
         cisin60=dcmplx(0.0d0,-sin60)
      end if

      if(la.eq.1) go to 32

      do 30 l=0,la-1

         ct1=ca(l+m)+ca(l+m*2)
         ct2=ca(l  )-0.5d0*ct1       
         ct3=cisin60*(ca(l+m)-ca(l+m*2))

         cb(l     )=ca(l)+ct1
         cb(l+la  )=(ct2+ct3)
         cb(l+la*2)=(ct2-ct3)

   30 continue

      do 31 kpla=la,m-la,la
      do 31 l=kpla,kpla+la-1

         j0=l+kpla*2         

         ct1=ca(l+m)+ca(l+m*2)
         ct2=ca(l  )-0.5d0*ct1       
         ct3=cisin60*(ca(l+m)-ca(l+m*2))

         cb(j0     )=ca(l)+ct1
         cb(j0+la  )=(ct2+ct3)*cw(kpla  )
         cb(j0+la*2)=(ct2-ct3)*cw(kpla*2)

   31 continue

      return

   32 ct1=ca(m)+ca(m*2)
      ct2=ca(0)-0.5d0*ct1       
      ct3=cisin60*(ca(m)-ca(m*2))

      cb(0)=ca(0)+ct1
      cb(1)=(ct2+ct3)
      cb(2)=(ct2-ct3)

      do 33 l=1,m-1

         j0=l*3         

         ct1=ca(l+m)+ca(l+m*2)
         ct2=ca(l  )-0.5d0*ct1       
         ct3=cisin60*(ca(l+m)-ca(l+m*2))

         cb(j0  )=ca(l)+ct1
         cb(j0+1)=(ct2+ct3)*cw(l  )
         cb(j0+2)=(ct2-ct3)*cw(l*2)

   33 continue

      return
      end
c---------------------------------------------
      subroutine fsb_lp4(ca,cb,cw,la,n,key)
c---------------------------------------------
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)
      dimension ca(0:n-1),cb(0:n-1)
      dimension cw(0:n-1)

      m=n/4

      if(key.gt.0) then
         ci=(0.0d0,+1.0d0)
      else
         ci=(0.0d0,-1.0d0)
      end if
   
      if(la.eq.1) go to 42

      do 40 l=0,la-1

         cc0 = ca(l  )+ca(l+m*2)
         cc1 = ca(l+m)+ca(l+m*3)
         cc2 = ca(l  )-ca(l+m*2)
         cc3 =(ca(l+m)-ca(l+m*3))*ci

         cb(l     )= cc0+cc1
         cb(l+la  )= cc2+cc3
         cb(l+la*2)= cc0-cc1
         cb(l+la*3)= cc2-cc3

   40 continue
   
      do 41 kpla=la,m-la,la
      do 41 l=kpla,kpla+la-1

         j0=l+kpla*3

         cc1 = ca(l  )+ca(l+m*2)
         cc2 = ca(l+m)+ca(l+m*3)
         cc3 = ca(l  )-ca(l+m*2)
         cc4 =(ca(l+m)-ca(l+m*3))*ci

         cb(j0     )= cc1+cc2
         cb(j0+la  )=(cc3+cc4)*cw(kpla  )
         cb(j0+la*2)=(cc1-cc2)*cw(kpla*2)
         cb(j0+la*3)=(cc3-cc4)*cw(kpla*3)

   41 continue
         
      return

   42 cc0 = ca(0)+ca(m*2)
      cc1 = ca(m)+ca(m*3)
      cc2 = ca(0)-ca(m*2)
      cc3 =(ca(m)-ca(m*3))*ci

      cb(0)= cc0+cc1
      cb(1)= cc2+cc3
      cb(2)= cc0-cc1
      cb(3)= cc2-cc3

      do 43 l=1,m-1

         j0=l*4

         cc1 = ca(l  )+ca(l+m*2)
         cc2 = ca(l+m)+ca(l+m*3)
         cc3 = ca(l  )-ca(l+m*2)
         cc4 =(ca(l+m)-ca(l+m*3))*ci

         cb(j0  )= cc1+cc2
         cb(j0+1)=(cc3+cc4)*cw(l  )
         cb(j0+2)=(cc1-cc2)*cw(l*2)
         cb(j0+3)=(cc3-cc4)*cw(l*3)

   43 continue
         
      return
      end
c---------------------------------------------
      subroutine fsb_lp5(ca,cb,cw,la,n,key)
c---------------------------------------------
      use m_Const_Parameters, only : PAI
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)
      dimension ca(0:n-1),cb(0:n-1)
      dimension cw(0:n-1)

      m=n/5

      r54=dsqrt(5.0d0)/4.0d0
c**      pai=3.14159265358979323846d0
      sin36=dsin(pai*0.2d0)
      sin72=dsin(pai*0.4d0)
      if(key.gt.0) then
         cisin36=dcmplx(0.0d0,+sin36)
         cisin72=dcmplx(0.0d0,+sin72)
      else
         cisin36=dcmplx(0.0d0,-sin36)
         cisin72=dcmplx(0.0d0,-sin72)
      end if

      if(la.eq.1) go to 52

      do 50 l=0,la-1

         ct1=ca(l+m  )+ca(l+m*4)
         ct2=ca(l+m*2)+ca(l+m*3)
         ct3=ca(l+m  )-ca(l+m*4)
         ct4=ca(l+m*2)-ca(l+m*3)
         ct5=     ct1+ct2
         ct6=r54*(ct1-ct2)
         ct7=ca(l   )-0.25d0*ct5
         ct8=ct7+ct6
         ct9=ct7-ct6
         ct10=cisin72*ct3+cisin36*ct4
         ct11=cisin36*ct3-cisin72*ct4

         cb(l     )=ca(l)+ct5
         cb(l+la  )=(ct8+ct10)
         cb(l+la*2)=(ct9+ct11)
         cb(l+la*3)=(ct9-ct11)
         cb(l+la*4)=(ct8-ct10)

   50 continue

      do 51 kpla=la,m-la,la
      do 51 l=kpla,kpla+la-1

         j0=l+kpla*4

         ct1=ca(l+m  )+ca(l+m*4)
         ct2=ca(l+m*2)+ca(l+m*3)
         ct3=ca(l+m  )-ca(l+m*4)
         ct4=ca(l+m*2)-ca(l+m*3)
         ct5=     ct1+ct2
         ct6=r54*(ct1-ct2)
         ct7=ca(l  )-0.25d0*ct5
         ct8=ct7+ct6
         ct9=ct7-ct6
         ct10=cisin72*ct3+cisin36*ct4
         ct11=cisin36*ct3-cisin72*ct4

         cb(j0     )=ca(l)+ct5
         cb(j0+la  )=(ct8+ct10)*cw(kpla  )
         cb(j0+la*2)=(ct9+ct11)*cw(kpla*2)
         cb(j0+la*3)=(ct9-ct11)*cw(kpla*3)
         cb(j0+la*4)=(ct8-ct10)*cw(kpla*4)

   51 continue

      return

   52 ct1=ca(m  )+ca(m*4)
      ct2=ca(m*2)+ca(m*3)
      ct3=ca(m  )-ca(m*4)
      ct4=ca(m*2)-ca(m*3)
      ct5=     ct1+ct2
      ct6=r54*(ct1-ct2)
      ct7=ca(0)-0.25d0*ct5
      ct8=ct7+ct6
      ct9=ct7-ct6
      ct10=cisin72*ct3+cisin36*ct4
      ct11=cisin36*ct3-cisin72*ct4

      cb(0)=ca(0)+ct5
      cb(1)=(ct8+ct10)
      cb(2)=(ct9+ct11)
      cb(3)=(ct9-ct11)
      cb(4)=(ct8-ct10)

      do 53 l=1,m-1

         j0=l*5

         ct1=ca(l+m  )+ca(l+m*4)
         ct2=ca(l+m*2)+ca(l+m*3)
         ct3=ca(l+m  )-ca(l+m*4)
         ct4=ca(l+m*2)-ca(l+m*3)
         ct5=     ct1+ct2
         ct6=r54*(ct1-ct2)
         ct7=ca(l  )-0.25d0*ct5
         ct8=ct7+ct6
         ct9=ct7-ct6
         ct10=cisin72*ct3+cisin36*ct4
         ct11=cisin36*ct3-cisin72*ct4

         cb(j0  )=ca(l)+ct5
         cb(j0+1)=(ct8+ct10)*cw(l  )
         cb(j0+2)=(ct9+ct11)*cw(l*2)
         cb(j0+3)=(ct9-ct11)*cw(l*3)
         cb(j0+4)=(ct8-ct10)*cw(l*4)

   53 continue

      return
      end
c---------------------------------------------
      subroutine fsb_lp6(ca,cb,cw,la,n,key)
c---------------------------------------------
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)
      dimension ca(0:n-1),cb(0:n-1)
      dimension cw(0:n-1)

      m=n/6

      sin60=0.5d0*dsqrt(3.0d0)
      if(key.gt.0) then
         cisin60=dcmplx(0.0d0,+sin60)
      else
         cisin60=dcmplx(0.0d0,-sin60)
      end if

      if(la.eq.1) go to 62

      do 60 l=0,la-1

         ct0=ca(l+m*2)+ca(l+m*4)
         ct1=ca(l+m  )+ca(l+m*5)
         ct2=ca(l    )-0.5d0*ct0
         ct3=ca(l+m*3)-0.5d0*ct1
         ct4=cisin60*(ca(l+m*2)-ca(l+m*4))
         ct5=cisin60*(ca(l+m*5)-ca(l+m  ))

         cd0 = ca(l    )+ct0
         cd1 = ct3-ct5
         cd2 = ct2+ct4
         cd3 = ca(l+m*3)+ct1
         cd4 = ct2-ct4
         cd5 = ct3+ct5

         cb(l     )=(cd0+cd3)
         cb(l+la  )=(cd2-cd5)
         cb(l+la*2)=(cd4+cd1)
         cb(l+la*3)=(cd0-cd3)
         cb(l+la*4)=(cd2+cd5)
         cb(l+la*5)=(cd4-cd1)

   60 continue

      do 61 kpla=la,m-la,la
      do 61 l=kpla,kpla+la-1

         j0=l+kpla*5

         ct0=ca(l+m*2)+ca(l+m*4)
         ct1=ca(l+m  )+ca(l+m*5)
         ct2=ca(l    )-0.5d0*ct0
         ct3=ca(l+m*3)-0.5d0*ct1
         ct4=cisin60*(ca(l+m*2)-ca(l+m*4))
         ct5=cisin60*(ca(l+m*5)-ca(l+m  ))

         cd0 = ca(l    )+ct0
         cd1 = ct3-ct5
         cd2 = ct2+ct4
         cd3 = ca(l+m*3)+ct1
         cd4 = ct2-ct4
         cd5 = ct3+ct5

         cb(j0     )=(cd0+cd3)
         cb(j0+la  )=(cd2-cd5)*cw(kpla  )
         cb(j0+la*2)=(cd4+cd1)*cw(kpla*2)
         cb(j0+la*3)=(cd0-cd3)*cw(kpla*3)
         cb(j0+la*4)=(cd2+cd5)*cw(kpla*4)
         cb(j0+la*5)=(cd4-cd1)*cw(kpla*5)

   61 continue

      return

   62 ct0=ca(m*2)+ca(m*4)
      ct1=ca(m  )+ca(m*5)
      ct2=ca(0  )-0.5d0*ct0
      ct3=ca(m*3)-0.5d0*ct1
      ct4=cisin60*(ca(m*2)-ca(m*4))
      ct5=cisin60*(ca(m*5)-ca(m  ))

      cd0=ca(0  )+ct0
      cd1=ct3-ct5
      cd2=ct2+ct4
      cd3=ca(m*3)+ct1
      cd4=ct2-ct4
      cd5=ct3+ct5

      cb(0)=(cd0+cd3)
      cb(1)=(cd2-cd5)
      cb(2)=(cd4+cd1)
      cb(3)=(cd0-cd3)
      cb(4)=(cd2+cd5)
      cb(5)=(cd4-cd1)

      do 63 l=1,m-1

         j0=l*6

         ct0=ca(l+m*2)+ca(l+m*4)
         ct1=ca(l+m  )+ca(l+m*5)
         ct2=ca(l    )-0.5d0*ct0
         ct3=ca(l+m*3)-0.5d0*ct1
         ct4=cisin60*(ca(l+m*2)-ca(l+m*4))
         ct5=cisin60*(ca(l+m*5)-ca(l+m  ))

         cd0=ca(l    )+ct0
         cd1=ct3-ct5
         cd2=ct2+ct4
         cd3=ca(l+m*3)+ct1
         cd4=ct2-ct4
         cd5=ct3+ct5

         cb(j0  )=(cd0+cd3)
         cb(j0+1)=(cd2-cd5)*cw(l  )
         cb(j0+2)=(cd4+cd1)*cw(l*2)
         cb(j0+3)=(cd0-cd3)*cw(l*3)
         cb(j0+4)=(cd2+cd5)*cw(l*4)
         cb(j0+5)=(cd4-cd1)*cw(l*5)

   63 continue

      return
      end
c---------------------------------------------
      subroutine fsb_lp8(ca,cb,cw,la,n,key)
c---------------------------------------------
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)
      dimension ca(0:n-1),cb(0:n-1)
      dimension cw(0:n-1)
   
      m=n/8

      x=dsqrt(0.5d0)
      if(key.gt.0) then
         ci1=dcmplx(x, +x)
      else
         ci1=dcmplx(x, -x)
      end if
      ci2=ci1*ci1
      ci3=ci1*ci2

      if(la.eq.1) go to 82

      do 80 l=0,la-1

         cc0 = ca(l    )+ca(l+m*4)
         cc1 = ca(l+m  )+ca(l+m*5)
         cc2 = ca(l+m*2)+ca(l+m*6)
         cc3 = ca(l+m*3)+ca(l+m*7)
         cc4 = ca(l    )-ca(l+m*4)
         cc5 = ca(l+m  )-ca(l+m*5)
         cc6 =(ca(l+m*2)-ca(l+m*6))*ci2
         cc7 =(ca(l+m*3)-ca(l+m*7))*ci2

         cd0 = cc0+cc2
         cd1 = cc1+cc3
         cd2 = cc4+cc6
         cd3 =(cc5+cc7)*ci1
         cd4 = cc0-cc2
         cd5 =(cc1-cc3)*ci2
         cd6 = cc4-cc6
         cd7 =(cc5-cc7)*ci3

         cb(l     ) = cd0+cd1
         cb(l+la  ) = cd2+cd3
         cb(l+la*2) = cd4+cd5
         cb(l+la*3) = cd6+cd7
         cb(l+la*4) = cd0-cd1
         cb(l+la*5) = cd2-cd3
         cb(l+la*6) = cd4-cd5
         cb(l+la*7) = cd6-cd7

   80 continue
   
      do 81 kpla=la,m-la,la
      do 81 l=kpla,kpla+la-1

         j0=l+kpla*7

         cc0 = ca(l    )+ca(l+m*4)
         cc1 = ca(l+m  )+ca(l+m*5)
         cc2 = ca(l+m*2)+ca(l+m*6)
         cc3 = ca(l+m*3)+ca(l+m*7)
         cc4 = ca(l    )-ca(l+m*4)
         cc5 = ca(l+m  )-ca(l+m*5)
         cc6 =(ca(l+m*2)-ca(l+m*6))*ci2
         cc7 =(ca(l+m*3)-ca(l+m*7))*ci2

         cd0 = cc0+cc2
         cd1 = cc1+cc3
         cd2 = cc4+cc6
         cd3 =(cc5+cc7)*ci1
         cd4 = cc0-cc2
         cd5 =(cc1-cc3)*ci2
         cd6 = cc4-cc6
         cd7 =(cc5-cc7)*ci3

         cb(j0     ) = cd0+cd1
         cb(j0+la  ) =(cd2+cd3)*cw(kpla  )
         cb(j0+la*2) =(cd4+cd5)*cw(kpla*2)
         cb(j0+la*3) =(cd6+cd7)*cw(kpla*3)
         cb(j0+la*4) =(cd0-cd1)*cw(kpla*4)
         cb(j0+la*5) =(cd2-cd3)*cw(kpla*5)
         cb(j0+la*6) =(cd4-cd5)*cw(kpla*6)
         cb(j0+la*7) =(cd6-cd7)*cw(kpla*7)

   81 continue
    
      return

   82 cc0 = ca(0  )+ca(m*4)
      cc1 = ca(m  )+ca(m*5)
      cc2 = ca(m*2)+ca(m*6)
      cc3 = ca(m*3)+ca(m*7)
      cc4 = ca(0  )-ca(m*4)
      cc5 = ca(m  )-ca(m*5)
      cc6 =(ca(m*2)-ca(m*6))*ci2
      cc7 =(ca(m*3)-ca(m*7))*ci2

      cd0 = cc0+cc2
      cd1 = cc1+cc3
      cd2 = cc4+cc6
      cd3 =(cc5+cc7)*ci1
      cd4 = cc0-cc2
      cd5 =(cc1-cc3)*ci2
      cd6 = cc4-cc6
      cd7 =(cc5-cc7)*ci3

      cb(0) = cd0+cd1
      cb(1) = cd2+cd3
      cb(2) = cd4+cd5
      cb(3) = cd6+cd7
      cb(4) = cd0-cd1
      cb(5) = cd2-cd3
      cb(6) = cd4-cd5
      cb(7) = cd6-cd7
   
      do 83 l=1,m-1

         j0=l*8

         cc0 = ca(l    )+ca(l+m*4)
         cc1 = ca(l+m  )+ca(l+m*5)
         cc2 = ca(l+m*2)+ca(l+m*6)
         cc3 = ca(l+m*3)+ca(l+m*7)
         cc4 = ca(l    )-ca(l+m*4)
         cc5 = ca(l+m  )-ca(l+m*5)
         cc6 =(ca(l+m*2)-ca(l+m*6))*ci2
         cc7 =(ca(l+m*3)-ca(l+m*7))*ci2

         cd0 = cc0+cc2
         cd1 = cc1+cc3
         cd2 = cc4+cc6
         cd3 =(cc5+cc7)*ci1
         cd4 = cc0-cc2
         cd5 =(cc1-cc3)*ci2
         cd6 = cc4-cc6
         cd7 =(cc5-cc7)*ci3

         cb(j0  ) = cd0+cd1
         cb(j0+1) =(cd2+cd3)*cw(l  )
         cb(j0+2) =(cd4+cd5)*cw(l*2)
         cb(j0+3) =(cd6+cd7)*cw(l*3)
         cb(j0+4) =(cd0-cd1)*cw(l*4)
         cb(j0+5) =(cd2-cd3)*cw(l*5)
         cb(j0+6) =(cd4-cd5)*cw(l*6)
         cb(j0+7) =(cd6-cd7)*cw(l*7)

   83 continue

      return
      end
c---------------------------------------------
      subroutine fsb_lp235v(ca,cb,cw,n,ndim,nn1,nn2,lp,key)
c---------------------------------------------
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)
      dimension ca(ndim,0:n-1),cb(0:n-1,2)
      dimension cw(0:n-1)
      dimension lp(8)

      do 10 i=nn1,nn2

         do 20 j=0,n-1
            cb(j,1) = ca(i,j)
  20     continue
 
         call fsb_lp235(cb(0,1),cb(0,2),cw,n,lp,key)

         do 30 j=0,n-1
            ca(i,j) = cb(j,1)
  30     continue

  10  continue
  
c***
      return
      end
c---------------------------------------------
      subroutine setcw(cw,n)
c---------------------------------------------
      use m_Const_Parameters, only : PAI
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)
      complex*16 cw(0:n-1)

c**      pai=3.14159265358979323846d0
      pc =2.0d0*pai/dble(n)

      do 10 i=0,n-1
         tr=dcos(pc*dble(i))
         ti=dsin(pc*dble(i))
         cw(i)=dcmplx(tr,ti)
   10 continue

      return
      end
c---------------------------------------------
      subroutine setlp(n,lp)
c---------------------------------------------
      implicit real*8(a-b,d-h,o-z)
      integer lp(8)

         nn =n

         lp(8)=0
         lp(6)=0
         lp(5)=0
         lp(4)=0
         lp(3)=0
         lp(2)=0

      do while(mod(nn,8).eq.0)
         lp(8)=lp(8)+1
         nn=nn/8
      end do

c     do while(mod(nn,6).eq.0)
c        lp(6)=lp(6)+1
c        nn=nn/6
c     end do

      do while(mod(nn,5).eq.0)
         lp(5)=lp(5)+1
         nn=nn/5
      end do

      do while(mod(nn,4).eq.0)
         lp(4)=lp(4)+1
         nn=nn/4
      end do

      do while(mod(nn,3).eq.0)
         lp(3)=lp(3)+1
         nn=nn/3
      end do

      do while(mod(nn,2).eq.0)
         lp(2)=lp(2)+1
         nn=nn/2
      end do

      if(nn.ne.1) stop 'n is an improper number.'

      return
      end
c---------------------------------------------
      subroutine trans_b(ca,cw1,n1h,kc1)
c---------------------------------------------
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)
      dimension ca (0:n1h)
      dimension cw1(0:n1h-1)
c**
      r0 = dble(ca(  0))
      r1 = dble(ca(n1h))
       
      ca(  0)=dcmplx( r0 + r1 , r0 - r1 )
      ca(n1h)=(0.0d0,0.0d0)
c**   
      if(kc1.eq.0) then
         nn=n1h/2
      else
         nn=min(n1h/2,kc1) 
      end if

c**COMMENT.
c* When n1h is an even number and k1 is n1h/2, k1 is eqaul to n1h-k1.
c* This does not make any problem.
c**
      do 10 k1 = 1,nn

      cc1 = ca(k1)
      cc2 = ca(n1h-k1)
      cc1j = dconjg(cc1)
      cc2j = dconjg(cc2)

      ca(k1    ) = cc1+cc2j + (cc1-cc2j)*cw1(k1    )
      ca(n1h-k1) = cc2+cc1j + (cc2-cc1j)*cw1(n1h-k1)

  10  continue
c 10  end do
c**   
      return
      end

c---------------------------------------------
      subroutine trans_f(ca,cw1,n1h,kc1)
c---------------------------------------------
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)
      dimension ca (0:n1h)
      dimension cw1(0:n1h-1)
c**
      rr = dble (ca(0))
      ri = dimag(ca(0))
       
      ca(0  )=dcmplx(rr+ri,0.0d0)
      ca(n1h)=dcmplx(rr-ri,0.0d0)
c**   
      if(kc1.eq.0) then
         nn=n1h/2
      else
         nn=min(n1h/2,kc1) 
      end if

c**COMMENT.
c* When n1h is an even number and k1 is n1h/2, k1 is eqaul to n1h-k1.
c* This does not make any problem.
c**
      do 10 k1 = 1,nn

      cc1 = ca(k1)
      cc2 = ca(n1h-k1)
      cc1j = dconjg(cc1)
      cc2j = dconjg(cc2)

      ca(k1    ) = ( cc1+cc2j + (cc1-cc2j)*cw1(k1    ) )*0.5d0
      ca(n1h-k1) = ( cc2+cc1j + (cc2-cc1j)*cw1(n1h-k1) )*0.5d0

  10  continue
c 10  end do
c**   
      return
      end
#endif
