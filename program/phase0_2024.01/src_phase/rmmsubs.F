!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: rdecomp, rsolve, cdecomp, csolve, csolve2,  rg_eispack, 
!           balanc, balbak, cdiv, elmhes, eltran, hqr, hqr2, cg_eispack,
!           cbabk2, cbal, comqr, comqr2, corth, csroot, mltpha4, 
!          decide_alpha, rmm3_uda, rmm2_uda,  crmm2_uda, crmm3_uda,
!
!  AUTHOR(S): T. Uda, T. Yamasaki   August/20/2003
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
!  $Id: rmmsubs.F 570 2017-04-21 20:34:50Z yamasaki $
c---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine rdecomp(n,a,w,ip,ier)
c     $Id: rmmsubs.F 570 2017-04-21 20:34:50Z yamasaki $
c                       rmmsubs.f 1.2 99/06/04 21:17:36
c                       rmmsubs.f 1.3 03/02/20 18:35:31
c     rmm routine is available for systems
c                        with inversion symmetry or
c                        GAMMA point sampling  
c
c***********************************************
c*    LU decomposition of real matrix a(n,n)
c***********************************************
      real*8 a(n,n),w(n)
      real*8 al
      real*8 eps
      integer ip(n)
      integer ier
c     parameter(eps=1.0d-75)
      parameter(eps=1.0d-30)
#ifdef __TIMER_SUB__
      call timer_sta(1135)
#endif

      do k=1,n
         ip(k)=k
      end do

      do k=1,n
         l=k
         al=abs(a(ip(l),k))
#ifdef __TIMER_DO__
      call timer_sta(1182)
#endif
	      do i=k+1,n
		 if(abs(a(ip(i),k)).gt.al) then
		    l=i
		    al=abs(a(ip(l),k))
		 end if
	      end do
#ifdef __TIMER_DO__
      call timer_end(1182)
#endif

      if(l.ne.k) then
         lv=ip(k)
         ip(k)=ip(l)
         ip(l)=lv
      end if

      if(dabs(a(ip(k),k)).lt.eps) then ! This matrix is singular.
c        write(6,*) 'dabs(a(ip(k),k))',dabs(a(ip(k),k))
         ier = 1 
         return
      end if

c** Gaussian elimination **
      a(ip(k),k)=1.0d0/a(ip(k),k)

#ifdef __TIMER_DO__
      call timer_sta(1183)
#endif
	      do i=k+1,n

		 a(ip(i),k)=a(ip(i),k)*a(ip(k),k)

		      do j=k+1,n
			 w(j)=a(ip(i),j)-a(ip(i),k)*a(ip(k),j)
		      end do

		      do j=k+1,n
			 a(ip(i),j)=w(j)
		      end do
	      end do
#ifdef __TIMER_DO__
      call timer_end(1183)
#endif
      end do
#ifdef __TIMER_SUB__
      call timer_end(1135)
#endif
         ier = 0 ! normal return
      return
      end
c---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine rsolve(n,m,a,b,x,ip)
c                           @(#)rmmsubs.f 1.2 99/06/04 21:17:36
c***************************************
c*    Solve the linear equation.
c*    matrix a(n,n) must be LU-decomposed
c*    before the calling of this routine.
c***************************************

      real*8 a(n,n),b(n,m),x(n,m)
      integer ip(n)

      real*8 t
#ifdef __TIMER_SUB__
      call timer_sta(1136)
#endif

#ifdef __TIMER_DO__
      call timer_sta(1184)
#endif
      do k=1,m
         do i=1,n
            t=b(ip(i),k)
            do j=1,i-1
               t=t-a(ip(i),j)*x(j,k)
            end do
            x(i,k)=t
         end do
      end do
#ifdef __TIMER_DO__
      call timer_end(1184)
#endif

#ifdef __TIMER_DO__
      call timer_sta(1185)
#endif
      do k=1,m
         do i=n,1,-1
            t=x(i,k)
            do j=i+1,n
               t=t-a(ip(i),j)*x(j,k)
            end do
            x(i,k)=t*a(ip(i),i)
         end do
      end do
#ifdef __TIMER_DO__
      call timer_end(1185)
#endif

#ifdef __TIMER_SUB__
      call timer_end(1136)
#endif
      return
      end
c---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine cdecomp(n,a,w,ip,ier)
c                           @(#)rmmsubs.f 1.2 99/06/04 21:17:36
c***********************************************
c*    LU decomposition of complex matrix a(n,n)
c***********************************************
      complex*16 a(n,n),w(n)
      real*8 al
      real*8 eps
      integer ip(n)
      integer ier 
      parameter(eps=1.0d-10)

      do k=1,n
         ip(k)=k
      end do

      do k=1,n
         l=k
         al=cdabs(a(ip(l),k))
	      do i=k+1,n
		 if(cdabs(a(ip(i),k)).gt.al) then
		    l=i
		    al=cdabs(a(ip(l),k))
		 end if
	      end do

      if(l.ne.k) then
         lv=ip(k)
         ip(k)=ip(l)
         ip(l)=lv
      end if

      if(cdabs(a(ip(k),k)).lt.eps) then ! This matrix is singular.
         ier = 1 
         return
      end if

c** Gaussian elimination **
      a(ip(k),k)=1.0d0/a(ip(k),k)

	      do i=k+1,n

		 a(ip(i),k)=a(ip(i),k)*a(ip(k),k)

		      do j=k+1,n
			 w(j)=a(ip(i),j)-a(ip(i),k)*a(ip(k),j)
		      end do

		      do j=k+1,n
			 a(ip(i),j)=w(j)
		      end do
	      end do
      end do
         ier = 0 ! normal return
      return
      end
c---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine csolve(n,m,a,b,x,ip)
c                           @(#)rmmsubs.f 1.2 99/06/04 21:17:36
c******************************************************
c*    Solve the linear equation.
c*    The complex matrix a(n,n) must be LU-decomposed
c*    before calling this routine.
c******************************************************
      complex*16 a(n,n),b(n,m),x(n,m)
      integer ip(n)

      complex*16 t

      do k=1,m
         do i=1,n
            t=b(ip(i),k)
            do j=1,i-1
               t=t-a(ip(i),j)*x(j,k)
            end do
            x(i,k)=t
         end do
      end do

      do k=1,m
         do i=n,1,-1
            t=x(i,k)
            do j=i+1,n
               t=t-a(ip(i),j)*x(j,k)
            end do
            x(i,k)=t*a(ip(i),i)
         end do
      end do

      return
      end

c---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine csolve2(n,m,a,b,x,ip)
c                           @(#)rmmsubs.f 1.2 99/06/04 21:17:36
c******************************************************
c*    Solve the linear equation.
c*    The complex matrix a(n,n) must be LU-decomposed
c*    before calling this routine.
c******************************************************
      real*8 a(n,n,2),b(n,m,2),x(n,m,2)
      integer ip(n)

      complex*16 t, q

      do k=1,m
         do i=1,n
            t=dcmplx(b(ip(i),k,1),b(ip(i),k,2))
            do j=1,i-1
               t=t-dcmplx(a(ip(i),j,1),a(ip(i),j,2)) 
     &             *  dcmplx(x(j,k,1),x(j,k,2))
            end do
            x(i,k,1) = dreal(t)
            x(i,k,2) = dimag(t)
         end do
      end do

      do k=1,m
         do i=n,1,-1
            t=dcmplx(x(i,k,1),x(i,k,2))
            do j=i+1,n
               t=t - dcmplx(a(ip(i),j,1),a(ip(i),j,2))
     &              *dcmplx(x(j,k,1),    x(j,k,2))
            end do
            q = t*dcmplx(a(ip(i),i,1),a(ip(i),i,2))
            x(i,k,1) = dreal(q)
            x(i,k,2) = dimag(q)
         end do
      end do

      return
      end
c---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine rg_eispack(nm,n,a,wr,wi,matz,z,iv1,fv1,ierr)
      integer n,nm,is1,is2,ierr,matz
      double precision a(nm,n),wr(n),wi(n),z(nm,n),fv1(n)
      integer iv1(n)
c
c     this subroutine calls the recommended sequence of
c     subroutines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     of a real general matrix.
c
c     on input
c
c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.
c
c        n  is the order of the matrix  a.
c
c        a  contains the real general matrix.
c
c        matz  is an integer variable set equal to zero if
c        only eigenvalues are desired.  otherwise it is set to
c        any non-zero integer for both eigenvalues and eigenvectors.
c
c     on output
c
c        wr  and  wi  contain the real and imaginary parts,
c        respectively, of the eigenvalues.  complex conjugate
c        pairs of eigenvalues appear consecutively with the
c        eigenvalue having the positive imaginary part first.
c
c        z  contains the real and imaginary parts of the eigenvectors
c        if matz is not zero.  if the j-th eigenvalue is real, the
c        j-th column of  z  contains its eigenvector.  if the j-th
c        eigenvalue is complex with positive imaginary part, the
c        j-th and (j+1)-th columns of  z  contain the real and
c        imaginary parts of its eigenvector.  the conjugate of this
c        vector is the eigenvector for the conjugate eigenvalue.
c
c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for hqr
c           and hqr2.  the normal completion code is zero.
c
c        iv1  and  fv1  are temporary storage arrays.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
c
   10 call  balanc(nm,n,a,is1,is2,fv1)
      call  elmhes(nm,n,is1,is2,a,iv1)
      if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      call  hqr(nm,n,is1,is2,a,wr,wi,ierr)
      go to 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 call  eltran(nm,n,is1,is2,a,iv1,z)
      call  hqr2(nm,n,is1,is2,a,wr,wi,z,ierr)
      if (ierr .ne. 0) go to 50
      call  balbak(nm,n,is1,is2,fv1,n,z)
   50 return
      end
      subroutine balanc(nm,n,a,low,igh,scale)
c
      integer i,j,k,l,m,n,jj,nm,igh,low,iexc
      double precision a(nm,n),scale(n)
      double precision c,f,g,r,s,b2,radix
      logical noconv
c
c     this subroutine is a translation of the algol procedure balance,
c     num. math. 13, 293-304(1969) by parlett and reinsch.
c     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
c
c     this subroutine balances a real matrix and isolates
c     eigenvalues whenever possible.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        a contains the input matrix to be balanced.
c
c     on output
c
c        a contains the balanced matrix.
c
c        low and igh are two integers such that a(i,j)
c          is equal to zero if
c           (1) i is greater than j and
c           (2) j=1,...,low-1 or i=igh+1,...,n.
c
c        scale contains information determining the
c           permutations and scaling factors used.
c
c     suppose that the principal submatrix in rows low through igh
c     has been balanced, that p(j) denotes the index interchanged
c     with j during the permutation step, and that the elements
c     of the diagonal matrix used are denoted by d(i,j).  then
c        scale(j) = p(j),    for j = 1,...,low-1
c                 = d(j,j),      j = low,...,igh
c                 = p(j)         j = igh+1,...,n.
c     the order in which the interchanges are made is n to igh+1,
c     then 1 to low-1.
c
c     note that 1 is returned for igh if igh is zero formally.
c
c     the algol procedure exc contained in balance appears in
c     balanc  in line.  (note that the algol roles of identifiers
c     k,l have been reversed.)
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      radix = 16.0d0
c
      b2 = radix * radix
      k = 1
      l = n
      go to 100
c     .......... in-line procedure for row and
c                column exchange ..........
   20 scale(m) = j
      if (j .eq. m) go to 50
c
      do 30 i = 1, l
         f = a(i,j)
         a(i,j) = a(i,m)
         a(i,m) = f
   30 continue
c
      do 40 i = k, n
         f = a(j,i)
         a(j,i) = a(m,i)
         a(m,i) = f
   40 continue
c
   50 go to (80,130), iexc
c     .......... search for rows isolating an eigenvalue
c                and push them down ..........
   80 if (l .eq. 1) go to 280
      l = l - 1
c     .......... for j=l step -1 until 1 do -- ..........
  100 do 120 jj = 1, l
         j = l + 1 - jj
c
         do 110 i = 1, l
            if (i .eq. j) go to 110
            if (a(j,i) .ne. 0.0d0) go to 120
  110    continue
c
         m = l
         iexc = 1
         go to 20
  120 continue
c
      go to 140
c     .......... search for columns isolating an eigenvalue
c                and push them left ..........
  130 k = k + 1
c
  140 do 170 j = k, l
c
         do 150 i = k, l
            if (i .eq. j) go to 150
            if (a(i,j) .ne. 0.0d0) go to 170
  150    continue
c
         m = k
         iexc = 2
         go to 20
  170 continue
c     .......... now balance the submatrix in rows k to l ..........
      do 180 i = k, l
  180 scale(i) = 1.0d0
c     .......... iterative loop for norm reduction ..........
  190 noconv = .false.
c
      do 270 i = k, l
         c = 0.0d0
         r = 0.0d0
c
         do 200 j = k, l
            if (j .eq. i) go to 200
            c = c + dabs(a(j,i))
            r = r + dabs(a(i,j))
  200    continue
c     .......... guard against zero c or r due to underflow ..........
         if (c .eq. 0.0d0 .or. r .eq. 0.0d0) go to 270
         g = r / radix
         f = 1.0d0
         s = c + r
  210    if (c .ge. g) go to 220
         f = f * radix
         c = c * b2
         go to 210
  220    g = r * radix
  230    if (c .lt. g) go to 240
         f = f / radix
         c = c / b2
         go to 230
c     .......... now balance ..........
  240    if ((c + r) / f .ge. 0.95d0 * s) go to 270
         g = 1.0d0 / f
         scale(i) = scale(i) * f
         noconv = .true.
c
         do 250 j = k, n
  250    a(i,j) = a(i,j) * g
c
         do 260 j = 1, l
  260    a(j,i) = a(j,i) * f
c
  270 continue
c
      if (noconv) go to 190
c
  280 low = k
      igh = l
      return
      end
      subroutine balbak(nm,n,low,igh,scale,m,z)
c
      integer i,j,k,m,n,ii,nm,igh,low
      double precision scale(n),z(nm,m)
      double precision s
c
c     this subroutine is a translation of the algol procedure balbak,
c     num. math. 13, 293-304(1969) by parlett and reinsch.
c     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
c
c     this subroutine forms the eigenvectors of a real general
c     matrix by back transforming those of the corresponding
c     balanced matrix determined by  balanc.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by  balanc.
c
c        scale contains information determining the permutations
c          and scaling factors used by  balanc.
c
c        m is the number of columns of z to be back transformed.
c
c        z contains the real and imaginary parts of the eigen-
c          vectors to be back transformed in its first m columns.
c
c     on output
c
c        z contains the real and imaginary parts of the
c          transformed eigenvectors in its first m columns.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      if (m .eq. 0) go to 200
      if (igh .eq. low) go to 120
c
      do 110 i = low, igh
         s = scale(i)
c     .......... left hand eigenvectors are back transformed
c                if the foregoing statement is replaced by
c                s=1.0d0/scale(i). ..........
         do 100 j = 1, m
  100    z(i,j) = z(i,j) * s
c
  110 continue
c     ......... for i=low-1 step -1 until 1,
c               igh+1 step 1 until n do -- ..........
  120 do 140 ii = 1, n
         i = ii
         if (i .ge. low .and. i .le. igh) go to 140
         if (i .lt. low) i = low - ii
         k = scale(i)
         if (k .eq. i) go to 140
c
         do 130 j = 1, m
            s = z(i,j)
            z(i,j) = z(k,j)
            z(k,j) = s
  130    continue
c
  140 continue
c
  200 return
      end
      subroutine cdiv(ar,ai,br,bi,cr,ci)
      double precision ar,ai,br,bi,cr,ci
c
c     complex division, (cr,ci) = (ar,ai)/(br,bi)
c
      double precision s,ars,ais,brs,bis
      s = dabs(br) + dabs(bi)
      ars = ar/s
      ais = ai/s
      brs = br/s
      bis = bi/s
      s = brs**2 + bis**2
      cr = (ars*brs + ais*bis)/s
      ci = (ais*brs - ars*bis)/s
      return
      end
      subroutine elmhes(nm,n,low,igh,a,int)
c
      integer i,j,m,n,la,nm,igh,kp1,low,mm1,mp1
      double precision a(nm,n)
      double precision x,y
      integer int(igh)
c
c     this subroutine is a translation of the algol procedure elmhes,
c     num. math. 12, 349-368(1968) by martin and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
c
c     given a real general matrix, this subroutine
c     reduces a submatrix situated in rows and columns
c     low through igh to upper hessenberg form by
c     stabilized elementary similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by the balancing
c          subroutine  balanc.  if  balanc  has not been used,
c          set low=1, igh=n.
c
c        a contains the input matrix.
c
c     on output
c
c        a contains the hessenberg matrix.  the multipliers
c          which were used in the reduction are stored in the
c          remaining triangle under the hessenberg matrix.
c
c        int contains information on the rows and columns
c          interchanged in the reduction.
c          only elements low through igh are used.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
c
      do 180 m = kp1, la
         mm1 = m - 1
         x = 0.0d0
         i = m
c
         do 100 j = m, igh
            if (dabs(a(j,mm1)) .le. dabs(x)) go to 100
            x = a(j,mm1)
            i = j
  100    continue
c
         int(m) = i
         if (i .eq. m) go to 130
c     .......... interchange rows and columns of a ..........
         do 110 j = mm1, n
            y = a(i,j)
            a(i,j) = a(m,j)
            a(m,j) = y
  110    continue
c
         do 120 j = 1, igh
            y = a(j,i)
            a(j,i) = a(j,m)
            a(j,m) = y
  120    continue
c     .......... end interchange ..........
  130    if (x .eq. 0.0d0) go to 180
         mp1 = m + 1
c
         do 160 i = mp1, igh
            y = a(i,mm1)
            if (y .eq. 0.0d0) go to 160
            y = y / x
            a(i,mm1) = y
c
            do 140 j = m, n
  140       a(i,j) = a(i,j) - y * a(m,j)
c
            do 150 j = 1, igh
  150       a(j,m) = a(j,m) + y * a(j,i)
c
  160    continue
c
  180 continue
c
  200 return
      end
      subroutine eltran(nm,n,low,igh,a,int,z)
c
      integer i,j,n,kl,mm,mp,nm,igh,low,mp1
      double precision a(nm,igh),z(nm,n)
      integer int(igh)
c
c     this subroutine is a translation of the algol procedure elmtrans,
c     num. math. 16, 181-204(1970) by peters and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
c
c     this subroutine accumulates the stabilized elementary
c     similarity transformations used in the reduction of a
c     real general matrix to upper hessenberg form by  elmhes.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by the balancing
c          subroutine  balanc.  if  balanc  has not been used,
c          set low=1, igh=n.
c
c        a contains the multipliers which were used in the
c          reduction by  elmhes  in its lower triangle
c          below the subdiagonal.
c
c        int contains information on the rows and columns
c          interchanged in the reduction by  elmhes.
c          only elements low through igh are used.
c
c     on output
c
c        z contains the transformation matrix produced in the
c          reduction by  elmhes.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
c     .......... initialize z to identity matrix ..........
      do 80 j = 1, n
c
         do 60 i = 1, n
   60    z(i,j) = 0.0d0
c
         z(j,j) = 1.0d0
   80 continue
c
      kl = igh - low - 1
      if (kl .lt. 1) go to 200
c     .......... for mp=igh-1 step -1 until low+1 do -- ..........
      do 140 mm = 1, kl
         mp = igh - mm
         mp1 = mp + 1
c
         do 100 i = mp1, igh
  100    z(i,mp) = a(i,mp-1)
c
         i = int(mp)
         if (i .eq. mp) go to 140
c
         do 130 j = mp, igh
            z(mp,j) = z(i,j)
            z(i,j) = 0.0d0
  130    continue
c
         z(i,mp) = 1.0d0
  140 continue
c
  200 return
      end
      double precision function epslon (x)
      double precision x
c
c     estimate unit roundoff in quantities of size x.
c
      double precision a,b,c,eps
c
c     this program should function properly on all systems
c     satisfying the following two assumptions,
c        1.  the base used in representing floating point
c            numbers is not a power of three.
c        2.  the quantity  a  in statement 10 is represented to 
c            the accuracy used in floating point variables
c            that are stored in memory.
c     the statement number 10 and the go to 10 are intended to
c     force optimizing compilers to generate code satisfying 
c     assumption 2.
c     under these assumptions, it should be true that,
c            a  is not exactly equal to four-thirds,
c            b  has a zero for its last bit or digit,
c            c  is not exactly equal to one,
c            eps  measures the separation of 1.0 from
c                 the next larger floating point number.
c     the developers of eispack would appreciate being informed
c     about any systems where these assumptions do not hold.
c
c     this version dated 4/6/83.
c
      a = 4.0d0/3.0d0
   10 b = a - 1.0d0
      c = b + b + b
      eps = dabs(c-1.0d0)
      if (eps .eq. 0.0d0) go to 10
      epslon = eps*dabs(x)
      return
      end
      subroutine hqr(nm,n,low,igh,h,wr,wi,ierr)
C  RESTORED CORRECT INDICES OF LOOPS (200,210,230,240). (9/29/89 BSG)
c
      integer i,j,k,l,m,n,en,ll,mm,na,nm,igh,itn,its,low,mp2,enm2,ierr
      double precision h(nm,n),wr(n),wi(n)
      double precision p,q,r,s,t,w,x,y,zz,norm,tst1,tst2
      logical notlas
c
c     this subroutine is a translation of the algol procedure hqr,
c     num. math. 14, 219-231(1970) by martin, peters, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 359-371(1971).
c
c     this subroutine finds the eigenvalues of a real
c     upper hessenberg matrix by the qr method.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by the balancing
c          subroutine  balanc.  if  balanc  has not been used,
c          set low=1, igh=n.
c
c        h contains the upper hessenberg matrix.  information about
c          the transformations used in the reduction to hessenberg
c          form by  elmhes  or  orthes, if performed, is stored
c          in the remaining triangle under the hessenberg matrix.
c
c     on output
c
c        h has been destroyed.  therefore, it must be saved
c          before calling  hqr  if subsequent calculation and
c          back transformation of eigenvectors is to be performed.
c
c        wr and wi contain the real and imaginary parts,
c          respectively, of the eigenvalues.  the eigenvalues
c          are unordered except that complex conjugate pairs
c          of values appear consecutively with the eigenvalue
c          having the positive imaginary part first.  if an
c          error exit is made, the eigenvalues should be correct
c          for indices ierr+1,...,n.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the limit of 30*n iterations is exhausted
c                     while the j-th eigenvalue is being sought.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated september 1989.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      norm = 0.0d0
      k = 1
c     .......... store roots isolated by balanc
c                and compute matrix norm ..........
      do 50 i = 1, n
c
         do 40 j = k, n
   40    norm = norm + dabs(h(i,j))
c
         k = i
         if (i .ge. low .and. i .le. igh) go to 50
         wr(i) = h(i,i)
         wi(i) = 0.0d0
   50 continue
c
      en = igh
      t = 0.0d0
      itn = 30*n
c     .......... search for next eigenvalues ..........
   60 if (en .lt. low) go to 1001
      its = 0
      na = en - 1
      enm2 = na - 1
c     .......... look for single small sub-diagonal element
c                for l=en step -1 until low do -- ..........
   70 do 80 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 100
         s = dabs(h(l-1,l-1)) + dabs(h(l,l))
         if (s .eq. 0.0d0) s = norm
         tst1 = s
         tst2 = tst1 + dabs(h(l,l-1))
         if (tst2 .eq. tst1) go to 100
   80 continue
c     .......... form shift ..........
  100 x = h(en,en)
      if (l .eq. en) go to 270
      y = h(na,na)
      w = h(en,na) * h(na,en)
      if (l .eq. na) go to 280
      if (itn .eq. 0) go to 1000
      if (its .ne. 10 .and. its .ne. 20) go to 130
c     .......... form exceptional shift ..........
      t = t + x
c
      do 120 i = low, en
  120 h(i,i) = h(i,i) - x
c
      s = dabs(h(en,na)) + dabs(h(na,enm2))
      x = 0.75d0 * s
      y = x
      w = -0.4375d0 * s * s
  130 its = its + 1
      itn = itn - 1
c     .......... look for two consecutive small
c                sub-diagonal elements.
c                for m=en-2 step -1 until l do -- ..........
      do 140 mm = l, enm2
         m = enm2 + l - mm
         zz = h(m,m)
         r = x - zz
         s = y - zz
         p = (r * s - w) / h(m+1,m) + h(m,m+1)
         q = h(m+1,m+1) - zz - r - s
         r = h(m+2,m+1)
         s = dabs(p) + dabs(q) + dabs(r)
         p = p / s
         q = q / s
         r = r / s
         if (m .eq. l) go to 150
         tst1 = dabs(p)*(dabs(h(m-1,m-1)) + dabs(zz) + dabs(h(m+1,m+1)))
         tst2 = tst1 + dabs(h(m,m-1))*(dabs(q) + dabs(r))
         if (tst2 .eq. tst1) go to 150
  140 continue
c
  150 mp2 = m + 2
c
      do 160 i = mp2, en
         h(i,i-2) = 0.0d0
         if (i .eq. mp2) go to 160
         h(i,i-3) = 0.0d0
  160 continue
c     .......... double qr step involving rows l to en and
c                columns m to en ..........
      do 260 k = m, na
         notlas = k .ne. na
         if (k .eq. m) go to 170
         p = h(k,k-1)
         q = h(k+1,k-1)
         r = 0.0d0
         if (notlas) r = h(k+2,k-1)
         x = dabs(p) + dabs(q) + dabs(r)
         if (x .eq. 0.0d0) go to 260
         p = p / x
         q = q / x
         r = r / x
  170    s = dsign(dsqrt(p*p+q*q+r*r),p)
         if (k .eq. m) go to 180
         h(k,k-1) = -s * x
         go to 190
  180    if (l .ne. m) h(k,k-1) = -h(k,k-1)
  190    p = p + s
         x = p / s
         y = q / s
         zz = r / s
         q = q / p
         r = r / p
         if (notlas) go to 225
c     .......... row modification ..........
         do 200 j = k, EN
            p = h(k,j) + q * h(k+1,j)
            h(k,j) = h(k,j) - p * x
            h(k+1,j) = h(k+1,j) - p * y
  200    continue
c
         j = min0(en,k+3)
c     .......... column modification ..........
         do 210 i = L, j
            p = x * h(i,k) + y * h(i,k+1)
            h(i,k) = h(i,k) - p
            h(i,k+1) = h(i,k+1) - p * q
  210    continue
         go to 255
  225    continue
c     .......... row modification ..........
         do 230 j = k, EN
            p = h(k,j) + q * h(k+1,j) + r * h(k+2,j)
            h(k,j) = h(k,j) - p * x
            h(k+1,j) = h(k+1,j) - p * y
            h(k+2,j) = h(k+2,j) - p * zz
  230    continue
c
         j = min0(en,k+3)
c     .......... column modification ..........
         do 240 i = L, j
            p = x * h(i,k) + y * h(i,k+1) + zz * h(i,k+2)
            h(i,k) = h(i,k) - p
            h(i,k+1) = h(i,k+1) - p * q
            h(i,k+2) = h(i,k+2) - p * r
  240    continue
  255    continue
c
  260 continue
c
      go to 70
c     .......... one root found ..........
  270 wr(en) = x + t
      wi(en) = 0.0d0
      en = na
      go to 60
c     .......... two roots found ..........
  280 p = (y - x) / 2.0d0
      q = p * p + w
      zz = dsqrt(dabs(q))
      x = x + t
      if (q .lt. 0.0d0) go to 320
c     .......... real pair ..........
      zz = p + dsign(zz,p)
      wr(na) = x + zz
      wr(en) = wr(na)
      if (zz .ne. 0.0d0) wr(en) = x - w / zz
      wi(na) = 0.0d0
      wi(en) = 0.0d0
      go to 330
c     .......... complex pair ..........
  320 wr(na) = x + p
      wr(en) = x + p
      wi(na) = zz
      wi(en) = -zz
  330 en = enm2
      go to 60
c     .......... set error -- all eigenvalues have not
c                converged after 30*n iterations ..........
 1000 ierr = en
 1001 return
      end
      subroutine hqr2(nm,n,low,igh,h,wr,wi,z,ierr)
c
      integer i,j,k,l,m,n,en,ii,jj,ll,mm,na,nm,nn,
     x        igh,itn,its,low,mp2,enm2,ierr
      double precision h(nm,n),wr(n),wi(n),z(nm,n)
      double precision p,q,r,s,t,w,x,y,ra,sa,vi,vr,zz,norm,tst1,tst2
      logical notlas
c
c     this subroutine is a translation of the algol procedure hqr2,
c     num. math. 16, 181-204(1970) by peters and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
c
c     this subroutine finds the eigenvalues and eigenvectors
c     of a real upper hessenberg matrix by the qr method.  the
c     eigenvectors of a real general matrix can also be found
c     if  elmhes  and  eltran  or  orthes  and  ortran  have
c     been used to reduce this general matrix to hessenberg form
c     and to accumulate the similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by the balancing
c          subroutine  balanc.  if  balanc  has not been used,
c          set low=1, igh=n.
c
c        h contains the upper hessenberg matrix.
c
c        z contains the transformation matrix produced by  eltran
c          after the reduction by  elmhes, or by  ortran  after the
c          reduction by  orthes, if performed.  if the eigenvectors
c          of the hessenberg matrix are desired, z must contain the
c          identity matrix.
c
c     on output
c
c        h has been destroyed.
c
c        wr and wi contain the real and imaginary parts,
c          respectively, of the eigenvalues.  the eigenvalues
c          are unordered except that complex conjugate pairs
c          of values appear consecutively with the eigenvalue
c          having the positive imaginary part first.  if an
c          error exit is made, the eigenvalues should be correct
c          for indices ierr+1,...,n.
c
c        z contains the real and imaginary parts of the eigenvectors.
c          if the i-th eigenvalue is real, the i-th column of z
c          contains its eigenvector.  if the i-th eigenvalue is complex
c          with positive imaginary part, the i-th and (i+1)-th
c          columns of z contain the real and imaginary parts of its
c          eigenvector.  the eigenvectors are unnormalized.  if an
c          error exit is made, none of the eigenvectors has been found.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the limit of 30*n iterations is exhausted
c                     while the j-th eigenvalue is being sought.
c
c     calls cdiv for complex division.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      norm = 0.0d0
      k = 1
c     .......... store roots isolated by balanc
c                and compute matrix norm ..........
      do 50 i = 1, n
c
         do 40 j = k, n
   40    norm = norm + dabs(h(i,j))
c
         k = i
         if (i .ge. low .and. i .le. igh) go to 50
         wr(i) = h(i,i)
         wi(i) = 0.0d0
   50 continue
c
      en = igh
      t = 0.0d0
      itn = 30*n
c     .......... search for next eigenvalues ..........
   60 if (en .lt. low) go to 340
      its = 0
      na = en - 1
      enm2 = na - 1
c     .......... look for single small sub-diagonal element
c                for l=en step -1 until low do -- ..........
   70 do 80 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 100
         s = dabs(h(l-1,l-1)) + dabs(h(l,l))
         if (s .eq. 0.0d0) s = norm
         tst1 = s
         tst2 = tst1 + dabs(h(l,l-1))
         if (tst2 .eq. tst1) go to 100
   80 continue
c     .......... form shift ..........
  100 x = h(en,en)
      if (l .eq. en) go to 270
      y = h(na,na)
      w = h(en,na) * h(na,en)
      if (l .eq. na) go to 280
      if (itn .eq. 0) go to 1000
      if (its .ne. 10 .and. its .ne. 20) go to 130
c     .......... form exceptional shift ..........
      t = t + x
c
      do 120 i = low, en
  120 h(i,i) = h(i,i) - x
c
      s = dabs(h(en,na)) + dabs(h(na,enm2))
      x = 0.75d0 * s
      y = x
      w = -0.4375d0 * s * s
  130 its = its + 1
      itn = itn - 1
c     .......... look for two consecutive small
c                sub-diagonal elements.
c                for m=en-2 step -1 until l do -- ..........
      do 140 mm = l, enm2
         m = enm2 + l - mm
         zz = h(m,m)
         r = x - zz
         s = y - zz
         p = (r * s - w) / h(m+1,m) + h(m,m+1)
         q = h(m+1,m+1) - zz - r - s
         r = h(m+2,m+1)
         s = dabs(p) + dabs(q) + dabs(r)
         p = p / s
         q = q / s
         r = r / s
         if (m .eq. l) go to 150
         tst1 = dabs(p)*(dabs(h(m-1,m-1)) + dabs(zz) + dabs(h(m+1,m+1)))
         tst2 = tst1 + dabs(h(m,m-1))*(dabs(q) + dabs(r))
         if (tst2 .eq. tst1) go to 150
  140 continue
c
  150 mp2 = m + 2
c
      do 160 i = mp2, en
         h(i,i-2) = 0.0d0
         if (i .eq. mp2) go to 160
         h(i,i-3) = 0.0d0
  160 continue
c     .......... double qr step involving rows l to en and
c                columns m to en ..........
      do 260 k = m, na
         notlas = k .ne. na
         if (k .eq. m) go to 170
         p = h(k,k-1)
         q = h(k+1,k-1)
         r = 0.0d0
         if (notlas) r = h(k+2,k-1)
         x = dabs(p) + dabs(q) + dabs(r)
         if (x .eq. 0.0d0) go to 260
         p = p / x
         q = q / x
         r = r / x
  170    s = dsign(dsqrt(p*p+q*q+r*r),p)
         if (k .eq. m) go to 180
         h(k,k-1) = -s * x
         go to 190
  180    if (l .ne. m) h(k,k-1) = -h(k,k-1)
  190    p = p + s
         x = p / s
         y = q / s
         zz = r / s
         q = q / p
         r = r / p
         if (notlas) go to 225
c     .......... row modification ..........
         do 200 j = k, n
            p = h(k,j) + q * h(k+1,j)
            h(k,j) = h(k,j) - p * x
            h(k+1,j) = h(k+1,j) - p * y
  200    continue
c
         j = min0(en,k+3)
c     .......... column modification ..........
         do 210 i = 1, j
            p = x * h(i,k) + y * h(i,k+1)
            h(i,k) = h(i,k) - p
            h(i,k+1) = h(i,k+1) - p * q
  210    continue
c     .......... accumulate transformations ..........
         do 220 i = low, igh
            p = x * z(i,k) + y * z(i,k+1)
            z(i,k) = z(i,k) - p
            z(i,k+1) = z(i,k+1) - p * q
  220    continue
         go to 255
  225    continue
c     .......... row modification ..........
         do 230 j = k, n
            p = h(k,j) + q * h(k+1,j) + r * h(k+2,j)
            h(k,j) = h(k,j) - p * x
            h(k+1,j) = h(k+1,j) - p * y
            h(k+2,j) = h(k+2,j) - p * zz
  230    continue
c
         j = min0(en,k+3)
c     .......... column modification ..........
         do 240 i = 1, j
            p = x * h(i,k) + y * h(i,k+1) + zz * h(i,k+2)
            h(i,k) = h(i,k) - p
            h(i,k+1) = h(i,k+1) - p * q
            h(i,k+2) = h(i,k+2) - p * r
  240    continue
c     .......... accumulate transformations ..........
         do 250 i = low, igh
            p = x * z(i,k) + y * z(i,k+1) + zz * z(i,k+2)
            z(i,k) = z(i,k) - p
            z(i,k+1) = z(i,k+1) - p * q
            z(i,k+2) = z(i,k+2) - p * r
  250    continue
  255    continue
c
  260 continue
c
      go to 70
c     .......... one root found ..........
  270 h(en,en) = x + t
      wr(en) = h(en,en)
      wi(en) = 0.0d0
      en = na
      go to 60
c     .......... two roots found ..........
  280 p = (y - x) / 2.0d0
      q = p * p + w
      zz = dsqrt(dabs(q))
      h(en,en) = x + t
      x = h(en,en)
      h(na,na) = y + t
      if (q .lt. 0.0d0) go to 320
c     .......... real pair ..........
      zz = p + dsign(zz,p)
      wr(na) = x + zz
      wr(en) = wr(na)
      if (zz .ne. 0.0d0) wr(en) = x - w / zz
      wi(na) = 0.0d0
      wi(en) = 0.0d0
      x = h(en,na)
      s = dabs(x) + dabs(zz)
      p = x / s
      q = zz / s
      r = dsqrt(p*p+q*q)
      p = p / r
      q = q / r
c     .......... row modification ..........
      do 290 j = na, n
         zz = h(na,j)
         h(na,j) = q * zz + p * h(en,j)
         h(en,j) = q * h(en,j) - p * zz
  290 continue
c     .......... column modification ..........
      do 300 i = 1, en
         zz = h(i,na)
         h(i,na) = q * zz + p * h(i,en)
         h(i,en) = q * h(i,en) - p * zz
  300 continue
c     .......... accumulate transformations ..........
      do 310 i = low, igh
         zz = z(i,na)
         z(i,na) = q * zz + p * z(i,en)
         z(i,en) = q * z(i,en) - p * zz
  310 continue
c
      go to 330
c     .......... complex pair ..........
  320 wr(na) = x + p
      wr(en) = x + p
      wi(na) = zz
      wi(en) = -zz
  330 en = enm2
      go to 60
c     .......... all roots found.  backsubstitute to find
c                vectors of upper triangular form ..........
  340 if (norm .eq. 0.0d0) go to 1001
c     .......... for en=n step -1 until 1 do -- ..........
      do 800 nn = 1, n
         en = n + 1 - nn
         p = wr(en)
         q = wi(en)
         na = en - 1
         if (q) 710, 600, 800
c     .......... real vector ..........
  600    m = en
         h(en,en) = 1.0d0
         if (na .eq. 0) go to 800
c     .......... for i=en-1 step -1 until 1 do -- ..........
         do 700 ii = 1, na
            i = en - ii
            w = h(i,i) - p
            r = 0.0d0
c
            do 610 j = m, en
  610       r = r + h(i,j) * h(j,en)
c
            if (wi(i) .ge. 0.0d0) go to 630
            zz = w
            s = r
            go to 700
  630       m = i
            if (wi(i) .ne. 0.0d0) go to 640
            t = w
            if (t .ne. 0.0d0) go to 635
               tst1 = norm
               t = tst1
  632          t = 0.01d0 * t
               tst2 = norm + t
               if (tst2 .gt. tst1) go to 632
  635       h(i,en) = -r / t
            go to 680
c     .......... solve real equations ..........
  640       x = h(i,i+1)
            y = h(i+1,i)
            q = (wr(i) - p) * (wr(i) - p) + wi(i) * wi(i)
            t = (x * s - zz * r) / q
            h(i,en) = t
            if (dabs(x) .le. dabs(zz)) go to 650
            h(i+1,en) = (-r - w * t) / x
            go to 680
  650       h(i+1,en) = (-s - y * t) / zz
c
c     .......... overflow control ..........
  680       t = dabs(h(i,en))
            if (t .eq. 0.0d0) go to 700
            tst1 = t
            tst2 = tst1 + 1.0d0/tst1
            if (tst2 .gt. tst1) go to 700
            do 690 j = i, en
               h(j,en) = h(j,en)/t
  690       continue
c
  700    continue
c     .......... end real vector ..........
         go to 800
c     .......... complex vector ..........
  710    m = na
c     .......... last vector component chosen imaginary so that
c                eigenvector matrix is triangular ..........
         if (dabs(h(en,na)) .le. dabs(h(na,en))) go to 720
         h(na,na) = q / h(en,na)
         h(na,en) = -(h(en,en) - p) / h(en,na)
         go to 730
  720    call cdiv(0.0d0,-h(na,en),h(na,na)-p,q,h(na,na),h(na,en))
  730    h(en,na) = 0.0d0
         h(en,en) = 1.0d0
         enm2 = na - 1
         if (enm2 .eq. 0) go to 800
c     .......... for i=en-2 step -1 until 1 do -- ..........
         do 795 ii = 1, enm2
            i = na - ii
            w = h(i,i) - p
            ra = 0.0d0
            sa = 0.0d0
c
            do 760 j = m, en
               ra = ra + h(i,j) * h(j,na)
               sa = sa + h(i,j) * h(j,en)
  760       continue
c
            if (wi(i) .ge. 0.0d0) go to 770
            zz = w
            r = ra
            s = sa
            go to 795
  770       m = i
            if (wi(i) .ne. 0.0d0) go to 780
            call cdiv(-ra,-sa,w,q,h(i,na),h(i,en))
            go to 790
c     .......... solve complex equations ..........
  780       x = h(i,i+1)
            y = h(i+1,i)
            vr = (wr(i) - p) * (wr(i) - p) + wi(i) * wi(i) - q * q
            vi = (wr(i) - p) * 2.0d0 * q
            if (vr .ne. 0.0d0 .or. vi .ne. 0.0d0) go to 784
               tst1 = norm * (dabs(w) + dabs(q) + dabs(x)
     x                      + dabs(y) + dabs(zz))
               vr = tst1
  783          vr = 0.01d0 * vr
               tst2 = tst1 + vr
               if (tst2 .gt. tst1) go to 783
  784       call cdiv(x*r-zz*ra+q*sa,x*s-zz*sa-q*ra,vr,vi,
     x                h(i,na),h(i,en))
            if (dabs(x) .le. dabs(zz) + dabs(q)) go to 785
            h(i+1,na) = (-ra - w * h(i,na) + q * h(i,en)) / x
            h(i+1,en) = (-sa - w * h(i,en) - q * h(i,na)) / x
            go to 790
  785       call cdiv(-r-y*h(i,na),-s-y*h(i,en),zz,q,
     x                h(i+1,na),h(i+1,en))
c
c     .......... overflow control ..........
  790       t = dmax1(dabs(h(i,na)), dabs(h(i,en)))
            if (t .eq. 0.0d0) go to 795
            tst1 = t
            tst2 = tst1 + 1.0d0/tst1
            if (tst2 .gt. tst1) go to 795
            do 792 j = i, en
               h(j,na) = h(j,na)/t
               h(j,en) = h(j,en)/t
  792       continue
c
  795    continue
c     .......... end complex vector ..........
  800 continue
c     .......... end back substitution.
c                vectors of isolated roots ..........
      do 840 i = 1, n
         if (i .ge. low .and. i .le. igh) go to 840
c
         do 820 j = i, n
  820    z(i,j) = h(i,j)
c
  840 continue
c     .......... multiply by transformation matrix to give
c                vectors of original full matrix.
c                for j=n step -1 until low do -- ..........
      do 880 jj = low, n
         j = n + low - jj
         m = min0(j,igh)
c
         do 880 i = low, igh
            zz = 0.0d0
c
            do 860 k = low, m
  860       zz = zz + z(i,k) * h(k,j)
c
            z(i,j) = zz
  880 continue
c
      go to 1001
c     .......... set error -- all eigenvalues have not
c                converged after 30*n iterations ..........
 1000 ierr = en
 1001 return
      end
c---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine cg_eispack(nm,n,ar,ai,wr,wi,matz,zr,zi,fv1,fv2,fv3
     &     ,ierr)
      integer n,nm,is1,is2,ierr,matz
      double precision ar(nm,n),ai(nm,n),wr(n),wi(n),zr(nm,n),zi(nm,n),
     x       fv1(n),fv2(n),fv3(n)
c
c     this subroutine calls the recommended sequence of
c     subroutines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     of a complex general matrix.
c
c     on input
c
c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.
c
c        n  is the order of the matrix  a=(ar,ai).
c
c        ar  and  ai  contain the real and imaginary parts,
c        respectively, of the complex general matrix.
c
c        matz  is an integer variable set equal to zero if
c        only eigenvalues are desired.  otherwise it is set to
c        any non-zero integer for both eigenvalues and eigenvectors.
c
c     on output
c
c        wr  and  wi  contain the real and imaginary parts,
c        respectively, of the eigenvalues.
c
c        zr  and  zi  contain the real and imaginary parts,
c        respectively, of the eigenvectors if matz is not zero.
c
c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for comqr
c           and comqr2.  the normal completion code is zero.
c
c        fv1, fv2, and  fv3  are temporary storage arrays.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
c
   10 call  cbal(nm,n,ar,ai,is1,is2,fv1)
      call  corth(nm,n,is1,is2,ar,ai,fv2,fv3)
      if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      call  comqr(nm,n,is1,is2,ar,ai,wr,wi,ierr)
      go to 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 call  comqr2(nm,n,is1,is2,fv2,fv3,ar,ai,wr,wi,zr,zi,ierr)
      if (ierr .ne. 0) go to 50
      call  cbabk2(nm,n,is1,is2,fv1,n,zr,zi)
   50 return
      end
      subroutine cbabk2(nm,n,low,igh,scale,m,zr,zi)
c
      integer i,j,k,m,n,ii,nm,igh,low
      double precision scale(n),zr(nm,m),zi(nm,m)
      double precision s
c
c     this subroutine is a translation of the algol procedure
c     cbabk2, which is a complex version of balbak,
c     num. math. 13, 293-304(1969) by parlett and reinsch.
c     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
c
c     this subroutine forms the eigenvectors of a complex general
c     matrix by back transforming those of the corresponding
c     balanced matrix determined by  cbal.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by  cbal.
c
c        scale contains information determining the permutations
c          and scaling factors used by  cbal.
c
c        m is the number of eigenvectors to be back transformed.
c
c        zr and zi contain the real and imaginary parts,
c          respectively, of the eigenvectors to be
c          back transformed in their first m columns.
c
c     on output
c
c        zr and zi contain the real and imaginary parts,
c          respectively, of the transformed eigenvectors
c          in their first m columns.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      if (m .eq. 0) go to 200
      if (igh .eq. low) go to 120
c
      do 110 i = low, igh
         s = scale(i)
c     .......... left hand eigenvectors are back transformed
c                if the foregoing statement is replaced by
c                s=1.0d0/scale(i). ..........
         do 100 j = 1, m
            zr(i,j) = zr(i,j) * s
            zi(i,j) = zi(i,j) * s
  100    continue
c
  110 continue
c     .......... for i=low-1 step -1 until 1,
c                igh+1 step 1 until n do -- ..........
  120 do 140 ii = 1, n
         i = ii
         if (i .ge. low .and. i .le. igh) go to 140
         if (i .lt. low) i = low - ii
         k = scale(i)
         if (k .eq. i) go to 140
c
         do 130 j = 1, m
            s = zr(i,j)
            zr(i,j) = zr(k,j)
            zr(k,j) = s
            s = zi(i,j)
            zi(i,j) = zi(k,j)
            zi(k,j) = s
  130    continue
c
  140 continue
c
  200 return
      end
      subroutine cbal(nm,n,ar,ai,low,igh,scale)
c
      integer i,j,k,l,m,n,jj,nm,igh,low,iexc
      double precision ar(nm,n),ai(nm,n),scale(n)
      double precision c,f,g,r,s,b2,radix
      logical noconv
c
c     this subroutine is a translation of the algol procedure
c     cbalance, which is a complex version of balance,
c     num. math. 13, 293-304(1969) by parlett and reinsch.
c     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
c
c     this subroutine balances a complex matrix and isolates
c     eigenvalues whenever possible.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        ar and ai contain the real and imaginary parts,
c          respectively, of the complex matrix to be balanced.
c
c     on output
c
c        ar and ai contain the real and imaginary parts,
c          respectively, of the balanced matrix.
c
c        low and igh are two integers such that ar(i,j) and ai(i,j)
c          are equal to zero if
c           (1) i is greater than j and
c           (2) j=1,...,low-1 or i=igh+1,...,n.
c
c        scale contains information determining the
c           permutations and scaling factors used.
c
c     suppose that the principal submatrix in rows low through igh
c     has been balanced, that p(j) denotes the index interchanged
c     with j during the permutation step, and that the elements
c     of the diagonal matrix used are denoted by d(i,j).  then
c        scale(j) = p(j),    for j = 1,...,low-1
c                 = d(j,j)       j = low,...,igh
c                 = p(j)         j = igh+1,...,n.
c     the order in which the interchanges are made is n to igh+1,
c     then 1 to low-1.
c
c     note that 1 is returned for igh if igh is zero formally.
c
c     the algol procedure exc contained in cbalance appears in
c     cbal  in line.  (note that the algol roles of identifiers
c     k,l have been reversed.)
c
c     arithmetic is real throughout.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      radix = 16.0d0
c
      b2 = radix * radix
      k = 1
      l = n
      go to 100
c     .......... in-line procedure for row and
c                column exchange ..........
   20 scale(m) = j
      if (j .eq. m) go to 50
c
      do 30 i = 1, l
         f = ar(i,j)
         ar(i,j) = ar(i,m)
         ar(i,m) = f
         f = ai(i,j)
         ai(i,j) = ai(i,m)
         ai(i,m) = f
   30 continue
c
      do 40 i = k, n
         f = ar(j,i)
         ar(j,i) = ar(m,i)
         ar(m,i) = f
         f = ai(j,i)
         ai(j,i) = ai(m,i)
         ai(m,i) = f
   40 continue
c
   50 go to (80,130), iexc
c     .......... search for rows isolating an eigenvalue
c                and push them down ..........
   80 if (l .eq. 1) go to 280
      l = l - 1
c     .......... for j=l step -1 until 1 do -- ..........
  100 do 120 jj = 1, l
         j = l + 1 - jj
c
         do 110 i = 1, l
            if (i .eq. j) go to 110
            if (ar(j,i) .ne. 0.0d0 .or. ai(j,i) .ne. 0.0d0) go to 120
  110    continue
c
         m = l
         iexc = 1
         go to 20
  120 continue
c
      go to 140
c     .......... search for columns isolating an eigenvalue
c                and push them left ..........
  130 k = k + 1
c
  140 do 170 j = k, l
c
         do 150 i = k, l
            if (i .eq. j) go to 150
            if (ar(i,j) .ne. 0.0d0 .or. ai(i,j) .ne. 0.0d0) go to 170
  150    continue
c
         m = k
         iexc = 2
         go to 20
  170 continue
c     .......... now balance the submatrix in rows k to l ..........
      do 180 i = k, l
  180 scale(i) = 1.0d0
c     .......... iterative loop for norm reduction ..........
  190 noconv = .false.
c
      do 270 i = k, l
         c = 0.0d0
         r = 0.0d0
c
         do 200 j = k, l
            if (j .eq. i) go to 200
            c = c + dabs(ar(j,i)) + dabs(ai(j,i))
            r = r + dabs(ar(i,j)) + dabs(ai(i,j))
  200    continue
c     .......... guard against zero c or r due to underflow ..........
         if (c .eq. 0.0d0 .or. r .eq. 0.0d0) go to 270
         g = r / radix
         f = 1.0d0
         s = c + r
  210    if (c .ge. g) go to 220
         f = f * radix
         c = c * b2
         go to 210
  220    g = r * radix
  230    if (c .lt. g) go to 240
         f = f / radix
         c = c / b2
         go to 230
c     .......... now balance ..........
  240    if ((c + r) / f .ge. 0.95d0 * s) go to 270
         g = 1.0d0 / f
         scale(i) = scale(i) * f
         noconv = .true.
c
         do 250 j = k, n
            ar(i,j) = ar(i,j) * g
            ai(i,j) = ai(i,j) * g
  250    continue
c
         do 260 j = 1, l
            ar(j,i) = ar(j,i) * f
            ai(j,i) = ai(j,i) * f
  260    continue
c
  270 continue
c
      if (noconv) go to 190
c
  280 low = k
      igh = l
      return
      end
      subroutine comqr(nm,n,low,igh,hr,hi,wr,wi,ierr)
c
      integer i,j,l,n,en,ll,nm,igh,itn,its,low,lp1,enm1,ierr
      double precision hr(nm,n),hi(nm,n),wr(n),wi(n)
      double precision si,sr,ti,tr,xi,xr,yi,yr,zzi,zzr,norm,tst1,tst2,
     x       pythag
c
c     this subroutine is a translation of a unitary analogue of the
c     algol procedure  comlr, num. math. 12, 369-376(1968) by martin
c     and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 396-403(1971).
c     the unitary analogue substitutes the qr algorithm of francis
c     (comp. jour. 4, 332-345(1962)) for the lr algorithm.
c
c     this subroutine finds the eigenvalues of a complex
c     upper hessenberg matrix by the qr method.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by the balancing
c          subroutine  cbal.  if  cbal  has not been used,
c          set low=1, igh=n.
c
c        hr and hi contain the real and imaginary parts,
c          respectively, of the complex upper hessenberg matrix.
c          their lower triangles below the subdiagonal contain
c          information about the unitary transformations used in
c          the reduction by  corth, if performed.
c
c     on output
c
c        the upper hessenberg portions of hr and hi have been
c          destroyed.  therefore, they must be saved before
c          calling  comqr  if subsequent calculation of
c          eigenvectors is to be performed.
c
c        wr and wi contain the real and imaginary parts,
c          respectively, of the eigenvalues.  if an error
c          exit is made, the eigenvalues should be correct
c          for indices ierr+1,...,n.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the limit of 30*n iterations is exhausted
c                     while the j-th eigenvalue is being sought.
c
c     calls cdiv for complex division.
c     calls csroot for complex square root.
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      if (low .eq. igh) go to 180
c     .......... create real subdiagonal elements ..........
      l = low + 1
c
      do 170 i = l, igh
         ll = min0(i+1,igh)
         if (hi(i,i-1) .eq. 0.0d0) go to 170
         norm = pythag(hr(i,i-1),hi(i,i-1))
         yr = hr(i,i-1) / norm
         yi = hi(i,i-1) / norm
         hr(i,i-1) = norm
         hi(i,i-1) = 0.0d0
c
         do 155 j = i, igh
            si = yr * hi(i,j) - yi * hr(i,j)
            hr(i,j) = yr * hr(i,j) + yi * hi(i,j)
            hi(i,j) = si
  155    continue
c
         do 160 j = low, ll
            si = yr * hi(j,i) + yi * hr(j,i)
            hr(j,i) = yr * hr(j,i) - yi * hi(j,i)
            hi(j,i) = si
  160    continue
c
  170 continue
c     .......... store roots isolated by cbal ..........
  180 do 200 i = 1, n
         if (i .ge. low .and. i .le. igh) go to 200
         wr(i) = hr(i,i)
         wi(i) = hi(i,i)
  200 continue
c
      en = igh
      tr = 0.0d0
      ti = 0.0d0
      itn = 30*n
c     .......... search for next eigenvalue ..........
  220 if (en .lt. low) go to 1001
      its = 0
      enm1 = en - 1
c     .......... look for single small sub-diagonal element
c                for l=en step -1 until low d0 -- ..........
  240 do 260 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 300
         tst1 = dabs(hr(l-1,l-1)) + dabs(hi(l-1,l-1))
     x            + dabs(hr(l,l)) + dabs(hi(l,l))
         tst2 = tst1 + dabs(hr(l,l-1))
         if (tst2 .eq. tst1) go to 300
  260 continue
c     .......... form shift ..........
  300 if (l .eq. en) go to 660
      if (itn .eq. 0) go to 1000
      if (its .eq. 10 .or. its .eq. 20) go to 320
      sr = hr(en,en)
      si = hi(en,en)
      xr = hr(enm1,en) * hr(en,enm1)
      xi = hi(enm1,en) * hr(en,enm1)
      if (xr .eq. 0.0d0 .and. xi .eq. 0.0d0) go to 340
      yr = (hr(enm1,enm1) - sr) / 2.0d0
      yi = (hi(enm1,enm1) - si) / 2.0d0
      call csroot(yr**2-yi**2+xr,2.0d0*yr*yi+xi,zzr,zzi)
      if (yr * zzr + yi * zzi .ge. 0.0d0) go to 310
      zzr = -zzr
      zzi = -zzi
  310 call cdiv(xr,xi,yr+zzr,yi+zzi,xr,xi)
      sr = sr - xr
      si = si - xi
      go to 340
c     .......... form exceptional shift ..........
  320 sr = dabs(hr(en,enm1)) + dabs(hr(enm1,en-2))
      si = 0.0d0
c
  340 do 360 i = low, en
         hr(i,i) = hr(i,i) - sr
         hi(i,i) = hi(i,i) - si
  360 continue
c
      tr = tr + sr
      ti = ti + si
      its = its + 1
      itn = itn - 1
c     .......... reduce to triangle (rows) ..........
      lp1 = l + 1
c
      do 500 i = lp1, en
         sr = hr(i,i-1)
         hr(i,i-1) = 0.0d0
         norm = pythag(pythag(hr(i-1,i-1),hi(i-1,i-1)),sr)
         xr = hr(i-1,i-1) / norm
         wr(i-1) = xr
         xi = hi(i-1,i-1) / norm
         wi(i-1) = xi
         hr(i-1,i-1) = norm
         hi(i-1,i-1) = 0.0d0
         hi(i,i-1) = sr / norm
c
         do 490 j = i, en
            yr = hr(i-1,j)
            yi = hi(i-1,j)
            zzr = hr(i,j)
            zzi = hi(i,j)
            hr(i-1,j) = xr * yr + xi * yi + hi(i,i-1) * zzr
            hi(i-1,j) = xr * yi - xi * yr + hi(i,i-1) * zzi
            hr(i,j) = xr * zzr - xi * zzi - hi(i,i-1) * yr
            hi(i,j) = xr * zzi + xi * zzr - hi(i,i-1) * yi
  490    continue
c
  500 continue
c
      si = hi(en,en)
      if (si .eq. 0.0d0) go to 540
      norm = pythag(hr(en,en),si)
      sr = hr(en,en) / norm
      si = si / norm
      hr(en,en) = norm
      hi(en,en) = 0.0d0
c     .......... inverse operation (columns) ..........
  540 do 600 j = lp1, en
         xr = wr(j-1)
         xi = wi(j-1)
c
         do 580 i = l, j
            yr = hr(i,j-1)
            yi = 0.0d0
            zzr = hr(i,j)
            zzi = hi(i,j)
            if (i .eq. j) go to 560
            yi = hi(i,j-1)
            hi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
  560       hr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            hr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            hi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
  580    continue
c
  600 continue
c
      if (si .eq. 0.0d0) go to 240
c
      do 630 i = l, en
         yr = hr(i,en)
         yi = hi(i,en)
         hr(i,en) = sr * yr - si * yi
         hi(i,en) = sr * yi + si * yr
  630 continue
c
      go to 240
c     .......... a root found ..........
  660 wr(en) = hr(en,en) + tr
      wi(en) = hi(en,en) + ti
      en = enm1
      go to 220
c     .......... set error -- all eigenvalues have not
c                converged after 30*n iterations ..........
 1000 ierr = en
 1001 return
      end
      subroutine comqr2(nm,n,low,igh,ortr,orti,hr,hi,wr,wi,zr,zi,ierr)
C  MESHED overflow control WITH vectors of isolated roots (10/19/89 BSG)
C  MESHED overflow control WITH triangular multiply (10/30/89 BSG)
c
      integer i,j,k,l,m,n,en,ii,jj,ll,nm,nn,igh,ip1,
     x        itn,its,low,lp1,enm1,iend,ierr
      double precision hr(nm,n),hi(nm,n),wr(n),wi(n),zr(nm,n),zi(nm,n),
     x       ortr(igh),orti(igh)
      double precision si,sr,ti,tr,xi,xr,yi,yr,zzi,zzr,norm,tst1,tst2,
     x       pythag
c
c     this subroutine is a translation of a unitary analogue of the
c     algol procedure  comlr2, num. math. 16, 181-204(1970) by peters
c     and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
c     the unitary analogue substitutes the qr algorithm of francis
c     (comp. jour. 4, 332-345(1962)) for the lr algorithm.
c
c     this subroutine finds the eigenvalues and eigenvectors
c     of a complex upper hessenberg matrix by the qr
c     method.  the eigenvectors of a complex general matrix
c     can also be found if  corth  has been used to reduce
c     this general matrix to hessenberg form.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by the balancing
c          subroutine  cbal.  if  cbal  has not been used,
c          set low=1, igh=n.
c
c        ortr and orti contain information about the unitary trans-
c          formations used in the reduction by  corth, if performed.
c          only elements low through igh are used.  if the eigenvectors
c          of the hessenberg matrix are desired, set ortr(j) and
c          orti(j) to 0.0d0 for these elements.
c
c        hr and hi contain the real and imaginary parts,
c          respectively, of the complex upper hessenberg matrix.
c          their lower triangles below the subdiagonal contain further
c          information about the transformations which were used in the
c          reduction by  corth, if performed.  if the eigenvectors of
c          the hessenberg matrix are desired, these elements may be
c          arbitrary.
c
c     on output
c
c        ortr, orti, and the upper hessenberg portions of hr and hi
c          have been destroyed.
c
c        wr and wi contain the real and imaginary parts,
c          respectively, of the eigenvalues.  if an error
c          exit is made, the eigenvalues should be correct
c          for indices ierr+1,...,n.
c
c        zr and zi contain the real and imaginary parts,
c          respectively, of the eigenvectors.  the eigenvectors
c          are unnormalized.  if an error exit is made, none of
c          the eigenvectors has been found.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the limit of 30*n iterations is exhausted
c                     while the j-th eigenvalue is being sought.
c
c     calls cdiv for complex division.
c     calls csroot for complex square root.
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated october 1989.
c
c     ------------------------------------------------------------------
c
      ierr = 0
c     .......... initialize eigenvector matrix ..........
      do 101 j = 1, n
c
         do 100 i = 1, n
            zr(i,j) = 0.0d0
            zi(i,j) = 0.0d0
  100    continue
         zr(j,j) = 1.0d0
  101 continue
c     .......... form the matrix of accumulated transformations
c                from the information left by corth ..........
      iend = igh - low - 1
      if (iend) 180, 150, 105
c     .......... for i=igh-1 step -1 until low+1 do -- ..........
  105 do 140 ii = 1, iend
         i = igh - ii
         if (ortr(i) .eq. 0.0d0 .and. orti(i) .eq. 0.0d0) go to 140
         if (hr(i,i-1) .eq. 0.0d0 .and. hi(i,i-1) .eq. 0.0d0) go to 140
c     .......... norm below is negative of h formed in corth ..........
         norm = hr(i,i-1) * ortr(i) + hi(i,i-1) * orti(i)
         ip1 = i + 1
c
         do 110 k = ip1, igh
            ortr(k) = hr(k,i-1)
            orti(k) = hi(k,i-1)
  110    continue
c
         do 130 j = i, igh
            sr = 0.0d0
            si = 0.0d0
c
            do 115 k = i, igh
               sr = sr + ortr(k) * zr(k,j) + orti(k) * zi(k,j)
               si = si + ortr(k) * zi(k,j) - orti(k) * zr(k,j)
  115       continue
c
            sr = sr / norm
            si = si / norm
c
            do 120 k = i, igh
               zr(k,j) = zr(k,j) + sr * ortr(k) - si * orti(k)
               zi(k,j) = zi(k,j) + sr * orti(k) + si * ortr(k)
  120       continue
c
  130    continue
c
  140 continue
c     .......... create real subdiagonal elements ..........
  150 l = low + 1
c
      do 170 i = l, igh
         ll = min0(i+1,igh)
         if (hi(i,i-1) .eq. 0.0d0) go to 170
         norm = pythag(hr(i,i-1),hi(i,i-1))
         yr = hr(i,i-1) / norm
         yi = hi(i,i-1) / norm
         hr(i,i-1) = norm
         hi(i,i-1) = 0.0d0
c
         do 155 j = i, n
            si = yr * hi(i,j) - yi * hr(i,j)
            hr(i,j) = yr * hr(i,j) + yi * hi(i,j)
            hi(i,j) = si
  155    continue
c
         do 160 j = 1, ll
            si = yr * hi(j,i) + yi * hr(j,i)
            hr(j,i) = yr * hr(j,i) - yi * hi(j,i)
            hi(j,i) = si
  160    continue
c
         do 165 j = low, igh
            si = yr * zi(j,i) + yi * zr(j,i)
            zr(j,i) = yr * zr(j,i) - yi * zi(j,i)
            zi(j,i) = si
  165    continue
c
  170 continue
c     .......... store roots isolated by cbal ..........
  180 do 200 i = 1, n
         if (i .ge. low .and. i .le. igh) go to 200
         wr(i) = hr(i,i)
         wi(i) = hi(i,i)
  200 continue
c
      en = igh
      tr = 0.0d0
      ti = 0.0d0
      itn = 30*n
c     .......... search for next eigenvalue ..........
  220 if (en .lt. low) go to 680
      its = 0
      enm1 = en - 1
c     .......... look for single small sub-diagonal element
c                for l=en step -1 until low do -- ..........
  240 do 260 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 300
         tst1 = dabs(hr(l-1,l-1)) + dabs(hi(l-1,l-1))
     x            + dabs(hr(l,l)) + dabs(hi(l,l))
         tst2 = tst1 + dabs(hr(l,l-1))
         if (tst2 .eq. tst1) go to 300
  260 continue
c     .......... form shift ..........
  300 if (l .eq. en) go to 660
      if (itn .eq. 0) go to 1000
      if (its .eq. 10 .or. its .eq. 20) go to 320
      sr = hr(en,en)
      si = hi(en,en)
      xr = hr(enm1,en) * hr(en,enm1)
      xi = hi(enm1,en) * hr(en,enm1)
      if (xr .eq. 0.0d0 .and. xi .eq. 0.0d0) go to 340
      yr = (hr(enm1,enm1) - sr) / 2.0d0
      yi = (hi(enm1,enm1) - si) / 2.0d0
      call csroot(yr**2-yi**2+xr,2.0d0*yr*yi+xi,zzr,zzi)
      if (yr * zzr + yi * zzi .ge. 0.0d0) go to 310
      zzr = -zzr
      zzi = -zzi
  310 call cdiv(xr,xi,yr+zzr,yi+zzi,xr,xi)
      sr = sr - xr
      si = si - xi
      go to 340
c     .......... form exceptional shift ..........
  320 sr = dabs(hr(en,enm1)) + dabs(hr(enm1,en-2))
      si = 0.0d0
c
  340 do 360 i = low, en
         hr(i,i) = hr(i,i) - sr
         hi(i,i) = hi(i,i) - si
  360 continue
c
      tr = tr + sr
      ti = ti + si
      its = its + 1
      itn = itn - 1
c     .......... reduce to triangle (rows) ..........
      lp1 = l + 1
c
      do 500 i = lp1, en
         sr = hr(i,i-1)
         hr(i,i-1) = 0.0d0
         norm = pythag(pythag(hr(i-1,i-1),hi(i-1,i-1)),sr)
         xr = hr(i-1,i-1) / norm
         wr(i-1) = xr
         xi = hi(i-1,i-1) / norm
         wi(i-1) = xi
         hr(i-1,i-1) = norm
         hi(i-1,i-1) = 0.0d0
         hi(i,i-1) = sr / norm
c
         do 490 j = i, n
            yr = hr(i-1,j)
            yi = hi(i-1,j)
            zzr = hr(i,j)
            zzi = hi(i,j)
            hr(i-1,j) = xr * yr + xi * yi + hi(i,i-1) * zzr
            hi(i-1,j) = xr * yi - xi * yr + hi(i,i-1) * zzi
            hr(i,j) = xr * zzr - xi * zzi - hi(i,i-1) * yr
            hi(i,j) = xr * zzi + xi * zzr - hi(i,i-1) * yi
  490    continue
c
  500 continue
c
      si = hi(en,en)
      if (si .eq. 0.0d0) go to 540
      norm = pythag(hr(en,en),si)
      sr = hr(en,en) / norm
      si = si / norm
      hr(en,en) = norm
      hi(en,en) = 0.0d0
      if (en .eq. n) go to 540
      ip1 = en + 1
c
      do 520 j = ip1, n
         yr = hr(en,j)
         yi = hi(en,j)
         hr(en,j) = sr * yr + si * yi
         hi(en,j) = sr * yi - si * yr
  520 continue
c     .......... inverse operation (columns) ..........
  540 do 600 j = lp1, en
         xr = wr(j-1)
         xi = wi(j-1)
c
         do 580 i = 1, j
            yr = hr(i,j-1)
            yi = 0.0d0
            zzr = hr(i,j)
            zzi = hi(i,j)
            if (i .eq. j) go to 560
            yi = hi(i,j-1)
            hi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
  560       hr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            hr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            hi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
  580    continue
c
         do 590 i = low, igh
            yr = zr(i,j-1)
            yi = zi(i,j-1)
            zzr = zr(i,j)
            zzi = zi(i,j)
            zr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            zi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
            zr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            zi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
  590    continue
c
  600 continue
c
      if (si .eq. 0.0d0) go to 240
c
      do 630 i = 1, en
         yr = hr(i,en)
         yi = hi(i,en)
         hr(i,en) = sr * yr - si * yi
         hi(i,en) = sr * yi + si * yr
  630 continue
c
      do 640 i = low, igh
         yr = zr(i,en)
         yi = zi(i,en)
         zr(i,en) = sr * yr - si * yi
         zi(i,en) = sr * yi + si * yr
  640 continue
c
      go to 240
c     .......... a root found ..........
  660 hr(en,en) = hr(en,en) + tr
      wr(en) = hr(en,en)
      hi(en,en) = hi(en,en) + ti
      wi(en) = hi(en,en)
      en = enm1
      go to 220
c     .......... all roots found.  backsubstitute to find
c                vectors of upper triangular form ..........
  680 norm = 0.0d0
c
      do 720 i = 1, n
c
         do 720 j = i, n
            tr = dabs(hr(i,j)) + dabs(hi(i,j))
            if (tr .gt. norm) norm = tr
  720 continue
c
      if (n .eq. 1 .or. norm .eq. 0.0d0) go to 1001
c     .......... for en=n step -1 until 2 do -- ..........
      do 800 nn = 2, n
         en = n + 2 - nn
         xr = wr(en)
         xi = wi(en)
         hr(en,en) = 1.0d0
         hi(en,en) = 0.0d0
         enm1 = en - 1
c     .......... for i=en-1 step -1 until 1 do -- ..........
         do 780 ii = 1, enm1
            i = en - ii
            zzr = 0.0d0
            zzi = 0.0d0
            ip1 = i + 1
c
            do 740 j = ip1, en
               zzr = zzr + hr(i,j) * hr(j,en) - hi(i,j) * hi(j,en)
               zzi = zzi + hr(i,j) * hi(j,en) + hi(i,j) * hr(j,en)
  740       continue
c
            yr = xr - wr(i)
            yi = xi - wi(i)
            if (yr .ne. 0.0d0 .or. yi .ne. 0.0d0) go to 765
               tst1 = norm
               yr = tst1
  760          yr = 0.01d0 * yr
               tst2 = norm + yr
               if (tst2 .gt. tst1) go to 760
  765       continue
            call cdiv(zzr,zzi,yr,yi,hr(i,en),hi(i,en))
c     .......... overflow control ..........
            tr = dabs(hr(i,en)) + dabs(hi(i,en))
            if (tr .eq. 0.0d0) go to 780
            tst1 = tr
            tst2 = tst1 + 1.0d0/tst1
            if (tst2 .gt. tst1) go to 780
            do 770 j = i, en
               hr(j,en) = hr(j,en)/tr
               hi(j,en) = hi(j,en)/tr
  770       continue
c
  780    continue
c
  800 continue
c     .......... end backsubstitution ..........
c     .......... vectors of isolated roots ..........
      do  840 i = 1, N
         if (i .ge. low .and. i .le. igh) go to 840
c
         do 820 j = I, n
            zr(i,j) = hr(i,j)
            zi(i,j) = hi(i,j)
  820    continue
c
  840 continue
c     .......... multiply by transformation matrix to give
c                vectors of original full matrix.
c                for j=n step -1 until low do -- ..........
      do 880 jj = low, N
         j = n + low - jj
         m = min0(j,igh)
c
         do 880 i = low, igh
            zzr = 0.0d0
            zzi = 0.0d0
c
            do 860 k = low, m
               zzr = zzr + zr(i,k) * hr(k,j) - zi(i,k) * hi(k,j)
               zzi = zzi + zr(i,k) * hi(k,j) + zi(i,k) * hr(k,j)
  860       continue
c
            zr(i,j) = zzr
            zi(i,j) = zzi
  880 continue
c
      go to 1001
c     .......... set error -- all eigenvalues have not
c                converged after 30*n iterations ..........
 1000 ierr = en
 1001 return
      end
      subroutine corth(nm,n,low,igh,ar,ai,ortr,orti)
c
      integer i,j,m,n,ii,jj,la,mp,nm,igh,kp1,low
      double precision ar(nm,n),ai(nm,n),ortr(igh),orti(igh)
      double precision f,g,h,fi,fr,scale,pythag
c
c     this subroutine is a translation of a complex analogue of
c     the algol procedure orthes, num. math. 12, 349-368(1968)
c     by martin and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
c
c     given a complex general matrix, this subroutine
c     reduces a submatrix situated in rows and columns
c     low through igh to upper hessenberg form by
c     unitary similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by the balancing
c          subroutine  cbal.  if  cbal  has not been used,
c          set low=1, igh=n.
c
c        ar and ai contain the real and imaginary parts,
c          respectively, of the complex input matrix.
c
c     on output
c
c        ar and ai contain the real and imaginary parts,
c          respectively, of the hessenberg matrix.  information
c          about the unitary transformations used in the reduction
c          is stored in the remaining triangles under the
c          hessenberg matrix.
c
c        ortr and orti contain further information about the
c          transformations.  only elements low through igh are used.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
c
      do 180 m = kp1, la
         h = 0.0d0
         ortr(m) = 0.0d0
         orti(m) = 0.0d0
         scale = 0.0d0
c     .......... scale column (algol tol then not needed) ..........
         do 90 i = m, igh
   90    scale = scale + dabs(ar(i,m-1)) + dabs(ai(i,m-1))
c
         if (scale .eq. 0.0d0) go to 180
         mp = m + igh
c     .......... for i=igh step -1 until m do -- ..........
         do 100 ii = m, igh
            i = mp - ii
            ortr(i) = ar(i,m-1) / scale
            orti(i) = ai(i,m-1) / scale
            h = h + ortr(i) * ortr(i) + orti(i) * orti(i)
  100    continue
c
         g = dsqrt(h)
         f = pythag(ortr(m),orti(m))
         if (f .eq. 0.0d0) go to 103
         h = h + f * g
         g = g / f
         ortr(m) = (1.0d0 + g) * ortr(m)
         orti(m) = (1.0d0 + g) * orti(m)
         go to 105
c
  103    ortr(m) = g
         ar(m,m-1) = scale
c     .......... form (i-(u*ut)/h) * a ..........
  105    do 130 j = m, n
            fr = 0.0d0
            fi = 0.0d0
c     .......... for i=igh step -1 until m do -- ..........
            do 110 ii = m, igh
               i = mp - ii
               fr = fr + ortr(i) * ar(i,j) + orti(i) * ai(i,j)
               fi = fi + ortr(i) * ai(i,j) - orti(i) * ar(i,j)
  110       continue
c
            fr = fr / h
            fi = fi / h
c
            do 120 i = m, igh
               ar(i,j) = ar(i,j) - fr * ortr(i) + fi * orti(i)
               ai(i,j) = ai(i,j) - fr * orti(i) - fi * ortr(i)
  120       continue
c
  130    continue
c     .......... form (i-(u*ut)/h)*a*(i-(u*ut)/h) ..........
         do 160 i = 1, igh
            fr = 0.0d0
            fi = 0.0d0
c     .......... for j=igh step -1 until m do -- ..........
            do 140 jj = m, igh
               j = mp - jj
               fr = fr + ortr(j) * ar(i,j) - orti(j) * ai(i,j)
               fi = fi + ortr(j) * ai(i,j) + orti(j) * ar(i,j)
  140       continue
c
            fr = fr / h
            fi = fi / h
c
            do 150 j = m, igh
               ar(i,j) = ar(i,j) - fr * ortr(j) - fi * orti(j)
               ai(i,j) = ai(i,j) + fr * orti(j) - fi * ortr(j)
  150       continue
c
  160    continue
c
         ortr(m) = scale * ortr(m)
         orti(m) = scale * orti(m)
         ar(m,m-1) = -g * ar(m,m-1)
         ai(m,m-1) = -g * ai(m,m-1)
  180 continue
c
  200 return
      end
      subroutine csroot(xr,xi,yr,yi)
      double precision xr,xi,yr,yi
c
c     (yr,yi) = complex dsqrt(xr,xi) 
c     branch chosen so that yr .ge. 0.0 and sign(yi) .eq. sign(xi)
c
      double precision s,tr,ti,pythag
      tr = xr
      ti = xi
      s = dsqrt(0.5d0*(pythag(tr,ti) + dabs(tr)))
      if (tr .ge. 0.0d0) yr = s
      if (ti .lt. 0.0d0) s = -s
      if (tr .le. 0.0d0) yi = s
      if (tr .lt. 0.0d0) yr = 0.5d0*(ti/yi)
      if (tr .gt. 0.0d0) yi = 0.5d0*(ti/yr)
      return
      end
      double precision function pythag(a,b)
      double precision a,b
c
c     finds dsqrt(a**2+b**2) without overflow or destructive underflow
c
      double precision p,r,s,t,u
      p = dmax1(dabs(a),dabs(b))
      if (p .eq. 0.0d0) go to 20
      r = (dmin1(dabs(a),dabs(b))/p)**2
   10 continue
         t = 4.0d0 + r
         if (t .eq. 4.0d0) go to 20
         s = r/t
         u = 1.0d0 + 2.0d0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end
c ---------------------------
      subroutine mltpha4(katm,matm,pos,pai2,kngp,ngabc,kng1,nbase,iiba
     &     ,zmltc,zmlts)
c                           @(#)rmmsubs.f 1.2 99/06/04 21:17:36
      implicit real*8 (a-h,o-z)
      integer nbase(*),ngabc(kngp,3)
      real*8 pos(katm,3)
      real*8 zmltc(kng1,matm),zmlts(kng1,matm)
      do ia = 1, matm
         fx = pos(ia,1)*pai2
         fy = pos(ia,2)*pai2
         fz = pos(ia,3)*pai2
         do i = 1, iiba
            nb = nbase(i)
            ph = ngabc(nb,1)*fx + ngabc(nb,2)*fy + ngabc(nb,3)*fz
            zmltc(i,ia) = dcos(ph)
            zmlts(i,ia) = dsin(ph)
         enddo
      enddo
      return
      end

      subroutine decide_alpha(kimg,dtim,nrmm,eigenr,eigeni,vector,alpha)
      implicit none
      integer  kimg
      real*8   dtim
      integer  nrmm
      real*8   eigenr(0:nrmm-1), eigeni(0:nrmm-1)
      real*8   vector(0:nrmm-1,0:nrmm-1,kimg)
      real*8   alpha(nrmm-1,kimg)

      real*8 eee
      integer i, k

      do k = 1, kimg
         do i = 1, nrmm-1
            alpha(i,k) = 0.d0
         end do
      end do

      if(nrmm .eq. 2) then
         if(eigeni(0) .ne. 0.d0) then
c            write(6,*) 'Warning: Eigen value is not real.'
            alpha(1,1) = dtim   
         else if( min(eigenr(0),eigenr(1)) .lt. 0.0d0 ) then
c            write(6,*) 'Warning: Eigen value is less than zero.'
            alpha(1,1) = dtim   
         else if( eigenr(0) .lt. eigenr(1) ) then
            if(kimg .eq. 1) then
               alpha(1,1) = vector(1,0,1)/vector(0,0,1)
               if(alpha(1,1) .le. 0.0d0) alpha(1,1) = dtim
            else if(kimg .eq. 2) then
               alpha(1,1) = dcmplx(vector(1,0,1),vector(1,0,kimg))
     &              /dcmplx(vector(0,0,1),vector(0,0,kimg))
            end if
         else
            if(kimg .eq. 1) then
               alpha(1,1) = vector(1,1,1)/vector(0,1,1)
               if(alpha(1,1) .le. 0.0d0) alpha(1,1) = dtim
            else if(kimg .eq. 2) then
               alpha(1,1) = dcmplx(vector(1,1,1),vector(1,1,kimg))
     &               /dcmplx(vector(0,1,1),vector(0,1,kimg))
            end if
         end if
      else if(nrmm .eq. 3) then
         eee = abs(eigeni(0))+abs(eigeni(1))+abs(eigeni(nrmm-1))
c*
         if(eee.eq.0.0d0) then  ! All eigenvalues are real.
c*
            if     ( eigenr(0).le.min(eigenr(1),eigenr(2)) ) then
               if(kimg .eq. 1) then
                  alpha(1,1) = vector(1,0,1)/vector(0,0,1)
                  alpha(2,1) = vector(2,0,1)/vector(0,0,1)
               else if(kimg .eq. 2) then
                  alpha(1,1) = dcmplx(vector(1,0,1),vector(1,0,kimg))
     &                      /dcmplx(vector(0,0,1),vector(0,0,kimg))
                  alpha(2,1) = dcmplx(vector(2,0,1),vector(2,0,kimg))
     &                      /dcmplx(vector(0,0,1),vector(0,0,kimg))
               endif
            else if( eigenr(1).le.min(eigenr(2),eigenr(0)) ) then
               if(kimg .eq. 1) then
                  alpha(1,1) = vector(1,1,1)/vector(0,1,1)
                  alpha(2,1) = vector(2,1,1)/vector(0,1,1)
               else if(kimg .eq. 2) then
                  alpha(1,1) = dcmplx(vector(1,1,1),vector(1,1,kimg))
     &                      /dcmplx(vector(0,1,1),vector(0,1,kimg))
                  alpha(2,1) = dcmplx(vector(2,1,1),vector(2,1,kimg))
     &                      /dcmplx(vector(0,1,1),vector(0,1,kimg))
               endif
            else if( eigenr(2).le.min(eigenr(0),eigenr(1)) ) then
               if(kimg .eq. 1) then
                  alpha(1,1) = vector(1,2,1)/vector(0,2,1)
                  alpha(2,1) = vector(2,2,1)/vector(0,2,1)
               else
                  alpha(1,1) = dcmplx(vector(1,2,1),vector(1,2,kimg))
     &                      /dcmplx(vector(0,2,1),vector(0,2,kimg))
                  alpha(2,1) = dcmplx(vector(2,2,1),vector(2,2,kimg))
     &                      /dcmplx(vector(0,2,1),vector(0,2,kimg))
               endif
            else 
               alpha(1,1) = dtim
               alpha(2,1) = 0.0d0
            end if
c*
         else                   ! Just one eigenvalue is real.
c*
            if     ( eigeni(0).eq.0.0d0 ) then
               if(kimg .eq. 1) then
                  alpha(1,1) = vector(1,0,1)/vector(0,0,1)
                  alpha(2,1) = vector(2,0,1)/vector(0,0,1)
               else if(kimg .eq. 2) then
                  alpha(1,1) = dcmplx(vector(1,0,1),vector(1,0,kimg))
     &                      /dcmplx(vector(0,0,1),vector(0,0,kimg))
                  alpha(2,1) = dcmplx(vector(2,0,1),vector(2,0,kimg))
     &                      /dcmplx(vector(0,0,1),vector(0,0,kimg))
               endif
            else if( eigeni(1).eq.0.0d0 ) then
               if(kimg .eq. 1) then
                  alpha(1,1) = vector(1,1,1)/vector(0,1,1)
                  alpha(2,1) = vector(2,1,1)/vector(0,1,1)
               else
                  alpha(1,1) = dcmplx(vector(1,1,1),vector(1,1,kimg))
     &                      /dcmplx(vector(0,1,1),vector(0,1,kimg))
                  alpha(2,1) = dcmplx(vector(2,1,1),vector(2,1,kimg))
     &                      /dcmplx(vector(0,1,1),vector(0,1,kimg))
               endif
            else if( eigeni(2).eq.0.0d0 ) then
               if(kimg .eq. 1) then
                  alpha(1,1) = vector(1,2,1)/vector(0,2,1)
                  alpha(2,1) = vector(2,2,1)/vector(0,2,1)
               else if(kimg .eq. 2) then
                  alpha(1,1) = dcmplx(vector(1,2,1),vector(1,2,kimg))
     &                      /dcmplx(vector(0,2,1),vector(0,2,kimg))
                  alpha(2,1) = dcmplx(vector(2,2,1),vector(2,2,kimg))
     &                      /dcmplx(vector(0,2,1),vector(0,2,kimg))
               endif
            else
               alpha(1,1) = dtim
               alpha(2,1) = 0.0d0
            end if
         end if
      end if
      return
      end

      subroutine rmm3_uda(ipri,nrmm,am,sm,alpha)
      implicit real*8(a-h,o-z)
      integer ipri, nrmm
      dimension am(0:nrmm-1,0:nrmm-1),sm(0:nrmm-1,0:nrmm-1)
     &     ,alpha(nrmm-1)
      dimension c(3)
c         print *,'am=',am(0,0),am(1,0),am(2,0)
c         print *,'   ',am(1,1),am(2,1)
**********12/13******************
c      if(am(0,0)/sm(0,0).gt.1.0d-2) then
*      if(am(0,0).gt.1.0d0) then
c         print *,'calculation is made by SD'
c         print *,'am=',am(0,0),am(1,0),am(2,0)
c         print *,'   ',am(1,1),am(2,1)
c         print *,'   ',am(2,2)

c         print *,'sm=',sm(0,0),sm(1,0),sm(2,0)
c         print *,'   ',sm(1,1),sm(2,1)
c         print *,'   ',sm(2,2)

c         alpha(1)=-0.10d0
c         alpha(2)=0.0d0
c         return
c      end if
************12/13****************
      a11 = am(0,0)
      a12 = am(1,0)
      a13 = am(2,0)
      a22 = am(1,1)
      a23 = am(2,1)
      a33 = am(2,2)
c$$$      a11=am(1)
c$$$      a12=am(2)
c$$$      a13=am(3)
c$$$      a22=am(4)
c$$$      a23=am(5)
c$$$      a33=am(6)

      s11 = sm(0,0)
      s12 = sm(1,0)
      s13 = sm(2,0)
      s22 = sm(1,1)
      s23 = sm(2,1)
      s33 = sm(2,2)
c$$$      s11=sm(1)
c$$$      s12=sm(2)
c$$$      s13=sm(3)
c$$$      s22=sm(4)
c$$$      s23=sm(5)
c$$$      s33=sm(6)

c      print *,'-------------------------------------------------------'
c      print *,'a11=',a11
c      print *,'a12=',a12
c      print *,'a13=',a13
c      print *,'a22=',a22
c      print *,'a23=',a23
c      print *,'a33=',a33
c      print *,'s11=',s11
c      print *,'s12=',s12
c      print *,'s13=',s13
c      print *,'s22=',s22
c      print *,'s23=',s23
c      print *,'s33=',s33
c      print *,'a11/s11=',a11/s11
c      print *,'a22/s22=',a22/s22
c      print *,'a33/s33=',a33/s33
c

      a21=a12
      a32=a23   
      a31=a13   
      s21=s12
      s32=s23
      s31=s13

c$$$      write(6,'(" -- s(<i|S|j>) (rmm3_uda) --")')
c$$$      write(6,'(3d12.4, " : s11 s12 s13 ")') s11,s12,s13
c$$$      write(6,'(3d12.4, " : s21 s22 s23 ")') s21,s22,s23
c$$$      write(6,'(3d12.4, " : s31 s32 s33 ")') s31,s32,s33
c$$$      write(6,'(" -- a(rr=<Ri|Rj>) (rmm3_uda) --")')
c$$$      write(6,'(3d12.4, " : r11 r12 r13 ")') a11,a12,a13
c$$$      write(6,'(3d12.4, " : r21 r22 r23 ")') a21,a22,a23
c$$$      write(6,'(3d12.4, " : r31 r32 r33 ")') a31,a32,a33
      
      r10=a22*a33-a32*a23
      r20=a21*a33-a31*a23
      r30=a21*a32-a31*a22

      r11=a32*s23+a23*s32 - a22*s33-s22*a33
      r21=a31*s23+s31*a23 - a21*s33-s21*a33
      r31=a31*s22+s31*a22 - a21*s32-s21*a32

      r12=s22*s33 - s32*s23
      r22=s21*s33 - s31*s23
      r32=s21*s32 - s31*s22

      f0= a11*r10 - a12*r20 + a13*r30
      f1= a11*r11 - a12*r21 + a13*r31
     -    -s11*r10 + s12*r20 - s13*r30

      f2= a11*r12 - a12*r22 + a13*r32
     -    -s11*r11 + s12*r21 - s13*r31

      f3= -s11*r12 + s12*r22 - s13*r32
c   -------------- 20021223 -------------------- 
        f3= f3/f0
        f2= f2/f0
        f1= f1/f0
        f0= 1.0d0
c ===============================================
   
      e1=0.d0
      e2=a11/s11
c      print *,'e2=',e2
      e0=0.d0
c        s=e2*(e2*(f3*e2+f2)+f1)+f0 
c        if(s.le.0.d0) then
c           e3=e2
c           e2=0.5d0*e2
c           else 
c            go to 50
c            end if
      do 40 k=1,100
         s=e2*(e2*(f3*e2+f2)+f1)+f0 
        if(s.le.0.d0) then
           go to 50
        else
         s=e2*(3.0d0*f3*e2+2.0d0*f2)+f1 
c           et=e2+1.d-30
c           s1=et*(et*(f3*et+f2)+f1)+f0 
c           if((s1-s).ge.0.d0 ) then
           if(s.ge.0.d0 ) then
              e2=e0 + (e2-e0)*0.5d0 
           else
              e3=e2
              e2 = e2 + (e2-e0)*0.5d0
              e0 = e2
           end if
        end if
 40   continue
      if(ipri .ge. 1) then
      write(6,*)'******************************************************'
      write(6,*)'*  warning from sub rmm3    k=100 is insufficient    *'
      write(6,*)'******************************************************'
      end if
 50   continue
c       write(6,*)'k=',k
c       write(6,*)'          ---------------------------------'
	e1=e0

      do 10 j=1,250

        s=e2*(e2*(f3*e2+f2)+f1)+f0 
c$$$        if((e2-e1).lt.5.d-22) go to 20
        if((e2-e1)/e2.lt.1.d-10) go to 20
c        if((dabs(s).lt.1.d-90).or.((e2-e1)/e2.lt.1.d-17)) go to 20
c        if((qabs(s).lt.1.d-70).or.((e2-e1)/e2.lt.1.d-10)) go to 20

        if(s.lt.0.d0) then
           e2=e1+(e2-e1)*0.5d0
        else if (s. gt. 0.d0)then
           e3=e2
           e2=e2+(e2-e1)*0.5d0
           e1=e3
        else 
           go to 20
        end if
 10   continue
      if(ipri .ge. 1) then
      write(6,*)'******************************************************'
      write(6,*)'*  warning from sub rmm3    j=250 is insufficient    *'
      write(6,*) 'e2,e1=',e2,e1
      write(6,*)'******************************************************'
      end if
 20   continue
      if((e2/a11)*s11.ge.1.d0) then
      if(ipri .ge. 1 ) then
      write(6,*)'******************************************************'
      write(6,*)'*  warning from sub rmm3        vector not improved  *'
      write(6,*) 'a11/s11=',a11/s11
      write(6,*) 'e2,s=',e2,s
      write(6,*) 'e1=',e1
      write(6,*) 'j=',j
      write(6,*)'******************************************************'
      end if
      end if
c      write(6,*) 'a11=',a11
c      write(6,*) 'e2=',e2

c      e3=(e1+e2)*0.5d0
      c(1)=1.d0
      d=(a12-e2*s12)*(a23-e2*s23)- (a13-e2*s13)*(a22-e2*s22)
      c(2)=(-(a11-e2*s11)*(a23-e2*s23)+(a13-e2*s13)*(a21-e2*s21))/d
      c(3)=(-(a12-e2*s12)*(a21-e2*s21)+(a22-e2*s22)*(a11-e2*s11))/d

      ee=(c(1)*c(1)*a11+2.d0*c(1)*c(2)*a12+c(2)*c(2)*a22+
     -     2.d0*c(1)*c(3)*a13+2.d0*c(2)*c(3)*a23+c(3)*c(3)*a33)/
     -     (c(1)*c(1)*s11+2.d0*c(1)*c(2)*s12+c(2)*c(2)*s22+
     -     2.d0*c(1)*c(3)*s13+ 2.d0*c(2)*c(3)*s23+c(3)*c(3)*s33)

      if((dabs(ee-e2)/e2.ge.1.d-5).and.(e2.gt.1.d0-20)) then
c      if((dabs(ee-e2)/e2.ge.1.d-3).and.(e2.gt.1.d0-13)) then
c      if(qabs(ee-e2)/e2.ge.1.d-5) then
      if(ipri .ge. 1) then
      write(6,*)'******************************************************'
      write(6,*)'*  warning from sub rmm3        vector not reliable  *'
      write(6,*) 'ee=',ee
      write(6,*) 'e1=',e1
      write(6,*) 'e2=',e2
      write(6,*)'******************************************************'

cccccc12/12/2002cccccccccccccccc
      write(6,*) ' f3,f2=',f3,f2
      write(6,*) ' f1,f0=',f1,f0
      write(6,*) 'a(0,0),(0,1),(0,2)=',am(0,0),am(0,1),am(0,2)
      write(6,*) '      a(1,1),(1,2)=',am(1,1),am(1,2)
      write(6,*) '            a(2,2)=',am(2,2)
      write(6,*) 's(0,0),(0,1),(0,2)=',sm(0,0),sm(0,1),sm(0,2)
      write(6,*) '      s(1,1),(1,2)=',sm(1,1),sm(1,2)
      write(6,*) '            s(2,2)=',sm(2,2)
      end if
      end if
c      write(6,*) 'ee=',ee
c      write(6,*) 'e2=',e2

c      write(6,*)'   e1=',e1
c      write(6,*)'   e2=',e2
c      write(6,*)'   e3=',e3
c      write(6,*)' jj=',j,'   f(e1) = ',s
c       write(6,*)'eigen vector =',c(1),c(2),c(3)
c      write(6,*) 'recalculated eigen value by using c1,c2,c3'
c      write(6,*) ee
 
 30   continue
      alpha(1) = c(2)
      alpha(2) = c(3)
c      alpha(1) = c(2)*0.75d0
c      alpha(2) = c(3)*0.75d0
      return
      end

      subroutine rmm2_uda(ipri,nrmm,am,sm,alpha)
      implicit real*8(a-h,o-z)
      integer  ipri, nrmm
      dimension am(0:nrmm-1,0:nrmm-1),sm(0:nrmm-1,0:nrmm-1)
     &     ,alpha(nrmm-1)
      dimension c(2)
c-----------------------------------------------------------------
c   
c     Residual Minimization Method with 2 basis
c  
c           R=H-<H>   
c   
c           Minimization of Residual  <i|R*R|j>/ <i|S|j> (i,j=1,2)
c 
c                      1>=|0> : trial vector
c                      2>=R|0>
c           
c           Initial Residual: Ri= <1|R*R|1>/ <1|S|1> 
c     
c                When Ri<1.d-20 this routine is skipped.      
c   
c                                          2003/02/21
c                                          written by T.Uda   
c    
c-----------------------------------------------------------------
      a11 = am(0,0)
      a12 = am(1,0)
      a22 = am(1,1)
      s11 = sm(0,0)
      s12 = sm(1,0)
      s22 = sm(1,1)

c      print *,'-------------------------------------------------------'
c      print *,'a11=',a11
c      print *,'a12=',a12
c      print *,'a22=',a22
c      print *,'s11=',s11
c      print *,'s12=',s12
c      print *,'s22=',s22

      a21=a12
      s21=s12

      f1= s11*s22 - s12*s21

      f2= a12*s21 + a21*s12 
     -     -s11*a22 - a11*s22

      f3= a11*a22  - a21*a12
    
      e1=0.d0
      e2=a11/s11
c$$$      e2 = -f2/(2*f1)
      if (e2.lt.1.d-20)then
         alpha(1)=0.0d0
         return
         end if

c      do 40 k=1,100
c           go to 50
c        else
c           et=e2+1.d-30
c           s1=et*(et*f1+f2)+f3 
c           if((s1-s).ge.0.d0 ) then
c              e2=e0 + (e2-e0)*0.5d0 
c           else
c              e3=e2
c              e2 = e2 + (e2-e0)*0.5d0
c              e0 = e2
c           end if
c        end if
        s=e2*(e2*f1+f2)+f3 
        if(s.gt.0.d0) then
           if(ipri .ge. 1) then
      write(6,*)'******************************************************'
      write(6,*)'*  warning from sub rmm2  vector is not improved    *'
      write(6,*)'******************************************************'
           end if
        alpha(1)=0.0d0
        return
        end if
c 50   continue
c       write(6,*)'k=',k
c       write(6,*)'          ---------------------------------'

      e1start = e1
      e2start = e2
      do 10 j=1,250

         s=e2*(e2*f1+f2)+f3
c$$$         if((e2-e1).lt.5.d-22) go to 20
c$$$        if((qabs(s).lt.1.d-70).or.((e2-e1)/e2.lt.1.d-10)) go to 20
        if((e2-e1)/e2.lt.1.d-10) go to 20

        if(s.lt.0.d0) then
           e2=e1+(e2-e1)*0.5d0
        else if (s. gt. 0.d0)then
           e3=e2
           e2=e2+(e2-e1)*0.5d0
           e1=e3
        else 
           go to 20
        end if

 10   continue
      if(ipri .ge. 1) then
      write(6,*)'******************************************************'
      write(6,*)'*  warning from sub rmm2    j=250 is insufficient    *'
      write(6,*)'******************************************************'
      write(6,'(" s = ",d15.5, " e2,e1,e2-e1 = ",3d15.5)') s,e2,e1,e2-e1
      write(6,'(" e2start, e1start, e2start-e1start = ",3d15.5)')
     &     e2start,e1start,e2start-e1start
      end if
 20   continue

      if((e2/a11)*s11.ge.1.d0) then
         if( ipri. ge. 1) then
      write(6,*)'******************************************************'
      write(6,*)'*  warning from sub rmm2        vector not improved  *'
      write(6,*)'******************************************************'
      write(6,*) 'a11/s11=',a11/s11
      write(6,*) 'e2=',e2
         end if
      end if
c      write(6,*) 'a11=',a11
c      write(6,*) 'e2=',e2

c      e3=(e1+e2)*0.5d0
      c(1)=1.d0
      c(2)=-(a11-e2*s11)/(a12-e2*s12)

      ee=(c(1)*c(1)*a11+2.d0*c(1)*c(2)*a12+c(2)*c(2)*a22)/
     -   (c(1)*c(1)*s11+2.d0*c(1)*c(2)*s12+c(2)*c(2)*s22)

      if(dabs(ee-e2)/e2.ge.1.d-5) then
c      if(qabs(ee-e2)/e2.ge.1.d-5) then
         if(ipri .ge. 1) then
      write(6,*)'******************************************************'
      write(6,*)'*  warning from sub rmm2        vector not reliable  *'
      write(6,*)'******************************************************'
      write(6,*) 'ee=',ee
      write(6,*) 'e2=',e2
         end if
      end if
c      write(6,*) 'ee=',ee
c      write(6,*) 'e2=',e2

c      write(6,*)'   e1=',e1
c      write(6,*)'   e2=',e2
c      write(6,*)'   e3=',e3
c      write(6,*)' jj=',j,'   f(e1) = ',s
c       write(6,*)'eigen vector =',c(1),c(2)
c      write(6,*) 'recalculated eigen value by using c1,c2'
c      write(6,*) ee
 
 30   continue

      alpha(1) = c(2)
      return
      end
c-----------------------------------------------------------------
      subroutine crmm2_uda(ipri,nrmm,rcam,rcsm,alpha)
      implicit real*8(a,b,d-h,o-z)
      implicit complex*16 (c)
      integer ipri, nrmm
      dimension rcam(0:nrmm-1,0:nrmm-1,2),rcsm(0:nrmm-1,0:nrmm-1,2)
     &     ,alpha(nrmm-1,2)
      complex*16 cam(0:1,0:1),csm(0:1,0:1)
      dimension c(2)
c-----------------------------------------------------------------
c   
c     Residual Minimization Method with 2 basis
c  
c           R=H-<H>   
c   
c           Minimization of Residual  <i|R*R|j>/ <i|S|j> (i,j=1,2)
c 
c                      1>=|0> : trial vector
c                      2>=R|0>
c           
c           Initial Residual: Ri= <1|R*R|1>/ <1|S|1> 
c     
c                When Ri<1.d-20 this rootine is skipped.      
c   
c                                          2003/02/21
c                                          written by T.Uda   
c    
c    * Revised by T. Yamasaki,  14th Mar. 2003 
c        + calpha
c               --> alpha
c        + calpha(1) = c(2)
c               -->
c             alpha(1) = dreal(c(2))
c    * Revised by T. Uda and T. Yamasaki, 19th Mar. 2003
c        + alpha(nrmm-1)
c               --> alpha(nrmm-1,2)
c        + alpha(1) = dreal(c(2))
c               -->
c             alpha(1,1) = dreal(c(2))
c             alpha(1,2) = dimag(c(2))
c-----------------------------------------------------------------


      do i = 0, nrmm-1
         do j = 0, nrmm-1
            cam(i,j) = dcmplx(rcam(i,j,1),rcam(i,j,2))
            csm(i,j) = dcmplx(rcsm(i,j,1),rcam(i,j,2))
         end do
      end do
      ca11 = cam(0,0)
      ca21 = cam(1,0)
      ca12 = cam(0,1)
      ca22 = cam(1,1)

      cs11 = csm(0,0)
      cs21 = csm(1,0)
      cs12 = csm(0,1)
      cs22 = csm(1,1)

c      print *,'-------------------------------------------------------'
c      print *,'ca11=',ca11
c      print *,'ca12=',ca12
c      print *,'ca22=',ca22
c      print *,'cs11=',cs11
c      print *,'cs12=',cs12
c      print *,'cs22=',cs22


      f1= cs11*cs22 - cs12*cs21

      f2= ca12*cs21 + ca21*cs12 
     -     -cs11*ca22 - ca11*cs22

      f3= ca11*ca22  - ca21*ca12



      e1=0.d0
      e2=ca11/cs11
c$$$      e20=e2

      if (e2.lt.1.d-20) then
         alpha(1,1)=0.0d0
         alpha(1,2)=0.0d0
         return
      end if

c      do 40 k=1,100
c        if(s.le.0.d0) then
c           go to 50
c        else
c           s=2.0d0*e2*f1+f2
c           et=e2+1.d-30
c           s1=et*(et*f1+f2)+f3

        s=e2*(e2*f1+f2)+f3

           if(s.ge.0.d0 ) then
              if(ipri .ge. 1) then
      write(6,*)'******************************************************'
      write(6,*)'*  warning from sub crmm2   vector is not improved   *'
      write(6,*)'******************************************************'
              end if
           
           alpha(1,1)=0.0d0
           alpha(1,2)=0.0d0
           return
           end if

c           else
c              e3=e2
c              e2 = e2 + (e2-e0)*0.5d0
c              e0 = e2
c           end if
c        end if
c 40   continue
c 50   continue
c       write(6,*)'k=',k
c       write(6,*)'          ---------------------------------'

      do 10 j=1,250

        s=e2*(e2*f1+f2)+f3
c$$$        if((e2-e1).lt.5.d-22) go to 20
c        if((dabs(s).lt.1.d-90).or.((e2-e1)/e2.lt.1.d-10)) go to 20
        if((e2-e1)/e2.lt.1.d-10) go to 20

        if(s. lt. 0.d0) then
           e2=e1+(e2-e1)*0.5d0
        else if (s. gt. 0.d0)then
           e3=e2
           e2=e2+(e2-e1)*0.5d0
           e1=e3
        else 
           go to 20
        end if
 10   continue
      if(ipri .ge. 1) then
      write(6,*)'******************************************************'
      write(6,*)'*  warning from sub crmm2    j=250 is insufficient   *'
      write(6,*)'******************************************************'
      write(6,'(" s = ",d20.5)') s
      end if
 20   continue
      
      test = (e2/ca11)*cs11
      if(test.ge.1.d0) then
         if(ipri .ge. 1) then
      write(6,*)'******************************************************'
      write(6,*)'*  warning from sub crmm2        vector not improved *'
      write(6,*)'******************************************************'
      write(6,*) 'a11/s11=',ca11/cs11
      write(6,*) 'e2=',e2
         end if
      end if
c      write(6,*) 'a11=',ca11
c      write(6,*) 'e2=',e2

c      e3=(e1+e2)*0.5d0
      c(1)=1.d0
      c(2)=-(ca11-e2*cs11)/(ca12-e2*cs12)

      ee=real(c(1)*c(1)*ca11+2.d0*c(1)*c(2)*ca12+
     -                                     dconjg(c(2))*c(2)*ca22)/
     -   real(c(1)*c(1)*cs11+2.d0*c(1)*c(2)*cs12+
     -                                     dconjg(c(2))*c(2)*cs22)

      if(dabs(ee-e2)/e2.ge.1.d-5) then

         if(ipri .ge. 1) then
      write(6,*)'******************************************************'
      write(6,*)'*  warning from sub crmm2       vector not reliable  *'
      write(6,*)'******************************************************'
      write(6,*) 'ee=',ee
      write(6,*) 'e1,e2=',e1,e2
         end if
      end if
 
 30   continue

      alpha(1,1) = dreal(c(2))
      alpha(1,2) = dimag(c(2))
      return
      end
c------------------------------------------------------------------
      subroutine crmm3_uda(ipri,nrmm,rcam,rcsm,alpha)
      implicit real*8(a,b,d-h,o-z)
      implicit complex*16 (c)
      integer ipri,nrmm
      dimension rcam(0:nrmm-1,0:nrmm-1,2),rcsm(0:nrmm-1,0:nrmm-1,2)
     &     ,alpha(nrmm-1,2)
      complex*16 cam(0:2,0:2),csm(0:2,0:2),c(3)
c-----------------------------------------------------------------
c   
c     Residual Minimization Method with 3basis
c  
c           R=H-<H>   
c   
c           Minimization of Residual  <i|R*R|j>/ <i|S|j> (i,j=1,2,3)
c 
c                      |1>=|0> : trial vector
c                      |2>=R|0>
c                      |3>=R|2>=R*R|0>
c           
c           Initial Residual: Ri= <1|R*R|1>/ <1|S|1> 
c     
c                When Ri<1.d-20 this rootine is skipped.      
c   
c                                          2003/02/21
c                                          written by T.Uda   
c
c    * Revised by T. Yamasaki,  14th Mar. 2003 
c        + calpha
c               --> alpha
c        + calpha(1) = c(2)
c          calpha(2) = c(3)
c               -->
c             alpha(1) = dreal(c(2))
c             alpha(2) = dreal(c(3))
c    * Revised by T. Uda and T. Yamasaki, 19th Mar. 2003
c        + alpha(nrmm-1)
c               --> alpha(nrmm-1,2)
c        + alpha(1) = dreal(c(2))
c          alpha(2) = dreal(c(3))
c               -->
c             alpha(1,1) = dreal(c(2))
c             alpha(2,1) = dreal(c(3))
c             alpha(1,2) = dimag(c(2))
c             alpha(2,2) = dimag(c(3))
c-----------------------------------------------------------------

c$$$      write(6,'(" <<crmm3_uda>>")')

      do i = 0, nrmm-1
         do j = 0, nrmm-1
            cam(i,j) = dcmplx(rcam(i,j,1),rcam(i,j,2))
            csm(i,j) = dcmplx(rcsm(i,j,1),rcam(i,j,2))
         end do
      end do
c
      ca11 = cam(0,0)
      ca12 = cam(0,1)
      ca13 = cam(0,2)
      ca22 = cam(1,1)
      ca23 = cam(1,2)
      ca33 = cam(2,2)
      ca21 = cam(1,0)
      ca31 = cam(2,0)
      ca32 = cam(2,1)

      cs11 = csm(0,0)
      cs12 = csm(0,1)
      cs13 = csm(0,2)
      cs22 = csm(1,1)
      cs23 = csm(1,2)
      cs33 = csm(2,2)
      cs21 = csm(1,0)
      cs31 = csm(2,0)
      cs32 = csm(2,1)
      

c$$$      ex1=ca11/cs11
c$$$      ex2=ca22/cs22
c$$$      ex3=ca33/cs33

c      print *,'ex1=',ex1
c      print *,'ex2=',ex2
c      print *,'ex3=',ex3
      
c      isign=0
c      if (ex3.lt.ex2) then
c       ctest=ca12
c       ca12 =ca13
c       ca13 =ctest



c$$$      print *,'-------------------------------------------------------'
c$$$      print *,'ca11=',ca11
c$$$      print *,'ca12=',ca12
c$$$      print *,'ca21=',ca21
c$$$      print *,'ca13=',ca13
c$$$      print *,'ca31=',ca31
c$$$      print *,'ca22=',ca22
c$$$      print *,'ca23=',ca23
c$$$      print *,'ca32=',ca32
c$$$      print *,'ca33=',ca33

c      print *,'cs11=',cs11
c      print *,'cs12=',cs12
c      print *,'cs21=',cs21
c      print *,'cs13=',cs13
c      print *,'cs31=',cs31
c      print *,'cs22=',cs22
c      print *,'cs23=',cs23
c      print *,'cs32=',cs32
c      print *,'cs33=',cs33

      f1= cs11*cs22 - cs12*cs21

      f2= ca12*cs21 + ca21*cs12 
     -     -cs11*ca22 - ca11*cs22

      f3= ca11*ca22  - ca21*ca12



      e1=0.d0
      e2=ca11/cs11
      e20=e2

      if (e2.lt.1.d-20) then
         goto 1001
c$$$         alpha(1,1)=0.0d0
c$$$         alpha(2,1)=0.0d0
c$$$         alpha(1,2)=0.0d0
c$$$         alpha(2,2)=0.0d0
c$$$         return
      end if

        s=e2*(e2*f1+f2)+f3

           if(s.ge.0.d0 ) then
              if(ipri .ge. 1) then
      write(6,*)'******************************************************'
      write(6,*)'*  warning from sub crmm3   vector is not improved   *'
      write(6,*)'******************************************************'
              end if
           
              goto 1001

           return
           end if


      do 10 j1=1,250
c$$$         write(6,'(" j1 = ",i5)') j1

        s=e2*(e2*f1+f2)+f3
c$$$        if((e2-e1).lt.5.d-22) go to 20
        if((e2-e1)/e2.lt.1.d-10) go to 20

        if(s.lt.0.d0) then
           e2=e1+(e2-e1)*0.5d0
        else if (s. gt. 0.d0)then
           e3=e2
           e2=e2+(e2-e1)*0.5d0
           e1=e3
        else 
           go to 20
        end if
 10   continue
      if(ipri .ge. 1) then
      write(6,*)'******************************************************'
      write(6,*)'*  warning from sub crmm3   j1=250 is insufficient   *'
      write(6,*)' e20,e21=',e20,e2
      write(6,*)' e1 ,e2 =',e1,e2
      write(6,*)'******************************************************'
      end if
 20   continue

c$$$      write(6,'(" !! out of 20 continue")')

      cr10=ca22*ca33-ca32*ca23
      cr20=ca21*ca33-ca31*ca23
      cr30=ca21*ca32-ca31*ca22

      cr11=ca32*cs23+ca23*cs32 - ca22*cs33-cs22*ca33
      cr21=ca31*cs23+cs31*ca23 - ca21*cs33-cs21*ca33
      cr31=ca31*cs22+cs31*ca22 - ca21*cs32-cs21*ca32

      cr12=cs22*cs33 - cs32*cs23
      cr22=cs21*cs33 - cs31*cs23
      cr32=cs21*cs32 - cs31*cs22

      f0= ca11*cr10 - ca12*cr20 + ca13*cr30
      f1= ca11*cr11 - ca12*cr21 + ca13*cr31
     -    -cs11*cr10 + cs12*cr20 - cs13*cr30

      f2= ca11*cr12 - ca12*cr22 + ca13*cr32
     -    -cs11*cr11 + cs12*cr21 - cs13*cr31

      f3= -cs11*cr12 + cs12*cr22 - cs13*cr32
c --------------------------------20030219---------------------
             f3= f3/f0
             f2= f2/f0
             f1= f1/f0
             f0= 1.0d0
c =============================================================
c$$$      do j1 = 1, 50
c$$$         write(6,'(" !! out of devided by f0: f0 = ",d20.10)') f0
c$$$      end do
    
c      e1=0.d0
c      e2=ca11/cs11
c      e20=ca11/cs11
c      print *,'e2=',e2
c      e0=0.d0
c        s=e2*(e2*(f3*e2+f2)+f1)+f0 
c        if(s.le.0.d0) then
c           e3=e2
c           e2=0.5d0*e2
c           else 
c            go to 50
c            end if
c      do 40 k=1,100
c        s=e2*(e2*(f3*e2+f2)+f1)+f0 
c        if(s.le.0.d0) then
c           go to 50
c        else
c           s=e2*(3.0d0*f3*e2+f2)+f1
c           et=e2+1.d-30
c           s1=et*(et*(f3*et+f2)+f1)+f0 
c           if(s.ge.0.d0 ) then
c              e2=e0 + (e2-e0)*0.5d0 
c           else
c              e3=e2
c              e2 = e2 + (e2-e0)*0.5d0
c              e0 = e2
c           end if
c        end if
c 40   continue
c      write(6,*)'******************************************************'
c      write(6,*)'*  warning from sub crmm3    k=100 is insufficient   *'
c      write(6,*)'******************************************************'
c 50   continue
c       write(6,*)'k=',k
c       write(6,*)'          ---------------------------------'

      e1=0.0d0
      e21=e2

 
      do 110 j2=1,250
      
        s=e2*(e2*(f3*e2+f2)+f1)+f0 
c$$$        write(6,'(" j2 = ",i5, " e2 = ",d20.10," s= ",d20.10)') j2, e2,s
c$$$        if((e2-e1).lt.1.0d-25) go to 120
        if((e2-e1)/e2.lt.1.0d-10) go to 120

        if(s.lt.0.d0) then
           e2=e1+(e2-e1)*0.5d0
        else if (s. gt. 0.d0)then
           e3=e2
           e2=e2+(e2-e1)*0.5d0
           e1=e3
        else 
           go to 120
        end if
 110   continue
       if(ipri .ge. 1) then
      write(6,*)'******************************************************'
      write(6,*)'*  warning from sub crmm3  j2=250 is insufficient    *'
      write(6,*)' e20,e21,e2=',e20,e21,e2
      write(6,*)' e1 ,e2 =',e1,e2
      write(6,*) 'j =',j
      write(6,*)'******************************************************'
       end if
 120   continue

      d=f2*f2 - 3.d0*f1*f3
      if (d.lt.0.0d0 .or. dabs(f3).lt.1.d-50) then
         if(ipri .ge. 1)
     &   write(6,'(" !rmmsubs d = ",d20.10," f1,f2,f3 = ",3d20.10)')
     &        d, f1,f2,f3
         if(d.lt.0.0d0 .and. ipri.ge.1) 
     &        write(6,*) ' !crmm3 warning d < 0.0'
         if(dabs(f3).lt.1.d-50 .and. ipri.ge.1) 
     &        write(6,*) ' !crmm3 warning |f3| <delta'
         goto 1001
      end if

c$$$       zeropoint1 = (-f2-sqrt(d))/(3.0d0*f3)
       zeropoint2 = (-f2+sqrt(d))/(3.0d0*f3)

c$$$      write(6,'(" !rmmsubs zeropoint2 = ",d20.10)') zeropoint2

c$$$       s0=zeropoint2*(zeropoint2*(zeropoint2*f3+f2)+f1)+f0
       s1=e1*(e1*(e1*f3+f2)+f1)+f0
       s2=e2*(e2*(e2*f3+f2)+f1)+f0
c$$$       do j1 = 1, 200
c$$$          write(6,'(" !!s1,s2,zeropoint2 = ",3d20.10)') s1,s2,zeropoint2
c$$$       end do
c       write(6,*) 'zeropoint2, f(zeropoint2)=',zeropoint2,s2
c       write(6,*)' e20,e21,e2=',e20,e21,e2

      if(e2.gt.zeropoint2) then
         if(ipri .ge. 1) then
      write(6,*)'******************************************************'
      write(6,*)'*  warning from sub crmm3      e2 > zero point       *'
      write(6,*)'   e2,e21,e20=',e2,e21,e20
      write(6,*)'******************************************************'
         end if
      end if
      if(s1.lt.-1.0d-10.or.s2.gt.1.0d-10) then
         if(ipri .ge. 1) then
      write(6,*)'******************************************************'
      write(6,*)'*  warning from sub crmm3                            *'
      write(6,*)'   e1,e2,zeropoint2=',e1,e2,zeropoint2
      write(6,*)'   e2,e21,e20=',e1,e21,e20
      write(6,*)' f(e1)=',s1,'f(e2)=',s2
      write(6,*)'******************************************************'
         end if
      end if
 

c      e3=(e1+e2)*0.5d0
      c(1)=1.d0
      cd=(ca12-e2*cs12)*(ca23-e2*cs23)- (ca13-e2*cs13)*(ca22-e2*cs22)
      c(2)=(-(ca11-e2*cs11)*(ca23-e2*cs23)+(ca13-e2*cs13)*
     -    (ca21-e2*cs21))/cd
      c(3)=(-(ca12-e2*cs12)*(ca21-e2*cs21)+(ca22-e2*cs22)*
     -    (ca11-e2*cs11))/cd

      ee=real(c(1)*c(1)*ca11+2.d0*c(1)*c(2)*ca12
     -                      +conjg(c(2))*c(2)*ca22+
     -     2.d0*c(1)*c(3)*ca13+2.d0*conjg(c(2))*c(3)*ca23
     -                      +conjg(c(3))*c(3)*ca33)/
     -     real(c(1)*c(1)*cs11+2.d0*c(1)*c(2)*cs12
     -                      +conjg(c(2))*c(2)*cs22+
     -     2.d0*c(1)*c(3)*cs13+ 2.d0*conjg(c(2))*c(3)*cs23
     -                      +conjg(c(3))*c(3)*cs33)

      if((dabs(ee-e2)/e2).ge.1.d-5) then
         if(ipri .ge. 1) then
      write(6,*)'******************************************************'
      write(6,*)'*  warning from sub crmm3        vector not reliable *'
      write(6,*) 'e20,ee=',e20,ee
      write(6,*) 'e1 ,e2=',e1,e2
      write(6,*)'******************************************************'
         end if
      end if

 30   continue
      alpha(1,1) = dreal(c(2))
      alpha(2,1) = dreal(c(3))
      alpha(1,2) = dimag(c(2))
      alpha(2,2) = dimag(c(3))
      return
 1001 continue
      alpha(1,1) = 0.d0
      alpha(2,1) = 0.d0
      alpha(1,2) = 0.d0
      alpha(2,2) = 0.d0
      return
      end

