!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: fft_c3d fft_r3d fft_trans_yz fft_trans_zy mkd_c3d
!
!  AUTHOR(S): T. Yamasaki, T. Yamamoto   June/01/2007
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


!     ****************************************************************
!     AUTHORS: T.KOKUBO(NEC), K.TAKADA(NIS), T.FURUTA(NIS), R.SAKA(NIS)
!              Y.SAKAGUCHI(NES)
!     REGD. ON 03/15/02 BY T.F.
!     UPGD. ON 04/18/06 BY T.K. & Y.S.
!     UPGD. ON 09/08/06 BY T.K.
!
      SUBROUTINE FFT_C3D(NX,NY,NZ,NY_D,NZ_D,A,lx,ly,lz,ly_d,lz_d,
     &                   INDEX,ISW,TIME,B,NT,FFT_WORLD,IERR)
!
!     IMPLICIT NONE
#ifdef MKLFFT
      USE MKL_DFTI
#endif
      INCLUDE "mpif.h"
!
      INTEGER,PARAMETER :: PREC=8
      INTEGER,INTENT(IN) :: FFT_WORLD,NX,NY,NZ,ISW,NT
      INTEGER,INTENT(IN) :: INDEX,NY_D,NZ_D
      INTEGER,INTENT(IN) :: lx,ly,ly_d,lz,lz_d
      INTEGER,INTENT(INOUT) :: IERR(2)
      COMPLEX(KIND=PREC),INTENT(INOUT) :: A(0:*),B(0:*)
      REAL(KIND=PREC),INTENT(OUT) :: TIME(5)
!
      INTEGER :: IRANK,ISIZE
#ifdef ASLFFT
      INTEGER :: IFAX(20)
#elif SRFFT
      INTEGER :: ip(4),iopt(3),n,nwkw
      INTEGER,PARAMETER :: MAXFAX = 4
      INTEGER,DIMENSION(MAXFAX) :: NFAX = (/2,3,5,7/)
#elif FFTW3
      INTEGER(KIND=PREC) :: PLAN
      INTEGER,DIMENSION(1) :: NFFTW3
      INTEGER,DIMENSION(1) :: NEMBED
      INTEGER :: ISWFFTW
      INTEGER,PARAMETER :: FFTW_RANK = 1
      INTEGER,PARAMETER :: FFTW_ESTIMATE = 64
#elif MKLFFT
      INTEGER :: STRIDE(2), Status
      TYPE(DFTI_DESCRIPTOR), POINTER :: Desc_Handle
#endif
      INTEGER :: I,J,K,MAXSIZE,KERR,ll,ll_d
      REAL(KIND=PREC) :: T0,T1,T2
      COMPLEX(KIND=PREC),DIMENSION(:),ALLOCATABLE :: TRIGS
!
      KERR = 0
      IERR = 0
#ifdef FFTW3
#ifdef _OPENMP
      if(NT>1) then
         call dfftw_init_threads
         call dfftw_plan_with_nthreads(NT)
      end if
#endif
      ISWFFTW   = -ISW
#endif
!
! *** MPI INITIALIZE ***
      CALL MPI_COMM_RANK(FFT_WORLD,IRANK,IERR(2))
      IF( IERR(2) .NE. 0 ) THEN
         IERR(1) = 8000
         RETURN
      ENDIF
      CALL MPI_COMM_SIZE(FFT_WORLD,ISIZE,IERR(2))
      IF( IERR(2) .NE. 0 ) THEN
         IERR(1) = 8000
         RETURN
      ENDIF
!
! *** ERROR CHECK ***
! NX,NY,NZ
      IF( ( NX .LT. 2 ) .OR. ( NY .LT. 2 ) .OR. ( NZ .LT. 2 ) ) THEN
         IERR(1) = 3000
         RETURN
      ENDIF
      IF( LY .LT. LY_D*ISIZE ) THEN
         IERR(1) = 3010
         RETURN
      ENDIF
      IF( LY_D .LT. NY_D ) THEN
         IERR(1) = 3020
         RETURN
      ENDIF
      IF( LZ .LT. NZ_D*ISIZE ) THEN
         IERR(1) = 3030
         RETURN
      ENDIF
      IF( LZ_D .LT. NZ_D ) THEN
         IERR(1) = 3040
         RETURN
      ENDIF
! ISW
      IF( ( ISW .NE. -1 ) .AND. ( ISW .NE. 1 ) .AND. (ISW .NE. 0 ) )THEN
         IERR(1) = 3050
      ENDIF
! NT
      IF( NT .LT. 1 ) THEN
         IERR(1) = 3060
         RETURN
      ENDIF
! INDEX
      IF( ( INDEX .LT. 2 ) .OR. ( INDEX .GT. 3 ) ) THEN
         IERR(1) = 3100
         RETURN
      ENDIF
!
      TIME = 0.0_PREC
!
! *** ONLY DATA TRANS ***
      IF( ISW .EQ. 0 ) THEN
         T1 = MPI_WTIME()
         if(index.eq.2) then
            CALL FFT_TRANS_YZ(A,LX,LY,LY_D,NY_D,NZ_D,B,FFT_WORLD,IERR)
         elseif(index.eq.3) then
            CALL FFT_TRANS_ZY(A,LX,LY,LY_D,NY_D,NZ_D,B,FFT_WORLD,IERR)
         endif
         T2 = MPI_WTIME()
         TIME(4) = T2 - T1
         TIME(5) = T2 - T1
         RETURN
      ENDIF
!
! WORKING AREA
      MAXSIZE = MAX(NX,NY,NZ)
      ALLOCATE(TRIGS(MAXSIZE),STAT=IERR(2))
      IF( IERR(2) .NE. 0 ) THEN
         IERR(1) = 7000
         RETURN
      ENDIF

      T0 = MPI_WTIME()
      if(index.eq.2) then !!! *** Y==>Z ***
!
! X FFT
      T1 = MPI_WTIME()
#ifdef ASLFFT
      if (mod(lx,4).eq.0) then
         ll=nz
         if(mod(nz,2).ne.0 .and. (lz.gt.nz)) then
            ll=ll+1
         endif
         call ztransb(a,lx,b,ly_d*ll,nx,ly_d*nz,nt)
         CALL HFCMFB(NX,LY_D*NZ,B,LY_D*LL,1,ISW,IFAX,TRIGS,A,NT,KERR)
         call ztransb(b,ly_d*ll,a,lx,ly_d*nz,nx,nt)
      else
         CALL HFCMFB(NX,LY_D*NZ,A,1,LX,ISW,IFAX,TRIGS,B,NT,KERR)
      endif
#elif SRFFT
      n=NX
      ip=0
      IFCNT=0
      DO I=1,MAXFAX
        IFACT=NFAX(I)
        DO WHILE (MOD(N,IFACT).EQ.0)
          ip(i)=ip(i)+1
          N=N/IFACT
        ENDDO
        IF(N.EQ.1) EXIT
      ENDDO
      if(isw.eq. 1) iopt(1)=1
      if(isw.eq.-1) iopt(1)=2
      iopt(2)=1
      iopt(3)=1
      nwkw=1
      call hzft5s(A,LX,LY_D*NZ,NX,ip,LY_D*NZ,iopt,trigs,B,nwkw,kerr)
#elif FFTW3
      NFFTW3(1) = NX
      NEMBED(1) = LX
      call dfftw_plan_many_dft
     &     (PLAN,FFTW_RANK,NFFTW3,LY_D*NZ,
     &      A,NEMBED,1,LX,
     &      A,NEMBED,1,LX,
     &      ISWFFTW,FFTW_ESTIMATE)
      call dfftw_execute(PLAN)
      call dfftw_destroy_plan(PLAN)
#elif MKLFFT
      Status = DftiCreateDescriptor(Desc_Handle,DFTI_DOUBLE,
     &         DFTI_COMPLEX, 1, NX)
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_NUMBER_OF_TRANSFORMS, LY_D*NZ)
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_INPUT_DISTANCE, LX)
     !!$ Status = DftiSetValue(Desc_Handle,
     !!$&         DFTI_OUTPUT_DISTANCE, LX)
      Status = DftiCommitDescriptor(Desc_Handle)
      IF(ISW == 1) THEN
         Status = DftiComputeForward(Desc_Handle,A)
      ELSE
         Status = DftiComputeBackward(Desc_Handle,A)
      ENDIF
      Status = DftiFreeDescriptor(Desc_Handle)
#endif
      T2 = MPI_WTIME()
      TIME(1) = T2 - T1
      IF( KERR .GE. 4000 ) THEN
         IERR(1) = 4000
         RETURN
      ENDIF
!
! Z FFT
      T1 = MPI_WTIME()
#ifdef ASLFFT
      CALL HFCMFB(NZ,LX*NY_D,A,LX*LY_D,1,ISW,IFAX,TRIGS,B,NT,KERR)
#elif SRFFT
      n=NZ
      ip=0
      IFCNT=0
      DO I=1,MAXFAX
        IFACT=NFAX(I)
        DO WHILE (MOD(N,IFACT).EQ.0)
          ip(i)=ip(i)+1
          N=N/IFACT
        ENDDO
        IF(N.EQ.1) EXIT
      ENDDO
      if(isw.eq. 1) iopt(1)=1
      if(isw.eq.-1) iopt(1)=2
      iopt(2)=1
      iopt(3)=2
      nwkw=1
      call hzft5s(A,LX*LY_D,NZ,NZ,ip,LX*NY_D,iopt,trigs,B,nwkw,kerr)
#elif FFTW3
      NFFTW3(1) = NZ
      NEMBED(1) = LZ
      call dfftw_plan_many_dft
     &     (PLAN,FFTW_RANK,NFFTW3,LX*NY_D,
     &      A,NEMBED,LX*LY_D,1,
     &      A,NEMBED,LX*LY_D,1,
     &      ISWFFTW,FFTW_ESTIMATE)
      call dfftw_execute(PLAN)
      call dfftw_destroy_plan(PLAN)
#elif MKLFFT
      Status = DftiCreateDescriptor(Desc_Handle,DFTI_DOUBLE,
     &         DFTI_COMPLEX, 1, NZ)
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_NUMBER_OF_TRANSFORMS, LX*NY_D)
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_INPUT_DISTANCE, 1)
     !!$ Status = DftiSetValue(Desc_Handle,
     !!$&         DFTI_OUTPUT_DISTANCE, 1)
      Stride(1)=0
      Stride(2)=LX*LY_D
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_INPUT_STRIDES, Stride)
     !!$ Status = DftiSetValue(Desc_Handle,
     !!$&         DFTI_OUTPUT_STRIDES, Stride)
      Status = DftiCommitDescriptor(Desc_Handle)
      IF(ISW == 1) THEN
         Status = DftiComputeForward(Desc_Handle,A)
      ELSE
         Status = DftiComputeBackward(Desc_Handle,A)
      ENDIF
      Status = DftiFreeDescriptor(Desc_Handle)
#endif
      T2 = MPI_WTIME()
      TIME(3) = T2 - T1
      IF( KERR .GE. 4000 ) THEN
         IERR(1) = 4000
         RETURN
      ENDIF
!
! *** TRANSPORT Y==>Z ***
      T1 = MPI_WTIME()
      CALL FFT_TRANS_YZ(A,LX,LY,LY_D,NY_D,NZ_D,B,FFT_WORLD,IERR)
      T2 = MPI_WTIME()
      TIME(4) = T2 - T1
      IF( IERR(1) .NE. 0 ) RETURN
!
! Y FFT
      T1 = MPI_WTIME()
!CDIR PARALLEL DO PRIVATE(I,J,K)
      DO K=0,NZ_D-1
      DO J=0,NY-1
      DO I=0,NX-1
         B(I+K*NX+J*NX*NZ_D)=A(I+J*LX+K*LX*LY)
      ENDDO
      ENDDO
      ENDDO
#ifdef ASLFFT
      CALL HFCMFB(NY,NX*NZ_D,B,NX*NZ_D,1,ISW,IFAX,TRIGS,A,NT,KERR)
#elif SRFFT
      n=NY
      ip=0
      IFCNT=0
      DO I=1,MAXFAX
        IFACT=NFAX(I)
        DO WHILE (MOD(N,IFACT).EQ.0)
          ip(i)=ip(i)+1
          N=N/IFACT
        ENDDO
        IF(N.EQ.1) EXIT
      ENDDO
      if(isw.eq. 1) iopt(1)=1
      if(isw.eq.-1) iopt(1)=2
      iopt(2)=1
      iopt(3)=2
      nwkw=1
      call hzft5s(B,NX*NZ_D,NY,NY,ip,NX*NZ_D,iopt,trigs,A,nwkw,kerr)
#elif FFTW3
      NFFTW3(1) = NY
      NEMBED(1) = NY
      call dfftw_plan_many_dft
     &     (PLAN,FFTW_RANK,NFFTW3,NX*NZ_D,
     &      B,NEMBED,NX*NZ_D,1,
     &      B,NEMBED,NX*NZ_D,1,
     &      ISWFFTW,FFTW_ESTIMATE)
      call dfftw_execute(PLAN)
      call dfftw_destroy_plan(PLAN)
#elif MKLFFT
      Status = DftiCreateDescriptor(Desc_Handle,DFTI_DOUBLE,
     &         DFTI_COMPLEX, 1, NY)
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_NUMBER_OF_TRANSFORMS, NX*NZ_D)
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_INPUT_DISTANCE, 1)
     !!$ Status = DftiSetValue(Desc_Handle,
     !!$&         DFTI_OUTPUT_DISTANCE, 1)
      Stride(1)=0
      Stride(2)=NX*NZ_D
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_INPUT_STRIDES, Stride)
     !!$ Status = DftiSetValue(Desc_Handle,
     !!$&         DFTI_OUTPUT_STRIDES, Stride)
      Status = DftiCommitDescriptor(Desc_Handle)
      IF(ISW == 1) THEN
         Status = DftiComputeForward(Desc_Handle,B)
      ELSE
         Status = DftiComputeBackward(Desc_Handle,B)
      ENDIF
      Status = DftiFreeDescriptor(Desc_Handle)
#endif
!CDIR PARALLEL DO PRIVATE(I,J,K)
      DO K=0,NZ_D-1
      DO J=0,NY-1
      DO I=0,NX-1
         A(I+J*LX+K*LX*LY)=B(I+K*NX+J*NX*NZ_D)
      ENDDO
      ENDDO
      ENDDO
      T2 = MPI_WTIME()
      TIME(2) = T2 - T1
      IF( KERR .GE. 4000 ) THEN
         IERR(1) = 4000
         RETURN
      ENDIF

      elseif(index.eq.3) then !!! *** Z==>Y ***
!
! X FFT
      T1 = MPI_WTIME()
#ifdef ASLFFT
      if (mod(lx,4).eq.0) then
         ll_d=nz_d
         if(mod(nz_d,2).ne.0 .and. (lz_d.gt.nz_d)) then
            ll_d=ll_d+1
         endif
         call ztransb(a,lx,b,ly*ll_d,nx,ly*nz_d,nt)
         CALL HFCMFB(NX,LY*NZ_D,B,LY*LL_D,1,ISW,IFAX,TRIGS,A,NT,KERR)
         call ztransb(b,ly*ll_d,a,lx,ly*nz_d,nx,nt)
      else
         CALL HFCMFB(NX,LY*NZ_D,A,1,LX,ISW,IFAX,TRIGS,B,NT,KERR)
      endif
#elif SRFFT
      n=NX
      ip=0
      IFCNT=0
      DO I=1,MAXFAX
        IFACT=NFAX(I)
        DO WHILE (MOD(N,IFACT).EQ.0)
          ip(i)=ip(i)+1
          N=N/IFACT
        ENDDO
        IF(N.EQ.1) EXIT
      ENDDO
      if(isw.eq. 1) iopt(1)=1
      if(isw.eq.-1) iopt(1)=2
      iopt(2)=1
      iopt(3)=1
      nwkw=1
      call hzft5s(A,LX,LY*NZ_D,NX,ip,LY*NZ_D,iopt,trigs,B,nwkw,kerr)
#elif FFTW3
      NFFTW3(1) = NX
      NEMBED(1) = LX
      call dfftw_plan_many_dft
     &     (PLAN,FFTW_RANK,NFFTW3,LY*NZ_D,
     &      A,NEMBED,1,LX,
     &      A,NEMBED,1,LX,
     &      ISWFFTW,FFTW_ESTIMATE)
      call dfftw_execute(PLAN)
      call dfftw_destroy_plan(PLAN)
#elif MKLFFT
      Status = DftiCreateDescriptor(Desc_Handle,DFTI_DOUBLE,
     &         DFTI_COMPLEX, 1, NX)
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_NUMBER_OF_TRANSFORMS, LY*NZ_D)
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_INPUT_DISTANCE, LX)
     !!$ Status = DftiSetValue(Desc_Handle,
     !!$&         DFTI_OUTPUT_DISTANCE, LX)
      Status = DftiCommitDescriptor(Desc_Handle)
      IF(ISW == 1) THEN
         Status = DftiComputeForward(Desc_Handle,A)
      ELSE
         Status = DftiComputeBackward(Desc_Handle,A)
      ENDIF
      Status = DftiFreeDescriptor(Desc_Handle)
#endif
      T2 = MPI_WTIME()
      TIME(1) = T2 - T1
      IF( KERR .GE. 4000 ) THEN
         IERR(1) = 4000
         RETURN
      ENDIF
!
! Y FFT
      T1 = MPI_WTIME()
!CDIR PARALLEL DO PRIVATE(I,J,K)
      DO K=0,NZ_D-1
      DO J=0,NY-1
      DO I=0,NX-1
         B(I+K*NX+J*NX*NZ_D)=A(I+J*LX+K*LX*LY)
      ENDDO
      ENDDO
      ENDDO
#ifdef ASLFFT
      CALL HFCMFB(NY,NX*NZ_D,B,NX*NZ_D,1,ISW,IFAX,TRIGS,A,NT,KERR)
#elif SRFFT
      n=NY
      ip=0
      IFCNT=0
      DO I=1,MAXFAX
        IFACT=NFAX(I)
        DO WHILE (MOD(N,IFACT).EQ.0)
          ip(i)=ip(i)+1
          N=N/IFACT
        ENDDO
        IF(N.EQ.1) EXIT
      ENDDO
      if(isw.eq. 1) iopt(1)=1
      if(isw.eq.-1) iopt(1)=2
      iopt(2)=1
      iopt(3)=2
      nwkw=1
      call hzft5s(B,NX*NZ_D,NY,NY,ip,NX*NZ_D,iopt,trigs,A,nwkw,kerr)
#elif FFTW3
      NFFTW3(1) = NY
      NEMBED(1) = NY
      call dfftw_plan_many_dft
     &     (PLAN,FFTW_RANK,NFFTW3,NX*NZ_D,
     &      B,NEMBED,NX*NZ_D,1,
     &      B,NEMBED,NX*NZ_D,1,
     &      ISWFFTW,FFTW_ESTIMATE)
      call dfftw_execute(PLAN)
      call dfftw_destroy_plan(PLAN)
#elif MKLFFT
      Status = DftiCreateDescriptor(Desc_Handle,DFTI_DOUBLE,
     &         DFTI_COMPLEX, 1, NY)
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_NUMBER_OF_TRANSFORMS, NX*NZ_D)
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_INPUT_DISTANCE, 1)
     !!$ Status = DftiSetValue(Desc_Handle,
     !!$&         DFTI_OUTPUT_DISTANCE, 1)
      Stride(1)=0
      Stride(2)=NX*NZ_D
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_INPUT_STRIDES, Stride)
     !!$ Status = DftiSetValue(Desc_Handle,
     !!$&         DFTI_OUTPUT_STRIDES, Stride)
      Status = DftiCommitDescriptor(Desc_Handle)
      IF(ISW == 1) THEN
         Status = DftiComputeForward(Desc_Handle,B)
      ELSE
         Status = DftiComputeBackward(Desc_Handle,B)
      ENDIF
      Status = DftiFreeDescriptor(Desc_Handle)
#endif
!CDIR PARALLEL DO PRIVATE(I,J,K)
      DO K=0,NZ_D-1
      DO J=0,NY-1
      DO I=0,NX-1
         A(I+J*LX+K*LX*LY)=B(I+K*NX+J*NX*NZ_D)
      ENDDO
      ENDDO
      ENDDO
      T2 = MPI_WTIME()
      TIME(2) = T2 - T1
      IF( KERR .GE. 4000 ) THEN
         IERR(1) = 4000
         RETURN
      ENDIF
!
! *** TRANSPORT Z==>Y ***
      T1 = MPI_WTIME()
      CALL FFT_TRANS_ZY(A,LX,LY,LY_D,NY_D,NZ_D,B,FFT_WORLD,IERR)
      T2 = MPI_WTIME()
      TIME(4) = T2 - T1
      IF( IERR(1) .NE. 0 ) RETURN
!
! Z FFT
      T1 = MPI_WTIME()
#ifdef ASLFFT
      CALL HFCMFB(NZ,LX*NY_D,A,LX*LY_D,1,ISW,IFAX,TRIGS,B,NT,KERR)
#elif SRFFT
      n=NZ
      ip=0
      IFCNT=0
      DO I=1,MAXFAX
        IFACT=NFAX(I)
        DO WHILE (MOD(N,IFACT).EQ.0)
          ip(i)=ip(i)+1
          N=N/IFACT
        ENDDO
        IF(N.EQ.1) EXIT
      ENDDO
      if(isw.eq. 1) iopt(1)=1
      if(isw.eq.-1) iopt(1)=2
      iopt(2)=1
      iopt(3)=2
      nwkw=1
      call hzft5s(A,LX*LY_D,NZ,NZ,ip,LX*NY_D,iopt,trigs,B,nwkw,kerr)
#elif FFTW3
      NFFTW3(1) = NZ
      NEMBED(1) = LZ
      call dfftw_plan_many_dft
     &     (PLAN,FFTW_RANK,NFFTW3,LX*NY_D,
     &      A,NEMBED,LX*LY_D,1,
     &      A,NEMBED,LX*LY_D,1,
     &      ISWFFTW,FFTW_ESTIMATE)
      call dfftw_execute(PLAN)
      call dfftw_destroy_plan(PLAN)
#elif MKLFFT
      Status = DftiCreateDescriptor(Desc_Handle,DFTI_DOUBLE,
     &         DFTI_COMPLEX, 1, NZ)
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_NUMBER_OF_TRANSFORMS, LX*NY_D)
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_INPUT_DISTANCE, 1)
     !!$ Status = DftiSetValue(Desc_Handle,
     !!$&         DFTI_OUTPUT_DISTANCE, 1)
      Stride(1)=0
      Stride(2)=LX*LY_D
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_INPUT_STRIDES, Stride)
     !!$ Status = DftiSetValue(Desc_Handle,
     !!$&         DFTI_OUTPUT_STRIDES, Stride)
      Status = DftiCommitDescriptor(Desc_Handle)
      IF(ISW == 1) THEN
         Status = DftiComputeForward(Desc_Handle,A)
      ELSE
         Status = DftiComputeBackward(Desc_Handle,A)
      ENDIF
      Status = DftiFreeDescriptor(Desc_Handle)
#endif
      T2 = MPI_WTIME()
      TIME(3) = T2 - T1
      IF( KERR .GE. 4000 ) THEN
         IERR(1) = 4000
         RETURN
      ENDIF

      endif
!
! *** FINALIZE ***
      T1 = MPI_WTIME()
      TIME(5) = T1 - T0
      DEALLOCATE(TRIGS,STAT=IERR(2))
      IF( IERR(2) .NE. 0 ) THEN
         IERR(1) = 7100
         RETURN
      ENDIF

      RETURN

      END

!     ****************************************************************
!     AUTHORS: T.KOKUBO(NEC), K.TAKADA(NIS), T.FURUTA(NIS), R.SAKA(NIS)
!
!     REGD. ON 03/15/02 BY R.S.
!     UPDD. ON 09/08/06 BY T.K.
!
      SUBROUTINE FFT_R3D(NX,NY,NZ,NY_D,NZ_D,A,LX,LY,LZ,LY_D,LZ_D,
     &                   index,isw, TIME, B,NT,FFT_WORLD, IERR )

#ifdef MKLFFT
      USE MKL_DFTI
#endif

      INTEGER,PARAMETER :: NPREC=8
!
! ::: ARGUMENTS :::
!
      INTEGER,INTENT(IN) :: NX,NY,NZ,NY_D,NZ_D, LX,LY,LZ,LY_D,LZ_D,
     &                      index,NT, FFT_WORLD,isw
      INTEGER,INTENT(INOUT),DIMENSION(:) :: IERR(2)
      REAL(KIND=NPREC),DIMENSION(5),INTENT(OUT) :: TIME
      REAL(KIND=NPREC),DIMENSION(:),INTENT(INOUT) :: A(0:*),B(0:*)
!
#ifdef ASLFFT
      INTEGER,DIMENSION(:) :: IFAX(20)
      integer :: kerr
#elif SRFFT
      INTEGER :: ip(4),iopt(3),n,nwkw
      INTEGER,PARAMETER :: MAXFAX = 4
      INTEGER,DIMENSION(MAXFAX) :: NFAX = (/2,3,5,7/)
#elif FFTW3
      INTEGER(KIND=NPREC) :: PLAN
      INTEGER,DIMENSION(1) :: NFFTW3
      INTEGER,DIMENSION(1) :: NEMBED
      INTEGER,DIMENSION(1) :: NEREAL
      INTEGER,PARAMETER :: FFTW_RANK = 1
      INTEGER,PARAMETER :: FFTW_ESTIMATE = 64
      INTEGER :: ISWFFTW
#elif MKLFFT
      INTEGER :: STRIDE(2), Status, II
      TYPE(DFTI_DESCRIPTOR), POINTER :: Desc_Handle
      COMPLEX(KIND=NPREC), allocatable :: C(:)
#endif
      INTEGER :: I,J,K,NX2,LX2,ISIZE,ll,ll_d
      COMPLEX(KIND=NPREC),DIMENSION(:),ALLOCATABLE :: TRIGS
      REAL(KIND=NPREC) :: T0,T1
!
      INCLUDE "mpif.h"
!
#ifdef FFTW3
      ISWFFTW=-ISW
#ifdef _OPENMP
      if(NT>1) then
         call dfftw_init_threads
         call dfftw_plan_with_nthreads(NT)
      end if
!
#endif
#endif
      CALL MPI_COMM_SIZE(FFT_WORLD,ISIZE,IERR(2))
      IF( IERR(2) .NE. 0 ) THEN
         IERR(1) = 8000
         RETURN
      ENDIF
!
      DO I=1,5
         TIME(I) = 0.D0
      ENDDO
      T0      = MPI_WTIME()
      IERR(1) = 0
      IERR(2) = 0

      NX2=(NX+2)/2
      LX2=LX/2

! *** ERROR CHECK
! ++ (1) NX,NY,NZ
      IF(  NX .LT. 2 .OR.
     &     NY .LT. 2 .OR. 
     &     NZ .LT. 2     )THEN
         IERR(1)=3000
         GOTO 999
      ENDIF

! ++ (2) LX
      IF( MOD(NX,2) .EQ. 1 )THEN
         IF( LX .LT. NX+1 )THEN
            IERR(1)=3010
            GOTO 999
         ENDIF
      ELSE
         IF( LX .LT. NX+2 )THEN
            IERR(1)=3010
            GOTO 999
         ENDIF
      ENDIF

! ++ (3) LY,LZ
      IF(  LY .LT. LY_D*ISIZE .OR.
     &     LZ .LT. LZ_D*ISIZE )THEN
         IERR(1)=3020
         GOTO 999
      ENDIF

! ++ (4) LY_d,LZ_d
      IF(  LY_d .LT. NY_d .OR.
     &     LZ_d .LT. NZ_d      )THEN
         IERR(1)=3030
         GOTO 999
      ENDIF

! ++ (5) MOD(LX,2) 
      IF( MOD(LX,2) .EQ. 1 )THEN
         IERR(1)=3040
         GOTO 999
      ENDIF

! ++ (6) NT
      IF( NT .LT. 1 )THEN
         IERR(1)=3050
         GOTO 999
      ENDIF

! ++ (7) index
      IF(  index .LT. 2 .OR.
     &     index .GT. 3      )THEN
         IERR(1)=3100
         GOTO 999
      ENDIF

! *** ONLY DATA TRANS ***
      IF( ISW .EQ. 0 ) THEN
         T1 = MPI_WTIME()
         if(index.eq.2) then
            CALL FFT_TRANS_YZ(A,LX2,LY,LY_D,NY_D,NZ_D,B,FFT_WORLD,IERR)
         elseif(index.eq.3) then
            CALL FFT_TRANS_ZY(A,LX2,LY,LY_D,NY_D,NZ_D,B,FFT_WORLD,IERR)
         endif
         T2 = MPI_WTIME()
         TIME(4) = T2 - T1
         TIME(5) = T2 - T1
         RETURN
      ENDIF

! *** MEMORY ALLOCATE ***
      NNMAX=MAX(NX,NY,NZ)
      ALLOCATE( TRIGS(NNMAX),STAT=IERR(2) )
      IF( IERR(2) .NE. 0 )THEN
         IERR(1)=7000
         GOTO 999
      ENDIF
#if MKLFFT
      NNMAX=MAX(LX*LY*LZ_D,LX*LY_D*LZ)
      ALLOCATE( C(0:NNMAX-1),STAT=IERR(2) )
      IF( IERR(2) .NE. 0 )THEN
         IERR(1)=7001
         GOTO 999
      ENDIF
#endif

      if(isw.eq.1) then
!!!!! forward(real to half complex FFT)

      if(index.eq.2) then !!! *** Y==>Z ***

! *** X - FFT ( REAL ) *** A( NX, NY_D, NZ )
      T1=MPI_WTIME()
#ifdef ASLFFT
      if (mod(lx,4).eq.0) then
        ll=nz
        if(mod(nz,2).ne.0 .and. (lz.gt.nz)) then
          ll=ll+1
        endif
        call transb(a,lx,b,ly_d*ll,nx+2,ly_d*nz,nt)
        CALL QFRMFB(
     &       NX, LY_D*NZ, B, ly_d*ll, 1,
     &       ISW,  IFAX, TRIGS,
     &       A, NT, IERR(1) )
        call transb(b,ly_d*ll,a,lx,ly_d*nz,nx+2,nt)
      else
        CALL QFRMFB(
     &       NX, LY_D*NZ, A, 1, LX,
     &       ISW,  IFAX, TRIGS,
     &       B, NT, IERR(1) )
      endif
#elif SRFFT
      n=NX
      ip=0
      IFCNT=0
      DO I=1,MAXFAX
        IFACT=NFAX(I)
        DO WHILE (MOD(N,IFACT).EQ.0)
          ip(i)=ip(i)+1
          N=N/IFACT
        ENDDO
        IF(N.EQ.1) EXIT
      ENDDO
      iopt(1)=1
      iopt(2)=1
      iopt(3)=1
      nwkw=1
      call hdfz5m(A,LX,LY_D*NZ,NX,ip,LY_D*NZ,iopt,trigs,B,nwkw,ierr(1))
#elif FFTW3
      NFFTW3(1) = NX
      NEREAL(1) = LX
      NEMBED(1) = LX2
      call dfftw_plan_many_dft_r2c
     &     (PLAN,FFTW_RANK,NFFTW3,LY_D*NZ,
     &      A,NEREAL,1,LX,
     &      A,NEMBED,1,LX2,
     &      FFTW_ESTIMATE)
      call dfftw_execute(PLAN)
      call dfftw_destroy_plan(PLAN)
#elif MKLFFT
      Status = DftiCreateDescriptor(Desc_Handle,DFTI_DOUBLE,
     &         DFTI_REAL, 1, NX)
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_NUMBER_OF_TRANSFORMS, LY_D*NZ)
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_INPUT_DISTANCE, LX)
      Status = DftiCommitDescriptor(Desc_Handle)
      Status = DftiComputeForward(Desc_Handle,A)
      Status = DftiFreeDescriptor(Desc_Handle)
#endif
      TIME(1)=MPI_WTIME() - T1
!     write(*,*)'xfft ierr=',ierr(1)
      IF( IERR(1) .GE. 4000 )THEN
         IERR(1)=4000
         GOTO 999
      ENDIF

! *** Z - FFT ( COMPLEX ) *** A( NX, NY_D, NZ )
      T1=MPI_WTIME()
#ifdef ASLFFT
      CALL HFCMFB(
     &     NZ, LX2*NY_D, A, LX2*LY_D, 1,
     &     ISW, IFAX, TRIGS,
     &     B, NT, IERR(1) )
#elif SRFFT
      n=NZ
      ip=0
      IFCNT=0
      DO I=1,MAXFAX
        IFACT=NFAX(I)
        DO WHILE (MOD(N,IFACT).EQ.0)
          ip(i)=ip(i)+1
          N=N/IFACT
        ENDDO
        IF(N.EQ.1) EXIT
      ENDDO
      iopt(1)=1
      iopt(2)=1
      iopt(3)=2
      nwkw=1
      call hzft5s(A,LX2*LY_D,NZ,NZ,ip,LX2*NY_D,
     &            iopt,trigs,B,nwkw,IERR(1))
#elif FFTW3
      NFFTW3(1) = NZ
      NEMBED(1) = LZ
      call dfftw_plan_many_dft
     &     (PLAN,FFTW_RANK,NFFTW3,LX2*NY_D,
     &      A,NEMBED,LX2*LY_D,1,
     &      A,NEMBED,LX2*LY_D,1,
     &      ISWFFTW,FFTW_ESTIMATE)
      call dfftw_execute(PLAN)
      call dfftw_destroy_plan(PLAN)
#elif MKLFFT
      Status = DftiCreateDescriptor(Desc_Handle,DFTI_DOUBLE,
     &         DFTI_COMPLEX, 1, NZ)
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_NUMBER_OF_TRANSFORMS, LX2*NY_D)
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_INPUT_DISTANCE, 1)
      Stride(1)=0
      Stride(2)=LX2*LY_D
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_INPUT_STRIDES, Stride)
      Status = DftiCommitDescriptor(Desc_Handle)
!!!!  next line may not be working, because A and B are real array....
      call DCOPY(LZ*LX2*LY_D*2,A,1,C,1)
      Status = DftiComputeForward(Desc_Handle,C)
      call DCOPY(LZ*LX2*LY_D*2,C,1,A,1)
      Status = DftiFreeDescriptor(Desc_Handle)
#endif
      TIME(3)=MPI_WTIME() - T1
!     write(*,*)'zfft ierr=',ierr(1)
      IF( IERR(1) .GE. 4000 )THEN
         IERR(1)=4000
         GOTO 999
      ENDIF

!     write(*,*) 'fft_rc3d'

! *** DATA TRANSPOSE *** Y==>Z
      T1=MPI_WTIME()
      CALL FFT_TRANS_YZ(A,LX2,LY,LY_D,NY_D,NZ_D,B,FFT_WORLD,IERR)
      TIME(4)=MPI_WTIME() - T1
      IF( IERR(1) .GE. 4000 )THEN
         IERR(1)=4000
         GOTO 999
      ENDIF

! *** Y FFT  ( COMPLEX ) *** A( NX2, NY, NZ_D )
      T1=MPI_WTIME()
!CDIR PARALLEL DO PRIVATE(I,J,K)
      DO K=0,NZ_D-1
      DO J=0,NY-1
      DO I=0,NX2*2-1
         B(I+K*NX2*2+J*NX2*2*NZ_D)=A(I+J*LX+K*LX*LY)
      ENDDO
      ENDDO
      ENDDO
#ifdef ASLFFT
      CALL HFCMFB(NY,NX2*NZ_D,B,NX2*NZ_D,1,ISW,IFAX,TRIGS,A,NT,KERR)
#elif SRFFT
      n=NY
      ip=0
      IFCNT=0
      DO I=1,MAXFAX
        IFACT=NFAX(I)
        DO WHILE (MOD(N,IFACT).EQ.0)
          ip(i)=ip(i)+1
          N=N/IFACT
        ENDDO
        IF(N.EQ.1) EXIT
      ENDDO
      iopt(1)=1
      iopt(2)=1
      iopt(3)=2
      nwkw=1
      call hzft5s(B,NX2*NZ_D,NY,NY,ip,NX2*NZ_D,
     &            iopt,trigs,A,nwkw,IERR(1))
#elif FFTW3
      NFFTW3(1) = NY
      NEMBED(1) = NY
      call dfftw_plan_many_dft
     &     (PLAN,FFTW_RANK,NFFTW3,NX2*NZ_D,
     &      B,NEMBED,NX2*NZ_D,1,
     &      B,NEMBED,NX2*NZ_D,1,
     &      ISWFFTW,FFTW_ESTIMATE)
      call dfftw_execute(PLAN)
      call dfftw_destroy_plan(PLAN)
#elif MKLFFT
      Status = DftiCreateDescriptor(Desc_Handle,DFTI_DOUBLE,
     &         DFTI_COMPLEX, 1, NY)
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_NUMBER_OF_TRANSFORMS, NX2*NZ_D)
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_INPUT_DISTANCE, 1)
      Stride(1)=0
      Stride(2)=NX2*NZ_D
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_INPUT_STRIDES, Stride)
      Status = DftiCommitDescriptor(Desc_Handle)
!!!!  next line may not be working, because A and B are real array....
      call DCOPY(NY*NX2*NZ_D*2,B,1,C,1)
      Status = DftiComputeForward(Desc_Handle,C)
      call DCOPY(NY*NX2*NZ_D*2,C,1,B,1)
      Status = DftiFreeDescriptor(Desc_Handle)
#endif
!CDIR PARALLEL DO PRIVATE(I,J,K)
      DO K=0,NZ_D-1
      DO J=0,NY-1
      DO I=0,NX2*2-1
         A(I+J*LX+K*LX*LY)=B(I+K*NX2*2+J*NX2*2*NZ_D)
      ENDDO
      ENDDO
      ENDDO
      IF( IERR(1) .GE. 4000 )THEN
         IERR(1)=4000
         TIME(2)=MPI_WTIME() - T1
         GOTO 999
      ENDIF
      TIME(2)=MPI_WTIME() - T1
!     write(*,*)'yfft ierr=',ierr(1)

      elseif(index.eq.3) then !!! *** Z==>Y ***

! *** X - FFT ( REAL ) *** A( NX, NY, NZ_D )
      T1=MPI_WTIME()
#ifdef ASLFFT
      if (mod(lx,4).eq.0) then
        ll_d=nz_d
        if(mod(nz_d,2).ne.0 .and. (lz_d.gt.nz_d)) then
          ll_d=ll_d+1
        endif
        call transb(a,lx,b,ly*ll_d,nx+2,ly*nz_d,nt)
        CALL QFRMFB(
     &       NX, LY*NZ_D, B, ly*ll_d, 1,
     &       ISW,  IFAX, TRIGS,
     &       A, NT, IERR(1) )
        call transb(b,ly*ll_d,a,lx,ly*nz_d,nx+2,nt)
      else
        CALL QFRMFB(
     &       NX, LY*NZ_D, A, 1, LX,
     &       ISW,  IFAX, TRIGS,
     &       B, NT, IERR(1) )
      endif
#elif SRFFT
      n=NX
      ip=0
      IFCNT=0
      DO I=1,MAXFAX
        IFACT=NFAX(I)
        DO WHILE (MOD(N,IFACT).EQ.0)
          ip(i)=ip(i)+1
          N=N/IFACT
        ENDDO
        IF(N.EQ.1) EXIT
      ENDDO
      iopt(1)=1
      iopt(2)=1
      iopt(3)=1
      nwkw=1
      call hdfz5m(A,LX,LY*NZ_D,NX,ip,LY*NZ_D,iopt,trigs,B,nwkw,ierr(1))
#elif FFTW3
      NFFTW3(1) = NX
      NEREAL(1) = LX
      NEMBED(1) = LX2
      call dfftw_plan_many_dft_r2c
     &     (PLAN,FFTW_RANK,NFFTW3,LY*NZ_D,
     &      A,NEREAL,1,LX,
     &      A,NEMBED,1,LX2,
     &      FFTW_ESTIMATE)
      call dfftw_execute(PLAN)
      call dfftw_destroy_plan(PLAN)
#elif MKLFFT
      Status = DftiCreateDescriptor(Desc_Handle,DFTI_DOUBLE,
     &         DFTI_REAL, 1, NX)
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_NUMBER_OF_TRANSFORMS, LY*NZ_D)
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_INPUT_DISTANCE, LX)
      Status = DftiCommitDescriptor(Desc_Handle)
      Status = DftiComputeForward(Desc_Handle,A)
      Status = DftiFreeDescriptor(Desc_Handle)
#endif
      TIME(1)=MPI_WTIME() - T1
      IF( IERR(1) .GE. 4000 )THEN
         IERR(1)=4000
         GOTO 999
      ENDIF

! *** Y - FFT ( COMPLEX ) *** A( NX2, NY, NZ_D )
      T1=MPI_WTIME()
!CDIR PARALLEL DO PRIVATE(I,J,K)
      DO K=0,NZ_D-1
      DO J=0,NY-1
      DO I=0,NX2*2-1
         B(I+K*NX2*2+J*NX2*2*NZ_D)=A(I+J*LX+K*LX*LY)
      ENDDO
      ENDDO
      ENDDO
#ifdef ASLFFT
      CALL HFCMFB(NY,NX2*NZ_D,B,NX2*NZ_D,1,ISW,IFAX,TRIGS,A,NT,KERR)
#elif SRFFT
      n=NY
      ip=0
      IFCNT=0
      DO I=1,MAXFAX
        IFACT=NFAX(I)
        DO WHILE (MOD(N,IFACT).EQ.0)
          ip(i)=ip(i)+1
          N=N/IFACT
        ENDDO
        IF(N.EQ.1) EXIT
      ENDDO
      iopt(1)=1
      iopt(2)=1
      iopt(3)=2
      nwkw=1
      call hzft5s(B,NX2*NZ_D,NY,NY,ip,NX2*NZ_D,
     &            iopt,trigs,A,nwkw,IERR(1))
#elif FFTW3
      NFFTW3(1) = NY
      NEMBED(1) = NY
      call dfftw_plan_many_dft
     &     (PLAN,FFTW_RANK,NFFTW3,NX2*NZ_D,
     &      B,NEMBED,NX2*NZ_D,1,
     &      B,NEMBED,NX2*NZ_D,1,
     &      ISWFFTW,FFTW_ESTIMATE)
      call dfftw_execute(PLAN)
      call dfftw_destroy_plan(PLAN)
#elif MKLFFT
      Status = DftiCreateDescriptor(Desc_Handle,DFTI_DOUBLE,
     &         DFTI_COMPLEX, 1, NY)
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_NUMBER_OF_TRANSFORMS, NX2*NZ_D)
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_INPUT_DISTANCE, 1)
      Stride(1)=0
      Stride(2)=NX2*NZ_D
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_INPUT_STRIDES, Stride)
      Status = DftiCommitDescriptor(Desc_Handle)
!!!!  next line may not be working, because A and B are real array....
      call DCOPY(NY*NX2*NZ_D*2,B,1,C,1)
      Status = DftiComputeForward(Desc_Handle,C)
      call DCOPY(NY*NX2*NZ_D*2,C,1,B,1)
      Status = DftiFreeDescriptor(Desc_Handle)
#endif
!CDIR PARALLEL DO PRIVATE(I,J,K)
      DO K=0,NZ_D-1
      DO J=0,NY-1
      DO I=0,NX2*2-1
         A(I+J*LX+K*LX*LY)=B(I+K*NX2*2+J*NX2*2*NZ_D)
      ENDDO
      ENDDO
      ENDDO
      IF( IERR(1) .GE. 4000 )THEN
         IERR(1)=4000
         TIME(2)=MPI_WTIME() - T1
         GOTO 999
      ENDIF
      TIME(2)=MPI_WTIME() - T1

!     write(*,*) 'fft_rc3d'

! *** MPI TRANSPOSE *** Z==Y
      T1=MPI_WTIME()
      CALL FFT_TRANS_ZY(A,LX2,LY,LY_D,NY_D,NZ_D,B,FFT_WORLD,IERR)
      TIME(4)=MPI_WTIME() - T1
      IF( IERR(1) .GE. 4000 )THEN
         IERR(1)=4000
         GOTO 999
      ENDIF

! *** Z FFT ( COMPLEX ) *** A( NX, NY_D, NZ )
      T1=MPI_WTIME()
#ifdef ASLFFT
      CALL HFCMFB(
     &     NZ, LX2*NY_D, A, LX2*LY_D, 1,
     &     ISW,  IFAX, TRIGS,
     &     B, NT, IERR(1) )
#elif SRFFT
      n=NZ
      ip=0
      IFCNT=0
      DO I=1,MAXFAX
        IFACT=NFAX(I)
        DO WHILE (MOD(N,IFACT).EQ.0)
          ip(i)=ip(i)+1
          N=N/IFACT
        ENDDO
        IF(N.EQ.1) EXIT
      ENDDO
      iopt(1)=1
      iopt(2)=1
      iopt(3)=2
      nwkw=1
      call hzft5s(A,LX2*LY_D,NZ,NZ,ip,LX2*NY_D,
     &            iopt,trigs,B,nwkw,IERR(1))
#elif FFTW3
      NFFTW3(1) = NZ
      NEMBED(1) = LZ
      call dfftw_plan_many_dft
     &     (PLAN,FFTW_RANK,NFFTW3,LX2*NY_D,
     &      A,NEMBED,LX2*LY_D,1,
     &      A,NEMBED,LX2*LY_D,1,
     &      ISWFFTW,FFTW_ESTIMATE)
      call dfftw_execute(PLAN)
      call dfftw_destroy_plan(PLAN)
#elif MKLFFT
      Status = DftiCreateDescriptor(Desc_Handle,DFTI_DOUBLE,
     &         DFTI_COMPLEX, 1, NZ)
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_NUMBER_OF_TRANSFORMS, LX2*NY_D)
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_INPUT_DISTANCE, 1)
      Stride(1)=0
      Stride(2)=LX2*LY_D
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_INPUT_STRIDES, Stride)
      Status = DftiCommitDescriptor(Desc_Handle)
!!!!  next line may not be working, because A and B are real array....
      call DCOPY(LZ*LX2*LY_D*2,A,1,C,1)
      Status = DftiComputeForward(Desc_Handle,C)
      call DCOPY(LZ*LX2*LY_D*2,C,1,A,1)
      Status = DftiFreeDescriptor(Desc_Handle)
#endif
      TIME(3)=MPI_WTIME() - T1
      IF( IERR(1) .GE. 4000 )THEN
         IERR(1)=4000
         GOTO 999
      ENDIF

      endif

      else
!!!!! backward(half complex to real FFT)

      if(index.eq.2) then !!! *** Y==>Z ***

! *** Z FFT ( COMPLEX ) *** A( NX2, NY_D, NZ )

!     write(*,*) 'cr3d zfft yz'
      T1=MPI_WTIME()
#ifdef ASLFFT
      CALL HFCMFB(
     &     NZ, LX2*NY_D, A, LX2*LY_D, 1,
     &     ISW,  IFAX, TRIGS,
     &     B, NT, IERR(1) )
#elif SRFFT
      n=NZ
      ip=0
      IFCNT=0
      DO I=1,MAXFAX
        IFACT=NFAX(I)
        DO WHILE (MOD(N,IFACT).EQ.0)
          ip(i)=ip(i)+1
          N=N/IFACT
        ENDDO
        IF(N.EQ.1) EXIT
      ENDDO
      iopt(1)=2
      iopt(2)=1
      iopt(3)=2
      nwkw=1
      call hzft5s(A,LX2*LY_D,NZ,NZ,ip,LX2*NY_D,
     &            iopt,trigs,B,nwkw,IERR(1))
#elif FFTW3
      NFFTW3(1) = NZ
      NEMBED(1) = LZ
      call dfftw_plan_many_dft
     &     (PLAN,FFTW_RANK,NFFTW3,LX2*NY_D,
     &      A,NEMBED,LX2*LY_D,1,
     &      A,NEMBED,LX2*LY_D,1,
     &      ISWFFTW,FFTW_ESTIMATE)
      call dfftw_execute(PLAN)
      call dfftw_destroy_plan(PLAN)
#elif MKLFFT
      Status = DftiCreateDescriptor(Desc_Handle,DFTI_DOUBLE,
     &         DFTI_COMPLEX, 1, NZ)
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_NUMBER_OF_TRANSFORMS, LX2*NY_D)
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_INPUT_DISTANCE, 1)
      Stride(1)=0
      Stride(2)=LX2*LY_D
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_INPUT_STRIDES, Stride)
      Status = DftiCommitDescriptor(Desc_Handle)
!!!!  next line may not be working, because A and B are real array....
      call DCOPY(LZ*LX2*LY_D*2,A,1,C,1)
      Status = DftiComputeBackward(Desc_Handle,C)
      call DCOPY(LZ*LX2*LY_D*2,C,1,A,1)
      Status = DftiFreeDescriptor(Desc_Handle)
#endif
      TIME(3)=MPI_WTIME() - T1
      IF( IERR(1) .GE. 4000 )THEN
         IERR(1)=4000
         GOTO 999
      ENDIF

!     write(*,*) 'cr3d trans yz'

! *** DATA TRANSPOSE *** Y==>Z
      T1=MPI_WTIME()
      CALL FFT_TRANS_YZ(A,LX2,LY,LY_D,NY_D,NZ_D,B,FFT_WORLD,IERR)
      TIME(4)=MPI_WTIME() - T1
      IF( IERR(1) .NE. 0 )THEN
         GOTO 999
      ENDIF

! *** Y - FFT ( COMPLEX ) *** A( NX2, NY, NZ_D )
!     write(*,*) 'cr3d yfft yz'
      T1=MPI_WTIME()
!CDIR PARALLEL DO PRIVATE(I,J,K)
      DO K=0,NZ_D-1
      DO J=0,NY-1
      DO I=0,NX2*2-1
         B(I+K*NX2*2+J*NX2*2*NZ_D)=A(I+J*LX+K*LX*LY)
      ENDDO
      ENDDO
      ENDDO
#ifdef ASLFFT
      CALL HFCMFB(NY,NX2*NZ_D,B,NX2*NZ_D,1,ISW,IFAX,TRIGS,A,NT,KERR)
#elif SRFFT
      n=NY
      ip=0
      IFCNT=0
      DO I=1,MAXFAX
        IFACT=NFAX(I)
        DO WHILE (MOD(N,IFACT).EQ.0)
          ip(i)=ip(i)+1
          N=N/IFACT
        ENDDO
        IF(N.EQ.1) EXIT
      ENDDO
      iopt(1)=2
      iopt(2)=1
      iopt(3)=2
      nwkw=1
      call hzft5s(B,NX2*NZ_D,NY,NY,ip,NX2*NZ_D,
     &            iopt,trigs,A,nwkw,IERR(1))
#elif FFTW3
      NFFTW3(1) = NY
      NEMBED(1) = NY
      call dfftw_plan_many_dft
     &     (PLAN,FFTW_RANK,NFFTW3,NX2*NZ_D,
     &      B,NEMBED,NX2*NZ_D,1,
     &      B,NEMBED,NX2*NZ_D,1,
     &      ISWFFTW,FFTW_ESTIMATE)
      call dfftw_execute(PLAN)
      call dfftw_destroy_plan(PLAN)
#elif MKLFFT
      Status = DftiCreateDescriptor(Desc_Handle,DFTI_DOUBLE,
     &         DFTI_COMPLEX, 1, NY)
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_NUMBER_OF_TRANSFORMS, NX2*NZ_D)
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_INPUT_DISTANCE, 1)
      Stride(1)=0
      Stride(2)=NX2*NZ_D
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_INPUT_STRIDES, Stride)
      Status = DftiCommitDescriptor(Desc_Handle)
!!!!  next line may not be working, because A and B are real array....
      call DCOPY(NY*NX2*NZ_D*2,B,1,C,1)
      Status = DftiComputeBackward(Desc_Handle,C)
      call DCOPY(NY*NX2*NZ_D*2,C,1,B,1)
      Status = DftiFreeDescriptor(Desc_Handle)
#endif
!CDIR PARALLEL DO PRIVATE(I,J,K)
      DO K=0,NZ_D-1
      DO J=0,NY-1
      DO I=0,NX2*2-1
         A(I+J*LX+K*LX*LY)=B(I+K*NX2*2+J*NX2*2*NZ_D)
      ENDDO
      ENDDO
      ENDDO
      IF( IERR(1) .GE. 4000 )THEN
         IERR(1)=4000
         TIME(2)=MPI_WTIME() - T1
         GOTO 999
      ENDIF
      TIME(2)=MPI_WTIME() - T1

! *** X - FFT ( REAL ) *** A( NX, NY, NZ_D )
!     write(*,*) 'cr3d xfft yz'
      T1=MPI_WTIME()
#ifdef ASLFFT
      if (mod(lx,4).eq.0) then
        ll_d=nz_d
        if(mod(nz_d,2).ne.0 .and. (lz_d.gt.nz_d)) then
          ll_d=ll_d+1
        endif
        call transb(A,lx,B,ly*ll_d,nx+2,ly*nz_d,nt)
        CALL QFRMFB(
     &       NX, LY*NZ_D, B, ly*ll_d, 1,
     &       ISW,  IFAX, TRIGS,
     &       A, NT, IERR(1) )
        call transb(B,ly*ll_d,A,lx,ly*nz_d,nx+2,nt)
      else
        CALL QFRMFB(
     &       NX, LY*NZ_D, A, 1, LX,
     &       ISW,  IFAX, TRIGS,
     &       B, NT, IERR(1) )
      endif
#elif SRFFT
      n=NX
      ip=0
      IFCNT=0
      DO I=1,MAXFAX
        IFACT=NFAX(I)
        DO WHILE (MOD(N,IFACT).EQ.0)
          ip(i)=ip(i)+1
          N=N/IFACT
        ENDDO
        IF(N.EQ.1) EXIT
      ENDDO
      iopt(1)=2
      iopt(2)=1
      iopt(3)=1
      nwkw=1
      call hzfd5m(A,LX,LY*NZ_D,NX,ip,LY*NZ_D,iopt,trigs,B,nwkw,ierr(1))
#elif FFTW3
      NFFTW3(1) = NX
      NEMBED(1) = LX2
      NEREAL(1) = LX
      call dfftw_plan_many_dft_c2r
     &     (PLAN,FFTW_RANK,NFFTW3,LY*NZ_D,
     &      A,NEMBED,1,LX2,
     &      A,NEREAL,1,LX,
     &      FFTW_ESTIMATE)
      call dfftw_execute(PLAN)
      call dfftw_destroy_plan(PLAN)
#elif MKLFFT
      Status = DftiCreateDescriptor(Desc_Handle,DFTI_DOUBLE,
     &         DFTI_REAL, 1, NX)
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_NUMBER_OF_TRANSFORMS, LY*NZ_D)
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_INPUT_DISTANCE, LX)
      Status = DftiCommitDescriptor(Desc_Handle)
      Status = DftiComputeBackward(Desc_Handle,A)
      Status = DftiFreeDescriptor(Desc_Handle)
#endif
      TIME(1)=MPI_WTIME() - T1
      IF( IERR(1) .GE. 4000 )THEN
         IERR(1)=4000
         GOTO 999
      ENDIF

      elseif(index.eq.3) then !!! *** Z==>Y ***

! *** Y FFT  ( COMPLEX ) *** A( NX2, NY, NZ_D )
!     write(*,*) 'cr3d yfft zy'
      T1=MPI_WTIME()
!CDIR PARALLEL DO PRIVATE(I,J,K)
      DO K=0,NZ_D-1
      DO J=0,NY-1
      DO I=0,NX2*2-1
         B(I+K*NX2*2+J*NX2*2*NZ_D)=A(I+J*LX+K*LX*LY)
      ENDDO
      ENDDO
      ENDDO
#ifdef ASLFFT
      CALL HFCMFB(NY,NX2*NZ_D,B,NX2*NZ_D,1,ISW,IFAX,TRIGS,A,NT,KERR)
#elif SRFFT
      n=NY
      ip=0
      IFCNT=0
      DO I=1,MAXFAX
        IFACT=NFAX(I)
        DO WHILE (MOD(N,IFACT).EQ.0)
          ip(i)=ip(i)+1
          N=N/IFACT
        ENDDO
        IF(N.EQ.1) EXIT
      ENDDO
      iopt(1)=2
      iopt(2)=1
      iopt(3)=2
      nwkw=1
      call hzft5s(B,NX2*NZ_D,NY,NY,ip,NX2*NZ_D,
     &            iopt,trigs,A,nwkw,IERR(1))
#elif FFTW3
      NFFTW3(1) = NY
      NEMBED(1) = NY
      call dfftw_plan_many_dft
     &     (PLAN,FFTW_RANK,NFFTW3,NX2*NZ_D,
     &      B,NEMBED,NX2*NZ_D,1,
     &      B,NEMBED,NX2*NZ_D,1,
     &      ISWFFTW,FFTW_ESTIMATE)
      call dfftw_execute(PLAN)
      call dfftw_destroy_plan(PLAN)
#elif MKLFFT
      Status = DftiCreateDescriptor(Desc_Handle,DFTI_DOUBLE,
     &         DFTI_COMPLEX, 1, NY)
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_NUMBER_OF_TRANSFORMS, NX2*NZ_D)
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_INPUT_DISTANCE, 1)
      Stride(1)=0
      Stride(2)=NX2*NZ_D
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_INPUT_STRIDES, Stride)
      Status = DftiCommitDescriptor(Desc_Handle)
!!!!  next line may not be working, because A and B are real array....
      call DCOPY(NY*NX2*NZ_D*2,B,1,C,1)
      Status = DftiComputeBackward(Desc_Handle,C)
      call DCOPY(NY*NX2*NZ_D*2,C,1,B,1)
      Status = DftiFreeDescriptor(Desc_Handle)
#endif
!CDIR PARALLEL DO PRIVATE(I,J,K)
      DO K=0,NZ_D-1
      DO J=0,NY-1
      DO I=0,NX2*2-1
         A(I+J*LX+K*LX*LY)=B(I+K*NX2*2+J*NX2*2*NZ_D)
      ENDDO
      ENDDO
      ENDDO
      IF( IERR(1) .GE. 4000 )THEN
         IERR(1)=4000
         GOTO 999
      ENDIF
      TIME(2)=MPI_WTIME() - T1

!     write(*,*) 'cr3d trans zy'

! *** DATA TRANSPOSE *** Z==>Y
      T1=MPI_WTIME()
      CALL FFT_TRANS_ZY(A,LX2,LY,LY_D,NY_D,NZ_D,B,FFT_WORLD,IERR)
      TIME(4)=MPI_WTIME() - T1
      IF( IERR(1) .NE. 0 )THEN
         GOTO 999
      ENDIF

! *** Z - FFT ( COMPLEX ) *** A( NX2, NY_D, NZ )
!     write(*,*) 'cr3d zfft zy'
      T1=MPI_WTIME()
#ifdef ASLFFT
      CALL HFCMFB(
     &     NZ, LX2*NY_D, A, LX2*LY_D, 1,
     &     ISW, IFAX, TRIGS,
     &     B, NT, IERR(1) )
#elif SRFFT
      n=NZ
      ip=0
      IFCNT=0
      DO I=1,MAXFAX
        IFACT=NFAX(I)
        DO WHILE (MOD(N,IFACT).EQ.0)
          ip(i)=ip(i)+1
          N=N/IFACT
        ENDDO
        IF(N.EQ.1) EXIT
      ENDDO
      iopt(1)=2
      iopt(2)=1
      iopt(3)=2
      nwkw=1
      call hzft5s(A,LX2*LY_D,NZ,NZ,ip,LX2*NY_D,
     &            iopt,trigs,B,nwkw,IERR(1))
#elif FFTW3
      NFFTW3(1) = NZ
      NEMBED(1) = LZ
      call dfftw_plan_many_dft
     &     (PLAN,FFTW_RANK,NFFTW3,LX2*NY_D,
     &      A,NEMBED,LX2*LY_D,1,
     &      A,NEMBED,LX2*LY_D,1,
     &      ISWFFTW,FFTW_ESTIMATE)
      call dfftw_execute(PLAN)
      call dfftw_destroy_plan(PLAN)
#elif MKLFFT
      Status = DftiCreateDescriptor(Desc_Handle,DFTI_DOUBLE,
     &         DFTI_COMPLEX, 1, NZ)
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_NUMBER_OF_TRANSFORMS, LX2*NY_D)
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_INPUT_DISTANCE, 1)
      Stride(1)=0
      Stride(2)=LX2*LY_D
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_INPUT_STRIDES, Stride)
      Status = DftiCommitDescriptor(Desc_Handle)
!!!!  next line may not be working, because A and B are real array....
      call DCOPY(LZ*LX2*LY_D*2,A,1,C,1)
      Status = DftiComputeBackward(Desc_Handle,C)
      call DCOPY(LZ*LX2*LY_D*2,C,1,A,1)
      Status = DftiFreeDescriptor(Desc_Handle)
#endif
      TIME(3)=MPI_WTIME() - T1
      IF( IERR(1) .GE. 4000 )THEN
         IERR(1)=4000
         GOTO 999
      ENDIF

! *** X - FFT ( REAL ) *** A( NX, NY_D, NZ )
!     write(*,*) 'cr3d xfft zy'
      T1=MPI_WTIME()
#ifdef ASLFFT
      if (mod(lx,4).eq.0) then
        ll=nz
        if(mod(nz,2).ne.0 .and. (lz.gt.nz)) then
          ll=ll+1
        endif
        call transb(a,lx,b,ly_d*ll,nx+2,ly_d*nz,nt)
        CALL QFRMFB(
     &       NX, LY_D*NZ, B, ly_d*ll, 1,
     &       ISW,  IFAX, TRIGS,
     &       A, NT, IERR(1) )
        call transb(b,ly_d*ll,a,lx,ly_d*nz,nx+2,nt)
      else
        CALL QFRMFB(
     &       NX, LY_D*NZ, A, 1, LX,
     &       ISW,  IFAX, TRIGS,
     &       B, NT, IERR(1) )
      endif
#elif SRFFT
      n=NX
      ip=0
      IFCNT=0
      DO I=1,MAXFAX
        IFACT=NFAX(I)
        DO WHILE (MOD(N,IFACT).EQ.0)
          ip(i)=ip(i)+1
          N=N/IFACT
        ENDDO
        IF(N.EQ.1) EXIT
      ENDDO
      iopt(1)=2
      iopt(2)=1
      iopt(3)=1
      nwkw=1
      call hzfd5m(A,LX,LY_D*NZ,NX,ip,LY_D*NZ,iopt,trigs,B,nwkw,ierr(1))
#elif FFTW3
      NFFTW3(1) = NX
      NEMBED(1) = LX2
      NEREAL(1) = LX
      call dfftw_plan_many_dft_c2r
     &     (PLAN,FFTW_RANK,NFFTW3,LY_D*NZ,
     &      A,NEMBED,1,LX2,
     &      A,NEREAL,1,LX,
     &      FFTW_ESTIMATE)
      call dfftw_execute(PLAN)
      call dfftw_destroy_plan(PLAN)
#elif MKLFFT
      Status = DftiCreateDescriptor(Desc_Handle,DFTI_DOUBLE,
     &         DFTI_REAL, 1, NX)
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_NUMBER_OF_TRANSFORMS, LY_D*NZ)
      Status = DftiSetValue(Desc_Handle,
     &         DFTI_INPUT_DISTANCE, LX)
      Status = DftiCommitDescriptor(Desc_Handle)
      Status = DftiComputeBackward(Desc_Handle,A)
      Status = DftiFreeDescriptor(Desc_Handle)
#endif
      TIME(1)=MPI_WTIME() - T1
      IF( IERR(1) .GE. 4000 )THEN
         IERR(1)=4000
         GOTO 999
      ENDIF

      endif

      endif

      DEALLOCATE( TRIGS, STAT=IERR(2) )
      IF( IERR(2) .NE. 0 )THEN
         IERR(1)=7100
         GOTO 999
      ENDIF

  999 CONTINUE
      TIME(5)=MPI_WTIME() - T0
      RETURN
      END

!     ****************************************************************
!     AUTHORS: T.KOKUBO(NEC), K.TAKADA(NIS), T.FURUTA(NIS), R.SAKA(NIS)
!
!     REGD. ON 03/19/02 BY T.K.
!
      subroutine fftflop(nn,m,flops,ierr)
      IMPLICIT NONE
      INTEGER,PARAMETER :: PREC=8
                                !
      INTEGER,INTENT(IN) :: nn,m
      INTEGER,INTENT(OUT) :: ierr
      REAL(KIND=PREC),INTENT(OUT) :: flops
                                !
      INTEGER :: IFAX(0:63)
      INTEGER :: I,N,IFACT,IFCNT
      INTEGER,PARAMETER :: MAXFAX = 4
      INTEGER,DIMENSION(MAXFAX) :: NFAX = (/2,4,3,5/)
                                !
      n=nn
      ierr=0
                                !
      IFCNT=0
                                !
      DO I=MAXFAX,1,-1
         IFACT=NFAX(I)
         DO WHILE (MOD(N,IFACT).EQ.0)
            IFCNT=IFCNT+1
            IFAX(IFCNT)=IFACT
            N=N/IFACT
         ENDDO
         IF(N.EQ.1) EXIT
      ENDDO
                                !
      IF(N.NE.1) THEN
         ierr=n
         return
      ENDIF
                                !
      IFAX(0)=IFCNT

      if(ifax(1).eq.2) then
         flops=flops+4.d0*dble(nn)*dble(m)/2.d0
      elseif(ifax(1).eq.3) then
         flops=flops+16.d0*dble(nn)*dble(m)/3.d0
      elseif(ifax(1).eq.4) then
         flops=flops+16.d0*dble(nn)*dble(m)/4.d0
      elseif(ifax(1).eq.5) then
         flops=flops+44.d0*dble(nn)*dble(m)/5.d0
      endif
      do i=2,ifax(0)
         if(ifax(i).eq.2) then
            flops=flops+10.d0*dble(nn)*dble(m)/2.d0
         elseif(ifax(i).eq.3) then
            flops=flops+28.d0*dble(nn)*dble(m)/3.d0
         elseif(ifax(i).eq.4) then
            flops=flops+34.d0*dble(nn)*dble(m)/4.d0
         elseif(ifax(i).eq.5) then
            flops=flops+66.d0*dble(nn)*dble(m)/5.d0
         endif
      enddo

      return
      end

!     ****************************************************************
!     AUTHORS: T.KOKUBO(NEC), K.TAKADA(NIS), T.FURUTA(NIS), R.SAKA(NIS)
!
!  REGD. ON 03/15/02 BY T.F.
!  UPDD. ON 08/31/06 BY T.K.
!
      SUBROUTINE FFT_TRANS_YZ(A,lx,ly,ly_d,ny_d,nz_d,
     &                        B,my_comm_world,IERR)
!
      IMPLICIT NONE
      INCLUDE "mpif.h"
!
      INTEGER,PARAMETER :: PREC=8
      INTEGER,INTENT(IN) :: my_comm_world
      INTEGER,INTENT(IN) :: lx,ly_d,ny_d,nz_d,ly
      INTEGER,INTENT(INOUT) :: IERR(2)
      COMPLEX(KIND=PREC),INTENT(INOUT) :: A(0:*)
      COMPLEX(KIND=PREC),INTENT(INOUT) :: B(0:*)
      COMPLEX(KIND=PREC),DIMENSION(:),ALLOCATABLE :: C
      INTEGER :: I,K,LXYZ_D,IS,IDIST
      INTEGER :: ITAG,IRECV,ISEND,IRANK,ISIZE
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER,DIMENSION(:),ALLOCATABLE :: IREQS,IREQR
!
! *** INITIALIZE ***
!
! *** MPI INITIALIZE ***
      CALL MPI_COMM_RANK(MY_COMM_WORLD,IRANK,IERR(2))
      IF( IERR(2) .NE. 0 ) THEN
         IERR(1) = 8000
         RETURN
      ENDIF
      CALL MPI_COMM_SIZE(MY_COMM_WORLD,ISIZE,IERR(2))
      IF( IERR(2) .NE. 0 ) THEN
         IERR(1) = 8000
         RETURN
      ENDIF

      LXYZ_D=LX*LY_D*NZ_D

      ALLOCATE(C(0:LXYZ_D-1),STAT=IERR(2))

      ALLOCATE(IREQS(0:ISIZE-1),STAT=IERR(2))
      IF( IERR(2) .NE. 0 ) THEN
         IERR(1) = 7000
         RETURN
      ENDIF

      ALLOCATE(IREQR(0:ISIZE-1),STAT=IERR(2))
      IF( IERR(2) .NE. 0 ) THEN
         IERR(1) = 7000
         RETURN
      ENDIF

! *** MPI RECIEVE ***
      DO IS=1,ISIZE-1
         IRECV=MOD(ISIZE+IRANK-IS,ISIZE)
         IDIST=LXYZ_D * IRECV
         ITAG =0
         CALL MPI_IRECV
     &        (B(IDIST),LXYZ_D,MPI_DOUBLE_COMPLEX,
     &         IRECV,ITAG,MY_COMM_WORLD,IREQR(IS),IERR(2))
         IF( IERR(2) .NE. 0 ) THEN
            IERR(1) = 8000
            RETURN
         ENDIF
      ENDDO

      CALL MPI_BARRIER(MY_COMM_WORLD,IERR(2))
      IF( IERR(2) .NE. 0 ) THEN
         IERR(1) = 8000
         RETURN
      ENDIF

! *** MPI SEND ***
      DO IS=1,ISIZE-1
         ISEND=MOD(IRANK+IS,ISIZE)
         IDIST=LXYZ_D * ISEND
         ITAG = 0
         CALL MPI_ISEND
     &        (A(IDIST),LXYZ_D,MPI_DOUBLE_COMPLEX,
     &         ISEND,ITAG,MY_COMM_WORLD,IREQS(IS),IERR(2))
         IF( IERR(2) .NE. 0 ) THEN
            IERR(1) = 8000
            RETURN
         ENDIF
      ENDDO

! *** DATA COPY IN MEMORY ON SAME RANK  ***
!CDIR PARALLEL DO PRIVATE(I,K)
      DO K=0,NZ_D-1
      DO I=0,LX*NY_D-1
         C(I+K*LX*LY_D) = A(I+K*LX*LY_D+LXYZ_D*IRANK)
      ENDDO
      ENDDO

! *** MPI WAIT ***
      DO IS=1,ISIZE-1
         CALL MPI_WAIT(IREQS(IS), STATUS, IERR(2))
         IF( IERR(2) .NE. 0 ) THEN
            IERR(1) = 8000
            RETURN
         ENDIF
         CALL MPI_WAIT(IREQR(IS), STATUS, IERR(2))
         IF( IERR(2) .NE. 0 ) THEN
            IERR(1) = 8000
            RETURN
         ENDIF
      ENDDO

! *** DATA COPY ***
!CDIR PARALLEL DO PRIVATE(I,K,IS,IRECV)
         DO K=0,NZ_D-1
      DO IS=1,ISIZE-1
         IRECV=MOD(ISIZE+IRANK-IS,ISIZE)
         DO I=0,LX*NY_D-1
            A(I+K*LX*LY+LX*NY_D*IRECV) = B(I+K*LX*LY_D+LXYZ_D*IRECV)
         ENDDO
      ENDDO
         ENDDO

! *** DATA COPY IN MEMORY ON SAME RANK  ***
!CDIR PARALLEL DO PRIVATE(I,K)
      DO K=0,NZ_D-1
      DO I=0,LX*NY_D-1
         A(I+K*LX*LY+LX*NY_D*IRANK) = C(I+K*LX*LY_D)
      ENDDO
      ENDDO
!edit_s(end  : sendrecv NO TYPE)------------------------------
 
      DEALLOCATE(C,STAT=IERR(2))
      DEALLOCATE(IREQS,STAT=IERR(2))
      IF( IERR(2) .NE. 0 ) THEN
         IERR(1) = 7100
         RETURN
      ENDIF
      DEALLOCATE(IREQR,STAT=IERR(2))
      IF( IERR(2) .NE. 0 ) THEN
         IERR(1) = 7100
         RETURN
      ENDIF
!     write(*,*) 'done fft_trans_yz'
!     write(*,*) 'done fft_trans_yz',irank
!     CALL MPI_BCAST(i,1,mpi_integer,0,MY_COMM_WORLD,IERR(2))
      RETURN
      END SUBROUTINE FFT_TRANS_YZ
!
      SUBROUTINE FFT_TRANS_ZY(A,lx,ly,ly_d,ny_d,nz_d,
     &                        B,my_comm_world,IERR)
!
      IMPLICIT NONE
      INCLUDE "mpif.h"
!
      INTEGER,PARAMETER :: PREC=8
      INTEGER,INTENT(IN) :: my_comm_world
      INTEGER,INTENT(IN) :: lx,ly_d,ny_d,nz_d,ly
      INTEGER,INTENT(INOUT) :: IERR(2)
      COMPLEX(KIND=PREC),INTENT(INOUT) :: A(0:*)
      COMPLEX(KIND=PREC),INTENT(INOUT) :: B(0:*)
      COMPLEX(KIND=PREC),DIMENSION(:),ALLOCATABLE :: C
      INTEGER :: I,K,LXYZ_D,IS,IDIST
      INTEGER :: ITAG,IRECV,ISEND,IRANK,ISIZE
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER,DIMENSION(:),ALLOCATABLE :: IREQS,IREQR
!
! *** INITIALIZE ***
!
! *** MPI INITIALIZE ***
      CALL MPI_COMM_RANK(MY_COMM_WORLD,IRANK,IERR(2))
      IF( IERR(2) .NE. 0 ) THEN
         IERR(1) = 8000
         RETURN
      ENDIF
      CALL MPI_COMM_SIZE(MY_COMM_WORLD,ISIZE,IERR(2))
      IF( IERR(2) .NE. 0 ) THEN
         IERR(1) = 8000
         RETURN
      ENDIF

      LXYZ_D=LX*LY_D*NZ_D
      ALLOCATE(C(0:LXYZ_D-1),STAT=IERR(2))

      ALLOCATE(IREQS(0:ISIZE-1),STAT=IERR(2))
      IF( IERR(2) .NE. 0 ) THEN
         IERR(1) = 7000
         RETURN
      ENDIF

      ALLOCATE(IREQR(0:ISIZE-1),STAT=IERR(2))
      IF( IERR(2) .NE. 0 ) THEN
         IERR(1) = 7000
         RETURN
      ENDIF

!CDIR PARALLEL DO PRIVATE(I,K,IS,ISEND)
         DO K=0,NZ_D-1
      DO IS=1,ISIZE-1
         ISEND=MOD(IRANK+IS,ISIZE)
         DO I=0,LX*NY_D-1
            B(I+K*LX*LY_D+LXYZ_D*ISEND) = A(I+K*LX*LY+LX*NY_D*ISEND)
         ENDDO
      ENDDO
         ENDDO

! *** DATA COPY ***
!CDIR PARALLEL DO PRIVATE(I,K)
      DO K=0,NZ_D-1
      DO I=0,LX*NY_D-1
         C(I+K*LX*LY_D) = A(I+K*LX*LY+LX*NY_D*IRANK)
      ENDDO
      ENDDO

! *** MPI RECIEVE ***
      DO IS=1,ISIZE-1
         IRECV=MOD(ISIZE+IRANK-IS,ISIZE)
         IDIST=LXYZ_D*IRECV
         ITAG =0
         CALL MPI_IRECV
     &        (A(IDIST),LXYZ_D,MPI_DOUBLE_COMPLEX,
     &         IRECV,ITAG,MY_COMM_WORLD,IREQR(IS),IERR(2))
         IF( IERR(2) .NE. 0 ) THEN
            IERR(1) = 8000
            RETURN
         ENDIF
      ENDDO
      CALL MPI_BARRIER(MY_COMM_WORLD,IERR(2))

! *** MPI SEND ***
      IF( IERR(2) .NE. 0 ) THEN
         IERR(1) = 8000
         RETURN
      ENDIF

      DO IS=1,ISIZE-1
         ISEND=MOD(IRANK+IS,ISIZE)
         IDIST=LXYZ_D*ISEND
         ITAG =0
         CALL MPI_ISEND
     &        (B(IDIST),LXYZ_D,MPI_DOUBLE_COMPLEX,
     &         ISEND,ITAG,MY_COMM_WORLD,IREQS(IS),IERR(2))
!        CALL MPI_SEND
!    &        (B(IDIST),LXYZ_D,MPI_DOUBLE_COMPLEX,
!    &         ISEND,ITAG,MY_COMM_WORLD,IERR(2))
         IF( IERR(2) .NE. 0 ) THEN
            IERR(1) = 8000
            RETURN
         ENDIF
      ENDDO

! *** MPI WAIT ***
      DO IS=1,ISIZE-1
         CALL MPI_WAIT(IREQS(IS), STATUS, IERR(2))
         IF( IERR(2) .NE. 0 ) THEN
            IERR(1) = 8000
            RETURN
         ENDIF
         CALL MPI_WAIT(IREQR(IS), STATUS, IERR(2))
         IF( IERR(2) .NE. 0 ) THEN
            IERR(1) = 8000
            RETURN
         ENDIF
      ENDDO

! *** DATA COPY ***
!!CDIR PARALLEL DO PRIVATE(I,K)
      DO K=0,NZ_D-1
      DO I=0,LX*NY_D-1
         A(I+K*LX*LY_D+LXYZ_D*IRANK) = C(I+K*LX*LY_D) 
      ENDDO
      ENDDO
!edit_s(end  : sendrecv NO TYPE)------------------------------
      DEALLOCATE(C,STAT=IERR(2))
      DEALLOCATE(IREQS,STAT=IERR(2))
      IF( IERR(2) .NE. 0 ) THEN
         IERR(1) = 7100
         RETURN
      ENDIF
      DEALLOCATE(IREQR,STAT=IERR(2))
      IF( IERR(2) .NE. 0 ) THEN
         IERR(1) = 7100
         RETURN
      ENDIF
!     write(*,*) 'done fft_trans_zy',irank
      RETURN
      END SUBROUTINE FFT_TRANS_ZY
      subroutine transb(a,la,b,lb,n1,n2,nt)
      INTEGER,PARAMETER :: PREC=8
      INTEGER,INTENT(IN) :: LA,LB,N1,N2,NT
      REAL(KIND=PREC),INTENT(INOUT) :: A(LA,N2),B(LB,N1)
      INTEGER :: I,J
      IF (N2 .GE. N1) THEN
!cdir pardo for=NT
        DO 70 I=0,N1-1
!cdir noloopchg
          DO 60 J=1,N1-I
            B(J,I+J)=A(I+J,J)
   60     CONTINUE
   70   CONTINUE
!cdir pardo for=NT
        DO 90 I=1,N2-N1
!cdir noloopchg
          DO 80 J=1,N1
            B(I+J,J)=A(J,I+J)
   80     CONTINUE
   90   CONTINUE
!cdir pardo for=NT
        DO 110 I=N2-N1+1,N2-1
!cdir noloopchg
          DO 100 J=1,N2-I
            B(I+J,J)=A(J,I+J)
  100     CONTINUE
  110   CONTINUE
      ELSE
!cdir pardo for=NT
        DO 130 I=0,N2-1
!cdir noloopchg
          DO 120 J=1,N2-I
            B(I+J,J)=A(J,I+J)
  120     CONTINUE
  130   CONTINUE
!cdir pardo for=NT
        DO 150 I=1,N1-N2
!cdir noloopchg
          DO 140 J=1,N2
            B(J,I+J)=A(I+J,J)
  140     CONTINUE
  150   CONTINUE
!cdir pardo for=NT
        DO 170 I=N1-N2+1,N1-1
!cdir noloopchg
          DO 160 J=1,N1-I
            B(J,I+J)=A(I+J,J)
  160     CONTINUE
  170   CONTINUE
      END IF
      RETURN
      END subroutine transb
      subroutine ztransb(a,la,b,lb,n1,n2,nt)
      INTEGER,PARAMETER :: PREC=8
      INTEGER,INTENT(IN) :: LA,LB,N1,N2,NT
      COMPLEX(KIND=PREC),INTENT(INOUT) :: A(LA,N2),B(LB,N1)
      INTEGER :: I,J
      IF (N2 .GE. N1) THEN
!cdir pardo for=NT
        DO 70 I=0,N1-1
!cdir noloopchg
          DO 60 J=1,N1-I
            B(J,I+J)=A(I+J,J)
   60     CONTINUE
   70   CONTINUE
!cdir pardo for=NT
        DO 90 I=1,N2-N1
!cdir noloopchg
          DO 80 J=1,N1
            B(I+J,J)=A(J,I+J)
   80     CONTINUE
   90   CONTINUE
!cdir pardo for=NT
        DO 110 I=N2-N1+1,N2-1
!cdir noloopchg
          DO 100 J=1,N2-I
            B(I+J,J)=A(J,I+J)
  100     CONTINUE
  110   CONTINUE
      ELSE
!cdir pardo for=NT
        DO 130 I=0,N2-1
!cdir noloopchg
          DO 120 J=1,N2-I
            B(I+J,J)=A(J,I+J)
  120     CONTINUE
  130   CONTINUE
!cdir pardo for=NT
        DO 150 I=1,N1-N2
!cdir noloopchg
          DO 140 J=1,N2
            B(J,I+J)=A(I+J,J)
  140     CONTINUE
  150   CONTINUE
!cdir pardo for=NT
        DO 170 I=N1-N2+1,N1-1
!cdir noloopchg
          DO 160 J=1,N1-I
            B(J,I+J)=A(I+J,J)
  160     CONTINUE
  170   CONTINUE
      END IF
      RETURN
      END subroutine ztransb

!     ****************************************************************
!     AUTHORS: T.KOKUBO(NEC), K.TAKADA(NIS), T.FURUTA(NIS), R.SAKA(NIS)
!
!     REGD. ON 02/18/02 BY K.T.
!     UPDD. ON 02/18/02 BY K.T.
!     UPDD. ON 02/25/02 BY K.T.
!     UPDD. ON 03/12/02 BY K.T.
!
!      SUBROUTINE MKDATA_ZFC3FB
      SUBROUTINE MKD_C3D
     &(A,LX,LY,NX,NY,NZ,ID_IN,N_IN,ISW,PARAM,IRANK,NPROC,IERR)
!
      USE m_Const_Parameters, only : PAI
      IMPLICIT NONE
      INTEGER,PARAMETER :: PREC=8
!
! ARGUMENTS
!
      INTEGER,INTENT(IN) :: LX,LY,NX,NY,NZ,ID_IN,N_IN,ISW,IRANK,NPROC
      INTEGER,INTENT(INOUT) :: IERR
      REAL(KIND=PREC),INTENT(IN) :: PARAM(*)
      COMPLEX(KIND=PREC),INTENT(OUT) :: A(*)
!
! FUNCTION:
!
! ALGORITHM:
!
! ARGUMENTS:
!
      REAL(KIND=PREC) :: SW1(20)
      REAL(KIND=PREC),PARAMETER :: PI=PAI
      INTEGER :: NW,I,J,K,II,JJ,KK,NN(3),NSTART(3),NEND(3)
      REAL(KIND=PREC) :: ARG,ARGI,ARGJ,ARGK
      REAL(KIND=PREC),PARAMETER :: C0=0.0D0,C1=1.0D0,C2=2.0D0
      COMPLEX(KIND=PREC) :: WZI,WZJ,WZK,BP,BQ,BR
!F-debug      COMPLEX(KIND=PREC),PARAMETER :: ZC1=(C1,C0)
      COMPLEX(KIND=PREC),PARAMETER :: ZC1=(1.0D0,0.0D0)
      INTEGER :: LXY,IS,JS,KS,IFLAG,JFLAG,KFLAG
!
      IERR = 0
!
      NN(1)=NX
      NN(2)=NY
      NN(3)=NZ
      IF(NN(ID_IN).LT.N_IN) THEN
        IERR=3000
        RETURN
      ENDIF
      IF(N_IN.LE.0) THEN
        IERR=3010
        RETURN
      ENDIF
      IF((IRANK.LT.0).OR.(IRANK.GE.NPROC)) THEN
        IERR=3100
        RETURN
      ENDIF
      IF((ID_IN.LT.1).OR.(ID_IN.GT.3)) THEN
        IERR=3110
        RETURN
      ENDIF
      IF((ISW.LT.-1).OR.(ISW.GT.1)) THEN
        IERR=3200
        RETURN
      ENDIF
!
      NN(1) = NX
      NN(2) = NY
      NN(3) = NZ
      NW=NN(ID_IN)
      NN(ID_IN)=N_IN
      NSTART(1)=1
      NSTART(2)=1
      NSTART(3)=1
      NEND(1)=NN(1)
      NEND(2)=NN(2)
      NEND(3)=NN(3)
      NSTART(ID_IN)=IRANK*NN(ID_IN)+1
      NEND(ID_IN)=MIN(NW,(IRANK+1)*NN(ID_IN))
      LXY=LX*LY
      IF(ISW.EQ.1) THEN
        SW1(1)=-PI*PARAM(1)/REAL(NX,PREC)
        SW1(2)=-PI*PARAM(2)/REAL(NY,PREC)
        SW1(3)=-PI*PARAM(3)/REAL(NZ,PREC)
        DO I=NSTART(1),NEND(1)
          II=I-NSTART(1)+1
          DO J=NSTART(2),NEND(2)
            JJ=J-NSTART(2)+1
            DO K=NSTART(3),NEND(3)
              KK=K-NSTART(3)+1
              ARG=SW1(1)*REAL(I-1,PREC)
     &            +SW1(2)*REAL(J-1,PREC)
     &            +SW1(3)*REAL(K-1,PREC)
              A(II+(JJ-1)*LX+(KK-1)*LXY)=CMPLX(COS(ARG),SIN(ARG))
     &        /(REAL(NX,PREC)*REAL(NY,PREC)*REAL(NZ,PREC))
            ENDDO
          ENDDO
        ENDDO
      ELSE IF(ISW.EQ.-1) THEN
        IF(ABS(PARAM(1)-REAL(INT(PARAM(1)),PREC)).LE.PARAM(4)) THEN
          IFLAG=1
          IS=PARAM(1)
        ELSE
          IFLAG=0
        ENDIF
        IF(ABS(PARAM(2)-REAL(INT(PARAM(2)),PREC)).LE.PARAM(4)) THEN
          JFLAG=1
          JS=PARAM(2)
        ELSE
          JFLAG=0
        ENDIF
        IF(ABS(PARAM(3)-REAL(INT(PARAM(3)),PREC)).LE.PARAM(4)) THEN
          KFLAG=1
          KS=PARAM(3)
        ELSE
          KFLAG=0
        ENDIF
        DO I=NSTART(1),NEND(1)
          II=I-NSTART(1)+1
          ARGI=-PI*(C2*REAL(I-1,PREC)+PARAM(1))
          WZI=ZC1-EXP(CMPLX(C0,ARGI/REAL(NX,PREC)))
          IF((IFLAG.EQ.1).AND.(MOD(IS+2*(I-1),2*NX).EQ.0)) THEN
            BP=REAL(NX,PREC)
          ELSE
            BP=(ZC1-EXP(CMPLX(C0,ARGI,PREC)))/WZI
          ENDIF
          DO J=NSTART(2),NEND(2)
            JJ=J-NSTART(2)+1
            ARGJ=-PI*(C2*REAL(J-1,PREC)+PARAM(2))
            WZJ=ZC1-EXP(CMPLX(C0,ARGJ/REAL(NY,PREC)))
            IF((JFLAG.EQ.1).AND.(MOD(JS+2*(J-1),2*NY).EQ.0)) THEN
              BQ=REAL(NY,PREC)
            ELSE
              BQ=(ZC1-EXP(CMPLX(C0,ARGJ,PREC)))/WZJ
            ENDIF
            DO K=NSTART(3),NEND(3)
              KK=K-NSTART(3)+1
              ARGK=-PI*(C2*REAL(K-1,PREC)+PARAM(3))
              WZK=ZC1-EXP(CMPLX(C0,ARGK/REAL(NZ,PREC)))
              IF((KFLAG.EQ.1).AND.(MOD(KS+2*(K-1),2*NZ).EQ.0)) THEN
                BR=REAL(NZ,PREC)
              ELSE
                BR=(ZC1-EXP(CMPLX(C0,ARGK,PREC)))/WZK
              ENDIF
              A(II+(JJ-1)*LX+(KK-1)*LXY)=BP*BQ*BR
     &        /(REAL(NX,PREC)*REAL(NY,PREC)*REAL(NZ,PREC))
            ENDDO
          ENDDO
        ENDDO
      ENDIF
!
      RETURN
      END

!     ****************************************************************
!     AUTHORS: T.KOKUBO(NEC), K.TAKADA(NIS), T.FURUTA(NIS), R.SAKA(NIS)
!
!     REGD. ON 03/01/02 BY K.T.
!     UPDD. ON 03/12/02 BY K.T.
!     UPDD. ON 03/14/02 BY K.T.
!
!      SUBROUTINE MKDATA_DFR3FB
      SUBROUTINE MKD_R3D
     &(A,LXW,LY,NX,NY,NZ,ID_IN,N_INW,ISW,PARAM,IRANK,NPROC,IERR)
!
      USE m_Const_Parameters, only : PAI
      IMPLICIT  NONE
      INTEGER,PARAMETER :: PREC=8
!
! ARGUMENTS
!
      INTEGER,INTENT(IN) :: LXW,LY,NX,NY,NZ,ID_IN,N_INW
      INTEGER,INTENT(IN) :: ISW,IRANK,NPROC
      INTEGER,INTENT(INOUT) :: IERR
      REAL(KIND=PREC),INTENT(IN) :: PARAM(*)
      REAL(KIND=PREC),INTENT(OUT) :: A(*)
!
! FUNCTION:
!
! ALGORITHM:
!
! ARGUMENTS:
!
      REAL(KIND=PREC) :: ARGX,ARGY,ARG,SW1(20)
      REAL(KIND=PREC),PARAMETER :: C0=0.0D0,C1=1.0D0,C2=2.0D0
      COMPLEX(KIND=PREC) :: WZI,WZI2,WZJ,WZJ2,WZK,WZK2,BP,BQ,BR
!F-debug      COMPLEX(KIND=PREC),PARAMETER :: ZC1=(C1,C0)
      COMPLEX(KIND=PREC),PARAMETER :: ZC1=(1.0D0,0.0D0)
      REAL(KIND=PREC),PARAMETER :: PI=PAI
      INTEGER :: LXY,NW,I,J,K,II,JJ,KK,NN(3),NSTART(3),NEND(3),LX
      INTEGER :: N_IN
!
      IERR = 0
!
      NN(1)=NX
      NN(2)=NY
      NN(3)=NZ
      N_IN=N_INW
      IF(ISW.GT.0) THEN
        LX=LXW
      ELSE
        LX=2*LXW
!       IF(ID_IN.EQ.1) THEN
!         N_IN=2*N_INW
!       ENDIF
      ENDIF
      IF(NN(ID_IN).LT.N_IN) THEN
        IERR=3000
        RETURN
      ENDIF
      IF(N_IN.LE.0) THEN
        IERR=3010
        RETURN
      ENDIF
      IF((IRANK.LT.0).OR.(IRANK.GE.NPROC)) THEN
        IERR=3100
        RETURN
      ENDIF
      IF((ID_IN.LT.1).OR.(ID_IN.GT.3)) THEN
        IERR=3110
        RETURN
      ENDIF
      IF((ISW.LT.-1).OR.(ISW.GT.1)) THEN
        IERR=3200
        RETURN
      ENDIF
!
      IF(ISW.GT.0) THEN
        NN(1) = NX
      ELSE
        NN(1) = (NX+2)/2
      ENDIF
      NN(2) = NY
      NN(3) = NZ
      NW=NN(ID_IN)
      NN(ID_IN)=N_IN
      NSTART(1)=1
      NSTART(2)=1
      NSTART(3)=1
      NEND(1)=NN(1)
      NEND(2)=NN(2)
      NEND(3)=NN(3)
      NSTART(ID_IN)=IRANK*NN(ID_IN)+1
      NEND(ID_IN)=MIN(NW,(IRANK+1)*NN(ID_IN))
      LXY=LX*LY
      IF(ISW.EQ.1) THEN
        SW1(1)=PARAM(1)/REAL(NX,PREC)
        SW1(2)=PARAM(2)/REAL(NY,PREC)
        SW1(3)=PARAM(3)/REAL(NZ,PREC)
        DO I=NSTART(1),NEND(1)
          II=I-NSTART(1)+1
          DO J=NSTART(2),NEND(2)
            JJ=J-NSTART(2)+1
            DO K=NSTART(3),NEND(3)
              KK=K-NSTART(3)+1
              ARG=SW1(1)*REAL(I-1,PREC)
     &            +SW1(2)*REAL(J-1,PREC)
     &            +SW1(3)*REAL(K-1,PREC)
              A(II+(JJ-1)*LX+(KK-1)*LXY)=EXP(ARG)
     &        /(REAL(NX,PREC)*REAL(NY,PREC)*REAL(NZ,PREC))
            ENDDO
          ENDDO
        ENDDO
      ELSE IF(ISW.EQ.-1) THEN
        DO I=NSTART(1),NEND(1)
          II=I-NSTART(1)+1
          ARGX=PARAM(1)
          ARGY=-PI*C2*REAL(I-1,PREC)
          WZI=ZC1-EXP(ARGX/REAL(NX,PREC))*
     &        CMPLX(COS(ARGY/REAL(NX,PREC)),SIN(ARGY/REAL(NX,PREC)))
          WZI2=ZC1-EXP(ARGX)*CMPLX(COS(ARGY),SIN(ARGY))
          IF(ABS(WZI).LT.PARAM(4)) THEN
            BP=ZC1
          ELSE
            BP=WZI2/(WZI*REAL(NX,PREC))
          ENDIF
          DO J=NSTART(2),NEND(2)
            JJ=J-NSTART(2)+1
            ARGX=PARAM(2)
            ARGY=-PI*C2*REAL(J-1,PREC)
            WZJ=ZC1-EXP(ARGX/REAL(NY,PREC))*
     &        CMPLX(COS(ARGY/REAL(NY,PREC)),SIN(ARGY/REAL(NY,PREC)))
            WZJ2=ZC1-EXP(ARGX)*CMPLX(COS(ARGY),SIN(ARGY))
            IF(ABS(WZJ).LT.PARAM(4)) THEN
              BQ=ZC1
            ELSE
              BQ=WZJ2/(WZJ*REAL(NY,PREC))
            ENDIF
            DO K=NSTART(3),NEND(3)
              KK=K-NSTART(3)+1
              ARGX=PARAM(3)
              ARGY=-PI*C2*REAL(K-1,PREC)
              WZK=ZC1-EXP(ARGX/REAL(NZ,PREC))*
     &        CMPLX(COS(ARGY/REAL(NZ,PREC)),SIN(ARGY/REAL(NZ,PREC)))
              WZK2=ZC1-EXP(ARGX)*CMPLX(COS(ARGY),SIN(ARGY))
              IF(ABS(WZK).LT.PARAM(4)) THEN
                BR=ZC1
              ELSE
                BR=WZK2/(WZK*REAL(NZ,PREC))
              ENDIF
              A(2*II-1+(JJ-1)*LX+(KK-1)*LXY)=REAL(BP*BQ*BR,PREC)
              A(2*II+(JJ-1)*LX+(KK-1)*LXY)=DIMAG(BP*BQ*BR)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
!
      RETURN
      END
