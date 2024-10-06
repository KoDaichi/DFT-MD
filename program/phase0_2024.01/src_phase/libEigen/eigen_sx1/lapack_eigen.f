!#include "f2c.h"
!#include <malloc.h>
!
!
!extern void
!dsyevd_(char *jobz, char *uplo, integer *n, doublereal *
!a, integer *lda, doublereal *w, doublereal *work, integer *lwork, 
!integer *iwork, integer *liwork, integer *info);
!
!int lapack_eigen2_(int *N, int *n, int *hbw, int *id, double *D, double *E, int *lde, double *Q, int *LDQ) {

      SUBROUTINE LAPACK_EIGEN2(N_GLOBAL, N, HBW, ID, D, E, LDE, Q, LDQ)
      implicit NONE

      INTEGER :: N_GLOBAL, N, HBW, ID, LDE, LDQ
      REAL(8) :: D(*), E(LDE,*), Q(LDQ,*)

      INTEGER :: i, j, info, lwork, liwork
      REAL(8) :: temp
      REAL(8), POINTER :: work(:)
      INTEGER, POINTER :: iwork(:)
      CHARACTER*1 :: JOBU, JOBVT

      DO J=1,N
      DO I=1,N
         Q(I,J)=0.0D0
      ENDDO
      ENDDO

      DO I=1,N
         Q(I,I)=D(I+ID-1)
      ENDDO

      DO J=1,HBW
      DO I=1,N-J
         Q(I,I+J)=E(I+ID-1,J)
         Q(I+J,I)=E(I+ID-1,J)
      ENDDO
      ENDDO

      JOBU = 'V'; JOBVT= 'U'

      lwork = -1
      liwork = -1

      call dsyevd(JOBU, JOBVT, N, Q, LDQ, D(1+ID-1),
     &            temp, lwork, i, liwork, info)
!      call DSYEV (JOBU, JOBVT, N, Q, LDQ, D(1+ID-1),
!     &            temp, lwork, info)

      lwork  = int(temp)
      liwork = i;

      allocate(work(lwork), iwork(liwork))
!      allocate(work(lwork))

      call dsyevd(JOBU, JOBVT, N, Q, LDQ, D(1+ID-1),
     &             work, lwork, iwork, liwork, info)
!      call DSYEV (JOBU, JOBVT, N, Q, LDQ, D(1+ID-1),
!     &            work, lwork, info)

      deallocate(work,iwork)
!      deallocate(work)

      END SUBROUTINE

