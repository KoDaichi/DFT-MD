      SUBROUTINE LAPACK_EIGEN2(N_GLOBAL, N, HBW, ID, D, E, LDE, Q, LDQ)
      implicit NONE

      INTEGER, intent(IN)    :: N_GLOBAL, N, HBW, ID, LDE, LDQ
      REAL(8), intent(INOUT) :: D(*), E(LDE,*), Q(LDQ,*)

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
         Q(I,I)=D(I)
      ENDDO

      DO J=1,HBW
      DO I=1,N-J
         Q(I,I+J)=E(I,J)
         Q(I+J,I)=E(I,J)
      ENDDO
      ENDDO

      JOBU = 'V'; JOBVT= 'U'

      lwork = -1
      liwork = -1

      call dsyevd(JOBU, JOBVT, N, Q, LDQ, D(1),
     &            temp, lwork, i, liwork, info)

      lwork  = int(temp)
      liwork = i;

      allocate ( work(lwork), iwork(liwork) )

      call dsyevd(JOBU, JOBVT, N, Q, LDQ, D(1),
     &             work, lwork, iwork, liwork, info)

      deallocate ( work, iwork )

      END SUBROUTINE LAPACK_EIGEN2

