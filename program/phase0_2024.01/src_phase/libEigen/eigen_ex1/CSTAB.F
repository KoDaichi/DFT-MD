       subroutine CSTAB_get_optdim(n_min, n_unroll,
     &                             delta_L1, delta_L2, n_opt)
       implicit none

       integer                :: n_min, n_unroll, delta_L1, delta_L2
       integer                :: n_opt, n_delta
       integer                :: n_opt2, n_delta2

       include 'CSTAB.h'

       integer                :: i,j,k


       n_opt=n_min
       do

          n_opt   = (n_opt-1)/L1_WINDOW+1
          n_opt   = (n_opt/2)*2+1
          n_opt   = n_opt*L1_WINDOW

          n_delta = 0

          do i=1,INT((n_unroll*1.2-1.0)/L1_WAY+1)

             k=mod(i*n_opt+L1_LSIZE/2,L1_LSIZE)-L1_LSIZE/2
             if(abs(k)<=delta_L1/2)then
                n_delta=(delta_L1/2-k-1)/i+1
                goto 10000
             end if

          end do

          do i=1,INT((n_unroll*1.2-1.0)/L2_WAY+1)

             k=mod(i*n_opt+L2_LSIZE/2,L2_LSIZE)-L2_LSIZE/2
             if(abs(k)<=delta_L2/2)then
                n_delta=(delta_L2/2-k-1)/i+1
                goto 10000
             end if

          end do

10000     continue

          if(n_delta==0)exit
          n_opt = n_opt + n_delta

       end do

       return
       end subroutine

       subroutine CSTAB_adjust_base(a,b,offset)
       implicit none

       real(8)                :: a(*), b(*)
       integer                :: offset

       include 'CSTAB.h'


       call get_delta(a(1),b(1),offset)
       offset=(offset/8)

       if(offset>0)then
         offset=MOD(L2_LSIZE-MOD(+offset,L2_LSIZE),L2_LSIZE)
       else
         offset=MOD(L2_LSIZE+MOD(-offset,L2_LSIZE),L2_LSIZE)
       endif

       return
       end subroutine
!
       subroutine CSTAB_adjust_page(a,b,offset)
       implicit none

       real(8)                :: a(*), b(*)
       integer                :: offset

       include 'CSTAB.h'


       call get_delta(a(1),b(1),offset)
       offset=(offset/8)
       if(offset>0)then
         offset=MOD(PAGE_LSIZE-MOD(+offset,PAGE_LSIZE),PAGE_LSIZE)
       else
         offset=MOD(PAGE_LSIZE+MOD(-offset,PAGE_LSIZE),PAGE_LSIZE)
       endif

       return
       end subroutine
!
       subroutine CSTAB_round_offset(offset)
       implicit none

       integer                :: offset

       include 'CSTAB.h'

       offset=MOD(offset,L2_LSIZE)

       return
       end subroutine
!
       subroutine CSTAB_pre_format(nx,a0,nm0,a,nm,n_,a_,nm_)
       implicit none
*
*   A0(NM0,NX) -> A_(NM_,N_) + A(NM,NX-N_)
*
       integer                :: nx, nm0, nm, n_, nm_
       real(8)                :: a0(nm0,*)
       real(8)                :: a(nm,*), a_(nm_,*)

       include 'CSTAB.h'

       integer, parameter     :: n_thread=8
       integer                :: i_thread

       integer, parameter     :: LX = (L2_SIZE/8/4)
       real(8), pointer       :: buff(:)

       allocate(buff(LX*n_thread))
       call CSTAB_pre_format_body(nx,a0,nm0,a,nm,n_,a_,nm_,buff)
       deallocate(buff)

       return
       end subroutine

       subroutine CSTAB_pre_format_body(nx,a0,nm0,a,nm,n_,a_,nm_,buff)
       implicit none
*
*   A0(NM0,NX) -> A_(NM_,N_) + A(NM,NX-N_)
*
       integer                :: nx, nm0, nm, n_, nm_
       real(8)                :: a0(nm0,*)
       real(8)                :: a(nm,*), a_(nm_,*)

       include 'CSTAB.h'

       integer, parameter     :: n_thread=8
       integer                :: i_thread

       integer, parameter     :: LX = (L2_SIZE/8/4)
       real(8)                :: buff(LX,n_thread)

       integer                ::  i_1, i_2, i_3
       integer                ::  j_1, j_2, j_3
       integer                ::  k_1, k_2, k_3, k_4
       integer                ::  kk_1
       integer                ::  MX, NN

       NN=MAX(1,MIN(nm0,nm))
       MX=(NN-1)/LX+1

       do  i_1=1,n_
          j_2=1
          j_3=MIN(nm0,nm_)
*voption vec
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          do  j_1=j_2,j_3
             buff(j_1,1)=a0(j_1,i_1)
          end do! j_1
*voption vec
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          do  j_1=j_2,j_3
             a_(j_1,i_1)=buff(j_1,1)
          end do! j_1
       end do! i_1

       k_2=1
       k_3=nx*MX
       k_4=MOD(k_3-k_2+1,n_thread)+k_2
       do k_1=k_2,k_4-1
             i_1=(k_1-1)/MX+1
             j_2=MOD(k_1-1,MX)*LX+1
             j_3=MIN(j_2+LX-1,NN)
*voption vec
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
             do  j_1=j_2,j_3
                a(j_1,i_1)
     &          =a0(j_1,i_1+n_)
             end do! j_1
       end do! k_1

       do kk_1=k_4,k_3,n_thread

!$OMP PARALLEL DO PRIVATE(k_1,i_1,j_1,j_2,j_3)
          do  i_thread=1,n_thread
             k_1=kk_1+i_thread-1
             i_1=(k_1-1)/MX+1
             j_2=MOD(k_1-1,MX)*LX+1
             j_3=MIN(j_2+LX-1,NN)
*voption vec
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
             do  j_1=j_2,j_3
                buff(j_1-j_2+1,i_thread)
     &          =a0(j_1,i_1+n_)
             end do! j_1
          end do! i_thread

!$OMP PARALLEL DO PRIVATE(k_1,i_1,j_1,j_2,j_3)
          do  i_thread=1,n_thread
             k_1=kk_1+i_thread-1
             i_1=(k_1-1)/MX+1
             j_2=MOD(k_1-1,MX)*LX+1
             j_3=MIN(j_2+LX-1,NN)
*voption vec
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
             do  j_1=j_2,j_3
                a(j_1,i_1)
     &          =buff(j_1-j_2+1,i_thread)
             end do! j_1
          end do! i_thread

       end do! k_1

       return
       end subroutine

       subroutine CSTAB_post_format(nx,a0,nm0,a,nm,n_,a_,nm_)
       implicit none
*
*   A_(NM_,N_) + A(NM,NX-N_) -> A0(NM0,NX)
*
       integer                :: nx, nm0, nm, n_, nm_
       real(8)                :: a0(*)
       real(8)                :: a(nm,*), a_(nm_,*)

       include 'CSTAB.h'

       integer, parameter     ::  n_thread=8
       integer                ::  i_thread

       integer, parameter     :: LX = (L2_SIZE/8/4)
       real(8), pointer       :: buff(:)

       allocate(buff(LX*n_thread))
       call CSTAB_post_format_body(nx,a0,nm0,a,nm,n_,a_,nm_,buff)
       deallocate(buff)

       return
       end subroutine

       subroutine CSTAB_post_format_body(nx,a0,nm0,a,nm,n_,a_,nm_,buff)
       implicit none
*
*   A_(NM_,N_) + A(NM,NX-N_) -> A0(NM0,NX)
*
       integer                :: nx, nm0, nm, n_, nm_
       real(8)                :: a0(nm0,*)
       real(8)                :: a(nm,*), a_(nm_,*)

       include 'CSTAB.h'

       integer, parameter     :: n_thread=8
       integer                :: i_thread

       integer, parameter     :: LX = (L2_SIZE/8/4)
       real(8)                :: buff(LX,n_thread)

       integer                :: i_1,i_2,i_3
       integer                :: j_1,j_2,j_3
       integer                :: k_1,k_2,k_3
       integer                :: kk_1
       integer                :: MX, NN
       integer                :: i


       NN=MAX(1,MIN(nm0,nm))
       MX=(NN-1)/LX+1

       k_2=1
       k_3=nx*MX
       do kk_1=k_3,k_2,-n_thread

!$OMP PARALLEL DO PRIVATE(k_1,i_1,j_1,j_2,j_3)
          do  i_thread=1,MIN(n_thread,kk_1-k_2+1)
             k_1=kk_1-i_thread+1
             i_1=(k_1-1)/MX+1
             j_2=MOD(k_1-1,MX)*LX+1
             j_3=MIN(j_2+LX-1,NN)
*voption vec
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
             do j_1=j_2,j_3
                buff(j_1-j_2+1,i_thread)
     &          =a(j_1,i_1)
             enddo! j_1
          enddo! i_thread

!$OMP PARALLEL DO PRIVATE(k_1,i_1,j_1,j_2,j_3)
          do  i_thread=1,MIN(n_thread,kk_1-k_2+1)
             k_1=kk_1-i_thread+1
             i_1=(k_1-1)/MX+1
             j_2=MOD(k_1-1,MX)*LX+1
             j_3=MIN(j_2+LX-1,NN)
*voption vec
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
             do  j_1=j_2,j_3
                a0(j_1,i_1+n_)
     &          =buff(j_1-j_2+1,i_thread)
             end do! j_1
          end do! i_thread
       end do! i_1

       do  i_1=n_,1,-1
          j_2=1
          j_3=MIN(nm0,nm_)
*voption vec
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          do  j_1=j_2,j_3
             buff(j_1,1)=a_(j_1,i_1)
          enddo
*voption vec
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          do  j_1=j_2,j_3
             a0(j_1,i_1)=buff(j_1,1)
          end do! j_1
       end do! i_1


       return
       end subroutine

