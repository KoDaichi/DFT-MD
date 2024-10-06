       subroutine tred1(n, a, nma0,
     &                  d_out, e_out, nme0, m_orig)
       implicit NONE

       include 'mpif.h'

       integer                :: n, nma0, nme0, m_orig
       integer                :: nm, ne, m
       real(8)                :: a(*)
       real(8)                :: d_out(*), e_out(*)

       include 'trd.h'
       include 'CSTAB.h'

       integer                :: nm0, nx0
       integer                :: i_1,j_1, i, j, k, j_2, j_3
       integer                :: nx, ierr

       integer, parameter     :: nm_max_L1 = 16*4
       integer, parameter     :: nm_max_L2 = 16*4*2


       real(8)                :: d1, d2
       real(8)                :: c1, c2

       integer :: lda, ldz, nz

!Kro   nm = nma0
       ne = nme0
       m  = m_orig

       call datacast_init(2)

       lda = nma0
!Kro   ldz = (n-1)/x_nnod+1
!Kro   call CSTAB_get_optdim(ldz,9,nm_max_L1,nm_max_L2,nx)
!Kro   ldz = nx
!Kro   nz = (n-1)/y_nnod+1

       call tred1_(a, lda, d_out, e_out, ne, n, m)

       call datacast_free(1)

       return
       end subroutine ! tred1

       subroutine tred1_(a, nm, d_out, e_out, ne, n, m_orig)
       implicit NONE

       include 'mpif.h'

       integer                :: n, nm, ne, nmx, nv, nvx, m_orig
       real(8)                :: a(1:nm,*)
       real(8)                :: d_out(1:n)
       real(8)                :: e_out(1:ne,*)

       real(8), pointer       :: u_t(:), v_t(:), d_t(:)
       real(8), pointer       :: w(:)
       real(8), pointer       :: ux_vx(:)
       real(8), pointer       :: uy_vy(:)

       include 'trd.h'
       include 'CSTAB.h'
       include 'param.h'

       real(8)                :: d1,d2
       real(8)                :: d3,d4
       integer                :: nm_orig, nx, n1, n2, n3
       integer                :: i,j,k

       integer                :: ierr, kx

       integer                :: offset1, offset2, offset3
       integer                :: offset4, offset5, offset6
       integer                :: offset7, offset8

       integer, parameter     :: nm_max_L1 = 16*4
       integer, parameter     :: nm_max_L2 = 16*6


!Kro
       integer                :: ny, nv1, nv2
           nx = (n-1)/x_nnod+1 +2
!          call CSTAB_get_optdim(nx, 6, nm_max_L1, nm_max_L2, nv)
           call CSTAB_get_optdim(nx, 6, nm_max_L1, nm_max_L2, nv1)
!Kro
           ny = (n-1)/y_nnod+1 +2
           call CSTAB_get_optdim(ny, 6, nm_max_L1, nm_max_L2, nv2)
           nv=max(nv1,nv2)
           if(TRD_inod==1)then
             PRINT*,"NV=",nv
           endif

       allocate(
     &           w(1:nm*m_orig+n_columns),
     &           v_t(1:MAX(2*nv,2*(n+3),2*8*m_orig)+n_columns),
     &           u_t(1:MAX(2*nv,2*(n+3),2*8*m_orig)+n_columns),
     &           d_t(1:nv+n_columns),
     &           ux_vx(1:nv*2*m_orig+2*n_columns),
     &           uy_vy(1:nv*2*m_orig+2*n_columns),
     &           stat=ierr)
       if ( ierr /= 0 ) then
          if ( TRD_inod == 1 ) print*,"Memory allocation error."
          call MPI_Abort(TRD_COMM_WORLD, 1, ierr)
       end if

       w = ZERO
       v_t = ZERO
       u_t = ZERO
       d_t = ZERO
       ux_vx = ZERO
       uy_vy = ZERO

       kx = nv*m_orig+n_columns
       call CSTAB_adjust_base(ux_vx(1), a(1,1),offset1)
       call CSTAB_adjust_base(ux_vx(1+kx), a(1,1),offset3)
       call CSTAB_adjust_base(uy_vy(1), a(1,1),offset2)
       call CSTAB_adjust_base(uy_vy(1+kx), a(1,1),offset4)
       call CSTAB_adjust_base(u_t(1),a(1,1),offset5)
       call CSTAB_adjust_base(v_t(1),a(1,1),offset6)
       call CSTAB_adjust_base(w(1),a(1,1),offset7)
       kx = (L1_WINDOW/8)
!     &           +(L1_WINDOW)
!     &           +(L1_LSIZE/8)
     &           +(L1_LSIZE)
     &           +(L2_LSIZE/8)
       offset1 = offset1+kx*1
       offset2 = offset2+kx*2
       offset3 = offset3+kx*3
       offset4 = offset4+kx*4
       offset5 = offset5+kx*5
       offset6 = offset6+kx*6
       offset7 = offset7+kx*7
       call CSTAB_round_offset(offset1)
       call CSTAB_round_offset(offset2)
       call CSTAB_round_offset(offset3)
       call CSTAB_round_offset(offset4)
       call CSTAB_round_offset(offset5)
       call CSTAB_round_offset(offset6)
       call CSTAB_round_offset(offset7)

!!!
       kx = nv*m_orig+n_columns

!$OMP PARALLEL
       call trd_body(a, nm, d_out, e_out, ne,
     &               n, m_orig,
     &               w(1+offset7),
     &               ux_vx(1   +offset1),    ! u_x(1:nv,m)
     &               uy_vy(1   +offset2),    ! u_y(1:nv,m)
     &               ux_vx(1+kx+offset3),    ! v_x(1:nv,m)
     &               uy_vy(1+kx+offset4),    ! v_y(1:nv,m)
     &               u_t(1+offset5),
     &               v_t(1+offset6),
     &               d_t(1),
     &               nv)
!$OMP END PARALLEL

       deallocate(w)
       deallocate(v_t)
       deallocate(u_t)
       deallocate(d_t)
       deallocate(ux_vx)
       deallocate(uy_vy)

       return
       end subroutine ! tred1_

       subroutine trd_body(a, nm, d_out, e_out, ne,
     &               n, m_orig,
     &               w, u_x, u_y, v_x, v_y, u_t, v_t, d_t, nv)
       implicit NONE

       integer, intent(in)    :: n, nm, ne, nv, m_orig
       real(8)                :: a(1:nm,*)
       real(8)                :: d_out(1:n), e_out(1:ne,*)
       real(8)                :: w(1:nm,*)
       real(8)                :: u_x(1:nv,*), u_y(1:nv,*)
       real(8)                :: v_x(1:nv,*), v_y(1:nv,*)
       real(8)                :: u_t(1:2*nv), v_t(1:2*nv)
       real(8)                :: d_t(1:nv)

       include 'trd.h'

       real(8)                :: c(MBAND,MBAND), t1, t2, t3
       integer                :: i, j, k, l, m

       integer                :: i_1, i_2, i_3, i_4
       integer                :: j_1, j_2, j_3, j_4
       integer                :: k_1, k_2, k_3, k_4
       integer                :: l_1, l_2, l_3, l_4

       integer                :: m0, mm
       integer                :: i_block, i_base

       integer                :: ii

       include 'mpif.h'
!$     include 'omp_lib.h'

       real(8)                :: d1,d2,dd(100)

!
! initialization
!
       dd = 0.0D+00
       m  = m_orig

!$OMP MASTER

       if ( TRD_inod == 1 ) then
          print*,"NUM.OF.PROCESS=",TRD_nnod,"(",x_nnod,y_nnod,")"
!$        print*,"NUM.OF.THREADS=",omp_get_num_threads()
       endif
       call flush(6)

       call tred1_init(a(1,1), nm, n,
     &         d_out(1), e_out(1,1), ne,
     &         u_t(1), v_t(1), nv)

!$OMP END MASTER

       if ( m > n ) then
          m = n
       end if
!
! main iteration to reduce the input matrix into a penta-diagonal
!
       mm = ((n-(MBAND+MOD(n,MBAND)))-1)/m+1+1

!$OMP BARRIER

       do i_block=mm,2,-1

!$OMP BARRIER

          i_base = (i_block-2)*m+(MBAND+MOD(n,MBAND))
          m0     = MIN(m,n-i_base)

!$OMP MASTER

!          print*, i_block, i_base, m0

          call tred1_load(a(1,1), w(1,1), nm,
     &            d_t(1),
     &            u_x(1,1), u_y(1,1), v_x(1,1), v_y(1,1), nv,
     &            m0, i_base, i_block)

!$OMP END MASTER

          k_2 = m0
          k_3 = MAX(1,3*(2-i_block))

!$OMP BARRIER

          do k_1=k_2,k_3,-MBAND; k_4=k_1-MBAND+1

!$OMP BARRIER

             i = i_base+k_1
!
! u=...
!
!$OMP MASTER

             d1=MPI_Wtime()
             call tred1_u(
     &               w(1,k_4), nm,
     &               u_x(1,k_4), u_y(1,k_4), nv,
     &               u_t(1), v_t(1), i,
     &               c(1,1), e_out(1,1), ne)
             d2=MPI_Wtime()
             dd(1)=dd(1)+(d2-d1)

!$OMP END MASTER

!$OMP BARRIER

!$OMP MASTER
             d1=MPI_Wtime()
!$OMP END MASTER
!
! v:=Au
!
             call tred1_au(
     &               a(1,1), nm,
     &               u_x(1,k_4), u_y(1,k_4), v_x(1,k_4), nv,
     &               u_t(1), v_t(1), d_t(1),
     &               i, i_base, m0)

!$OMP BARRIER

!$OMP MASTER
             d2=MPI_Wtime()
             dd(2)=dd(2)+(d2-d1)
!
! v=v-(UV+VU)u
!
             d1=MPI_Wtime()
             call tred1_vorhto(
     &               u_x(1,1), v_x(1,1), nv,
     &               u_t(1), v_t(1),
     &               i, i_base, m0)
             d2=MPI_Wtime()
             dd(5)=dd(5)+(d2-d1)

!
! v':= v-((u,v)/2|u|^2)u
!
             d1=MPI_Wtime()
             call tred1_v(
     &               u_x(1,k_4), v_x(1,k_4), v_y(1,k_4), nv,
     &               u_t(1), v_t(1), c(1,1), i)
             d2=MPI_Wtime()
             dd(4)=dd(4)+(d2-d1)
!
! w':= w-uv-vu
!
             d1=MPI_Wtime()
             call tred1_local_2update(
     &               w(1,1), nm,
     &               u_x(1,1), u_y(1,1), v_x(1,1), v_y(1,1), nv,
     &               i_base, i)
             d2=MPI_Wtime()
             dd(6)=dd(6)+(d2-d1)

!$OMP END MASTER

!$OMP BARRIER

          end do! k_1

!$OMP BARRIER

!!$OMP MASTER
!             print*,"A="
!             print*,a(1:3,1)
!             print*,a(1:3,2)
!             print*,a(1:3,3)
!!$OMP END MASTER

!$OMP MASTER

          call tred1_store(a(1,1), w(1,1), nm,
     &            d_t(1),
     &            m0, i_base)
!$OMP END MASTER

!!$OMP MASTER
!             print*,"A'="
!             print*,a(1:3,1)
!             print*,a(1:3,2)
!             print*,a(1:3,3)
!!$OMP END MASTER


!$OMP BARRIER

!$OMP MASTER
          d1=MPI_Wtime()
!$OMP END MASTER
!
! A:=A-v^Tu-uT^v
!
!          if ( i_block >= 1 ) then
!!$OMP MASTER
!             print*,"A''="
!             print*,a(1:3,1)
!             print*,a(1:3,2)
!             print*,a(1:3,3)
!!$OMP END MASTER
             call tred1_2update(
     &               a(1,1), nm,
     &               u_x(1,1), u_y(1,1), v_x(1,1), v_y(1,1), nv,
     &               m0, i_base)
!!$OMP MASTER
!             print*,"A'''="
!             print*,a(1:3,1)
!             print*,a(1:3,2)
!             print*,a(1:3,3)
!!$OMP END MASTER
!          end if

!$OMP BARRIER

!$OMP MASTER
          d2=MPI_Wtime()
          dd(3)=dd(3)+(d2-d1)

!$OMP END MASTER

!$OMP BARRIER

       end do! i_block

!$OMP BARRIER

!$OMP MASTER

!       print*,"D=",d_out(1:n)
!       print*,"E=",e_out(1:n,1)
!       print*,"F=",e_out(1:n,2)


       call tred1_final(a(1,1), nm, n,
     &         d_out(1), e_out(1,1), ne, u_t(1))

!       print*,"D=",d_out(1:n)
!       print*,"E=",e_out(1:n,1)
!       print*,"F=",e_out(1:n,2)

       if(TRD_inod==1)then
            print*,"calc (u,beta)    ",dd(1)
            print*,"mat-vec (Au)     ",dd(2),dble(n)**3*2d-9/3/dd(2)
            print*,"2update (A-uv-vu)",dd(3),dble(n)**3*2d-9/3/dd(3)
            print*,"calc v           ",dd(4)
            print*,"v=v-(UV+VU)u     ",dd(5)
            print*,"UV post reduction",dd(6)
       endif
       call flush(6)

!$OMP END MASTER

       return
       end subroutine ! trd_body

