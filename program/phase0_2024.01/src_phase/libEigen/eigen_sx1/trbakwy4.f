       subroutine trbakwy(n, a, nma0, z, nmz0, e, m0)!, iblk)
       implicit double precision(a-h,o-z),integer(i-n)
       real(8) :: a(*)
       real(8) :: z(*)
       real(8) :: e(*)

       real(8) , pointer :: w(:),d(:)
       real(8) , pointer :: v(:)
       real(8) , pointer :: tt(:), ss(:)

       real(8) :: d1, d2

       include 'mpif.h'
!$     include 'omp_lib.h'
       include 'trd.h'
       include 'CSTAB.h'
*
       iblk = 2

       nma  = nma0
       nmz  = nmz0
       m    = MIN(nsm,m0)
       if ( m < 1 ) m = 1
*
       call datacast_init(2)
*
       na   = (n-1)/y_nnod+1
       na   = na  +MOD(na  -1,2)
       call CSTAB_get_optdim(nma,9,16*4,16*6,nm)
*
       allocate(
     &          d(1:n),
     &          v(1:MAX(nm*m,n)+n_columns),
     &          ss(1:na*m+nsm*nsm+n_columns),
     &          tt(1:na*m+nsm*nsm+n_columns),
     &          stat=i_stat)
       if(i_stat/=0)then
         if(TRD_inod==1)print*,"Memory allocation error."
         call mpi_abort(TRD_COMM_WORLD,1,ierr)
       endif
*
       call CSTAB_adjust_base(v(1), z(1), i_v)
       call CSTAB_adjust_base(ss(1), z(1), i_s)
       call CSTAB_adjust_base(tt(1), z(1), i_t)
          kx = (L1_WINDOW/4)
!     &           +(L1_WINDOW)
!     &           +(L1_LSIZE/8)
     &           +(L1_LSIZE)
     &           +(L2_LSIZE/4)
       i_v = i_v + kx*1
       i_s = i_s + kx*2
       i_t = i_t + kx*3
       call CSTAB_round_offset(i_v)
       call CSTAB_round_offset(i_s)
       call CSTAB_round_offset(i_t)
*
       call MPI_Barrier(TRD_COMM_WORLD,ierr)
       d1 = MPI_Wtime()
!$OMP PARALLEL
       call trbakwy_body( n,
     $      a(1), nma,
     $      z(1), nmz,
     $      e(1), d(1), v(1+i_v), nm, m, iblk,
     $      ss(1+i_s), tt(1+i_t), na )
!$OMP END PARALLEL
       call MPI_Barrier(TRD_COMM_WORLD,ierr)
       d2 = MPI_Wtime()
*
       deallocate( d )
       deallocate( v )
       deallocate( ss )
       deallocate( tt )
*
       call datacast_free(1)

       IF(TRD_inod==1)print*,"TRBAKWY",d2-d1
*
       return
       end
*
*
*
       subroutine trbakwy_body(n, a, nma, z, nmz, e, d,
     &             v, nm, m, iblk, ss, tt, nss)

       implicit double precision(a-h,o-z),integer(i-n)

       real(8) :: a(1:nma,*)
       real(8) :: z(1:nmz,*)
       real(8) :: e(1:n)
       real(8) :: d(1:n)

       real(8) :: v(1:nm,*)

       real(8) :: ss(*), tt(*)

       include 'mpif.h'
!$     include 'omp_lib.h'
       include 'trd.h'
*
       real(8) :: vv(nsx,nsm)
       integer :: nodes(0:nsm-1)

       real(8), pointer :: wk(:)

*
!$OMP MASTER
       call MPI_Barrier(TRD_COMM_WORLD,ierr)
       d1=MPI_Wtime()
       dcom=0.0D0+00
       dx = 0.0; dy = 0.0

       lwk=((m-1)/y_nnod+1)*((n-1)/x_nnod+1)
       allocate(wk(lwk))
!$OMP END MASTER
*
       i_2 = loop_c_start(1+iblk, y_nnod,y_inod)
       i_3 = loop_c_end  (n,      y_nnod,y_inod)
*
!$OMP MASTER
       d(1:n)=0.0
       do i_1=i_2,i_3
!          i = loop_c_l2g_depth(i_1, y_nnod,y_inod)
          i = (i_1-1)*y_nnod+y_inod
          l = i-iblk
!          x_root = loop_c_node(l, x_nnod,x_inod)
          x_root = MOD(l-1,x_nnod)+1
          if ( x_root == x_inod ) then
!             j_1 = loop_c_g2l_depth(l, x_nnod,x_inod)
             j_1  = (l-1)/x_nnod+1
             d(i) = a(j_1,i_1)
          end if
       end do! i_1
       call reduce_dbl( d, v, n, 1, TRD_COMM_WORLD )
       do i=1+iblk,n
          if ( e(i)*d(i) == 0.0D+00 ) then
             s0 = 0.0
          else
             s0 = (1.0D0/d(i))/e(i)
          end if
          d(i) = s0
       end do! i
!$OMP END MASTER
*
       i_2 = loop_c_start(1, y_nnod,y_inod)
       i_3 = loop_c_end  (n, y_nnod,y_inod)
*
       nx = MIN(MOD(n-(1+iblk)+1,m)+(1+iblk)-1,n)
*
       do i=(1+iblk),nx
*
!$OMP BARRIER
*
          if ( e(i) == 0.0 ) cycle
*
          j_2 = loop_c_start(1,   x_nnod,x_inod)
          j_3 = loop_c_end  (i-iblk, x_nnod,x_inod)
          i_4=MOD(i_3-i_2+1,4)+i_2
*
!$OMP MASTER
          ds=MPI_Wtime()
          nodes(0) = loop_c_node(i, y_nnod, y_inod)
          if ( nodes(0) == y_inod ) then
             i_1 = loop_c_g2l_depth(i, y_nnod, y_inod)
             do j_1=j_2,j_3
                v(j_1,1) = a(j_1,i_1)
             end do! j_1
          end if
          call bcast_dbl(v(j_2,1), j_3-j_2+1, nodes(0), y_COMM_WORLD)
          de=MPI_Wtime()
          dcom=dcom+(de-ds)
!$OMP END MASTER
*
!$OMP BARRIER
*
!$OMP MASTER
          if ( i_4 == i_2 + 1 ) then
             i_1 = i_2
             s0=0.0D+00
             do j_1=j_2,j_3
                s0=s0+v(j_1,1)*z(j_1,i_1+0)
             end do! j_1
             ss(i_1+0)=s0
          end if
          if ( i_4 == i_2 + 2 ) then
             i_1 = i_2
             s0=0.0D+00
             s1=0.0D+00
             do j_1=j_2,j_3
                s0=s0+v(j_1,1)*z(j_1,i_1+0)
                s1=s1+v(j_1,1)*z(j_1,i_1+1)
             end do! j_1
             ss(i_1+0)=s0
             ss(i_1+1)=s1
          end if
          if ( i_4 == i_2 + 3 ) then
             i_1 = i_2
             s0=0.0D+00
             s1=0.0D+00
             s2=0.0D+00
             do j_1=j_2,j_3
                s0=s0+v(j_1,1)*z(j_1,i_1+0)
                s1=s1+v(j_1,1)*z(j_1,i_1+1)
                s2=s2+v(j_1,1)*z(j_1,i_1+2)
             end do! j_1
             ss(i_1+0)=s0
             ss(i_1+1)=s1
             ss(i_1+2)=s2
          end if
!$OMP END MASTER
!$OMP DO
          do i_1=i_4,i_3,4
             s0=0.0D+00
             s1=0.0D+00
             s2=0.0D+00
             s3=0.0D+00
             do j_1=j_2,j_3
                s0=s0+v(j_1,1)*z(j_1,i_1+0)
                s1=s1+v(j_1,1)*z(j_1,i_1+1)
                s2=s2+v(j_1,1)*z(j_1,i_1+2)
                s3=s3+v(j_1,1)*z(j_1,i_1+3)
             end do! j_1
             ss(i_1+0)=s0
             ss(i_1+1)=s1
             ss(i_1+2)=s2
             ss(i_1+3)=s3
          end do! i_1
!$OMP END DO
*
!$OMP BARRIER
*
!$OMP MASTER
          ds=MPI_Wtime()
          call reduce_dbl(ss(i_2),tt, i_3-i_2+1, 1, x_COMM_WORLD)
          de=MPI_Wtime()
          dcom=dcom+(de-ds)

          s0 = d(i)
          do i_1=i_2,i_3
             ss(i_1) = ss(i_1) * s0
          end do! i_1
!$OMP END MASTER

!$OMP BARRIER

          j_2 = loop_c_start(1,      x_nnod,x_inod)
          j_3 = loop_c_end  (i-iblk, x_nnod,x_inod)
          i_4=MOD(i_3-i_2+1,4)+i_2

!$OMP MASTER
          if ( i_4 == i_2 + 1 ) then
             i_1 = i_2
             s0 = ss(i_1+0)
             do j_1=j_2,j_3
                z(j_1,i_1+0) = z(j_1,i_1+0) + s0 * v(j_1,1)
             end do! j_1
          end if
          if ( i_4 == i_2 + 2 ) then
             i_1 = i_2
             s0 = ss(i_1+0)
             s1 = ss(i_1+1)
             do j_1=j_2,j_3
                z(j_1,i_1+0) = z(j_1,i_1+0) + s0 * v(j_1,1)
                z(j_1,i_1+1) = z(j_1,i_1+1) + s1 * v(j_1,1)
             end do! j_1
          end if
          if ( i_4 == i_2 + 3 ) then
             i_1 = i_2
             s0 = ss(i_1+0)
             s1 = ss(i_1+1)
             s2 = ss(i_1+2)
             do j_1=j_2,j_3
                z(j_1,i_1+0) = z(j_1,i_1+0) + s0 * v(j_1,1)
                z(j_1,i_1+1) = z(j_1,i_1+1) + s1 * v(j_1,1)
                z(j_1,i_1+2) = z(j_1,i_1+2) + s2 * v(j_1,1)
             end do! j_1
          end if
!$OMP END MASTER
!$OMP DO
          do i_1=i_4,i_3,4
             s0 = ss(i_1+0)
             s1 = ss(i_1+1)
             s2 = ss(i_1+2)
             s3 = ss(i_1+3)
             do j_1=j_2,j_3
                z(j_1,i_1+0) = z(j_1,i_1+0) + s0 * v(j_1,1)
                z(j_1,i_1+1) = z(j_1,i_1+1) + s1 * v(j_1,1)
                z(j_1,i_1+2) = z(j_1,i_1+2) + s2 * v(j_1,1)
                z(j_1,i_1+3) = z(j_1,i_1+3) + s3 * v(j_1,1)
             end do! j_1
          end do! i_1
!$OMP ENDDO
*
!$OMP BARRIER
*
       end do
*
       do i=nx+1, n, m

*
!$OMP BARRIER
*
!$OMP MASTER
          ds=MPI_Wtime()

        IF ( m > y_nnod ) THEN

          do j=0,m-1
             nodes(j) = loop_c_node(i+j, y_nnod, y_inod)
          enddo

          do iy=1,y_nnod

             k0=0
             do j=0,m-1
             if ( nodes(j) == iy ) then
                i_1 = loop_c_g2l_depth(i+j, y_nnod, iy)
                j_2 = loop_c_start(1,    x_nnod,x_inod)
                j_3 = loop_c_end  (i+m-1-iblk,x_nnod,x_inod)
                if ( y_inod == iy ) then
                   do j_1=j_2,j_3
                      wk(k0+j_1) = a(j_1,i_1)
                   end do! k
                endif
                k0=k0+(j_3-j_2+1)
             endif
             enddo

             call bcast_dbl(wk, k0, iy, y_COMM_WORLD)

             k0=0
             do j=0,m-1
             if ( nodes(j) == iy ) then
                i_1 = loop_c_g2l_depth(i+j, y_nnod, iy)
                j_2 = loop_c_start(1,    x_nnod,x_inod)
                j_3 = loop_c_end  (i+m-1-iblk,x_nnod,x_inod)
                do j_1=j_2,j_3
                   v(j_1,j+1) = wk(k0+j_1)
                end do! k
                k0=k0+(j_3-j_2+1)
                j_2 = loop_c_start(i+j,  x_nnod,x_inod)
                j_3 = loop_c_end  (i+m-1-iblk,x_nnod,x_inod)
                do j_1=j_2,j_3
                   v(j_1,j+1) = 0.0D+00
                end do
             endif
             enddo

          enddo

        ELSE

          do j=0,m-1
             nodes(j) = loop_c_node(i+j, y_nnod, y_inod)
             if ( nodes(j) == y_inod ) then
                i_1 = loop_c_g2l_depth(i+j, y_nnod, y_inod)
                j_2 = loop_c_start(1,    x_nnod,x_inod)
                j_3 = loop_c_end  (i+m-1-iblk,x_nnod,x_inod)
                do j_1=j_2,j_3
                   v(j_1,j+1) = a(j_1,i_1)
                end do! k
                j_2 = loop_c_start(i+j,  x_nnod,x_inod)
                j_3 = loop_c_end  (i+m-1-iblk,x_nnod,x_inod)
!CDIR NOVECTOR
                do j_1=j_2,j_3
                   v(j_1,j+1) = 0.0D+00
                end do
             end if
          end do
          j_2 = loop_c_start(1,     x_nnod,x_inod)
          j_3 = loop_c_end  (i+m-1-iblk, x_nnod,x_inod)
          do j=0,m-1
             call bcast_dbl(v(1,j+1), j_3-j_2+1,
     &                      nodes(j), y_COMM_WORLD)
          enddo

        ENDIF

          de=MPI_Wtime()
          dcom=dcom+(de-ds)
!$OMP END MASTER
*
!$OMP BARRIER
*
          call trbakwy_body2(i_3, z, nmz,
     $         d, e, v, vv, nm, m, iblk, i, ss, tt, nss, dcom, dx, dy)
*
!$OMP BARRIER
*
       end do
*
!$OMP MASTER
       deallocate(wk)
       call MPI_Barrier(TRD_COMM_WORLD,ierr)
       d2=MPI_Wtime()
       if ( 1 == TRD_inod ) then
          print*,"TRBAK=",(d2-d1)
          print*,"COMM=",dcom
          print*,"   ",(2d0*n*n*n)/(d2-d1)*1d-9,"GFLOPS"
          print*,"   ",(1d0*n*n*n)/(dx)*1d-9,"GFLOPS"
          print*,"   ",(1d0*n*n*n)/(dy)*1d-9,"GFLOPS"
       end if
!$OMP END MASTER
*
*
       return
       end subroutine

