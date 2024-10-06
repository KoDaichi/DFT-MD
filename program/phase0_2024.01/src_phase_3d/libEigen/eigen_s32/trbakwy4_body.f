       subroutine trbakwy_body2(n, z, nmz,
     $            d, e, v, vv, nm, m, i, ss, tt, nss, dcom, dx, dy)
       implicit double precision(a-h,o-z),integer(i-n)
*
       integer :: n, nmz, nm, m0, i, nss
       real(8) :: z(nmz,*),d(*),e(*),v(nm,*)
       real(8) :: ss(*), tt(*), dcom, dx, dy, ds, de
*
       include 'mpif.h'
!$     include 'omp_lib.h'
       include 'trd.h'
*
       real(8) :: vv(nsx,*)
*
*-
          i_2 = 1
          i_3 = n
*
          j_2 = loop_c_start(1,    x_nnod,x_inod)
          j_3 = loop_c_end  (i+m-2,x_nnod,x_inod)
*
!$OMP MASTER
! SM:= 0
! SS:= 0
          ss(1:(i_3-i_2+1)*m+nsm*nsm) = 0.0D+0
!$OMP END MASTER
*-
!$OMP BARRIER
! SM:= SM+lower(V^TV)
! SS(1:m,1:n):= SS+V(j_2:j_3,1:m)^T*Z(j_2:j_3,1:n)
          ds=MPI_Wtime()
          call trbk1(z, nmz,
     $               v, vv, nm, m, ss(1+nsm*nsm), ss(1),
     $               i_2,i_3,j_2,j_3)
!$OMP BARRIER
          de=MPI_Wtime()

!$OMP MASTER
          dx = dx + (de-ds)
*-
          ds=MPI_Wtime()
          call reduce_dbl(ss, tt, (i_3-i_2+1)*m+nsm*nsm,
     $         1, x_COMM_WORLD)
          de=MPI_Wtime()
          dcom=dcom+(de-ds)
*-
! SM:= D*SM
          do m_0=1,m-1
!DIR$ IVDEP
             do i_0=m_0+1,m
                ss(i_0+(m_0-1)*nsm) = ss(i_0+(m_0-1)*nsm) * d(i+i_0-1)
             end do! i_0
          end do! m_0
*-
! SS:= D*SS
          i_4 = MOD(i_3-i_2+1,4)+i_2
          do i_1=i_2,i_4-1                       ! 0
!DIR$ IVDEP
             do m_0=1,m
                ss((i_1+0-i_2)*m+m_0+nsm*nsm) =
     $          ss((i_1+0-i_2)*m+m_0+nsm*nsm) * d(i+m_0-1)
             end do! m_0
          end do! i_1
          do i_1=i_4,i_3,4                  ! 3
!DIR$ IVDEP
             do m_0=1,m
                ss((i_1+0-i_2)*m+m_0+nsm*nsm) =
     $          ss((i_1+0-i_2)*m+m_0+nsm*nsm) * d(i+m_0-1)
                ss((i_1+1-i_2)*m+m_0+nsm*nsm) =
     $          ss((i_1+1-i_2)*m+m_0+nsm*nsm) * d(i+m_0-1)
                ss((i_1+2-i_2)*m+m_0+nsm*nsm) =
     $          ss((i_1+2-i_2)*m+m_0+nsm*nsm) * d(i+m_0-1)
                ss((i_1+3-i_2)*m+m_0+nsm*nsm) =
     $          ss((i_1+3-i_2)*m+m_0+nsm*nsm) * d(i+m_0-1)
             end do! m_0
          end do! i_1
!$OMP END MASTER
*-
!$OMP BARRIER
! V:= V*(I-SM)^{-1}
! Z(j_2:j_3,1:n):= Z + V(j_2:j_3,1:m)*SS(1:m,1:n)
          ds=MPI_Wtime()
          call trbk2( z,
     $                nmz, v, vv, nm, m,
     $                ss(1+nsm*nsm), ss(1),
     $                i_2,i_3,j_2,j_3 )

!$OMP BARRIER
          de=MPI_Wtime()

!$OMP MASTER
          dy = dy + (de-ds)
!$OMP END MASTER
*-

       return
       end subroutine  trbakwy_body2
*-
       subroutine trbk1(z, nmz,
     $            v, vv, nm, m, ss, sm, i_2,i_3,j_2,j_3)
       implicit double precision(a-h,o-z),integer(i-n)
*
       integer :: nmz, nm, m
       real(8) :: z(nmz,*),v(nm,*)
       real(8) :: ss(m,*)
*
       include 'mpif.h'
!$     include 'omp_lib.h'
       include 'trd.h'
*
       real(8) :: vv(nsx,*), sm(nsm,nsm)
       integer :: local_rank, local_size
*
       local_size = 1
       local_rank = 0
!$     local_size = omp_get_num_threads()
!$     local_rank = omp_get_thread_num()
*
*-
          do j_0=j_2,j_3,nsx; j_4=MIN(j_0+nsx-1,j_3)

             do m_0=1+local_rank,m-1,local_size

                i_0 = m_0 + 1
                i_4 = MOD(m-m_0-1+1,4)+i_0
                if ( i_0 == i_4 - 1 ) then
                   s0 = sm(i_0+0,m_0)
!DIR$ IVDEP
!DIR$ VECTOR ALIGNED
                   do j_1=j_0,j_4
                      t0 = v(j_1,m_0)
                      v0 = v(j_1,i_0+0)
                      s0 = s0+v0*t0
                   end do! k
                   sm(i_0+0,m_0) = s0
                end if
                if ( i_0 == i_4 - 2 ) then
                   s0 = sm(i_0+0,m_0)
                   s1 = sm(i_0+1,m_0)
!DIR$ IVDEP
!DIR$ VECTOR ALIGNED
                   do j_1=j_0,j_4
                      t0 = v(j_1,m_0)
                      v0 = v(j_1,i_0+0)
                      v1 = v(j_1,i_0+1)
                      s0 = s0+v0*t0
                      s1 = s1+v1*t0
                   end do! k
                   sm(i_0+0,m_0) = s0
                   sm(i_0+1,m_0) = s1
                end if
                if ( i_0 == i_4 - 3 ) then
                   s0 = sm(i_0+0,m_0)
                   s1 = sm(i_0+1,m_0)
                   s2 = sm(i_0+2,m_0)
!DIR$ IVDEP
!DIR$ VECTOR ALIGNED
                   do j_1=j_0,j_4
                      t0 = v(j_1,m_0)
                      v0 = v(j_1,i_0+0)
                      v1 = v(j_1,i_0+1)
                      v2 = v(j_1,i_0+2)
                      s0 = s0+v0*t0
                      s1 = s1+v1*t0
                      s2 = s2+v2*t0
                   end do! k
                   sm(i_0+0,m_0) = s0
                   sm(i_0+1,m_0) = s1
                   sm(i_0+2,m_0) = s2
                end if
                do i_0=i_4,m,4              ! 3
                   s0 = sm(i_0+0,m_0)
                   s1 = sm(i_0+1,m_0)
                   s2 = sm(i_0+2,m_0)
                   s3 = sm(i_0+3,m_0)
!DIR$ IVDEP
!DIR$ VECTOR ALIGNED
                   do j_1=j_0,j_4
                      t0 = v(j_1,m_0)
                      v0 = v(j_1,i_0+0)
                      s0 = s0+v0*t0
                      v1 = v(j_1,i_0+1)
                      s1 = s1+v1*t0
                      v2 = v(j_1,i_0+2)
                      s2 = s2+v2*t0
                      v3 = v(j_1,i_0+3)
                      s3 = s3+v3*t0
                   end do! k
                   sm(i_0+0,m_0) = s0
                   sm(i_0+1,m_0) = s1
                   sm(i_0+2,m_0) = s2
                   sm(i_0+3,m_0) = s3
                end do! i_0

             end do! m_0

          end do! j_0
*-
          j_5 = j_3 - j_2 + 1

        do ii_2 = i_2, i_3, (1024*local_size)

          ii_3 = MIN(ii_2+(1024*local_size)-1,i_3)
          ii_4 = ii_2-i_2+1

          i_0 = (ii_3-ii_2) / local_size + 1
          i_4 = i_0 * local_rank
          i_5 = MIN(i_0, ii_3-(ii_2+i_4)+1)

          if ( j_5 > 0 .AND. i_5 > 0 ) then
             CALL DGEMM('T','N',
     &               m, i_5, j_5,
     &               1.0D+00, v (j_2    ,1       ), nm,
     &                        z (j_2    ,ii_2+i_4), nmz,
     &               1.0D+00, ss(1      ,ii_4+i_4), m)
          endif
        end do

*-
       return
       end subroutine  trbk1
*-
       subroutine trbk2( z,
     $            nmz,
     $            v, vv, nm, m, ss, sm, i_2,i_3, j_2,j_3 )
       implicit double precision(a-h,o-z),integer(i-n)
*
       integer :: nmz, nm, m
       real(8) :: z(nmz,*)
       real(8) :: v(nm,*)
       real(8) :: ss(m,*)
       real(8) :: dcom, ds, de
*
       include 'mpif.h'
!$     include 'omp_lib.h'
       include 'trd.h'
*
       real(8) :: vv(nsx,*), sm(nsm,*)
       integer :: local_rank, local_size

       local_size = 1
       local_rank = 0
!$     local_size = omp_get_num_threads()
!$     local_rank = omp_get_thread_num()
*
*-
          j_5 = (j_3-j_2) / local_size + 1
          j_5 = ((j_5-1)/2+1)*2
          j_4 = j_5 * local_rank
          j_5 = MIN(j_5, j_3-(j_2+j_4)+1)
          j_6 = j_2+j_4
          j_7 = j_6+j_5-1
         
          do j_0=j_6,j_7,nsx; j_8=MIN(j_0+nsx-1,j_7)

             do m_0=m,1,-2

                i_0 = m
                i_4 = MOD(m-m_0, 3)
                if ( i_4 == 1 ) then
!DIR$ IVDEP
!DIR$ VECTOR ALIGNED
                   do j_1=j_0,j_8
                      v0 = v(j_1,i_0-0)
                      u0 =
     $                    + sm(i_0-0,m_0  ) * v0
                      u1 =
     $                    + sm(i_0-0,m_0-1) * v0
                      v(j_1,m_0  ) = v(j_1,m_0  ) + u0
                      v(j_1,m_0-1) = v(j_1,m_0-1) + u1
                   end do! j_1
                end if
                if ( i_4 == 2 ) then
!DIR$ IVDEP
!DIR$ VECTOR ALIGNED
                   do j_1=j_0,j_8
                      v0 = v(j_1,i_0-0)
                      v1 = v(j_1,i_0-1)
                      u0 =
     $                    + sm(i_0-0,m_0  ) * v0
     $                    + sm(i_0-1,m_0  ) * v1
                      u1 =
     $                    + sm(i_0-0,m_0-1) * v0
     $                    + sm(i_0-1,m_0-1) * v1
                      v(j_1,m_0  ) = v(j_1,m_0  ) + u0
                      v(j_1,m_0-1) = v(j_1,m_0-1) + u1
                   end do! j_1
                end if

                do i_0=m-i_4,m_0+1,-3
!DIR$ IVDEP
!DIR$ VECTOR ALIGNED
                   do j_1=j_0,j_8
                      v0 = v(j_1,i_0-0)
                      v1 = v(j_1,i_0-1)
                      v2 = v(j_1,i_0-2)
                      u0 =
     $                    + sm(i_0-0,m_0  ) * v0
     $                    + sm(i_0-1,m_0  ) * v1
     $                    + sm(i_0-2,m_0  ) * v2
                      u1 =
     $                    + sm(i_0-0,m_0-1) * v0
     $                    + sm(i_0-1,m_0-1) * v1
     $                    + sm(i_0-2,m_0-1) * v2
                      v(j_1,m_0  ) = v(j_1,m_0  ) + u0
                      v(j_1,m_0-1) = v(j_1,m_0-1) + u1
                   end do! j_1
                end do! i_0

                i_0=m_0
!DIR$ IVDEP
!DIR$ VECTOR ALIGNED
                   do j_1=j_0,j_8
                      v0 = v(j_1,i_0-0)
                      u1 =
     $                    + sm(i_0-0,m_0-1) * v0
                      v(j_1,m_0-1) = v(j_1,m_0-1) + u1
                   end do! j_1

             end do! m_0

          end do! j_0

!$OMP BARRIER
*-

          j_5 = j_3 - j_2 + 1

          i_5 = (i_3-i_2) / local_size + 1
          i_4 = i_5 * local_rank
          i_5 = MIN(i_5, i_3-(i_2+i_4)+1)

          if ( j_5 > 0 .AND. i_5 > 0 ) then
             CALL DGEMM('N','N',
     &               j_5, i_5, m,
     &               1.0D+00, v (j_2    ,1      ), nm,
     &                        ss(1      ,1  +i_4), m,
     &               1.0D+00, z (j_2    ,i_2+i_4), nmz)
          end if
*-


       return
       end subroutine  trbk2

