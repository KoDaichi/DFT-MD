      subroutine eigen_hbk(ar,ai,lda, zr,zi,ldz, tau,ldtau, n)

      use communication_h , only : eigen_init
     &                , eigen_free
     &                , cyc1d_cyc2d
     &                , cyc2d_cyc1d

      implicit double precision(a-h,o-z),integer(i-n)

      real(8) :: ar(*), ai(*)
      real(8) :: zr(*), zi(*)
      real(8) :: tau(*)

      real(8) , allocatable :: wkr(:), wki(:)
      real(8) , pointer :: vr(:),vi(:), wr(:),wi(:)
      real(8) , pointer :: ssr(:),ssi(:)
      real(8) , pointer :: tt(:)

      include 'mpif.h'
!$    include 'omp_lib.h'
      include 'trd.h'

      call eigen_init(2)
      m = 128; mm = m*(m-1)/2/1

      s0=mpi_wtime()
      lda2 = (n-1)/size_of_col+1
      na   = (n-1)/size_of_row+1
      lda2 = lda2+mod(lda2-1,2)
      na   = na  +mod(na  -1,2)
      allocate( wkr(lda2*na) )
      allocate( wki(lda2*na) )

      call cyc1d_cyc2d(n,zr(1),ldz,wkr(1),lda2,size_of_col,size_of_row,
     &                 ierr)
      call cyc1d_cyc2d(n,zi(1),ldz,wki(1),lda2,size_of_col,size_of_row,
     &                 ierr)
      call cyc1d_cyc2d(n,ar(1),lda,zr(1),lda2,size_of_col,size_of_row,
     &                 ierr)
      call cyc1d_cyc2d(n,ai(1),lda,zi(1),lda2,size_of_col,size_of_row,
     &                 ierr)

      nm = lda2
      allocate(
     &         vr(1:nm*m), vi(1:nm*m),
     &         wr(1:nm*m), wi(1:nm*m),
     &         ssr(1:max(na*m+mm,n)),
     &         ssi(1:max(na*m+mm,n)),
     &         tt(1:max(na*m+mm,n)),
     &         stat=i_stat)
      if(i_stat/=0)then
         if(myrank==1)print*,"memory allocation error."
         call mpi_abort(mpi_comm_eigen,1,ierr)
      endif

      call mpi_barrier(mpi_comm_eigen,ierr)
      call eigen_hbk_main1( 
     &       n, zr(1),zi(1), lda2,
     &       wkr(1),wki(1), lda2,
     &       tau(1),ldtau, m, mm,
     &       vr(1),vi(1), wr(1),wi(1), nm, ssr(1),ssi(1),tt(1)
     &     )
      call mpi_barrier(mpi_comm_eigen,ierr)

      call cyc2d_cyc1d(n,zr(1),ldz,wkr(1),lda2,size_of_col,size_of_row,
     &                 ierr)
      call cyc2d_cyc1d(n,zi(1),ldz,wki(1),lda2,size_of_col,size_of_row,
     &                 ierr)

      deallocate( vr,vi )
      deallocate( wr,wi )
      deallocate( ssr,ssi )
      deallocate( tt )
      deallocate( wkr,wki )

      call eigen_free(3)

#ifdef TIMER
      s1=mpi_wtime()
      if(myrank==1)then
          print*,"Exectime of \"eigen_hbk\" routine =",s1-s0,"(sec)"
!----     s1-s0,8e-9*dble(n)**3/(s1-s0),"GFLOPS"
      endif
#endif

      return
      end subroutine
!----
       subroutine eigen_hbk_main1(n,ar,ai,nma,zr,zi,nmz,tau,nmtau,
     &                          m,mm, vr,vi,wr,wi,nm,ssr,ssi,tt)
!----
      use communication_h, only : get_loop_start
     &               , get_loop_end
     &               , get_owner_node
     &               , translate_l2g
     &               , translate_g2l
     &               , bcast_dbl
     &               , reduce_dbl
!----
       implicit double precision(a-h,o-z),integer(i-n)

       real(8) :: ar(1:nma,*), ai(1:nma,*)
       real(8) :: zr(1:nmz,*), zi(1:nmz,*)
       real(8) :: tau(1:nmtau,2)

       real(8) :: vr(1:nm,*), vi(1:nm,*)
       real(8) :: wr(1:nm,*), wi(1:nm,*)

       real(8) :: ssr(*),ssi(*), tt(*)

       integer :: nodes(0:m-1)

       include 'mpif.h'
!$     include 'omp_lib.h'
       include 'trd.h'
!----
       call mpi_barrier(mpi_comm_eigen,ierr)
       d1=mpi_wtime()
       dcom=0.0d0+00
!----
       i_2 = get_loop_start(1, size_of_row,my_row)
       i_3 = get_loop_end  (n, size_of_row,my_row)
       j_2 = get_loop_start(1, size_of_col,my_col)
       j_3 = get_loop_end  (n, size_of_col,my_col)

       do i_1=i_2,i_3
          do j_1=j_2,j_3
             k = (j_1-1)*size_of_col+my_col
             zi(j_1,i_1) =  zr(j_1,i_1) * tau(k,2)
             zr(j_1,i_1) =  zr(j_1,i_1) * tau(k,1)
          end do! j_1
       end do! i_1
       tau(1:n,1) = 0.0
       do i_1=i_2,i_3
          k = (i_1-1)*size_of_row+my_row
          x_root = mod(k-1,size_of_col)+1
          if ( x_root == my_col ) then
             j_1  = (k-1)/size_of_col+1
             tau(k,1) = ai(j_1,i_1)
          end if
       end do! i_1
       call reduce_dbl(tau(1,1), tau(1,2), n, mpi_comm_eigen)

       i=1
1000   i=i+1; if ( i > n ) go to 9000

       if ( tau(i,1) == 0.0d+00 ) go to 1000

1100   if ( (i+1<=n) .and. (tau(i+1,1)/=0.0d+00) ) go to 1200

1111   continue

          j_2 = get_loop_start(1,   size_of_col,my_col)
          j_3 = get_loop_end  (i-1, size_of_col,my_col)
!----
          ds=mpi_wtime()
          nodes(0) = get_owner_node(i, size_of_row, my_row)
          if ( nodes(0) == my_row ) then
             i_1 = translate_g2l(i, size_of_row, my_row)
             do j_1=j_2,j_3
                vr(j_1,1) = ar(j_1,i_1)
                vi(j_1,1) = ai(j_1,i_1)
             end do! j_1
          end if
          call bcast_dbl(vr(j_2,1), j_3-j_2+1, nodes(0), mpi_comm_row)
          call bcast_dbl(vi(j_2,1), j_3-j_2+1, nodes(0), mpi_comm_row)
          de=mpi_wtime()
          dcom=dcom+(de-ds)
!----
          i_4=mod(i_3-i_2+1,8)+i_2
          do i_1=i_2,i_4-1 !                  0
             sr0=0.0d+00
             si0=0.0d+00
             do j_1=j_2,j_3
                sr0=sr0+vr(j_1,1)*zr(j_1,i_1+0)
     &                       +vi(j_1,1)*zi(j_1,i_1+0)
                si0=si0+vr(j_1,1)*zi(j_1,i_1+0)
     &                       -vi(j_1,1)*zr(j_1,i_1+0)
             end do! j_1
             ssr(i_1+0)=sr0
             ssi(i_1+0)=si0
          end do! i_1
          do i_1=i_4,i_3,8 !             7
             sr0=0.0d+00
             si0=0.0d+00
             sr1=0.0d+00
             si1=0.0d+00
             sr2=0.0d+00
             si2=0.0d+00
             sr3=0.0d+00
             si3=0.0d+00
             sr4=0.0d+00
             si4=0.0d+00
             sr5=0.0d+00
             si5=0.0d+00
             sr6=0.0d+00
             si6=0.0d+00
             sr7=0.0d+00
             si7=0.0d+00
             do j_1=j_2,j_3
                sr0=sr0+vr(j_1,1)*zr(j_1,i_1+0)
     &                       +vi(j_1,1)*zi(j_1,i_1+0)
                si0=si0+vr(j_1,1)*zi(j_1,i_1+0)
     &                       -vi(j_1,1)*zr(j_1,i_1+0)
                sr1=sr1+vr(j_1,1)*zr(j_1,i_1+1)
     &                       +vi(j_1,1)*zi(j_1,i_1+1)
                si1=si1+vr(j_1,1)*zi(j_1,i_1+1)
     &                       -vi(j_1,1)*zr(j_1,i_1+1)
                sr2=sr2+vr(j_1,1)*zr(j_1,i_1+2)
     &                       +vi(j_1,1)*zi(j_1,i_1+2)
                si2=si2+vr(j_1,1)*zi(j_1,i_1+2)
     &                       -vi(j_1,1)*zr(j_1,i_1+2)
                sr3=sr3+vr(j_1,1)*zr(j_1,i_1+3)
     &                       +vi(j_1,1)*zi(j_1,i_1+3)
                si3=si3+vr(j_1,1)*zi(j_1,i_1+3)
     &                       -vi(j_1,1)*zr(j_1,i_1+3)
                sr4=sr4+vr(j_1,1)*zr(j_1,i_1+4)
     &                       +vi(j_1,1)*zi(j_1,i_1+4)
                si4=si4+vr(j_1,1)*zi(j_1,i_1+4)
     &                       -vi(j_1,1)*zr(j_1,i_1+4)
                sr5=sr5+vr(j_1,1)*zr(j_1,i_1+5)
     &                       +vi(j_1,1)*zi(j_1,i_1+5)
                si5=si5+vr(j_1,1)*zi(j_1,i_1+5)
     &                       -vi(j_1,1)*zr(j_1,i_1+5)
                sr6=sr6+vr(j_1,1)*zr(j_1,i_1+6)
     &                       +vi(j_1,1)*zi(j_1,i_1+6)
                si6=si6+vr(j_1,1)*zi(j_1,i_1+6)
     &                       -vi(j_1,1)*zr(j_1,i_1+6)
                sr7=sr7+vr(j_1,1)*zr(j_1,i_1+7)
     &                       +vi(j_1,1)*zi(j_1,i_1+7)
                si7=si7+vr(j_1,1)*zi(j_1,i_1+7)
     &                       -vi(j_1,1)*zr(j_1,i_1+7)
             end do! j_1
             ssr(i_1+0)=sr0
             ssi(i_1+0)=si0
             ssr(i_1+1)=sr1
             ssi(i_1+1)=si1
             ssr(i_1+2)=sr2
             ssi(i_1+2)=si2
             ssr(i_1+3)=sr3
             ssi(i_1+3)=si3
             ssr(i_1+4)=sr4
             ssi(i_1+4)=si4
             ssr(i_1+5)=sr5
             ssi(i_1+5)=si5
             ssr(i_1+6)=sr6
             ssi(i_1+6)=si6
             ssr(i_1+7)=sr7
             ssi(i_1+7)=si7
          end do! i_1

          ds=mpi_wtime()
          call reduce_dbl(ssr(i_2), tt, i_3-i_2+1, mpi_comm_col)
          call reduce_dbl(ssi(i_2), tt, i_3-i_2+1, mpi_comm_col)
          de=mpi_wtime()
          dcom=dcom+(de-ds)

          do i_1=i_2,i_3
             h0 = tau(i,1)
             if ( h0 /= 0.0d+00 ) then
                ssr(i_1) = (ssr(i_1)/h0)/h0
                ssi(i_1) = (ssi(i_1)/h0)/h0
             else
                ssr(i_1) = 0.0d+00
                ssi(i_1) = 0.0d+00
             endif
          end do! i_1

          do i_1=i_2,i_4-1 !                  0
             sr0=ssr(i_1+0)
             si0=ssi(i_1+0)
             do j_1=j_2,j_3
                zr(j_1,i_1+0)=zr(j_1,i_1+0)
     &                                           -sr0*vr(j_1,1)
     &                                           +si0*vi(j_1,1)
                zi(j_1,i_1+0)=zi(j_1,i_1+0)
     &                                           -si0*vr(j_1,1)
     &                                           -sr0*vi(j_1,1)
             end do! j_1
          end do! i_1
          do i_1=i_4,i_3,8 !             7
             sr0=ssr(i_1+0)
             si0=ssi(i_1+0)
             sr1=ssr(i_1+1)
             si1=ssi(i_1+1)
             sr2=ssr(i_1+2)
             si2=ssi(i_1+2)
             sr3=ssr(i_1+3)
             si3=ssi(i_1+3)
             sr4=ssr(i_1+4)
             si4=ssi(i_1+4)
             sr5=ssr(i_1+5)
             si5=ssi(i_1+5)
             sr6=ssr(i_1+6)
             si6=ssi(i_1+6)
             sr7=ssr(i_1+7)
             si7=ssi(i_1+7)
             do j_1=j_2,j_3
                zr(j_1,i_1+0)=zr(j_1,i_1+0)
     &                                           -sr0*vr(j_1,1)
     &                                           +si0*vi(j_1,1)
                zi(j_1,i_1+0)=zi(j_1,i_1+0)
     &                                           -si0*vr(j_1,1)
     &                                           -sr0*vi(j_1,1)
                zr(j_1,i_1+1)=zr(j_1,i_1+1)
     &                                           -sr1*vr(j_1,1)
     &                                           +si1*vi(j_1,1)
                zi(j_1,i_1+1)=zi(j_1,i_1+1)
     &                                           -si1*vr(j_1,1)
     &                                           -sr1*vi(j_1,1)
                zr(j_1,i_1+2)=zr(j_1,i_1+2)
     &                                           -sr2*vr(j_1,1)
     &                                           +si2*vi(j_1,1)
                zi(j_1,i_1+2)=zi(j_1,i_1+2)
     &                                           -si2*vr(j_1,1)
     &                                           -sr2*vi(j_1,1)
                zr(j_1,i_1+3)=zr(j_1,i_1+3)
     &                                           -sr3*vr(j_1,1)
     &                                           +si3*vi(j_1,1)
                zi(j_1,i_1+3)=zi(j_1,i_1+3)
     &                                           -si3*vr(j_1,1)
     &                                           -sr3*vi(j_1,1)
                zr(j_1,i_1+4)=zr(j_1,i_1+4)
     &                                           -sr4*vr(j_1,1)
     &                                           +si4*vi(j_1,1)
                zi(j_1,i_1+4)=zi(j_1,i_1+4)
     &                                           -si4*vr(j_1,1)
     &                                           -sr4*vi(j_1,1)
                zr(j_1,i_1+5)=zr(j_1,i_1+5)
     &                                           -sr5*vr(j_1,1)
     &                                           +si5*vi(j_1,1)
                zi(j_1,i_1+5)=zi(j_1,i_1+5)
     &                                           -si5*vr(j_1,1)
     &                                           -sr5*vi(j_1,1)
                zr(j_1,i_1+6)=zr(j_1,i_1+6)
     &                                           -sr6*vr(j_1,1)
     &                                           +si6*vi(j_1,1)
                zi(j_1,i_1+6)=zi(j_1,i_1+6)
     &                                           -si6*vr(j_1,1)
     &                                           -sr6*vi(j_1,1)
                zr(j_1,i_1+7)=zr(j_1,i_1+7)
     &                                           -sr7*vr(j_1,1)
     &                                           +si7*vi(j_1,1)
                zi(j_1,i_1+7)=zi(j_1,i_1+7)
     &                                           -si7*vr(j_1,1)
     &                                           -sr7*vi(j_1,1)
             end do! j_1
          end do! i_1
          i=i+0; goto 1000

1200      continue

            i_ng=0
            do i_1=2,m-1
               if ( (i+i_1<=n) .and. (tau(i+i_1,1)/=0.0) ) then
                  continue
               else
                  i_ng=i_ng+1
                  exit
               endif
            enddo
            if (i_ng/=0) goto 1111

          ds=mpi_wtime()
          do j=0,m-1
             nodes(j) = get_owner_node(i+j, size_of_row, my_row)
             if ( nodes(j) == my_row ) then
                i_1 = translate_g2l(i+j, size_of_row, my_row)
                j_2 = get_loop_start(1,  size_of_col,my_col)
                j_3 = get_loop_end  (i+m-2,size_of_col,my_col)
                do j_1=j_2,j_3
                   vr(j_1,j+1)=ar(j_1,i_1)
                   vi(j_1,j+1)=ai(j_1,i_1)
                end do! k
                j_2 = get_loop_start(i+j,size_of_col,my_col)
                j_3 = get_loop_end  (i+m-2,size_of_col,my_col)
!cdir novector
                do j_1=j_2,j_3
                   vr(j_1,j+1)=0.0d+00
                   vi(j_1,j+1)=0.0d+00
                end do
             end if
          end do
          j_2 = get_loop_start(1,     size_of_col,my_col)
          j_3 = get_loop_end  (i+m-2, size_of_col,my_col)
          do j=0,m-1
             call bcast_dbl(vr(1,j+1), j_3-j_2+1,
     &                      nodes(j), mpi_comm_row)
             call bcast_dbl(vi(1,j+1), j_3-j_2+1,
     &                      nodes(j), mpi_comm_row)
          enddo
          de=mpi_wtime()
          dcom=dcom+(de-ds)

          call eigen_hbk_main2(
     &           i_3, zr, zi, nmz,
     &           tau, nmtau, m, mm,
     &           vr,vi, wr,wi, nm, i, ssr,ssi,tt, dcom
     &         )

          i=i+m-1; goto 1000

       continue

9000   continue

       call mpi_barrier(mpi_comm_eigen,ierr)
#ifdef DETAIL
       d2=mpi_wtime()
       if (myrank==1) then
          print*," "
          print*,"detail of exectime in eigen_hbk "
          print*,"   time of eigen_hbk_main1=",(d2-d1),"(sec)"
          print*,"   communication in eigen_hbk_main1=",dcom,"(sec)"
!----     print*,"   ",(8d0*n*n*n)/(d2-d1)*1d-9,"gflops"
       end if
#endif

       return
       end subroutine ! eigen_hbk_main1

       subroutine eigen_hbk_main2(n, zr, zi, nmz,
     &            tau, nmtau, m, mm,
     &            vr,vi, wr,wi, nm, i, ssr,ssi, tt, dcom)
!----
      use communication_h, only : get_loop_start
     &               , get_loop_end
     &               , reduce_dbl
!----
       implicit double precision(a-h,o-z),integer(i-n)

       integer :: n, nmz, nm, i
       real(8) :: zr(nmz,*),zi(nmz,*),tau(nmtau,2)
       real(8) :: vr(nm,*),vi(nm,*),wr(nm,*),wi(nm,*)
       real(8) :: ssr(*),ssi(*), tt(*), dcom, ds, de

       include 'mpif.h'
!$     include 'omp_lib.h'
       include 'trd.h'

          i_2=1
          i_3=n
          j_2 = get_loop_start(1,    size_of_col,my_col)
          j_3 = get_loop_end  (i+m-2,size_of_col,my_col)

          do iss=1,m*(m-1)/2
             ssr(iss)=0.0d0
             ssi(iss)=0.0d0
          enddo
          do j_10=j_2,j_3,128; j_4=j_10; j_5=min(j_3,j_10+128-1)
          iss=m*(m-1)/2
          do l_1=m,2,-1
          do k_1=l_1-1,1,-1
             sr0=ssr(iss)
             si0=ssi(iss)
             do j_1=j_4,j_5
             vr1=vr(j_1,l_1)
             vi1=vi(j_1,l_1)
             vr0=vr(j_1,k_1)
             vi0=vi(j_1,k_1)
             sr0=sr0+vr1*vr0+vi1*vi0
             si0=si0-vi1*vr0+vr1*vi0
             enddo
             ssr(iss)=sr0
             ssi(iss)=si0
             iss=iss-1
          enddo
          enddo
          enddo

          i_5 = i_3 - i_2 + 1
          j_5 = j_3 - j_2 + 1
          if ( j_5 > 0 .and. i_5 > 0 ) then
                 call dgemm('t','n',
     &               m, i_5, j_5,
     &               1.0d+00, vr(j_2    ,1       ), nm,
     &                        zr(j_2    ,i_2     ), nmz,
     &               0.0d+00, ssr(1     +mm     ), m)
                 call dgemm('t','n',
     &               m, i_5, j_5,
     &               1.0d+00, vi(j_2    ,1       ), nm,
     &                        zi(j_2    ,i_2     ), nmz,
     &               1.0d+00, ssr(1     +mm     ), m)

                 call dgemm('t','n',
     &               m, i_5, j_5,
     &              -1.0d+00, vi(j_2    ,1       ), nm,
     &                        zr(j_2    ,i_2     ), nmz,
     &               0.0d+00, ssi(1     +mm     ), m)
                 call dgemm('t','n',
     &               m, i_5, j_5,
     &               1.0d+00, vr(j_2    ,1       ), nm,
     &                        zi(j_2    ,i_2     ), nmz,
     &               1.0d+00, ssi(1     +mm     ), m)
          endif
!----
          ds=mpi_wtime()
          call reduce_dbl(ssr, tt, (i_3-i_2+1)*m+mm, mpi_comm_col)
          call reduce_dbl(ssi, tt, (i_3-i_2+1)*m+mm, mpi_comm_col)
          de=mpi_wtime()
          dcom=dcom+(de-ds)
!----
             iss=m*(m-1)/2
          do k_1=m,1,-1
             h0=tau(i+k_1-1,1)
          do i_1=k_1-1,1,-1
             ssr(iss)=(ssr(iss)/h0)/h0
             ssi(iss)=(ssi(iss)/h0)/h0
             iss=iss-1
          enddo
          do i_1=i_2,i_3
             ssr((i_1-i_2)*m+k_1+mm)=
     $          (ssr((i_1-i_2)*m+k_1+mm)/h0)/h0
             ssi((i_1-i_2)*m+k_1+mm)=
     $          (ssi((i_1-i_2)*m+k_1+mm)/h0)/h0
          end do! i_1
          end do! k_1
!----
          do j_10=j_2,j_3,128; j_4=j_10; j_5=min(j_3,j_10+128-1)
          do l_1=m,1,-1
             do j_1=j_4,j_5
                wr(j_1,l_1)=vr(j_1,l_1)
                wi(j_1,l_1)=vi(j_1,l_1)
             enddo
          enddo
             iss=m*(m-1)/2
          do l_1=m,2,-1
          do k_1=l_1-1,1,-1
             do j_1=j_4,j_5
                wr(j_1,k_1)=wr(j_1,k_1)
     &                     -wr(j_1,l_1)*ssr(iss)
     &                     +wi(j_1,l_1)*ssi(iss)
                wi(j_1,k_1)=wi(j_1,k_1)
     &                     -wi(j_1,l_1)*ssr(iss)
     &                     -wr(j_1,l_1)*ssi(iss)
             enddo
             iss=iss-1
          enddo
          enddo
          enddo
!----
          i_5=i_3-i_2+1
          j_5=j_3-j_2+1
          if ( j_5 > 0 .and. i_5 > 0 ) then
             call dgemm('n','n',
     &               j_5, i_5, m,
     &              -1.0d+00, wr(j_2 ,1  ), nm,
     &                        ssr(1   +mm), m,
     &               1.0d+00, zr(j_2 ,i_2), nmz)
             call dgemm('n','n',
     &               j_5, i_5, m,
     &              +1.0d+00, wi(j_2 ,1  ), nm,
     &                        ssi(1   +mm), m,
     &               1.0d+00, zr(j_2 ,i_2), nmz)
             call dgemm('n','n',
     &               j_5, i_5, m,
     &              -1.0d+00, wr(j_2 ,1  ), nm,
     &                        ssi(1   +mm), m,
     &               1.0d+00, zi(j_2 ,i_2), nmz)
             call dgemm('n','n',
     &               j_5, i_5, m,
     &              -1.0d+00, wi(j_2 ,1  ), nm,
     &                        ssr(1   +mm), m,
     &               1.0d+00, zi(j_2 ,i_2), nmz)
          endif

       return
       end subroutine ! eigen_hbk_main2
