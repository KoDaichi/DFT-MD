      subroutine eigen_hrd_2update(
     &             ar, ai, nm,
     &             ur_x, ui_x,
     &             ur_y, ui_y,
     &             vr_x, vi_x,
     &             vr_y, vi_y,
     &             nv, m, n,
     &             ar_,ai_,n_,nm_
     &           )

!======================================================================!
!$    use omp_lib
      implicit none
!----
      real(8), intent(inout) :: ar (nm,*),  ai (nm,*)
      real(8), intent(in)    :: ur_x(nv,*), ui_x(nv,*)
      real(8), intent(in)    :: ur_y(nv,*), vi_y(nv,*)
      real(8), intent(in)    :: vr_x(nv,*), vi_x(nv,*)
      real(8), intent(in)    :: vr_y(nv,*), ui_y(nv,*)
      integer, intent(in)    :: nm, nv, m, n
      real(8) :: ar_, ai_
      integer :: n_, nm_
!----
      if ( n <= 0 ) return
!----
      call eigen_hrd_2update_body(
     &       ar(1,1), ai(1,1),
     &       nm,
     &       ur_x(1,1), ui_x(1,1), ur_y(1,1), ui_y(1,1),
     &       vr_x(1,1), vi_x(1,1), vr_y(1,1), vi_y(1,1),
     &       nv, m, n
     &     )
!----
      end subroutine
!======================================================================!
      subroutine eigen_hrd_2update_body(
     &             ar_alias_0, ai_alias_0,
     &             nma,
     &             ur, ui, uyr, uyi, vr, vi, vyr, vyi, nmv,
     &             m_size, i_base
     &           )
!----
      use communication_h, only : get_loop_start
     &               , get_loop_end
     &               , translate_l2g
!----
!$    use omp_lib
      implicit none
!----
      real(8), intent(inout) :: ar_alias_0(nma,*), ai_alias_0(nma,*)
      real(8), intent(in)    :: ur (nmv,*), ui (nmv,*)
      real(8), intent(in)    :: uyr(nmv,*), uyi(nmv,*)
      real(8), intent(in)    :: vr (nmv,*), vi (nmv,*)
      real(8), intent(in)    :: vyr(nmv,*), vyi(nmv,*)
      integer, intent(in)    :: nma, nmv, m_size, i_base
!----
      integer :: j, k0, n

      integer, parameter :: nvec = 2128
      integer, parameter :: blas_chunk = 24
      integer, parameter :: blas_nvect = 1024

      include 'trd.h'

      integer :: i_1,i_2,i_3,i_4
      integer :: j_1,j_2,j_3,j_4
      integer :: ii_1,ii_2,ii_3,ii_4
      integer :: jjj_1,jjj_2,jjj_3,jjj_4
      integer :: ii_step
!----
      integer :: local_rank, local_size
!----
      local_rank = 0
      local_size = 1
!$    local_rank = omp_get_thread_num()
!$    local_size = omp_get_num_threads()
!----
      n = get_loop_end  (i_base, size_of_row,my_row)
      k0 = translate_l2g(n, size_of_row,my_row)

      jjj_2 = 1 ! get_loop_start(1,      size_of_col, my_col)
      jjj_3 = get_loop_end  (k0, size_of_col,my_col)
      ii_step = 0
      do jjj_1 = jjj_2, jjj_3, blas_nvect
         jjj_4 = min(jjj_1+blas_nvect-1, jjj_3)
         k0   = translate_l2g(jjj_1, size_of_col, my_col)
         ii_2 = get_loop_start(k0,     size_of_row, my_row)
         ii_2 = max(jjj_2,ii_2)
         ii_3 = n

         do i_1=ii_2,ii_3,blas_chunk
            i_4=min(ii_3-i_1+1,blas_chunk)
            i_3 = i_1+i_4-1
            j   = translate_l2g(i_3, size_of_row, my_row)
            j_4 = get_loop_end(j, size_of_col, my_col)
            j_4 = min(j_4,jjj_4)-jjj_1+1

            ii_step = ii_step+1
            if(mod(ii_step-1,local_size)==local_rank)then

               call dgemm('n','t',
     &                    j_4, i_4, m_size,
     &                    -1.0d+00, ur (jjj_1,1), nmv,
     &                    vyr(i_1,1), nmv,
     &                    1.0d+00, ar_alias_0(jjj_1, i_1), nma)
               call dgemm('n','t',
     &                    j_4, i_4, m_size,
     &                    -1.0d+00, ui (jjj_1,1), nmv,
     &                    vyi(i_1,1), nmv,
     &                    1.0d+00, ar_alias_0(jjj_1, i_1), nma)
               call dgemm('n','t',
     &                    j_4, i_4, m_size,
     &                    -1.0d+00, vr (jjj_1,1), nmv,
     &                    uyr(i_1,1), nmv,
     &                    1.0d+00, ar_alias_0(jjj_1,i_1), nma)
               call dgemm('n','t',
     &                    j_4, i_4, m_size,
     &                    -1.0d+00, vi (jjj_1,1), nmv,
     &                    uyi(i_1,1), nmv,
     &                    1.0d+00, ar_alias_0(jjj_1,i_1), nma)

            end if

            ii_step = ii_step+1
            if(mod(ii_step-1,local_size)==local_rank)then

               call dgemm('n','t',
     &                    j_4, i_4, m_size,
     &                    -1.0d+00, ui (jjj_1,1), nmv,
     &                    vyr(i_1,1), nmv,
     &                    1.0d+00, ai_alias_0(jjj_1,i_1), nma)
               call dgemm('n','t',
     &                    j_4, i_4, m_size,
     &                    1.0d+00, ur (jjj_1,1), nmv,
     &                    vyi(i_1,1), nmv,
     &                    1.0d+00, ai_alias_0(jjj_1,i_1), nma)
               call dgemm('n','t',
     &                    j_4, i_4, m_size,
     &                    -1.0d+00, vi (jjj_1,1), nmv,
     &                    uyr(i_1,1), nmv,
     &                    1.0d+00, ai_alias_0(jjj_1,i_1), nma)
               call dgemm('n','t',
     &                    j_4, i_4, m_size,
     &                    1.0d+00, vr (jjj_1,1), nmv,
     &                    uyi(i_1,1), nmv,
     &                    1.0d+00, ai_alias_0(jjj_1,i_1), nma)

            end if

         end do

      end do
!----
      return
      end subroutine
!======================================================================!
