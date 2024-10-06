       subroutine matrix_adjust_h(a,n,b,lda)

      use communication_h, only : loop_c
     &               , translate_l2g
     &               , eigen_init
     &               , eigen_free

       implicit double precision(a-h,o-z),integer(i-n)
       real(8)  :: a(n,n)
       real(8)  :: b(1:lda,*)

       include 'mpif.h'
       include 'trd.h'

       call eigen_init(1)

       i_2=1; i_3=n
       call loop_c(i_2,i_3,i_4,n,size_of_row,my_row)
       j_2=1; j_3=n
       call loop_c(j_2,j_3,j_4,n,size_of_col,my_col)

       do i_1=i_2,i_3
          i = translate_l2g(i_1,size_of_row,my_row)
          do j_1=j_2,j_3
             j = translate_l2g(j_1,size_of_col,my_col)
             b(j_1,i_1) = a(j,i)
          end do
       end do

       call eigen_free(0)

       return
       end

       subroutine matrix_adjust_hi(a,n,b,lda)

      use communication_h, only : loop_c
     &               , translate_l2g
     &               , eigen_init
     &               , eigen_free

       implicit double precision(a-h,o-z),integer(i-n)
       real(8)  :: a(n,n)
       real(8)  :: b(1:lda,*)

       include 'mpif.h'
       include 'trd.h'

       call eigen_init(1)

       i_2=1; i_3=n
       call loop_c(i_2,i_3,i_4,n,size_of_row,my_row)
       j_2=1; j_3=n
       call loop_c(j_2,j_3,j_4,n,size_of_col,my_col)

       do i_1=i_2,i_3
          i = translate_l2g(i_1,size_of_row,my_row)
          do j_1=j_2,j_3
             j = translate_l2g(j_1,size_of_col,my_col)
             if(j<i)then
               b(j_1,i_1) = a(j,i)
             else if(j==i)then
               b(j_1,i_1) = a(j,i)
             else
               b(j_1,i_1) = a(j,i)
             endif
          end do
       end do

       call eigen_free(0)

       return
       end

