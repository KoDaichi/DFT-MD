       subroutine matrix_set(a,n)
       implicit double precision(a-h,o-z),integer(i-n)

       real(8)   :: a(n,n)

       do i=1,n
          do j=1,n
             a(j,i) = (n+1-Max(n+1-i,n+1-j))*1.0D+00
          end do
       end do

       return
       end

       subroutine matrix_seti(a,n)
       implicit double precision(a-h,o-z),integer(i-n)

       real(8)   :: a(n,n)

       do i=1,n
          do j=1,n
             if(j<i)then
               a(j,i) = (n+1-Max(n+1-i,n+1-j))*1.0D+00
             else if(j==i)then
               a(j,i) = 0.0D+00
             else
               a(j,i) =-(n+1-Max(n+1-i,n+1-j))*1.0D+00
             endif
          end do
       end do

       return
       end

