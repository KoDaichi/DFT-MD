*
*
* NOTICE : THIS SOURCE IS LIMITTED FOR THE MICROPROCESSORS
*          WHICH ADOPT THE IEEE754 FLOATING POINT ARITHMETIC FORMAT
*          AND OPERATIONS
*
*
       real(8) function machine_epsilon()

       implicit NONE
       real(8) :: eps
       data eps /z'3CB0000000000000'/

       machine_epsilon = eps

       return
       end function  machine_epsilon
