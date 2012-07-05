program main
        use globals
        use initialise 
        implicit none

        real :: ukz=0.1
        call init_wn(ukz)
        call print_fft_matrix(kx,N)
end program main

subroutine print_sq_matrix(A,n)
        implicit none
        integer, intent(in) :: n
        real, intent(in) :: A(n,n)
        integer :: i
        do i=1,n
                print*, A(i,:)
        end do
end 

subroutine print_fft_matrix(A,n)
        implicit none
        integer, intent(in) :: n
        real, intent(in) :: A(n,n)
        integer :: i,j
        do i=1,n
                print*, A(:,i)
        end do
end 
