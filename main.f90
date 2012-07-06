program main
        use globals
        use initialise 
        implicit none
        integer :: i,j
        real :: ukz=0.1
        call alloc_matrices()
        call init_wn(ukz)
        call init_grid()
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
        real, intent(in) :: A(n,(n/2+1))
        integer :: i,j
        do i=1,n
                print*, A(i,:)
        end do
end 
