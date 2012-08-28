module display
        use globals
        use fft_routines
        implicit none
        
contains 
        subroutine print_fft(A)
                complex(C_DOUBLE_COMPLEX), intent(in) :: A(:,:)
                print *, A
        end subroutine print_fft
        
        subroutine print_mat(A)
                real(kind=8), intent(in) :: A(:,:)
                integer :: i
                do i=1,N
                        print '(50f8.4)', A(i,1:N) 
                end do 
        end subroutine print_mat
        subroutine print_fft_matrix(A)
                implicit none
                complex(C_DOUBLE_COMPLEX), intent(in) :: A((N/2+1),N)
                integer :: i,j
                do i=1,Nr
                        print '(50f8.4)', A(i,1:N) 
                end do
        end subroutine print_fft_matrix
end module display
