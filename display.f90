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

        subroutine mat_w2f(mat,fname,N)
                integer, intent(in) :: N
                real(kind=8), intent(in), dimension(N,N) :: mat 
                character (len=*), intent(in) :: fname
                integer :: i,j
                open(unit=1,file=fname)
                do i=1,N
                        do j=1,N
                        write(1,*), mat(i,j)
                        end do
                end do
                close(1)
        end subroutine mat_w2f  

        subroutine mat_w2f_c(mat,fname,N)
                integer, intent(in) :: N
                complex(C_DOUBLE_COMPLEX), intent(in), dimension((N/2+1),N) :: mat 
                character (len=*), intent(in) :: fname
                integer :: i,j 
                open(unit=1,file=fname)
                do i=1,(N/2+1)
                        do j=1,N
                                write(1,*), real(mat(i,j))
                                write(1,*), imag(mat(i,j))
                        end do
                end do
                close(1)
        end subroutine mat_w2f_c

        subroutine mat_w2f_fft(mat,fname,N)
                integer, intent(in) :: N
                real(kind=8), intent(in), dimension(N/2+1,N) :: mat 
                character (len=*), intent(in) :: fname
                integer :: i,j
                open(unit=1,file=fname)
                do i=1,(N/2+1)
                        do j=1,N
                        write(1,*), mat(i,j)
                        end do
                end do
                close(1)
        end subroutine mat_w2f_fft
end module display
