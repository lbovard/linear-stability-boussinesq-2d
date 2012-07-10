program main
        use globals
        use initialise
        use fft_routines
        use display
        implicit none
        integer :: i,j
        real :: ukz=0.1 
        complex(C_DOUBLE_COMPLEX), dimension(Nr,N) :: u_hat
        call alloc_matrices()
        call init_fft()
        call init_wn(4.0)
        call init_grid()
        call init_projection()
        forall(i=1:N,j=1:N) u(i,j)=1.0_8/real(i+j-1,8)  
        do i=1,N-1
                u(i,i+1)=real(i+1,8)
        end do
!        call print_fft_matrix(p11,N)
!        call print_fft_matrix(p31,N)
!        call print_fft_matrix(cut,N)
!        call afft2(u,u_hat)
!        call print_fft(u_hat)
!        u_hat=u_hat*kx*ii
!        call ifft2(u_hat,u)
!        call normalise(u)
        call dealloc_matrices()

end program main
