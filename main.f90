program main
        use globals
        use initialise
        use fft_routines
        use display
        implicit none
        integer :: i,j
        call alloc_matrices()
        call alloc_fft()
        call init_fft()
        call init_wn()
        call init_grid()
        call init_projection()
        call init_velocities()

        call afft2(uu,uu_hat)
        call afft2(vv,vv_hat)
        call afft2(ww,ww_hat)
        call afft2(rho,rho_hat)
        uu_hat=uu_hat*p11+vv_hat*p12+ww_hat*p13
        vv_hat=uu_hat*p21+vv_hat*p22+ww_hat*p23
        ww_hat=uu_hat*p31+vv_hat*p32+ww_hat*p33
        irho_hat=rho_hat*exp(k_sq*t/Re/Sc)
        iuu_hat=uu_hat*exp(k_sq*t/Re)
        ivv_hat=vv_hat*exp(k_sq*t/Re)
        iww_hat=ww_hat*exp(k_sq*t/Re)
!        call print_fft_matrix(p11,N)

!        call print_fft_matrix(p31,N)
!        call print_fft_matrix(cut,N)
!        w_hat=w_hat*kx*ii
!        call ifft2(w_hat,w)
!        call normalise(w)
        print '()'
!        call print_mat(w)
        call dealloc_matrices()
        call dealloc_fft()

end program main
