program main
        use globals
        use initialise
        use fft_routines
        use display
        use solver
        implicit none
        integer :: i,j
       ! num_steps=floor(t_final/dt)
        num_steps=1
        call alloc_matrices()
        call alloc_fft()
        call init_fft()
        call init_wn()
        call init_grid()
        call init_projection()
        call init_velocities()
        call init_vorticity()

        call afft2(uu,uu_hat)
        call afft2(vv,vv_hat)
        call afft2(ww,ww_hat)
        call afft2(rho,rho_hat)

        print *, 'w/o div'
        call print_fft_matrix(uu_hat)
        
        call mat_w2f(u_0,"u_0.dat",N)
        call mat_w2f(v_0,"v_0.dat",N)
        call mat_w2f(om_i,"om_i.dat",N)
        call mat_w2f(uu,"uu.dat",N)
        call mat_w2f(vv,"vv.dat",N)
        call mat_w2f(ww,"ww.dat",N)
        ! ensures divergence free
        uu_hat=uu_hat*p11+vv_hat*p12+ww_hat*p13
        print *, 'w/ div'
        call print_fft_matrix(uu_hat)
        vv_hat=uu_hat*p21+vv_hat*p22+ww_hat*p23
        ww_hat=uu_hat*p31+vv_hat*p32+ww_hat*p33
        irho_hat=rho_hat*exp(k_sq*t/Re/Sc)
        iuu_hat=uu_hat*exp(k_sq*t/Re)
        ivv_hat=vv_hat*exp(k_sq*t/Re)
        iww_hat=ww_hat*exp(k_sq*t/Re)
        call rho_right() 
        irho_hat_new=irho_hat+dt*rr
        call vel_right() 
        iuu_hat_new=iuu_hat+dt*ur
        ivv_hat_new=ivv_hat+dt*vr
        iww_hat_new=iww_hat+dt*wr


        ! update scheme
        iuu_hat=iuu_hat_new
        ivv_hat=ivv_hat_new
        iww_hat=iww_hat_new
        irho_hat=irho_hat_new
        rr_old=rr
        ur_old=ur
        vr_old=vr
        wr_old=wr
!
!        do i=1,num_steps
!                t=dt*real(i,8)
!
!                call rho_right()
!                irho_hat_new=irho_hat+1.5_8*dt*rr-0.5_8*dt*rr_old
!                call vel_right()
!                iuu_hat_new=iuu_hat+1.5_8*dt*ur-0.5_8*dt*ur_old
!                ivv_hat_new=ivv_hat+1.5_8*dt*vr-0.5_8*dt*vr_old
!                iww_hat_new=iww_hat+1.5_8*dt*wr-0.5_8*dt*wr_old
!        
!                !update scheme
!                iuu_hat=iuu_hat_new
!                ivv_hat=ivv_hat_new
!                iww_hat=iww_hat_new
!                irho_hat=irho_hat_new
!                rr_old=rr
!                ur_old=ur
!                vr_old=vr
!                wr_old=wr
!                call energy()
!                growth_rate(i)=en
!                print *, en
!        end do 
        call dealloc_matrices()
        call dealloc_fft()

end program main


subroutine mat_w2f(mat,fname,N)
        integer, intent(in) :: N
        real(kind=8), intent(in), dimension(N,N) :: mat 
        character (len=*), intent(in) :: fname
        integer :: i
        open(unit=1,file=fname)
        do i=1,N
                do j=1,N
                write(1,*), mat(i,j)
                end do
        end do
        close(1)
end subroutine mat_w2f    
