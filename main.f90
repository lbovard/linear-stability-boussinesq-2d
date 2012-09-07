program main
        use globals
        use initialise
        use fft_routines
        use display
        use solver
        implicit none
        integer :: i,j
       ! num_steps=floor(t_final/dt)
        call alloc_matrices()
        call alloc_fft()
        call init_fft()
        call init_wn()
        call init_grid()
        call init_projection()
        call init_velocities()
        call init_vorticity()

        call mat_w2f_c(uu,"uu_init.dat",N)
        call mat_w2f_c(vv,"vv_init.dat",N)
        call mat_w2f_c(ww,"ww_init.dat",N)

        call afft2(uu,uu_hat)
        call afft2(vv,vv_hat)
        call afft2(ww,ww_hat)
        
        call init_fft()


        call afft2(rho,rho_hat)
        uu_hat_temp=uu_hat
        vv_hat_temp=vv_hat
        ww_hat_temp=ww_hat
        ! ensures divergence free
        uu_hat=uu_hat_temp*p11+vv_hat_temp*p12+ww_hat_temp*p13
        vv_hat=uu_hat_temp*p21+vv_hat_temp*p22+ww_hat_temp*p23
        ww_hat=uu_hat_temp*p31+vv_hat_temp*p32+ww_hat_temp*p33
        !convert to integrating factors
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

!
        ! update scheme
        iuu_hat=iuu_hat_new
        ivv_hat=ivv_hat_new
        iww_hat=iww_hat_new
        irho_hat=irho_hat_new
        rr_old=rr
        ur_old=ur
        vr_old=vr
        wr_old=wr

        num_steps=floor(t_final/dt)
!        num_steps=10
        print *, num_steps
        do i=1,num_steps
                print *, i
                t=dt*cmplx(i,0,8)

                call rho_right()
                irho_hat_new=irho_hat+1.5_8*dt*rr-0.5_8*dt*rr_old
                call vel_right()
                iuu_hat_new=iuu_hat+1.5_8*dt*ur-0.5_8*dt*ur_old
                ivv_hat_new=ivv_hat+1.5_8*dt*vr-0.5_8*dt*vr_old
                iww_hat_new=iww_hat+1.5_8*dt*wr-0.5_8*dt*wr_old
        
                !update scheme
                iuu_hat=iuu_hat_new
                ivv_hat=ivv_hat_new
                iww_hat=iww_hat_new
                irho_hat=irho_hat_new
                rr_old=rr
                ur_old=ur
                vr_old=vr
                wr_old=wr
                call energy()
                growth_rate(i)=en
        end do 
!        call dealloc_matrices()
!        call dealloc_fft()

        call mat_w2f_c(irho_hat,"irho_hat.dat",N)        
        call mat_w2f_c(iuu_hat,"iuu_hat.dat",N)        
        call mat_w2f_c(ivv_hat,"ivv_hat.dat",N)        
        call mat_w2f_c(iww_hat,"iww_hat.dat",N)        
        call w2f(growth_rate,"energy.dat",num_steps)
end program main

