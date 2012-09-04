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

        call mat_w2f(uu,"uu_init.dat",N)
        call mat_w2f(vv,"vv_init.dat",N)
        call mat_w2f(ww,"ww_init.dat",N)

        call afft2(uu,uu_hat)
        call afft2(vv,vv_hat)
        call afft2(ww,ww_hat)
        
!        call mat_w2f_c(uu_hat,"uu_hat.dat",N)
!        call mat_w2f_c(vv_hat,"vv_hat.dat",N)
!        call mat_w2f_c(ww_hat,"ww_hat.dat",N)


!        call ifft2(uu_hat,r1)
!        call ifft2(vv_hat,r2)
!        call ifft2(ww_hat,r3)
!        
!        call normalise(r1)
!        call normalise(r2)
!        call normalise(r3)
        call init_fft()
!        call mat_w2f(r1,"uu.dat",N)
!        call mat_w2f(r2,"vv.dat",N)
!        call mat_w2f(r3,"ww.dat",N)


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

!        call mat_w2f_c(ivv_hat_new,"vv_hat.dat",N)        
!        call mat_w2f_c(iww_hat_new,"ww_hat.dat",N)        
!
        ! update scheme
        iuu_hat=iuu_hat_new
        ivv_hat=ivv_hat_new
        iww_hat=iww_hat_new
        irho_hat=irho_hat_new
!        call mat_w2f_c(irho_hat,"pp_hat.dat",N)        
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

