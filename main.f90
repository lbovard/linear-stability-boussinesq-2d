program main
        use globals
        use initialise
        use fft_routines
        use display
        use solver
        integer :: i
        character(len=32) :: arg
        character(len=20) :: outputname 

        i=command_argument_count()
        if(i /= 1) then
                print *, "ERROR: insufficient cl arguments"
                call EXIT(0)
        end if
        call get_command_argument(i,arg)
        read(arg,'(g5.2)') ukz
        print *, outputname

        num_steps=floor(t_final/dt)
        call alloc_matrices()
        call alloc_fft()
        call init_fft()
        call init_wn()
        call init_grid()
        call init_projection()
        call init_velocities()
        call init_vorticity()

!        call mat_w2f_c(uu,"uu_init.dat",N)
!        call mat_w2f_c(vv,"vv_init.dat",N)
!        call mat_w2f_c(ww,"ww_init.dat",N)
!        call mat_w2f_c(rho,"rho_init.dat",N)

        call afft2(uu,uu_hat)
        call afft2(vv,vv_hat)
        call afft2(ww,ww_hat)
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


        ! update scheme
        iuu_hat=iuu_hat_new
        ivv_hat=ivv_hat_new
        iww_hat=iww_hat_new
        irho_hat=irho_hat_new
        rr_old=rr
        ur_old=ur
        vr_old=vr
        wr_old=wr
        print *, num_steps
        print *, ukz
        print *, dt
        print *, dx
        prev_en=0._8 
        do i=1,num_steps
                print *,i 
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
                tote(i)=en
                growth_rate(i)=(en-prev_en)/dt 
                prev_en=en  
        end do 
        r1_hat=exp(-k_sq*t/Re)*iuu_hat
        r2_hat=exp(-k_sq*t/Re)*ivv_hat
        r3_hat=exp(-k_sq*t/Re)*iww_hat
        r4_hat=exp(-k_sq*t/Re/Sc)*irho_hat
        call ifft2(r1_hat,r1)
        call ifft2(r2_hat,r2)
        call ifft2(r3_hat,r3)
        call ifft2(r4_hat,r4)
        outputname='k_z.'//trim(arg)//'.u.dat'
        call mat_w2f_c(r1,outputname,N)
        outputname='k_z.'//trim(arg)//'.v.dat'
        call mat_w2f_c(r2,outputname,N)
        outputname='k_z.'//trim(arg)//'.w.dat'
        call mat_w2f_c(r3,outputname,N)
        outputname='k_z.'//trim(arg)//'.rho.dat'
        call mat_w2f_c(r4,outputname,N)
        outputname='k_z.'//trim(arg)//'.totE.dat'
        call w2f(tote,outputname,num_steps)
        outputname='k_z.'//trim(arg)//'.sigma.dat'
        call w2f(growth_rate,outputname,num_steps)
        call dealloc_matrices()
        call dealloc_fft()
end program main

