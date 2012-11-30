
program main
        use globals
        use initialise
        use fft_routines
        use display
        use solver
        integer :: i,j,fulldump, current_time, num_dumps, time_dump
        integer :: nanspres
        !netcdf variables
        integer :: mdata_id, x_dimid,y_dimid,t_id,cmplx_id
        integer :: uid,vid,wid,rhoid
        integer :: time_rep
        integer, dimension(4) :: dimids,curr_dim, count_dim
        character(len=8) :: ct
        character(len=32) :: kzs,Ns, res, fhs
        character(len=80) :: outputname 
        
        !get kz,Fh,Re,N from command line 
        i=command_argument_count()
        if(i /= 4) then
                print *, "ERROR: insufficient cl arguments"
                call EXIT(0)
        end if
        call get_command_argument(1,kzs)
        read(kzs,'(g5.2)') ukz
        call get_command_argument(2,fhs)
        read(fhs,'(g5.2)') Fh
        call get_command_argument(3,res)
        read(res,'(g7.2)') Re
        call get_command_argument(4,Ns)
        read(Ns,'(i4)') N

        !define some important information
        dx=L/real(N,8)
        dy=dx
        num_steps=floor(t_final/dt)
        !number of times I want to dump the data
        num_dumps=50
        !how often I dump the data, in terms of num_steps
        fulldump=floor((t_final/num_dumps)/dt)
        print *, fulldump
        print *, num_steps

        !allocate
        call alloc_matrices()
        call alloc_fft()
        call init_fft()
        call init_wn()
        call init_grid()
        call init_projection()
        call init_conditions(cont,kzs)
        call init_vorticity()

        if  (hypervis == 1) then 
                Rev=Re
                Re=Re*(n_k**2)
        end if
      
        ifactor=k_sq/Re
        call afft2(uu,uu_hat)
        call afft2(vv,vv_hat)
        call afft2(ww,ww_hat)
        call afft2(rho,rho_hat)
       
        ! do NETCDF allocation
        outputname='kz.'//trim(kzs)//'.'//trim(Ns)//'.re.'//trim(res)//'.fh.'//trim(fhs)//'.nc'
        call check(nf90_create(outputname,nf90_64bit_offset,mdata_id))
        call check(nf90_def_dim(mdata_id,"X",N,x_dimid))
        call check(nf90_def_dim(mdata_id,"Y",N,y_dimid))
        !num_dumps + 2 to include first and last
        call check(nf90_def_dim(mdata_id,"t",num_dumps+2,t_id))
        call check(nf90_def_dim(mdata_id,"cmplx",2,cmplx_id))
        
        dimids=(/ x_dimid,y_dimid,t_id,cmplx_id /)
        call check(nf90_def_var(mdata_id,"u",nf90_double,dimids,uid))
        call check(nf90_def_var(mdata_id,"v",nf90_double,dimids,vid))
        call check(nf90_def_var(mdata_id,"w",nf90_double,dimids,wid))
        call check(nf90_def_var(mdata_id,"rho",nf90_double,dimids,rhoid))

        call check(nf90_enddef(mdata_id))

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

        !Euler method for first step        
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

        !dump initial data
        curr_dim=(/1,1,1,1/)
        count_dim=(/N,N,1,1/)
        call check(nf90_put_var(mdata_id,uid,real(ur),curr_dim,count_dim))
        call check(nf90_put_var(mdata_id,vid,real(vr),curr_dim,count_dim))
        call check(nf90_put_var(mdata_id,wid,real(wr),curr_dim,count_dim))
        call check(nf90_put_var(mdata_id,rhoid,real(rr),curr_dim,count_dim))
!
        curr_dim=(/1,1,1,2/)
        call check(nf90_put_var(mdata_id,uid,imag(ur),curr_dim,count_dim))
        call check(nf90_put_var(mdata_id,vid,imag(vr),curr_dim,count_dim))
        call check(nf90_put_var(mdata_id,wid,imag(wr),curr_dim,count_dim))
        call check(nf90_put_var(mdata_id,rhoid,imag(rr),curr_dim,count_dim))
        time_dump=1
        j=1
        prev_en=0._8 
        do i=1,num_steps
                tstep=i
                
                t=dt*cmplx(i,0,8) 
                time_rep=dt*cmplx(j,0,8)
                if (mod(j,1000)==0) then
                     j=0 
                end if
                j=j+1
                t_ifactor=ifactor*time_rep
                call isnan_matrix(t_ifactor,nanspres,N) 
                if(nanspres==1) then
                    print *, 'problem is in t_ifactor', tstep
                end if 
                exp_ifactor_p=exp(ifactor)
                call isnan_matrix(exp_ifactor_p,nanspres,N) 
                if(nanspres==1) then
                    print *, 'problem is in exp(t_ifactor)', tstep
                end if 
                exp_ifactor_n=exp(-ifactor)
                call isnan_matrix(exp_ifactor_n,nanspres,N) 
                if(nanspres==1) then
                    print *, 'problem is in exp(-t_ifactor)', tstep
                end if 
                !apply Adams-Basforth 2nd order time-stepping
                call rho_right()
                call vel_right()
!                call isnan_matrix(rr,nanspres,N) 
!                if(nanspres==1) then
!                    print *, 'problem in before updating rho at', tstep
!                end if 
!                call isnan_matrix(ur,nanspres,N) 
!                if(nanspres==1) then
!                    print *, 'problem in before updating ur at', tstep
!                end if 
!                call isnan_matrix(vr,nanspres,N) 
!                if(nanspres==1) then
!                    print *, 'problem in before updating vr at', tstep
!                end if 
!                call isnan_matrix(wr,nanspres,N) 
!                if(nanspres==1) then
!                    print *, 'problem in before updating wr at', tstep
!                end if 
                irho_hat_new=irho_hat+1.5_8*dt*rr-0.5_8*dt*rr_old
                iuu_hat_new=iuu_hat+1.5_8*dt*ur-0.5_8*dt*ur_old
                ivv_hat_new=ivv_hat+1.5_8*dt*vr-0.5_8*dt*vr_old
                iww_hat_new=iww_hat+1.5_8*dt*wr-0.5_8*dt*wr_old
!
!                call isnan_matrix(irho_hat_new,nanspres,N) 
!                if(nanspres==1) then
!                    print *, 'problem in time stepping irho_hat at', tstep
!                end if 
!
!                call isnan_matrix(iuu_hat_new,nanspres,N) 
!                if(nanspres==1) then
!                    print *, 'problem in time stepping iuu_hat at', tstep
!                end if
!        
!                call isnan_matrix(ivv_hat_new,nanspres,N) 
!                if(nanspres==1) then
!                    print *, 'problem in time stepping ivv_hat at', tstep
!                end if
!        
!                call isnan_matrix(iww_hat_new,nanspres,N) 
!                if(nanspres==1) then
!                    print *, 'problem in time stepping iww_hat at', tstep
!                end if
        
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
                if (isnan(en)) then
                        print *, "ERROR: NaNs detected at timestep"
                        print *, i
                        print *, 'time dump', time_dump, 'at time step', i
                        time_dump=time_dump+1
                        curr_dim=(/1,1,time_dump,1/)
                        count_dim=(/N,N,1,1/)
                        call check(nf90_put_var(mdata_id,uid,real(r1),curr_dim,count_dim))
                        call check(nf90_put_var(mdata_id,vid,real(r2),curr_dim,count_dim))
                        call check(nf90_put_var(mdata_id,wid,real(r3),curr_dim,count_dim))
                        call check(nf90_put_var(mdata_id,rhoid,real(r4),curr_dim,count_dim))
                        curr_dim=(/1,1,time_dump,2/)
                        call check(nf90_put_var(mdata_id,uid,imag(r1),curr_dim,count_dim))
                        call check(nf90_put_var(mdata_id,vid,imag(r2),curr_dim,count_dim))
                        call check(nf90_put_var(mdata_id,wid,imag(r3),curr_dim,count_dim))
                        call check(nf90_put_var(mdata_id,rhoid,imag(r4),curr_dim,count_dim))
                        outputname='kz.'//trim(kzs)//'.totE.'//trim(Ns)//'.re.'//trim(res)//'.fh.'//trim(fhs)//'.dat'
                        call arr_w2f(tote,outputname,num_steps)
                        outputname='kz.'//trim(kzs)//'.sigma.'//trim(Ns)//'.re.'//trim(res)//'.fh.'//trim(fhs)//'.dat'
                        call arr_w2f(growth_rate,outputname,num_steps)
                        ! close netcdf
                        call check(nf90_close(mdata_id))
                        ! deallocate everything
                        call dealloc_matrices()
                        call dealloc_fft()
                        call EXIT(0)
                end if

                ! every 100 time steps dump some info 
                if (mod(i,100)==0) then
                        outputname='kz.'//trim(kzs)//'_data_'//trim(Ns)//'.dat'
                        call data_w2f(growth_rate(i),outputname,i)
                end if

               if(mod(i,fulldump)==0)  then
                        time_dump=time_dump+1
                        curr_dim=(/1,1,time_dump,1/)
                        count_dim=(/N,N,1,1/)
                        call check(nf90_put_var(mdata_id,uid,real(r1),curr_dim,count_dim))
                        call check(nf90_put_var(mdata_id,vid,real(r2),curr_dim,count_dim))
                        call check(nf90_put_var(mdata_id,wid,real(r3),curr_dim,count_dim))
                        call check(nf90_put_var(mdata_id,rhoid,real(r4),curr_dim,count_dim))
                        curr_dim=(/1,1,time_dump,2/)
                        call check(nf90_put_var(mdata_id,uid,imag(r1),curr_dim,count_dim))
                        call check(nf90_put_var(mdata_id,vid,imag(r2),curr_dim,count_dim))
                        call check(nf90_put_var(mdata_id,wid,imag(r3),curr_dim,count_dim))
                        call check(nf90_put_var(mdata_id,rhoid,imag(r4),curr_dim,count_dim))
                        outputname='kz.'//trim(kzs)//'.totE.'//trim(Ns)//'.re.'//trim(res)//'.fh.'//trim(fhs)//'.dat'
                        call arr_w2f(tote,outputname,num_steps)
                        outputname='kz.'//trim(kzs)//'.sigma.'//trim(Ns)//'.re.'//trim(res)//'.fh.'//trim(fhs)//'.dat'
                        call arr_w2f(growth_rate,outputname,num_steps)
                end if
        end do 

        !return to the real space for plotting 
!        r1_hat=exp(-k_sq*t/Re)*iuu_hat
!        r2_hat=exp(-k_sq*t/Re)*ivv_hat
!        r3_hat=exp(-k_sq*t/Re)*iww_hat
!        r4_hat=exp(-k_sq*t/Re/Sc)*irho_hat
        r1_hat=exp_i_factor_n*iuu_hat
        r2_hat=exp_i_factor_n*ivv_hat
        r3_hat=exp_i_factor_n*iww_hat
        r4_hat=exp_i_factor_n*irho_hat
        call ifft2(r1_hat,r1)
        call ifft2(r2_hat,r2)
        call ifft2(r3_hat,r3)
        call ifft2(r4_hat,r4)
        curr_dim=(/1,1,num_dumps+2,1/)
        count_dim=(/N,N,1,1/)
        call check(nf90_put_var(mdata_id,uid,real(r1),curr_dim,count_dim))
        call check(nf90_put_var(mdata_id,vid,real(r2),curr_dim,count_dim))
        call check(nf90_put_var(mdata_id,wid,real(r3),curr_dim,count_dim))
        call check(nf90_put_var(mdata_id,rhoid,real(r4),curr_dim,count_dim))
        curr_dim=(/1,1,num_dumps+2,2/)
        call check(nf90_put_var(mdata_id,uid,imag(r1),curr_dim,count_dim))
        call check(nf90_put_var(mdata_id,vid,imag(r2),curr_dim,count_dim))
        call check(nf90_put_var(mdata_id,wid,imag(r3),curr_dim,count_dim))
        call check(nf90_put_var(mdata_id,rhoid,imag(r4),curr_dim,count_dim))
        outputname='kz.'//trim(kzs)//'.totE.'//trim(Ns)//'.re.'//trim(res)//'.fh.'//trim(fhs)//'.dat'
        call arr_w2f(tote,outputname,num_steps)
        outputname='kz.'//trim(kzs)//'.sigma.'//trim(Ns)//'.re.'//trim(res)//'.fh.'//trim(fhs)//'.dat'
        call arr_w2f(growth_rate,outputname,num_steps)
        ! close netcdf
        call check(nf90_close(mdata_id))
        ! deallocate everything
        call dealloc_matrices()
        call dealloc_fft()
contains
        subroutine check(status)
        integer, intent(in) :: status
        if(status/=nf90_noerr) then
                print *, trim(nf90_strerror(status))
                stop "Stopped. Error in netcdf."
        end if
        end subroutine check
end program main


