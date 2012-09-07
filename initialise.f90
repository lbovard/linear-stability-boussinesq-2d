! Module contains routines involved in initialising, allocating, deallocating
! all the various matrices and arrays that are needed for simulation
!

module initialise       
        use globals
        use display
        implicit none
       ! include "fftw3.f03"
contains
        ! allocate the memory  
        subroutine alloc_matrices()
                allocate(kx(N,N),ky(N,N),kz(N,N))
                allocate(k_sq(N,N))
                allocate(kinv_sq(N,N))
                allocate(cut(N,N))
                !underlying mesh
                allocate(X(N,N),Y(N,N))
                ! initial omega and velocity values
                allocate(om_i(N,N),u_0(N,N),v_0(N,N))
                !projection tensors
                allocate(p11(N,N),p12(N,N),p13(N,N),p21(N,N),p22(N,N),p23(N,N),p31(N,N),p32(N,N),p33(N,N))
                ! integrating factor 
                allocate(iuu_hat(N,N),ivv_hat(N,N),iww_hat(N,N),irho_hat(N,N))
                ! time step variables
                allocate(iuu_hat_new(N,N),ivv_hat_new(N,N),iww_hat_new(N,N),irho_hat_new(N,N))
                ! right side variables
               ! allocate(rr(N,N),ur(N,N),vr(N,N),wr(N,N))
                allocate(rr(N,N),ur(N,N),vr(N,N),wr(N,N))
                allocate(rr_old(N,N),ur_old(N,N),vr_old(N,N),wr_old(N,N))
                allocate(growth_rate(N))
        end subroutine alloc_matrices

        ! iniitialise wavenumber matrices
        subroutine init_wn()
                integer :: i,j 
                kz=ukz 
                forall(i=1:N,j=1:N/2+1) kx(i,j)=cmplx(tpiL,0,8)*cmplx(j-1,0,8)
                forall(i=1:N,j=N/2+2:N) kx(i,j)=cmplx(tpiL,0,8)*cmplx(j-N-1,0,8)
                ky=transpose(kx)
                k_sq=kx*kx+ky*ky+kz*kz
                kinv_sq=1.0_8/k_sq
                !dealiasing information using 2/3s rule
                n_k=ceiling(2.0/3.0*N)
                if(mod(n_k,2)==1) then
                        n_k=n_k-1
                end if
                n_k=(N-n_k)/2
                cut=cmplx(1.0,0,8)
                forall(i=1:N,j=(N/2-n_k+1):N/2) cut(i,j)=0
                forall(i=1:N,j=N/2+2:(N/2+n_k+1)) cut(i,j)=0
                cut=cut*transpose(cut)
        end subroutine init_wn

        !initialise the grid 
        subroutine init_grid()
                real, dimension(N) :: xx
                integer ::  i,j
                do i=1,N
                        xx(i)=-L/2+dx*real(i)
                end do
                forall(i=1:N,j=1:N) X(i,j)=xx(j)
                Y=transpose(X)
        end subroutine init_grid

        !initialise the projection tensor        
        subroutine init_projection() 
                p11=1-kx*kx*kinv_sq
                p22=1-ky*ky*kinv_sq
                p33=1-kz*kz*kinv_sq
                p12=-kx*ky*kinv_sq
                p13=-kx*kz*kinv_sq
                p23=-ky*kz*kinv_sq
                p32=p23
                p21=p12
                p31=p13
        end subroutine init_projection
       
        subroutine init_velocities()
                real(kind=8), dimension(:,:), allocatable  ::uur,vvr,wwr
                allocate(uur(N,N),vvr(N,N),wwr(N,N))
                call random_number(uur)
                call random_number(vvr)
                call random_number(wwr)
                uu=cmplx(uur,0,8)
                vv=cmplx(vvr,0,8)
                ww=cmplx(wwr,0,8)
        end subroutine init_velocities 
        
        !deallocate memory
        subroutine dealloc_matrices()
                deallocate(kx,ky,kz)
                deallocate(k_sq)
                deallocate(kinv_sq)
                deallocate(cut)
                deallocate(X,Y)
                deallocate(om_i,u_0,v_0)
                deallocate(p11,p12,p13,p21,p22,p23,p31,p32,p33)
                deallocate(iuu_hat,ivv_hat,iww_hat,irho_hat)
                deallocate(iuu_hat_new,ivv_hat_new,iww_hat_new,irho_hat_new)
                deallocate(rr,ur,vr,wr)
                deallocate(rr_old,ur_old,vr_old,wr_old)
                deallocate(growth_rate)
        end subroutine dealloc_matrices 
       
        
end module initialise
