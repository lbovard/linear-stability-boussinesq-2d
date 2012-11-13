! Module contains routines involved in initialising, allocating, deallocating
! all the various matrices and arrays that are needed for simulation
!
#define HYPERVIS 1 
module initialise       
        use globals
        use display
        implicit none
       ! include "fftw3.f03"
contains
        ! allocate the memory  
        subroutine alloc_matrices()
                allocate(kx(N,N),ky(N,N),kz(N,N),k_sq(N,N),kinv_sq(N,N),cut(N,N))
                !underlying mesh
                allocate(X(N,N),Y(N,N))
                ! initial omega and velocity values
                allocate(om_i(N,N),u_0(N,N),v_0(N,N))
                allocate(uu_hat_temp(N,N),vv_hat_temp(N,N),ww_hat_temp(N,N)) 
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
                allocate(growth_rate(num_steps))
                allocate(tote(num_steps))
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
                if (hypervis==1) then
                        k_sq=k_sq*k_sq
                end if
                !dealiasing information using 2/3s rule
                n_k=ceiling(2.0/3.0*N)
                if(mod(n_k,2)==1) then
                        n_k=n_k-1
                end if
                n_k=(N-n_k)/2
                cut=cmplx(1.0,0,8)
                forall(i=1:N,j=(N/2-n_k+1):(N/2+n_k+1)) cut(i,j)=0
!                forall(i=1:N,j=N/2+2:(N/2+n_k+1)) cut(i,j)=0
                cut=cut*transpose(cut)
        end subroutine init_wn

        !initialise the grid 
        subroutine init_grid()
                real(kind=8), dimension(N) :: xx
                integer ::  i,j
                do i=1,N
                        xx(i)=-L/2.0_8+dx*real(i,8)
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
       
        subroutine init_conditions(cont,kzs)
                integer, intent(in) :: cont
                real(kind=8), dimension(:,:), allocatable  ::uur,vvr,wwr
                character(len=30) :: fname
                character(len=30), intent(in) :: kzs
                !if new run generate new initial conditions
                if(cont==0) then
                        allocate(uur(N,N),vvr(N,N),wwr(N,N))
                        call random_number(uur)
                        call random_number(vvr)
                        call random_number(wwr)
                        uu=cmplx(uur,0,8)
                        vv=cmplx(vvr,0,8)
                        ww=cmplx(wwr,0,8)
                        rho=cmplx(0,0,8)
                else 
                        fname='k_z.'//trim(kzs)//'.u.dat'
                        call mat_r2f_c(uu,fname,N)
                        fname='k_z.'//trim(kzs)//'.v.dat'
                        call mat_r2f_c(vv,fname,N)
                        fname='k_z.'//trim(kzs)//'.w.dat'
                        call mat_r2f_c(ww,fname,N)
                        fname='k_z.'//trim(kzs)//'.rho.dat'
                        call mat_r2f_c(rho,fname,N)
                end if
        end subroutine init_conditions
        
        !deallocate memory
        subroutine dealloc_matrices()
        end subroutine dealloc_matrices 
       
        
end module initialise
