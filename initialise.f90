! Module contains routines involved in initialising, allocating, deallocating
! all the various matrices and arrays that are needed for simulation
!

module initialise       
        use globals
        implicit none
        include "fftw3.f03"
contains
        ! allocate the memory  
        subroutine alloc_matrices()
                ! see fftw documentation for allocation sizes
                allocate(kx(Nr,N),ky(Nr,N),kz(Nr,N))
                allocate(k_sq(Nr,N))
                allocate(kinv_sq(Nr,N))
                allocate(cut(Nr,N))
                !underlying mesh
                allocate(X(N,N),Y(N,N))
                ! initial omega and velocity values
                allocate(om_i(N,N),u_0(N,N),v_0(N,N))
                !projection tensors
                allocate(p11(Nr,N),p12(Nr,N),p13(Nr,N),p21(Nr,N),p22(Nr,N),p23(Nr,N),p31(Nr,N),p32(Nr,N),p33(Nr,N))
                ! integrating factor 
                allocate(iuu_hat(Nr,N),ivv_hat(Nr,N),iww_hat(Nr,N),irho_hat(Nr,N))
                ! time step variables
                allocate(iuu_hat_new(Nr,N),ivv_hat_new(Nr,N),iww_hat_new(Nr,N),irho_hat_new(Nr,N))
                ! right side variables
               ! allocate(rr(N,N),ur(N,N),vr(N,N),wr(N,N))
                allocate(rr(Nr,N),ur(Nr,N),vr(Nr,N),wr(Nr,N))
                allocate(rr_old(N,N),ur_old(N,N),vr_old(N,N),wr_old(N,N))
                allocate(growth_rate(N))
        end subroutine alloc_matrices

        ! iniitialise wavenumber matrices
        subroutine init_wn()
                integer :: i,j 
                ! better on parallel machines
                ! initialises kx,ky matrices according to fftw storage
                kz=ukz 
                forall(i=1:Nr,j=1:N/2) kx(i,j)=tpiL*real(j-1,8)
                forall(i=1:Nr,j=N/2+2:N) kx(i,j)=tpiL*real(j-N-1,8)
                forall(i=1:Nr-1,j=1:N) ky(i,j)=tpiL*real(i-1,8)
                ky(Nr,:)=0.0_8
                k_sq=kx*kx+ky*ky+kz*kz
                kinv_sq=1.0_8/k_sq

                !dealiasing information using 2/3s rule
                n_k=ceiling(2.0/3.0*N)
                if(mod(n_k,2)==1) then
                        n_k=n_k-1
                end if
                n_k=(N-n_k)/2
                cut=1.0_8
                forall(i=1:Nr,j=(N/2-n_k+1):N/2) cut(i,j)=0.0_8
                forall(i=1:Nr,j=(N/2+2):(N/2+n_k+1)) cut(i,j)=0.0_8
                forall(i=(Nr-n_k):Nr-1, j=1:N) cut(i,j)=0.0_8
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
                call random_number( harvest=uu(1:N,1:N))
                call random_number( harvest=vv(1:N,1:N))
                call random_number( harvest=ww(1:N,1:N))
        end subroutine init_velocities 

        subroutine dealloc_matrices()
                deallocate(kx,ky,kz)
                deallocate(X,Y)
                deallocate(p11,p12,p13,p21,p22,p23,p31,p32,p33)
        end subroutine dealloc_matrices 
       
        ! test matrix u for debugging purposes 
!        subroutine test(A)
!                real(C_DOUBLE), intent(inout):: A(:,:)
!                integer :: i,j
!                forall(i=1:N,j=1:N) A(i,j)=1.0_8/real(i+j-1,8)
!                do i=1,N-1
!                        A(i,i+1)=real(i+1,8)
!                end do
!        end subroutine test
        
end module initialise
