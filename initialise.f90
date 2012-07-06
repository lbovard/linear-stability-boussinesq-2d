! Module contains routines involved in initialising all the various 
! matrices and arrays that are needed
!
! Groups everything together
!

module initialise       
        use globals
        implicit none
contains
        ! allocate the memory  
        subroutine alloc_matrices()
                ! see fftw documentation for allocation sizes
                allocate(kx(N,Nr))
                allocate(ky(N,Nr))
                allocate(kz(N,Nr))
                !underlying mesh
                allocate(X(N,N))
                allocate(Y(N,N))
        end subroutine alloc_matrices
        ! iniitialise kz matrices
        subroutine init_wn(ukz)
                real, intent(in) :: ukz
                integer :: i,j 
                ! better on parallel machines
                ! initialises kx,ky matrices according to fftw storage
                forall(j=1:Nr,i=1:N) kx(i,j)=2*pi/L*real(j-1)
                forall(j=1:Nr,i=1:Nr-1) ky(i,j)=2*pi/L*real(i-1) !possible oneliner
                forall(j=1:Nr,i=Nr:N) ky(i,j)=2*pi/L*real(i-N-1)        
                kz=ukz !this vs forall, which is faster? investigate at some point
        end subroutine init_wn
        
        subroutine init_grid()
                real, dimension(N) :: x
                integer ::  i,j
                do i=1,N
                        x(i)=-L/2+dx*real(i)
                        print *, x(i)
                end do
        end subroutine init_grid
        subroutine dealloc_matrices()
                deallocate(kx)
                deallocate(ky)
                deallocate(kz)
                deallocate(X)
                deallocate(Y)
        end subroutine dealloc_matrices 
                
end module initialise
