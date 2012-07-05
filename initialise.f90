!initialise variables
module initialise       
        use globals
        implicit none
contains 
        ! iniitialise kz matrices
        subroutine init_wn(ukz)
                real, intent(in) :: ukz
                integer i,j 
                allocate(kx(1:(N/2+1),1:N))
                do j=1,(N/2+1)
                        do i=1,N
                                kx(i,j)=real(i)
                        end do
                end do
        end subroutine init_wn
                
end module initialise
