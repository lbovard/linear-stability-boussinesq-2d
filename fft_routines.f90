module fft_routines
        use globals
        use, intrinsic :: iso_c_binding
        implicit none
        include "fftw3.f03"
        type(C_PTR) :: forward, inverse

        real(kind=8),  dimension(N,N) :: func 
        complex(C_DOUBLE_COMPLEX), dimension(Nr,N) :: func_hat
contains
        ! initialises fft routines
        subroutine init_fft()
                forward=fftw_plan_dft_r2c_2d(N,N,func,func_hat,FFTW_ESTIMATE)
                inverse=fftw_plan_dft_c2r_2d(N,N,func_hat,func,FFTW_ESTIMATE)
        end subroutine init_fft
        
        ! 2d unaliased fft
        subroutine fft2(A,A_hat) 
                real(kind=8), intent(inout) :: A(:,:)
                complex(C_DOUBLE_COMPLEX), intent(inout) :: A_hat(:,:)
                call fftw_execute_dft_r2c(forward,A,A_hat)
        end subroutine fft2
       
        ! 2d aliased fft 
        subroutine afft2(A,A_hat) 
                real(kind=8), intent(inout) :: A(:,:)
                complex(C_DOUBLE_COMPLEX), intent(inout) :: A_hat(:,:)
                call fftw_execute_dft_r2c(forward,A,A_hat)
                A_hat=A_hat*cut
        end subroutine afft2
        
        ! 2d ifft
        subroutine ifft2(A_hat,A) 
                real(kind=8), intent(inout) :: A(:,:)
                complex(C_DOUBLE_COMPLEX), intent(inout) :: A_hat(:,:)
                call fftw_execute_dft_c2r(inverse,A_hat,A)
        end subroutine ifft2

        ! normalise fft since FFTW does not
        subroutine normalise(A)
                real(kind=8), intent(inout) :: A(:,:)
                A=A/real(N*N,8) 
        end subroutine normalise 
                
end module fft_routines
