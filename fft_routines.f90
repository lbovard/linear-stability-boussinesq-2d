module fft_routines
        use globals
        use, intrinsic :: iso_c_binding
        implicit none
        include "fftw3.f03"

        ! fftw variables
        type(C_PTR) :: forward, inverse
        real(C_DOUBLE), pointer :: func(:,:)
        complex(C_DOUBLE_COMPLEX), pointer :: func_hat(:,:)

contains
        ! allocate fftw variables
        subroutine alloc_fft()
                up=fftw_alloc_real(int(N*N,C_SIZE_T))
                vp=fftw_alloc_real(int(N*N,C_SIZE_T))
                wp=fftw_alloc_real(int(N*N,C_SIZE_T))
                rhop=fftw_alloc_real(int(N*N,C_SIZE_T))
                funcp=fftw_alloc_real(int(N*N,C_SIZE_T))
                uq=fftw_alloc_complex(int(Nr*N,C_SIZE_T))
                vq=fftw_alloc_complex(int(Nr*N,C_SIZE_T))
                wq=fftw_alloc_complex(int(Nr*N,C_SIZE_T))
                rhoq=fftw_alloc_complex(int(Nr*N,C_SIZE_T))
                funcq=fftw_alloc_complex(int(Nr*N,C_SIZE_T))
                call c_f_pointer(up,ww,[n,n])
                call c_f_pointer(vp,vv,[N,N])
                call c_f_pointer(wp,uu,[N,N])
                call c_f_pointer(rhop,rho,[N,N])
                call c_f_pointer(funcp,func,[N,N])
                call c_f_pointer(uq,uu_hat,[Nr,N])
                call c_f_pointer(vq,vv_hat,[Nr,N])
                call c_f_pointer(wq,ww_hat,[Nr,N])
                call c_f_pointer(rhoq,rho_hat,[Nr,N])
                call c_f_pointer(funcq,rho_hat,[Nr,N])
        end subroutine alloc_fft

        ! initialises fft routines
        subroutine init_fft()
                forward=fftw_plan_dft_r2c_2d(N,N,func,func_hat,FFTW_ESTIMATE)
                inverse=fftw_plan_dft_c2r_2d(N,N,func_hat,func,FFTW_ESTIMATE)
        end subroutine init_fft
        
        ! 2d unaliased fft
        subroutine fft2(A,A_hat) 
                real(C_DOUBLE), intent(inout) :: A(:,:)
                complex(C_DOUBLE_COMPLEX), intent(inout) :: A_hat(:,:)
                call fftw_execute_dft_r2c(forward,A,A_hat)
        end subroutine fft2
       
        ! 2d aliased fft 
        subroutine afft2(A,A_hat) 
                real(C_DOUBLE), intent(inout) :: A(:,:)
                complex(C_DOUBLE_COMPLEX), intent(inout) :: A_hat(:,:)
                call fftw_execute_dft_r2c(forward,A,A_hat)
                A_hat=A_hat*cut
        end subroutine afft2
        
        ! 2d ifft
        subroutine ifft2(A_hat,A) 
                real(C_DOUBLE), intent(inout) :: A(:,:)
                complex(C_DOUBLE_COMPLEX), intent(inout) :: A_hat(:,:)
                call fftw_execute_dft_c2r(inverse,A_hat,A)
        end subroutine ifft2

        ! normalise fft since FFTW does not
        subroutine normalise(A)
                real(C_DOUBLE), intent(inout) :: A(:,:)
                A=A/real(N*N,8) 
        end subroutine normalise 

        !clean up and free memory
        subroutine dealloc_fft()
                call fftw_destroy_plan(forward)
                call fftw_destroy_plan(inverse)
        end subroutine dealloc_fft

                
end module fft_routines
