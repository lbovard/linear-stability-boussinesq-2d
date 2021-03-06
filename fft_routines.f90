module fft_routines
        use globals
        use, intrinsic :: iso_c_binding
        implicit none
        include "fftw3.f03"

        ! fftw variables
        type(C_PTR) :: forward, inverse

contains
        ! allocate fftw variables
        subroutine alloc_fft()
                up=fftw_alloc_complex(int(N*N,C_SIZE_T))
                vp=fftw_alloc_complex(int(N*N,C_SIZE_T))
                wp=fftw_alloc_complex(int(N*N,C_SIZE_T))
                rhop=fftw_alloc_complex(int(N*N,C_SIZE_T))
                om1p=fftw_alloc_complex(int(N*N,C_SIZE_T))
                om2p=fftw_alloc_complex(int(N*N,C_SIZE_T))
                om3p=fftw_alloc_complex(int(N*N,C_SIZE_T))
                rp1=fftw_alloc_complex(int(N*N,C_SIZE_T))
                rp2=fftw_alloc_complex(int(N*N,C_SIZE_T))
                rp3=fftw_alloc_complex(int(N*N,C_SIZE_T))
                rp4=fftw_alloc_complex(int(N*N,C_SIZE_T))
                A1p=fftw_alloc_complex(int(N*N,C_SIZE_T))
                B1p=fftw_alloc_complex(int(N*N,C_SIZE_T))
                C1p=fftw_alloc_complex(int(N*N,C_SIZE_T))
                uq=fftw_alloc_complex(int(N*N,C_SIZE_T))
                vq=fftw_alloc_complex(int(N*N,C_SIZE_T))
                wq=fftw_alloc_complex(int(N*N,C_SIZE_T))
                rhoq=fftw_alloc_complex(int(N*N,C_SIZE_T))
                om1hp=fftw_alloc_complex(int(N*N,C_SIZE_T))
                om2hp=fftw_alloc_complex(int(N*N,C_SIZE_T))
                om3hp=fftw_alloc_complex(int(N*N,C_SIZE_T))
                rp1_hat=fftw_alloc_complex(int(N*N,C_SIZE_T))
                rp2_hat=fftw_alloc_complex(int(N*N,C_SIZE_T))
                rp3_hat=fftw_alloc_complex(int(N*N,C_SIZE_T))
                rp4_hat=fftw_alloc_complex(int(N*N,C_SIZE_T))
                A1p_hat=fftw_alloc_complex(int(N*N,C_SIZE_T))
                B1p_hat=fftw_alloc_complex(int(N*N,C_SIZE_T))
                C1p_hat=fftw_alloc_complex(int(N*N,C_SIZE_T))
                call c_f_pointer(up,ww,[N,N])
                call c_f_pointer(vp,vv,[N,N])
                call c_f_pointer(wp,uu,[N,N])
                call c_f_pointer(rhop,rho,[N,N])
                call c_f_pointer(om1p,om1,[N,N])
                call c_f_pointer(om2p,om2,[N,N])
                call c_f_pointer(om3p,om3,[N,N])
                call c_f_pointer(rp1,r1,[N,N])
                call c_f_pointer(rp2,r2,[N,N])
                call c_f_pointer(rp3,r3,[N,N])
                call c_f_pointer(rp4,r4,[N,N])
                call c_f_pointer(A1p,A_1,[N,N])
                call c_f_pointer(B1p,B_1,[N,N])
                call c_f_pointer(C1p,C_1,[N,N])
                call c_f_pointer(uq,uu_hat,[N,N])
                call c_f_pointer(vq,vv_hat,[N,N])
                call c_f_pointer(wq,ww_hat,[N,N])
                call c_f_pointer(rhoq,rho_hat,[N,N])
                call c_f_pointer(om1hp,om1_hat,[N,N])
                call c_f_pointer(om2hp,om2_hat,[N,N])
                call c_f_pointer(om3hp,om3_hat,[N,N])
                call c_f_pointer(rp1_hat,r1_hat,[N,N])
                call c_f_pointer(rp2_hat,r2_hat,[N,N])
                call c_f_pointer(rp3_hat,r3_hat,[N,N])
                call c_f_pointer(rp4_hat,r4_hat,[N,N])
                call c_f_pointer(A1p_hat,A_1_hat,[N,N])
                call c_f_pointer(B1p_hat,B_1_hat,[N,N])
                call c_f_pointer(C1p_hat,C_1_hat,[N,N])
        end subroutine alloc_fft

        ! initialises fft routines
        subroutine init_fft()
                forward=fftw_plan_dft_2d(N,N,uu,uu_hat,FFTW_FORWARD,FFTW_ESTIMATE)
                inverse=fftw_plan_dft_2d(N,N,uu_hat,uu,FFTW_BACKWARD,FFTW_ESTIMATE)
        end subroutine init_fft
        
        ! 2d unaliased fft
        subroutine fft2(A,A_hat) 
                complex(C_DOUBLE_COMPLEX), pointer, intent(inout) :: A(:,:)
                complex(C_DOUBLE_COMPLEX), pointer,  intent(inout) :: A_hat(:,:)
                call fftw_execute_dft(forward,A,A_hat)
        end subroutine fft2
       
        ! 2d aliased fft 
        subroutine afft2(A,A_hat) 
                complex(C_DOUBLE_COMPLEX), pointer, intent(inout) :: A(:,:)
                complex(C_DOUBLE_COMPLEX), pointer,  intent(inout) :: A_hat(:,:)
                call fftw_execute_dft(forward,A,A_hat)
                A_hat=A_hat*cut
        end subroutine afft2
        
        ! 2d ifft
        subroutine ifft2(A_hat,A) 
                complex(C_DOUBLE_COMPLEX), pointer, intent(inout) :: A(:,:)
                complex(C_DOUBLE_COMPLEX), pointer, intent(inout) :: A_hat(:,:)
                call fftw_execute_dft(inverse,A_hat,A)
                A=A/cmplx(N*N,0,8) !normalise
        end subroutine ifft2

        !clean up and free memory
        subroutine dealloc_fft()
                call fftw_destroy_plan(forward)
                call fftw_destroy_plan(inverse)
        end subroutine dealloc_fft

                
end module fft_routines
