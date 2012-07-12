module solver
        use globals
        use fft_routines
        implicit none 

contains 
        !compute the vorticity
        subroutine omega(uh,vh,wh)
                complex(C_DOUBLE_COMPLEX), intent(in) :: uh(:,:),vh(:,:),wh(:,:)
                om1=ii*ky*wh-ii*kz*vh
                om2=ii*kz*uh-ii*kx*wh
                om3=ii*kx*vh-ii*ky*uh
        end subroutine omega

        subroutine rho_right(wh,irho_hat,rr)
                complex(C_DOUBLE_COMPLEX), intent(in) :: wh(:,:),irho_hat(:,:)
                complex(C_DOUBLE_COMPLEX), intent(out) :: rr(:,:)
        end subroutine rho_right
end module solver
