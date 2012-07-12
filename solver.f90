module solver
        use globals
        use fft_routines
        use display
        implicit none 

contains 
        !compute the vorticity
        subroutine omega(uh,vh,wh)
                complex(C_DOUBLE_COMPLEX), intent(in) :: uh(:,:),vh(:,:),wh(:,:)
                om1_hat=ii*ky*wh-ii*kz*vh
                om2_hat=ii*kz*uh-ii*kx*wh
                om3_hat=ii*kx*vh-ii*ky*uh
                call ifft2(om1_hat,om1) 
                call ifft2(om2_hat,om2) 
                call ifft2(om3_hat,om3) 
        end subroutine omega

        subroutine rho_right()
                r1_hat=exp(-k_sq*t/Re/Sc)*irho_hat
                call ifft2(r1_hat,r1)
                r2=u_0*r1
                r3=v_0*r2
                call afft2(r2,r2_hat)
                call afft2(r3,r3_hat)
                rr=exp(k_sq*t/Re/Sc)*(-ii*kx*r2_hat-ii*ky*r3_hat)
        end subroutine rho_right

        subroutine vel_right()
                r1_hat=exp(-k_sq*t/Re)*iuu_hat
                r2_hat=exp(-k_sq*t/Re)*ivv_hat
                r3_hat=exp(-k_sq*t/Re)*iww_hat
                r4_hat=exp(-k_sq*t/Re/Sc)*irho_hat
                call ifft2(r1_hat,r1)
                call ifft2(r2_hat,r2)
                call ifft2(r3_hat,r3)
                call ifft2(r4_hat,r4)
                call omega(r1_hat,r2_hat,r3_hat)        
                A_1=r2*om_i+v_0*om3
                B_1=-r1*om_i-v_0*om3
                C_1=u_0*om1-v_0*om1

                call afft2(A_1,A_1_hat)
                call afft2(B_1,B_1_hat)
                call afft2(C_1,C_1_hat)

                ur=exp(k_sq*t/Re)*(p11*A_1_hat+p12*B_1_hat+p13*C_1_hat)
                vr=exp(k_sq*t/Re)*(p21*A_1_hat+p22*B_1_hat+p23*C_1_hat)
                wr=exp(k_sq*t/Re)*(p31*A_1_hat+p32*B_1_hat+p33*C_1_hat)
        end subroutine vel_right

        subroutine energy()
                integer :: i,j
                en=0.0
                do i=1,N
                        do j=1,N
                                en=en+abs(r1(i,j))**2+abs(r2(i,j))**2+abs(r3(i,j))**2
                        end do
                end do

                en=0.5_8*log(en)
         end subroutine energy

        subroutine init_vorticity()      
                real(kind=8), dimension(:,:), allocatable :: rr,rr_sq,theta
                real(kind=8), parameter :: AA=-0.402759395702553_8, BB=3.83170597020751_8
                integer :: i,j
                allocate(rr(N,N),rr_sq(N,N),theta(N,N))
                rr=X*X+Y*Y
                rr_sq=sqrt(rr)
                theta=atan2(Y,X)
                do i=1,N
                        do j=1,N
                                if(rr_sq(i,j)< 1.0_8 .and. rr_sq(i,j) /= 0.0) then
                                        u_0(i,j)=2.0/(AA*BB*rr_sq(i,j)**3)*((X(i,j)**2-Y(i,j)**2)*bessel_j1(BB*rr_sq(i,j))+ & 
                                        BB*Y(i,j)**2*rr_sq(i,j)*bessel_j0(BB*rr_sq(i,j)))
                                        v_0(i,j)=-2*X(i,j)*Y(i,j)/(rr_sq(i,j)**3*BB*AA)*(-2*bessel_j1(BB*rr_sq(i,j))+ &
                                                        BB*rr_sq(i,j)*bessel_j0(BB*rr_sq(i,j)))
                                        om_i(i,j)=2*BB*bessel_j1(BB*rr_sq(i,j))*sin(theta(i,j))/AA
                                else if(rr_sq(i,j) == 0.0) then
                                        u_0(i,j)=1/AA
                                        v_0(i,j)=0.0
                                else 
                                        u_0(i,j)=(1.0+1.0/rr(i,j))-2*X(i,j)*X(i,j)/(rr(i,j)*rr(i,j))
                                        v_0(i,j)=-2.0*X(i,j)*Y(i,j)/(rr(i,j)*rr(i,j))
                                        om_i(i,j)=0.0
                                end if
                        end do
                end do
        end subroutine init_vorticity              


end module solver
