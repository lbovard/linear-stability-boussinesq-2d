! Module contains various global constants
!
! Eventually make N, Fh, Re input parameters from command line

module globals
        use, intrinsic :: iso_c_binding
        implicit none
        ! N must be even
        integer, parameter :: N=16
        integer, parameter :: Nr=(N/2)+1
        real(kind=8) :: t_final=40.0
        integer :: n_k
        real, parameter :: pi=3.14159265
        real, parameter :: Fh=0.2, Re=10000, Sc=1.0_8
        real, parameter :: L=9.0
        real, parameter :: tpiL=2.0_8*pi/L
        real, parameter :: dx=L/real(N), dy=dx
        real, parameter :: dt=1e-2
        real, parameter :: ukz=4.0
        
        integer :: num_steps
        real(C_DOUBLE), pointer :: uu(:,:), vv(:,:), ww(:,:), rho(:,:)
        real(C_DOUBLE), pointer :: om1(:,:),om2(:,:),om3(:,:)
        real(C_DOUBLE), pointer :: r1(:,:),r2(:,:),r3(:,:),r4(:,:)
        real(C_DOUBLE), pointer :: A_1(:,:),B_1(:,:),C_1(:,:)
        real(kind=8), dimension(:,:), allocatable :: kx,ky,kz,k_sq,kinv_sq,X,Y,cut
        real(kind=8), dimension(:,:), allocatable :: p11,p12,p13,p21,p22,p23,p31,p32,p33
        real(kind=8), dimension(:,:), allocatable :: om_i,u_0,v_0
        real(kind=8) :: t=0.0_8, en
        real(kind=8), dimension(:), allocatable :: growth_rate(:)
        complex, parameter :: ii=(0.0_8,1.0_8)
        complex(C_DOUBLE_COMPLEX), pointer :: uu_hat(:,:), vv_hat(:,:), ww_hat(:,:), rho_hat(:,:)
        !integrating factor versions
        complex(C_DOUBLE_COMPLEX), dimension(:,:), allocatable :: iuu_hat(:,:),ivv_hat(:,:), iww_hat(:,:), irho_hat(:,:)
        !time step variables
        complex(C_DOUBLE_COMPLEX), dimension(:,:), allocatable :: iuu_hat_new(:,:), ivv_hat_new(:,:), iww_hat_new(:,:), &
                 irho_hat_new(:,:)
        ! new variables
        complex(C_DOUBLE_COMPLEX), dimension(:,:), allocatable :: rr(:,:),ur(:,:),vr(:,:),wr(:,:)
        complex(C_DOUBLE_COMPLEX), dimension(:,:), allocatable :: rr_old(:,:),ur_old(:,:),vr_old(:,:),wr_old(:,:)
        complex(C_DOUBLE_COMPLEX), pointer :: om1_hat(:,:),om2_hat(:,:),om3_hat(:,:)
        complex(C_DOUBLE_COMPLEX), pointer :: r1_hat(:,:), r2_hat(:,:),r3_hat(:,:),r4_hat(:,:)
        complex(C_DOUBLE_COMPLEX), pointer :: A_1_hat(:,:), B_1_hat(:,:),C_1_hat(:,:)
        type(C_PTR) :: up,vp,wp,rhop,uq,vq,wq,rhoq,funcp,funcq,om1hp,om2hp,om3hp,om1p,om2p,om3p
        type(C_PTR) :: rp1,rp2,rp3,rp4,rp1_hat,rp2_hat,rp3_hat,rp4_hat
        type(C_PTR) :: A1p,B1p,C1p,A1p_hat,B1p_hat,C1p_hat
end module globals
