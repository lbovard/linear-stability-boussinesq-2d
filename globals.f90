! Module contains various global constants
!
! Eventually make N, Fh, Re input parameters from command line

module globals
        use, intrinsic :: iso_c_binding
        implicit none
        ! N must be even
        integer, parameter :: N=16
        integer, parameter :: Nr=(N/2)+1
        integer :: n_k
        real, parameter :: pi=3.14159265
        real, parameter :: Fh=0.2, Re=10000, Sc=1.0_8
        real, parameter :: L=30.0
        real, parameter :: tpiL=2.0_8*pi/L
        real, parameter :: dx=L/real(N), dy=dx
        real, parameter :: dt=1e-2
        real, parameter :: ukz=4.0
        real(C_DOUBLE), pointer :: uu(:,:), vv(:,:), ww(:,:), rho(:,:)
        real(kind=8), dimension(:,:), allocatable :: kx,ky,kz,k_sq,kinv_sq,X,Y,cut
        real(kind=8), dimension(:,:), allocatable :: p11,p12,p13,p21,p22,p23,p31,p32,p33
        real(kind=8), dimension(:,:), allocatable :: om1,om2,om3
        real(kind=8), dimension(:,:), allocatable :: om_i,u_0,v_0
        real(kind=8) :: t=0.0_8
        complex, parameter :: ii=(0.0_8,1.0_8)
        complex(C_DOUBLE_COMPLEX), pointer :: uu_hat(:,:), vv_hat(:,:), ww_hat(:,:), rho_hat(:,:)
        !integrating factor versions
        complex(C_DOUBLE_COMPLEX), dimension(:,:), allocatable :: iuu_hat(:,:),ivv_hat(:,:), iww_hat(:,:), irho_hat(:,:)
        !time step variables
        complex(C_DOUBLE_COMPLEX), dimension(:,:), allocatable :: uu_hat_new(:,:), vv_hat_new(:,:), ww_hat_new(:,:), &
                 rho_hat_new(:,:)
        type(C_PTR) :: up,vp,wp,rhop,uq,vq,wq,rhoq,funcp,funcq
end module globals
