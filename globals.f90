! Module contains various global constants
! Eventually make N, Fh, Re input parameters from command line

module globals
        use netcdf
        use, intrinsic :: iso_c_binding
        implicit none
        ! N must be even
        integer :: N
        integer :: n_k, num_steps       
        integer :: cont=0
        real(kind=8) :: t_final=100.0_8

        real(kind=8), parameter :: pi=3.14159265358979323846264338328_8
!        real(kind=8), parameter :: Fh=0.2_8, Re=10000_8, Sc=1.0_8
        real(kind=8) :: Fh, Re, Sc=1.0_8
        real(kind=8), parameter :: L=9.0_8
        real(kind=8), parameter :: tpiL=2.0_8*pi/L
!        real(kind=8) :: dx=L/real(N,8), dy=dx
        real(kind=8) :: dx, dy
        real(kind=8), parameter :: dt=0.00075_8/2.0_8
        real(kind=8) :: t=0.0_8, en, prev_en, ukz
        complex(kind=8), parameter :: ii=(0.0_8,1.0_8)

        real(kind=8), dimension(:), allocatable :: growth_rate(:),tote(:)
        real(kind=8), dimension(:,:), allocatable :: X,Y
        complex(kind=8), dimension(:,:), allocatable :: kx,ky,kz,k_sq,kinv_sq,cut
        complex(kind=8), dimension(:,:), allocatable :: p11,p12,p13,p21,p22,p23,p31,p32,p33
        complex(kind=8), dimension(:,:), allocatable :: om_i,u_0,v_0
        complex(kind=8), dimension(:,:), allocatable :: uu_hat_temp,vv_hat_temp,ww_hat_temp
        complex(kind=8), dimension(:,:), allocatable :: iuu_hat_new,ivv_hat_new,iww_hat_new,irho_hat_new
        complex(kind=8), dimension(:,:), allocatable :: rr,ur,vr,wr
        complex(kind=8), dimension(:,:), allocatable :: rr_old,ur_old,vr_old,wr_old
        complex(C_DOUBLE_COMPLEX), pointer :: uu(:,:), vv(:,:), ww(:,:), rho(:,:)
        complex(C_DOUBLE_COMPLEX), pointer :: om1(:,:),om2(:,:),om3(:,:)
        complex(C_DOUBLE_COMPLEX), pointer :: r1(:,:),r2(:,:),r3(:,:),r4(:,:)
        complex(C_DOUBLE_COMPLEX), pointer :: A_1(:,:),B_1(:,:),C_1(:,:)
        complex(C_DOUBLE_COMPLEX), pointer :: uu_hat(:,:), vv_hat(:,:), ww_hat(:,:), rho_hat(:,:)
        complex(C_DOUBLE_COMPLEX), dimension(:,:), allocatable :: iuu_hat(:,:),ivv_hat(:,:), iww_hat(:,:), irho_hat(:,:)
        complex(C_DOUBLE_COMPLEX), pointer :: om1_hat(:,:),om2_hat(:,:),om3_hat(:,:)
        complex(C_DOUBLE_COMPLEX), pointer :: r1_hat(:,:), r2_hat(:,:),r3_hat(:,:),r4_hat(:,:)
        complex(C_DOUBLE_COMPLEX), pointer :: A_1_hat(:,:), B_1_hat(:,:),C_1_hat(:,:)
        type(C_PTR) :: up,vp,wp,rhop,uq,vq,wq,rhoq,funcp,funcq,om1hp,om2hp,om3hp,om1p,om2p,om3p
        type(C_PTR) :: rp1,rp2,rp3,rp4,rp1_hat,rp2_hat,rp3_hat,rp4_hat
        type(C_PTR) :: A1p,B1p,C1p,A1p_hat,B1p_hat,C1p_hat


end module globals
