! Module contains various global constants
!
! Eventually make N, Fh, Re input parameters from command line
!
!
!
!

module globals
        implicit none
        ! N must be even
        integer, parameter :: N=8
        integer, parameter :: Nr=(N/2)+1
        integer :: n_k
        real, parameter :: pi=3.14159265
        real, parameter :: Fh=0.2, Re=10000, Sc=1
        real, parameter :: L=30.0
        real, parameter :: tpiL=2.0_8*pi/L
        real, parameter :: dx=L/real(N), dy=dx
        real, parameter :: dt=1e-2
        real(kind=8), dimension(:,:), allocatable :: u,v,w,rho
        real(kind=8), dimension(:,:), allocatable :: kx,ky,kz,k_sq,kinv_sq,X,Y,omega,cut
        real(kind=8), dimension(:,:), allocatable :: p11,p12,p13,p21,p22,p23,p31,p32,p33
!        real, dimension(:,:), allocatable :: kx,ky,kz,k_sq,kinv_sq,X,Y,omega
        complex, parameter :: ii=(0.0_8,1.0_8)
end module globals
