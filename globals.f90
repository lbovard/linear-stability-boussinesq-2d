! Module contains various global constants
!
! Eventually make N, Fh, Re input parameters from command line
!
!
!
!

module globals
        implicit none
        integer, parameter :: N=6
        integer, parameter :: Nr=(N/2)+1
        real, parameter :: pi=3.14159265
        real, parameter :: Fh=0.2, Re=10000, Sc=1
        real, parameter :: L=30.0
        real, parameter :: dx=L/real(N), dy=dx
        real, parameter :: dt=1e-2
        real, dimension(:,:), allocatable :: kx,ky,kz,k_sq,kinv_sq,X,Y,omega
end module globals
