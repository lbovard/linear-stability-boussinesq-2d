! contains some global parameters
module globals
        implicit none
        integer, parameter :: N=6
        real, parameter :: Fh=0.2, Re=10000, Sc=1
        real, parameter :: L=30.0
        real, parameter :: dx=L/real(N), dy=dx
        real, parameter :: dt=1e-2
        real, dimension(:,:), allocatable :: kx,ky,kz,k_sq,kinv_sq,X,Y,omega
end module globals
