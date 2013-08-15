module display
        use globals
        implicit none
        
contains 

        !write real matrix to file
        subroutine mat_w2f(mat,fname,N)
                integer, intent(in) :: N
                real(kind=8), intent(in), dimension(N,N) :: mat 
                character (len=*), intent(in) :: fname
                integer :: i,j
                open(unit=1,file=fname)
                do i=1,N
                        do j=1,N
                        write(1,*), mat(i,j)
                        end do
                end do
                close(1)
        end subroutine mat_w2f  

        !write array to file
        subroutine arr_w2f(mat,fname,N)
                integer, intent(in) :: N
                real(kind=8), intent(in), dimension(N) :: mat 
                character (len=*), intent(in) :: fname
                integer :: i
                open(unit=1,file=fname)
                do i=1,N
                        write(1,*), mat(i)
                end do
                close(1)
        end subroutine arr_w2f  
        
        !write complex matrix to file
        subroutine mat_w2f_c(mat,fname,N)
                integer, intent(in) :: N
                complex(C_DOUBLE_COMPLEX), intent(in), dimension(N,N) :: mat 
                character (len=*), intent(in) :: fname
                integer :: i,j 
                open(unit=1,file=fname)
                do i=1,N
                        do j=1,N
                                write(1,*), real(mat(i,j))
                                write(1,*), imag(mat(i,j))
                        end do
                end do
                close(1)
        end subroutine mat_w2f_c

        !write some data to file
        subroutine data_w2f(mat,fname,time_step)
                integer, intent(in) :: time_step
                real(kind=8), intent(in) :: mat 
                character (len=*), intent(in) :: fname
                open(unit=1,file=fname)
                        write(1,*), 'Fh=', Fh
                        write(1,*), 'N=', N
                        write(1,*), 'Re=', Re
                        if (hypervis==1) then
                          write(1,*), 'ReV=',Rev
                        end if
                        write(1,*), 'kz=', ukz
                        write(1,*), 'dt=', dt
                        write(1,*), 'dealias=',dealias_coeff
                        write(1,*), 'L=' ,L
                        write(1,*), 'growth rate=', mat
                        write(1,*), 'current time step=', time_step
                        write(1,*), 'total time steps=', num_steps
                        write(1,*), 'current time=', time_step*dt
                close(1)
        end subroutine data_w2f  

        subroutine mat_r2f_c(mat,fname,N)
                integer, intent(in) :: N
                complex(C_DOUBLE_COMPLEX), intent(inout), dimension(N,N) :: mat 
                character (len=*), intent(in) :: fname
                real(kind=8) :: re,im
                integer :: i,j 
                open(unit=1,file=fname)
                do i=1,N
                        do j=1,N
                                read(1,*), re
                                read(1,*), im
                                mat(i,j)=cmplx(re,im,8)
                        end do
                end do
                close(1)
        end subroutine mat_r2f_c 
  
        subroutine isnan_matrix(inmat,hasnan,N)
                integer, intent(in) :: N 
                integer, intent(inout) :: hasnan 
                integer :: i,j
                complex(C_DOUBLE_COMPLEX), intent(inout), dimension(N,N) :: inmat  
                real :: remat,immat
                hasnan=0
                do i=1,N
                  do j=1,N
                    if(isnan(real(inmat(i,j)))) then
                        hasnan=1  
                    end if
                    if(isnan(imag(inmat(i,j)))) then
                        hasnan=1  
                    end if
                  end do
                end do
        end subroutine isnan_matrix        
end module display
