program main
use guarda
implicit none
integer, parameter :: ninp=5,nres=6
real(8), allocatable, dimension(:) :: x,y
real(8) a,b,h,sm
integer n,i
open(unit=nres, file='result',status='replace')
open(unit=ninp, file='input')
!write(*,*) 'Введите границы а и b'
read(ninp,*) a,b
!write(*,*) 'Введите количество значений функций'
read(ninp,*) n
allocate(x(n),y(n))
h=(b-a)/n
do i=2,n
x(i)=a+h*(i-1)
y(i)=f(x(i))
enddo
write(nres,'(a,8x,a,14x,a)') '#','x','y'
do i=2,n
    write(nres,'(1x,e15.7,e15.7)') x(i),y(i)
enddo
sm=trap(y,a,b,n)
write(nres,'(9x,a,17x,a)') '#a','b'
write(nres,'(a,e15.7,5x,e15.7)') '#',a,b
write(nres,'(a,i8)') '#n=',n
write(nres,'(a,e15.7)') '#S=',sm
deallocate(x,y);
close(nres)
close(ninp)
end