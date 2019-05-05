program main
use guarda
implicit none
real(8) x(1000000),y(1000000),a,b,h,sm
integer n,i
write(*,*) 'Введите границы а и b'
read(*,*) a,b
write(*,*) 'Введите количество значений функций'
read(*,*) n
h=(b-a)/n
do i=1,n
x(i)=a+h*(i-0.5)
y(i)=f(x(i))
enddo
sm=rectan(y,a,b,n)
write(*,'(9x,a,17x,a)') 'a','b'
write(*,'(e15.7,5x,e15.7)') a,b
write(*,'(a,i8)') 'n=',n
write(*,*) 'S=',sm
end