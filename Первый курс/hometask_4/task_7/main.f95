program main
use sort
implicit none
real, allocatable, dimension(:) :: a,b,c
integer res,i,n,m
write(*,*) 'Введите кол-во элементов a вектора'
read(*,*) n
write(*,*) 'Введите кол-во элементов b вектора'
read(*,*) m
allocate(a(n),b(m),c(n+m))
write(*,*) 'Введите значение вектора a'
read(*,*) (a(i),i=1,n)
write(*,*) 'Введите значение вектора b'
read(*,*) (b(i),i=1,m)
write(*,*) 'Вектор a'
write(*,*) a
write(*,*) 'Вектор b'
write(*,*) b
c=unite2(a,b,n,m)
write(*,*) 'Вектор а+b'
write(*,*) c
deallocate(a,b,c)
end