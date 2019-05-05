program main
use sort
implicit none
real, allocatable, dimension(:) :: a,b
integer k,i,n
write(*,*) 'Введите кол-во элементов a вектора'
read(*,*) n
allocate(a(n),b(n))
write(*,*) 'Введите значение вектора a'
read(*,*) (a(i),i=1,n)
write(*,*) 'Вектор a'
write(*,*) a
write(*,*) 'Нерекурсивный вариант инверсии'
b=invert(a,n)
write(*,*) b
write(*,*) 'Рекурсивный вариант инверсии'
b=0; k=1;
b=invertr(a,k,n)
write(*,*) b
deallocate(a,b)
end