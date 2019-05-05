program main
use sort
implicit none
real, allocatable, dimension(:) :: a,b,c
integer k,i,n,m
write(*,*) 'Введите кол-во элементов a и b вектора'
read(*,*) n,m
k=(m+n)
allocate(a(n),b(m),c(k))
write(*,*) 'Введите значение вектора a'
read(*,*) (a(i),i=1,n)
write(*,*) 'Введите значение вектора b'
read(*,*) (b(i),i=1,m)
write(*,*) 'a='
write(*,'(e15.7)') a
write(*,*) 'b='
write(*,'(e15.7)') b
call unite(a,b,m,n,c)
write(*,'(i8,e15.7)') (i,c(i),i=1,n+m)
deallocate(a,b,c)
end