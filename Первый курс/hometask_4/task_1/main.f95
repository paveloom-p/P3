program main
use sort
implicit none
integer, allocatable, dimension(:) :: a
integer k,i,n
write(*,*) 'Введите кол-во элементов вектора'
read(*,*) n
allocate(a(n))
write(*,*) 'Введите значение вектора'
read(*,*) (a(i),i=1,n)
write(*,*) a
k=0
k=chet(a,n)
deallocate(a)
write(*,'(a,i8)') 'Кол-во нечетных=',k
end