program main
use sort
implicit none
real, allocatable, dimension(:) :: a
integer ier,i,n
write(*,*) 'Введите кол-во элементов a вектора'
read(*,*) n
allocate(a(n))
write(*,*) 'Введите значение вектора a'
read(*,*) (a(i),i=1,n)
write(*,*) 'Вектор a'
write(*,*) a
ier=qsort(a,n)
  select case(ier)
    case(0); write(*,*) 'Упорядоченность по возр'
    case(1); write(*,*) 'Упорядоченность по убыв'
    case(2); write(*,*) 'Упорядоченность не убыв'
    case(3); write(*,*) 'Упорядоченность не возр'
    case(4); write(*,*) 'Упорядоченности нет'
  end select
deallocate(a)
end