program main
use sort
implicit none
real, allocatable, dimension(:) :: a
integer res,i,n,res1
real o
res=0; res1=0;
write(*,*) 'Введите кол-во элементов a вектора'
read(*,*) n
allocate(a(n))
write(*,*) 'Введите значение вектора a'
read(*,*) (a(i),i=1,n)
write(*,*) 'Введите образец'
read(*,*) o
write(*,*) 'Вектор a'
write(*,*) a
write(*,'(a,e15.2)') 'o=',o
res=search(a,n,o)
call searchr(1,n,a,o,res1)
if (res.ne.0.and.res1.ne.0) then
write(*,'(a,i8)') 'Номер элемента массива',res
write(*,'(a,i8)') 'Recursive version',res1
else
write(*,*) res,res1
write(*,*) 'Нет такого элемента'
endif
deallocate(a)
end