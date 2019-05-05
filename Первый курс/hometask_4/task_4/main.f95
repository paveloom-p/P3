program main
use sort
implicit none
real, allocatable, dimension(:) :: a
integer ier,i,n
character(5) c
write(*,*) 'Введите кол-во элементов a вектора'
read(*,*) n
allocate(a(n))
write(*,*) 'Введите значение вектора a'
read(*,*) (a(i),i=1,n)
write(*,*) 'Введите направление сдвига'
read(*,*) c
if (c.eq.'right') then
    ier=0;
                  else
    if (c.eq.'left ') then
       ier=1;
                      else
       ier=2;
       write(*,*) 'некорректный ввод'
                      endif
                  endif
write(*,*) 'Вектор a'
write(*,*) a
a=shift(a,n,ier)
write(*,*) a
deallocate(a)
end