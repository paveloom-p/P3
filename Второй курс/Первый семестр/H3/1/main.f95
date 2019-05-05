program main
use function
implicit none

integer a
read(*,*) a

write(*,'(/,a,15x,i1)') 'Однозначное целое число:', a
write(*,'(a,5x,a,/)') 'Обозначающая его десятичная цифра:', int_to_digit(a)

end
