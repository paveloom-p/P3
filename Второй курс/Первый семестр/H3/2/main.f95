program main
use function
implicit none

character b
read(*,*) b

write(*,'(/,a,20x,a)') 'Десятичная цифра:', b
write(*,'(a,5x,i1,/)') 'Числовой эквивалент целого типа:', digit_to_int(b)

end
