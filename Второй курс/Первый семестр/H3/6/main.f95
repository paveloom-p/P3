program main
use function
implicit none

character(50) a
character(50) b

read(*,*) a, b

write(*,'(/,5x,a50,/,5x,a,/,5x,a50,/,5x,a,/,1x,a50,/)') a, ' + ', b, ' = ', summ(a,b)

end
