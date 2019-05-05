program main
use function
implicit none

character(50) b

read(*,*) b

write(*,'(/,a,7x,a)') 'Прямой вывод:', trim(b)
write(*,'(a,3x,a,/)') 'Обратный порядок:', reverse(trim(b))

end
