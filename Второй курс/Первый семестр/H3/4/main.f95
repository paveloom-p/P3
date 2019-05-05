program main
use function
implicit none

integer a
read(*,*) a

write(*,'(/,a,t52,i10)') 'Исходное число:', a
write(*,'(a,2x,a,/)') 'Его шестнадцатеричное представление:', hex(a)

end
