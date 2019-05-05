program main
use function
implicit none

integer a
integer z

read(*,*) a
read(*,*) z

write(*,'(/,a,t65,i10)') 'Исходное число:', a
write(*,'(a,i2,a,t93,a/)') 'Его представление в ', z, '-ичной системе счисления:', z_hex(a,z)

end
