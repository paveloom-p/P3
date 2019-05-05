program main
use function
implicit none

character(100) a
character(100) b

character(3) ta, tb, tsumm

read(*,*) a, b

write(ta,'(i3)') 5+dim(len(trim(b)),len(trim(a)))
write(tb,'(i3)') 5+dim(len(trim(a)),len(trim(b)))
write(tsumm,'(i3)') max(len(trim(a)),len(trim(b)))+1

write(*,'(/,'//ta//'x,a100,/,5x,a,/,'//tb//'x,a100,/,5x,a,/,4x,a'//tsumm//',/)') a, ' + ', b, ' = ', summ(trim(a),trim(b))

end
