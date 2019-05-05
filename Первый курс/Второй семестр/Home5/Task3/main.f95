program main
use functions
implicit none

integer b(8) / 1,2,3,4, 10,20,30,40 /
integer ord(2) / 2,1 / !изменять здесь!

!рассмотрим такой вариант: 8=4*2
!то есть s(1)=4, s(2)=2

integer shape1(2,4) !изменять здесь!

write(*,*) 'Исходный вектор: ', b

write(*,*) myresh2(ord,shape1,b)

end program main
