program main
use subprograms
implicit none

real(8) x, eps, res, a, b, h, erfc_x
integer k, n, i, ier

read(*,*) a, b, n, eps

h=(b-a)/n

write(*,'(/,2x,a,/)') 'Моделирование функции erf(x) по разложению &
&в ряд Маклорена и по её асимптотическому разложению.'
write(*,'(2x,a,e15.7,/)') 'Выбранная относительная погрешность ограничения: ', eps

write(*,'(2x,a,15x,a,5x,a,3x,a,4x,a,19x,a,22x,a)') 'x', 'res (apx. to erfc)', 'ier', 'k', 'erfc(x)', 'aer', 'rer (%)'

do i=0, n, 1

x=a+h*i
call erfc_as(x,eps,res,k,ier)

erfc_x=erfc(x)
write(*,'(e15.7,e25.16,2x,i1,2x,i3,1x,e25.16,1x,e25.16,e25.16)') x, res, ier, &
&k, erfc_x, abs(res-erfc_x), abs(res-erfc_x)/erfc_x*100

enddo

write(*,*) ' '

end program main
