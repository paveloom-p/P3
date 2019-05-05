program main
use subprograms
implicit none

real(8) x, eps, res, umax, a, b, h, erf_x
integer k, n, i

read(*,*) a, b, n, eps

h=(b-a)/n

write(*,'(/,2x,a,/)') 'Моделирование функции erf(x) по разложению в ряд Маклорена.'
write(*,'(2x,a,e15.7,/)') 'Выбранная относительная погрешность ограничения: ', eps

write(*,'(2x,a,15x,a,22x,a,19x,a,22x,a,18x,a,21x,a)') 'x', 'res', 'erf(x)', 'aer', 'rer (%)', 'umax', 'k'

do i=0, n, 1

x=a+h*i
call erf_my(x,eps,res,umax,k)

erf_x=erf(x)
write(*,'(e15.7,e25.16,e25.16,e25.16,e25.16,e25.16,2x,i3)') x, res, erf_x, abs(res-erf_x), abs(res-erf_x)/erf_x*100, umax, k

enddo

write(*,*) ' '

end program main
