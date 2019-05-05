program main
use subprograms
implicit none

real(8) x, eps, res, umax, a, b, h, gamma
integer k, n, i

gamma=0.57721566490153286061d0

read(*,*) a, b, n, eps

h=(b-a)/n

write(*,'(/,2x,a,/)') 'Моделирование функции E_1(x) + gamma + ln(x)&
& по разложению в ряд Маклорена.'
write(*,'(2x,a,e15.7,/)') 'Выбранная относительная погрешность ограничения: ', eps

write(*,'(2x,a,15x,a,22x,a,19x,a,21x,a)') 'x', 'res', 'E_1(x)', 'umax', 'k' 

do i=0, n, 1

x=a+h*i
call E1_my(x,eps,res,umax,k)

write(*,'(e15.7,e25.16,e25.16,e25.16,2x,i3)') x, res, res-gamma-log(x), umax, k

enddo

write(*,*) ' '

end program main
