program main
use subprograms
implicit none

real(8) x, eps, res_mac, res_as, E_1, a, b, h, gamma
integer k_mac, k_as, n, i, ier

gamma=0.57721566490153286061d0

read(*,*) a, b, n, eps

h=(b-a)/n

write(*,'(/,2x,a,/)') 'Моделирование функции E_1(x) + gamma + ln(x) по разложению & 
&в ряд Маклорена и по её асимптотическому разложению.'
write(*,'(2x,a,e15.7,/)') 'Выбранная относительная погрешность ограничения: ', eps

write(*,'(2x,a,15x,a,2x,a,2x,a,4x,a,16x,a,12x,a,2x,a)') 'x', 'res (as) (apx. to E_1)', &
& 'k (as)', 'ier', 'res (mac)', 'E_1(x) (mac)', 'k (mac)', 'methods aer'

do i=0, n, 1

x=a+h*i
call E1_my(x,eps,res_mac,k_mac)
call E1_as(x,eps,res_as,k_as,ier)

E_1=res_mac-gamma-log(x)
write(*,'(e15.7,e25.16,t43,i3,6x,i1,2x,e25.16,e25.16,3x,i3,2x,e25.16)') x, res_as, &
& k_as, ier, res_mac, E_1, k_mac, abs(res_as-E_1)

enddo

write(*,*) ' '

end program main
