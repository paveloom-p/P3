program main
use real_type
use actual
use my
implicit none

real(mp) x, eps, res_ac, res_my, a, b, h
integer k, n, i, j, ier
character(2) rtype
character(1) rn

include "wolfram.f" ! Точные результаты, полученные в Wolfram Alpha

read(*,*) eps ! Относительная погрешность ограничения

rtype=' ' ! Проверка, сколько символов займет запись числа mp
write(rtype,'(i2)') mp
if (rtype(1:1) .eq. ' ') then; rn='1'; else; rn='2'; endif

! Вывод заголовка
write(*,'(/,2x,a)') 'Моделирование функции &
& f(x) = (E_1(x) - E_1(x/2) + ln(2) - x/2 * (1 - 3 * x^2/8)) * 144/(7x^3)&
& по разложению в ряд Маклорена.'
write(*,'(2x,a,/)') 'Сравнение с разложением по полиномам Чебышева и& 
& приблизительно точными результатами, вычисленными с помощью Wolfram Alpha.'
write(*,'(2x,a,i'//rn//',a)') 'Выбранный тип вещественных чисел: REAL(', mp, ')'
write(*,'(2x,a,e15.7,/)') 'Выбранная относительная погрешность ограничения: ', eps
write(*,'(2x,a,/)') '(прим.: значение ier = 1 означает, что указанная точность&
& не может быть достигнута при данной разновидности REAL)'

do j=1, 3 ! Начало внешнего цикла (последнее значение - число рассматриваемых диапазонов)
read(*,*) a, b, n

h=(b-a)/n

write(*,'(2x,a,15x,a,22x,a,15x,a,22x,a,19x,a,2x,a,2x,a)') 'x', 'res', 'actual res', 'aer', 'rer (%)', 'k', 'ier', 'wolfram (≈)'

do i=0, n, 1 ! Начало внутреннего цикла

x=a+h*i

res_ac=actual_function(x)             ! Вычисление по формуле "в лоб"
call my_function(x,eps,res_my,k,ier)  ! Вычисление по разложению в ряд Маклорена

! Вывод заголовка таблицы для текущего диапазона
write(*,'(e15.7,e25.16,e25.16,e25.16,e25.16,2x,i3,3x,i1,e33.24)') x, res_my, res_ac, abs(res_my-res_ac), &
& abs(res_my-res_ac)/abs(res_ac)*100, k, ier, wolfram(i+j+(j-1)*10)

enddo ! Окончание внутреннего цикла

write(*,*) ' '

enddo ! Окончание внешнего цикла

end program main
