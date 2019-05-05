! Подсчет числа всех звезд в каталоге.

program prog_2
use hipmain ! Подключение модуля HipMain
implicit none

integer(4) nstar   ! Счетчик числа звезд
type(THipparcos) s ! Описание переменной s типа THipparcos.

call OpenHipparcosMain

nstar=0 ! Число звезд до чтения

do while (ReadHipparcosMain(s)) ! Бесконечный цикл:
nstar=nstar+1 ! Учет в счётчике очередной звезы
enddo

call CloseHipparcosMain

write(*,'(/,1x,a,i6,a,i6,a,/)') 'Число звезд (nstar) = ', nstar, ' (действительно, ', HipNumOfStars, ')'
stop 0000

end program prog_2
