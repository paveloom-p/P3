program list_spectral ! Вывод списка звёзд заданного спектрального класса. Файл distrib.f95
use Hipmain
implicit none

character(10) Want_Sp

type(THipparcos) s ! Описание переменной s типа THipparcos.

read(*,*) Want_Sp
write(*,'(/,a,5x,A10)') ' Выбранный спектральный класс:', Want_Sp

!call get_out

call OpenHipparcosMain

do while (ReadHipparcosMain(s,Want_Sp)) ! Пока не прочли весь каталог
enddo

call CloseHipparcosMain

write(*,'(a,/)') ' Далее будут выведены (если найдены) HIP звёзд:'
call system('grep "[1-9]" output')
write(*,*) ' '

end program list_spectral
