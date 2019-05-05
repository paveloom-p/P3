program list_spectral ! Вывод списка звёзд заданного спектрального класса. Файл distrib.f95
use Hipmain
implicit none

! Расположение полной версии каталога Hipparcos
character(*), parameter :: HipparcosName ='/home/paveloom/Hipparcos/1/version_cd/cats/hip_main.dat'

character(10) Want_Sp

type(THipparcos) s ! Описание переменной s типа THipparcos.

read(*,*) Want_Sp
write(*,'(/,a,5x,A10)') ' Выбранный спектральный класс:', Want_Sp


write(*,'(a,/)') ' Далее будут выведены (если найдены) HIP звёзд:'
open(3, file = HipparcosName); open(4, file = 'result')

do while (ReadHipparcosMain(s)) ! Пока не прочли весь каталог
if (s%Sp .eq. Want_Sp) write(*,*) s%HIP
enddo

close(4); close(3)

write(*,*) ' '

end program list_spectral
