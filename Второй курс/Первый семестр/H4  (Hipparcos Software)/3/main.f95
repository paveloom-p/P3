program coordinates ! Задача 3
use Hipmain
use subprograms
implicit none

character(*), parameter :: HipparcosName ='/home/paveloom/Hipparcos/1/version_cd/cats/hip_main.dat'

TYPE(THipparcos) s

integer i, j, m
integer(4) input_HIP(100)

integer order(100) ! "Маршрут" соединения точек
real(8) top_era ! Интересующая эпоха (текущая приравнивается равной нулю)
integer era_incr ! Дискретезация отрезка времени

do i=1, 100, 1 ! Считывание HIP исследуемых звёзд (i-1 - число звёзд)
read(*,*) input_HIP(i)
if (input_HIP(i) .eq. 0) exit
enddo

order=0
do m=1, 100, 1 ! Считывание "маршрута" (m-1 - число точек в маршруте)
read(*,*) order(m)
if (order(m) .eq. 0) exit
enddo

read(*,*) top_era, era_incr ! Считывание эры

! Открытие каталога Hipparcos (чтение), файла для вывода на дисплей, data-файла, файла для рассчёта
open(3, file=HipparcosName); open(4, file='result'); open(5, file='data_1'); open(6, file='prmotions')

write(4,'(/,7x,a,5x,a,21x,a,21x,a,21x,a)') 'HIP', 'RAdeg', 'DEdeg', 'pmRA', 'pmDE'

do while (ReadHipparcosMain(s)) ! (HIP-ы считываются в порядке возрастания)
        do j=1, i-1, 1          ! Вывести текущие координаты и собственные движения
        if (s%HIP .eq. input_HIP(j)) then
                write(4,*) s%HIP, s%RAdeg, s%DEdeg, s%pmRA, s%pmDE ! На дисплей
                write(5,*) "'", s%RAdeg, s%DEdeg, "'" ! В data-файл
                write(6,*) s%pmRA; write(6,*) s%pmDE; write(6,*) s%RAdeg; write(6,*) s%DEdeg ! В файл для рассчёта
        endif
        enddo
enddo

write(4,*) ' '

! Закрытие файлов (кроме result)
close(6); close(5); close(3)
! Переход к обработке полученных данных

write(4,'(7x,a,2x,i3)') 'Число указанных звёзд:', i-1
write(4,*) '      Выбранная эпоха:', top_era
write(4,*) '      Дискретизация:', era_incr

write(4,'(/,7x,a,/,7x,a,/)') 'Перерасчёт координат на выбранную эпоху:', &
& '(группа строчек -- звезда, строчка -- координаты в определенной эре)'
call calculation(i,top_era,era_incr) ! Перерасчёт координат на выбранную эпоху и их вывод

call plotting_data(i,order(1:m-1),era_incr) ! Подготовка данных к выводу

close(4)

end program coordinates

!Функции
