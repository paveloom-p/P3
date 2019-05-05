!Листинг E 3.2. Подсчет звезд без данных о координатах,
!собственных движениях и параллаксах

program prog_1
use HipMain

integer :: NoCoord = 0 ! Счетчик звезд без точных координат
integer :: NoProp = 0 ! Счетчик звезд без собств. движений
integer :: NoPar = 0 ! Счетчик звезд без параллаксов

type(THipparcos) :: s
call OpenHipparcosMain


do while (ReadHipparcosMain(s))
! Сравнение логических переменных
if (s%NoRADE) NoCoord = NoCoord+1
if (s%Nopm) NoProp = NoProp+1
if (s%NoPlx) NoPar = NoPar+1
end do

call CloseHipparcosMain

write(*,'(/,a,t15,i6)') ' No coord', NoCoord
write(*,'(a,t15,i6)') ' No PM', NoProp
write(*,'(a,t15,i6,/)') ' No Par ', NoPar

end program
