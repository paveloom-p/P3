!Листинг E 3.1. Модуль чтения каталога Hipparcos

module HipMain
implicit none





type THipparcos
sequence
integer(4) :: HIP ! Номер звезды по Hipparocs
character(10) Sp  ! Развернутый спектральный класс

end type THipparcos


contains


logical function ReadHipparcosMain(s) ! Чтение данных о звезде
type(THipparcos), intent(out) :: s
character(450) hs ! Запись строки каталога

integer ier

read(3,'(A450)',iostat=ier) hs ! Чтение одной строки каталога
ReadHipparcosMain = (ier .eq. 0)

!if (ier .eq. -1) then; print *, 'достигнут маркер окончания файла'; return; endif
if (ier .gt.  0) then; print *, 'в данном обнаружена ошибка'; stop 222; endif

! Интерпретация с 12 байт, начиная с 3-го - это номер HIP
read(hs(3:14),*) s%HIP

s%Sp=hs(436:445) ! Чтение данных о спектральном классе

end function ReadHipparcosMain

end module HipMain
