!Листинг E 3.1. Модуль чтения каталога Hipparcos

module HipMain
implicit none

type THipparcos
sequence
integer(4) :: HIP ! Номер звезды по Hipparocs

! Астрометрическая информация
real(8) :: RAdeg,DEdeg ! экваториальные координаты в градусах
real(8) :: pmRa,pmDE ! собственные движения ma*cos(d) и md

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
read(hs(3:14),*) s%hip

! Чтение координат: по 12 байт с 52 и с 65 позиции
! Функция TRIM удаляет из строки пробелы, а LEN возвращает длину
! строки, соответственно, если это 0, то в строке только пробелы

if (len(trim(hs(52:63))) .ne. 0) then
read(hs(52:63),*) s%RAdeg
end if

if (len(trim(hs(65:76))) .ne. 0) then
read(hs(65:76),*) s%DEdeg
end if

! Чтение собственных движений: по 8 байт с 88 и с 97 позиции

if (len(trim(hs(88:95))) .ne. 0) then
read(hs(80:86),*) s%pmRA
end if

if (len(trim(hs(88:95))) .ne. 0) then
read(hs(97:104),*) s%pmDE
end if

end function ReadHipparcosMain
end module HipMain
