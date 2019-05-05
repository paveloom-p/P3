!Листинг E 3.1. Модуль чтения каталога Hipparcos

module HipMain
implicit none


! Расположение полной версии каталога Hipparcos
character(*), parameter :: HipparcosName ='/home/paveloom/Hipparcos/1/version_cd/cats/hip_main.dat'

integer, parameter :: HipNumOfStars = 118218 ! Число звезд
integer, parameter :: u = 10 ! Номер файла


type THipparcos
sequence
integer(4) :: HIP ! Номер звезды по Hipparocs

! Астрометрическая информация
real(8) :: RAdeg,DEdeg ! экваториальные координаты в градусах
real(8) :: Plx ! тригонометрический параллакс в mas
real(8) :: pmRa,pmDE ! собственные движения ma*cos(d) и md
character(1) :: AstroRef ! Флаг для кратных систем

! Фотометрическая информация
real(4) :: VMag ! Звездная вел. по шкале Джонсона
real(4) :: B_V ! Показатель цвета B-V по шкале Джонсона

! Ошибки соответствующих величин
real(8) :: sigma_RAdeg,sigma_DEdeg
real(8) :: sigma_Plx
real(8) :: sigma_pmRa,sigma_pmDE

character(10) Sp ! Развернутый спектральный класс
logical NoRaDe ! Нет данных о точных координатах
logical NoPlx ! Нет данных о параллаксе
logical Nopm ! Нет данных о собственных движениях
logical NoVMag ! Нет данных о звездной величине
logical NoB_V ! Нет данных о показателе цвета

end type THipparcos


contains



subroutine OpenHipparcosMain ! Открытие файла каталога
open(u, file = HipparcosName)
end subroutine OpenHipparcosMain


subroutine CloseHipparcosMain ! Закрытие файла каталога
CLOSE(u)
end subroutine CloseHipparcosMain


logical function ReadHipparcosMain(s) ! Чтение данных о звезде
type(THipparcos), intent(out) :: s
character(450) hs ! Запись строки каталога

integer ier

!if (eof(u)) then
!ReadHipparcosMain=.false.
!return
!else
!ReadHipparcosMain=.true.
!end if

read(u,'(A450)',iostat=ier) hs ! Чтение одной строки каталога
ReadHipparcosMain = (ier .eq. 0)

!if (ier .eq. -1) then; print *, 'достигнут маркер окончания файла'; return; endif
if (ier .gt.  0) then; print *, 'в данном обнаружена ошибка'; stop 222; endif

! Сбрасываем флаги событий
s%NoRaDe = .False.
s%NoPlx =.False.
s%Nopm = .False.
s%NoVMag = .False.
s%NoB_V = .False.

! Интерпретация с 12 байт, начиная с 3-го - это номер HIP
read(hs(3:14),*) s%hip

! Чтение координат: по 12 байт с 52 и с 65 позиции
! Функция TRIM удаляет из строки пробелы, а LEN возвращает длину
! строки, соответственно, если это 0, то в строке только пробелы

if (len(trim(hs(52:63))) .eq. 0) then
s%NoRaDe = .true.
s%RADeg = 0.0 ! на всякий случай записываем 0
else
read(hs(52:63),*) s%RAdeg
end if

if (len(trim(hs(65:76))) .eq. 0) then
s%NoRaDe = .true.
s%DEDeg = 0.0
else
read(hs(65:76),*) s%DEdeg
end if

! Чтение параллакса - 7 байт с 80-й позиции
if (len(trim(hs(80:86))) .eq. 0) then
s%NoPlx = .true.
s%Plx = 0.0
else
read(hs(80:86),*) s%Plx
end if

! Чтение собственных движений: по 8 байт с 88 и с 97 позиции

if (len(trim(hs(88:95))) .eq. 0) then
s%NoPM = .true.
s%pmRA = 0.0
else
read(hs(80:86),*) s%pmRA
end if

if (len(trim(hs(88:95))) .eq. 0) then
s%NoPM = .true.
s%pmDE = 0.0
else
read(hs(97:104),*) s%pmDE
end if


s%AstroRef=hs(78:78) ! Флаг кратной звезды


! Чтение зв.величины и показателя цвета B-V по шкале Джонсона
if (len(trim(hs(42:46))) .eq. 0) then
s%NoVMag = .true.
s%VMag = 0.0
else
read(hs(42:46),*) s%VMag
end if

if (len(trim(hs(246:251))) .eq. 0) then
s%NoB_V = .true.
s%B_V = 0.0
else
read(hs(246:251),*) s%B_V
end if


! Данные об ошибках всегда есть, если присутствуют сами величины

if (.not. s%NoRADE) then
read(hs(106:111),*) s%sigma_RAdeg
read(hs(113:118),*) s%sigma_DEdeg
end if

if (.not. s%NoPlx) then
read(hs(120:125),*) s%sigma_Plx
end if

if (.not. s%Nopm) then
read(hs(127:132),*) s%sigma_pmRA
read(hs(134:139),*) s%sigma_pmDE
end if


s%Sp=hs(436:445) ! Чтение данных о спектральном классе


end function ReadHipparcosMain
end module HipMain
