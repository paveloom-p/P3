program distrib ! Вычисление распределения звезд Файл distrib.f95
use Hipmain ! Hipparcos по абсолютной величине.
implicit none

integer, parameter :: nplt=10
integer, parameter :: Mmin=-13
integer, parameter :: Mmax= 15

integer(4) nstar, nstar1 ! Счетчик числа звезд

TYPE(THipparcos) s ! Описание переменной s типа THipparcos.

integer ia(Mmin:Mmax) ! Статистика
integer round ! Тип функции округления
integer i ! Вспомогательная переменная
real(8) r ! Расстояние
real(8) M ! Абсолютная звёздная величина

ia=0 ! Обнуление статистики


call OpenHipparcosMain


nstar=0


do while (ReadHipparcosMain(s)) ! Пока не прочли весь каталог

if ( s%NoPlx ) cycle ! Нет данных о параллаксе
if ( s%Plx .le. 0d0) cycle ! Неположительный параллакс
if (s%sigma_Plx/s%Plx .gt. 0.5d0) cycle ! Точность параллакса хуже 50%

r=1000d0/s%Plx ! Расчёт расстояния в парсеках
M=s%Vmag-5d0*dlog10(r)+5d0 ! Расчёт абс.зв.величины
i=round(M)

if (i .ge. Mmin .and. i .le. Mmax) ia(i)=ia(i)+1

nstar=nstar+1

enddo


call CloseHipparcosMain


write(*,'(i3,3x,i7)') (i,ia(i),i=Mmin,Mmax)

nstar1=0

do i=Mmin, Mmax
nstar1=nstar1+ia(i)
enddo

print *, ' # nstar=',nstar
print *, ' # nstar1=', nstar1, sum(ia)

open (nplt,file='dist.plt',status='replace')
call plot(ia,Mmin, Mmax,nplt)
close(nplt)

stop 0000
end program distrib

!Функции

function round(x)
implicit none
integer round
real(8) x; if (x .lt. 0d0) then; round=idint(x-0.5d0)
else; round=idint(x+0.5d0)
endif
end function round

function str(nn)
implicit none
character(*) str
character(10) w
character(1) sf
integer nn, n, f

n=nn
str=''

do while (n/=0)
f=mod(n,10)
sf=achar(f+iachar('0'))
w=str
str=sf//w
n=n/10
enddo

end function str

!Субрутина

subroutine plot(ia,Mmin,Mmax,nplt)
implicit none
integer Mmin, Mmax
integer ia(Mmin:Mmax)
character(1) z; character(3) sx
character(5) str, sy, sp; character(128) txt(38)
integer s1, s2, ii, i, nplt

txt=' #'
txt(1)='set terminal postscript eps 32'
txt(2)='set output ''distrib.eps'' '
txt(3)='set nokey'
txt(4)='set border 15 lw 3'
txt(5)='set size 1.5, 1'
txt(6)='set xtics -8, 1, 16'

do i=Mmin, Mmax
z=' '

if (i .lt. 0) z='-'

ii=iabs(i)
s1=ii/10; s2=mod(ii,10)

if(s1 .ne. 0) then
sx=z //achar(s1+iachar('0')) //achar(s2+iachar('0'))
else
sx=z //achar(s2+iachar('0'))
endif

sy=str(ia(i))
sp=str(ia(i)+iabs(5000-ii*100))
txt(20+i)='set label "'//sy//'" at '//sx//' , '//sp//' center rotate by 90'

enddo

txt(36)='plot [-8:16] [0:30000] ''result'' u 1:2 w boxes lt -1 lw 10'
txt(37)='unset label'
txt(38)='reset'

do i=1,38
write(nplt,'(a)') txt(i)
enddo

end subroutine plot
