module subroutines
implicit none

!interface rdr
!module procedure rdr_1, rdr_1l
!end interface rdr

contains

! Подсчет значения
function your_function(n,x) result(y)
integer n
real(8) x, y

y=x**n

end function

! Ввод коээфициентов

subroutine rdr_1(p,n) ! Ввод в обычный массив
integer n
real(8) p(0:n)

read(*,*) p

end subroutine rdr_1

! Вывод коэффициентов

subroutine wrt_1(p,n)
integer n
real(8) p(0:n)

character(3) sk
character(60) sf

write(sk,'(i3)') n+1
sf='(" Массив коэффициентов:",'//sk//'e15.7)'

write(*,sf) p

end subroutine wrt_1

! Ввод первой строки

subroutine initial_1(d,a,n)
integer n
real(8) a, d(-1:n)

character(3) sk
character(60) sf

d(-1)=a; d(0)=your_function(n,a)
d(1:n)=(/15,50,60,24/)

write(*,*) ' '
write(*,'(17x,a,2x,a)') '============================', '------------------&
&----------------------------------------'

write(sk,'(i3)') n+2
sf='(" Первая строка:",'//sk//'e15.7)'

write(*,sf) d

end subroutine initial_1

! Две субрутины на вывод i-ой строки

subroutine babbage_search(d,n,i) ! II
integer n, i, j
real(8) d(-1:n) 
real(8) find(i), result_

find=d(0:i-1)
write(*,*) find

do while(i .gt. 1)

do j=1, i
find(j)=find(j+1)+find(j)
enddo
i=i-1
write(*,*) find(1:i)

enddo

result_=find(1)
write(*,'(/,a,e15.7)') 'Для третьего аргумента был получен результат', result_

end subroutine babbage_search

subroutine babbage(d,n) ! I
integer n
real(8) d(-1:n)
integer i

read(*,*) i
write(*,'(/,a,i2,a,/)') ' Будет произведен поиск', i, '-ой строки и &
&воспроизведена часть таблицы.'

call babbage_search(d,n,i)

end subroutine babbage


end module subroutines
