program main
use subroutines
implicit none

integer n
real(8) p(0:100)  ! nmax=100
real(8) d(-1:100) ! обращение будет к первым n элементам
real(8) a, b
integer m

! Ввод степени полинома
read(*,*) n
write(*,'(/,a,i3,/)') ' Была введена степень полинома:', n

! Ввод коэффициентов
call rdr_1(p,n)

! Вывод коэффициентов
call wrt_1(p,n)

! Задание промежутка с дискретизацией по m
read(*,*) a, b, m

write(*,'(a,e15.7,a,e15.7,a,i3,a)') ' Задан промежуток от  ', a, '  до', b, '  с &
& делением на', m, '  промежутков,'
write(*,'(25x,a,e15.7)') ' то есть с шагом', (b-a)/m

! Ввод первой строки
call initial_1(d,a,n)

! Поиск i-ой строки

call babbage(d,n)

end
