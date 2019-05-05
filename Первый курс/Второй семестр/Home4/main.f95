program main
use sort
use quadra
implicit none

!--------------SORT--------------

integer n(5), a(5), b(2), rev(5), rev_r(size(rev)), cyc(6), ord(5), find(7), findex(3), find_z(size(find)), sort1(5), sort2(7)
integer length
integer j,k
character direction

!-----------QUADRA---------------

integer n_r, i
real(8) a_r, b_r, h_r
real(8), dimension(10000)::codomain
real(8), dimension(size(codomain))::domain

integer n_t
real(8) a_t, b_t, h_t
real(8), dimension(10000)::codomain_t
real(8), dimension(size(codomain))::domain_t

integer n_s
real(8) a_s, b_s, h_s
real(8), dimension(10000)::codomain_s
real(8), dimension(size(codomain))::domain_s

!--------------SORT--------------

write(*,*) ' '
write(*,*) ' SORT '
write(*,*) ' '

write(*,*) '[ Извиняюсь за лишний вывод после каждого выполнения функции, &
так и не удосужился его убрать.       ]'
write(*,*) '[ (особенно в третьем задании (инверсия вектора) у меня &
 много лишнего выходит, будьте внимательны!) ]'
write(*,*) ' '

read(*,*) n
write (*,*) 'Введён вектор n:', n

write(*,*) 'Количество нечетных чисел:', task1(n)
write(*,*) ' '; write(*,*) ' '

read(*,*) a
write(*,*) 'Введён вектор a: ', a
read(*,*) b
write(*,*) 'Введён вектор b: ', b

if (size(a) .eq. size(b)) then
	length=size(a)
elseif (size(a) .gt. size(b)) then
	length=size(a)
else
	length=size(b)
endif

write(*,*) 'Вектор с чередованиями из векторов a и b:'
write(*,*) 'Полученный с помощью функции:', task2(a, b, length)
call subtask2(a, b, length)
!write(*,*) 'Тест: '
!write(*,*) test(a, b, stat)
write(*,*) ' '; write(*,*) ' '

read(*,*) rev
write(*,*) 'Введен вектор rev: ', rev
write(*,*) 'Его инверсия, полученная нерекурсивным алгоритмом: ', task3(rev)
i=1
write(*,*) 'Его инверсия, полученная рекурсивным алгоритмом:', rectask3(rev,i,rev_r)
write(*,*) ' '; write(*,*) ' '

read(*,*) cyc
write(*,*) 'Введен вектор cyc: ', cyc

read(*,*) direction
if (direction .eq. 'r') then
write(*,*) 'Выбран сдвиг вправо.'
elseif (direction .eq. 'l') then 
write(*,*) 'Выбран сдвиг влево.'
else 
stop 'Неправильно введено направление сдвига.'
endif

write(*,*) 'Циклический сдвиг вектора cyc: ', task4(cyc, direction)
write(*,*) ' '; write(*,*) ' '

write(*,*) 'Проверка на упорядоченность (строгое возрастание).'
read(*,*) ord
write(*,*) 'Введен вектор ord: ', ord
write(*,*) task5(ord)
write(*,*) ' '; write(*,*) ' '

write(*,*) 'Поиск по образцу.'
read(*,*) find
write(*,*) 'Введен вектор для поиска: ', find
read(*,*) findex
write(*,*) 'Введен образец поиска', findex
write(*,*) ' '
write(*,*) 'Результаты поиска по нерекурсивному алгоритму: '
write(*,*) task6(find, findex)
write(*,*) ' '

!>

do i=1, size(find), 1
	find_z(i)=0
enddo

do i=1, size(find), 1
	if (findex(1) .eq. find(i)) then
	find_z(i)=i
	endif
enddo

write(*,*) 'Маска первого элемента образца (индексы): ', find_z

k=0

i=1
j=0

write(*,*) 'Результаты поиска по рекурсивному алгоритму: '
call rectask6(find, findex, find_z, i,j,k)

!<

write(*,*) ' '; write(*,*) ' '

read(*,*) sort1, sort2
write(*,*) 'Были введены два упорядоченных массива:'
write(*,*) sort1
write(*,*) sort2

write(*,*) task7(sort1,sort2)


!-------------QUADRA-------------

write(*,*) ' '
write(*,*) ' QUADRA '
write(*,*) ' '

write(*,*) ' Формула средних прямоугольников '
write(*,*) ' '

read(*,*) a_r, b_r, n_r
write(*,*) 'Был задан промежуток: от ', a_r, 'до ', b_r, '.'
write(*,*) 'Будут рассмотрены первые ', n_r, 'значений(-я).'

h_r=(b_r-a_r)/n_r

do i=1, size(codomain), 1
	domain(i)=0
	codomain(i)=0
enddo

do i=1, n_r, 1
	domain(i)=a_r+(i-0.5)*h_r
	codomain(i)=your_function(domain(i))
enddo

write(*,*) 'Таблица подинтегральной функции выведена в файл results.dat'
open (unit=101, file='results.dat')
do i=1, n_r, 1
	write(101, '(e15.7,e15.7)') domain(i), codomain(i)
enddo
close(101)

write(*,*) 'Rectan:'
write(*,*) rectan(codomain, a_r, b_r, n_r) 

!--------------------------------

write(*,*) ' '
write(*,*) ' Формула трапеций '
write(*,*) ' '

read(*,*) a_t, b_t, n_t
write(*,*) 'Был задан промежуток: от ', a_t, 'до ', b_t, '.'
write(*,*) 'Будут рассмотрены первые ', n_t+1, 'значений(-я).'

h_t=(b_t-a_t)/n_t

do i=1, size(codomain_t), 1
	domain_t(i)=0
	codomain_t(i)=0
enddo

do i=1, n_t+1, 1
	domain_t(i)=a_t+(i-1)*h_t
	codomain_t(i)=your_function(domain_t(i))
enddo

write(*,*) 'Таблица подинтегральной функции выведена в файл results_t.dat'
open (unit=102, file='results_t.dat')
do i=1, n_t+1, 1
	write(102, '(e15.7,e15.7)') domain_t(i), codomain_t(i)
enddo
close(102)

write(*,*) 'Trap:'
write(*,*) trap(codomain_t, a_t, b_t, n_t) 

!--------------------------------

write(*,*) ' '
write(*,*) ' Формула Симпсона (парабол) '
write(*,*) ' '

read(*,*) a_s, b_s, n_s
write(*,*) 'Был задан промежуток: от ', a_s, 'до ', b_s, '.'
write(*,*) 'Будут рассмотрены первые ', n_s+1, 'значений(-я).'

h_s=(b_s-a_s)/n_s

do i=1, size(codomain_s), 1
	domain_s(i)=0
	codomain_s(i)=0
enddo

do i=1, n_s+1, 1
	domain_s(i)=a_s+(i-1)*h_s
	codomain_s(i)=your_function(domain_s(i))
enddo

write(*,*) 'Таблица подинтегральной функции выведена в файл results_s.dat'
open (unit=103, file='results_s.dat')
do i=1, n_s+1, 1
	write(103, '(e15.7,e15.7)') domain_s(i), codomain_s(i)
enddo
close(103)

write(*,*) 'Sim:'
write(*,*) sim(codomain_s, a_s, b_s, n_s)

end
