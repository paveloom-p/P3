program main
use sort
use quadra
implicit none
integer ier
!--------------SORT-------------------------
integer k                            ! task1
integer, allocatable :: n(:)         ! task1
character( 3) sk                     ! task1
character(60) sf                     ! task1
! ------------------------------------------
!Оформляйте каждую задачу в отдельной директории!!!

! И ошибки искать легче, и не отвлекает мусор от других задач.
! И программа короче


integer, allocatable :: a(:)         ! task2
integer, allocatable :: b(:)         ! task2
integer, allocatable :: u(:)         ! task2


integer, allocatable :: rev(:)       ! task3
integer, allocatable :: rev_r(:)     ! task3
integer, allocatable :: w(:)         ! task3


integer, allocatable :: cyc(:)       ! task4
character direction                  ! task4


integer, allocatable :: ord(:)       ! task5


integer, allocatable :: find(:)      ! task6
!integer j                        
!integer findex(3), find_z(size(find))


integer, allocatable :: sort1(:)     ! task7
integer, allocatable :: sort2(:)     ! task7

!-----------QUADRA---------------

integer n_r, i                                   ! task8  (из пред. дом. задания)
real(8) a_r, b_r, h_r                            ! task8
real(8), allocatable :: codomain(:)              ! task8
real(8), allocatable :: domain(:)                ! task8

integer n_t                                      ! task9  (8 на листе заданий)
real(8) a_t, b_t, h_t                            ! task9
real(8), allocatable :: codomain_t(:)            ! task9
real(8), allocatable :: domain_t(:)              ! task9

integer n_s                                      ! task10 (9 на листе заданий)
real(8) a_s, b_s, h_s                            ! task10
real(8), allocatable :: codomain_s(:)            ! task10
real(8), allocatable :: domain_s(:)              ! task10

!--------------SORT--------------

! Задание 1

write(*,'(/a/)')  ' SORT ' ! Вывод тот же, а строк в исходном тексте одна

write(*,'(/,a,/)') ' Задание 1'

read(*,*) k; write(*,*) 'Число элементов:',k
allocate(n(k), stat=ier)

if (ier/=0) then; write(*,*) 'ier=',ier,' do not allocate array!'
                  stop 1
endif

read(*,*) n
!write(*,'(a,5i3)') 'Введён вектор n:', n

write(sk,'(i3)') k!; write(*,*) ' sk=',sk
sf='(" Введён вектор n:",'//sk//'i3)'!; write(*,*) ' sf=',trim(sf)     

! П: Но подождите, результат функции
!    task1(n) -- число. Ему и стандартного
!    вывода хватит. Давайте я использую этот
!    прием для вывода моего вектора, как вы и
!    указали в readme_4.

write(*,sf) n
write(*,'(a,i3,a,i7//)') ' Количество нечетных чисел task1(n)=', task1(n)
!write(*,sf) task1(n)
deallocate(n, stat=ier)   ! П: Буду высвобождать память после размещаемых массивов.
write(*,'(//)')           ! П: Буду отделять задания тремя строками в выводе.

!-----------------------------------------------------------------------------

! Задание 2

write(*,'(/,a,/)') ' Задание 2'

allocate(a(5), stat=ier)
if (ier/=0) then; write(*,*) 'ier=',ier,' do not allocate array!'
                  stop 1
endif

read(*,*) a
write(*,*) 'Введён вектор a: ', a


allocate(b(2), stat=ier)
if (ier/=0) then; write(*,*) 'ier=',ier,' do not allocate array!'
                  stop 1
endif

read(*,*) b
write(*,*) 'Введён вектор b: ', b

!if (size(a) .eq. size(b)) then      ! А не проще
!	length=size(a)               !            length=max(size(a),size(b))
!elseif (size(a) .gt. size(b)) then  ! Кроме того, предполагалось, что
!	length=size(a)               ! в векторе результата должны быть 
!else                                ! все элементы обоих векторов, т.е. 
!	length=size(b)               ! nc=na+nb Элементы вектора, не участвующеие
!endif                               ! в чередовании, просто дописываются в хвост
!write(*,*) '1) length=', length
!mlength=min(size(a),size(b))        ! П: Угу, исправил. Плюс, оказалось, что мне полезнее
                                     !    знать минимум, да и вычислять его не здесь, а
                                     !    внутри функции.
!write(*,*) '2) mlength=', mlength


write(*,'(/,a)') ' Вектор с чередованиями из векторов a и b:'
write(*,*) 'Полученный с помощью функции:  ', task2(a,b)   
                                     ! П: Жалко выделять и передавать размещаемый массив функции.

allocate (u(size(a)+size(b)), stat=ier)
        if (ier/=0) then; write(*,*) 'ier=',ier,' do not allocate array!'
                  stop 1
        endif
                                     ! П: А вот для субрутины выделить не жалко, да.
call subtask2(a,b,u)

write(*,*) 'Полученный с помощью субрутины:', u

deallocate(a, stat=ier)
deallocate(b, stat=ier)
deallocate(u, stat=ier)
write(*,'(//)')


!-----------------------------------------------------------------------------

! Задание 3

write(*,'(/,a,/)') ' Задание 3'

allocate (rev(5), stat=ier)
        if (ier/=0) then; write(*,*) 'ier=',ier,' do not allocate array!'
                  stop 1
        endif
allocate (rev_r(size(rev)), stat=ier)
        if (ier/=0) then; write(*,*) 'ier=',ier,' do not allocate array!'
                  stop 1
        endif
allocate (w(size(rev)), stat=ier)
        if (ier/=0) then; write(*,*) 'ier=',ier,' do not allocate array!'
                  stop 1
        endif

read(*,*) rev
write(*,*) 'Введен вектор rev: ', rev
write(*,*) 'Его инверсия, полученная нерекурсивным алгоритмом: ', task3(rev)

i=1
w=rectask3(rev,i,rev_r)
write(*,*) 'Его инверсия, полученная рекурсивным алгоритмом:   ', w

deallocate(rev, stat=ier)
deallocate(rev_r, stat=ier)
deallocate(w, stat=ier)
write(*,'(//)')


!-----------------------------------------------------------------------------

! Задание 4

write(*,'(/,a,/)') ' Задание 4'

allocate (cyc(6), stat=ier)
        if (ier/=0) then; write(*,*) 'ier=',ier,' do not allocate array!'
                  stop 1
        endif


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
deallocate(cyc, stat=ier)
write(*,'(//)')

!-----------------------------------------------------------------------------

! Задание 5

write(*,'(/,a,/)') ' Задание 5'

write(*,*) 'Проверка на упорядоченность (строгое возрастание).'

allocate(ord(5), stat=ier)
        if (ier/=0) then; write(*,*) 'ier=',ier,' do not allocate array!'
                  stop 1
        endif
 
read(*,*) ord
write(*,*) 'Введен вектор ord:  ', ord
write(*,*) '"0" - да, "1" - нет:', task5(ord)

deallocate(ord, stat=ier)
write(*,'(//)')

!-----------------------------------------------------------------------------

! Задание 6

write(*,'(/,a,/)') ' Задание 6'

write(*,*) 'Поиск по образцу.'

allocate(find(7), stat=ier)
        if (ier/=0) then; write(*,*) 'ier=',ier,' do not allocate array!'
                  stop 1
        endif

read(*,*) find
write(*,*) 'Введен вектор для поиска: ', find
!read(*,*) findex
!write(*,*) 'Введен образец поиска', findex
!write(*,*) ' '
read(*,*) k
write(*,*) 'Искомый элемент:', k
write(*,*) 'Результаты поиска по нерекурсивному алгоритму: ', task6(find,k)
!write(*,*) task6(find, findex)

!>

!do i=1, size(find), 1
!	find_z(i)=0
!enddo

!do i=1, size(find), 1
!	if (findex(1) .eq. find(i)) then
!	find_z(i)=i
!	endif
!enddo

!write(*,*) 'Маска первого элемента образца (индексы): ', find_z

!k=0

!i=1
!j=0

i=0
write(*,*) 'Результаты поиска по рекурсивному алгоритму:   ', rectask6(find,k,i)
!call rectask6(find, findex, find_z, i,j,k)

!<

deallocate(find, stat=ier)
write(*,'(//)')

!-----------------------------------------------------------------------------

! Задание 7

write(*,'(/,a,/)') ' Задание 7'

allocate(sort1(5), stat=ier)
        if (ier/=0) then; write(*,*) 'ier=',ier,' do not allocate array!'
                  stop 1
        endif

allocate(sort2(7), stat=ier)
        if (ier/=0) then; write(*,*) 'ier=',ier,' do not allocate array!'
                  stop 1
        endif

read(*,*) sort1, sort2
write(*,*) 'Были введены два упорядоченных массива:'
write(*,*) sort1
write(*,*) sort2
write(*,*) 'Результат объединения с упорядочиванием:'
write(*,*) task7(sort1,sort2)

deallocate(sort1, stat=ier)
deallocate(sort2, stat=ier)
write(*,'(//)')


!-------------QUADRA-------------

write(*,'(/,a,/)') ' QUADRA '

!-----------------------------------------------------------------------------

! П: Объяснил причину существования массива domain в письме.

! Задание 8

write(*,'(/,a,/)') ' Задание 8'

write(*,'(a,/)') ' Формула средних прямоугольников '

read(*,*) a_r, b_r, n_r
write(*,*) 'Был задан промежуток: от ', a_r, 'до ', b_r, '.'
write(*,*) 'Будут рассмотрены первые ', n_r, 'значений(-я).'

h_r=(b_r-a_r)/n_r

allocate(domain(10000), stat=ier)
        if (ier/=0) then; write(*,*) 'ier=',ier,' do not allocate array!'
                  stop 1
        endif

allocate(codomain(size(domain)), stat=ier)
        if (ier/=0) then; write(*,*) 'ier=',ier,' do not allocate array!'
                  stop 1
        endif

domain=0; codomain=0

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
write(*,*) 'Integral=~', rectan(codomain, a_r, b_r, n_r) 

deallocate(domain, stat=ier)
deallocate(codomain, stat=ier)
write(*,'(//)')
!
! Т.е. 
!
!
!
!
!
!
!
!
!


!-----------------------------------------------------------------------------

! Задание 9

write(*,'(/,a)') ' Задание 9'

write(*,'(/,a,/)') ' Формула трапеций '

read(*,*) a_t, b_t, n_t
write(*,*) 'Был задан промежуток: от ', a_t, 'до ', b_t, '.'
write(*,*) 'Будут рассмотрены первые ', n_t+1, 'значений(-я).'

h_t=(b_t-a_t)/n_t

allocate(domain_t(10000), stat=ier)
        if (ier/=0) then; write(*,*) 'ier=',ier,' do not allocate array!'
                  stop 1
        endif

allocate(codomain_t(size(domain_t)), stat=ier)
        if (ier/=0) then; write(*,*) 'ier=',ier,' do not allocate array!'
                  stop 1
        endif


domain_t=0; codomain_t=0

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
write(*,*) 'Integral=~', trap(codomain_t, a_t, b_t, n_t) 

deallocate(domain_t, stat=ier)
deallocate(codomain_t, stat=ier)
write(*,'(//)')

!!-----------------------------------------------------------------------------

! Задание 10

write(*,'(/,a)') ' Задание 10'

write(*,'(/,a,/)') ' Формула Симпсона (парабол) '

read(*,*) a_s, b_s, n_s
write(*,*) 'Был задан промежуток: от ', a_s, 'до ', b_s, '.'
write(*,*) 'Будут рассмотрены первые ', n_s+1, 'значений(-я).'

h_s=(b_s-a_s)/n_s

allocate(domain_s(10000), stat=ier)
        if (ier/=0) then; write(*,*) 'ier=',ier,' do not allocate array!'
                  stop 1
        endif

allocate(codomain_s(size(domain_s)), stat=ier)
        if (ier/=0) then; write(*,*) 'ier=',ier,' do not allocate array!'
                  stop 1
        endif


domain_s=0; codomain_s=0

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
write(*,*) 'Integral=~', sim(codomain_s, a_s, b_s, n_s)

deallocate(domain_s, stat=ier)
deallocate(codomain_s, stat=ier)
write(*,'(//)')

end
