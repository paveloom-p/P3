module sort
contains

!--------------------------------

function task1(v)     result(r)     !  Функция должна возвращать результат
integer v(:), i, r                  !  через своё имя. Каким именно манером
r=0                                 ! --- Ваше дело: хотите так 
do i=1, size(v), 1                  ! но можно и по старинке 
	if (mod(v(i), 2) .eq. 1) then 
	r=r+1
	endif
enddo
                                    ! <--= тут дописав task1=r, или даже
                                    ! вообще, используя вместо r имя
                                    ! task1 в теле цикла do. Однако,
                                    ! наиболее оптимален приведённый вариант
                                    ! с result(r).
end function task1

!-----------------------------------------------------------------------------








!-----------------------------------------------------------------------------

function task2(v,w) result(u)                        ! П: Исправил вывод.
integer v(:), w(:)!, i, c
integer u(size(v)+size(w))                           ! П: Насколько понял, это не этично, да
integer mlength, length                              !    и алгоритмически не очень разумно (реализуемо)
                                                     !    выделять динамические массивы внутри функции
                                                     !    (а выделять в главной и лишний раз передавать сюда не хочу),
                                                     !    поэтому оставлю здесь и далее для функций автоматические массивы.

!c=0                                                 ! Почему так сложно?
                                                     ! см. so.f90 
                                                     ! ==================
                                                     ! У нас же тема операции
                                                     ! над массивами:
                                                     ! индексный триплет и т.д.

mlength=min(size(v),size(w))
                                                     ! П: Исправил. Так получше будет.
u(1:mlength:2)=v(1:mlength)
u(2:mlength:2)=w(1:mlength)

if (size(v) .gt. size(w)) then
u(mlength*2+1:size(u))=v(mlength+1:size(v))
else
u(mlength*2+1:size(u))=w(mlength+1:size(w))
endif


!write(*,*) u

                                                     
!       do i=1, size(u), 1
!              if (mod(i, 2) .eq. 1) then
!                     if (i-c .le. size(v)) then
!                        u(i)=v(i-c)
!                        c=c+1
!                        else
!                        u(i)=0
!                        endif
!              else 
!                     if (i-i/2 .le. size(w)) then
!                        u(i)=w(i-i/2)
!                        else
!                       u(i)=0
!                        endif
!                     endif
!       enddo

end function task2



subroutine subtask2(v,w,u)

integer v(:), w(:)!, i, c
integer u(size(v)+size(w))                       
integer mlength, length                         



mlength=min(size(v),size(w))
                    
u(1:mlength:2)=v(1:mlength)
u(2:mlength:2)=w(1:mlength)

if (size(v) .gt. size(w)) then
u(mlength*2+1:size(u))=v(mlength+1:size(v))
else
u(mlength*2+1:size(u))=w(mlength+1:size(w))
endif

end subroutine subtask2

!-----------------------------------------------------------------------------








!-----------------------------------------------------------------------------

function task3(v) result(r)           ! Вообще-то имелось ввиду, что
integer v(:), r(size(v))              ! все перестановки элементов
                                      ! должны происходить не выходя за 
                                      ! пределы одного единственного
r(1:size(v))=v(size(v):1:-1)          ! вектора, подаваемого на вход. 
                                      ! Он же и на выход
!write(*,*) r                         ! Так что выгодно оформлять:
                                      ! subroutine invert(v)
                                      !                         или 
                                      ! subroutine invert(v,n)

                                      ! П: В задании указана разработка функций --
                                      !    пускай они и будут (ведь там и рекурсия
                                      !    в придачу). А для вывода мне
                                      !    понадобится дополнительный вектор,
                                      !    что поделать, ничего страшного. А вот
                                      !    индексные триплеты облегчают жизнь, да.
end function task3

!-----------------

! П: Понятно, что субрутина лучше. Но что поделать.

!subroutine task3(v)
!integer v(:)
!r(1:size(v))=r(size(v):1:-1)
!end subroutine task3

recursive function rectask3(v,i,r) result(w)    ! П: Исправил вывод.
integer v(:), r(:), w(size(r))
integer i

if (i .le. size(v)) then                        ! П: Здесь придется оставить итерацию,
        r(i)=v((size(v)+1-i))                   !    а то в триплеты её как-то не вставить.
	
        i=i+1
        w=rectask3(v,i,r)
else
w=r
endif

end function rectask3

!-----------------------------------------------------------------------------








!-----------------------------------------------------------------------------

function task4(v, direction) result(r)     ! Используйте индексные триплеты!!!!
integer v(:), r(size(v)), ch               ! И опять все подвижки только
character direction                        ! в пределах вектора поданного на вход 
                                           ! Можно использовать одну простую рабочую  
!if (direction .eq. 'r') then              ! переменную типа integer.
!              ch=v(size(v))
!	       do i=size(v), 2, -1         ! Никаких циклов !!!
!              r(i)=v(i-1)
!              enddo
!              r(1)=ch                     ! П: Переделал в триплетах.
!elseif (direction .eq. 'l') then
!              ch=v(1)
!              do i=1, size(v)-1, 1
!              r(i)=v(i+1)
!              enddo
!              r(size(v))=ch
!endif

if (direction .eq. 'r') then 
        ch=v(size(v))
        r(2:size(v))=v(1:size(v)-1)
        r(1)=ch
elseif (direction .eq. 'l') then
        ch=v(1)
        r(1:size(v)-1)=v(2:size(v))
        r(size(v))=ch
endif

!write(*,*) r                           ! subroutine вход и выход через v

                                        ! П: Опять же, вся домашняя работа была
                                        !    на разработку функций, судя по описанию
                                        !    заданий.
                                        
                                        !    Ну да, исчез бы формальный аргумент, да
                                        !    изменилась бы буквально пара символов.

!        ch=v(size(v))
!        v(2:size(v))=v(1:size(v)-1)
!        v(1)=ch

end function task4

!-----------------------------------------------------------------------------








!-----------------------------------------------------------------------------

function task5(v) result(stat) !упорядочен по возрастанию
integer v(:), i, stat

stat=0

do i=2, size(v)   !, 1     зачем, лишнее засорят обзор
if (v(i-1) .lt. v(i)) then

cycle                               ! Проще! Зачем две альтернативы?
else                                !======= 
stat=1                              ! if (v(i-1) .gt. v(i)) then
exit                                   !         stat=1; exit
                                    ! endif
endif                               ! stat --- это результат, возвращаемый
enddo                               !          функцией, с тем, чтобы
                                    ! где-то в вызывающей её программе
                                    ! последняя могла выбрать ветвь
                                    ! дальнейшего хода своего алгоритма:
                                    !  if (task5(v)==0) then; ....
                                    !                   else; ....
                                    !  endif
!if (stat .eq. 0) then              !
!write(*,*) 'Цикл упорядочен.'      ! Вообще-то для вектора типа integer,
                                    ! естественно, результат сделать типа 
!else                               !  logical. Да, уж ладно! 
!write(*,*) 'Цикл не упорядочен.'
!endif                              ! А вот, если тип элементов real, то
                                    ! не всё так просто. См.

                                    ! П: Зато красиво было.



end function task5


!-----------------------------------------------------------------------------








!-----------------------------------------------------------------------------

!function task6(v,w) !кажется, я случайно написал алгоритм поиска образца для 
                     !любого массива
function task6(v,k) result(i)
! =============================================================================
! Не усдожняйте задачу! Запутаться можно и в гораздо более простой.
!
! На вход подаётся один вектор типа integer (причём обязательно строго
! упорядоченный) и образец в виде простой переменной типа integer
! Функция же должна вернуть либо 0, если образец не найден, либо индекс
! элемента равного образцу.

!  Никакими дополнительными массивами внутри  алгоритма пользоваться нельзя 
!  (максимум 2-3 простые переменные целого типа, если не одна)
! =============================================================================

! П: А жаль, хороший алгоритм был, универсальный! Ладно, переделал.

!integer v(:), w(:), z(size(v)), i, j, k

!do i=1, size(z), 1
!	z(i)=0
!enddo

!do i=1, size(v), 1
!	if (w(1) .eq. v(i)) then
!	z(i)=i
!	endif
!enddo

!write(*,*) 'Маска первого элемента образца (индексы): ', z

!k=0
!do i=1, size(z), 1
!	if (z(i) .ne. 0) then
!		do j=0, size(w)-1, 1
!			if (v(i+j) .eq. w(j+1)) then
!			k=k+1
!			else
!			k=0
!			cycle
!			endif
!		enddo
!		
!		if (k .ne. 0) then 
!		write(*,*) 'Значение счетчика совпадений k:', k
!		endif		
!		if (k .eq. size(w)) then
!		write(*,*) 'Найден образец на промежутке с индексами от ', i, ' до ', i+size(w)-1
!		k=0
!		endif
!	endif
!enddo

!!do i=1, size(w), 1
!!v(k+i)=w(k+i)
!!enddo

integer v(1:), k, i

i=0
do i=1, size(v)
        if (v(i) .eq. k) then
        exit
	endif
enddo

! П: Вот как-то грустно даже стало.

end function task6

!-----------------

!recursive subroutine rectask6(v,w,z, i,j,k)
recursive function rectask6(v,k,i) result(w)
!integer v(:), w(:), z(size(v)), i, j, k
!integer us, uf

!if (i .le. size(z)) then
!	if (z(i) .ne. 0) then
!		if (j .le. size(w)-1) then
!			if (v(i+j) .eq. w(j+1)) then
!			k=k+1
!			j=j+1
!			call rectask6(v,w,z, i,j,k)
!			else
!			k=0
!			j=0
!			i=i+1
!			call rectask6(v,w,z, i,j,k)
!			endif
!		
!		else
!
!		write(*,*) 'Значение счетчика совпадений k:', k
!		if (k .eq. size(w)) then
!		us=i
!		uf=i+size(w)-1
!		write(*,*) 'Найден образец на промежутке с индексами от ', us, ' до ', uf
!		k=0
!		j=0
!		endif
!		endif
!	endif
!	
!	i=i+1
!	call rectask6(v,w,z, i,j,k)
!endif

integer v(:), k, i, w

if (i .le. size(v)) then
        if (v(i) .eq. k) then
	w=i
	else
        i=i+1
        w=rectask6(v,k,i)
        endif
endif


end function rectask6

!--------------------------------

function task7(a,b)  result(c)    !упорядоченность по нестрогому возрастанию
integer a(:), b(:), c(size(a)+size(b))
!integer i, j, k, c, sw
                                  ! Нужно, имея на входе два упорядоченных
!do i=1, size(v), 1               ! массива так сразу вставлять их элементы в  
!	z(i)=v(i)                 ! в результирующий, чтобы  получился
!enddo                            ! тоже упорядоченный массив. 
                                  !
                                  ! Не нужно упорядочивать результирующий массив
                                  ! Он сразу должен компоноваться упорядоченным.
                                  ! 
                                  ! na=size(a); nb=size(b); nc=na+nb
                                  ! Пусть i --- текущий индекс массива a
                                  !       j --- текущий индекс массива b
                                  !       k --- текущий индекс компонуемого массива с
                                  !  Сначала: i=1  j=1  k=1
                                  !  Далее  : Пока (i<=na) и (j<=nb)
                                  !         :    Если a(i)< b(i), то с(k)=a(i); i=i+1
!do i=size(v)+1, size(z), 1       !         :                  иначе c(k)=b(j); j=j+1
!       z(i)=w(i-size(v))         !         :    илсЕ
!enddo                            !         :    k=k+1
                                  !         : акоП
                                  !  Таким образом, как только будет исчерпан
                                  !  наиболее короткий  массив, то в хвост массива С
                                  !  останется дописатьоставшиеся элементы более
                                  !  длинного
                                  !  Пока i<=na дописываем в хвост С остатки A
                                  !  Пока j<=nb дописываем в хвос  С остатки B
                                  !  И всё!
                                  
                                  !  П: Исправил.

integer na, nb, nc
integer i,j,k
na=size(a); nb=size(b); nc=na+nb

i=1; j=1; k=1

do while (i .le. na .and. j .le. nb)
        if (a(i) .lt. b(j)) then
                c(k)=a(i)
                i=i+1
        else
                c(k)=b(j)
                j=j+1
        endif
        k=k+1
enddo

do while (i .le. na)
        c(k)=a(i); i=i+1; k=k+1
enddo

do while (j .le. nb)
        c(k)=b(j); j=j+1; k=k+1
enddo


!write(*,*) 'Смешанный массив:', z

!c=z(1)
!do i=1, size(z), 1
!       if (z(i) .lt. c) then
!       c=z(i)
!       endif
!enddo

!do i=1, size(z), 1
!       if (z(i) .eq. c) then
!       sw=z(1)
!       z(1)=z(i)
!       z(i)=sw
!       endif
!enddo

!do i=2, size(z), 1
!       k=z(i)
!       z(1)=k
!       j=i
!              do while (k .lt. z(j-1))
!              z(j)=z(j-1)
!              j=j-1
!              enddo
!       z(j)=k
!enddo
!z(1)=c

!write(*,*) 'Упорядоченный массив:', z

end function task7

!--------------------------------

end module
