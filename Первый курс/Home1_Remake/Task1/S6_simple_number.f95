subroutine simple_number

integer(8) i
integer, parameter :: jmax=1000000

integer(8) a_f, k_f
real(8) t21, t22, t23, t24
real(8) c


write(*,*) '6. Проверка на простое натуральное число.'

!a_f=7
read(*,*) a_f

write(*,*) 'Было введено натуральное число:', a_f
write(*,*) ' '

call cpu_time(t21)
do j=1,jmax
call simple_number_nonrec(a_f,k_f)
enddo

write(*,*) 'Результат нерекурсивного метода:'
if (k_f .eq. 1) then
write(*,*) 'Введенное число - простое.'
else
write(*,*) 'Введенное число не является простым.'
endif

call cpu_time(t22)

write(*,*) 'Время, затраченное на нерекурсивный алгоритм:', t22-t21
write(*,*) ' '

call cpu_time(t23)
do j=1, jmax

k_f=0
i=1
c=a_f

call simple_number_rec(a_f,k_f,i,c)

enddo

write(*,*) 'Результат рекурсивного метода:'
if (k_f .eq. 1) then
write(*,*) 'Введенное число - простое.'
else
write(*,*) 'Введенное число не является простым.'
endif
call cpu_time(t24)

write(*,*) 'Время, затраченное на рекурсивный алгоритм:  ', t24-t23
write(*,*) ' '


!=======
contains

subroutine simple_number_nonrec(a,k)
integer(8) a, i, k			! П: Перевел числа в тип integer(8).
					!    Не так много: 170001 ещё считает,
					!    180001 уже нет.
real(8) c

c=a

!write(*,*) c
!write(*,*) sqrt(c)
!write(*,*) int(sqrt(c))+1
!write(*,*) ' '

k=0
do i=1, int(sqrt(c))+1, 2
	if (mod(a,i) .eq. 0) then	! П: Теперь алгоритм проверяет ~(sqrt(n)+1)/2
	k=k+1				!    претендентов на делитель простого числа.
	endif
enddo

end subroutine simple_number_nonrec

!----------------------------------

recursive subroutine simple_number_rec(a,k,i,c)
integer(8) a, i, k
real(8) c

if (i .le. int(sqrt(c))+1) then
	if (mod(a,i) .eq. 0) then
	k=k+1
	endif

	i=i+1
	call simple_number_rec(a,k,i,c)
endif

end subroutine simple_number_rec

end subroutine simple_number
