subroutine less_than_1_binary

integer i, j
integer, parameter :: jmax=1000000

real(8) a_e, c_e
integer w_e(32)
real(8) t17, t18, t19, t20


write(*,*) '5. Перевод дробного неотрицательного числа &
в двоичную систему счисления.'

!a_e=0.125 (по умолчанию)
read (*,*) a_e

write(*,*) 'Было введено число:', a_e
write(*,*) ' '

call cpu_time(t17)
do j=1, jmax
call less_than_1_binary_nonrec(a_e,w_e)
enddo
write(*,*) 'Двоичное представление дробной части по нерекурсивному методу:'
write(*,'(" ",32i1)') w_e			! П: Сделал поприличнее вывод.
call cpu_time(t18)

write(*,*) 'Время, затраченное на нерекурсивный алгоритм:', t18-t17
write(*,*) ' '

call cpu_time(t19)
do j=1, jmax

do i=1, 32, 1
	w_e(i)=0
enddo

c_e=a_e
i=1

call less_than_1_binary_rec(a_e, c_e, w_e, i)

enddo

write(*,*) 'Двоичное представление дробной части по рекурсивному методу:'
write(*,'(" ",32i1)') w_e
call cpu_time(t20)

write(*,*) 'Время, затраченное на рекурсивный алгоритм:  ', t20-t19
write(*,*) ' '

!=======
contains

subroutine less_than_1_binary_nonrec(a,w)
real(8) a, c
integer b, w(32), i

do i=1, 32, 1
	w(i)=0
enddo

c=a
do i=1, 32, 1
	c=c*2
	b=int(c)
	
	w(i)=b			! П: Упростил заполнение массива w.
	c=c-b	

enddo

end subroutine less_than_1_binary_nonrec

!----------------------------------

recursive subroutine less_than_1_binary_rec(a,c,w,i)
real(8) a, c
integer b, w(32), i

if (i .le. 32) then
	c=c*2
	b=int(c)
	
	w(i)=b
	c=c-b
	
	i=i+1

	call less_than_1_binary_rec(a,c,w,i)
endif

end subroutine less_than_1_binary_rec

end subroutine less_than_1_binary
