program main
use subroutines
implicit none

integer a_def, b_def, a,b, c,d
real(8) t1, t2, t3, t4

integer a2_def, b2_def, a2,b2, c2,d2
real(8) t5, t6, t7, t8

integer v, i, c_, k_
real(8) t9, t10, t11, t12
integer a_(32)

integer a_d, n_d, r_d
real(8) t13, t14, t15, t16

real(8) a_e, c_e
integer w_e(32)
real(8) t17, t18, t19, t20

integer a_f, k_f
real(8) t21, t22, t23, t24

!числа для целочисленного деления
!a_def=12 (по умолчанию)
!b_def=4 (по умолчанию)
read(*,*) a_def, b_def

write(*,*) ' '
write(*,*) '1. Целочисленное деление.'
write(*,*) 'Введены числа', a_def, '         и', b_def
write(*,*) ' '

if (a_def .ge. b_def) then
	c=1
	if (mod(a_def,b_def) .ne. 0) then
	stop 'Выберите целые числа для целочисленного деления.'
	endif
else
	c=2
	if (mod(b_def,a_def) .ne. 0) then
	stop 'Выберите целые числа для целочисленного деления.'
	endif
endif

a=a_def
b=b_def
call cpu_time(t1)
call division_nonrec(a,b,c)
call cpu_time(t2)
write(*,*) 'Время, затраченное на нерекурсивный алгоритм:', t2-t1
write(*,*) ' '

d=0
a=a_def
b=b_def
call cpu_time(t3)
call division_rec(a,b,c,d)
write(*,*) 'Результат деления рекурсивного алгоритма:  ', d
call cpu_time(t4)
write(*,*) 'Время, затраченное на рекурсивный алгоритм:  ', t4-t3

!----------------------------------------------------------------------------

!числа для нецелочисленного деления
!a2_def=11 (по умолчанию)
!b2_def=4 (по умолчанию)
read(*,*) a2_def, b2_def

write(*,*) ' '
write(*,*) '2. Деление с остатком.'
write(*,*) 'Введены числа', a2_def, '         и', b2_def
write(*,*) ' '

if (a2_def .ge. b2_def) then
	c2=1
else
	c2=2
endif

a2=a2_def
b2=b2_def
call cpu_time(t5)
call division_left_nonrec(a2,b2,c2)
call cpu_time(t6)
write(*,*) 'Время, затраченное на нерекурсивный алгоритм:', t6-t5
write(*,*) ' '

a2=a2_def
b2=b2_def
call cpu_time(t7)
d2=0
call division_left_rec(a2,b2,c2,d2)
	write(*,*) 'Результат деления рекурсивного алгоритма:  ', d2
	if (c2 .eq. 1) then
	write(*,*) 'Найденный остаток:', a2
	elseif (c2 .eq. 2) then
	write(*,*) 'Найденный остаток:', a2 
	endif 
	call cpu_time(t8)
write(*,*) 'Время, затраченное на рекурсивный алгоритм:  ', t8-t7

write(*,*) ' '

!----------------------------------------------------------------------------

write(*,*) '3. Перевод числа в двоичную систему счисления.'

!программа рассчитана на целое число
!v=4 (по умолчанию)
read (*,*) v

write(*,*) 'Было введено целое число:', v
write(*,*) ' '

call cpu_time(t9)
call binary_nonrec(v)
call cpu_time(t10)
write(*,*) 'Время, затраченное на нерекурсивный алгоритм:  ', t10-t9
write(*,*) ' '

call cpu_time(t11)
i=32
k_=abs(v)
a_(1)=0
call binary_rec(i,k_,a_)
	if(v .lt. 0) then
		a_(1)=1
	endif
	write(*,*) 'Двоичное представление по рекурсивному методу:'
	write(*,1000) a_
	1000 format (/5(32i15/))
call cpu_time(t12)
write(*,*) 'Время, затраченное на рекурсивный алгоритм:  ', t12-t11
write(*,*) ' '

!----------------------------------------------------------------------------

write(*,*) '4. Возведение целого числа в &
неотрицательную целочисленную степень.'

!a_d=2 (по умолчанию)
!n_d=7 (по умолчанию) (>0, целое)
read(*,*) a_d, n_d

write(*,*) 'Было введено основание', a_d, '          и степень', n_d
write(*,*)

call cpu_time(t13)
call degree_nonrec(a_d, n_d)
call cpu_time(t14)

write(*,*) 'Время, затраченное на нерекурсивный алгоритм:', t14-t13
write(*,*) ' '

call cpu_time(t15)
r_d=a_d
i=0
call degree_rec(a_d,n_d,i,r_d)
write(*,*) 'Результат возведения в степень по рекурсивному методу:  ', r_d
call cpu_time(t16)

write(*,*) 'Время, затраченное на рекурсивный алгоритм:  ', t16-t15
write(*,*) ' '

!----------------------------------------------------------------------------

write(*,*) '5. Перевод дробного неотрицательного числа &
в двоичную систему счисления.'

!a_e=0.125 (по умолчанию)
read (*,*) a_e

write(*,*) 'Было введено число:', a_e
write(*,*) ' '

call cpu_time(t17)
call less_then_1_binary_nonrec(a_e)
call cpu_time(t18)

write(*,*) 'Время, затраченное на нерекурсивный алгоритм:', t18-t17
write(*,*) ' '

call cpu_time(t19)
do i=1, 32, 1
	w_e(i)=0
enddo

c_e=a_e
i=1

call less_then_1_binary_rec(a_e, c_e, w_e, i)

write(*,*) 'Двоичное представление дробной части по рекурсивному методу:'
write(*,1001) w_e
1001 format (/5(32i15/))
call cpu_time(t20)

write(*,*) 'Время, затраченное на рекурсивный алгоритм:  ', t20-t19
write(*,*) ' '

!----------------------------------------------------------------------------

write(*,*) '6. Проверка на простое натуральное число.'

!a_f=7
read(*,*) a_f

write(*,*) 'Было введено натуральное число:', a_f
write(*,*) ' '

call cpu_time(t21)
call simple_number_nonrec(a_f)
call cpu_time(t22)

write(*,*) 'Время, затраченное на нерекурсивный алгоритм:', t22-t21
write(*,*) ' '

call cpu_time(t23)
k_f=0
i=2

call simple_number_rec(a_f,k_f,i)

write(*,*) 'Результат рекурсивного метода:'
if (k_f .eq. 1) then
write(*,*) 'Введенное число - простое.'
else
write(*,*) 'Введенное число не является простым.'
endif
call cpu_time(t24)

write(*,*) ' '

end
