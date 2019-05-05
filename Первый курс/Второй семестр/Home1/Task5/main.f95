program main
use subroutines
implicit none

integer a_d, n_d, r_d, i
real(8) t13, t14, t15, t16, ts1, ts2

write(*,*) ' '
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
write(*,*) 'Результат возведения в степень по рекурсивному методу:  	  ', r_d
call cpu_time(t16)

write(*,*) 'Время, затраченное на рекурсивный алгоритм:  ', t16-t15
write(*,*) ' '

call cpu_time(ts1)
r_d=a_d
i=0
call degree_tail_rec(a_d,n_d,i,r_d)
write(*,*) 'Результат возведения в степень по методу &
 c хвостовой рекурсией:  ', r_d
call cpu_time(ts2)

write(*,*) 'Время, затраченное на рекурсивный алгоритм:  ', ts2-ts1
write(*,*) ' '

write(*,*) 'Действительно, не обнаружил разницу между моим &
рекурсивным алгоритмом и '
write(*,*) 'алгоритмом с хвостовой рекурсией, так как у меня, &
кажется, это и одно и'
write(*,*) 'то же.'
write(*,*) ' '

end
