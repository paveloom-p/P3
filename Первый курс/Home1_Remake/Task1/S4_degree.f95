subroutine degree
implicit none

integer i
integer, parameter :: imax=1000000

integer a, a_d, n, n_d, r_d
integer k,q,s
real(8) t13, t14, t15, t16

write(*,*) '4. Возведение целого числа в &
неотрицательную целочисленную степень.'

!a_d=2 (по умолчанию)
!n_d=7 (по умолчанию) (>0, целое)
read(*,*) a, n

write(*,*) 'Было введено основание', a, '          и степень', n
write(*,*)

a_d=a
n_d=n
call cpu_time(t13)
	r_d=1
	do i=1, imax
	call degree_nonrec(a_d, n_d, r_d)
	enddo
	write(*,*) 'Результат возведения в степень по нерекурсивному методу:', r_d
call cpu_time(t14)

write(*,*) 'Время, затраченное на нерекурсивный алгоритм:', t14-t13
write(*,*) ' '

a_d=a
n_d=n
call cpu_time(t15)
	r_d=1
	do i=1, imax
	call degree_rec(a_d,n_d,r_d)
	enddo
write(*,*) 'Результат возведения в степень по рекурсивному методу:  ', r_d
call cpu_time(t16)

write(*,*) 'Время, затраченное на рекурсивный алгоритм:  ', t16-t15
write(*,*) ' '

!=======
contains

subroutine degree_nonrec(a,n,r)
integer a, n
integer r

do while (n .ne. 0)
	do while (mod(n,2) .eq. 0)	! П: Переписал оба алгоритма с использованием
		a=a*a			!    метода половинного деления.
		n=n/2
	enddo
	r=r*a
	n=n-1
enddo

end subroutine degree_nonrec

!----------------------------------

recursive subroutine degree_rec(a,n,r)
integer a, n
integer r

if (n .gt. 0) then
	if (mod(n,2) .eq. 0) then
		a=a*a
		n=n/2
		call degree_rec(a,n,r)
	endif

	if (n .gt. 0) then		! П: Ещё один условный оператор, потому что
	r=r*a				!    без него программа зацикливала именно эту часть,
	n=n-1				!    несмотря на то, что условия, активирующие её,
	call degree_rec(a,n,r)		!    уже не соблюдались.
	endif
endif


end subroutine degree_rec

!-------------------------------------------------

end subroutine degree
