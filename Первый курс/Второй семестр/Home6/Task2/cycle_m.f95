module cycle_m
implicit none
contains

function bin_max(a)
integer a(:), b(size(a))
integer bin_max, i, j, c, k

b=a
c=b(1)
do j=2, size(a), 1
	c=c*2+b(j)
enddo
write(*,*) 'Число в десятичной системе:', c
k=c

c=0
do i=1, size(a)-1, 1
b=cshift(a, i)
write(*,*) 'Массив после сдвига:',b
	c=b(1)
	do j=2, size(a), 1
		c=c*2+b(j)
	enddo
	write(*,*) 'Число в десятичной системе:', c
	if (c .gt. k) then
		k=c
	endif
enddo

write(*,*) 'Максимальное число, полученное в результате сдвигов:', k

end function bin_max

end module cycle_m
