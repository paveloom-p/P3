module binrad
implicit none
contains

!-----------

function binrad5(a) !будем считать, что числа всегда положительны
integer a(:), i, k
integer binrad5

k=a(1)
do i=2, size(a), 1
	k=k*2+a(i)
enddo

write(*,*) 'Число в десятичной системе счисления:', k

end function binrad5

!-----------

function binrad6(a) !будем считать, что числа всегда положительны
integer a(:), i, k, c, b(size(a))
integer binrad6

k=a(1)
do i=2, size(a), 1
	k=k*2+a(i)
enddo

do i=1, size(b), 1
	c=mod(k,8)
	k=(k-c)/8
	b(size(b)+1-i)=c
enddo

write(*,*) 'Восьмеричная система счисления:', b

end function binrad6

!-----------

end module binrad
