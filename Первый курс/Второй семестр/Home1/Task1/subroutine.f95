module subroutines
implicit none
contains

!----------------------------------------------------------------------------
subroutine division_nonrec(a,b,c) !целочисленное деление
integer a,b,c
integer d

! > нерекурсивный метод

d=0
if (c .eq. 1) then
	do while (a .ne. 0)
	d=d+1
	a=a-b
	end do
elseif (c .eq. 2) then
	do while (b .ne. 0)
	d=d+1
	b=b-a
	end do
endif
write(*,*) 'Результат деления нерекурсивного алгоритма:', d

! < нерекурсивный метод
end subroutine division_nonrec

!----------------------------------

recursive subroutine division_rec(a,b,c,d) !целочисленное деление
integer a,b,c,d

! > рекурсивный метод
if (c .eq. 1) then
	if (a .ne. 0) then
		d=d+1
		a=a-b
		call division_rec(a,b,c,d)
	endif
elseif (c .eq. 2) then
	if (b .ne. 0) then
		d=d+1
		b=b-a
		call division_rec(a,b,c,d)
	endif
endif
! < рекурсивный метод

end subroutine division_rec

!----------------------------------------------------------------------------

subroutine division_left_nonrec(a,b,c) !поиск остатка

integer a,b,c,d

! > нерекурсивный метод

d=0
if (c .eq. 1) then
	do while (a .gt. b )
	d=d+1
	a=a-b
	end do
elseif (c .eq. 2) then
	do while (b .gt. a)
	d=d+1
	b=b-a
	end do
endif
write(*,*) 'Результат деления нерекурсивного алгоритма:', d
if (c .eq. 1) then
write(*,*) 'Найденный остаток:', a
elseif (c .eq. 2) then
write(*,*) 'Найденный остаток:', a 
endif

! < нерекурсивный метод

end subroutine division_left_nonrec

!----------------------------------

recursive subroutine division_left_rec(a,b,c,d) !поиск остатка
integer a,b,c
integer d

! > рекурсивный метод

if (c .eq. 1) then
	if (a .gt. b ) then
	d=d+1
	a=a-b
	call division_left_rec(a,b,c,d)
	endif
elseif (c .eq. 2) then
	if (b .gt. a)  then
	d=d+1
	b=b-a
	call division_left_rec(a,b,c,d)
	endif
endif

! < рекурсивный метод

end subroutine division_left_rec

!----------------------------------------------------------------------------

subroutine binary_nonrec(v)
integer v, k, c, i
integer a(32)

k=abs(v)

do i=32, 1, -1
	c=mod(k,2)
	k=(k-c)/2
	a(i)=c
enddo

if(v .lt. 0) then
	a(1)=1
endif

write(*,*) 'Двоичное представление по нерекурсивному методу:'
write(*,1000) a
1000 format (/5(32i15/))

end subroutine binary_nonrec

!----------------------------------

recursive subroutine binary_rec(i,k,a)
integer k, c, i
integer a(32)

if (i .ge. 1) then
	c=mod(k,2)
	k=(k-c)/2
	a(i)=c
	i=i-1
	call binary_rec(i,k,a)
endif

end subroutine binary_rec

!----------------------------------------------------------------------------

subroutine degree_nonrec(a,n)
integer a, n
integer r, i

if (n .ge. 1) then

r=a
do i=1, n-1, 1
	r=r*a
enddo

elseif (n .eq. 0) then
r=1
endif

write(*,*) 'Результат возведения в степень по нерекурсивному методу:', r

end subroutine degree_nonrec

!----------------------------------

recursive subroutine degree_rec(a,n,i,r)
integer a, n
integer r, i

if (n .ge. 1) then

if (i .lt. n-1) then
	r=r*a
	i=i+1
	call degree_rec(a,n,i,r)
endif

elseif (n .eq. 0) then
r=1
endif

end subroutine degree_rec

!----------------------------------------------------------------------------

subroutine less_then_1_binary_nonrec(a)
real(8) a, c
integer b, w(32), i

do i=1, 32, 1
	w(i)=0
enddo

c=a
do i=1, 32, 1
	c=c*2
	b=int(c)
	
	if (mod(b,2) .eq. 0) then
	w(i)=0
	elseif (mod(b,2) .eq. 1) then
	w(i)=1
	endif
enddo

write(*,*) 'Двоичное представление дробной части по нерекурсивному методу:'
write(*,1000) w
1000 format (/5(32i15/))

end subroutine less_then_1_binary_nonrec

!----------------------------------

recursive subroutine less_then_1_binary_rec(a,c,w,i)
real(8) a, c
integer b, w(32), i

if (i .le. 32) then
	c=c*2
	b=int(c)
	
	if (mod(b,2) .eq. 0) then
	w(i)=0
	elseif (mod(b,2) .eq. 1) then
	w(i)=1
	endif
	
	i=i+1

	call less_then_1_binary_rec(a,c,w,i)
endif

end subroutine less_then_1_binary_rec

!----------------------------------------------------------------------------

subroutine simple_number_nonrec(a)
integer a, i, k

k=0
do i=2, a, 1
	if (mod(a,i) .eq. 0) then
	k=k+1
	endif
enddo

write(*,*) 'Результат нерекурсивного метода:'
if (k .eq. 1) then
write(*,*) 'Введенное число - простое.'
else
write(*,*) 'Введенное число не является простым.'
endif

end subroutine simple_number_nonrec

!----------------------------------

recursive subroutine simple_number_rec(a,k,i)
integer a, i, k

if (i .le. a) then
	if (mod(a,i) .eq. 0) then
	k=k+1
	endif
	
	i=i+1	

	call simple_number_rec(a,k,i)
endif

end subroutine simple_number_rec

end module subroutines
