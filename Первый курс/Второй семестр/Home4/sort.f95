module sort
contains

!--------------------------------

function task1(v)
integer v(:), i, r
r=0
do i=1, size(v), 1
	if (mod(v(i), 2) .eq. 1) then 
	r=r+1
	endif
enddo
write(*,*) r

end function task1

!--------------------------------

function task2(v,w, length)
integer v(:), w(:), i, c
integer u(length*2)

c=0
	do i=1, size(u), 1
		if (mod(i, 2) .eq. 1) then
			if (i-c .le. size(v)) then
			u(i)=v(i-c)
			c=c+1
			else
			u(i)=0
			endif
		else 
			if (i-i/2 .le. size(w)) then
			u(i)=w(i-i/2)
			else
			u(i)=0
			endif
		endif
	enddo

write(*,*) u

end function task2

!function test(v,w,stat)
!integer v(:), w(:), stat, u(2,size(v))

!if (stat .eq. 1) then
!u=reshape(v,(/2,size(v)/),order=(/1,2/))
!endif

!write(*,*) u 

!end function test

subroutine subtask2(v, w, length)

integer v(:), w(:), i, c
integer u(length*2)

c=0
	do i=1, size(u), 1
		if (mod(i, 2) .eq. 1) then
			if (i-c .le. size(v)) then
			u(i)=v(i-c)
			c=c+1
			else
			u(i)=0
			endif
		else 
			if (i-i/2 .le. size(w)) then
			u(i)=w(i-i/2)
			else
			u(i)=0
			endif
		endif
	enddo

write(*,*) 'Полученный с помощью субрутины:', u

end subroutine subtask2

!--------------------------------

function task3(v)
integer v(:), r(size(v))

do i=1, size(v), 1
	r(i)=v((size(v)+1-i))
enddo
write(*,*) r

end function task3

!-----------------

recursive function rectask3(v,i,r) result(w)
integer v(:), r(size(v)), w(size(v))
integer i

if (i .le. size(v)) then
	r(i)=v((size(v)+1-i))
	
	i=i+1
	write(*,*) rectask3(v,i,r)
elseif (i .eq. size(v)+1) then
w=r
endif

end function rectask3

!--------------------------------

function task4(v, direction)
integer v(:), r(size(v)), ch
character direction

if (direction .eq. 'r') then
		ch=v(size(v))
		do i=size(v), 2, -1
		r(i)=v(i-1)
		enddo
		r(1)=ch
elseif (direction .eq. 'l') then
		ch=v(1)
		do i=1, size(v)-1, 1
		r(i)=v(i+1)
		enddo
		r(size(v))=ch
endif

write(*,*) r

end function task4

!--------------------------------

function task5(v) !упорядочен по возрастанию
integer v(:), i, stat

stat=0

do i=2, size(v), 1
if (v(i-1) .lt. v(i)) then
cycle
else
stat=1
exit
endif
enddo

if (stat .eq. 0) then
write(*,*) 'Цикл упорядочен.'
else
write(*,*) 'Цикл не упорядочен.'
endif

end function task5

!--------------------------------

function task6(v,w) !кажется, я случайно написал алгоритм поиска образца для любого массива
integer v(:), w(:), z(size(v)), i, j, k

do i=1, size(z), 1
	z(i)=0
enddo

do i=1, size(v), 1
	if (w(1) .eq. v(i)) then
	z(i)=i
	endif
enddo

write(*,*) 'Маска первого элемента образца (индексы): ', z

k=0
do i=1, size(z), 1
	if (z(i) .ne. 0) then
		do j=0, size(w)-1, 1
			if (v(i+j) .eq. w(j+1)) then
			k=k+1
			else
			k=0
			cycle
			endif
		enddo
		
		if (k .ne. 0) then 
		write(*,*) 'Значение счетчика совпадений k:', k
		endif		
		if (k .eq. size(w)) then
		write(*,*) 'Найден образец на промежутке с индексами от ', i, ' до ', i+size(w)-1
		k=0
		endif
	endif
enddo

!do i=1, size(w), 1
!v(k+i)=w(k+i)
!enddo

end function task6

!-----------------

recursive subroutine rectask6(v,w,z, i,j,k)
integer v(:), w(:), z(size(v)), i, j, k
integer us, uf

if (i .le. size(z)) then
	if (z(i) .ne. 0) then
		if (j .le. size(w)-1) then
			if (v(i+j) .eq. w(j+1)) then
			k=k+1
			j=j+1
			call rectask6(v,w,z, i,j,k)
			else
			k=0
			j=0
			i=i+1
			call rectask6(v,w,z, i,j,k)
			endif
		
		else

		write(*,*) 'Значение счетчика совпадений k:', k
		if (k .eq. size(w)) then
		us=i
		uf=i+size(w)-1
		write(*,*) 'Найден образец на промежутке с индексами от ', us, ' до ', uf
		k=0
		j=0
		endif
		endif
	endif
	
	i=i+1
	call rectask6(v,w,z, i,j,k)
endif

end subroutine rectask6

!--------------------------------

function task7(v,w) !упорядоченность по нестрогому возрастанию
integer v(:), w(:), z(size(v)+size(w))
integer i, j, k, c, sw

do i=1, size(v), 1
	z(i)=v(i)
enddo

do i=size(v)+1, size(z), 1
	z(i)=w(i-size(v))
enddo

write(*,*) 'Смешанный массив:', z

c=z(1)
do i=1, size(z), 1
	if (z(i) .lt. c) then
	c=z(i)
	endif
enddo

do i=1, size(z), 1
	if (z(i) .eq. c) then
	sw=z(1)
	z(1)=z(i)
	z(i)=sw
	endif
enddo

do i=2, size(z), 1
	k=z(i)
	z(1)=k
	j=i
		do while (k .lt. z(j-1))
		z(j)=z(j-1)
		j=j-1
		enddo
	z(j)=k
enddo
z(1)=c

write(*,*) 'Упорядоченный массив:', z

end function task7

!--------------------------------

end module
