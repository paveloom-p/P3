module triangle
implicit none
contains

function tri(a, s1, s2)
real(8) a(:,:), switch(s1), nullifier(s2-1)
integer i, j, s1, s2, imp_i, imp_j, k, w, z, pr, t
real(8) tri
real(8) imp 

! Преобразование I - меняет строки местами
! Преобразование II - делит на скаляр
! Преобразование III - складывает строки с домножением на константу
! (в данном случае, с целью обнуления)

do k=0, 2, 1

if (k .gt. 0) then
	write(*,*) ' '
	write(*,*) 'Рассмотрим меньшую матрицу.'
endif

!---------------------------

!проверка на уже ступенчатую форму

w=0
do i=1, s1-k, 1
	do j=1, s2-k, 1
		if (j .gt. i .and. a(i+k,j+k) .eq. 0) then
		w=w+1
		endif
	enddo
enddo

pr=0
if (mod(s1-k,s2-k) .ge. 1 .and. s1-k .lt. s2-k) then
pr=s1-k
elseif (mod(s1-k,s2-k) .eq. 0) then
pr=s1-1-k
elseif(mod(s1-k,s2-k) .ge. 1 .and. s1-k .gt. s2-k) then
pr=s1-mod(s1-k,s2-k)-1-k
endif  

z=0

do i=1, pr, 1
	z=z+i
enddo

if (w .eq. z) then
stop 'Матрица уже имеет ступенчатый вид. &
Преобразования не необходимы.'
endif

!---------------------------

!поиск главного элемента

imp=0
imp_i=0
imp_j=0

do i=1, s1-k, 1
	if (imp .ne. 0) then
	exit
	endif
	do j=1, s2-k, 1
		if (a(i+k,j+k) .ne. 0) then
		imp=a(i+k,j+k)
		imp_i=i; imp_j=j
		exit
		endif
	enddo
enddo

if (imp .eq. 0) then
stop 'Главного элемента не нашлось. &
Матрица приведена к ступенчатому виду.'
endif

write(*,*) 'Главный элемент:', imp
write(*,*) 'Его индексы:', imp_i, imp_j

!---------------------------

if (imp_j .ne. 1) then
	do i=1, s1-k, 1
		switch(i)=a(i+k,imp_j+k)
	enddo

	do i=1, s1-k, 1
		a(i+k,imp_j+k)=a(i+k,1+k)
	enddo

	do i=1, s1-k, 1
		a(i+k,1+k)=switch(i)
	enddo

	imp_j=1

write(*,*) ''; write(*,*) 'Преобразование I: '; write(*, 1001) a

endif

!---------------------------

if (imp .ne. 1) then
do i=1, s1-k, 1
	a(i+k,1+k)=a(i+k,1+k)/imp
enddo

write(*,*) ''; write(*,*) 'Преобразование II: '; write(*, 1001) a
endif

!---------------------------

do i=1, size(nullifier), 1
	nullifier(i)=0
enddo

do j=2, s2-k, 1
	nullifier(j-1)=-a(1+k,j+k)
enddo

t=0
do i=1, size(nullifier), 1
	if (nullifier(i) .eq. 0) then
	t=t+1
	endif
enddo

if (t .ne. size(nullifier)) then

do j=2, s2-k, 1
	do i=1, s1-k, 1
		a(i+k,j+k)=a(i+k,j+k)+nullifier(j-1)*a(i+k,1+k)
	enddo
enddo

write(*,*) ''; write(*,*) 'Преобразование III: '; write(*, 1001) a

endif

!---------------------------

enddo

1001 format (/5(10x,3e15.7/))

end function tri

end module triangle
