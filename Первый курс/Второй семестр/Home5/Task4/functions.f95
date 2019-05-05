module functions
contains

!-----------------------------

function myresh2(ord,shape1,b)
integer b(:), ord(3), s1, s2, s3, i, j, k, a
integer shape1(:,:,:), test(size(shape1(:,1,1)),size(shape1(1,:,1)),size(shape1(1,1,:)))

s1=size(shape1(:,1,1))
s2=size(shape1(1,:,1))
s3=size(shape1(1,1,:))

write(*,*) 's1=', s1
write(*,*) 's2=', s2
write(*,*) 's3=', s3

!------------

if (ord(1) .eq. 1 .and. ord(2) .eq. 2 .and. ord(3) .eq. 3) then
do k=1, s3, 1
	do j=1, s2, 1
		do i=1, s1, 1
		shape1(i,j,k)=b(i+(j-1)*s1+(k-1)*s1*s2)
		enddo
	enddo
enddo

	elseif (ord(1) .eq. 2 .and. ord(2) .eq. 1 .and. ord(3) .eq. 3) then

do k=1, s3, 1
	do j=1, s1, 1
		do i=1, s2, 1
		shape1(j,i,k)=b(i+(j-1)*s2+(k-1)*s1*s2)
		enddo
	enddo
enddo

	elseif (ord(1) .eq. 3 .and. ord(2) .eq. 2 .and. ord(3) .eq. 1) then

do k=1, s1, 1
	do j=1, s3, 1
		do i=1, s2, 1
		shape1(k,i,j)=b(j+(i-1)*s3+(k-1)*s3*s2)
		enddo
	enddo
enddo

	elseif (ord(1) .eq. 3 .and. ord(2) .eq. 1 .and. ord(3) .eq. 2) then

do k=1, s1, 1
	do j=1, s2, 1
		do i=1, s3, 1
		shape1(k,j,i)=b(i+(k-1)*s3+(j-1)*s1*s3)
		enddo
	enddo
enddo

	elseif (ord(1) .eq. 1 .and. ord(2) .eq. 3 .and. ord(3) .eq. 2) then

do k=1, s1, 1
	do j=1, s2, 1
		do i=1, s3, 1
		shape1(k,j,i)=b(k+(i-1)*s1+(j-1)*s1*s3)
		enddo
	enddo
enddo

	elseif (ord(1) .eq. 2 .and. ord(2) .eq. 3 .and. ord(3) .eq. 1) then

do k=1, s1, 1
	do j=1, s2, 1
		do i=1, s3, 1
		shape1(k,j,i)=b(j+(i-1)*s2+(k-1)*s2*s3)
		enddo
	enddo
enddo

endif

!------------

write(*,*) 'Результат функции myresh2:', shape1

!проверка

if (ord(1) .eq. 1 .and. ord(2) .eq. 2 .and. ord(3) .eq. 3) then
test=reshape(b, (/s1,s2,s3/), order=(/1,2,3/))
write(*,*) 'Результат функции reshape:', test

elseif (ord(1) .eq. 2 .and. ord(2) .eq. 1 .and. ord(3) .eq. 3) then
test=reshape(b, (/s1,s2,s3/), order=(/2,1,3/))
write(*,*) 'Результат функции reshape:', test

elseif (ord(1) .eq. 3 .and. ord(2) .eq. 2 .and. ord(3) .eq. 1) then
test=reshape(b, (/s1,s2,s3/), order=(/3,2,1/))
write(*,*) 'Результат функции reshape:', test

elseif (ord(1) .eq. 3 .and. ord(2) .eq. 1 .and. ord(3) .eq. 2) then
test=reshape(b, (/s1,s2,s3/), order=(/3,1,2/))
write(*,*) 'Результат функции reshape:', test

elseif (ord(1) .eq. 1 .and. ord(2) .eq. 3 .and. ord(3) .eq. 2) then
test=reshape(b, (/s1,s2,s3/), order=(/1,3,2/))
write(*,*) 'Результат функции reshape:', test

elseif (ord(1) .eq. 2 .and. ord(2) .eq. 3 .and. ord(3) .eq. 1) then
test=reshape(b, (/s1,s2,s3/), order=(/2,3,1/))
write(*,*) 'Результат функции reshape:', test

endif


end function myresh2

!-----------------------------

end module
