module functions
contains

!-----------------------------

function myresh2(ord,shape1,b)
integer b(:), ord(2), s1, s2, i, j, a
integer shape1(:,:), test(size(shape1(:,1)),size(shape1(1,:)))

s1=size(shape1(:,1))
s2=size(shape1(1,:))

write(*,*) 's1=', s1
write(*,*) 's2=', s2

!------------

if (ord(1) .eq. 2 .and. ord(2) .eq. 1) then

	do j=1, s1, 1
		do i=1, s2, 1
		shape1(j,i)=b(i+(j-1)*s2)
		enddo
	enddo

	else

	do j=1, s2, 1
		do i=1, s1, 1
		shape1(i,j)=b(i+(j-1)*s1)
		enddo
	enddo
	
endif
!------------

write(*,*) 'Результат функции myresh2:', shape1

!проверка

if (ord(1) .eq. 2 .and. ord(2) .eq. 1) then
test=reshape(b, (/s1,s2/), order=(/2,1/))
write(*,*) 'Результат функции reshape:', test
else
test=reshape(b, (/s1,s2/), order=(/1,2/))
write(*,*) 'Результат функции reshape:', test
endif

end function myresh2

!-----------------------------

end module
