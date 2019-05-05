module functions 
contains

!-----------------------------

function myresh2(shape1,b,s1,s2,ord) result(shape2)
integer b(:), s1, s2, i, j
integer, optional :: ord(2)
integer o
integer shape1(:,:), shape2(s1,s2)

!write(*,*) 's1=', s1
!write(*,*) 's2=', s2

!------------

o=0
if (present(ord)) o=1     ! П: Странно, но сработало только так. Хотел, конечно,
                          !    обойтись без лишних проверки и цикла, ну да ладно.

if (o .eq. 1) then
if (ord(1) .eq. 1 .and. ord(2) .eq. 2) then
        do j=1, s2, 1
                do i=1, s1, 1
                shape1(i,j)=b(i+(j-1)*s1)
                enddo
        enddo

elseif (ord(1) .eq. 2 .and. ord(2) .eq. 1) then

        do j=1, s1, 1
                do i=1, s2, 1
                shape1(j,i)=b(i+(j-1)*s2)
                enddo
        enddo
endif
else
        do j=1, s2, 1
                do i=1, s1, 1
                shape1(i,j)=b(i+(j-1)*s1)
                enddo
        enddo
endif

shape2=shape1

end function myresh2

!-----------------------------

end module
