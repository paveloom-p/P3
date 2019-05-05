program main
use functions
implicit none

integer b(8) / 1,2,3,4, 10,20,30,40 /
integer i, j, s1, s2

! Рассмотрим такой вариант: 8=4*2,
! то есть s(1)=4, s(2)=2.

integer shape1(4,2) ! Изменять здесь.
integer shape_f(2)

integer A(size(shape1(:,1)),size(shape1(1,:))), AA(size(A(:,1)),size(A(1,:)))
integer ord(2) / 2,1 /  ! Изменять здесь.

s1=size(shape1(:,1))
s2=size(shape1(1,:))
shape_f(1)=s1; shape_f(2)=s2

write(*,*) 'Исходный вектор: ', b
write(*,*) 'Заданная форма:  ', shape_f
write(*,*) 'Заданный порядок:', ord

A=myresh2(ord,shape1,b,s1,s2)
AA=reshape(b,shape=shape_f,order=ord)

! П: В формате указывается количество столбцов,
!    стоит помнить при смене формы (указать s(2)).
!    (проверил на (4,2) и (2,4))

write(*,*) 'Матрица A:'
write(*,'(2i7)') ((A(i,j),j=1,s2),i=1,s1)
write(*,*) 'Матрица AA:'
write(*,'(2i7)') ((AA(i,j),j=1,s2),i=1,s1)


end program main
