program main
use myarray
implicit none

integer a(5), b(5)
!изменять матрицы здесь (my_matmul)
!если захочется вектор, положите первое измерение равным 1
integer :: c(4,3)=reshape((/1,2,3, 4,5,6, 7,8,9, 10,11,12/), shape=(/size(c(:,1)),size(c(1,:))/))
integer :: d(3,4)=reshape((/1,2,3,4, 5,6,7,8, 9,10,11,12/), shape=(/size(d(:,1)),size(d(1,:))/))
integer sc1, sc2, sd1, sd2
!тоже изменять здесь (my_maxval)
!параметр dim изменять в input
integer :: e(3,4)=reshape((/1,2,3, 4,-4,-3, -2,-1,5, 0,-5,2/), shape=(/size(e(:,1)),size(e(1,:))/))
integer se1, se2, dim, i, j
!опять же (my_transpose)
integer :: f(4,4)=reshape((/1,2,3,4, -8,2,5,6, -8,-8,3,7, -8,-8,-8,4/), shape=(/size(f(:,1)),size(f(1,:))/))
integer sf1, sf2
!(my_minloc)
integer :: g(3,4)=reshape((/1,2,3, 4,-20,-3, -2,-1,15, 0,-5,2/), shape=(/size(g(:,1)),size(g(1,:))/))
integer mask, dim2
integer sg1, sg2


!-----------------------------

read(*,*) a, b
write(*,*) 'Были введены массивы:'
write(*,*) a
write(*,*) b

write(*,*) my_dot_production(a,b)

!-----------------------------

sc1=size(c(:,1))
sc2=size(c(1,:))
sd1=size(d(:,1))
sd2=size(d(1,:))

write(*,*) ' '
write(*,*) 'Были введены две матрицы:'
write(*,*) ' '
write(*,*) 'Измерения c:', sc1, sc2
write(*,1000) c 
!write(*,*) 'c:',c
write(*,*) ' '
write(*,*) 'Измерения d:', sd1, sd2
write(*,1001) d
!write(*,*) 'd:',d

write(*,*) my_matmul(c,d,sc2,sd1, sc1,sd2)

!для второго задания (matmul)
1000 format (1x, 'Матрица c:',(/5(10x,4i15/)))
1001 format (1x, 'Матрица d:',(/5(10x,4i15/)))

!-----------------------------

sf1=size(f(:,1))
sf2=size(f(1,:))

write(*,*) ' '
write(*,*) 'Была введена матрица f:'
write(*,*) ' '
write(*,*) 'Измерения f:', sf1, sf2
write(*,1005) f

write(*,*) ' '
write(*,1006) my_transpose(f,sf1,sf2)

1005 format (1x, 'Матрица f:', (/5(10x,4i15/)))
1006 format (1x, (/5(10x,4i15/)))

!-----------------------------

se1=size(e(:,1))
se2=size(e(1,:))

write(*,*) ' '
write(*,*) 'Была введена матрица.'
write(*,*) ' '
write(*,*) 'Измерения e:', se1, se2
write(*,1002) ((e(i,j),j=1,se2),i=1,se1)
!write(*,*) 'e:',e

read(*,*) dim
write(*,*) 'Параметр dim:', dim

write(*,*) ' '
write(*,*) 'Функция my_maxval расчитана только на двумерные массивы &
(матрицы).'
write(*,*) 'Для задания одномерного массива (вектора) положите размерность &
строк равной 1.'
write(*,*) 'Возможно использование параметра dim. За примером использования &
маски обратитесь к моей функции my_minloc'
write(*,*) ' '

write(*,1003) my_maxval(e, se1, se2, dim)

1002 format (1x, 'Матрица e:',(/5(10x,4i15/)))
1003 format (1x, 'my_maxval:',(/5(10x,4i15/)))

!-----------------------------

sg1=size(g(:,1))
sg2=size(g(1,:))

write(*,*) ' '
write(*,*) 'Была введена матрица.'
write(*,*) ' '
write(*,*) 'Измерения g:', se1, se2
write(*,1007) ((g(i,j),j=1,sg2),i=1,sg1)

write(*,*) ' '
write(*,*) 'Функция my_minloc расчитана только на двумерные массивы &
(матрицы).'
write(*,*) 'Для задания одномерного массива (вектора) положите размерность &
строк равной 1.'
write(*,*) 'Возможно использование параметров dim и mask. Маска задается &
вручную в тексте программы.'
write(*,*) ' '

! 0 - не использовать маску
! 1 - использовать маску
read(*,*) mask
! 0 - не использовать параметр dim
! 1-2 - номер измерения (только матрицы)
read(*,*) dim2

write(*,*) ' '
write(*,*) my_minloc(g,sg1,sg2,mask,dim2)

1007 format (1x, 'Матрица g:',(/5(10x,4i15/)))

end program
