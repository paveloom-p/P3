module myarray
implicit none
contains

function my_dot_production(a, b)
integer a(:), b(:), test, i, sc, pr
integer my_dot_production
real(8) t1, t2, t3, t4

write(*,*) ' '
write(*,*) 'Скалярное произведение.'
write(*,*) ' '

call cpu_time(t1)

sc=0
do i=1, size(a), 1
	sc=sc+a(i)*b(i)
enddo

call cpu_time(t2)

write(*,*) 'Результат функции my_dot_production:', sc

call cpu_time(t3)
test=dot_product(a,b)
call cpu_time(t4)
write(*,*) 'Результат функции dot_product:      ', test

write(*,*) ' '
write(*,1001) t2-t1
write(*,1002) t4-t3

1001 format(1x,'Время работы my_dot_production:', e15.7)
1002 format(1x,'Время работы dot_product:', e15.7)

end function my_dot_production

!-----------------------------------

function my_matmul(a,b, s1, s2, pr1, pr2)
integer a(:,:), b(:,:), res(pr1,pr2), res_m_v(pr1), res_v_m(pr2)
integer s1, s2, pr1, pr2, i, j, k, l, m, summ
integer my_matmul
real(8) t1, t2, t3

integer vector(pr2), vector_b(s1)

if (pr1 .gt. 1 .and. s2 .gt. 1) then

call cpu_time(t1)

if (s1 .ne. s2) then
stop 'Матрицы не согласованы.'
endif

write(*,*) ' '
write(*,*) 'Перемножение согласованных матриц.'
write(*,*) ' '

do i=1, pr1, 1
	do j=1, pr2, 1
		res(i,j)=0
		k=i
		l=j
			do m=1, s1, 1
				res(i,j)=res(i,j)+a(k,m)*b(m,l)
			enddo

	enddo
enddo

!write(*,*) res
write(*,1003) res

call cpu_time(t2)

!write(*,*) matmul(a,b)
write(*,1002) matmul(a,b)

call cpu_time(t3)

1002 format (1x, 'Результат matmul:', (/5(10x,3i15/)))
1003 format (1x, 'Результат my_matmul:', (/5(10x,3i15/)))

elseif (pr1 .gt. 1 .and. s2 .eq. 1) then

call cpu_time(t1)

if (s1 .ne. pr2) then
stop 'Перемножение матрицы на вектор невозможно, &
	так как различны количества столбцов.'
endif

do i=1, pr2, 1
	vector(i)=b(1,i)
enddo 

write(*,*) ' '
write(*,*) 'В виде вектора:', vector

write(*,*) ' '
write(*,*) 'Перемножение матрицы на вектор.'
write(*,*) ' '


do i=1, pr1, 1
	do j=1, pr2, 1
		res_m_v(i)=0
		k=i
			do m=1, s1, 1
				res_m_v(i)=res_m_v(i)+a(k,m)*b(m,1)
			enddo

	enddo
enddo

!write(*,*) res_m_v
write(*,1003) res_m_v

call cpu_time(t2)

!write(*,*) matmul(a,b)
write(*,1002) matmul(a,vector)

call cpu_time(t3)

elseif (pr1 .eq. 1 .and. s2 .gt. 1) then

call cpu_time(t1)

if (s1 .ne. s2) then
stop 'Перемножение вектора на матрицу невозможно, &
	так как различны количества столбцов в векторе и строк в матрице.'
endif

do i=1, s1, 1
	vector_b(i)=a(1,i)
enddo 

write(*,*) ' '
write(*,*) 'Вектор c в виде вектора:', vector_b

write(*,*) ' '
write(*,*) 'Перемножение вектора на матрицу.'
write(*,*) ' '

do i=1, pr1, 1
	do j=1, pr2, 1
		res_v_m(j)=0
		l=j
			do m=1, s1, 1
				res_v_m(j)=res_v_m(j)+a(1,m)*b(m,l)
			enddo

	enddo
enddo

!write(*,*) res_v_m
write(*,1003) res_v_m

call cpu_time(t2)

!write(*,*) matmul(vector_b,b)
write(*,1002) matmul(vector_b,b)

endif

write(*,*) ' '
write(*,1004) t2-t1
write(*,1005) t3-t2

1004 format(1x,'Время работы my_matmul:', e15.7)
1005 format(1x,'Время работы matmul:', e15.7)

end function my_matmul

!-----------------------------------

function my_transpose(a, sf1, sf2)
integer a(:,:), b(sf1, sf2), sf1, sf2
integer changer, i, j, k
integer my_transpose

b=a

k=0
do i=1, sf1-1, 1
	do j=2+k, sf2, 1
		changer=b(i,j)
		b(i,j)=b(j,i)
		b(j,i)=changer
	enddo
	k=k+1
enddo

write(*,1005) b
write(*,1006) transpose(a)

1005 format (1x, 'Результат my_transpose:',(/5(10x,4i15/)))
1006 format (1x, 'Результат transpose:',(/5(10x,4i15/)))

end function my_transpose

!-----------------------------------

function my_maxval(a,se1,se2,dim)
integer a(:,:), maxarr1(se1), maxarr2(se2)
integer dims(size(shape(a)))
integer i, j, se1, se2, maxv, dim, act
integer my_maxval

if (dim .eq. 0) then

maxv=-9999999
do i=1, se1, 1
	do j=1, se2, 1
		if (a(i,j) .gt. maxv) then
		maxv=a(i,j)
		endif	
	enddo
enddo

write(*,*) 'Максимальное значение.'
write(*,*) 'Результат my_maxval:', maxv
write(*,*) 'Результат maxval:', maxval(a)

elseif (dim .gt. 0) then

do i=1, se1, 1
maxarr1(i)=-9999999
enddo

do i=1, se2, 1
maxarr2(i)=-9999999
enddo

if (dim .eq. 1) then

do j=1, se2, 1
	do i=1, se1, 1
		if (a(i,j) .gt. maxarr2(j)) then
		maxarr2(j)=a(i,j)
		endif
	enddo
enddo

elseif (dim .eq. 2) then

do i=1, se1, 1
	do j=1, se2, 1
		if (a(i,j) .gt. maxarr1(i)) then
		maxarr1(i)=a(i,j)
		endif
	enddo
enddo

endif

write(*,*) ' '
write(*,*) 'Максимальное значение по измерению.'

if (dim .eq. 1) then
write(*,*) 'Результат my_maxval:', maxarr2
elseif (dim .eq. 2) then
write(*,*) 'Результат my_maxval:', maxarr1
endif
write(*,*) 'Результат maxval:', maxval(a, dim=dim)

endif

write(*,*) ' '

!write(*,*) my_maxval(a(1,:), se1)

end function my_maxval

!-----------------------------------

function mask_maker(a,sg1,sg2,dim2)
integer a(:,:), mask(sg1,sg2), sg1, sg2, dim2
integer minarr1(sg1), minarr2(sg2)
integer j_minarr1(sg1), i_minarr2(sg2)
integer i, j, minv, i_minv, j_minv
integer mask_maker

!Изменяйте ваше условие маски в данной конструкции:
mask=a>0
!а также в вызове встроенной функции minloc ниже.

write(*,*) 'Было введено условие для маски. Ищите &
его описание в тексте программы.'
write(*,*) ' '
write(*,1007) ((mask(i,j),j=1,sg2),i=1,sg1)

if (dim2 .eq. 0) then
minv=99999999
do i=1, sg1, 1
	do j=1, sg2, 1
		if (a(i,j) .lt. minv .and. mask(i,j) .eq. 1) then
		minv=a(i,j)
		i_minv=i
		j_minv=j
		endif
	enddo
enddo

write(*,*) 'Найденное минимальное значение:', minv
write(*,*) ' '
write(*,*) 'Его индексы:'
write(*,*) 'Результат my_minloc:', i_minv, j_minv
write(*,*) 'Результат minloc:', minloc(a,mask=a>0)

elseif (dim2 .gt. 0) then

do i=1, sg1, 1
minarr1(i)=99999999
enddo

do i=1, sg2, 1
minarr2(i)=99999999
enddo

if (dim2 .eq. 1) then

do j=1, sg2, 1
	do i=1, sg1, 1
		if (a(i,j) .lt. minarr2(j) .and. mask(i,j) .eq. 1) then
		minarr2(j)=a(i,j)
		i_minarr2(j)=i
		endif
	enddo
enddo

write(*,*) 'Поиск по измерению:', dim2
write(*,*) 'Найденные минимальные значения:', minarr2
write(*,*) 'Их индексы вдоль этого измерения:'
write(*,*) 'Результат my_minloc:', i_minarr2
write(*,*) 'Результат minloc:   ', minloc(a,dim=1,mask=a>0)

elseif (dim2 .eq. 2) then

do i=1, sg1, 1
	do j=1, sg2, 1
		if (a(i,j) .lt. minarr1(i) .and. mask(i,j) .eq. 1) then
		minarr1(i)=a(i,j)
		j_minarr1(i)=j
		endif
	enddo
enddo

write(*,*) 'Поиск по измерению:', dim2
write(*,*) 'Найденные минимальные значения:', minarr1
write(*,*) 'Их индексы вдоль этого измерения:'
write(*,*) 'Результат my_minloc:', j_minarr1
write(*,*) 'Результат minloc:   ', minloc(a,dim=2,mask=a>0)


endif
endif

1007 format (1x, 'Маска:',(/5(10x,4i15/)))

write(*,*) ' '

end function mask_maker

function my_minloc(a,sg1,sg2, mask, dim2)
integer a(:,:), sg1, sg2, minv, i_minv, j_minv, dim2
integer minarr1(sg1), minarr2(sg2)
integer j_minarr1(sg1), i_minarr2(sg2)
integer my_minloc, i, j
integer mask

if (mask .eq. 0) then
if (dim2 .eq. 0) then

minv=99999999
do i=1, sg1, 1
	do j=1, sg2, 1
		if (a(i,j) .lt. minv) then
		minv=a(i,j)
		i_minv=i
		j_minv=j
		endif
	enddo
enddo

write(*,*) 'Найденное минимальное значение:', minv
write(*,*) ' '
write(*,*) 'Его индексы:'
write(*,*) 'Результат my_minloc:', i_minv, j_minv
write(*,*) 'Результат minloc:', minloc(a)

elseif (dim2 .gt. 0) then

do i=1, sg1, 1
minarr1(i)=99999999
enddo

do i=1, sg2, 1
minarr2(i)=99999999
enddo

if (dim2 .eq. 1) then

do j=1, sg2, 1
	do i=1, sg1, 1
		if (a(i,j) .lt. minarr2(j)) then
		minarr2(j)=a(i,j)
		i_minarr2(j)=i
		endif
	enddo
enddo

write(*,*) 'Поиск по измерению:', dim2
write(*,*) 'Найденные минимальные значения:', minarr2
write(*,*) 'Их индексы вдоль этого измерения:'
write(*,*) 'Результат my_minloc:', i_minarr2
write(*,*) 'Результат minloc:   ', minloc(a,dim=1)

elseif (dim2 .eq. 2) then

do i=1, sg1, 1
	do j=1, sg2, 1
		if (a(i,j) .lt. minarr1(i)) then
		minarr1(i)=a(i,j)
		j_minarr1(i)=j
		endif
	enddo
enddo

write(*,*) 'Поиск по измерению:', dim2
write(*,*) 'Найденные минимальные значения:', minarr1
write(*,*) 'Их индексы вдоль этого измерения:'
write(*,*) 'Результат my_minloc:', j_minarr1
write(*,*) 'Результат minloc:   ', minloc(a,dim=2)

endif
endif

elseif (mask .ne. 0) then

write(*,*) mask_maker(a,sg1,sg2,dim2)

endif

end function my_minloc

end module myarray
