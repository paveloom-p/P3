program h4_4
implicit none
integer n, i, k, x
real y
integer,dimension(1)::a
data a / 0 /
write(*,*) 'Введите степень.'
read(*,*) n
write(*,*) 'Введите коэффициенты, начиная со старшего члена.'

do i=1, n+1, 1
read(*,*) k           !от младшей степени до старшей
a(i)=k
enddo

write(*,*) 'Введите точку, в которой хотите узнать значение.'
read(*,*) x

y=a(1)
do i=1, n, 1
y=y*x+a(i+1)
enddo

write(*,*) 'Значение многочлена в этой точке: ', y
end
