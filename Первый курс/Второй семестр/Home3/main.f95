program main
use poly
use binrad
use quadra

integer n, i, x, a, b, v, n_r
real y, y2, y4, yn4
real(8) a_r, b_r, x_r, h_r
real,dimension(100)::m
integer, dimension(100)::t
integer, dimension(8)::binar
real(8), dimension(1000)::domain
read(*,*) n
write(*,*) 'Введенная степень: ', n

read(*,*) (m(i),i=1,n+1)
write(*,*) 'Введенные коэффициенты, начиная со старшей степени: '
write(*,'(i3,e15.7)') (i,m(i),i=1,n+1)

read(*,*) x
write(*,*) 'Введенная точка, в которой хотите узнать значение: ', x

write(*,*) 'Poly1: '
write(*,*) poly1(n,x,m)

write(*,*) 'Poly2: '
write(*,*) poly2(x, m, n)

read(*,*) (t(i),i=1,n)
write(*,*) 'Введен массив: '
write(*,'(i3,i15)') (i,t(i),i=1,n)

write(*,*) 'Poly3: '
write(*,*) poly3(n, t)

read(*,*) v
write(*,*) 'Было введено целое число v: ', v

write(*,*) 'Poly4: '
write(*,*) poly4(v)

read(*,*) binar
write(*,*) 'Был введен набор нулей и единиц:', binar

write(*,*) 'Binrad5:'
write(*,*) binrad5(binar)

write(*,*) 'Binrad6:'
write(*,*) binrad6(binar)

!quadra

read(*,*) a_r, b_r, n_r
write(*,*) 'Был задан промежуток: от ', a_r, 'до ', b_r, '.'
write(*,*) 'Будут рассмотрены первые ', n_r, 'значений(-я).'

h_r=(b_r-a_r)/n_r

do i=1, size(domain), 1
	domain(i)=0
enddo

do i=1, n_r, 1
	x_r=a_r+(i-1)*h_r
	domain(i)=your_function(x_r)
enddo

write(*,*) 'Был получен массив значений domain:'
write(*,'(i3,e15.7, e15.7)') (i,domain(i),i=1,n_r)

write(*,*) 'Rectan:'
write(*,*) rectan(domain, a_r, b_r, n_r) 


end
