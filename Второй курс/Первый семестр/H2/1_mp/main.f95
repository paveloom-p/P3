program main
use real_type
implicit none

!-----------QUADRA---------------

!ИНТЕРФЕЙС

interface 

function rectan(f,a,b,n) result(integral)
use real_type
implicit none
real(mp), external :: f
real(mp) a, b, h, x, summa, integral
integer n, i

end function rectan

function trap(f,a,b,n) result(integral)
use real_type
implicit none
real(mp), external :: f
real(mp) a, b, h, x, summa, integral
integer n, i

end function trap

function sim(f,a,b,n) result(integral)
use real_type
implicit none
real(mp), external :: f
real(mp) a, b, h, x, summa_ch, summa_nech, summa, integral
integer n, i

end function sim

function f(x) result(y)
use real_type
implicit none
real(mp) p, q, r
common / parameters / p, q, r
real(mp) x, y
                              
end function f

end interface

!УЧАСТКИ ИНТЕГРИРОВАНИЯ

integer n_r                             
real(mp) a_r, b_r

integer n_t
real(mp) a_t, b_t

integer n_s
real(mp) a_s, b_s

!ВРЕМЯ

real(8) t1, t2, t_rectan, t3, t4, t_trap, t5, t6, t_sim

!ПАРАМЕТРЫ

!Пусть на подинтегральную функцию действуют смещения p, q и r;

real(mp) p, q, r
common / parameters / p, q, r

!-----------QUADRA---------------

write(*,'(/,a,/)') ' QUADRA '

!-----------------------------------------------------------------------------

!УЧАСТКИ ИНТЕГРИРОВАНИЯ

read(*,*) a_r, b_r, n_r
read(*,*) a_t, b_t, n_t
read(*,*) a_s, b_s, n_s

!ПАРАМЕТРЫ

read(*,*) p, q, r
write(*,*) 'Действующие параметры:'
write(*,*) 'p = ', p
write(*,*) 'q = ', q
write(*,*) 'r = ', r
write(*,'(/)')

!ОПИСАНИЕ И ВЫЧИСЛЕНИЕ

write(*,'(a,/)') ' Формула средних прямоугольников '

write(*,*) 'Был задан промежуток: от ', a_r, 'до ', b_r, '.'
write(*,*) 'Будут рассмотрены первые ', n_r, 'значений(-я).'

write(*,*) 'Rectan:'
call cpu_time(t1)
write(*,*) 'Integral=~', rectan(f, a_r, b_r, n_r)
call cpu_time(t2)

write(*,'(//)')

!-----------------------------------------------------------------------------

write(*,'(/,a,/)') ' Формула трапеций '

write(*,*) 'Был задан промежуток: от ', a_t, 'до ', b_t, '.'
write(*,*) 'Будут рассмотрены первые ', n_t+1, 'значений(-я).'

write(*,*) 'Trap:'
call cpu_time(t3)
write(*,*) 'Integral=~', trap(f, a_t, b_t, n_t) 
call cpu_time(t4)

write(*,'(//)')

!-----------------------------------------------------------------------------

write(*,'(/,a,/)') ' Формула Симпсона (парабол) '

write(*,*) 'Был задан промежуток: от ', a_s, 'до ', b_s, '.'
write(*,*) 'Будут рассмотрены первые ', n_s+1, 'значений(-я).'

write(*,*) 'Sim:'
call cpu_time(t5)
write(*,*) 'Integral=~', sim(f, a_s, b_s, n_s)
call cpu_time(t6)

write(*,'(//)')

!-----------------------------------------------------------------------------

!ВРЕМЯ

t_rectan=t2-t1
t_trap=(t4-t3)/t_rectan
t_sim=(t6-t5)/t_rectan

write(*,'(a,/)') ' Временные затраты (в единицах временных затрат rectan)'

write(*,*) '------------------------------------------------------------'
write(*,'(a,4x,a,5x,a,35x,a)') ' Время', ':', 'Задача N1', ':'
write(*,*) '------------------------------------------------------------'
write(*,'(a,6x,a,5x,a,8x,a,14x,a,9x,a)') ' Тип', ':', 'rectan', 'trap', 'sim', ':'
write(*,*) '------------------------------------------------------------'
write(*,'(a,i2,a,1x,a,7x,a,5x,e15.7,3x,e15.7,3x,a)') ' real(', mp, ')', ':', '1', t_trap, t_sim, ':'
write(*,'(a,///)') ' ------------------------------------------------------------'

end
