program main
implicit none

!-----------QUADRA---------------

!ИНТЕРФЕЙСЫ

interface rectan

        function rectan_real_4(f_real_4,a,b,n) result(integral)
        real(4), external :: f_real_4
        real(4) a, b, h, x, summa, integral
        integer n, i

        end function rectan_real_4

        function rectan_real_8(f_real_8,a,b,n) result(integral)
        real(8), external :: f_real_8
        real(8) a, b, h, x, summa, integral
        integer n, i

        end function rectan_real_8

        function rectan_real_10(f_real_10,a,b,n) result(integral)
        real(10), external :: f_real_10
        real(10) a, b, h, x, summa, integral
        integer n, i

        end function rectan_real_10

        function rectan_real_16(f_real_16,a,b,n) result(integral)
        real(16), external :: f_real_16
        real(16) a, b, h, x, summa, integral
        integer n, i

        end function rectan_real_16

end interface rectan

interface trap

        function trap_real_4(f_real_4,a,b,n) result(integral)
        real(4), external :: f_real_4
        real(4) a, b, h, x, summa, integral
        integer n, i

        end function trap_real_4

        function trap_real_8(f_real_8,a,b,n) result(integral)
        real(8), external :: f_real_8
        real(8) a, b, h, x, summa, integral
        integer n, i

        end function trap_real_8

        function trap_real_10(f_real_10,a,b,n) result(integral)
        real(10), external :: f_real_10
        real(10) a, b, h, x, summa, integral
        integer n, i

        end function trap_real_10

        function trap_real_16(f_real_16,a,b,n) result(integral)
        real(16), external :: f_real_16
        real(16) a, b, h, x, summa, integral
        integer n, i

        end function trap_real_16

end interface trap

interface sim

        function sim_real_4(f_real_4,a,b,n) result(integral)
        real(4), external :: f_real_4
        real(4) a, b, h, x, summa, integral
        integer n, i

        end function sim_real_4

        function sim_real_8(f_real_8,a,b,n) result(integral)
        real(8), external :: f_real_8
        real(8) a, b, h, x, summa, integral
        integer n, i

        end function sim_real_8

        function sim_real_10(f_real_10,a,b,n) result(integral)
        real(10), external :: f_real_10
        real(10) a, b, h, x, summa, integral
        integer n, i

        end function sim_real_10

        function sim_real_16(f_real_16,a,b,n) result(integral)
        real(16), external :: f_real_16
        real(16) a, b, h, x, summa, integral
        integer n, i

        end function sim_real_16

end interface sim

!ПОДИНТЕГРАЛЬНЫЕ ФУНКЦИИ

external f_real_4
external f_real_8
external f_real_10
external f_real_16
real(4) f_real_4
real(8) f_real_8
real(10) f_real_10
real(16) f_real_16

!УЧАСТКИ ИНТЕГРИРОВАНИЯ

integer n_r
                          
real(4) a_r_r4, b_r_r4
real(8) a_r_r8, b_r_r8
real(10) a_r_r10, b_r_r10
real(16) a_r_r16, b_r_r16


integer n_t 
                  
real(4) a_t_r4, b_t_r4
real(8) a_t_r8, b_t_r8
real(10) a_t_r10, b_t_r10
real(16) a_t_r16, b_t_r16

integer n_s 
                                
real(4) a_s_r4, b_s_r4
real(8) a_s_r8, b_s_r8
real(10) a_s_r10, b_s_r10
real(16) a_s_r16, b_s_r16

!ПАРАМЕТРЫ

!Пусть на подинтегральную функцию действуют смещения p, q и r;

real(4) p_r4, q_r4, r_r4
real(8) p_r8, q_r8, r_r8
real(10) p_r10, q_r10, r_r10
real(16) p_r16, q_r16, r_r16
common / parameters_real_4 / p_r4, q_r4, r_r4
common / parameters_real_8 / p_r8, q_r8, r_r8
common / parameters_real_10 / p_r10, q_r10, r_r10
common / parameters_real_16 / p_r16, q_r16, r_r16

!ВРЕМЯ

real(8) t1_r4, t2_r4, t_rectan_r4, t3_r4, t4_r4, t_trap_r4, t5_r4, t6_r4, t_sim_r4
real(8) t1_r8, t2_r8, t_rectan_r8, t3_r8, t4_r8, t_trap_r8, t5_r8, t6_r8, t_sim_r8
real(8) t1_r10, t2_r10, t_rectan_r10, t3_r10, t4_r10, t_trap_r10, t5_r10, t6_r10, t_sim_r10
real(8) t1_r16, t2_r16, t_rectan_r16, t3_r16, t4_r16, t_trap_r16, t5_r16, t6_r16, t_sim_r16

!ОПРЕДЕЛИТЕЛЬ ТИПА

! ВНИМАНИЕ

! Для смены типа переменных в вычислениях измените первое значение (real_type)
! в файле input. Возможные значения:
!          4 - вычислять по типу REAL 4
!          8 - вычислять по типу REAL 8
!          10 - вычислять по типу REAL 10
!          16 - вычислять по типу REAL 16
!          любое другое целое - вычислять по ВСЕМ типам REAL

integer real_type
    
!-----------QUADRA---------------

write(*,'(/,a,/)') ' QUADRA '

!-----------------------------------------------------------------------------

read(*,*) real_type

if (real_type .eq. 4 .or. real_type .eq. 8 .or. real_type .eq. 10 .or. real_type .eq. 16) then
        write(*,'(a,i2,/)') ' ВЫБРАННЫЙ ТИП - REAL', real_type
else
        write(*,'(a,/)') ' ВЫБРАННЫЙ ТИП - ВЫЧИСЛЕНИЕ ПО ВСЕМ ТИПАМ'
endif

read(*,*) p_r4, q_r4, r_r4 
read(*,*) p_r8, q_r8, r_r8   
read(*,*) p_r10, q_r10, r_r10
read(*,*) p_r16, q_r16, r_r16

write(*,'(a,/)') ' Действующие параметры (смещение по x):'
write(*,'(a,3x,a,9x,a,9x,a,9x,a)') ' Тип', 'real_4', 'real_8', 'real_10', 'real_16'
write(*,'(a,e15.7,e15.7,e15.7,e15.7)') ' p = ', p_r4, p_r8, p_r10, p_r16
write(*,'(a,e15.7,e15.7,e15.7,e15.7)') ' q = ', q_r4, q_r8, q_r10, q_r16
write(*,'(a,e15.7,e15.7,e15.7,e15.7,/,/)') ' r = ', r_r4, r_r8, r_r10, r_r16

read(*,*) n_r
read(*,*) a_r_r4, b_r_r4
read(*,*) a_r_r8, b_r_r8
read(*,*) a_r_r10, b_r_r10
read(*,*) a_r_r16, b_r_r16

read(*,*) n_t
read(*,*) a_t_r4, b_t_r4
read(*,*) a_t_r8, b_t_r8
read(*,*) a_t_r10, b_t_r10
read(*,*) a_t_r16, b_t_r16

read(*,*) n_s
read(*,*) a_s_r4, b_s_r4
read(*,*) a_s_r8, b_s_r8
read(*,*) a_s_r10, b_s_r10
read(*,*) a_s_r16, b_s_r16



write(*,'(a,/)') ' Формула средних прямоугольников '

if (real_type .eq. 4) then
        write(*,*) 'Был задан промежуток: от ', a_r_r4, 'до ', b_r_r4, '.'
elseif (real_type .eq. 8) then
        write(*,*) 'Был задан промежуток: от ', a_r_r8, 'до ', b_r_r8, '.'
elseif (real_type .eq. 10) then
        write(*,*) 'Был задан промежуток: от ', a_r_r10, 'до ', b_r_r10, '.'
elseif (real_type .eq. 16) then
        write(*,*) 'Был задан промежуток: от ', a_r_r16, 'до ', b_r_r16, '.'
else
        write(*,'(a,/)') ' Заданные промежутки: '
        write(*,'(a,5x,a,9x,a,9x,a,9x,a)') ' Тип', 'real_4', 'real_8', 'real_10', 'real_16'
        write(*,'(a,e15.7,e15.7,e15.7,e15.7)') ' a_r = ', a_r_r4, a_r_r8, a_r_r10, a_r_r16
        write(*,'(a,e15.7,e15.7,e15.7,e15.7,/)') ' b_r = ', b_r_r4, b_r_r8, b_r_r10, b_r_r16
endif

write(*,'(a,i3,3x,a,/)') ' Будут рассмотрены первые ', n_r, 'значений(-я).'

write(*,*) 'Rectan:'

if (real_type .eq. 4) then
        call cpu_time(t1_r4)
        write(*,*) 'Integral=~', rectan(f_real_4, a_r_r4, b_r_r4, n_r)
        call cpu_time(t2_r4)
elseif (real_type .eq. 8) then
        call cpu_time(t1_r8)
        write(*,*) 'Integral=~', rectan(f_real_8, a_r_r8, b_r_r8, n_r)
        call cpu_time(t2_r8)
elseif (real_type .eq. 10) then
        call cpu_time(t1_r10)
        write(*,*) 'Integral=~', rectan(f_real_10, a_r_r10, b_r_r10, n_r)
        call cpu_time(t2_r10)
elseif (real_type .eq. 16) then
        call cpu_time(t1_r16)
        write(*,*) 'Integral=~', rectan(f_real_16, a_r_r16, b_r_r16, n_r)
        call cpu_time(t2_r16)
else
        call cpu_time(t1_r4)
        write(*,*) 'Integral (REAL 4) =~', rectan(f_real_4, a_r_r4, b_r_r4, n_r)
        call cpu_time(t2_r4)

        call cpu_time(t1_r8)
        write(*,*) 'Integral (REAL 8) =~', rectan(f_real_8, a_r_r8, b_r_r8, n_r)
        call cpu_time(t2_r8)

        call cpu_time(t1_r10)
        write(*,*) 'Integral (REAL 10)=~', rectan(f_real_10, a_r_r10, b_r_r10, n_r)
        call cpu_time(t2_r10)

        call cpu_time(t1_r16)
        write(*,*) 'Integral (REAL 16)=~', rectan(f_real_16, a_r_r16, b_r_r16, n_r)
        call cpu_time(t2_r16)
endif

write(*,'(/)')

!-----------------------------------------------------------------------------

write(*,'(/,a,/)') ' Формула трапеций '

if (real_type .eq. 4) then
        write(*,*) 'Был задан промежуток: от ', a_t_r4, 'до ', b_t_r4, '.'
elseif (real_type .eq. 8) then
        write(*,*) 'Был задан промежуток: от ', a_t_r8, 'до ', b_t_r8, '.'
elseif (real_type .eq. 10) then
        write(*,*) 'Был задан промежуток: от ', a_t_r10, 'до ', b_t_r10, '.'
elseif (real_type .eq. 16) then
        write(*,*) 'Был задан промежуток: от ', a_t_r16, 'до ', b_t_r16, '.'
else
        write(*,'(a,/)') ' Заданные промежутки: '
        write(*,'(a,5x,a,9x,a,9x,a,9x,a)') ' Тип', 'real_4', 'real_8', 'real_10', 'real_16'
        write(*,'(a,e15.7,e15.7,e15.7,e15.7)') ' a_t = ', a_t_r4, a_t_r8, a_t_r10, a_t_r16
        write(*,'(a,e15.7,e15.7,e15.7,e15.7,/)') ' b_t = ', b_t_r4, b_t_r8, b_t_r10, b_t_r16
endif

write(*,'(a,i3,3x,a,/)') ' Будут рассмотрены первые ', n_t+1, 'значений(-я).'

write(*,*) 'Trap:'

if (real_type .eq. 4) then
        call cpu_time(t3_r4)
        write(*,*) 'Integral=~', trap(f_real_4,a_t_r4,b_t_r4,n_t) 
        call cpu_time(t4_r4)
elseif (real_type .eq. 8) then
        call cpu_time(t3_r8)
        write(*,*) 'Integral=~', trap(f_real_8,a_t_r8,b_t_r8,n_t) 
        call cpu_time(t4_r8)
elseif (real_type .eq. 10) then
        call cpu_time(t3_r10)
        write(*,*) 'Integral=~', trap(f_real_10,a_t_r10,b_t_r10,n_t)
        call cpu_time(t4_r10)
elseif (real_type .eq. 16) then
        call cpu_time(t3_r16)
        write(*,*) 'Integral=~', trap(f_real_16,a_t_r16,b_t_r16,n_t) 
        call cpu_time(t4_r16)
else
        call cpu_time(t3_r4)
        write(*,*) 'Integral (REAL 4) =~', trap(f_real_4, a_t_r4, b_t_r4, n_t)
        call cpu_time(t4_r4)

        call cpu_time(t3_r8)
        write(*,*) 'Integral (REAL 8) =~', trap(f_real_8, a_t_r8, b_t_r8, n_t)
        call cpu_time(t4_r8)

        call cpu_time(t3_r10)
        write(*,*) 'Integral (REAL 10)=~', trap(f_real_10, a_t_r10, b_t_r10, n_t)
        call cpu_time(t4_r10)

        call cpu_time(t3_r16)
        write(*,*) 'Integral (REAL 16)=~', trap(f_real_16, a_t_r16, b_t_r16, n_t)
        call cpu_time(t4_r16)
endif

write(*,'(/)')

!-----------------------------------------------------------------------------

write(*,'(/,a,/)') ' Формула Симпсона (парабол) '

if (real_type .eq. 4) then
        write(*,*) 'Был задан промежуток: от ', a_s_r4, 'до ', b_s_r4, '.'
elseif (real_type .eq. 8) then
        write(*,*) 'Был задан промежуток: от ', a_s_r8, 'до ', b_s_r8, '.'
elseif (real_type .eq. 10) then
        write(*,*) 'Был задан промежуток: от ', a_s_r10, 'до ', b_s_r10, '.'
elseif (real_type .eq. 16) then
        write(*,*) 'Был задан промежуток: от ', a_s_r16, 'до ', b_s_r16, '.'
else
        write(*,'(a,/)') ' Заданные промежутки: '
        write(*,'(a,5x,a,9x,a,9x,a,9x,a)') ' Тип', 'real_4', 'real_8', 'real_10', 'real_16'
        write(*,'(a,e15.7,e15.7,e15.7,e15.7)') ' a_s = ', a_s_r4, a_s_r8, a_s_r10, a_s_r16
        write(*,'(a,e15.7,e15.7,e15.7,e15.7,/)') ' b_s = ', b_s_r4, b_s_r8, b_s_r10, b_s_r16
endif

write(*,'(a,i3,3x,a,/)') ' Будут рассмотрены первые ', n_s+1, 'значений(-я).'

write(*,*) 'Sim:'

if (real_type .eq. 4) then
        call cpu_time(t5_r4)
        write(*,*) 'Integral=~', sim(f_real_4,a_s_r4,b_s_r4,n_s) 
        call cpu_time(t6_r4)
elseif (real_type .eq. 8) then
        call cpu_time(t5_r8)
        write(*,*) 'Integral=~', sim(f_real_8,a_s_r8,b_s_r8,n_s) 
        call cpu_time(t6_r8)
elseif (real_type .eq. 10) then
        call cpu_time(t5_r10)
        write(*,*) 'Integral=~', sim(f_real_10,a_s_r10,b_s_r10,n_s)
        call cpu_time(t6_r10)
elseif (real_type .eq. 16) then
        call cpu_time(t5_r16)
        write(*,*) 'Integral=~', sim(f_real_16,a_s_r16,b_s_r16,n_s) 
        call cpu_time(t6_r16)
else
        call cpu_time(t5_r4)
        write(*,*) 'Integral (REAL 4) =~', sim(f_real_4, a_s_r4, b_s_r4, n_s)
        call cpu_time(t6_r4)

        call cpu_time(t5_r8)
        write(*,*) 'Integral (REAL 8) =~', sim(f_real_8, a_s_r8, b_s_r8, n_s)
        call cpu_time(t6_r8)

        call cpu_time(t5_r10)
        write(*,*) 'Integral (REAL 10)=~', sim(f_real_10, a_s_r10, b_s_r10, n_s)
        call cpu_time(t6_r10)

        call cpu_time(t5_r16)
        write(*,*) 'Integral (REAL 16)=~', sim(f_real_16, a_s_r16, b_s_r16, n_s)
        call cpu_time(t6_r16)
endif

write(*,'(/)')

!-----------------------------------------------------------------------------

write(*,'(a,/)') ' Временные затраты (в единицах временных затрат rectan)'

write(*,*) '------------------------------------------------------------'
write(*,'(a,4x,a,5x,a,35x,a)') ' Время', ':', 'Задача N1', ':'
write(*,*) '------------------------------------------------------------'
write(*,'(a,6x,a,5x,a,8x,a,14x,a,9x,a)') ' Тип', ':', 'rectan', 'trap', 'sim', ':'
write(*,*) '------------------------------------------------------------'


if (real_type .eq. 4) then
        t_rectan_r4=t2_r4-t1_r4
        t_trap_r4=(t4_r4-t3_r4)/t_rectan_r4
        t_sim_r4=(t6_r4-t5_r4)/t_rectan_r4
        write(*,'(a,2x,a,4x,i3,6x,e15.7,3x,e15.7,3x,a)') ' real(4)', ':', 1, t_trap_r4, t_sim_r4, ':'
elseif (real_type .eq. 8) then
        t_rectan_r8=t2_r8-t1_r8
        t_trap_r8=(t4_r8-t3_r8)/t_rectan_r8
        t_sim_r8=(t6_r8-t5_r8)/t_rectan_r8
        write(*,'(a,2x,a,4x,i3,6x,e15.7,3x,e15.7,3x,a)') ' real(8)', ':', 1, t_trap_r8, t_sim_r8, ':'
elseif (real_type .eq. 10) then
        t_rectan_r10=t2_r10-t1_r10
        t_trap_r10=(t4_r10-t3_r10)/t_rectan_r10
        t_sim_r10=(t6_r10-t5_r10)/t_rectan_r10
        write(*,'(a,1x,a,4x,i3,6x,e15.7,3x,e15.7,3x,a)') ' real(10)', ':', 1, t_trap_r10, t_sim_r10, ':'
elseif (real_type .eq. 16) then
        t_rectan_r16=t2_r16-t1_r16
        t_trap_r16=(t4_r16-t3_r16)/t_rectan_r16
        t_sim_r16=(t6_r16-t5_r16)/t_rectan_r16
        write(*,'(a,1x,a,4x,i3,6x,e15.7,3x,e15.7,3x,a)') ' real(16)', ':', 1, t_trap_r16, t_sim_r16, ':'
else
        t_rectan_r4=t2_r4-t1_r4
        t_rectan_r8=t2_r8-t1_r8
        t_rectan_r10=t2_r10-t1_r10
        t_rectan_r16=t2_r16-t1_r16

        t_trap_r4=(t4_r4-t3_r4)/t_rectan_r4
        t_trap_r8=(t4_r8-t3_r8)/t_rectan_r8
        t_trap_r10=(t4_r10-t3_r10)/t_rectan_r10
        t_trap_r16=(t4_r16-t3_r16)/t_rectan_r16

        t_sim_r4=(t6_r4-t5_r4)/t_rectan_r4
        t_sim_r8=(t6_r8-t5_r8)/t_rectan_r8
        t_sim_r10=(t6_r10-t5_r10)/t_rectan_r10
        t_sim_r16=(t6_r16-t5_r16)/t_rectan_r16

        write(*,'(a,2x,a,4x,i3,6x,e15.7,3x,e15.7,3x,a)') ' real(4)', ':', 1, t_trap_r4, t_sim_r4, ':'
        write(*,'(a,2x,a,4x,i3,6x,e15.7,3x,e15.7,3x,a)') ' real(8)', ':', 1, t_trap_r8, t_sim_r8, ':'
        write(*,'(a,1x,a,4x,i3,6x,e15.7,3x,e15.7,3x,a)') ' real(10)', ':', 1, t_trap_r10, t_sim_r10, ':'
        write(*,'(a,1x,a,4x,i3,6x,e15.7,3x,e15.7,3x,a)') ' real(16)', ':', 1, t_trap_r16, t_sim_r16, ':'
endif

write(*,'(//)')

end
