program main
use quadr
use functions
implicit none

!-----------QUADRA---------------

!ИНТЕРФЕЙС

interface quadra
module procedure quadra_real_4, quadra_real_8, quadra_real_10, quadra_real_16
end interface

!УЧАСТОК ИНТЕГРИРОВАНИЯ

integer n
                       
real(4) a_r4, b_r4
real(8) a_r8, b_r8
real(10) a_r10, b_r10
real(16) a_r16, b_r16

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

!ОПРЕДЕЛИТЕЛЬ ФУНКЦИИ

! Для смены типа функции, используемой в вычислениях, измените второе значение
! (function_type) в файле input. Возможные значения:
!          rectan - вычислять, используя rectan
!          trap - вычислять, используя trap
!          sim - вычислять, используя sim
!          любая другая символьная переменная - вычислить, используя ВСЕ типы функций

character(6) function_type

!-----------QUADRA---------------

write(*,'(/,a,/)') ' QUADRA '

!-----------------------------------------------------------------------------

read(*,*) real_type

if (real_type .eq. 4 .or. real_type .eq. 8 .or. real_type .eq. 10 .or. real_type .eq. 16) then
        write(*,'(a,i2)') ' ВЫБРАННЫЙ ТИП ДАННЫХ - REAL', real_type
else
        write(*,'(a)') ' ВЫБРАННЫЙ ТИП ДАННЫХ - ВЫЧИСЛЕНИЕ ПО ВСЕМ ТИПАМ'
endif

read(*,*) function_type

if (function_type .eq. 'rectan' .or. function_type .eq. 'trap' .or. function_type .eq. 'sim') then
        write(*,'(a,a,/)') ' ВЫБРАННЫЙ ТИП ФУНКЦИИ - ', function_type
else
        write(*,'(a,/)') ' ВЫБРАННЫЙ ТИП ФУНКЦИИ - ВЫЧИСЛЕНИЕ ПО ВСЕМ ТИПАМ'
endif

read(*,*) p_r4, q_r4, r_r4 
read(*,*) p_r8, q_r8, r_r8   
read(*,*) p_r10, q_r10, r_r10
read(*,*) p_r16, q_r16, r_r16

write(*,'(a,/)') ' Действующие параметры (смещение по x):'
write(*,'(a,3x,a,9x,a,9x,a,8x,a)') ' Тип', 'real_4', 'real_8', 'real_10', 'real_16'
write(*,'(a,e15.7,e15.7,e15.7,e15.7)') ' p = ', p_r4, p_r8, p_r10, p_r16
write(*,'(a,e15.7,e15.7,e15.7,e15.7)') ' q = ', q_r4, q_r8, q_r10, q_r16
write(*,'(a,e15.7,e15.7,e15.7,e15.7,/,/)') ' r = ', r_r4, r_r8, r_r10, r_r16

read(*,*) n

read(*,*) a_r4, b_r4
read(*,*) a_r8, b_r8
read(*,*) a_r10, b_r10
read(*,*) a_r16, b_r16

if (real_type .eq. 4) then
        write(*,*) 'Был задан промежуток: от ', a_r4, 'до ', b_r4, '.'
elseif (real_type .eq. 8) then
        write(*,*) 'Был задан промежуток: от ', a_r8, 'до ', b_r8, '.'
elseif (real_type .eq. 10) then
        write(*,*) 'Был задан промежуток: от ', a_r10, 'до ', b_r10, '.'
elseif (real_type .eq. 16) then
        write(*,*) 'Был задан промежуток: от ', a_r16, 'до ', b_r16, '.'
else
        write(*,'(a,/)') ' Заданные промежутки: '
        write(*,'(a,3x,a,9x,a,9x,a,8x,a)') ' Тип', 'real_4', 'real_8', 'real_10', 'real_16'
        write(*,'(a,e15.7,e15.7,e15.7,e15.7)') ' a = ', a_r4, a_r8, a_r10, a_r16
        write(*,'(a,e15.7,e15.7,e15.7,e15.7,/)') ' b = ', b_r4, b_r8, b_r10, b_r16
endif

write(*,'(a,i3,2x,a,i3,2x,a,/)') ' Будут рассмотрены первые ', n, &
'значений(-я) (', n+1, ', если trap или sim).'

write(*,'(a,/)') ' Обратите внимание: для оценки временных затрат функция rectan &
&будет вычисляться ВСЕГДА.'


if (real_type .eq. 4) then
        call cpu_time(t1_r4)
        write(*,*) 'Integral (RECTAN) =~', quadra(rectan_real_4, f_real_4, a_r4, b_r4, n)
        call cpu_time(t2_r4)
elseif (real_type .eq. 8) then
        call cpu_time(t1_r8)
        write(*,*) 'Integral (RECTAN) =~', quadra(rectan_real_8, f_real_8, a_r8, b_r8, n)
        call cpu_time(t2_r8)
elseif (real_type .eq. 10) then
        call cpu_time(t1_r10)
        write(*,*) 'Integral (RECTAN) =~', quadra(rectan_real_10, f_real_10, a_r10, b_r10, n)
        call cpu_time(t2_r10)
elseif (real_type .eq. 16) then
        call cpu_time(t1_r16)
        write(*,*) 'Integral (RECTAN) =~', quadra(rectan_real_16, f_real_16, a_r16, b_r16, n)
        call cpu_time(t2_r16)
else
        call cpu_time(t1_r4)
        write(*,*) 'Integral (RECTAN) (REAL 4) =~', quadra(rectan_real_4, f_real_4, a_r4, b_r4, n)
        call cpu_time(t2_r4)

        call cpu_time(t1_r8)
        write(*,*) 'Integral (RECTAN) (REAL 8) =~', quadra(rectan_real_8, f_real_8, a_r8, b_r8, n)
        call cpu_time(t2_r8)

        call cpu_time(t1_r10)
        write(*,*) 'Integral (RECTAN) (REAL 10)=~', quadra(rectan_real_10, f_real_10, a_r10, b_r10, n)
        call cpu_time(t2_r10)

        call cpu_time(t1_r16)
        write(*,*) 'Integral (RECTAN) (REAL 16)=~', quadra(rectan_real_16, f_real_16, a_r16, b_r16, n)
        call cpu_time(t2_r16)
endif

if (function_type .eq. 'trap') then

if (real_type .eq. 4) then
        call cpu_time(t3_r4)
        write(*,*) 'Integral (TRAP) =~', quadra(trap_real_4, f_real_4, a_r4, b_r4, n)
        call cpu_time(t4_r4)
elseif (real_type .eq. 8) then
        call cpu_time(t3_r8)
        write(*,*) 'Integral (TRAP) =~', quadra(trap_real_8, f_real_8, a_r8, b_r8, n)
        call cpu_time(t4_r8)
elseif (real_type .eq. 10) then
        call cpu_time(t3_r10)
        write(*,*) 'Integral (TRAP) =~', quadra(trap_real_10, f_real_10, a_r10, b_r10, n)
        call cpu_time(t4_r10)
elseif (real_type .eq. 16) then
        call cpu_time(t3_r16)
        write(*,*) 'Integral (TRAP) =~', quadra(trap_real_16, f_real_16, a_r16, b_r16, n)
        call cpu_time(t4_r16)
else
        call cpu_time(t3_r4)
        write(*,*) 'Integral (TRAP) (REAL 4) =~', quadra(trap_real_4, f_real_4, a_r4, b_r4, n)
        call cpu_time(t4_r4)

        call cpu_time(t3_r8)
        write(*,*) 'Integral (TRAP) (REAL 8) =~', quadra(trap_real_8, f_real_8, a_r8, b_r8, n)
        call cpu_time(t4_r8)

        call cpu_time(t3_r10)
        write(*,*) 'Integral (TRAP) (REAL 10)=~', quadra(trap_real_10, f_real_10, a_r10, b_r10, n)
        call cpu_time(t4_r10)

        call cpu_time(t3_r16)
        write(*,*) 'Integral (TRAP) (REAL 16)=~', quadra(trap_real_16, f_real_16, a_r16, b_r16, n)
        call cpu_time(t4_r16)
endif

elseif (function_type .eq. 'sim') then

if (real_type .eq. 4) then
        call cpu_time(t5_r4)
        write(*,*) 'Integral (SIM) =~', quadra(sim_real_4, f_real_4, a_r4, b_r4, n)
        call cpu_time(t6_r4)
elseif (real_type .eq. 8) then
        call cpu_time(t5_r8)
        write(*,*) 'Integral (SIM) =~', quadra(sim_real_8, f_real_8, a_r8, b_r8, n)
        call cpu_time(t6_r8)
elseif (real_type .eq. 10) then
        call cpu_time(t5_r10)
        write(*,*) 'Integral (SIM) =~', quadra(sim_real_10, f_real_10, a_r10, b_r10, n)
        call cpu_time(t6_r10)
elseif (real_type .eq. 16) then
        call cpu_time(t5_r16)
        write(*,*) 'Integral (SIM) =~', quadra(sim_real_16, f_real_16, a_r16, b_r16, n)
        call cpu_time(t6_r16)
else
        call cpu_time(t5_r4)
        write(*,*) 'Integral (SIM) (REAL 4) =~', quadra(sim_real_4, f_real_4, a_r4, b_r4, n)
        call cpu_time(t6_r4)

        call cpu_time(t5_r8)
        write(*,*) 'Integral (SIM) (REAL 8) =~', quadra(sim_real_8, f_real_8, a_r8, b_r8, n)
        call cpu_time(t6_r8)

        call cpu_time(t5_r10)
        write(*,*) 'Integral (SIM) (REAL 10)=~', quadra(sim_real_10, f_real_10, a_r10, b_r10, n)
        call cpu_time(t6_r10)

        call cpu_time(t5_r16)
        write(*,*) 'Integral (SIM) (REAL 16)=~', quadra(sim_real_16, f_real_16, a_r16, b_r16, n)
        call cpu_time(t6_r16)
endif

elseif (function_type .ne. 'rectan') then

if (real_type .eq. 4) then
        call cpu_time(t3_r4)
        write(*,*) 'Integral (TRAP) =~', quadra(trap_real_4, f_real_4, a_r4, b_r4, n)
        call cpu_time(t4_r4)
        call cpu_time(t5_r4)
        write(*,*) 'Integral (SIM) =~', quadra(sim_real_4, f_real_4, a_r4, b_r4, n)
        call cpu_time(t6_r4)
elseif (real_type .eq. 8) then
        call cpu_time(t3_r8)
        write(*,*) 'Integral (TRAP) =~', quadra(trap_real_8, f_real_8, a_r8, b_r8, n)
        call cpu_time(t4_r8)
        call cpu_time(t5_r8)
        write(*,*) 'Integral (SIM) =~', quadra(sim_real_8, f_real_8, a_r8, b_r8, n)
        call cpu_time(t6_r8)
elseif (real_type .eq. 10) then
        call cpu_time(t3_r10)
        write(*,*) 'Integral (TRAP) =~', quadra(trap_real_10, f_real_10, a_r10, b_r10, n)
        call cpu_time(t4_r10)
        call cpu_time(t5_r10)
        write(*,*) 'Integral (SIM) =~', quadra(sim_real_10, f_real_10, a_r10, b_r10, n)
        call cpu_time(t6_r10)
elseif (real_type .eq. 16) then
        call cpu_time(t3_r16)
        write(*,*) 'Integral (TRAP)=~', quadra(trap_real_16, f_real_16, a_r16, b_r16, n)
        call cpu_time(t4_r16)
        call cpu_time(t5_r16)
        write(*,*) 'Integral (SIM) =~', quadra(sim_real_16, f_real_16, a_r16, b_r16, n)
        call cpu_time(t6_r16)
else
        call cpu_time(t3_r4)
        write(*,*) 'Integral (TRAP) (REAL 4) =~', quadra(trap_real_4, f_real_4, a_r4, b_r4, n)
        call cpu_time(t4_r4)
        call cpu_time(t3_r8)
        write(*,*) 'Integral (TRAP) (REAL 8) =~', quadra(trap_real_8, f_real_8, a_r8, b_r8, n)
        call cpu_time(t4_r8)
        call cpu_time(t3_r10)
        write(*,*) 'Integral (TRAP) (REAL 10)=~', quadra(trap_real_10, f_real_10, a_r10, b_r10, n)
        call cpu_time(t4_r10)
        call cpu_time(t3_r16)
        write(*,*) 'Integral (TRAP) (REAL 16)=~', quadra(trap_real_16, f_real_16, a_r16, b_r16, n)
        call cpu_time(t4_r16)

        call cpu_time(t5_r4)
        write(*,*) 'Integral (SIM) (REAL 4) =~', quadra(sim_real_4, f_real_4, a_r4, b_r4, n)
        call cpu_time(t6_r4)
        call cpu_time(t5_r8)
        write(*,*) 'Integral (SIM) (REAL 8) =~', quadra(sim_real_8, f_real_8, a_r8, b_r8, n)
        call cpu_time(t6_r8)
        call cpu_time(t5_r10)
        write(*,*) 'Integral (SIM) (REAL 10)=~', quadra(sim_real_10, f_real_10, a_r10, b_r10, n)
        call cpu_time(t6_r10)
        call cpu_time(t5_r16)
        write(*,*) 'Integral (SIM) (REAL 16)=~', quadra(sim_real_16, f_real_16, a_r16, b_r16, n)
        call cpu_time(t6_r16)
endif

endif

write(*,*) ' '

!-----------------------------------------------------------------------------

write(*,'(a,/)') ' Временные затраты (в единицах временных затрат rectan)'

write(*,*) '------------------------------------------------------------'
write(*,'(a,4x,a,5x,a,35x,a)') ' Время', ':', 'Задача N3', ':'
write(*,*) '------------------------------------------------------------'
write(*,'(a,6x,a,5x,a,8x,a,14x,a,9x,a)') ' Тип', ':', 'rectan', 'trap', 'sim', ':'
write(*,*) '------------------------------------------------------------'

if (function_type .eq. 'rectan') then

if (real_type .eq. 4) then
        write(*,'(a,2x,a,4x,i3,42x,a)') ' real(4)', ':', 1, ':'
elseif (real_type .eq. 8) then
        write(*,'(a,2x,a,4x,i3,42x,a)') ' real(8)', ':', 1, ':'
elseif (real_type .eq. 10) then
        write(*,'(a,1x,a,4x,i3,42x,a)') ' real(10)', ':', 1, ':'
elseif (real_type .eq. 16) then
        write(*,'(a,1x,a,4x,i3,42x,a)') ' real(16)', ':', 1, ':'
else
        write(*,'(a,2x,a,4x,i3,42x,a)') ' real(4)', ':', 1, ':'
        write(*,'(a,2x,a,4x,i3,42x,a)') ' real(8)', ':', 1, ':'
        write(*,'(a,1x,a,4x,i3,42x,a)') ' real(10)', ':', 1, ':'
        write(*,'(a,1x,a,4x,i3,42x,a)') ' real(16)', ':', 1, ':'
endif

elseif (function_type .eq. 'trap') then

if (real_type .eq. 4) then
        t_rectan_r4=t2_r4-t1_r4
        t_trap_r4=(t4_r4-t3_r4)/t_rectan_r4
        write(*,'(a,2x,a,4x,i3,6x,e15.7,21x,a)') ' real(4)', ':', 1, t_trap_r4, ':'
elseif (real_type .eq. 8) then
        t_rectan_r8=t2_r8-t1_r8
        t_trap_r8=(t4_r8-t3_r8)/t_rectan_r8
        write(*,'(a,2x,a,4x,i3,6x,e15.7,21x,a)') ' real(8)', ':', 1, t_trap_r8, ':'
elseif (real_type .eq. 10) then
        t_rectan_r10=t2_r10-t1_r10
        t_trap_r10=(t4_r10-t3_r10)/t_rectan_r10
        write(*,'(a,1x,a,4x,i3,6x,e15.7,21x,a)') ' real(10)', ':', 1, t_trap_r10, ':'
elseif (real_type .eq. 16) then
        t_rectan_r16=t2_r16-t1_r16
        t_trap_r16=(t4_r16-t3_r16)/t_rectan_r16
        write(*,'(a,1x,a,4x,i3,6x,e15.7,21x,a)') ' real(16)', ':', 1, t_trap_r16, ':'
else
        t_rectan_r4=t2_r4-t1_r4
        t_rectan_r8=t2_r8-t1_r8
        t_rectan_r10=t2_r10-t1_r10
        t_rectan_r16=t2_r16-t1_r16

        t_trap_r4=(t4_r4-t3_r4)/t_rectan_r4
        t_trap_r8=(t4_r8-t3_r8)/t_rectan_r8
        t_trap_r10=(t4_r10-t3_r10)/t_rectan_r10
        t_trap_r16=(t4_r16-t3_r16)/t_rectan_r16

        write(*,'(a,2x,a,4x,i3,6x,e15.7,21x,a)') ' real(4)', ':', 1, t_trap_r4, ':'
        write(*,'(a,2x,a,4x,i3,6x,e15.7,21x,a)') ' real(8)', ':', 1, t_trap_r8, ':'
        write(*,'(a,1x,a,4x,i3,6x,e15.7,21x,a)') ' real(10)', ':', 1, t_trap_r10, ':'
        write(*,'(a,1x,a,4x,i3,6x,e15.7,21x,a)') ' real(16)', ':', 1, t_trap_r16, ':'
endif

elseif (function_type .eq. 'sim') then

if (real_type .eq. 4) then
        t_rectan_r4=t2_r4-t1_r4
        t_sim_r4=(t6_r4-t5_r4)/t_rectan_r4
        write(*,'(a,2x,a,4x,i3,24x,e15.7,3x,a)') ' real(4)', ':', 1, t_sim_r4, ':'
elseif (real_type .eq. 8) then
        t_rectan_r8=t2_r8-t1_r8
        t_sim_r8=(t6_r8-t5_r8)/t_rectan_r8
        write(*,'(a,2x,a,4x,i3,24x,e15.7,3x,a)') ' real(8)', ':', 1, t_sim_r8, ':'
elseif (real_type .eq. 10) then
        t_rectan_r10=t2_r10-t1_r10
        t_sim_r10=(t6_r10-t5_r10)/t_rectan_r10
        write(*,'(a,1x,a,4x,i3,24x,e15.7,3x,a)') ' real(10)', ':', 1, t_sim_r10, ':'
elseif (real_type .eq. 16) then
        t_rectan_r16=t2_r16-t1_r16
        t_sim_r16=(t6_r16-t5_r16)/t_rectan_r16
        write(*,'(a,1x,a,4x,i3,24x,e15.7,3x,a)') ' real(16)', ':', 1, t_sim_r16, ':'
else
        t_rectan_r4=t2_r4-t1_r4
        t_rectan_r8=t2_r8-t1_r8
        t_rectan_r10=t2_r10-t1_r10
        t_rectan_r16=t2_r16-t1_r16

        t_sim_r4=(t6_r4-t5_r4)/t_rectan_r4
        t_sim_r8=(t6_r8-t5_r8)/t_rectan_r8
        t_sim_r10=(t6_r10-t5_r10)/t_rectan_r10
        t_sim_r16=(t6_r16-t5_r16)/t_rectan_r16

        write(*,'(a,2x,a,4x,i3,24x,e15.7,3x,a)') ' real(4)', ':', 1, t_sim_r4, ':'
        write(*,'(a,2x,a,4x,i3,24x,e15.7,3x,a)') ' real(8)', ':', 1, t_sim_r8, ':'
        write(*,'(a,1x,a,4x,i3,24x,e15.7,3x,a)') ' real(10)', ':', 1, t_sim_r10, ':'
        write(*,'(a,1x,a,4x,i3,24x,e15.7,3x,a)') ' real(16)', ':', 1, t_sim_r16, ':'
endif

else

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

endif

write(*,'(//)')

end
