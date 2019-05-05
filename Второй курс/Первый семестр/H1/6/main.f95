program main
use quadr
use functions
implicit none

!-----------QUADRA---------------

integer n, ier, k                                
real(8) a, b, eps                                    

!-----------QUADRA---------------

write(*,'(/,a,/)') ' QUADRA '

!-----------------------------------------------------------------------------

read(*,*) a, b, k, eps

ier=0

write(*,*) 'Был задан промежуток: от ', a, 'до ', b, '.'
write(*,'(a,i3)') '   k   = ', k
write(*,'(a,e16.8,/)') ' eps   =', eps
write(*,*) 'Значение интеграла:', quadrauto(trap,f,a,b,n,eps,k,ier)

write(*,*) ' '

if (ier .eq. 1) then
write(*,'(a,/,a,/,a,/)') '        Необходимая точность не была достигнута.', &
                                   '        Попробуйте увеличить k (число удвоений),', &
                                   '        либо понизить требуемую точность eps.'
else
write(*,'(a,i3,/)') ' Найдено решение при k   = ', k
endif

end
