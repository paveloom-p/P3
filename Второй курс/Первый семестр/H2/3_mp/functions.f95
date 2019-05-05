module functions
use real_type
implicit none

!ПАРАМЕТРЫ

!Пусть на подинтегральную функцию действуют смещения p, q и r;

real(mp) p, q, r

contains

function f(x) result(y)
use real_type
implicit none
real(mp) x, y                              
y=x+p+q+r                         
end function f

end module functions
