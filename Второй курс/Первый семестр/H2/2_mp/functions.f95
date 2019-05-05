module functions
implicit none
contains

function f(x) result(y)
use real_type
implicit none
real(mp) p, q, r
common / parameters / p, q, r 
real(mp) x, y                              
y=x+p+q+r                         
end function f

end module functions
