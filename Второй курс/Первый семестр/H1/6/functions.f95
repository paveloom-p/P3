module functions
implicit none
contains

function f(x) result(y)     
real(8) x, y                              
y=x**3                                 
end function f

end module functions
