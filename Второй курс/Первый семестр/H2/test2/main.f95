program test
use my_prec
implicit none

interface
function wr(x) result(y); use my_prec; real(mp) x,y; end function wr
end interface

real(mp) x    !, wr
read(*,*) x
write(*,*) ' x=',x
write(*,*) wr(x)

end
