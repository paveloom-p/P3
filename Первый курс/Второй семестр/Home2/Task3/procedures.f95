module procedures
use TypeCoord
implicit none
contains

subroutine rdr8(p)
 type (coord8) p
 read(*,*) p%x, p%y, p%z
end subroutine rdr8
 
subroutine wrt8(p)
 type (coord8) P
 write (*,*) p%x, p%y, p%z
end subroutine wrt8

subroutine dist8(a,b,d)
 type (coord8) a, b
 real(8) d
 d=sqrt( ((a%x-b%x)*(a%x-b%x)) + ((a%y-b%y)*(a%y-b%y)) + ((a%z-b%z)*(a%z-b%z)) )
 return
end subroutine dist8
end module procedures


 

 





