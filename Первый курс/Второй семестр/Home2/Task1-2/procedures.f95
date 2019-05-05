module procedures
use TypeCoord
implicit none
contains

subroutine rdr4(p)
 type (coord4) p
 read(*,*) p%x, p%y, p%z
end subroutine rdr4
 
subroutine wrt4(p)
 type (coord4) P
 write (*,*) p%x, p%y, p%z
end subroutine wrt4

subroutine dist4(a,b,d)
 type (coord4) a, b
 real(4) d
 d=sqrt( ((a%x-b%x)*(a%x-b%x)) + ((a%y-b%y)*(a%y-b%y)) + ((a%z-b%z)*(a%z-b%z)) )
 return
end subroutine dist4
end module procedures


 

 





