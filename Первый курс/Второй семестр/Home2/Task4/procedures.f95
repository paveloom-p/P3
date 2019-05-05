module procedures
use TypeCoord
implicit none

interface rdr
module procedure rdr4, rdr8
end interface rdr

interface wrt
module procedure wrt4, wrt8
end interface wrt

interface dist
module procedure dist4, dist8
end interface dist

contains

subroutine rdr4(p)
 type (coord4) p
 read(*,*) p%x, p%y, p%z
end subroutine rdr4

subroutine rdr8(p)
 type (coord8) p
 read(*,*) p%x, p%y, p%z
end subroutine rdr8

subroutine wrt4(p)
 type (coord4) P
 write (*,*) p%x, p%y, p%z
end subroutine wrt4
 
subroutine wrt8(p)
 type (coord8) P
 write (*,*) p%x, p%y, p%z
end subroutine wrt8

subroutine dist4(a,b,d4)
 type (coord4) a, b
 real(4) d4
 d4=sqrt( ((a%x-b%x)*(a%x-b%x)) + ((a%y-b%y)*(a%y-b%y)) + ((a%z-b%z)*(a%z-b%z)) )
 return
end subroutine dist4

subroutine dist8(a,b,d8)
 type (coord8) a, b
 real(8) d8
 d8=sqrt( ((a%x-b%x)*(a%x-b%x)) + ((a%y-b%y)*(a%y-b%y)) + ((a%z-b%z)*(a%z-b%z)) )
 return
end subroutine dist8

end module procedures


 

 





