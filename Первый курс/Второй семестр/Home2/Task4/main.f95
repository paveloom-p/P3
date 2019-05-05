program program1
use TypeCoord
use procedures
implicit none
type (coord4) A4, B4
type (coord8) A8, B8
real(4) d4
real(8) d8

write(*,*) 'Enter coordinates of A:'
call rdr(A4)
call rdr(A8)
write(*,*) 'Coordinates of A:'
call wrt(A4)
call wrt(A8)

write(*,*) 'Enter coordinates of B:'
call rdr(B4)
call rdr(B8)
write(*,*) 'Coordinates of B:'
call wrt(B4)
call wrt(B8)

write(*,*) 'Distance between A and B:'
call dist(A4,B4,d4)
write(*,*) d4
call dist(A8,B8,d8)
write(*,*) d8
end





 







