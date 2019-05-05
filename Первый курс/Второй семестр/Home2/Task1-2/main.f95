program program1
use TypeCoord
use procedures
implicit none
type (coord4) A, B
real d

call rdr4(A)
write(*,*) 'Coordinates of A:'
call wrt4(A)

call rdr4(B)
write(*,*) 'Coordinates of B:'
call wrt4(B)

write(*,*) 'Distance between A and B:'
call dist4(A,B,d)
write(*,*) d
end





 







