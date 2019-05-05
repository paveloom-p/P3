program program1
use TypeCoord
use procedures
implicit none
type (coord8) A, B
real(8) d

write(*,*) 'Enter coordinates of A:'
call rdr8(A)
write(*,*) 'Coordinates of A:'
call wrt8(A)

write(*,*) 'Enter coordinates of B:'
call rdr8(B)
write(*,*) 'Coordinates of B:'
call wrt8(B)

write(*,*) 'Distance between A and B:'
call dist8(A,B,d)
write(*,*) d
end





 







