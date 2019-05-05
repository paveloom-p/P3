program test
implicit none

real(8) A, B
integer i

open(100,file='data_decart'); open(101,file='data_decart_out')

do i=1, 7, 1
read(100,*) B, A
write(101,*) sin(90-A)*cos(B), sin(90-A)*sin(B), cos(90-A)
enddo

close(101); close(100)

end program test
