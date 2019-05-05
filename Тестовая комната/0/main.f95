program test
implicit none

character(1) a(6)
character b(6)
integer i

read(*,*) a

do i=1, 6, 1
write(*,*) ichar(a(i)), achar(ichar(a(i)))
enddo

write(*,*) a

write(*,*) transfer(a,b)

end program test
