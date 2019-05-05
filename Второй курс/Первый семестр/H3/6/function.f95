module function
implicit none
contains

function summ(a,b) result(c)
character(50) a, b, c
integer(16) ia, ib, ic

read(a,*) ia
read(b,*) ib

ic=ia+ib

write(c,*) ic

end function summ

end module
