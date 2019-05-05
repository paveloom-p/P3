module function
implicit none
contains

function hex(a) result(c)
integer a
character(10) c

write(c,'(z10)') a

end function hex

end module
