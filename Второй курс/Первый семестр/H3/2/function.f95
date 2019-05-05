module function
implicit none
contains

function digit_to_int(b) result(a)
integer a
character b

a=ichar(b)-ichar('0')

end function digit_to_int

end module
