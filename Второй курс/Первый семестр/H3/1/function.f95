module function
implicit none
contains

function int_to_digit(a) result(b)
integer a
character b

b=achar(ichar('0')+a)

end function int_to_digit

end module
