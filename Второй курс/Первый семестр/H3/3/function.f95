module function
implicit none
contains

function reverse(d) result(f)
character(*) d
character, dimension(len(d)) :: c
character(len(d)) f

c=transfer(d,c)

c=c(len(d):1:-1)

f=transfer(c,f)

end function reverse

end module
