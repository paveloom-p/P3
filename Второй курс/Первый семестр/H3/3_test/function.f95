module function
implicit none
contains

character(len(d)) function reverse(d)
character(*) d
character, dimension(len(d)) :: c

c=transfer(d,c)

c=c(len(d):1:-1)

reverse=transfer(c,reverse)

end function reverse

character(len(d)) function reverse_2(d)
character(*) d

reverse_2=d

end function reverse_2

end module
