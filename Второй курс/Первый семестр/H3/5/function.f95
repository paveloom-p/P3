module function
implicit none

integer c, i

character digit(0:31)  / '0','1','2','3','4','5','6','7','8','9',&
&'A','B','C','D','E','F','G','H','I','J',&
&'K','L','M','N','O','P','Q','R','S','T',&
&'U','V' /

contains

function z_hex(a,z) result(f)
integer a, z
character(50) f

f=' '
i=len(f)

do while (a .gt. z)
c=mod(a,z)
f(i:i)=digit(c)
a=a/z
i=i-1

enddo

f(i:i)=digit(a)

f=adjustl(f)

end function z_hex

end module
