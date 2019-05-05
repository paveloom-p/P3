module functions
implicit none
contains

function fibon(n) result(f3)
implicit none
integer f1, f2, f3, n, i 

f1=1d0
f2=1d0
f3=1d0
 
do i=3,n
f3=f1+f2
f1=f2
f2=f3
enddo

end function fibon

recursive integer function fibor1(n) result(w)
implicit none
integer n

if (n<3) then
w=1
else
w=fibor1(n-1)+fibor1(n-2)
endif

end function fibor1

recursive integer function fibor2(k,n,fk1,fk) result(w)
implicit none
integer n, k
integer fk1, fk

if (k.eq.n) then
w=fk
else
w=fibor2(k+1,n,fk, fk+fk1)
endif

end function fibor2

recursive integer function fibor3(n,f) result(w)
implicit none
integer n
integer f(100)

if ((n.eq.1).or.(n.eq.2)) then
f(n)=1; w=1
endif

if (f(n).ne.0) then
w=f(n)
else
f(n)=fibor3(n-1,f)+fibor3(n-2,f)
w=f(n)
endif

end function fibor3

real(8) function fibof(n) result(w)
implicit none
integer n
real(8), parameter :: c5=sqrt(5d0) 
real(8), parameter :: x=(1d0+c5)*0.5d0
real(8), parameter :: y=(1d0-c5)*0.5d0

w=(x**n-y**n)/c5

end function fibof

end module functions
