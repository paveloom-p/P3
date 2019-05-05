module subprograms
implicit none

real(8), parameter :: pi=4.D0*datan(1.D0)

contains

subroutine erfc_as(x,eps,res,k,ier)
implicit none

real(8), intent(in) :: x, eps
real(8), intent(out) :: res
integer, intent(out) :: k, ier

real(8) a, b, koef, xx

xx=x*x
koef=sqrt(pi)*x*exp(xx)

b=1d0/xx/2d0
res=1d0
a=res

k=1

do

if (b .lt. abs(res)*eps) then 
ier=0; exit
endif

a=b
res=res+(-1)**(k)*b

k=k+1
b=b/2d0/xx*(2d0*k-1d0)

if (abs(b) .gt. abs(a)) then
ier=2; exit
endif

enddo

res=res/koef

end subroutine erfc_as

end module subprograms
