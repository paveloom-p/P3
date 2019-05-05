module subprograms
implicit none

contains

subroutine erf_my(x,eps,res,umax,k)
implicit none

real(8), intent(in) :: x, eps
real(8), intent(out) :: res, umax
integer, intent(out) :: k

real(8) b, koef, xx, pi

pi=4.D0*datan(1.D0)

koef=2/sqrt(pi)
xx=x*x

umax=x*xx/3
b=umax
res=x

if (x .gt. umax) umax=x

k=1

do while (b .gt. abs(res)*eps)

res=res+(-1)**(k)*b

k=k+1
b=b*xx*(2*(k-1)+1)/(k*(2*k+1))

if (b .gt. umax) umax=b

enddo

res=koef*res

end subroutine erf_my

end module subprograms
