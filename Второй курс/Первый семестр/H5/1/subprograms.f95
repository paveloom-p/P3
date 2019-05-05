module subprograms
implicit none

contains

subroutine E1_my(x,eps,res,umax,k)
implicit none

real(8), intent(in) :: x, eps
real(8), intent(out) :: res, umax
integer, intent(out) :: k

real(8) b

umax=x*x/4
b=umax
res=x

if (x .gt. umax) umax=x

k=2

do while (b .gt. abs(res)*eps)

res=res+(-1)**(k+1)*b

k=k+1
b=b*x*(k-1)/(k*k)

if (b .gt. umax) umax=b

enddo

end subroutine E1_my

end module subprograms
