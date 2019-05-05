module subprograms
implicit none

contains

subroutine E1_as(x,eps,res,k,ier)
implicit none

real(8), intent(in) :: x, eps
real(8), intent(out) :: res
integer, intent(out) :: k, ier

real(8) a, b, koef

koef=exp(-x)/x

b=1d0/x
res=1d0
a=res

k=1

do while (1 .gt. 0)

if (b .lt. abs(res)*eps) then 
ier=0; exit
endif

a=b
res=res+(-1)**(k)*b

k=k+1
b=b*k/x

if (abs(b) .gt. abs(a)) then
ier=2; exit
endif

enddo

res=koef*res

end subroutine E1_as

subroutine E1_my(x,eps,res,k)
implicit none

real(8), intent(in) :: x, eps
real(8), intent(out) :: res
integer, intent(out) :: k

real(8) b

b=x*x/4
res=x

k=2

do while (b .gt. abs(res)*eps)

res=res+(-1)**(k+1)*b

k=k+1
b=b*x*(k-1)/(k*k)

enddo

end subroutine E1_my

end module subprograms
