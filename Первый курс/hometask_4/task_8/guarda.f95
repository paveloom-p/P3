module guarda
implicit none
contains
function rectan(y,a,b,n) result(s)
implicit none
    real(8) y(n),a,b,h,summ,s
    integer i,n
h=(b-a)/n; summ=0;
do i=1,n
summ=summ+y(i);
enddo
s=summ*h
end function
function trap(y,a,b,n) result(s)
implicit none
    real(8) y(n),a,b,h,s
    integer n,i
h=(b-a)/n;
s=((f(a)+f(b))/2+sum(y(2:n)))*h
end function
function f(x) result(w)
implicit none
real(8) w,x
w=2*x-1
end function
end module