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
function f(x) result(w)
implicit none
real(8) w,x
w=3*x*x-1
end function
end module