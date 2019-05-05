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
function sim(y,a,b,n) result(s)
    implicit none
       real(8) a,b,s,y(n+1),s1,s2,h,sp
       integer i,n
s1=0; s2=0; h=(b-a)/n
do i=3,n-1,2
    s2=s2+y(i)
enddo
do i=2,n,2
    s1=s1+y(i)
enddo
sp=y(1)+y(n+1)+4*s1+2*s2
s=(h/3)*sp
end function
function f(x) result(w)
implicit none
real(8) w,x
w=x*x*x
end function
end module