module sort
implicit none
contains
function shift(a,n,ir) result(w)
integer n,i,ir
real a(n),k,w(n),l
 if (ir.eq.0) then
    l=a(n)
    do i=n,2,-1
       k=a(i-1)
       a(i)=k
    enddo
    a(1)=l
               else
    l=a(1)
    do i=1,n-1
       a(i)=a(i+1)
    enddo
    a(n)=l
               endif
w=a
end function
function invert(a,n) result(w)
integer n,i
real a(n),w(n)
  do i=1,n
    w(i)=a(n-i+1)
  enddo
end function
recursive function invertr(a,k,n) result(w)
integer n,k
real a(n),w(n)
 if (k<=n) then
    w(k)=a(n-k+1)
    k=k+1
    w=invertr(a,k,n)
 endif
end function
subroutine unite(a,b,k,n,w)
integer k,n,i,j
real a(n),b(k),w(k+n)
   if (n.eq.k) then
do i=1,k+n
    if (mod(i,2).eq.0) then
       w(i)=b(i/2)
                       else
       w(i)=a((i+1)/2)
                       endif
enddo
             endif
    if (n<k) then
i=1;
do while (i<=n)
  w(2*i-1)=a(i); w(2*i)=b(i); i=i+1;
enddo
  j=2*(i-1)+1;
do while (i<=k)
  w(j)=b(i); i=i+1; j=j+1
enddo
             endif
    if (n>k) then
i=1;
do while (i<=k)
   w(2*i-1)=a(i); w(2*i)=b(i); i=i+1;
enddo
   j=2*(i-1)+1;
do while (i<=n)
  w(j)=a(i); i=i+1; j=j+1;
enddo
             endif
end subroutine
function chet(a,n) result(w)
integer a(n),w,i,n
w=0;
do i=1,n
    if (mod(a(i),2).eq.1) then
       w=w+1
    endif
enddo
end function
end module
