module sort
implicit none
contains
function search(a,n,o) result(i)
real a(n),o,a1,b1,c1
integer i1,n,c,i2,i
i1=1; i2=n; i=0;
c=n/2; a1=a(i1); b1=a(i2); c1=a(c);
do while (abs(i2-i1)>=0)
    if ((o<c1).or.(o>c1)) then
            if (o<c1) then
               i2=c; b1=a(i2); c=i2/2; c1=a(c)
                      else
               i1=c; a1=a(i1); c=(i2+i1)/2; c1=a(c)
                      endif
    endif
    if ((o.eq.a1).or.(o.eq.b1)) then
       if (o.eq.a1) then
          i=i1; exit
                    else
          i=i2; exit
                    endif
    endif
    if (o.eq.c1) then
       i=c; exit
    endif
    if ((i2-i1).eq.2) then
       if (a(i2-1).eq.o) then
          i=i2-1; exit
                         else
          exit
       endif
    endif
enddo
end function
recursive subroutine searchr(a1,b1,a,o,i)
implicit none
real a(:),o,ac
integer i,n,c,a1,b1
n=size(a)
if (a(a1).eq.o) then
    i=a1; return
endif
if (a(b1).eq.o) then
    i=b1; return
endif
if (a(b1)-a(a1)>=2) then
    i=(a1+b1)/2; c=i; ac=a(c)
    if (ac.eq.o) then
    return
    else
    c=(a1+b1)/2; ac=a(c)
    if (o>ac) then
       call searchr(c,b1,a,o,i)
    else
       call searchr(a1,c,a,o,i)
    endif
    endif
    endif
end subroutine
function qsort(a,n) result(ier)
integer n,ier,i,k,l
real a(n)
k=0; ier=4; l=0;
  do i=1,n-1
     if (a(i)<a(i+1)) then
        k=k+1
                     else
        if (a(i)<=a(i+1)) then
           l=l+1;
                         endif
                     endif
  enddo
    if (k.eq.n-1) then
       ier=0;
                  endif
    if (((l+k).eq.(n-1)).and.(l.ne.0)) then
       ier=2;
                    endif
k=0; l=0;
  do i=1,n-1
    if (a(i)>a(i+1)) then
       k=k+1
                    else
       if (a(i)>=a(i+1)) then
          l=l+1
                        endif
                     endif
    enddo
     if (k.eq.n-1) then
        ier=1;
                   endif
     if (((k+l).eq.n-1).and.(l.ne.0)) then
        ier=3;
                       endif
end function
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
