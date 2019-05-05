module poly
implicit none
contains
subroutine tentotwo1(a,w)
implicit none
    integer a,w(0:31),i,k
w(0:31)=0
if (a<0) then
w(31)=1
  do i=0,30
    if (mod(a,2).eq.0) then
w(i)=1
                       else
w(i)=0
                       endif
a=a/2
  enddo
  i=0
    if (w(i).eq.1) then
       w(i)=0; k=1;
       do while (k.ne.0)
          i=i+1
          if (w(i).eq.0) then; w(i)=1
             k=0
          else; w(i)=0
          endif
          enddo
          else
          w(i)=1
          endif
endif
if (a>0) then
  do i=0,30
w(i)=mod(a,2)
a=a/2
  enddo
endif
end subroutine
subroutine pol3(a,n,r)
integer n,i
integer a(n),r(-16:16)
r=0
do i=1,n
r(a(i))=r(a(i))+1;
enddo
end subroutine
subroutine pol1(x,p,n,w3,w4)
integer n,i
real p(0:n),x,w3,w4,w1,w2,x2
i=2;w1=p(0); x2=x*x;
 if (n.eq.0) then; w3=p(0); w4=p(0); return; endif
 if (n.eq.1) then; w3=p(1)+p(0)*x; w4=p(1)-p(0)*x; return; endif
 if (n>1) then
    do while (i<=n)
             w1=w1*x2+p(i)
             i=i+2
    enddo
 w2=p(1); i=3;
    do while (i<=n)
             w2=w2*x2+p(i)
             i=i+2
    enddo
    if (mod(n,2).eq.0) then; w3=w1+x*w2; w4=w1-x*w2
                       else; w3=w1*x+w2; w4=-w1*x+w2
    endif
    endif
end subroutine
subroutine pol(x,p,n,w3)
integer n,i
real p(0:n),x,w3,w1,w2,x2
i=2;w1=p(0); x2=x*x;
 if (n.eq.0) then; w3=p(0); return; endif
 if (n.eq.1) then; w3=p(1)+p(0)*x; return; endif
 if (n>1) then
    do while (i<=n)
             w1=w1*x2+p(i)
             i=i+2
    enddo
 w2=p(1); i=3;
    do while (i<=n)
             w2=w2*x2+p(i)
             i=i+2
    enddo
    if (mod(n,2).eq.0) then; w3=w1+x*w2;
                       else; w3=w1*x+w2;
    endif
    endif
end subroutine
end module
