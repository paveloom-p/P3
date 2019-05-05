module binrad
implicit none
contains
subroutine tentotwo2(a,w)
implicit none
    integer w(0:10),a(0:31),i,k
w=0
k=0
!if (a(31).eq.1) then
!  i=0
!    if (a(i).eq.0) then
!       a(i)=1; k=1;
!       do while (k.ne.0)
!          i=i+1
!          if (a(i).eq.1) then
!             a(i)=0
!             k=0
!          else
!             a(i)=1
!          endif
!       enddo
!    else
!          a(i)=0
!    endif
! do i=0,31
!  if (a(i).eq.0) then
!a(i)=1;
! else
!a(i)=0
!  endif
! enddo
!k=2;
!endif
do i=0,9
w(i)=a(3*i)*1+a(3*i+1)*2+a(3*i+2)*4
enddo
w(10)=a(30)*1+a(31)*2
end subroutine
subroutine tentotwo1(a,w)
implicit none
    integer w,a(0:31),i,k,r
w=0
k=0
if (a(31).eq.1) then
  i=0
    if (a(i).eq.0) then
       a(i)=1; k=1;
       do while (k.ne.0)
          i=i+1
          if (a(i).eq.1) then
             a(i)=0
             k=0
          else
             a(i)=1
          endif
       enddo
    else
          a(i)=0
    endif
 do i=0,31
  if (a(i).eq.0) then
a(i)=1;
  else
a(i)=0
  endif
 enddo
k=2;
endif
r=1
do i=0,31
w=w+a(i)*r
r=r*2
enddo
if (k.eq.2) then; w=-w; endif
end subroutine
end module
