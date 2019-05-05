program main
use poly
implicit none
real t1,t2,t12
integer n,i,a(10000),res(-16:16)
read(*,*) n
do i=1,n
read(*,*) a(i)
enddo
write(*,*) n
write(*,*) (i,a(i),i=1,n)
res=0;
call second(t1)
call pol3(a,n,res)
call second(t2)
write(*,*) 'Time',t2-t1
write(*,'(2i3)') (i,res(i),i=-16,16)
end
