program main
use poly
implicit none
integer a,r(0:31),i
read(*,*) a
write(*,*) 'a=',a
write(*,'(b32.32)') a
call tentotwo1(a,r)
write(*,'(32i1)') (r(i),i=31,0,-1)

end
