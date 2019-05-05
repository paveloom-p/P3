program main
use binrad
implicit none
integer r(0:10),a(0:31),i
read(*,'(32i1)') (a(i),i=31,0,-1)
write(*,'(a,32i1)') 'a=',(a(i),i=31,0,-1)
call tentotwo2(a,r)
write(*,'(a,11i1)') 'a=',(r(i),i=10,0,-1)
end
