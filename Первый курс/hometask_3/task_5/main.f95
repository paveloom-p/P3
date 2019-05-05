program main
use binrad
implicit none
integer r,a(0:31),i
read(*,'(32i1)') (a(i),i=31,0,-1)
write(*,'(a,32i1)') 'a=',(a(i),i=31,0,-1)
call tentotwo1(a,r)
write(*,'(a,i5)') 'a=',r
end
