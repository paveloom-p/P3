program tsthanoi
use thanoi
implicit none
integer n

write(*,*) 'Введите число колец:'
read (*,*) n; 
write(*,*) ' Число колец =', n
write(*,*) ' '

call hanoi(n,1,3)
write(*,*) ' '

end
