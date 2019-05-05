program main
use cycle_m
implicit none

integer a(6)

read(*,*) a
write(*,*) 'Был введен массив a:', a

write(*,*) bin_max(a)



end program
