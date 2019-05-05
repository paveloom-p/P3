program main
implicit none

external f, quadra
external trap, sim, rectan

!-----------QUADRA---------------

real(8) f, trap, sim, rectan, quadra

integer n                                  
real(8) a, b                                      

!-----------QUADRA---------------

write(*,'(/,a,/)') ' QUADRA '

!-----------------------------------------------------------------------------

read(*,*) a, b, n

write(*,*) 'Был задан промежуток: от ', a, 'до ', b, '.'
write(*,'(a,i3,2x,a,i3,2x,a)') ' Будут рассмотрены первые ', n, &
'значений(-я) (', n+1, ', если trap или sim).'
write(*,*) 'Integral=~', quadra(sim, f, a, b, n)

write(*,*) ' '

end
