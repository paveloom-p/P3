program main
use quadr
use functions
implicit none

!-----------QUADRA---------------

integer n, i, ier                           
real(8) a, b, h

real(8), allocatable :: domain(:) 

real(8), allocatable :: domain_r(:) 

!-----------QUADRA---------------

write(*,'(/,a,/)') ' QUADRA '

!-----------------------------------------------------------------------------

read(*,*) a, b, n

h=(b-a)/n

allocate(domain_r(10000), stat=ier)
        if (ier/=0) then; write(*,*) 'ier=',ier,' do not allocate array!'
                  stop 1
        endif

allocate(domain(10000), stat=ier)
        if (ier/=0) then; write(*,*) 'ier=',ier,' do not allocate array!'
                  stop 1
        endif

domain=0; domain_r=0;

do i=1, n+1, 1
        domain_r(i)=a+(i-0.5)*h

        domain(i)=a+(i-1)*h
enddo


write(*,*) 'Был задан промежуток: от ', a, 'до ', b, '.'
write(*,'(a,i3,2x,a,i3,2x,a)') ' Будут рассмотрены первые ', n, &
'значений(-я) (', n+1, ', если trap или sim).'
write(*,*) 'Integral=~', quadra(sim_function, f, a, b, n)
!write(*,*) 'Integral=~', quadra(rectan_array, domain_r, a, b, n)
write(*,*) 'Integral=~', quadra(trap_array, domain, a, b, n)

write(*,*) ' '


deallocate(domain, stat=ier)
deallocate(domain_r, stat=ier)

end
