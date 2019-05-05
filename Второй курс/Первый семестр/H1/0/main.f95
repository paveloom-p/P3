program main
use quadra
implicit none

!-----------QUADRA---------------

integer ier, i

integer n                                  
real(8) a, b, h 
                          
real(8), allocatable :: codomain(:)              
real(8), allocatable :: domain(:)               
                           
real(8), allocatable :: codomain_t(:)            
real(8), allocatable :: domain_t(:)              
                           
real(8), allocatable :: codomain_s(:)            
real(8), allocatable :: domain_s(:)              

!-----------QUADRA---------------

write(*,'(/,a,/)') ' QUADRA '

!-----------------------------------------------------------------------------

write(*,'(a,/)') ' Формула средних прямоугольников '

read(*,*) a, b, n
write(*,*) 'Был задан промежуток: от ', a, 'до ', b, '.'
write(*,*) 'Будут рассмотрены первые ', n, 'значений(-я).'

h=(b-a)/n

allocate(domain(10000), stat=ier)
        if (ier/=0) then; write(*,*) 'ier=',ier,' do not allocate array!'
                  stop 1
        endif

allocate(codomain(size(domain)), stat=ier)
        if (ier/=0) then; write(*,*) 'ier=',ier,' do not allocate array!'
                  stop 1
        endif

domain=0; codomain=0

do i=1, n, 1
        domain(i)=a+(i-0.5)*h
        codomain(i)=f(domain(i))
enddo

write(*,*) 'Таблица подинтегральной функции выведена в файл results.dat'
open (unit=101, file='results.dat')
do i=1, n, 1
        write(101, '(e15.7,e15.7)') domain(i), codomain(i)
enddo
close(101)

write(*,*) 'Rectan:'
write(*,*) 'Integral=~', rectan(codomain, a, b, n) 

deallocate(domain, stat=ier)
deallocate(codomain, stat=ier)
write(*,'(//)')

!-----------------------------------------------------------------------------

write(*,'(/,a,/)') ' Формула трапеций '

write(*,*) 'Был задан промежуток: от ', a, 'до ', b, '.'
write(*,*) 'Будут рассмотрены первые ', n+1, 'значений(-я).'

h=(b-a)/n

allocate(domain_t(10000), stat=ier)
        if (ier/=0) then; write(*,*) 'ier=',ier,' do not allocate array!'
                  stop 1
        endif

allocate(codomain_t(size(domain_t)), stat=ier)
        if (ier/=0) then; write(*,*) 'ier=',ier,' do not allocate array!'
                  stop 1
        endif


domain_t=0; codomain_t=0

do i=1, n+1, 1
        domain_t(i)=a+(i-1)*h
        codomain_t(i)=f(domain_t(i))
enddo

write(*,*) 'Таблица подинтегральной функции выведена в файл results_t.dat'
open (unit=102, file='results_t.dat')
do i=1, n+1, 1
        write(102, '(e15.7,e15.7)') domain_t(i), codomain_t(i)
enddo
close(102)

write(*,*) 'Trap:'
write(*,*) 'Integral=~', trap(codomain_t, a, b, n) 

deallocate(domain_t, stat=ier)
deallocate(codomain_t, stat=ier)
write(*,'(//)')

!-----------------------------------------------------------------------------

write(*,'(/,a,/)') ' Формула Симпсона (парабол) '

write(*,*) 'Был задан промежуток: от ', a, 'до ', b, '.'
write(*,*) 'Будут рассмотрены первые ', n+1, 'значений(-я).'

h=(b-a)/n

allocate(domain_s(10000), stat=ier)
        if (ier/=0) then; write(*,*) 'ier=',ier,' do not allocate array!'
                  stop 1
        endif

allocate(codomain_s(size(domain_s)), stat=ier)
        if (ier/=0) then; write(*,*) 'ier=',ier,' do not allocate array!'
                  stop 1
        endif


domain_s=0; codomain_s=0

do i=1, n+1, 1
        domain_s(i)=a+(i-1)*h
        codomain_s(i)=f(domain_s(i))
enddo

write(*,*) 'Таблица подинтегральной функции выведена в файл results_s.dat'
open (unit=103, file='results_s.dat')
do i=1, n+1, 1
        write(103, '(e15.7,e15.7)') domain_s(i), codomain_s(i)
enddo
close(103)

write(*,*) 'Sim:'
write(*,*) 'Integral=~', sim(codomain_s, a, b, n)

deallocate(domain_s, stat=ier)
deallocate(codomain_s, stat=ier)
write(*,'(//)')

end
