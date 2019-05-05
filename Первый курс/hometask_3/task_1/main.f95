program main
use poly
implicit none
real p(0:10000),res,x1
real t1,t2,t12
integer n,i
write(*,*) 'Введите степень полинома'
read(*,*) n
write(*,*) 'Введите коэффиценты от нулевой степени до n'
read(*,*) (p(i),i=0,n)
write(*,*) 'Введите аргумент полинома'
read(*,*) x1
write(*,*) 'Значение коэффицентов соответствующих степеней'
write(*,'(i3,e15.7)') (i,p(i),i=0,n)
call second(t1)
do i=1,10000
call pol(x1,p,n,res)
enddo
call second(t2)
write(*,1000)
write(*,1001) t2-t1,x1,res
1000 format(8x,'Time',14x,'x',18x,'P(x)')
1001 format(1x,e15.7,2x,e15.7,5x,e15.7,6x,e15.7)
end