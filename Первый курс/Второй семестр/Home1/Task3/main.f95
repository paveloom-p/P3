program tstfibo4
use functions
implicit none

integer, parameter :: krep=10000000
integer i, n, fn, fr1, fr2, fr3, fk1, fk, j
real t1,t2,t3,t4,t5,t6, t12, t23, t34, t45, t56
real*8 ff
integer work(100)

!interface

!	integer function fibon(n)
!	integer n 
!	end function fibon
	
!	recursive integer function fibor1(n) result(w)
!	integer n
!	end function fibor1

!	recursive integer function fibor2(k,n,fk1,fk) result(w)
!	integer k, n, fk1, fk
!	end function fibor2
	
!	recursive integer function fibor3(n,f) result(w)
!	integer n
!	integer f(100)
!	end function fibor3

!	real(8) function fibof(n)
!	integer n
!	end function fibof

!end interface

work=0
write(*,*) ' Число нерекурсивных вызовов (krep)=', krep
write(*,*) ' Введите n:'
read(*,*) n
write(*,*) ' n=',n

call second(t1)
do i=1,krep
fn=fibon(n)
enddo

call second(t2)
do i=1,krep
fr1=fibor1(n)
enddo

fk1=0d0
fk=1d0
call second(t3)
do i=1,krep
fr2=fibor2(1,n,fk1,fk)
enddo

call second(t4)
do i=1,krep
fr3=fibor3(n,work)
enddo

call second(t5)
do i=1,krep
ff=fibof(n)
enddo

call second(t6)
write(*,1000)
write(*,1001) fn, fr1, fr2, fr3, ff

t12=t2-t1
t23=t3-t2
t34=t4-t3
t45=t5-t4
t56=t6-t5

write(*,1002) t12, t23, t34, t45, t56
write(*,*) ' '

1000 format( 1x,'Функция',5x,'fibon',7x,'fibor1',7x,'fibor2',&
7x,'fibor3',7x,'fiborf')
1001 format(1x,'Значение',i10,2x,i10,2x,i10,2x,i10,6x,e11.4)
1002 format(1x,'Время',3x,e11.4,2x,e11.4,2x,e11.4,2x,e11.4,2x,e11.4)

end
