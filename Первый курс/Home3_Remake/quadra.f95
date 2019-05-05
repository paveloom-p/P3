module quadra
implicit none
contains

!-----------

function your_function(x) result(y)
real(8) x, y
y=x
end function your_function

!-----------

function rectan(y,a,b,n) result(integral)
real(8) y(:), a, b, h, alpha, summa, integral
!real(8) rectan
integer n, i

h=(b-a)/n
alpha=a+h/2

summa=0
do i=1, n,1
	summa=summa+y(i)
enddo

integral=h*summa

!write(*,*) 'Summa=',summa

end function rectan

!-----------

end module quadra
