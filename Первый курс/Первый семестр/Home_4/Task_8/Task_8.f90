program h4_8
implicit none
real(8) x0, x1, eps
integer i
write(*,*) 'Укажите число и точность измерений.'
read(*,*) x0, eps

do while (abs(x0-x1) .LT. eps)
x1=2*x0+
enddo

end
