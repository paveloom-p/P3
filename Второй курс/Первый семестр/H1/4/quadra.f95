module quadr
implicit none
contains

function quadra(fqua,f,a,b,n) result(i)
real(8) fqua, i
real(8), external :: f
integer n                                  
real(8) a, b

i=fqua(f,a,b,n)

end function quadra

function rectan(f,a,b,n) result(integral)
real(8) f, a, b, h, x, summa, integral
integer n, i       

h=(b-a)/n

summa=0
do i=1, n, 1
        x=a+(i-0.5)*h
        summa=summa+f(x)
enddo

integral=h*summa

end function rectan

!-----------

function trap(f,a,b,n) result(integral)
real(8) f, a, b, h, x, summa, integral
integer n, i

h=(b-a)/n

summa=(f(a)+f(b))/2
do i=2, n, 1
        x=a+(i-1)*h
        summa=summa+f(x)
enddo

integral=h*summa

end function trap

!-----------

function sim(f,a,b,n) result(integral)
real(8) f, a, b, h, x, summa_ch, summa_nech, summa, integral
integer n, i

h=(b-a)/n

summa_ch=0
do i=2, n, 2
        x=a+(i-1)*h
        summa_ch=summa_ch+f(x)
enddo

summa_nech=0
do i=3, n-1, 2
        x=a+(i-1)*h
        summa_nech=summa_nech+f(x)
enddo

summa=f(a)+f(b)+4*summa_ch+2*summa_nech

integral=(h/3)*summa

end function sim

end module quadr
