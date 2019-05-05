module quadr
implicit none
contains

!QUADRA_FUNCTION

function quadra_real_4(fqua,f,a,b,n) result(i)
implicit none
real(4) fqua, i
real(4), external :: f
integer n                                  
real(4) a, b

i=fqua(f,a,b,n)

end function quadra_real_4

function quadra_real_8(fqua,f,a,b,n) result(i)
implicit none
real(8) fqua, i
real(8), external :: f
integer n                                  
real(8) a, b

i=fqua(f,a,b,n)

end function quadra_real_8

function quadra_real_10(fqua,f,a,b,n) result(i)
implicit none
real(10) fqua, i
real(10), external :: f
integer n                                  
real(10) a, b

i=fqua(f,a,b,n)

end function quadra_real_10

function quadra_real_16(fqua,f,a,b,n) result(i)
implicit none
real(16) fqua, i
real(16), external :: f
integer n                                  
real(16) a, b

i=fqua(f,a,b,n)

end function quadra_real_16

!-----------

!RECTAN

function rectan_real_4(f,a,b,n) result(integral)
implicit none
real(4) f, a, b, h, x, summa, integral
integer n, i       

h=(b-a)/n

summa=0
do i=1, n, 1
        x=a+(i-0.5)*h
        summa=summa+f(x)
enddo

integral=h*summa

end function rectan_real_4

function rectan_real_8(f,a,b,n) result(integral)
implicit none
real(8) f, a, b, h, x, summa, integral
integer n, i       

h=(b-a)/n

summa=0
do i=1, n, 1
        x=a+(i-0.5)*h
        summa=summa+f(x)
enddo

integral=h*summa

end function rectan_real_8

function rectan_real_10(f,a,b,n) result(integral)
implicit none
real(10) f, a, b, h, x, summa, integral
integer n, i       

h=(b-a)/n

summa=0
do i=1, n, 1
        x=a+(i-0.5)*h
        summa=summa+f(x)
enddo

integral=h*summa

end function rectan_real_10

function rectan_real_16(f,a,b,n) result(integral)
implicit none
real(16) f, a, b, h, x, summa, integral
integer n, i       

h=(b-a)/n

summa=0
do i=1, n, 1
        x=a+(i-0.5)*h
        summa=summa+f(x)
enddo

integral=h*summa

end function rectan_real_16

!-----------

!TRAP

function trap_real_4(f,a,b,n) result(integral)
implicit none
real(4) f, a, b, h, x, summa, integral
integer n, i

h=(b-a)/n

summa=(f(a)+f(b))/2
do i=2, n, 1
        x=a+(i-1)*h
        summa=summa+f(x)
enddo

integral=h*summa

end function trap_real_4

function trap_real_8(f,a,b,n) result(integral)
implicit none
real(8) f, a, b, h, x, summa, integral
integer n, i

h=(b-a)/n

summa=(f(a)+f(b))/2
do i=2, n, 1
        x=a+(i-1)*h
        summa=summa+f(x)
enddo

integral=h*summa

end function trap_real_8

function trap_real_10(f,a,b,n) result(integral)
implicit none
real(10) f, a, b, h, x, summa, integral
integer n, i

h=(b-a)/n

summa=(f(a)+f(b))/2
do i=2, n, 1
        x=a+(i-1)*h
        summa=summa+f(x)
enddo

integral=h*summa

end function trap_real_10

function trap_real_16(f,a,b,n) result(integral)
implicit none
real(16) f, a, b, h, x, summa, integral
integer n, i

h=(b-a)/n

summa=(f(a)+f(b))/2
do i=2, n, 1
        x=a+(i-1)*h
        summa=summa+f(x)
enddo

integral=h*summa

end function trap_real_16

!-----------

!SIM

function sim_real_4(f,a,b,n) result(integral)
implicit none
real(4) f, a, b, h, x, summa_ch, summa_nech, summa, integral
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

end function sim_real_4

function sim_real_8(f,a,b,n) result(integral)
implicit none
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

end function sim_real_8

function sim_real_10(f,a,b,n) result(integral)
implicit none
real(10) f, a, b, h, x, summa_ch, summa_nech, summa, integral
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

end function sim_real_10

function sim_real_16(f,a,b,n) result(integral)
implicit none
real(16) f, a, b, h, x, summa_ch, summa_nech, summa, integral
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

end function sim_real_16

end module quadr
