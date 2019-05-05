module quadr
implicit none

real(4), private :: h_r4, x_r4, summa_r4
real(8), private :: h_r8, x_r8, summa_r8
real(10), private :: h_r10, x_r10, summa_r10
real(16), private :: h_r16, x_r16, summa_r16
integer, private :: i

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
real(4) f, a, b, integral
integer n    

h_r4=(b-a)/n

summa_r4=0
do i=1, n, 1
        x_r4=a+(i-0.5)*h_r4
        summa_r4=summa_r4+f(x_r4)
enddo

integral=h_r4*summa_r4

end function rectan_real_4

function rectan_real_8(f,a,b,n) result(integral)
implicit none
real(8) f, a, b, integral
integer n    

h_r8=(b-a)/n

summa_r8=0
do i=1, n, 1
        x_r8=a+(i-0.5)*h_r8
        summa_r8=summa_r8+f(x_r8)
enddo

integral=h_r8*summa_r8

end function rectan_real_8

function rectan_real_10(f,a,b,n) result(integral)
implicit none
real(10) f, a, b, integral
integer n     

h_r10=(b-a)/n

summa_r10=0
do i=1, n, 1
        x_r10=a+(i-0.5)*h_r10
        summa_r10=summa_r10+f(x_r10)
enddo

integral=h_r10*summa_r10

end function rectan_real_10

function rectan_real_16(f,a,b,n) result(integral)
implicit none
real(16) f, a, b, integral
integer n  

h_r16=(b-a)/n

summa_r16=0
do i=1, n, 1
        x_r16=a+(i-0.5)*h_r16
        summa_r16=summa_r16+f(x_r16)
enddo

integral=h_r16*summa_r16

end function rectan_real_16

!-----------

!TRAP

function trap_real_4(f,a,b,n) result(integral)
implicit none
real(4) f, a, b, integral
integer n

h_r4=(b-a)/n

summa_r4=(f(a)+f(b))/2
do i=2, n, 1
        x_r4=a+(i-1)*h_r4
        summa_r4=summa_r4+f(x_r4)
enddo

integral=h_r4*summa_r4

end function trap_real_4

function trap_real_8(f,a,b,n) result(integral)
implicit none
real(8) f, a, b, integral
integer n

h_r8=(b-a)/n

summa_r8=(f(a)+f(b))/2
do i=2, n, 1
        x_r8=a+(i-1)*h_r8
        summa_r8=summa_r8+f(x_r8)
enddo

integral=h_r8*summa_r8

end function trap_real_8

function trap_real_10(f,a,b,n) result(integral)
implicit none
real(10) f, a, b, integral
integer n

h_r10=(b-a)/n

summa_r10=(f(a)+f(b))/2
do i=2, n, 1
        x_r10=a+(i-1)*h_r10
        summa_r10=summa_r10+f(x_r10)
enddo

integral=h_r10*summa_r10

end function trap_real_10

function trap_real_16(f,a,b,n) result(integral)
implicit none
real(16) f, a, b, integral
integer n

h_r16=(b-a)/n

summa_r16=(f(a)+f(b))/2
do i=2, n, 1
        x_r16=a+(i-1)*h_r16
        summa_r16=summa_r16+f(x_r16)
enddo

integral=h_r16*summa_r16

end function trap_real_16

!-----------

!SIM

function sim_real_4(f,a,b,n) result(integral)
implicit none
real(4) f, a, b, summa_ch, summa_nech, integral
integer n

h_r4=(b-a)/n

summa_ch=0
do i=2, n, 2
        x_r4=a+(i-1)*h_r4
        summa_ch=summa_ch+f(x_r4)
enddo

summa_nech=0
do i=3, n-1, 2
        x_r4=a+(i-1)*h_r4
        summa_nech=summa_nech+f(x_r4)
enddo

summa_r4=f(a)+f(b)+4*summa_ch+2*summa_nech

integral=(h_r4/3)*summa_r4

end function sim_real_4

function sim_real_8(f,a,b,n) result(integral)
implicit none
real(8) f, a, b, summa_ch, summa_nech, integral
integer n

h_r8=(b-a)/n

summa_ch=0
do i=2, n, 2
        x_r8=a+(i-1)*h_r8
        summa_ch=summa_ch+f(x_r8)
enddo

summa_nech=0
do i=3, n-1, 2
        x_r8=a+(i-1)*h_r8
        summa_nech=summa_nech+f(x_r8)
enddo

summa_r8=f(a)+f(b)+4*summa_ch+2*summa_nech

integral=(h_r8/3)*summa_r8

end function sim_real_8

function sim_real_10(f,a,b,n) result(integral)
implicit none
real(10) f, a, b, summa_ch, summa_nech, integral
integer n

h_r10=(b-a)/n

summa_ch=0
do i=2, n, 2
        x_r10=a+(i-1)*h_r10
        summa_ch=summa_ch+f(x_r10)
enddo

summa_nech=0
do i=3, n-1, 2
        x_r10=a+(i-1)*h_r10
        summa_nech=summa_nech+f(x_r10)
enddo

summa_r10=f(a)+f(b)+4*summa_ch+2*summa_nech

integral=(h_r10/3)*summa_r10

end function sim_real_10

function sim_real_16(f,a,b,n) result(integral)
implicit none
real(16) f, a, b, summa_ch, summa_nech, integral
integer n

h_r16=(b-a)/n

summa_ch=0
do i=2, n, 2
        x_r16=a+(i-1)*h_r16
        summa_ch=summa_ch+f(x_r16)
enddo

summa_nech=0
do i=3, n-1, 2
        x_r16=a+(i-1)*h_r16
        summa_nech=summa_nech+f(x_r16)
enddo

summa_r16=f(a)+f(b)+4*summa_ch+2*summa_nech

integral=(h_r16/3)*summa_r16

end function sim_real_16

end module quadr
