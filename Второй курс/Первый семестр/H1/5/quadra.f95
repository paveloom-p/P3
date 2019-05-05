module quadr
implicit none

interface quadra
module procedure quadra_function, quadra_array
end interface

interface rectan
module procedure rectan_function, rectan_array
end interface

interface sim
module procedure sim_function, sim_array
end interface

interface trap
module procedure trap_function, trap_array
end interface

contains



function quadra_function(fqua,f,a,b,n) result(i)
real(8) fqua, i
real(8), external :: f
integer n                                  
real(8) a, b

i=fqua(f,a,b,n)

end function quadra_function

function quadra_array(fqua,arr,a,b,n) result(i)
real(8) fqua, i
real(8) arr(:)
integer n                                  
real(8) a, b

i=fqua(arr,a,b,n)

end function quadra_array

!------------------------------------------

function rectan_function(f,a,b,n) result(integral)
real(8) f, a, b, h, x, summa, integral
integer n, i       

h=(b-a)/n

summa=0
do i=1, n, 1
        x=a+(i-0.5)*h
        summa=summa+f(x)
enddo

integral=h*summa

end function rectan_function

function rectan_array(arr,a,b,n) result(integral)

real(8) arr(n), a, b, h, summa, integral
integer n, i

h=(b-a)/n

summa=0
do i=1, n,1
        summa=summa+arr(i)
enddo

integral=h*summa

end function rectan_array

!-----------

function trap_function(f,a,b,n) result(integral)
real(8) f, a, b, h, x, summa, integral
integer n, i

h=(b-a)/n

summa=(f(a)+f(b))/2
do i=2, n, 1
        x=a+(i-1)*h
        summa=summa+f(x)
enddo

integral=h*summa

end function trap_function

function trap_array(arr,a,b,n) result(integral)
real(8) arr(n), a, b, h, summa, integral
integer n, i

h=(b-a)/n

summa=(arr(1)+arr(n+1))/2
do i=2, n, 1
        summa=summa+arr(i)
enddo

integral=h*summa

end function trap_array

!-----------

function sim_function(f,a,b,n) result(integral)
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

end function sim_function

function sim_array(arr,a,b,n) result(integral)
real(8) arr(n+1), a, b, h, summa_ch, summa_nech, summa, integral
integer n, i

h=(b-a)/n

summa_ch=0
do i=2, n, 2
        summa_ch=summa_ch+arr(i)
enddo

summa_nech=0
do i=3, n-1, 2
        summa_nech=summa_nech+arr(i)
enddo

summa=arr(1)+arr(n+1)+4*summa_ch+2*summa_nech

integral=(h/3)*summa

end function sim_array

end module quadr
