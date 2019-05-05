module quadr
implicit none
contains

function quadrauto(fqua,f,a,b,n,eps,k,ier) result(i)
real(8) fqua, i
real(8), external :: f
integer n, ier, k                          
real(8) a, b, eps

i=fqua(f,a,b,n,eps,k,ier)

end function quadrauto

function rectan(f,a,b,n,eps,k,ier) result(integral)
real(8) f, a, b, h, x, summa, integral, check, eps
integer n, ier, k
integer i, m  

check=0

do m=1, k, 1

n=2**m
!write(*,*) n

h=(b-a)/n

summa=0
do i=1, n, 1
        x=a+(i-0.5)*h
        summa=summa+f(x)
enddo

integral=h*summa
!write(*,*) integral
!write(*,*) abs(integral-check)
!write(*,*) m, k

if (abs(integral-check) .gt. eps) then
        if (m .eq. k) then
        ier=1
        endif
        continue
else
        ier=0
        k=m
        exit
endif

check=integral

enddo

end function rectan

!-----------

function trap(f,a,b,n,eps,k,ier) result(integral)
real(8) f, a, b, h, x, summa, integral, check, eps
integer n, ier, k
integer i, m

check=0

do m=1, k, 1

n=2**m

h=(b-a)/n

summa=(f(a)+f(b))/2
do i=2, n, 1
        x=a+(i-1)*h
        summa=summa+f(x)
enddo

integral=h*summa

if (abs(integral-check) .gt. eps) then
        if (m .eq. k) then
        ier=1
        endif
        continue
else
        ier=0
        k=m
        exit
endif

check=integral

enddo

end function trap

!-----------

function sim(f,a,b,n,eps,k,ier) result(integral)
real(8) f, a, b, h, x, summa_ch, summa_nech, summa, integral, check, eps
integer n, ier, k
integer i, m

check=0

do m=1, k, 1

n=2**m

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

if (abs(integral-check) .gt. eps) then
        if (m .eq. k) then
        ier=1
        endif
        continue
else
        ier=0
        k=m
        exit
endif

check=integral

enddo

end function sim

end module quadr
