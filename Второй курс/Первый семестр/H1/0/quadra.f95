module quadra
implicit none
contains

!-----------

function f(x) result(y)     
real(8) x, y                              
y=x   !**3                                  
end function f      
                                        
!-----------                            

function rectan(y,a,b,n) result(integral)
real(8) y(:), a, b, h, summa, integral
integer n, i

h=(b-a)/n

summa=0
do i=1, n,1
        summa=summa+y(i)
enddo

integral=h*summa

!write(*,*) 'Summa=',summa
!write(*,*) 'Integral=~', integral

end function rectan

!-----------

!-----------

function trap(y,a,b,n) result(integral)
real(8) y(:), a, b, h, summa, integral
integer n, i

h=(b-a)/n

summa=(y(1)+y(n+1))/2
do i=2, n, 1
        summa=summa+y(i)
enddo

integral=h*summa

!write(*,*) 'Summa=',summa
!write(*,*) 'Integral=~', integral

end function trap

!-----------

function sim(y,a,b,n) result(integral)
real(8) y(:), a, b, h, summa_ch, summa_nech, summa, integral
integer n, i

h=(b-a)/n

summa_ch=0
do i=2, n, 2
        summa_ch=summa_ch+y(i)
enddo

summa_nech=0
do i=3, n-1, 2
        summa_nech=summa_nech+y(i)
enddo

summa=y(1)+y(n+1)+4*summa_ch+2*summa_nech

integral=(h/3)*summa

!write(*,*) 'Summa=',summa
!write(*,*) 'Integral=~', integral

end function sim

!-----------

end module quadra
