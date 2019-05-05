program main
implicit none

!-----------QUADRA---------------

interface
        function trap(f,a,b,n) result(integral)
                real(8) f, a, b, h, summa, integral
                integer n, i, ier
        end function trap
        function sim(f,a,b,n) result(integral)
                real(8) f, a, b, h, summa_ch, summa_nech, summa, integral
                integer n, i, ier
        end function sim
        function rectan(f,a,b,n) result(integral)
                real(8) f, a, b, h, summa, integral
                integer n, i, ier       
        end function rectan
        function f(x) result(y)     
                real(8) x, y
        end function f
        function quadra(fqua,f,a,b,n) result(i)
                real(8), external :: fqua, f
                integer n                                  
                real(8) a, b, i
        end function quadra
end interface

integer n                                  
real(8) a, b                                      

!-----------QUADRA---------------

write(*,'(/,a,/)') ' QUADRA '

!-----------------------------------------------------------------------------

read(*,*) a, b, n

write(*,*) 'Был задан промежуток: от ', a, 'до ', b, '.'
write(*,'(a,i3,2x,a,i3,2x,a)') ' Будут рассмотрены первые ', n, &
'значений(-я) (', n+1, ', если trap или sim).'
write(*,*) 'Integral=~', quadra(sim, f, a, b, n)

write(*,*) ' '

end
