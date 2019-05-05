program h3_7
implicit none
integer n5, n3, x, s
integer, dimension (1:5) :: tab
tab = (/(0, 0), (-1, 2), (-2,4), (0, 1), (-1, 3)/)
read(*,*) x
s = tab(mod(x,5))
write(*,*) s
n5=x/5+s
write(*,*) n5, n3
end


!koeffs = [(0, 0), (-1, 2), (-2,4), (0, 1), (-1, 3)]
