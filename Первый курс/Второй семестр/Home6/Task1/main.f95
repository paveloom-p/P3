program main
use subroutine_
implicit none

integer a(8) / 1,2,3,4,5,6,7,8 /

integer l 

write(*,*) ' '
read(*,*) l
write(*,*) 'Был введен начальный индекс массива:', l

write(*,*) 'Текущий диапазон индексов:	     ', l, l+size(a)-1
write(*,*) ' '

call array(a,l)
end
