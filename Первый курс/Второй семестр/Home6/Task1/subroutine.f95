module subroutine_
implicit none
contains

!----------------------------------
subroutine lbound_(a,l)
integer a(l:), b(1)
integer l

b=lbound(a)
write(*,*) 'Нижняя граница измерения: ', b

end subroutine lbound_
!----------------------------------
subroutine ubound_(a,l)
integer a(l:), b(1)
integer l

b=ubound(a)
write(*,*) 'Верхняя граница измерения:', b

end subroutine ubound_
!----------------------------------
subroutine array(a,l)
integer l
integer a(l:)

read(*,*) a
write(*,*) 'Был введен вектор:', a

call lbound_(a,l)
call ubound_(a,l)

end subroutine array

end module subroutine_
