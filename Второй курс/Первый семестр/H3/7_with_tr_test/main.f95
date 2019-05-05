program main
implicit none

integer i, j
character(50) a(100)

a=' '
write(*,*) ' '

do i=1, 100, 1
read(*,*) a(i)
if (a(i) .eq. 'stop reading') then; exit; endif;
enddo

call system("tr абвгдежзийклмнопрстуфхцчшщъыьэюя &
&евгджёзкойлмнпурстфыхцчшщбъэьюяа < input > input_tr")

do j=1, i-1, 1
write(*,*) a(j)
enddo

end
