program main
implicit none

character(50) a
read(*,*) a

write(*,'(/,a,1x,a,/)') ' Chosen letter:', a

write(*,*) 'Char:'
write(*,*) char(ichar(a(1:1))), char(ichar(a(2:2)))

write(*,*) 'Achar:'
write(*,*) achar(iachar(a(1:1))), achar(iachar(a(2:2)))

write(*,*) ' '

end
