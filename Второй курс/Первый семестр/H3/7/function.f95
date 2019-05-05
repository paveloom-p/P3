module function
implicit none

integer, private :: i
character(1) glas(168:344)

contains

function coder(a) result(b)
character(*) a
character(len(a)) b

b=a

!Проверить корректность замены:
!do i=168, 255, 1
!write(*,*) i, glas(i), glas(i+89), i+89
!enddo

do i=1, len_trim(a), 1
        if (ichar(a(i:i)) .ge. 168 .and. ichar(a(i:i)) .le. 255) then
        b(i:i)=glas((ichar(a(i:i))+89))
        endif
enddo

end function coder

end module
