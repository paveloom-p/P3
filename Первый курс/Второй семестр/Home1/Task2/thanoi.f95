module thanoi
implicit none
contains

subroutine ring(from,where)
implicit none
integer from, where

write(*,*) 'С ',from,' на ',where

end subroutine ring



recursive subroutine hanoi(n,p,q)
implicit none
integer n, p, q

if (n.eq.1) then
call ring(p,q)
else
call hanoi(n-1,p,6-p-q)
call ring(p,q);
call hanoi(n-1,6-p-q,q)
endif
end subroutine hanoi

end module thanoi
