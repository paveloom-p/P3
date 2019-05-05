module subroutines
implicit none
contains

subroutine degree_nonrec(a,n)
integer a, n
integer r, i

if (n .ge. 1) then

r=a
do i=1, n-1, 1
	r=r*a
enddo

elseif (n .eq. 0) then
r=1
endif

write(*,*) 'Результат возведения в степень по нерекурсивному методу:         ', r

end subroutine degree_nonrec

!----------------------------------

recursive subroutine degree_rec(a,n,i,r)
integer a, n
integer r, i

if (n .ge. 1) then

if (i .lt. n-1) then
	r=r*a
	i=i+1
	call degree_rec(a,n,i,r)
endif

elseif (n .eq. 0) then
r=1
endif

end subroutine degree_rec

!----------------------------------

recursive subroutine degree_tail_rec(a,n,i,r)
integer a, n
integer r, i

if (n .ge. 1) then

if (i .lt. n-1) then
	r=r*a
	i=i+1
	call degree_tail_rec(a,n,i,r)
endif

elseif (n .eq. 0) then
r=1
endif

end subroutine degree_tail_rec


end module subroutines
