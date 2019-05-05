module function
implicit none
contains

function summ(a,b) result(c)
character(*) a, b
character(max(len(a),len(b))+1) c

integer(1) d(0:max(len(a),len(b)))

integer(1) ib
integer i

d=0
c=' '

do i=len(a), 1, -1
read(a(i:i),*) d(i+dim(len(b),len(a)))
enddo

do i=len(b), 1, -1
read(b(i:i),*) ib
if (d(i+dim(len(a),len(b)))+ib .gt. 9) d(i+dim(len(a),len(b))-1)=d(i+dim(len(a),len(b))-1)+1_1
write(c(i+1+dim(len(a),len(b)):i+1+dim(len(a),len(b))),'(i1)') mod(d(i+dim(len(a),len(b)))+ib,10)
enddo

do i=dim(len(a),len(b)), 1, -1
if (d(i) .eq. 10) d(i-1)=d(i-1)+1_1; d(i)=0
write(c(i+1:i+1),'(i1)') d(i)
enddo

if (d(0) .ne. 0) write(c(1:1),'(i1)') d(0)

end function summ

end module
