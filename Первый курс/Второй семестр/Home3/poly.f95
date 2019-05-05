module poly
implicit none
contains

!-----------

function poly1(n ,x ,m) result(y)
integer n, i, x
real y
real,dimension(100)::m

y=m(1)
do i=1, n, 1
y=y*x+m(i+1)
enddo

end function poly1

!-----------

function poly2(x, m, n)

integer i, n, x
real poly2(2), m(0:n), s1, s2, ax, x2
ax=abs(x)
select case (n)
	case(0); poly2(1)=m(0); poly2(2)=m(0); return
	case(1); poly2(1)=m(1)+ax*m(0); poly2(2)=m(1)-ax*m(0); return
	case default; x2=x*x; s1=m(0); do i=2,n,2; s1=s1*x2+m(i); enddo
			      s2=m(1); do i=3,n,2; s2=s2*x2+m(i); enddo
	if (mod(n,2)==0) then; poly2(1)=s1+ax*s2; poly2(2)=s1-ax*s2
			 else; poly2(1)=s2+ax*s1; poly2(2)=s2-ax*s1
	endif
endselect

end function poly2

!-----------

function poly3(n, t)
integer i, n, t(:)
integer poly3
integer cm16, cm15, cm14, cm13, cm12, cm11, cm10, cm9, cm8, cm7, cm6, cm5, cm4, cm3, cm2, cm1, c0, c1, c2, c3, c4, c5, c6, c7, &
&c8, c9, c10, c11, c12, c13, c14, c15, c16
integer, dimension(33)::c

do i=1, 33, 1
	c(i)=0
enddo

do i=1, n, 1
	select case (t(i))
		case(-16); c(1)=c(1)+1
		case(-15); c(2)=c(2)+1
		case(-14); c(3)=c(3)+1
		case(-13); c(4)=c(4)+1
		case(-12); c(5)=c(5)+1
		case(-11); c(6)=c(6)+1
		case(-10); c(7)=c(7)+1
		case(-9); c(8)=c(8)+1
		case(-8); c(9)=c(9)+1
		case(-7); c(10)=c(10)+1
		case(-6); c(11)=c(11)+1
		case(-5); c(12)=c(12)+1
		case(-4); c(13)=c(13)+1
		case(-3); c(14)=c(14)+1
		case(-2); c(15)=c(15)+1
		case(-1); c(16)=c(16)+1
		case(0); c(17)=c(17)+1
		case(1); c(18)=c(18)+1
		case(2); c(19)=c(19)+1
		case(3); c(20)=c(20)+1
		case(4); c(21)=c(21)+1
		case(5); c(22)=c(22)+1
		case(6); c(23)=c(23)+1
		case(7); c(24)=c(24)+1
		case(8); c(25)=c(25)+1
		case(9); c(26)=c(26)+1
		case(10); c(27)=c(27)+1
		case(11); c(28)=c(28)+1
		case(12); c(29)=c(29)+1
		case(13); c(30)=c(30)+1
		case(14); c(31)=c(31)+1
		case(15); c(32)=c(32)+1
		case(16); c(33)=c(33)+1
	endselect
enddo

write(*,'(i3,i15)') (i-17,c(i),i=1,33)

end function poly3

!-----------

function poly4(v)
integer poly4
integer v, k, c, i
integer a(32)

k=abs(v)

do i=1, 32, 1
	c=mod(k,2)
	k=(k-c)/2
	a(i)=c
enddo

if(v .lt. 0) then
	a(32)=1
endif

write(*,*) 'Слева - номер разряда, справа - его содержимое.'
write(*,'(i3,i15)') (i,a(i),i=1,32)

end function poly4

!-----------

end module poly
