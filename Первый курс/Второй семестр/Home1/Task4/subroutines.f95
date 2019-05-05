function f(x) result(w)
implicit none
real(8) w, x
w=3*x-1
end

subroutine bisecn0(aa,bb,epsx,epsy,res,k,ier)
implicit none
real(8) aa, bb, epsx, epsy, res, f, a, b, c, fa, fb, fc
integer k, ier, kmax / 1000 /

a=aa
b=bb
fa=f(a)
fb=f(b)

if ( (fa*fb).le.0d0 ) then 
k=0

do

	c=(a+b)*0.5d0
	fc=f(c)

	if ( fa*fc.gt.0d0) then
		a=c
		fa=fc
	else 
		b=c
		fb=fc
	endif

	k=k+1

	if ((k.gt.kmax).or.& 
	& (dabs(b-a).lt.epsx.and.&
	& dabs(fc ).lt.epsy)) exit

enddo

res=(a+b)*0.5d0
if (k.le.kmax) then
	ier=0 
else
	ier=1
endif

else
ier=2
endif

end subroutine bisecn0



recursive subroutine bisecr0(a,b,epsx,epsy,res,k,ier)
implicit none
real(8) aa, bb, epsx, epsy, res, f, a, b, c, fa, fb, fc
integer k, ier, kmax / 1000 /
! a=aa; b=bb;

fa=f(a)
fb=f(b)

if ( (fa*fb).le.0d0 ) then
	res=(a+b)*0.5d0
	c=res
	fc=f(c)
	
	if ( (dabs(b-a).lt.epsx) .and.&
	& (dabs( fc).lt.epsy) ) then
	ier=0 
	else

		if (k.gt.kmax) then
		ier=1
		else
		c=(a+b)*0.5d0
		fc=f(c)
		k=k+1
	
			if (f(a)*fc.gt.0) then
			call bisecr0(c,b,epsx,epsy,res,k,ier)
			else
			call bisecr0(a,c,epsx,epsy,res,k,ier)
			endif
		endif
	endif

else
ier=2
endif

end subroutine bisecr0
