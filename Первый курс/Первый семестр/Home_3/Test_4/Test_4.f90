program h3_4
implicit none
real x, y
integer i,j
integer, dimension (-1:1, -1:1) :: tab
tab = transpose(reshape ((/-3,-1,4,-2,0,2,-4,1,3/), shape(tab)))
write (*,*) 'Coordinates '
read (*,*) x, y
if (x/=0) then; i=x/abs(x);
else; i=0;
end if
if (y/=0) then; j=y/abs(y);
else; j=0;
end if
write (*,*) 'Key', tab(i,j);
end program
