program h3_6
implicit none
integer x, r, n5, n3
read(*,*) x
r=mod(x,5)
if (r==0) then 
n5=x/5; n3=0
endif

if (r==1) then 
n5=x/5-1; n3=2
endif

if (r==2) then 
n5=x/5-2; n3=4
endif

if (r==3) then 
n5=x/5; n3=1
endif

if (r==4) then 
n5=x/5-1; n3=3
endif

write(*,*) n5, n3
end
