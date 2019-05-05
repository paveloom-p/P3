program tsfs2p18
data ninp / 5 /, nres / 6 /
open(unit=ninp, file='input')
open(unit=nres, file='result', status='replace')
read(ninp,100) t0, tn
read(ninp,101) n
write(nres, *) ' #  t0=',t0,'  tn=',tn,'   n=',n
ht=(tn-t0)/n
write(nres,1100)
do i=0, n
t=t0+i*ht
x=exp(-t*t)
r0=fun0(x)
write(nres,1001) i, t, r0
enddo
close(nres)
100 format(e15.7)
101 format(i15)
1100 format(1x,' #',2x,'i',12x,'t',11x,'fun0')
1001 format(1x,i5,2x,2x,e15.7,e15.7)
end
