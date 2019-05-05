program tsfs2p18
integer, parameter :: ninp=5, nres=6
integer n, i
real fun0, t0, tn, ht, t, x
open(unit=ninp, file='input')
open(unit=nres, file='result', status='replace')
read(ninp,'(e15.7)') t0, tn
read(ninp,'(i15)') n; close(ninp)
write(nres, *) ' #  t0=',t0,'  tn=',tn,'   n=',n
ht=(tn-t0)/n
write(nres,1100)
do i=0, n
t=t0+i*ht; x=exp(-t*t); write(nres,1001) i, t, fun0(x)
enddo
close(nres)
1100 format(1x,' #',2x,'i',12x,'t',11x,'fun0')
1001 format(1x,i5,2x,2x,e15.7,e15.7)
end
