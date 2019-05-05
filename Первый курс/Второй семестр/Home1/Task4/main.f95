program tstbis
implicit none
real(8) f, aa, bb, epsx, epsy, resn, resr
real t1, t2, t3, t4
character(8) mess(2) /'bisecn0','bisecr0'/
integer itern, iterr, iern, ierr, i
integer ninp / 5 /, nres / 6 /, krep / 1000000 /

write(*,*) 'Результат будет выведен в файл result.'

open (unit=ninp,file='input')
open (unit=nres,file='result',status='replace')

read (ninp,*) aa, bb, epsx, epsy
write(nres,1100) aa, bb, epsx, epsy, krep

call second(t1)
do i=1,krep
call bisecn0(aa,bb,epsx,epsy,resn,itern,iern)
enddo

call second(t2)
write(nres,1101) resn, f(resn), itern, iern, t2-t1, mess(1)

call second(t3)
do i=1,krep
iterr=0
call bisecr0(aa,bb,epsx,epsy,resr,iterr,ierr)
enddo

call second(t4)
write(nres,1101) resr, f(resr), iterr, ierr, t4-t3, mess(2)
close(nres)

1100 format(1x,' a=',d15.7,' b=',d15.7/ 1x,'epsx=',d15.7,&
& ' epsy=',d15.7,' krep=',i10/13x,'res',14x,'f(res)',3x,'iter',3x,'ier',&
& 3x,'Время',3x,'Подпрограмма')
1101 format(1x,d25.17,2x,d9.2,2x,i5,2x,i2,2x,f7.3,6x,a)

end
