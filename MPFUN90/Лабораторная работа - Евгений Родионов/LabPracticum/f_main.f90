program f_main
use mpmodule
use RungeKutta
implicit none
integer :: i,k,mp,mpwr
integer :: Q,L
type(mp_real):: T,eps,h
type(mp_real),allocatable :: m(:)
type(mp_real),allocatable :: R0(:,:)
character*3 body
character*1 En,Hn,Gn

 call getarg(1,En)
 call getarg(2,Hn)
 call getarg(3,Gn)

 open(2000, file='./initial/mp.dat')

    read(2000,*) mp
    read(2000,*) mpwr

 close(2000)

 call mpinit (mp)                                 !Желаемая точность
 mpoud = mpwr                                     !Вывод

 open(1000,file='./initial/parameters.dat')

  read(1000,*) Q
  read(1000,*) L
  call mpread(1000,T)
  call mpread(1000,eps)
  call mpread(1000,h)

 close(1000)

 allocate(R0(Q,L*2),m(Q))

 open(1001,file='./initial/masses.dat')
 do i=1,Q
    call mpread(1001,m(i))
 enddo
 close(1001)

 do i=1,Q
    write(body,'(I3.3)') i
    open(1001+i,file='./initial/body'//trim(body)//'.dat')
    do k=1,2*L
	call mpread(1001+i,R0(i,k))
    enddo
    close(1001+i)
 enddo

 call RK(Q,L,T,eps,h,m,R0,En,Hn,Gn)

end program f_main
