module RungeKutta
use mpmodule
implicit none

contains

function F(Q,L,m,G,R,k)
integer :: Q,L,k,i,s,j
type(mp_real) :: m(Q)
type(mp_real) :: R(Q,2*L),Ro,F(2*L),G

 do i=1,2*L
    F(i)='0'
 enddo

 F(1:L)=R(k,L+1:2*L)

 do i=L+1,L*2
    do s=1,Q
	if (s/=k) then
	    Ro=0
	    do j=1,L
		Ro=Ro+(R(k,j)-R(s,j))**2
	    enddo
	    F(i)=F(i)+G*m(s)*(R(s,i-L)-R(k,i-L))/Ro**(3.d0/2.d0)
	endif
    enddo
 enddo

end function F


subroutine RK(Q,L,T,eps,h,m,R0,En,Hn,Gn)
integer :: Q,L
integer :: i,s,k
type(mp_real) :: m(Q),T,eps
type(mp_real) :: R(Q,2*L),R0(Q,2*L),error(Q,L*2)
type(mp_real) :: E1,Elast,error_s
type(mp_real) :: k1(Q,L*2),k2(Q,L*2),k3(Q,L*2),k4(Q,L*2),k5(Q,L*2),k6(Q,L*2),k_add(Q,L*2)
type(mp_real) :: h,sum_h,G
character*3 body
character*1 En,Hn,Gn

 G=4*(4*atan(mpreal(1.d0)))**2

 if (En=='G' .or. Hn=='G' .or. Gn=='G') then
 G='1'
 endif

 sum_h='0'

 if (En=='E' .or. Hn=='E' .or. Gn=='E') then
    E1=energy(Q,L,m,G,R0)
 endif

 ! do k=1,Q
!    write(body,'(I3.3)') k
!    open(k+10,file='./results/X_of_body_'//trim(body)//'.dat')
!    open(Q+k+10,file='./results/Y_of_body_'//trim(body)//'.dat')
!    open(2*Q+k+10,file='./results/Z_of_body_'//trim(body)//'.dat')
!    open(3*Q+k+10,file='./results/V_x_of_body_'//trim(body)//'.dat')
!    open(4*Q+k+10,file='./results/V_y_of_body_'//trim(body)//'.dat')
!    open(5*Q+k+10,file='./results/V_z_of_body_'//trim(body)//'.dat')
 ! enddo

 do k=1,Q
    write(body,'(I3.3)') k
    open(k+10,file='./results/X_of_body_'//trim(body)//'.dat')
    open(Q+k+10,file='./results/V_x_of_body_'//trim(body)//'.dat')
    open(2*Q+k+10,file='./results/Y_of_body_'//trim(body)//'.dat')
    open(3*Q+k+10,file='./results/V_y_of_body_'//trim(body)//'.dat')
    if (L==3) then
     open(4*Q+k+10,file='./results/Z_of_body_'//trim(body)//'.dat')
     open(5*Q+k+10,file='./results/V_z_of_body_'//trim(body)//'.dat')
    endif
 enddo


10 do while (sum_h <= T)

					

    do s=1,Q
	k1(s,:)=F(Q,L,m,G,R0(:,:),s)
    enddo
    k1=multiply_mp(k1,h)

					

    do s=1,Q
	k2(s,:) = F(Q,L,m,G,add_array_mp(R0,divide_int(k1,4)),s)
    enddo
    k2=multiply_mp(k2,h)

					

    do s=1,Q
	k3(s,:) = F(Q,L,m,G,add_array_mp(add_array_mp(R0,divide_int(multiply_int(k1,3),32)),divide_int(multiply_int(k2,9),32)),s)
    enddo
    k3=multiply_mp(k3,h)

					

    k_add=divide_int(multiply_int(k1,1932),2197)
    k_add=deduct_array_mp(k_add,divide_int(multiply_int(k2,7200),2197))
    k_add=add_array_mp(k_add,divide_int(multiply_int(k2,7296),2197))
    do s=1,Q
	k4(s,:) = F(Q,L,m,G,add_array_mp(R0,k_add),s)
    enddo
    k4=multiply_mp(k4,h)

					

    k_add=divide_int(multiply_int(k1,439),216)
    k_add=deduct_array_mp(k_add,multiply_int(k2,8))
    k_add=add_array_mp(k_add,divide_int(multiply_int(k3,3680),513))
    k_add=deduct_array_mp(k_add,divide_int(multiply_int(k4,845),4104))
    do s=1,Q
	k5(s,:) = F(Q,L,m,G,add_array_mp(R0,k_add),s)
    enddo
    k5=multiply_mp(k5,h)

					

    k_add=multiply_int(k2,2)
    k_add=deduct_array_mp(k_add,divide_int(multiply_int(k1,8),27))
    k_add=deduct_array_mp(k_add,divide_int(multiply_int(k3,3544),2565))
    k_add=add_array_mp(k_add,divide_int(multiply_int(k4,1859),4104))
    k_add=deduct_array_mp(k_add,divide_int(multiply_int(k5,11),40))
    do s=1,Q
	k6(s,:) = F(Q,L,m,G,add_array_mp(R0,k_add),s)
    enddo
    k6=multiply_mp(k6,h)

					

    if (En=='H' .or. Hn=='H' .or. Gn=='H') then

    error=divide_int(k1,360)
    error=deduct_array_mp(error,divide_int(multiply_int(k3,128),4275))
    error=deduct_array_mp(error,divide_int(multiply_int(k4,2197),75240))
    error=add_array_mp(error,divide_int(k5,50))
    error=add_array_mp(error,divide_int(multiply_int(k6,2),55))

					

    do s=1,Q
	error_s=0
	do i=1,2*L
	    error_s=error_s+error(s,i)**2
	enddo
	if (sqrt(error_s)>eps) then
	    h=h/2
	    goto 10
	endif
    enddo

    endif

    R=add_array_mp(R0,divide_int(multiply_int(k1,16),135))
    R=add_array_mp(R,divide_int(multiply_int(k3,6656),12825))
    R=add_array_mp(R,divide_int(multiply_int(k4,28561),56430))
    R=deduct_array_mp(R,divide_int(multiply_int(k5,9),50))
    R=add_array_mp(R,divide_int(multiply_int(k6,2),55))

    R0=R

    do k=1,Q
	do i=1,L
	    call mpwrite(k+10,R(k,i))
	enddo
    enddo

    do k=1,Q
	do i=1,2*L
	    call mpwrite(Q*(i-1)+k+10,R(k,i))
	enddo
    enddo

    sum_h = sum_h + h
    enddo

 do k=1,6*Q+10
    close(k+10)
 enddo

 if(En=='E' .or. Hn=='E' .or. Gn=='E') then
    Elast=energy(Q,L,m,G,R)

    write(*,*) 'Полная энергия в начальный момент:'
    call mpwrite(6,E1)
    write(*,*) 'Полная энергия в конечный момент:'
    call mpwrite(6,Elast)
    write(*,*) 'Относительная ошибка:'
    call mpwrite(6,abs(Elast/E1-1))
 endif

end subroutine RK

function add_array_mp(a,b)
integer :: i,k
type(mp_real), dimension(:,:), intent(in) :: a,b
type(mp_real) :: add_array_mp(ubound(a,1),ubound(a,2))

 do i=lbound(a,1),ubound(a,1)
    do k=lbound(a,2),ubound(a,2)
	add_array_mp(i,k)=a(i,k)+b(i,k)
    enddo
 enddo

end function add_array_mp

function deduct_array_mp(a,b)
integer :: i,k
type(mp_real), dimension(:,:), intent(in) :: a,b
type(mp_real) :: deduct_array_mp(ubound(a,1),ubound(a,2))

 do i=lbound(a,1),ubound(a,1)
    do k=lbound(a,2),ubound(a,2)
	deduct_array_mp(i,k)=a(i,k)-b(i,k)
    enddo
 enddo

end function deduct_array_mp



function multiply_int(a,c)
integer :: i,k,c
type(mp_real),dimension(:,:) :: a
type(mp_real) :: multiply_int(ubound(a,1),ubound(a,2))

 do i=lbound(a,1),ubound(a,1)
    do k=lbound(a,2),ubound(a,2)
	multiply_int(i,k)=a(i,k)*c
    enddo
 enddo

end function multiply_int

function multiply_mp(a,c)
integer :: i,k
type(mp_real) :: c
type(mp_real),dimension(:,:) :: a
type(mp_real) :: multiply_mp(ubound(a,1),ubound(a,2))

 do i=lbound(a,1),ubound(a,1)
    do k=lbound(a,2),ubound(a,2)
	multiply_mp(i,k)=a(i,k)*c
    enddo
 enddo

end function multiply_mp


function divide_int(a,c)
integer :: i,k,c
type(mp_real),dimension(:,:) :: a
type(mp_real) :: divide_int(ubound(a,1),ubound(a,2))

 do i=lbound(a,1),ubound(a,1)
    do k=lbound(a,2),ubound(a,2)
	divide_int(i,k)=a(i,k)/mpreal(c)
    enddo
 enddo

end function divide_int

function energy(Q,L,m,G,R)
integer :: i,s,Q,L
type(mp_real) :: m(Q),R(Q,2*L),energy,G

 energy='0'

 do i=1,Q
    energy=energy+m(i)*(R(i,4)**2+R(i,5)**2+R(i,6)**2)/2.d0
    do s=1,Q
	if (s/=i) then
	    energy=energy-G/2*m(i)*m(s)/(sqrt((R(i,1)-R(s,1))**2 + (R(i,2)-R(s,2))**2 + (R(i,3)-R(s,3))**2))
	endif
    enddo
 enddo

end function energy


end module RungeKutta