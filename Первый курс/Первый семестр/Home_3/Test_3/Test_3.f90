program h3_3
implicit none
real a, pi
character what
read(*,*) what
pi=3.14;

select case (what)
	case ('C')
		read(*,*) a
		a=a/2/pi
		write(*,*) a
	case ('G')
		read(*,*) a
		write(*,*) a
	case ('M')
		read(*,*) a
		
		write(*,*) a
	case ('S')
		read(*,*) a
			do while (a>60)
			a=a-60
			enddo
		write(*,*) a
	case ('d')
		read(*,*) a
		
		write(*,*) a
	case ('h')
	
	case ('m')
	
	case ('s')
	
end select
end
