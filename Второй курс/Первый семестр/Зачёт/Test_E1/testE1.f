      program teste1               ! Программа тестирует алгоритмы расчета
      implicit none                ! первой интегрально-показательной функции
      integer  nres                ! E1(x). Первый -- функция e1mac(x) -- ведет
      real*8 rel                   ! расчет через разложение в ряд Маклорена.
      real*8 x,r1,r2,a,b,h,gamma   ! выражения
      real*8 e1mac, e141a, e141c   !                E1(x)+gamma+ln(x) 
      integer n, i                 !
      real*8 e141
                                   ! Второй -- функция e141a(x) -- через   
      gamma=.57721566490153286061d0!
      a=0.1d0; b=1.1d0; n=10       ! разложение по смещённым полиномам Чебышева
      h=(b-a)/n                    ! при рабочем диапазоне аргумента:
      write(*,1200)                ! 
      do i=0,n                     !            0 <= x <= 8.
        x=a+i*h; r1=e1mac(x)*x     !            
                 r2=e141a(x)       ! Третья -- функция e141c(x) -- через 
        rel=dabs(r1-r2)/r2         ! разложение по смещённым полиномам
        write(*,1201) x,r1,r2, rel ! Чебышева функции E1(x)*x*exp(x) при
      enddo                        ! рабочем диапазоне аргумента: x >= 5.
      a=1; b=9; n=8                 
      h=(b-a)/n; write(*,'(/)')                         
      write(*,1200)               
      do i=0,n                      
        x=a+i*h; r1=e1mac(x)*x      
                 r2=e141a(x) 
        rel=dabs(r1-r2)/r2		       
        write(*,1201) x,r1, r2, rel  
      enddo                        
      a=5d0; b=8d0; n=3            ! Сравнение работы e141a и e141с в
      h=(b-a)/n; write(*,1300)     ! области перекрытия рабочих диапазонов:
      do i=0,n                     !            5 <=  x <= 8.
        x=a+i*h; r1=e141a(x)                            
        r2=e141c(x)*dexp(-x)/x + gamma+dlog(x)
        rel=dabs(r1-r2)/r2
        write(*,1201) x, r1, r2, rel
      enddo

 1200 format(1x,' #  x',12x,'e1mac(x)*x     = ',4x,
     >                          'e141a=E1+g+ln(x)',6x,'relerr')
 1201 format(1x,2x,f5.2, d25.17, d25.17, d10.2) 
 1300 format(//1x,' #  x',7x,'e141a=E1(x)+g+ln(x) = ',
     >                     'e141c(x)exp(-x)/x+g+ln(x)',2x,'reller')
      end
