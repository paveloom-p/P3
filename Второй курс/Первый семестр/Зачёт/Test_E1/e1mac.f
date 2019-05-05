      function e1mac(x)           ! Функция e1mac(x) вычисляет для заданного
      implicit none               ! аргумента ( x ) значение выражения
      real*8 x,e1mac,a,s,s0       !          
      integer k                   !         E1(x)+gamma+ln(x)
      a=1                         !         -----------------     ,  где
      s=0d0                       !                 x
      k=1                         !                
   10 continue                    !  gamma - постоянная Эйлера, по     
        s0=s                      !  разложению в ряд Маклорена:
        s=s+a                     ! 
                                  !                    
        a=-a*x/(k+1)/(k+1)*k      !               (-1)^(k+1) x^(k-1)
        k=k+1                     !        [k=1,) ------------------     
      if (s0.ne.s) goto 10        !                  k * k! 
      e1mac=s
      end
