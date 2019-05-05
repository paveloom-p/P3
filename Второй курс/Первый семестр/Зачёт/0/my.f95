module my
use real_type
implicit none

contains

subroutine my_function(x,eps,res,k,ier)
use real_type
implicit none

real(mp), intent(in) :: x, eps
real(mp), intent(out) :: res
integer, intent(out) :: k        ! Число слагаемых

! Вспомогательные переменные
real(mp) b, koef_1, koef_2, koef_3, koef_4, epsi
real(mp) kk, mm
integer m, ier

epsi=epsilon(1.0_mp) ! Погрешность ошибки округления

koef_1=144.0_mp/7.0_mp ! Коэффициент домножения ряда
koef_2=34.0_mp/7.0_mp-27.0_mp/7.0_mp/x+(log(2.0_mp)-log(x)+log(x/2.0_mp))*koef_1 ! Сумма слагаемых минус первого и нулевого порядка,
                                                                                 ! домноженная на koef_1
res=-15.0_mp*x/1536.0_mp      ! Слагаемое первого порядка
koef_3=31.0_mp/32.0_mp        ! Начальные значения коэффициентов перед
koef_4=1.0_mp/32.0_mp         ! применением реккурентной формулы для коэффициентов
b=-res*x*64.0_mp/375.0_mp     ! Слагаемое второго порядка (рекуррентная формула)

k=2 ! На данном этапе в сумме правильной части ряда два слагаемых

do while (abs(koef_3*b) .gt. abs(res)*eps) ! Главный цикл

! Проверка на достижение суммой погрешности ошибки округления
if (abs(res)*eps .lt. abs(res)*epsi) then; ier=1; exit; endif

res=res+koef_3*b ! Прибавление вычисленного члена ряда

k=k+1; m=k+3; kk=k; mm=m*m; ! Вычисление следующего члена ряда
koef_4=-koef_4/2.0_mp ! Рекуррентные формулы для коэффициентов
koef_3=-koef_3+koef_4
b=b*x*(kk+2.0_mp)/mm  ! Рекуррентная формула для модульной части ряда

enddo

res=koef_1*res+koef_2 ! Домножение на общий коэффициент и финальное суммирование

end subroutine my_function

end module my
