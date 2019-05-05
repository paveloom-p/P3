module functions
implicit none

!ПАРАМЕТРЫ

!Пусть на подинтегральную функцию действуют смещения p, q и r;

real(4) p_r4, q_r4, r_r4
real(8) p_r8, q_r8, r_r8
real(10) p_r10, q_r10, r_r10
real(16) p_r16, q_r16, r_r16

contains


function f_real_4(x) result(y)
implicit none   
real(4) x, y               
y=x+p_r4+q_r4+r_r4  
end function f_real_4

function f_real_8(x) result(y)
implicit none
real(8) x, y               
y=x+p_r8+q_r8+r_r8  
end function f_real_8

function f_real_10(x) result(y)
implicit none
real(10) x, y                
y=x+p_r10+q_r10+r_r10  
end function f_real_10

function f_real_16(x) result(y)
implicit none
real(16) x, y                
y=x+p_r16+q_r16+r_r16
end function f_real_16

end module functions
