function fun0(x)
use my_prec
implicit none
real(mp) fun0, x
fun0=(1.7_mp-(1.7_mp**3-x)**(1.0_mp/3))/x
end
