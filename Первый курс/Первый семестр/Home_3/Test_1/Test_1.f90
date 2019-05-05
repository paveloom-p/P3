program h3_1
implicit none
character(5),dimension(10)::a 
data a /'zero', 'one', 'two', 'three', 'four', 'five', 'six', 'seven', 'eight', 'nine'/
integer p
p=1
p=mod(p/10)
write(*,*) a(b)
end
