set terminal postscript eps enhanced 32
set output 'result.eps'
set border 15 lw 3
set key left box
set label 't'           at 6.3,0.025 center
set label 'x=exp(-t*t)' at 4 , 0.15 center
plot 'result' u 2:3 title ' fun0(x) '\
                     w lp lt 3 lw 6 pt 70 ps 3,\
     ''       u 2:4 title 'fun1(x)' w l lt -1 lw 3
reset 
