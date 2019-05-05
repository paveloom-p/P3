set terminal wxt size 1024,700 enhanced font 'Verdana,10' persist
set output "output.eps"

splot "data_decart_out" with lines
