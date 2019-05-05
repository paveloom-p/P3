set terminal wxt dashed size 700,700 enhanced font 'Verdana,10' persist
# set terminal eps dashed size 25,25 enhanced font 'Helvetica,14'
# set output "output.eps"

set style line 1 linecolor rgb '#0060ad' linetype 1 linewidth 2 pointtype 7 pointsize 1.5
set style line 2 linecolor rgb '#dd181f' linetype 5 linewidth 2 pointtype 5 pointsize 1.5
plot "plotting_data" index 0 with linespoints linestyle 1, "" index 1 with linespoints linestyle 2
