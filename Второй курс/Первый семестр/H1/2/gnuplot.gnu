set terminal eps
set output 'graph.eps'
plot 'results.dat' with dots

set output 'graph_t.eps'
plot 'results_t.dat' with dots

set output 'graph_s.eps'
plot 'results_s.dat' with dots 
