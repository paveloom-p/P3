ext=f95
comp=gfortran
files=$(wildcard *.f95)
obj=$(patsubst %.$(ext), %.o, $(files))
name=Run

$(name): sort.f95 quadra.f95 main.f95
	$(comp) sort.f95 quadra.f95 main.f95 -o $@
Go:$(name)
	./$(name) < input
Plot:
	gnuplot gnuplot.gnu; gv graph.eps
Plot2:
	gnuplot gnuplot.gnu; gv graph_t.eps
Plot3:
	gnuplot gnuplot.gnu; gv graph_s.eps
clean:
	rm -f $(obj) $(wildcard *.eps)
exterminate:
	rm -f $(obj) $(name) $(wildcard *.o) $(wildcard *.mod) $(wildcard *.eps)
