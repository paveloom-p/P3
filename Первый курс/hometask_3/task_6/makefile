prog: binrad.o main.o
	gfortran $^ -o $@
main.o: main.f95 binrad.mod
	gfortran -c $<
binrad.o binrad.mod: binrad.f95
	gfortran -c $<
clean:
	rm -f *.mod *.o prog
result: prog
	./prog < input > result