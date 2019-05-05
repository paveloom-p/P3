prog: poly.o main.o
	gfortran $^ -o $@
main.o:main.f95 poly.mod
	gfortran -c $<
poly.o poly.mod: poly.f95
	gfortran -c $<
clean:
	rm -f *.mod *.o prog
result: prog
	./prog < input > result