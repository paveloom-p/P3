prog: sort.o main.o
	gfortran $^ -o $@
main.o: main.f95 sort.mod
	gfortran -c $<
sort.o sort.mod: sort.f95
	gfortran -c $<
result: prog
	./prog
clean:
	rm -f *.mod *.o prog
