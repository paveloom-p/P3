prog: main.o guarda.o
	gfortran $^ -o $@
main.o: main.f95 guarda.mod
	gfortran -c $<
guarda.o guarda.mod:guarda.f95
	gfortran -c $<
clean:
	rm -f *.mod *.o prog
result: prog
	./prog < input > result