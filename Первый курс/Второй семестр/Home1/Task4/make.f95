ext=f95
comp=gfortran
files=$(wildcard *.f95)
obj=$(patsubst %.$(ext), %.o, $(files))
name=Run

$(name): subroutines.f95 main.f95
	$(comp) subroutines.f95 main.f95 -o $@
Go:$(name)
	./$(name)
clean:
	rm -f $(obj)
exterminate:
	rm -f $(obj) $(name) $(wildcard *.o) $(wildcard *.mod)
