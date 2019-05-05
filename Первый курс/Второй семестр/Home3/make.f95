ext=f95
comp=gfortran
files=$(wildcard *.f95)
obj=$(patsubst %.$(ext), %.o, $(files))
name=Run

$(name): poly.f95 binrad.f95 quadra.f95 main.f95
	$(comp) poly.f95 binrad.f95 quadra.f95 main.f95 -o $@
Go:$(name)
	./$(name) < input
clean:
	rm -f $(obj)
exterminate:
	rm -f $(obj) $(name) $(wildcard *.o) $(wildcard *.mod)
