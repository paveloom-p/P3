ext=f95
comp=gfortran
files=$(wildcard *.f95)
obj=$(patsubst %.$(ext), %.o, $(files))
name=Run

$(name): thanoi.f95 main.f95
	$(comp) thanoi.f95 main.f95 -o $@
Go:$(name)
	time -p ./$(name)
clean:
	rm -f $(obj)
exterminate:
	rm -f $(obj) $(name) $(wildcard *.o) $(wildcard *.mod)
