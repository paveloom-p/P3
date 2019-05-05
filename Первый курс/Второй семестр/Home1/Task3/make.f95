ext=f95
comp=gfortran
files=$(wildcard *.f95)
obj=$(patsubst %.$(ext), %.o, $(files))
name=Run

$(name): functions.f95 main.f95
	$(comp) functions.f95 main.f95 -o $@
Go:$(name)
	time -p ./$(name)
clean:
	rm -f $(obj)
exterminate:
	rm -f $(obj) $(name) $(wildcard *.o) $(wildcard *.mod)
