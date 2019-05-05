ext=f95
comp=gfortran
files=$(wildcard *.f95)
obj=$(patsubst %.$(ext), %.o, $(files))
name=Run

$(name): S4_degree.f95 S5_less_than_1_binary.f95 S6_simple_number.f95 main.f95
	$(comp) S4_degree.f95 S5_less_than_1_binary.f95 S6_simple_number.f95 main.f95 -o $@
Go:$(name)
	time -p ./$(name) < input
clean:
	rm -f $(obj)
exterminate:
	rm -f $(obj) $(name) $(wildcard *.o) $(wildcard *.mod)
