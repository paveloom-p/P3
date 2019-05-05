SOURCE:=$(wildcard *.f)
OBJECT:=$(patsubst %.f, %.o, $(SOURCE))
main : $(OBJECT) 
	 gfortran $^ -o $@
%.o  :  %.f
	gfortran -c $^
result  : main
	./main  > result
result-r  : main
	make result
	cat result
clear :
	 rm -f *.o *.mod main
