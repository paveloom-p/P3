   comp:=gfortran
    opt:=-c -Wall -Wtabs
pattern:=*.f95
source :=$(wildcard $(pattern))
   obj :=$(patsubst %.f95, %.o, $(source))
      main : $(obj)
	     $(comp) $^ -o $@
       %.o : %.f95
		$(comp) $(opt) $< -o $@
     %.mod : %.f95
		$(comp) $(opt) $<

       main.o : real_type.mod
     quadra.o : real_type.mod
  functions.o : real_type.mod

       result : main input
	        time ./$<  < input > result
       result-r : main input
		rm result
		make result
		cat result
        clear :
	 rm -f *.o *.mod main

        clean-all :
	 rm -f *.o *.mod main *.eps *.dat result
	
