   comp:=gfortran
    opt:=-c -Wall -Wtabs
pattern:=*.f95
source :=$(wildcard $(pattern))
   obj :=$(patsubst %.f95, %.o, $(source))
      main : $(obj)
	     $(comp) $^ -o $@
     %.o %.mod : %.f95
	        $(comp) $(opt) $< -o $@
	        touch $@
       result : main
	        time ./$< > result
       result-r : main
		rm result
		make result
		cat result
       result-a : main
		rm result
		make result
		cat result
		gnuplot dist.plt; gv distrib.eps
        clear :
	 rm -f *.o *.mod main

        clean-all :
	 rm -f *.o *.mod main *.eps *.dat result
	
