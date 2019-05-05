   comp:=gfortran
    opt:=-c -Wall -Wtabs
pattern:=*.f95
source :=$(wildcard $(pattern))
   obj :=$(patsubst %.f95, %.o, $(source))
      main : $(obj)
	     $(comp) $^ -o $@

       %.o : %.f95
		$(comp) $(opt) $< -o $@
		touch $@
     %.mod : %.f95
		$(comp) $(opt) $<
		touch $@

    main.o : subprograms.mod

       result : main input
	        time ./$< < input > result
       result-r : main input
		rm result
		make result
		cat result
       result-a : main input
		rm result
		make result
		cat result
		gnuplot script.gnu
        clear :
	 rm -f *.o *.mod main fort.*

        clean-all :
	 rm -f *.o *.mod main *.eps *.dat result
	
