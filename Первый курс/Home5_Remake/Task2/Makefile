   comp:=gfortran
    opt:=-c -Wall -Wtabs
pattern:=*.f95
source :=$(wildcard $(pattern))
   obj :=$(patsubst %.f95, %.o, $(source))
      main : $(obj)
	     $(comp) $^ -o $@
     %.o %.mod : %.f95
	        $(comp) $(opt) $<
	        touch $@
       result : main
	        time ./$< > result
       result-c : main
	          time ./$< > result
	          cat result
        clear :
	 rm -f *.o *.mod main
        noresult : result
	 rm -f result
