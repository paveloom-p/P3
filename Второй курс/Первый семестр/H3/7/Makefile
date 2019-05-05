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

 input_cp1251 : input
		iconv -f utf-8 -t cp1251 input -o input_cp1251 
result_cp1251 : main input_cp1251
	        time ./$<  < input_cp1251 > result_cp1251
       result : main result_cp1251
	        iconv -f cp1251 -t utf-8 result_cp1251 -o result
       result-r : main input_cp1251
		rm result
		make result
		cat result
        clear :
	 rm -f *.o *.mod main

        clean-all :
	 rm -f *.o *.mod main *.eps *.dat result
	
