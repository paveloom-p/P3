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

   input_code : input
		iconv -f utf-8 -t cp1251 input -o input_code

     input_tr : input_code
		tr абвгдежзийклмнопрстуфхцчшщъыьэюя евгджёзкойлмнпурстфыхцчшщбъэьюяа < input > input_tr

       result : result_code
		iconv -f cp1251 -t utf-8 result_code -o result

  result_code : main input_tr
	        time ./$<  < input_tr > result_code
     result-r : main input_tr
		rm result
		make result
		cat result
        clear :
	 rm -f *.o *.mod main

        clean-all :
	 rm -f *.o *.mod main *.eps *.dat result
	
