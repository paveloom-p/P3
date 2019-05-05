# vpath %.f /home/aw/nagaw8/my_module/
   comp:=gfortran
    opt:=-c -Wall -Wtabs
  mymod:=poly.o binrad.o quadra.o
  input:=input
pattern:=*.f95
source :=$(wildcard $(pattern))
   obj :=$(patsubst %.f95, %.o, $(source))
   main:=main
      $(main) : $(mymod) $(obj)
	        $(comp) $^ -o $@
     %.o %.mod : %.f95
	        $(comp) $(opt) $<
	        touch $@
    $(main).o : poly.mod  binrad.mod quadra.mod
#   my_const.o : my_prec.mod
#    my_besk.o : my_prec.mod my_const.mod
       result : $(main) $(input)
	        time ./$<  < $(input) #> result
        clear :
	 rm -f *.o *.mod $(main)
