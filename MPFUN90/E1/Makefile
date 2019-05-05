##############################################################################
#  example/Makefile
#
#  This software is provided for research purposes only.
#  Commercial usage requires license agreement.
#
#  This work was supported by the Director, Office of Science, Division
#  of Mathematical, Information, and Computational Sciences of the
#  U.S. Department of Energy under contract number DE-AC03-76SF00098.
# 
#  Copyright (c) 2001-2002-2003-2004
#  
#  Modified:       2008-08-25
#  Version:        1.0
#
##############################################################################

###### Fortran Test Programs ######
F90 = gfortran

F90FLAGS = -O2 -ffree-form

LD = gfortran

F90MOD = p_mpfun_e1.mod p_e1.mod
F90O =  mpfun90.o mpmod90.o p_main.o second.o p_mpfun_e1.o p_e1.o

main:  $(F90O)
	$(LD) -o $@ $(F90O)

###### General Rules ######

p_main.o : $(F90MOD)

.f.o:
	$(F90) $(F90FLAGS) -c $<

%.mod : %.f
	$(F90) -c $(F90FLAGS) $<

result: main input
	touch mp_output_copy
	rm mp_output_copy
	./main > result < input

result-r: input
	make result
	cat result

###### Clean Up ######

clear: 
	rm -f *.o *.mod *~ *.d work.pc work.pcl fort.* main


