.PHONY: all
all: main.out oucd.out

#source=main
#object=a.out

#${object}: ${source}.f90
#	gfortran -g ${source} -o ${object} -Wall -Wtabs -Wextra -Wconversion -fbacktrace -fbounds-check -ffpe-trap=zero,overflow,underflow -fcheck=all -fmax-errors=7
#	gfortran -g ${source}.f90 -o ${object} -Wall -Wno-tabs -Wextra -Wconversion -fbacktrace -fbounds-check -ffpe-trap=zero,overflow,underflow -fcheck=all -fmax-errors=1 #
#	gfortran ${source} -o ${object} -Ofast #PRODUCTION

%.out: %.f90
	gfortran standard.f90 $< -o $@ -Ofast #PRODUCTION
	#gfortran -g $< -o $@ -Wall -Wno-tabs -Wextra -Wconversion -fbacktrace -fbounds-check -ffpe-trap=zero,overflow,underflow -fcheck=all -fmax-errors=1

.PHONY: clean
clean:
	rm -f *.mod

.PHONY: erase
erase: clean
	rm -f *.out

