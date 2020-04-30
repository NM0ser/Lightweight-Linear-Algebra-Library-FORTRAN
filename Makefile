run_verification: main_verification.o nmoser_linear_algebra.o
	gfortran -o run_verification main_verification.o nmoser_linear_algebra.o

main_verification.o: main_verification.f nmoser_linear_algebra.o
	gfortran -c main_verification.f

nmoser_linear_algebra.o: nmoser_linear_algebra.f
	gfortran -c nmoser_linear_algebra.f

clean:
	-rm *.o *.mod *.txt run_verification