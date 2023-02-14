

FC = gfortran
FCFLAGS = -O2
FCLIBS = -lblas -llapack

EXEC = hedis

PARAMS = noh-problem.f90

SOURCES = code/m_constants.f90 code/m_quantitys.f90 $(wildcard code/f_*.f90) $(wildcard code/s_*.f90) examples/$(PARAMS)
OBJECTS = $(SOURCES:.f90=.o)



$(EXEC): $(OBJECTS)
	$(FC) -o $(EXEC) $(OBJECTS) $(FCLIBS) $(FCFLAGS)


clean:
	rm -f $(EXEC) *.mod */*.o */*~ *~


allclean:
	make clean
	rm -f results/*


%.o: %.f90
	$(FC) -c $< -o $@
