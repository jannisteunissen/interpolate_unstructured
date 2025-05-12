F90 := gfortran
FFLAGS := -O2 -g -std=f2008 -Wall -Wextra
INCDIRS := kdtree2/build
LIBDIRS := kdtree2/build
LIBS := kdtree2
TESTS := test_triangle test_quad test_tetra
EXAMPLES := benchmark $(TESTS)


.PHONY:	all clean test

all: 	$(EXAMPLES)

test:	$(TESTS)
	@for test in $(TESTS); do \
		./$$test || (echo "FAILED $$test" && exit 1); \
		echo "PASSED: $$test"; \
	done

clean:
	$(RM) $(EXAMPLES) *.o *.mod

# Dependency information
$(EXAMPLES): m_interp_unstructured.o
$(EXAMPLES:%=%.o): m_interp_unstructured.o

# How to get .o object files from .f90 source files
%.o: %.f90
	$(F90) -c -o $@ $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

# How to get executables from .o object files
%: %.o
	$(F90) -o $@ $^ $(FFLAGS) $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))
