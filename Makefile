F90 := gfortran
FFLAGS := -O2 -g -std=f2008 -Wall -Wextra
ifeq ($(DEBUG), 1)
	FFLAGS += -O0 -fcheck=all
endif
TESTS := test_triangle test_quad test_tetra test_vtk test_trace_field
EXAMPLES := benchmark $(TESTS)
LIB := libinterp_unstructured.a
OBJECTS := m_interp_unstructured.o kdtree2_module.o m_binda.o m_vtk.o

# So that kdtree2_module.f90 can be found
vpath %.f90 kdtree2/src

.PHONY:	all clean test

all: 	$(EXAMPLES) $(LIB)

test:	$(TESTS)
	@for test in $(TESTS); do \
		./$$test || (echo "FAILED $$test" && exit 1); \
		echo "PASSED: $$test"; \
	done

clean:
	$(RM) $(EXAMPLES) *.o *.mod

$(LIB): $(OBJECTS)
	$(AR) rcs $@ $^

# Dependency information
$(EXAMPLES): $(OBJECTS)
$(EXAMPLES:%=%.o): $(OBJECTS)
m_interp_unstructured.o: kdtree2_module.o m_binda.o m_vtk.o

# How to get .o object files from .f90 source files
%.o: %.f90
	$(F90) -c -o $@ $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

# How to get executables from .o object files
%: %.o
	$(F90) -o $@ $^ $(FFLAGS) $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))
