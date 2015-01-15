.SUFFIXES:
.SUFFIXES: .c .f .f90 .F90 .o

F90            = ifort
DEBUG          = 1
PARALLEL       = 0

# Remove trailing whitespace
# ----------------------------------------------------------------
F90           := $(strip $(F90))
PROFILE       := $(strip $(PROFILE))

# Compiler flags
# ----------------------------------------------------------------
ifeq ($(F90),ifort)
OPT += -O3 -static-intel -ip -xHOST
ifeq ($(PARALLEL),1)
OPT += -openmp 
endif
ifeq ($(DEBUG),1)
OPT += -g -fpe0 -traceback -check all -check noarg_temp_created -debug extended -debug-parameters all
OPT += -warn all -warn nounused -debug all
endif
endif

# Code compilation
# ----------------------------------------------------------------
# List of executables to be built within the package
PROGRAMS = convert_to_single dump_level_cells

.PHONY: all clean distclean

all: $(PROGRAMS)

amr_utils.o: amr_utils.f90
	$(F90) $(OPT) -c amr_utils.f90

%: %.o amr_utils.o
	$(F90) $(OPT) -o $@ $^

%.o: %.F90 amr_utils.o
	$(F90) $(OPT) -c $<

%.o: %.f90 amr_utils.o
	$(F90) $(OPT) -c $<

clean:
	\rm -f *.o
	\rm -f *.mod
	\rm -f *.F90~

distclean: clean
	\rm -f $(PROGRAMS)

