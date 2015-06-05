.SUFFIXES:

F90            = gfortran
MPIF90         = mpif90    # Only if using MPI
DEBUG          = 2         # Either 0, 1 or 2
OPT            = 0         # Either 0, 1, 2, 3 or FAST
OPENMP         = 0         # Either 0 or 1
MPI            = 0         # Either 0 or 1

# List of executables to be built within the package
PROGRAMS = amr_utils.o convert_to_single dump_level_cells velocity_dispersion

# Compilation directories
SRCDIR = ${PWD}
EXEDIR = ${PWD}

VPATH = ${SRCDIR}
MODULEDIR = ${SRCDIR}
# HEADERS = ${SRCDIR}/headers
INCLUDE_DIR += ${MODULEDIR} #${HEADERS}

FFLAGS += $(foreach INCDIR, ${INCLUDE_DIR}, -I ${INCDIR})

# Remove trailing whitespace
# ----------------------------------------------------------------
F90           := $(strip ${F90})
PROFILE       := $(strip $(PROFILE))

# Compiler flags
# ----------------------------------------------------------------

# Debug flags
ifeq (${DEBUG},1)
    ifeq (${F90},ifort)
        FFLAGS += -g -fpe0 -traceback -debug extended -debug-parameters all -warn all -warn nounused -debug all
    else ifeq (${F90},gfortran)
        FFLAGS += -fimplicit-none -fbacktrace -g -pedantic -Wall -Wextra
    endif
    FFLAGS += -DDEBUG
else ifeq (${DEBUG},2)
    ifeq (${F90},ifort)
        FFLAGS += -g -fpe0 -traceback -check all -check noarg_temp_created -debug extended -debug-parameters all
        FFLAGS += -warn all -warn nounused -debug all
    else ifeq (${F90},gfortran)
        FFLAGS += -fimplicit-none -fbounds-check -fbacktrace -g
        FFLAGS += -fcheck=all -pedantic -Wall -Wextra -ffpe-trap=invalid,zero,overflow,underflow
    endif
    FFLAGS += -DDEBUG
endif

# Optimisation flags
ifeq (${OPT},0)
    FFLAGS += -O0
else ifeq (${OPT},1)
    FFLAGS += -O1
else ifeq (${OPT},2)
    FFLAGS += -O2
else ifeq (${OPT},3)
    FFLAGS += -O3
else ifeq (${OPT},FAST)
    ifeq (${F90},ifort)
        FFLAGS += -fast #-O3 -static-intel -ip -ipo -xHOST
    else ifeq (${F90},gfortran)
        FFLAGS += -Ofast -march=native
    else
        FFLAGS += -O3
    endif
endif

# Parallel flags
ifeq (${PARALLEL},1)
    ifeq (${F90},ifort)
        FFLAGS += -openmp -parallel
    else ifeq (${F90},gfortran)
        FFLAGS += -fopenmp
    endif
endif

# Code compilation
# ----------------------------------------------------------------

.PHONY: all clean depend distclean
.PRECIOUS: %.o
.DEFAULT: all

all: fort.dep ${PROGRAMS}

%: %.o
	${F90} ${FFLAGS} ${LDFLAGS} -o ${EXEDIR}/$@ $^

%.o: %.F90
	${F90} ${FFLAGS} -c $<

%.o: %.f90
	${F90} ${FFLAGS} -c $<

clean: depend
	\rm -f *.o
	\rm -f *.mod
	\rm -f *.F90~

distclean: clean
	\rm -f ${PROGRAMS}
	rm fort.dep

depend:
	python3 make_depends.py ${VPATH}

# Dependencies
# ----------------------------------------------------------------

fort.dep : depend

-include fort.dep
