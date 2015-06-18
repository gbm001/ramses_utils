.SUFFIXES:

export

F90            = gfortran
MPIF90         = mpif90    # Only if using MPI
DEBUG          = 2         # Either 0, 1 or 2
OPT            = FAST      # Either 0, 1, 2, 3 or FAST
OPENMP         = 0         # Either 0 or 1
MPI            = 0         # Either 0 or 1

# List of executables to be built within the package
PROGRAMS = amr_utils.o convert_to_single dump_level_cells velocity_dispersion \
           power_spectrum
TEST_PROGRAMS = test_binning

# Compilation directories
SRCDIR = ${PWD}
EXEDIR = ${PWD}

VPATH = ${SRCDIR}
MODULEDIR = ${SRCDIR}
# HEADERS = ${SRCDIR}/headers
INCLUDE_DIR += ${MODULEDIR} #${HEADERS}

FFLAGS += $(foreach INCDIR, ${INCLUDE_DIR}, -I ${INCDIR})
LDFLAGS +=
LDOBJ += libfftw3.a

# Remove trailing whitespace
# ----------------------------------------------------------------
F90           := $(strip ${F90})
MPIF90        := $(strip ${MPIF90})
DEBUG         := $(strip $(DEBUG))
OPT           := $(strip $(OPT))
OPENMP        := $(strip $(OPENMP))
MPI           := $(strip $(MPI))

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

# Standards flags
ifeq (${F90},ifort)
    #FFLAGS += -fpscomp logicals
    FFLAGS += -assume noold_ldout_format,noold_maxminloc,noold_unit_star,noold_xor,std_mod_proc_name,fpe_summary,byterecl
                    # minus0,protect_parens,realloc_lhs
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
        FFLAGS += -ipo -O3 -no-prec-div -static-intel -xHost #-fast
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

.PHONY: all clean fileclean depend distclean programs tests ${PROGRAMS}
.PRECIOUS: %.o
.DEFAULT_GOAL := programs

all: programs tests

programs: depend ${PROGRAMS}

tests: depend ${TEST_PROGRAMS}

${PROGRAMS}: depend
	@${MAKE} --no-print-directory -f makefile.inc $@

clean: fileclean depend

fileclean:
	\rm -f *.o
	\rm -f *.mod
	\rm -f *.f90~
	\rm -f *.F90~

distclean: fileclean
	\rm -f ${PROGRAMS}
	\rm -f fort.dep

depend:
	python3 make_depends.py ${VPATH}

Makefile: ;
