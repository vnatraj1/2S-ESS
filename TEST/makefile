#  PCAL2S tool.

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#               UV HISTORY

#  This Makefile, constructed 10/25/2011 by R. Spurr and V. Natraj
#  This Makefile, Updated     01/16/2012 by R. Spurr and V. Natraj
#  This Makefile, Updated     04/16/2012 by R. Spurr and V. Natraj
#  This Makefile, Updated     06/04/2012 by R. Spurr and V. Natraj

#  This Makefile, Updated     11/30/2012 by R. Spurr
#  This Makefile, Updated     04/04/2013 by R. Spurr
#  This Makefile, Updated     08/06/2013 by R. Spurr

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#  This Makefile, Updated     05/09/2015 by R. Spurr
#     - New versions of LIDORT and 2STREAM
#     - Re-engineered codes, OpenMP options

#  THIS IS THE UVVSWIR SECTION

# Path variables
# ==============

#  General sourcecode

GEN_PATH      =  ../sourcecode
2S_PATH       = $(GEN_PATH)/2stream_2p4_feb15
FO_PATH_1P4   = $(GEN_PATH)/FO_1p4_sourcecode

#  UVVSWIR Specific paths

UVVSWIR_TESTER_PATH    = .

MOD_PATH = mod
OBJ_PATH = obj

MOD_FILES = $(MOD_PATH)/*.mod
OBJ_FILES = $(OBJ_PATH)/*.o
SU_FILES  = $(OBJ_PATH)/*.su
GENMOD_FILES = $(MOD_PATH)/*genmod.f90

# Compiler options
# ================

#  Intel
#FC = ifort

#  gfortran
FC = gfortran

#  g95
#FC = g95

#  NAG
#FC = f95

# Define FORTRAN90 compiler flags
# ===============================

FFLAGS = -c

# Additional flags for Intel
ifeq ($(FC), ifort)
	FFLAGS := $(FFLAGS) -mcmodel=large -fPIC -i-dynamic -I$(MOD_PATH) -module $(MOD_PATH)
	FFLAGS_DEBUG = -g -warn all -check all -traceback -zero
	FFLAGS_OPT = -O3
	FFLAGS_OPENMP = -openmp
endif

# Additional flags for gfortran
ifeq ($(FC), gfortran)
	FFLAGS := $(FFLAGS) -mcmodel=large -I$(MOD_PATH) -J$(MOD_PATH) -O3
#	FFLAGS_DEBUG = -g -C -Wall -mcmodel=large -fbounds-check -frange-check -ffpe-trap=invalid,zero,overflow -fbacktrace
#	FFLAGS := $(FFLAGS) -I$(MOD_PATH) -J$(MOD_PATH)

#	FFLAGS_DEBUG = -fstack-usage
	FFLAGS_OPT = -O3
        FFLAGS_OPENMP = -fopenmp -frecursive
endif

# Additional flags for g95
ifeq ($(FC), g95)
#      older g95
	FFLAGS := $(FFLAGS) -I$(MOD_PATH) -fmod=$(MOD_PATH)
	FFLAGS_DEBUG = -g -O2 -Wall -Wuninitialized -fbounds-check
	FFLAGS_OPT = -O3
#      g95 v0.92
#	FFLAGS := $(FFLAGS) -I$(MOD_PATH) -fmod=$(MOD_PATH)
#	FFLAGS_DEBUG = -g -Wall -fbounds-check -ftrace=full
endif

# Additional flags for NAG
#ifeq ($(FC), f95)
ifeq ($(FC), $(filter $(FC), f95 nagfor))
	FFLAGS := $(FFLAGS) -mdir $(MOD_PATH) -I$(MOD_PATH)
#	FFLAGS_DEBUG = -g -C=all -C=undefined -gline -mtrace=all -nan
	FFLAGS_DEBUG = -w=obs -w=unused -C=all -C=undefined -gline
#	FFLAGS_DEBUG = -w=obs -w=unused -g -C=all -C=undefined -gline -mtrace=all -nan
	FFLAGS_OPT = -O3
        FFLAGS_OPENMP = -openmp
endif

# Additional flags for PGF90
ifeq ($(FC), pgf90)
#   PLACEHOLDER
endif


# For debug build, use "make DEBUG=t"
ifeq ($(DEBUG), t)
	FFLAGS := $(FFLAGS) $(FFLAGS_DEBUG)
endif

# For optimized build, use "make OPT=t"
ifeq ($(OPT), t)
	FFLAGS := $(FFLAGS) $(FFLAGS_OPT)
endif

# For parallel build using OpenMP, use "make OPENMP=t"
ifeq ($(OPENMP), t)
	FFLAGS := $(FFLAGS) $(FFLAGS_OPENMP)
endif


#
# Define list of source files
# %%%%%%%%%%%%%%%%%%%%%%%%%%%

# (Note: ordering is important because of dependencies)

#  2STREAM SOURCES
#  ===============

# Regular

# new version 2.4
SOURCES_2S =
SOURCES_2S +=   \
   $(2S_PATH)/2stream_Taylor.f90 \
   $(2S_PATH)/2stream_inputs.f90 \
   $(2S_PATH)/2stream_writemodules.f90 \
   $(2S_PATH)/2stream_miscsetups.f90 \
   $(2S_PATH)/2stream_solutions.f90  \
   $(2S_PATH)/2stream_bvproblem.f90  \
   $(2S_PATH)/2stream_thermalsup.f90 \
   $(2S_PATH)/2stream_intensity.f90

#  masters

SOURCES_2S_MASTER =
SOURCES_2S_MASTER +=   \
   $(SOURCES_2S)      \
   $(2S_PATH)/2stream_master.f90

#  FIRST_ORDER SOURCES
#  ===================

#  Version 1.4
#  ===========

#  Geometry and Regular

SOURCES_FO_1P4_G =
SOURCES_FO_1P4_G +=   \
   $(FO_PATH_1P4)/FO_geometry_Generic.f90  \
   $(FO_PATH_1P4)/FO_geometry_Routines.f90 \
   $(FO_PATH_1P4)/FO_SSgeometry_master.f90 \
   $(FO_PATH_1P4)/FO_ScalarSS_Spherfuncs.f90

SOURCES_FO_1P4_MASTER =
SOURCES_FO_1P4_MASTER +=   \
   $(SOURCES_FO_1P4_G)      \
   $(FO_PATH_1P4)/FO_ScalarSS_RTCalcs_I.f90

#  EXACT MASTERS
#  =============

SOURCES_UVVSWIR_EXACT =
SOURCES_UVVSWIR_EXACT +=     \
    $(SOURCES_2S_MASTER)     \
    $(SOURCES_FO_1P4_MASTER) \
    $(UVVSWIR_TESTER_PATH)/test_regular.f90

# Define pattern rules for creating object files:
# ===============================================

.SUFFIXES:

# For f90 source files

$(OBJ_PATH)/%.o : $(2S_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(FO_PATH_1P4)/%.f90
	$(FC) $(FFLAGS) $< -o $@

$(OBJ_PATH)/%.o : $(UVVSWIR_TESTER_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@

# Define object files
# ===================

#  Exact

F90SOURCES_UVVSWIR_EXACT := $(notdir $(filter %.f90, $(SOURCES_UVVSWIR_EXACT)))
F90OBJECTS_UVVSWIR_EXACT := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_UVVSWIR_EXACT)))

# Define desired targets
# ======================

APPLICATION_UVVSWIR_EXACT        = test_regular.exe

Main:  $(APPLICATION_UVVSWIR_EXACT)

test_regular.exe: $(F90OBJECTS_UVVSWIR_EXACT)
	$(FC) $^ -o $@

#  make clean
#  ==========

.PHONY: clean
clean:
	rm -f *.o $(OBJ_FILES) $(SU_FILES) *.mod $(MOD_FILES) $(GENMOD_FILES) *.log *.exe

