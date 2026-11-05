# Use gfortran unless already defined
F90 ?= gfortran

ifeq ($(F90), gfortran)
	FFLAGS	?= -O2 -g -Wall -Wextra
else ifeq ($(F90), ifort)
	FFLAGS	:= -O2 -stand f08 -warn all
endif

OBJS := obj/m_config.o obj/Mu_LaB_SWE.o
LIB := libconfig_fortran.a
EXECUTABLE := bin/main

.PHONY:	all test clean

all: 	$(LIB) $(EXECUTABLE)

$(LIB): $(OBJS)
	$(RM) $@
	$(AR) rcs $@ $^

clean:
	$(RM) $(EXECUTABLE) obj/m_config.o m_config.mod obj/Mu_LaB_SWE.o mu_lab_swe.mod $(LIB)

# Dependency information
$(EXECUTABLE): obj/m_config.o obj/Mu_LaB_SWE.o

# How to get .o object files from .f90 source files
obj/%.o: src/%.f90
	$(F90) -c -o $@ $< $(FFLAGS)

# How to get executables from .o object files
bin/%: obj/%.o
	$(F90) -o $@ $^ $(FFLAGS)
