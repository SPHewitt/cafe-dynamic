#$Id: makefile 2240 2017-04-18 16:36:34Z mexas $
#
# Use GNU make (gmake). Will not work with BSD make.
# Some programs will not build with Intel 15 compiler,
# so I use make conditionals based on the FC var.
# Adjust the path to CGPACK library!
#
# Based on p121 by Louise M. Lever
# -----------------------------------------------------------------------
# Compiles xx14 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build xx14 executable as RELEASE edition
#   debug: Build xx14 executable as DEBUG edition
#   install: Copy xx14 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local xx14 executable
#   relink: Rebuild xx14 from object files and libraries
# -----------------------------------------------------------------------
# Notes:
# -----------------------------------------------------------------------

include ../../../../build/mk_defs.inc
include ../../../../build/$(MACHINE).inc

MODULE=	../../../../include/mpi
LIB=	../../../../lib
BIN_DIR=../../../../bin

CGPACK=	$(HOME)/cgpack-code/lib/libcgpack.a
CGPACKMOD=$(HOME)/cgpack-code/mod
CGPACKLIB=$(HOME)/cgpack-code/lib

SRC:=	xx20std.f90 
ifeq ($(FC),ftn)
  SRC:=	$(SRC) xx20.f90 
endif

OBJ=	$(SRC:.f90=.o)
EXE=	$(SRC:.f90=.x)

.SUFFIXES: .o .f90 .x

all: release

release:
	@echo
	@echo "Building XX20 RELEASE"
	$(MAKE) $(EXE) \
	FFLAGS="-c -I$(MODULE) $(XX14_MPI_REL_FFLAGS) -I$(CGPACKMOD)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(XX14_MPI_REL_LDFLAGS) -L$(CGPACKLIB) -lcgpack"
	@echo "Done" $(EXE) "RELEASE"

debug:
	@echo
	@echo "Building xx20 DEBUG"
	$(MAKE) $(EXE) \
	FFLAGS="-c -I$(MODULE) $(XX14_MPI_DBG_FFLAGS)" -I$(CGPACK) \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(XX14_MPI_DBG_LDFLAGS)"
	@echo "Done" $(EXE) "DEBUG"

# always rebuild if the CGPACK module has updated
$(OBJ): $(CGPACK)

.f90.o:
	$(FC) $< $(FFLAGS)

.o.x:
	$(FC) $< -o $@ $(LDFLAGS)

install:
	cp $(EXE) $(BIN_DIR)

clean:
	rm -f $(OBJ)

execlean:
	rm -f $(BIN_DIR)/$(EXE)

relink: $(EXE)
