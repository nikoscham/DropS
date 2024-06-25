#*******************************************************************************
#*** DROPlet Simulation (DropS)
#*** BY NIKOLAOS CHAMAKOS (nikoscham@gmail.com)
#*** NATIONAL TECHNICAL UNIVERSITY OF ATHENS, GREECE
#*******************************************************************************

all: DropS

#*******************************************************************************
#*** DESCRIPTION: COMPILE DropS
#*******************************************************************************

#*** FORTRAN COMPILER
FCOMPILER = gfortran

#*** OPTIMIZATION FLAGS
OPTFLAGS = -O3

#*** VARIOUS FLAGS
VARFLAGS = -fbounds-check

#*** LIBMUMPS
topdir = /usr/local/MUMPS
libdir = $(topdir)/lib
include $(topdir)/Makefile.inc
LIBMUMPS_COMMON = $(libdir)/libmumps_common$(PLAT)$(LIBEXT)
LIBDMUMPS = $(libdir)/libdmumps$(PLAT)$(LIBEXT) $(LIBMUMPS_COMMON)

DropS: cvar.o auxfunct.o initsol.o nodnum.o coord.o defpar.o tsfun.o eikon.solv.o printres.o ew.solv.o ew.paramsolv.o ayl.paramsolv.o ayl.solv.o cyl.solv.o DropS.o 
	$(FCOMPILER) cvar.o auxfunct.o initsol.o  nodnum.o coord.o defpar.o tsfun.o eikon.solv.o printres.o ew.solv.o ew.paramsolv.o ayl.paramsolv.o ayl.solv.o cyl.solv.o DropS.o $(INCS) $(LIBDMUMPS) $(LORDERINGS) $(LIBS) $(LIBBLAS) $(LIBOTHERS) -I$(topdir)/include $(VARFLAGS) $(OPTFLAGS) -o DropS

cvar.o: cvar.f90
	$(FCOMPILER)  -c $(OPTFLAGS) cvar.f90
coord.o: coord.f90
	$(FCOMPILER)  -c $(OPTFLAGS) coord.f90
defpar.o: defpar.f90
	$(FCOMPILER)  -c $(OPTFLAGS) defpar.f90
ew.solv.o: ew.solv.f90
	$(FCOMPILER)  -c $(OPTFLAGS) ew.solv.f90
ew.paramsolv.o: ew.paramsolv.f90
	$(FCOMPILER)  -c $(OPTFLAGS) ew.paramsolv.f90
ayl.paramsolv.o: ayl.paramsolv.f90
	$(FCOMPILER)  -c $(OPTFLAGS) ayl.paramsolv.f90
ayl.solv.o: ayl.solv.f90
	$(FCOMPILER)  -c $(OPTFLAGS) ayl.solv.f90
cyl.solv.o: cyl.solv.f90
	$(FCOMPILER)  -c $(OPTFLAGS) cyl.solv.f90
initsol.o: initsol.f90
	$(FCOMPILER)  -c $(OPTFLAGS) initsol.f90
DropS.o: DropS.f90
	$(FCOMPILER)  -c $(OPTFLAGS) DropS.f90
nodnum.o: nodnum.f90
	$(FCOMPILER)  -c $(OPTFLAGS) nodnum.f90
tsfun.o: tsfun.f90
	$(FCOMPILER)  -c $(OPTFLAGS) tsfun.f90
eikon.solv.o: eikon.solv.f90
	$(FCOMPILER)  -c $(OPTFLAGS) eikon.solv.f90
auxfunct.o: auxfunct.f90
	$(FCOMPILER)  -c $(OPTFLAGS) auxfunct.f90 -I$(topdir)/include
printres.o: printres.f90
	$(FCOMPILER)  -c $(OPTFLAGS) printres.f90

clean:
	$(RM) *.o *.mod DropS

#*******************************************************************************
