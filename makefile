FC = ifort
#FC = gfortran

# for debug 
#FFLAGS = -fno-backtrace -Wall -g
#FFLAGS = -check all -fpe0 -traceback -qopenmp -g

# for performance
#FFLAGS = -O3 -traceback -qopenmp -g
FFLAGS = -O3 -traceback -g

LIB = -mkl
BINDIR = ../bin

MODULES = const.o \
          thermostat.o \
          pes_hd.o \
          ctrl.o

SUBROUTINES = initialize.o \
              NHC_VV.o \
              MD_VV.o  \
              pes_interface.o \
              simulation.o \
	      pes_ruFIXco.o

OBJECTS = $(MODULES) main.o $(SUBROUTINES)

MODS = $(MODULES:.o=.mod)

#@echo $(MODULES)
#@echo $(SUBROUTINES)

all:	xmd.x pmf2d.x gen.x

xmd.x: $(OBJECTS)
	$(FC) $(FFLAGS) $(LIB) -o xmd.x $(OBJECTS)

gen.x:	gen.f90
	$(FC) $(FFLAGS) $(LIB) $+ -o $@

pmf2d.x: pmf2d.f90
	$(FC) $(FFLAGS) $(LIB) $+ -o $@

## pmf2d.x 300 traj_co_300.dat PMF2 0 1 50 0 1 50

$(MODULES): %.o: %.f90
	$(FC) -c $(FFLAGS) $<

main.o: main.f90
	$(FC) -c $(FFLAGS) main.f90

$(SUBROUTINES): %.o: %.f90
	$(FC) -c $(FFLAGS) $<

# Modify the clean target to remove gen.x and pmf2d.x as well

clean:
	rm  *.mod *.o xmd.x gen.x pmf2d.x

install:
	install -d $(BINDIR)
	install -m 0755 xmd.x gen.x pmf2d.x $(BINDIR) # Install all binaries

uninstall:
	rm  $(BINDIR)/xmd.x $(BINDIR)/gen.x $(BINDIR)/pmf2d.x # Remove all installed binaries
